#!/usr/bin/env python3
r"""esm_embed.py — ESM-2 蛋白质嵌入提取(残基级 / 序列级)

从 FASTA 提取 ESM-2 表征:残基级 embedding(L×D)、可选 mean-pool 序列级 embedding。
下游用途:序列聚类、突变效应、口袋特征、机器学习输入、嵌入可视化。

运行环境(torch + fair-esm 独立 venv,CPU 版):
  & 'D:\bioai\venv-esm\Scripts\python.exe' esm_embed.py --fasta query.fasta \
      --model esm2_t6_8M_UR50D --repr-layers 6 --outdir D:\bioai\jobs\esm_embed
模型(dl.fbaipublicfiles.com 自动下载到 TORCH_HOME):
  esm2_t6_8M_UR50D   (34MB,  8M 参数,  D=320  — 默认,快速)
  esm2_t12_35M_UR50D (140MB, 35M,      D=480)
  esm2_t30_150M_UR50D(600MB, 150M,     D=640)
  esm2_t33_650M_UR50D(2.5GB, 650M,     D=1280 — 科研级,下载大)
  esm2_t36_3B_UR50D / esm2_t48_15B_UR50D (显存/内存要求高,CPU 慎用)

输出(--outdir):
  <id>_residue_<layer>.npy   残基级 (L, D)
  <id>_sequence_<layer>.npy  序列级 (D,)(--mean-pool 时)
  metadata.json              模型/层/序列/形状

编程调用:
  from esm_embed import embed_fasta
  result = embed_fasta("query.fasta", model="esm2_t6_8M_UR50D", repr_layers=[6])
"""
import argparse
import json
import os
import sys

import numpy as np


def _load_model(name):
    import esm
    return esm.pretrained.load_model_and_alphabet(name)


def embed_sequences(seqs, model="esm2_t6_8M_UR50D", repr_layers=None, mean_pool=True):
    """seqs: [(id, seq)] → dict {id: {layer_l: {"residue": np.array, "sequence": np.array|None}}}"""
    import torch
    import esm
    torch.set_grad_enabled(False)
    model_fn, alphabet = _load_model(model)
    model_fn.eval()
    batch_converter = alphabet.get_batch_converter()
    data = [("q%d" % i, s) for i, (_, s) in enumerate(seqs)]
    _, _, toks = batch_converter(data)
    repr_layers = [(repr_layers or [model_fn.num_layers])[-1]] \
        if isinstance(repr_layers, int) else list(repr_layers or [model_fn.num_layers])
    out = model_fn(toks, repr_layers=repr_layers, return_contacts=False)
    reps = out["representations"]
    result = {}
    for i, (sid, seq) in enumerate(seqs):
        entry = {}
        for layer in repr_layers:
            t = reps[layer][i]
            residue = t[1:1 + len(seq)].detach().cpu().numpy().astype(np.float32)
            entry[layer] = {
                "residue": residue,
                "sequence": residue.mean(0) if mean_pool else None,
            }
        result[sid] = entry
    return result, model


def embed_fasta(fasta_path, model="esm2_t6_8M_UR50D", repr_layers=None,
                mean_pool=True, outdir=None):
    seqs = []
    with open(fasta_path, encoding="utf-8") as f:
        sid, buf = None, []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None:
                    seqs.append((sid, "".join(buf)))
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
        if sid is not None:
            seqs.append((sid, "".join(buf)))
    seqs = [(sid, s.replace("*", "")) for sid, s in seqs if s]
    # 复合物 fasta(colabfold 格式:链用 ':' 连接)→ 按链拆分为独立序列
    split_seqs = []
    for sid, s in seqs:
        parts = [p for p in s.split(":") if p]
        if len(parts) > 1:
            for i, p in enumerate(parts):
                split_seqs.append(("%s_chain%d" % (sid, i + 1), p))
        else:
            split_seqs.append((sid, parts[0] if parts else ""))
    seqs = [x for x in split_seqs if x[1]]
    if not seqs:
        raise ValueError("no sequences parsed from %s" % fasta_path)
    result, model_name = embed_sequences(seqs, model=model, repr_layers=repr_layers,
                                         mean_pool=mean_pool)
    meta = {"model": model_name, "mean_pool": mean_pool, "records": []}
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        for sid, entry in result.items():
            rec = {"id": sid, "seq_len": None, "layers": {}}
            for layer, tensors in entry.items():
                resid = "%s_residue_L%d.npy" % (sid, layer)
                np.save(os.path.join(outdir, resid), tensors["residue"])
                layer_info = {"layer": layer, "residue_npy": resid,
                              "shape": list(tensors["residue"].shape)}
                if tensors["sequence"] is not None:
                    seqnpy = "%s_sequence_L%d.npy" % (sid, layer)
                    np.save(os.path.join(outdir, seqnpy), tensors["sequence"])
                    layer_info["sequence_npy"] = seqnpy
                rec["layers"][str(layer)] = layer_info
                rec["seq_len"] = len(tensors["residue"])
            meta["records"].append(rec)
        with open(os.path.join(outdir, "metadata.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f, ensure_ascii=False, indent=2)
    return result, meta


def main(argv=None):
    ap = argparse.ArgumentParser(description="ESM-2 embedding extractor (fair-esm)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--model", default="esm2_t6_8M_UR50D")
    ap.add_argument("--repr-layers", default="last",
                    help="comma-separated layer indices or 'last' (default)")
    ap.add_argument("--no-mean-pool", action="store_true")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)
    try:
        if args.repr_layers.strip().lower() == "last":
            layers = None
        else:
            layers = [int(x) for x in args.repr_layers.split(",")]
        result, meta = embed_fasta(args.fasta, model=args.model, repr_layers=layers,
                                   mean_pool=not args.no_mean_pool, outdir=args.outdir)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    for rec in meta["records"]:
        print("%s: len=%s  layers=%s" % (
            rec["id"], rec["seq_len"],
            ", ".join("%s(%s)" % (k, "x".join(map(str, v["shape"]))) for k, v in rec["layers"].items())))
    print("metadata: %s" % os.path.join(args.outdir, "metadata.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
