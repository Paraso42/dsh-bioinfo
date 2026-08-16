#!/usr/bin/env python3
r"""pp_interact.py — 蛋白-蛋白互作界面分析(接触对 / 界面残基 / 埋藏表面积 BSA)

纯 Biopython(Bio.PDB,1.87 实测 API),生信环境运行:
  $env:PYTHONPATH='D:\biopython'
  & 'C:\Program Files\Python313\python.exe' pp_interact.py --complex complex.pdb --chains A B --out interface.json

输出:
  1. 链间原子接触对(NeighborSearch,默认 5 Å 截断,残基对去重取最小距离)
  2. 界面残基(参与接触的残基;默认排除 HETATM,可用 --include-het 保留)
  3. BSA 埋藏表面积:ΔSASA = SASA(A 单独) + SASA(B 单独) - SASA(复合物)
     (ShrakeRupley,逐残基;阈值 --bsa-min,默认 1.0 Å²)
  4. JSON 报告 + 人类可读摘要

编程调用:
  from pp_interact import analyze_complex
  report = analyze_complex("complex.pdb", chains=("A","B"), cutoff=5.0)
"""
import argparse
import json
import math
import os
import sys
import tempfile
from collections import defaultdict

DEFAULT_CUTOFF = 5.0          # 接触截断(Å)
DEFAULT_BSA_MIN = 1.0         # 界面残基 BSA 阈值(Å²)


def _load(path, model_id=0):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("complex", path)
    return struct[model_id]


def _atoms(model, chain_id, include_het=False):
    chain = model[chain_id]
    atoms = []
    for res in chain:
        if not include_het and res.id[0] != " ":   # id[0]==' ' 为标准残基
            continue
        atoms.extend(res.get_atoms())
    return atoms


def _res_label(res):
    return {
        "chain": res.get_parent().id,
        "resid": res.id[1],
        "resname": res.get_resname().strip(),
    }


def _compute_sasa_map(model, include_het=False):
    from Bio.PDB.SASA import ShrakeRupley
    ShrakeRupley().compute(model, level="R")
    m = {}
    for chain in model:
        for res in chain:
            if not include_het and res.id[0] != " ":
                continue
            m[(chain.id, res.id[1])] = float(res.sasa)
    return m


def _extract_chain_sasa(struct, chain_id, tmpdir, include_het=False):
    """把单一链写出临时 PDB 再解析,得到隔离状态 SASA。"""
    from Bio.PDB import PDBIO, Select

    class ChainOnly(Select):
        def accept_chain(self, chain):
            return chain.id == chain_id

    io = PDBIO()
    io.set_structure(struct)
    out = os.path.join(tmpdir, "chain_%s.pdb" % chain_id)
    io.save(out, ChainOnly())
    model = _load(out)
    return _compute_sasa_map(model, include_het=include_het)


def analyze_complex(path, chains=("A", "B"), cutoff=DEFAULT_CUTOFF,
                    include_het=False, bsa_min=DEFAULT_BSA_MIN):
    """分析复合物 PDB 的两条链之间的互作界面。返回 dict 报告。"""
    from Bio.PDB import NeighborSearch

    model = _load(path)
    for cid in chains:
        if cid not in model:
            raise KeyError("chain %r not found in %s" % (cid, path))

    ca, cb = chains
    atoms_a = _atoms(model, ca, include_het)
    atoms_b = _atoms(model, cb, include_het)
    if not atoms_a or not atoms_b:
        raise ValueError("chain %s or %s has no atoms (HETATM excluded? try --include-het)" % (ca, cb))

    # 1) 原子接触 -> 残基对去重(取最小距离)
    ns = NeighborSearch(atoms_b)
    best = {}  # (label_a, label_b) -> min dist
    for a in atoms_a:
        for b in ns.search(a.coord, cutoff):
            la, lb = _res_label(a.get_parent()), _res_label(b.get_parent())
            key = (la["chain"], la["resid"], lb["chain"], lb["resid"])
            d = float(a - b)
            if key not in best or d < best[key][0]:
                best[key] = (d, la, lb)

    contacts = [
        {"a": la, "b": lb, "dist": round(d, 2)}
        for (_, _, _, _), (d, la, lb) in sorted(best.items(), key=lambda kv: kv[1][0])
    ]
    interface = defaultdict(list)
    for c in contacts:
        interface[c["a"]["chain"]].append(c["a"])
        interface[c["b"]["chain"]].append(c["b"])
    for cid in interface:
        seen = set()
        uniq = []
        for r in interface[cid]:
            if r["resid"] not in seen:
                seen.add(r["resid"])
                uniq.append(r)
        interface[cid] = sorted(uniq, key=lambda r: r["resid"])

    # 2) SASA:复合物 + 隔离链 -> BSA
    sasa_complex = _compute_sasa_map(model, include_het)
    with tempfile.TemporaryDirectory() as tmp:
        sasa_iso_a = _extract_chain_sasa(struct=model.get_parent() if hasattr(model, "get_parent") else model,
                                         chain_id=ca, tmpdir=tmp, include_het=include_het)
        sasa_iso_b = _extract_chain_sasa(model.get_parent() if hasattr(model, "get_parent") else model,
                                         chain_id=cb, tmpdir=tmp, include_het=include_het)

    resname_map = {}
    for chain in model:
        for res in chain:
            if not include_het and res.id[0] != " ":
                continue
            resname_map[(chain.id, res.id[1])] = res.get_resname().strip()

    bsa_rows = []
    for (cid, rid), sasa_iso in list(sasa_iso_a.items()) + list(sasa_iso_b.items()):
        sasa_c = sasa_complex.get((cid, rid), 0.0)
        delta = round(sasa_iso - sasa_c, 2)
        if delta >= bsa_min:
            bsa_rows.append({"chain": cid, "resid": rid, "resname": resname_map.get((cid, rid), "?"),
                             "bsa": delta,
                             "sasa_isolated": round(sasa_iso, 2), "sasa_complex": round(sasa_c, 2)})
    bsa_rows.sort(key=lambda r: -r["bsa"])

    report = {
        "complex": path,
        "chains": list(chains),
        "cutoff_angstrom": cutoff,
        "n_contacts": len(contacts),
        "contacts": contacts[:200],               # 防超大输出,完整列表可调
        "n_contacts_truncated": len(contacts) > 200,
        "interface_residues": {k: v for k, v in interface.items()},
        "bsa_total": round(sum(r["bsa"] for r in bsa_rows), 2),
        "bsa_min_threshold": bsa_min,
        "bsa_per_residue": bsa_rows,
    }
    return report


def _json_safe(obj):
    """把非有限浮点(NaN/Infinity,如空界面的 log10(0))清洗为 None,
    保证 JSON 报告可被严格解析器(含 DSH 工具结果层)无损读取。"""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    return obj


def _print_report(r):
    print("complex: %s  chains: %s" % (r["complex"], ",".join(r["chains"])))
    print("cutoff: %.1f A | contacts: %d%s" % (
        r["cutoff_angstrom"], r["n_contacts"],
        " (truncated)" if r.get("n_contacts_truncated") else ""))
    print("interface residues:")
    for cid, res_list in r["interface_residues"].items():
        print("  chain %s: %d residues -> %s" % (
            cid, len(res_list),
            ", ".join("%s%d" % (x["resname"], x["resid"]) for x in res_list[:40]) +
            (" ..." if len(res_list) > 40 else "")))
    print("BSA total: %.1f A^2 (top 10):" % r["bsa_total"])
    for row in r["bsa_per_residue"][:10]:
        print("  %s %s%d  dSASA=%6.2f" % (row["chain"], row.get("resname", "?"), row["resid"], row["bsa"]))


def main(argv=None):
    ap = argparse.ArgumentParser(description="protein-protein interface analysis (Bio.PDB)")
    ap.add_argument("--complex", required=True, help="complex PDB file")
    ap.add_argument("--chains", nargs="+", default=["A", "B"], help="chain ids (default: A B)")
    ap.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF, help="contact cutoff in Angstrom")
    ap.add_argument("--include-het", action="store_true", help="include HETATM residues")
    ap.add_argument("--bsa-min", type=float, default=DEFAULT_BSA_MIN, help="min BSA per residue")
    ap.add_argument("--out", help="write JSON report to file")
    args = ap.parse_args(argv)

    try:
        report = analyze_complex(args.complex, chains=tuple(args.chains),
                                 cutoff=args.cutoff, include_het=args.include_het,
                                 bsa_min=args.bsa_min)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    _print_report(report)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(_json_safe(report), f, ensure_ascii=False, indent=2,
                      allow_nan=False)
        print("JSON written: %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
