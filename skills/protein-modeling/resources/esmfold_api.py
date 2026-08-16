#!/usr/bin/env python3
r"""esmfold_api.py — ESM Metagenomic Atlas API 客户端(序列 → PDB,免费无 key)

纯标准库(urllib),浏览器 UA + 退避重试(纪律同 ncbi_blast.py)。
适合:快速粗模、网络兜底时的结构预测通道;重活走 LocalColabFold(AF2)。

用法:
  & 'C:\Program Files\Python313\python.exe' esmfold_api.py --fasta query.fasta --out D:\bioai\jobs\esm_out.pdb
  & 'C:\Program Files\Python313\python.exe' esmfold_api.py "MQIFVKTLTGKTITLEVEPS..." --out out.pdb

编程调用:
  from esmfold_api import fold_sequence
  pdb_text = fold_sequence("MQIFVKT...", retries=4)
"""
import argparse
import re
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

DEFAULT_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36")
SEQ_RE = re.compile(r"^[A-Za-z*\-]+$")


def fold_sequence(sequence, url=DEFAULT_URL, timeout=240, retries=4, quiet=False):
    """提交一条蛋白序列,返回 PDB 文本。失败自动退避重试,最终抛 RuntimeError。"""
    seq = "".join(sequence.split())            # 去空白/换行
    if not seq or not SEQ_RE.match(seq):
        raise ValueError("invalid protein sequence (letters/*/- only)")
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, data=seq.encode("ascii"),
                          headers={"User-Agent": UA, "Content-Type": "text/plain"})
            with urlopen(req, timeout=timeout) as resp:
                body = resp.read().decode("utf-8", "replace")
            if "ATOM" in body or "HETATM" in body:
                return body
            raise RuntimeError("unexpected response (%d bytes): %s" % (len(body), body[:200]))
        except (URLError, HTTPError, TimeoutError, OSError, RuntimeError) as e:
            last = e
            if not quiet:
                print("[attempt %d/%d] %s" % (attempt, retries, e), file=sys.stderr)
            if attempt < retries:
                time.sleep(min(20 * attempt, 120))
    raise RuntimeError(
        "ESM Atlas API failed after %d attempts: %s. "
        "NOTE: the Atlas has been intermittently unavailable (repeated HTTP 504). "
        "Local fallback (offline, GPU): run_colabfold.ps1 -MsaMode single_sequence "
        "(see protein-modeling skill / af2_predict tool)" % (retries, last))


def read_fasta(path):
    """解析 fasta,返回 [(header, seq), ...]。"""
    out, cur_id, cur = [], None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    out.append((cur_id, "".join(cur)))
                cur_id, cur = line[1:].split()[0] if line[1:].strip() else "seq", []
            else:
                cur.append(line)
    if cur_id is not None:
        out.append((cur_id, "".join(cur)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="ESM Atlas API protein folding client")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("sequence", nargs="?", help="protein sequence (no spaces)")
    src.add_argument("--fasta", help="fasta file (all sequences folded)")
    ap.add_argument("--out", help="output PDB path; with --fasta use {id} placeholder or a directory")
    ap.add_argument("--retries", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=240)
    ap.add_argument("--url", default=DEFAULT_URL)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    if args.fasta:
        records = read_fasta(args.fasta)
    else:
        records = [("seq", args.sequence)]
    if not records:
        print("ERROR: no sequences", file=sys.stderr)
        return 2

    ok = 0
    for idx, (seq_id, seq) in enumerate(records):
        try:
            pdb_text = fold_sequence(seq, url=args.url, timeout=args.timeout,
                                     retries=args.retries, quiet=args.quiet)
        except Exception as e:
            print("ERROR [%s]: %s" % (seq_id, e), file=sys.stderr)
            continue
        if args.out:
            if len(records) == 1:
                dest = args.out
            elif "{id}" in args.out:
                dest = args.out.replace("{id}", seq_id)
            else:
                import os
                os.makedirs(args.out, exist_ok=True)
                dest = os.path.join(args.out, seq_id + ".pdb")
            with open(dest, "w", encoding="utf-8") as f:
                f.write(pdb_text)
            print("wrote %s (%d bytes)" % (dest, len(pdb_text)))
        else:
            sys.stdout.write(pdb_text)
        ok += 1
    return 0 if ok == len(records) else 1


if __name__ == "__main__":
    sys.exit(main())
