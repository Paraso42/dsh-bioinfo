#!/usr/bin/env python3
r"""pdb_fetch.py — RCSB PDB REST 客户端(元数据 / 结构下载 / 序列 / 检索)

子命令:
  meta      pdb_fetch.py meta 1yph --json meta.json
            (含逐链聚合物残基数 chain_residues_polymer,来自结构文件 ATOM 计数)
  download  pdb_fetch.py download 1yph --format pdb --out 1yph.pdb   (或 --format cif)
  fasta     pdb_fetch.py fasta 1yph --out 1yph.fasta
  search    pdb_fetch.py search '{"query":{"type":"terminal","service":"text",
              "parameters":{"value":"trypsin","attribute":"struct.title.pdbx"}},"return_type":"entry"}' --out hits.json

基址 data.rcsb.org / files.rcsb.org;内置重试。
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

DATA = "https://data.rcsb.org/rest/v1"
FILES = "https://files.rcsb.org"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                     "(KHTML, like Gecko) Chrome/124.0 Safari/537.36 DSH-bioinfo/1.0"}


def _get(url, retries=5, timeout=90):
    for attempt in range(retries):
        req = urllib.request.Request(url, headers=UA)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("HTTP %s for %s" % (e.code, url))
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("network error for %s: %s" % (url, e))
    raise RuntimeError("giving up on %s" % url)


def _post_json(url, payload, retries=3, timeout=120):
    for attempt in range(retries):
        req = urllib.request.Request(url, data=json.dumps(payload).encode("utf-8"),
                                     headers={**UA, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("HTTP %s for search" % e.code)
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("network error in search: %s" % e)
    raise RuntimeError("search giving up")


def _parse_chain_stats_from_pdb(text):
    """从 PDB 文本提取:链 ID、非聚合配体、逐链聚合物残基数。

    残基计数只统计 ATOM 记录(标准聚合物残基;结晶水/配体为 HETATM 不计数),
    按 (chain, resseq, icode) 去重,只取第一个 MODEL(ENDMDL 即停)。
    """
    chains, ligs, resmap = set(), set(), {}
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM  "):
            c = line[21:22].strip()
            chains.add(c)
            resmap.setdefault(c, set()).add((line[22:27].strip(), line[26:27]))
        elif line.startswith("HETATM"):
            res = line[17:20].strip()
            if res not in ("HOH", "WAT"):
                ligs.add(res)
            chains.add(line[21:22].strip())
    return (sorted(c for c in chains if c), sorted(ligs),
            {c: len(keys) for c, keys in sorted(resmap.items())})


def fetch_metadata(pid):
    meta = json.loads(_get(DATA + "/core/entry/" + urllib.parse.quote(pid)))
    summary = _meta_summary(meta)
    try:
        text = download_structure(pid, fmt="pdb")
        chains, ligs, rescounts = _parse_chain_stats_from_pdb(text)
        if not summary["chains"]:
            summary["chains"], summary["ligands"] = chains, ligs
        summary["chain_residues_polymer"] = rescounts
        summary["source"] = "core/entry + structure file"
    except Exception:
        summary["source"] = "core/entry" if summary["chains"] else "core/entry (chains unknown)"
    return summary


def download_structure(pid, fmt="pdb"):
    url = FILES + "/download/%s.%s" % (urllib.parse.quote(pid), fmt)
    return _get(url, timeout=300)


def fetch_fasta(pid):
    return _get(FILES + "/fasta/%s" % urllib.parse.quote(pid))


def search_entries(payload):
    return _post_json(SEARCH, payload)


def _meta_summary(m):
    """core/entry 返回的字段在顶层('entry' 只是 ID 字符串)。"""
    info = m.get("rcsb_entry_info", {})
    res = info.get("resolution_combined")
    if res is None:
        res = (m.get("diffrn_resolution_high") or {}).get("value")
    method = None
    exptl = m.get("exptl") or []
    if exptl:
        method = exptl[0].get("method")
    return {
        "pdb_id": m.get("rcsb_id"),
        "title": (m.get("struct") or {}).get("title"),
        "resolution": res,
        "method": method,
        "deposit_date": info.get("deposit_date"),
        "polymer_count": info.get("polymer_entity_count_deposited"),
        "chains": _chain_ids(m),
        "ligands": _ligands(m),
    }


def _chain_ids(m):
    ids = []
    for inst in (m.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_instances") or []):
        ids.append(inst.get("asym_id"))
    return ids


def _ligands(m):
    out = []
    for grp in (m.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_instances") or []):
        out.append(grp.get("chem_comp_id"))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="RCSB PDB REST client")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_meta = sub.add_parser("meta")
    p_meta.add_argument("pdb_id")
    p_meta.add_argument("--json")
    p_dl = sub.add_parser("download")
    p_dl.add_argument("pdb_id")
    p_dl.add_argument("--format", default="pdb", choices=["pdb", "cif"])
    p_dl.add_argument("--out", required=True)
    p_fa = sub.add_parser("fasta")
    p_fa.add_argument("pdb_id")
    p_fa.add_argument("--out", required=True)
    p_s = sub.add_parser("search")
    p_s.add_argument("query_json", help="RCSB search API query object (JSON string or @file)")
    p_s.add_argument("--out")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "meta":
            summary = fetch_metadata(args.pdb_id)
            print(json.dumps(summary, ensure_ascii=False, indent=2))
            if args.json:
                with open(args.json, "w", encoding="utf-8") as f:
                    json.dump(summary, f, ensure_ascii=False, indent=2)
                print("JSON written: %s" % args.json)
        elif args.cmd == "download":
            text = download_structure(args.pdb_id, fmt=args.format)
            with open(args.out, "w", encoding="utf-8") as f:
                f.write(text)
            n = sum(1 for l in text.splitlines() if l.startswith(("ATOM", "HETATM")))
            print("wrote %s (%d coordinate records)" % (args.out, n))
        elif args.cmd == "fasta":
            text = fetch_fasta(args.pdb_id)
            with open(args.out, "w", encoding="utf-8") as f:
                f.write(text)
            print("wrote %s" % args.out)
        else:
            q = args.query_json
            if q.startswith("@"):
                with open(q[1:], encoding="utf-8") as f:
                    q = f.read()
            result = search_entries(json.loads(q))
            ids = [r.get("identifier") for r in result.get("result_set", [])]
            print("hits (%d): %s" % (len(ids), ", ".join(ids[:30])))
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    json.dump(result, f, ensure_ascii=False, indent=2)
                print("written: %s" % args.out)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
