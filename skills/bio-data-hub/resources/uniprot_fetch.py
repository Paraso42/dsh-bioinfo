#!/usr/bin/env python3
r"""uniprot_fetch.py — UniProt REST 客户端(纯标准库,浏览器 UA + 重试)

子命令:
  get     按 accession 取注释/序列:uniprot_fetch.py get P00698 --json out.json
  fasta   取 FASTA:uniprot_fetch.py fasta P00698 P00761 --out seqs.fasta
  search  检索:uniprot_fetch.py search "organism:9606 AND reviewed:true AND keyword:Kinase" \
             --size 10 --fields accession,id,gene_names,protein_name --format tsv

基址 rest.uniprot.org;自动处理 429(退避)/5xx(重试);结果可写文件。
"""
import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://rest.uniprot.org"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                     "(KHTML, like Gecko) Chrome/124.0 Safari/537.36 DSH-bioinfo/1.0"}


def _get(url, retries=5, timeout=60, accept=None):
    headers = dict(UA)
    if accept:
        headers["Accept"] = accept
    for attempt in range(retries):
        req = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("HTTP %s for %s: %s" % (e.code, url, e.read().decode("utf-8", "replace")[:300]))
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("network error for %s: %s" % (url, e))
    raise RuntimeError("giving up on %s" % url)


def fetch_entry(accession):
    url = BASE + "/uniprotkb/%s.json" % urllib.parse.quote(accession)
    return json.loads(_get(url, accept="application/json"))


def fetch_fasta(accessions):
    url = BASE + "/uniprotkb/stream?query=" + urllib.parse.quote(
        " OR ".join("accession:%s" % a for a in accessions)) + "&format=fasta"
    return _get(url)


def search(query, size=10, fields=None, fmt="tsv"):
    # 常见别名纠错:organism: → organism_id:
    query = re.sub(r"\borganism:", "organism_id:", query)
    url = BASE + "/uniprotkb/search?query=" + urllib.parse.quote(query) + "&size=%d" % size
    if fields:
        url += "&fields=" + urllib.parse.quote(",".join(fields))
    url += "&format=" + fmt
    return _get(url)


def _entry_summary(e):
    def _g(obj, *keys):
        cur = obj
        for k in keys:
            if isinstance(cur, dict) and k in cur:
                cur = cur[k]
            else:
                return None
        return cur
    return {
        "accession": e.get("primaryAccession"),
        "id": e.get("uniProtkbId"),
        "name": _g(e, "proteinDescription", "recommendedName", "fullName", "value"),
        "gene": (_g(e, "genes", 0, "geneName", "value") if e.get("genes") else None),
        "organism": _g(e, "organism", "scientificName"),
        "length": _g(e, "sequence", "length"),
        "reviewed": e.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
        "function": _g(e, "comments", 0, "texts", 0, "value") if e.get("comments") else None,
        "subcellular": [_g(c, "subcellularLocation", "location", "value")
                        for c in e.get("comments", []) if c.get("commentType") == "SUBCELLULAR LOCATION"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description="UniProt REST client")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_get = sub.add_parser("get")
    p_get.add_argument("accession")
    p_get.add_argument("--out")
    p_fa = sub.add_parser("fasta")
    p_fa.add_argument("accessions", nargs="+")
    p_fa.add_argument("--out", required=True)
    p_s = sub.add_parser("search")
    p_s.add_argument("query")
    p_s.add_argument("--size", type=int, default=10)
    p_s.add_argument("--fields", nargs="+")
    p_s.add_argument("--format", default="tsv", choices=["tsv", "json", "list", "fasta"])
    p_s.add_argument("--out")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "get":
            entry = fetch_entry(args.accession)
            summary = _entry_summary(entry)
            print(json.dumps(summary, ensure_ascii=False, indent=2))
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    json.dump(entry, f, ensure_ascii=False, indent=2)
                print("full JSON written: %s" % args.out)
        elif args.cmd == "fasta":
            text = fetch_fasta(args.accessions)
            with open(args.out, "w", encoding="utf-8") as f:
                f.write(text)
            n = sum(1 for l in text.splitlines() if l.startswith(">"))
            print("wrote %s (%d records)" % (args.out, n))
        else:
            text = search(args.query, size=args.size, fields=args.fields, fmt=args.format)
            print(text)
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    f.write(text)
                print("written: %s" % args.out)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
