#!/usr/bin/env python3
r"""kegg_fetch.py — KEGG REST 客户端(find / get / link,纯标准库)

子命令:
  find     kegg_fetch.py find genes kinase --org hsa --limit 20
           kegg_fetch.py find pathway apoptosis --org hsa
  get      kegg_fetch.py get hsa:207 hsa:5594 --out kegg.txt
  link     kegg_fetch.py link pathway hsa:207          (基因 → 通路)
           kegg_fetch.py link genes path:hsa04210      (通路 → 基因)

基址 rest.kegg.jp;单请求限 10 个 ID(get 自动分批)。
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://rest.kegg.jp"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                     "(KHTML, like Gecko) Chrome/124.0 Safari/537.36 DSH-bioinfo/1.0"}


def _get(path, retries=5, timeout=90):
    for attempt in range(retries):
        req = urllib.request.Request(BASE + path, headers=UA)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("HTTP %s for %s" % (e.code, path))
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < retries - 1:
                time.sleep(min(2 ** attempt, 30))
                continue
            raise RuntimeError("network error for %s: %s" % (path, e))
    raise RuntimeError("giving up on %s" % path)


def find(db, query, org=None, limit=None):
    q = urllib.parse.quote(query)
    if org:
        q += "+" + org                      # KEGG find 的 org 用 '+' 拼接
    path = "/find/%s/%s" % (db, q)
    lines = _get(path).splitlines()
    if limit:
        lines = lines[:limit]
    out = []
    for l in lines:
        if "\t" in l:
            k, v = l.split("\t", 1)
            out.append({"id": k, "description": v})
    return out


def get(ids, chunk=10):
    chunks = [ids[i:i + chunk] for i in range(0, len(ids), chunk)]
    out = []
    for c in chunks:
        time.sleep(0.6)  # KEGG 限速礼貌
        text = _get("/get/" + "+".join(c))
        out.append(text)
    return "".join(out)


def link(target_db, source_ids, chunk=10):
    chunks = [source_ids[i:i + chunk] for i in range(0, len(source_ids), chunk)]
    rows = []
    for c in chunks:
        time.sleep(0.6)
        text = _get("/link/%s/%s" % (target_db, "+".join(c)))
        for l in text.splitlines():
            if "\t" in l:
                k, v = l.split("\t", 1)
                rows.append({"source": k, "target": v})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description="KEGG REST client")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_find = sub.add_parser("find")
    p_find.add_argument("db", help="e.g. genes / pathway / compound / enzyme")
    p_find.add_argument("query")
    p_find.add_argument("--org", help="e.g. hsa")
    p_find.add_argument("--limit", type=int)
    p_find.add_argument("--out")
    p_get = sub.add_parser("get")
    p_get.add_argument("ids", nargs="+")
    p_get.add_argument("--out", required=True)
    p_link = sub.add_parser("link")
    p_link.add_argument("target_db", help="e.g. pathway / genes / enzyme")
    p_link.add_argument("ids", nargs="+")
    p_link.add_argument("--out")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "find":
            rows = find(args.db, args.query, org=args.org, limit=args.limit)
            for r in rows:
                print("%-12s %s" % (r["id"], r["description"][:80]))
            print("(%d hits)" % len(rows))
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    json.dump(rows, f, ensure_ascii=False, indent=2)
        elif args.cmd == "get":
            text = get(args.ids)
            with open(args.out, "w", encoding="utf-8") as f:
                f.write(text)
            print("wrote %s (%d entries)" % (args.out, text.count("ENTRY")))
        else:
            rows = link(args.target_db, args.ids)
            for r in rows:
                print("%-12s -> %s" % (r["source"], r["target"]))
            print("(%d links)" % len(rows))
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    json.dump(rows, f, ensure_ascii=False, indent=2)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
