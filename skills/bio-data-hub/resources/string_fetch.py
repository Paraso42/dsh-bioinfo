#!/usr/bin/env python3
r"""string_fetch.py — STRING-DB REST 客户端(互作网络 / ID 映射,纯标准库)

子命令:
  network  互作伙伴:  string_fetch.py network P00698 --species 9606 --score 400 --limit 50
  map      ID 映射:   string_fetch.py map P00698 ENSG00000130203 --species 9606
  image    网络图 PNG:string_fetch.py image P00698 --species 9606 --out net.png

基址 string-db.org/api/tsv(无 key;大任务限速,脚本已做节流)。
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://string-db.org/api"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                     "(KHTML, like Gecko) Chrome/124.0 Safari/537.36 DSH-bioinfo/1.0"}


def _get(url, retries=5, timeout=60):
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


def _parse_tsv(text):
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def network(proteins, species=9606, score=400, limit=50):
    q = urllib.parse.quote("%0d".join(proteins))
    url = ("%s/tsv/network?identifiers=%s&species=%s&required_score=%d&limit=%d"
           % (BASE, q, species, score, limit))
    rows = _parse_tsv(_get(url))
    out = []
    for r in rows:
        out.append({
            "protein_a": r.get("preferredName_A") or r.get("stringId_A"),
            "protein_b": r.get("preferredName_B") or r.get("stringId_B"),
            "score": float(r.get("score", 0)),
            "experiments": r.get("experiments"), "database": r.get("database"),
            "textmining": r.get("textmining"),
        })
    return out


def map_ids(identifiers, species=9606):
    # get_string_ids 仅接受单个 ID(多 ID 新行编码会 404)→ 逐个请求
    rows = []
    for ident in identifiers:
        q = urllib.parse.quote(ident)
        url = "%s/tsv/get_string_ids?identifiers=%s&species=%s" % (BASE, q, species)
        try:
            parsed = _parse_tsv(_get(url))
            for r in parsed:
                r["queryItem"] = ident        # v12 响应无 queryItem 列,补上输入 ID
            rows.extend(parsed)
        except Exception as e:
            rows.append({"queryItem": ident, "error": str(e)})
        time.sleep(1.0)   # STRING 限速礼貌
    return [{k: r.get(k) for k in ("queryItem", "stringId", "preferredName", "annotation")}
            for r in rows]


def network_image(proteins, species=9606, score=400, out=None):
    q = urllib.parse.quote("%0d".join(proteins))
    url = ("%s/image/network?identifiers=%s&species=%s&required_score=%d"
           % (BASE, q, species, score))
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    if out:
        with open(out, "wb") as f:
            f.write(data)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="STRING-DB REST client")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_net = sub.add_parser("network")
    p_net.add_argument("proteins", nargs="+")
    p_net.add_argument("--species", default="9606")
    p_net.add_argument("--score", type=int, default=400)
    p_net.add_argument("--limit", type=int, default=50)
    p_net.add_argument("--out")
    p_map = sub.add_parser("map")
    p_map.add_argument("identifiers", nargs="+")
    p_map.add_argument("--species", default="9606")
    p_map.add_argument("--out")
    p_img = sub.add_parser("image")
    p_img.add_argument("proteins", nargs="+")
    p_img.add_argument("--species", default="9606")
    p_img.add_argument("--score", type=int, default=400)
    p_img.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "network":
            rows = network(args.proteins, species=args.species, score=args.score, limit=args.limit)
            print("partners: %d" % len(rows))
            for r in rows[:20]:
                print("  %-12s - %-12s  score=%.2f (exp=%s db=%s)" % (
                    r["protein_a"], r["protein_b"], r["score"], r["experiments"], r["database"]))
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    json.dump(rows, f, ensure_ascii=False, indent=2)
                print("JSON written: %s" % args.out)
        elif args.cmd == "map":
            rows = map_ids(args.identifiers, species=args.species)
            for r in rows:
                print("  %-20s -> %-12s %s" % (r["queryItem"], r.get("stringId"), r.get("annotation", "")))
            if args.out:
                with open(args.out, "w", encoding="utf-8") as f:
                    json.dump(rows, f, ensure_ascii=False, indent=2)
        else:
            p = network_image(args.proteins, species=args.species, score=args.score, out=args.out)
            print("PNG written: %s" % p)
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
