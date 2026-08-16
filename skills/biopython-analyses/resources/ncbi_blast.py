# -*- coding: utf-8 -*-
"""
ncbi_blast.py — 带浏览器 UA 与重试的 NCBI BLAST URLAPI 客户端
（biopython-analyses 技能的随包资源；本机 2026-07 实测通过）

背景（重要）：
  NCBI BLAST 服务端会按 User-Agent 间歇性拒连：Bio.Blast.NCBIWWW.qblast
  （UA="BiopythonClient"）与默认 Python urllib 的 POST 常被直接断开
  （RemoteDisconnected: Remote end closed connection without response），
  对 nr 大库 + entrez 过滤的重请求几乎必现；换浏览器 UA 即稳定成功。
  DNS/TCP/GET 均正常——是服务端行为，不是防火墙问题。Entrez（eutils）不受影响。

用法（命令行）：
  python ncbi_blast.py <program> <database> <query.fasta>
      [--entrez-query "Viridiplantae[ORGN]"] [--expect 1e-10] [--hits 10]
      [-o result.xml] [--max-wait 480] [--email you@example.com]
  program : blastp | blastn | blastx | tblastn | tblastx
  database: nr | nt | refseq_protein | swissprot | ...
  输出    : BLAST XML（-o 存盘；否则打印到 stdout）

编程调用：
  import sys; sys.path.insert(0, r"C:\\deepseek-harness\\.dsh\\.agent-presets\\bioinfo\\skills\\biopython-analyses\\resources")
  from ncbi_blast import ncbi_blast
  xml = ncbi_blast("blastp", "nr", ">q\\nMSPQTETK...", entrez_query="Viridiplantae[ORGN]",
                   expect=1e-20, hitlist=10)
  from Bio import SearchIO
  from io import StringIO
  q = SearchIO.read(StringIO(xml), "blast-xml")

附注：
  - NCBI 官方限频：两次 Put 间隔 ≥10 s；同一 RID 轮询间隔 ≥60 s（本脚本已内置）。
  - NCBI FTP（ftp.ncbi.nlm.nih.gov）大文件被带宽整形（~1 KB/s 级），
    改用 Datasets API（https://api.ncbi.nlm.nih.gov/datasets/v2alpha/...）
    或 Entrez efetch 的 seq_start/seq_stop 分段下载。
"""
import argparse
import re
import sys
import time
from urllib.request import Request, urlopen
from urllib.parse import urlencode

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36")


def ncbi_blast(program, database, query, entrez_query="(none)", expect=1e-10,
               hitlist=10, max_wait=480, email="your@email.com",
               tool="dsh-ncbi-blast-helper", url=URL, quiet=False):
    """提交一次 BLAST 并轮询取回 XML 字符串；失败返回 None。

    query 为 FASTA 文本（含 > 头行）。
    协议：CMD=Put（POST，浏览器 UA，最多 6 次重试）→ CMD=Get 轮询（20 s 起、之后 60 s）
          → 返回最终 XML。
    """
    put = {
        "CMD": "Put", "PROGRAM": program, "DATABASE": database, "QUERY": query,
        "EXPECT": str(expect), "HITLIST_SIZE": str(hitlist),
        "ALIGNMENTS": str(hitlist), "DESCRIPTIONS": str(hitlist),
        "FORMAT_TYPE": "XML", "MATRIX_NAME": "BLOSUM62",
        "GAPCOSTS": "11 1", "FILTER": "F", "COMPOSITION_BASED_STATISTICS": "2",
        "ENTREZ_QUERY": entrez_query, "email": email, "tool": tool,
    }
    rid = None
    for attempt in range(6):
        try:
            req = Request(url, urlencode(put).encode(), {"User-Agent": UA})
            with urlopen(req, timeout=90) as r:
                page = r.read().decode("utf-8", "replace")
            m = re.search(r"RID = (\S+)", page)
            if m:
                rid = m.group(1)
                if not quiet:
                    print(f"[ncbi_blast] RID = {rid}", file=sys.stderr, flush=True)
                break
            if not quiet:
                print(f"[ncbi_blast] put attempt {attempt + 1}: no RID; "
                      f"first 200 chars: {page[:200]!r}", file=sys.stderr, flush=True)
        except Exception as e:
            if not quiet:
                print(f"[ncbi_blast] put attempt {attempt + 1}: "
                      f"{type(e).__name__} {repr(e)[:100]}", file=sys.stderr, flush=True)
        time.sleep(10 + attempt * 10)
    if not rid:
        return None

    get = {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid,
           "ALIGNMENTS": str(hitlist), "DESCRIPTIONS": str(hitlist)}
    t0 = time.time()
    delay = 20
    while time.time() - t0 < max_wait:
        time.sleep(delay)
        delay = 60
        try:
            req = Request(url + "?" + urlencode(get), headers={"User-Agent": UA})
            with urlopen(req, timeout=90) as r:
                body = r.read().decode("utf-8", "replace")
        except Exception as e:
            if not quiet:
                print(f"[ncbi_blast] get error: {type(e).__name__} "
                      f"{repr(e)[:100]}", file=sys.stderr, flush=True)
            continue
        if "Status=" not in body:
            return body  # 结果已就绪（XML 无 Status 标签）
        m = re.search(r"Status=(\S+)", body)
        st = m.group(1) if m else "?"
        if not quiet:
            print(f"[ncbi_blast] status={st} elapsed={time.time() - t0:.0f}s",
                  file=sys.stderr, flush=True)
        if st.upper() == "READY":
            time.sleep(5)
            req = Request(url + "?" + urlencode(get), headers={"User-Agent": UA})
            with urlopen(req, timeout=90) as r:
                return r.read().decode("utf-8", "replace")
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description="NCBI BLAST URLAPI 客户端（浏览器 UA + 重试）")
    ap.add_argument("program", help="blastp/blastn/blastx/tblastn/tblastx")
    ap.add_argument("database", help="nr/nt/refseq_protein/swissprot/...")
    ap.add_argument("query_fasta", help="query FASTA 文件")
    ap.add_argument("--entrez-query", default="(none)", help='如 "Viridiplantae[ORGN]"')
    ap.add_argument("--expect", default=1e-10, type=float)
    ap.add_argument("--hits", default=10, type=int,
                    help="hitlist/alignments/descriptions 数量")
    ap.add_argument("-o", "--output", help="输出 XML 文件（默认打印到 stdout）")
    ap.add_argument("--max-wait", default=480, type=int, help="轮询最长等待秒数")
    ap.add_argument("--email", default="your@email.com", help="NCBI 联系邮箱（建议填写真实地址）")
    a = ap.parse_args(argv)
    with open(a.query_fasta, encoding="utf-8") as f:
        query = f.read()
    xml = ncbi_blast(a.program, a.database, query, entrez_query=a.entrez_query,
                     expect=a.expect, hitlist=a.hits, max_wait=a.max_wait,
                     email=a.email)
    if xml is None:
        print("BLAST 失败：重试后仍无法取得结果", file=sys.stderr)
        return 1
    if a.output:
        with open(a.output, "w", encoding="utf-8") as f:
            f.write(xml)
        print(f"已保存 {a.output}（{len(xml)} bytes）")
    else:
        print(xml)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
