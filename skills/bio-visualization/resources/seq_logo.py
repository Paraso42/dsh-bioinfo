#!/usr/bin/env python3
r"""seq_logo.py — 序列 Logo 生成(Logomaker)

输入:比对好的 FASTA(等长)或多行文本;输出概率型或信息量型 Logo PNG。

  & 'D:\bioai\venv\Scripts\python.exe' seq_logo.py aligned.fasta --out logo.png
  & 'D:\bioai\venv\Scripts\python.exe' seq_logo.py aligned.fasta --info --out logo_bits.png \
      --title "barnase homologs logo"
  & 'D:\bioai\venv\Scripts\python.exe' seq_logo.py aligned.fasta --first 40 --stack-width 0.9

选项:
  --info  信息量(bits,基于 BLOSUM62 背景);默认概率频率
  --first 只画前 N 列(长 MSA 建议 40-80)
  --start 起始列(1-based)
"""
import argparse
import sys


def _setup():
    """Agg 后端 + 中文字体自动选择(雅黑/黑体/思源黑体/宋体 → DejaVu 兜底)。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    installed = {f.name for f in font_manager.fontManager.ttflist}
    for cand in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC",
                 "WenQuanYi Micro Hei", "SimSun"):
        if cand in installed:
            plt.rcParams["font.sans-serif"] = [cand, "DejaVu Sans"]
            break
    else:
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        print("WARN: no CJK font found; Chinese labels will render as boxes",
              file=sys.stderr)
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def load_alignment(path):
    seqs, sid = [], None
    buf = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if sid is not None:
                    seqs.append(("".join(buf)))
                sid = line[1:].split()[0]
                buf = []
            elif line.strip():
                buf.append(line.strip().upper().replace(" ", ""))
        if sid is not None:
            seqs.append("".join(buf))
    seqs = [s for s in seqs if s]
    if not seqs:
        raise ValueError("no sequences parsed from %s" % path)
    L = {len(s) for s in seqs}
    if len(L) > 1:
        raise ValueError("sequences are not aligned (lengths: %s)" % sorted(L))
    return seqs


def make_logo(seqs, out="logo.png", info=False, start=1, first=None, title=None,
              stack_width=0.9, dpi=150):
    plt = _setup()
    import logomaker
    import pandas as pd
    cols = list("ACDEFGHIKLMNPQRSTVWY")
    start = max(1, start)
    end = len(seqs[0]) if first is None else min(len(seqs[0]), start - 1 + first)
    sub = [s[start - 1:end] for s in seqs]
    import numpy as np
    mat = np.zeros((len(sub), len(cols)))
    for i, s in enumerate(sub):
        for j, c in enumerate(s):
            if c in cols:
                mat[i, cols.index(c)] += 1.0
    counts = pd.Series(mat.sum(axis=0), index=cols)
    counts = counts[counts > 0]
    probs = counts / counts.sum()
    if info:
        bg = pd.Series({
            "A": 0.078, "C": 0.019, "D": 0.053, "E": 0.063, "F": 0.039,
            "G": 0.072, "H": 0.023, "I": 0.053, "K": 0.059, "L": 0.091,
            "M": 0.023, "N": 0.043, "P": 0.052, "Q": 0.042, "R": 0.051,
            "S": 0.068, "T": 0.059, "V": 0.066, "W": 0.014, "Y": 0.032})
        import numpy as np
        h = probs * np.log2((probs / bg[probs.index]).clip(lower=1e-12))
        df = h.to_frame().T
        ylabel = "information (bits)"
    else:
        df = probs.to_frame().T
        ylabel = "frequency"
    fig, ax = plt.subplots(figsize=(max(4, 0.35 * len(df.columns)), 3))
    ww = logomaker.Logo(df, ax=ax, width=stack_width)
    ww.style_spines(visible=False)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("position (alignment %d..%d, n=%d)" % (start, end, len(seqs)))
    if title:
        ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return {"out": out, "n_sequences": len(seqs), "positions": [start, end],
            "conserved": [str(k) for k, v in probs.items() if v >= 0.9]}


def main(argv=None):
    ap = argparse.ArgumentParser(description="sequence logo via Logomaker")
    ap.add_argument("fasta")
    ap.add_argument("--out", default="logo.png")
    ap.add_argument("--info", action="store_true")
    ap.add_argument("--first", type=int)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--title")
    ap.add_argument("--stack-width", type=float, default=0.9)
    args = ap.parse_args(argv)
    try:
        seqs = load_alignment(args.fasta)
        r = make_logo(seqs, out=args.out, info=args.info, start=args.start,
                      first=args.first, title=args.title, stack_width=args.stack_width)
        print("n=%d positions %d-%d -> %s" % (r["n_sequences"], r["positions"][0],
                                              r["positions"][1], r["out"]))
        if r["conserved"]:
            print("fully conserved columns: %s" % ",".join(r["conserved"]))
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
