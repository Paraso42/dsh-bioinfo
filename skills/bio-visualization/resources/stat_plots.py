#!/usr/bin/env python3
r"""stat_plots.py — 生信统计可视化(火山图 / MA 图 / 热图 / 环形基因组图)

运行环境(D:\bioai\venv:matplotlib/seaborn/pycirclize/pandas):
  & 'D:\bioai\venv\Scripts\python.exe' stat_plots.py volcano  deg.csv --log2fc log2FoldChange \
      --pvalue padj --genes gene --thresholds 1 0.05 --out volcano.png
  & 'D:\bioai\venv\Scripts\python.exe' stat_plots.py ma  deg.csv --log2fc log2FoldChange \
      --base-mean baseMean --genes gene --out ma.png
  & 'D:\bioai\venv\Scripts\python.exe' stat_plots.py heatmap matrix.csv --zscore --cluster \
      --cmap RdBu_r --out heatmap.png
  & 'D:\bioai\venv\Scripts\python.exe' stat_plots.py circos sectors.csv --out circos.png

CSV 约定:
  volcano/ma:行 = 基因;列名可用 --log2fc/--pvalue/--genes 指定(默认自动识别常见列名)
  heatmap:首列行名,其余列为数值(样本 × 特征或特征 × 样本均可)
  circos:列 chrom,start,end[,name,value](BED 风格;value 用于颜色映射)
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


def _load_csv(path):
    import pandas as pd
    return pd.read_csv(path)


def _find_col(df, candidates, default=None):
    for c in candidates:
        if c in df.columns:
            return c
    return default


def volcano(csv_path, log2fc_col=None, pvalue_col=None, gene_col=None,
            lfc_thr=1.0, p_thr=0.05, top_n=15, out="volcano.png", dpi=150):
    import numpy as np
    import pandas as pd
    plt = _setup()
    df = _load_csv(csv_path)
    lfc = log2fc_col or _find_col(df, ["log2FoldChange", "log2FC", "logFC", "LFC", "log2(fc)"])
    pv = pvalue_col or _find_col(df, ["padj", "p_adjusted", "adj.P.Val", "FDR", "P.Value", "pvalue"])
    g = gene_col or _find_col(df, ["gene", "Gene", "symbol", "SYMBOL", "name"])
    if not lfc or not pv:
        raise ValueError("cannot auto-detect log2fc/pvalue columns; pass --log2fc/--pvalue")
    x = df[lfc].astype(float)
    y = -np.log10(df[pv].astype(float).clip(lower=1e-300))
    up = (x >= lfc_thr) & (y >= -np.log10(p_thr))
    down = (x <= -lfc_thr) & (y >= -np.log10(p_thr))
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x[~(up | down)], y[~(up | down)], s=8, c="#b0b0b0", alpha=0.6, label="ns")
    ax.scatter(x[up], y[up], s=8, c="#d62728", alpha=0.7, label="up")
    ax.scatter(x[down], y[down], s=8, c="#1f77b4", alpha=0.7, label="down")
    ax.axhline(-np.log10(p_thr), ls="--", c="k", lw=0.8)
    ax.axvline(lfc_thr, ls="--", c="k", lw=0.8)
    ax.axvline(-lfc_thr, ls="--", c="k", lw=0.8)
    if g:
        labels = df[g].astype(str).values
        sig = x[(up | down)].sort_values(ascending=False)
        for i in sig.index[:top_n]:
            ax.annotate(labels[i], (x[i], y[i]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(lfc)
    ax.set_ylabel("-log10(%s)" % pv)
    ax.set_title("volcano: up=%d down=%d" % (int(up.sum()), int(down.sum())))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return {"out": out, "up": int(up.sum()), "down": int(down.sum())}


def ma(csv_path, log2fc_col=None, base_mean_col=None, gene_col=None, out="ma.png", dpi=150):
    import numpy as np
    plt = _setup()
    df = _load_csv(csv_path)
    lfc = log2fc_col or _find_col(df, ["log2FoldChange", "log2FC", "logFC"])
    bm = base_mean_col or _find_col(df, ["baseMean", "AveExpr", "mean", "baseMean+1"])
    if not lfc or not bm:
        raise ValueError("cannot auto-detect log2fc/baseMean columns")
    x = np.log10(df[bm].astype(float).clip(lower=1e-6))
    y = df[lfc].astype(float)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y, s=6, c="#4c72b0", alpha=0.5)
    ax.axhline(0, c="k", lw=0.8)
    ax.set_xlabel("log10(%s)" % bm)
    ax.set_ylabel(lfc)
    ax.set_title("MA plot")
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return {"out": out}


def heatmap(csv_path, zscore=True, cluster=True, cmap="RdYlBu_r", figsize=None,
            out="heatmap.png", dpi=150):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    plt = _setup()
    df = _load_csv(csv_path)
    first = df.columns[0]
    idx = df[first].astype(str).tolist()
    mat = df.drop(columns=[first]).apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if mat.shape[0] == 0:
        raise ValueError("empty matrix")
    if zscore:
        mat = (mat - np.nanmean(mat, axis=1, keepdims=True)) / (np.nanstd(mat, axis=1, keepdims=True) + 1e-12)
    h, w = mat.shape
    fs = figsize or (max(6, w * 0.35), max(5, h * 0.28))
    cg = sns.clustermap(pd.DataFrame(mat, index=idx, columns=list(df.columns[1:])),
                        cmap=cmap, figsize=fs, row_cluster=cluster, col_cluster=cluster,
                        yticklabels=(h <= 60), xticklabels=True, cbar_kws={"label": "z-score" if zscore else "value"})
    cg.savefig(out, dpi=dpi)
    import matplotlib.pyplot as plt2
    plt2.close("all")
    return {"out": out, "shape": [h, w]}


def circos(csv_path, out="circos.png", dpi=150):
    from pycirclize import Circos
    import pandas as pd
    plt = _setup()
    df = _load_csv(csv_path)
    for col in ("chrom", "start", "end"):
        if col not in df.columns:
            raise ValueError("circos CSV needs columns chrom,start,end")
    chroms = df["chrom"].astype(str).unique()
    chrom_size = {c: int(df.loc[df["chrom"] == c, "end"].max()) for c in chroms}
    circos = Circos(sectors={c: s for c, s in chrom_size.items()}, space=3)
    has_value = "value" in df.columns
    vmin, vmax = (float(df["value"].min()), float(df["value"].max())) if has_value else (0, 1)
    for c in chroms:
        sector = circos.get_sector(c)
        track = sector.add_track((70, 90))
        sub = df[df["chrom"] == c]
        for _, row in sub.iterrows():
            v = float(row["value"]) if has_value else 1.0
            frac = (v - vmin) / max(1e-12, vmax - vmin) if has_value else 0.5
            track.rect(int(row["start"]), int(row["end"]), fc=plt.get_cmap("viridis")(frac), ec="none")
        if has_value:
            track.text(c, x=sector.size / 2, r=100, size=10)
    fig = circos.plotfig()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return {"out": out, "chromosomes": list(chroms)}


def main(argv=None):
    ap = argparse.ArgumentParser(description="bioinformatics statistical plots")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_v = sub.add_parser("volcano")
    p_v.add_argument("csv")
    p_v.add_argument("--log2fc"); p_v.add_argument("--pvalue"); p_v.add_argument("--genes")
    p_v.add_argument("--thresholds", nargs=2, type=float, default=[1.0, 0.05])
    p_v.add_argument("--top", type=int, default=15)
    p_v.add_argument("--out", default="volcano.png")
    p_m = sub.add_parser("ma")
    p_m.add_argument("csv")
    p_m.add_argument("--log2fc"); p_m.add_argument("--base-mean"); p_m.add_argument("--genes")
    p_m.add_argument("--out", default="ma.png")
    p_h = sub.add_parser("heatmap")
    p_h.add_argument("csv")
    p_h.add_argument("--zscore", action="store_true")
    p_h.add_argument("--cluster", action="store_true")
    p_h.add_argument("--cmap", default="RdYlBu_r")
    p_h.add_argument("--out", default="heatmap.png")
    p_c = sub.add_parser("circos")
    p_c.add_argument("csv")
    p_c.add_argument("--out", default="circos.png")
    args = ap.parse_args(argv)
    try:
        if args.cmd == "volcano":
            r = volcano(args.csv, log2fc_col=args.log2fc, pvalue_col=args.pvalue,
                        gene_col=args.genes, lfc_thr=args.thresholds[0],
                        p_thr=args.thresholds[1], top_n=args.top, out=args.out)
            print("up=%d down=%d -> %s" % (r["up"], r["down"], r["out"]))
        elif args.cmd == "ma":
            r = ma(args.csv, log2fc_col=args.log2fc, base_mean_col=args.base_mean,
                   gene_col=args.genes, out=args.out)
            print("written: %s" % r["out"])
        elif args.cmd == "heatmap":
            r = heatmap(args.csv, zscore=args.zscore, cluster=args.cluster,
                        cmap=args.cmap, out=args.out)
            print("shape %s -> %s" % (r["shape"], r["out"]))
        else:
            r = circos(args.csv, out=args.out)
            print("chromosomes %s -> %s" % (r["chromosomes"], r["out"]))
    except Exception as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
