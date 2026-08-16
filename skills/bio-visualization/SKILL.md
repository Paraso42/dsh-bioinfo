---
name: bio-visualization
description: "生信统计可视化(matplotlib/seaborn/pyCirclize/logomaker,环境 D:\\bioai\\venv):火山图/MA 图/聚类热图/环形基因组图(stat_plots.py)+ 序列 Logo 概率型与信息量型(seq_logo.py);CSV 约定与实测示例见文档"
---

# 生信统计可视化(bio-visualization)

一张图顶千行表:差异表达、样本聚类、基因组区段、序列保守性的标准出图。

## 〇、环境

- 解释器:`D:\bioai\venv\Scripts\python.exe`(matplotlib 3.11 + seaborn + pycirclize 1.10 + logomaker 0.8.7 + pandas 3)
- 脚本目录:`C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\bio-visualization\resources\`
- 全部 `Agg` 后端无头出图,PNG dpi=150;已实测四个统计图 + 两种 Logo

## 一、stat_plots.py — 统计图

### volcano 火山图

```powershell
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\stat_plots.py" volcano deg.csv `
    --log2fc log2FoldChange --pvalue padj --genes gene `
    --thresholds 1 0.05 --top 15 --out volcano.png
```
- 列名可省略(自动识别 log2FoldChange/log2FC/padj/pvalue 等常见名);输出 up/down 计数
- 阈值默认 |log2FC|≥1 且 padj<0.05;前 N 显著基因自动标注

### ma MA 图
```powershell
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\stat_plots.py" ma deg.csv --out ma.png
```

### heatmap 热图(可选 z-score + 双向聚类)
```powershell
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\stat_plots.py" heatmap matrix.csv --zscore --cluster --cmap RdYlBu_r --out heatmap.png
```
- CSV 首列 = 行名,其余列数值;`--cluster` 用 seaborn clustermap(行/列双向)

### circos 环形基因组图
```powershell
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\stat_plots.py" circos sectors.csv --out circos.png
```
- CSV 列:`chrom,start,end[,value]`(BED 风格);value 用于 viridis 颜色映射;已适配 pycirclize 1.10 API

## 二、seq_logo.py — 序列 Logo

输入**已比对**的等长 FASTA(不等长会明确报错):

```powershell
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\seq_logo.py" aligned.fasta --out logo.png
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\seq_logo.py" aligned.fasta --info --first 60 `
    --title "barnase-family" --out logo_bits.png
```
- 默认概率(频率)型;`--info` 信息量型(blosum62 背景的 bits)
- `--first N` 只画前 N 列(长 MSA 建议 40-80);`--start` 指定起始列
- 比对来源建议:MUSCLE/MAFFT(未本地部署)或 Biopython `PairwiseAligner` 星形比对(本技能实测配方:
  `res.coordinates` 逐段展开——注意 1.87 的 `alignment[1]` 返回原始序列而非比对行,勿直接 str())

## 三、实测数据

- `D:\bioai\jobs\acceptance\`:volcano.png(40 up/20 down)、ma.png、heatmap.png(24×6)、circos.png(3 染色体)、
  barnase_logo_bits.png / barnase_logo_prob.png(6 条 barnase 家族序列,157 列,含保守的 VINTFDGVADYL 基序)

## 四、常见错误

| 报错/现象 | 原因与解法 |
|---|---|
| `Failed to convert value(s) to axis units`(circos) | pycirclize 1.10 `text(text, x=, r=)` 参数序;脚本已适配,勿手改签名 |
| `Glyph.__init__() got an unexpected keyword argument 'stack_width'` | logomaker 0.8.7 用 `width=`;脚本已适配 |
| `iloc cannot enlarge its target object` | pandas 3.0 标量构造行为变化;脚本已改用 numpy 计数,勿回退 |
| 序列不等长报错 | seq_logo 要求比对好的输入;先比对再出图 |
| 中文字体方块 | 未配置中文字体;标题/标签用英文或 ASCII 规避 |
