---
name: biopython-analyses
description: "本机 Biopython 1.87 进阶分析技能：PairwiseAligner 比对、AlignIO、系统发育 Phylo/TreeConstruction、蛋白质结构 PDB、BLAST/SearchIO/Entrez、限制酶、模体、PAML 选择压力、图形输出"
---
# Biopython 进阶分析技能（比对 · 进化 · 结构 · 数据库 · 高级分析）

本技能教你熟练调用本机 Biopython 1.87 完成**序列比对、系统发育、蛋白质结构、BLAST/Entrez 数据库、限制酶、模体、选择压力(PAML)**等分析。所有 API 均在本机实测验证。

## 〇、环境（与 biopython skill 相同，速查）

- 解释器：`C:\Program Files\Python313\python.exe`（**不是**默认 python 3.14！）
- 库路径：脚本开头 `import sys; sys.path.insert(0, r"D:\biopython")`，或 `$env:PYTHONPATH="D:\biopython"`
- 运行：`& "C:\Program Files\Python313\python.exe" 脚本.py`
- 可选依赖已装：matplotlib（绘图）、reportlab（基因组图）、scipy、networkx、pandas
- 测试数据：`D:\biopython\code\Tests\`（含 `PDB\1A8O.pdb`、`PhyloXML`、`TreeConstruction\msa.phy`、`PAML\Trees\lysin.trees`、`Blast\` 等）

## 一、双序列比对：Bio.Align.PairwiseAligner（pairwise2 已弃用！）

**1.87 关键：** `Bio.pairwise2` 仍存在但已**弃用**（会告警，未来版本移除）；一律用新式 `Bio.Align.PairwiseAligner`。

```python
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner()
aligner.mode = "global"          # "global" 全局比对 | "local" 局部比对
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -2      # 缺口开放罚分（负值）
aligner.extend_gap_score = -0.5  # 缺口延伸罚分（负值）
# 也可以直接加载计分矩阵：
from Bio.Align import substitution_matrices
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

score = aligner.score("ACCGT", "ACG")      # 最优得分
alignments = aligner.align("ACCGT", "ACG") # 全部最优比对（迭代器/列表）
a = alignments[0]
print(a)                                   # 打印比对
print(a.aligned)                           # 坐标映射 ((target区间),(query区间))
# 蛋白质局部比对示例：
local = PairwiseAligner(); local.mode = "local"
local.substitution_matrix = substitution_matrices.load("BLOSUM62")
local.open_gap_score = -5; local.extend_gap_score = -1
hits = local.align("MGHQQLYWSHPRKFGQGSRS", "MERLVLKS")
print(hits[0]); print(hits[0].score)
```

## 二、多序列比对：Bio.AlignIO + MultipleSeqAlignment

```python
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 读取（支持格式：clustal, fasta, phylip, phylip-relaxed, nexus, stockholm,
#   emboss, msf, maf, mauve, fasta-m10 等——1.87 实测列表）
aln = AlignIO.read("msa.phy", "phylip")          # 单比对块
alns = AlignIO.parse("many.sth", "stockholm")    # 多比对块迭代器
print(len(aln), aln.get_alignment_length())      # 行数、列数

# 手工构造 MultipleSeqAlignment：
msa = MultipleSeqAlignment([SeqRecord(Seq("ACGT"), id="a"),
                            SeqRecord(Seq("ACGA"), id="b")])
print(msa[:, 2])      # 取第 2 列（列切片，返回该列字符）
print(msa[0].seq)     # 取第 0 行
# 新式 Alignment 对象（1.80+ 推荐）：
aln2 = msa.alignment                # MultipleSeqAlignment -> Alignment
print(aln2.shape)                   # (行, 列)
print(aln2[:, 1])                   # 列切片
# AlignInfo.SummaryInfo 已弃用（1.87 实测无共识方法），共识序列用列统计自己算：
from collections import Counter
consensus = "".join(Counter(aln2[:, i]).most_common(1)[0][0] for i in range(aln2.shape[1]))

# 写入
AlignIO.write(msa, "out.clustal", "clustal")
AlignIO.convert("in.phy", "phylip", "out.fasta", "fasta")
```

## 三、系统发育树：Bio.Phylo + Bio.Phylo.TreeConstruction

```python
import io
from Bio import Phylo

# 读树（newick / nexus / phyloxml）
tree = Phylo.read(io.StringIO("((A:0.1,B:0.2):0.3,C:0.4);"), "newick")
# 或 tree = Phylo.read("tree.nwk", "newick")
tree.count_terminals()            # 叶子数
tree.distance("A", "B")           # 两节点间分支长度
tree.common_ancestor("A", "B")    # 最近共同祖先
[c.name for c in tree.find_clades()]          # 遍历全部节点（内部节点 name=None）
tree.get_nonterminals(); tree.get_terminals()
# 剪枝、重命名、子树提取：
sub = tree.from_clade(tree.common_ancestor("A", "B"))
for c in tree.find_clades(): c.name = (c.name or "inner").upper()

# 从多序列比对构建树（距离法 NJ/UPGMA）：
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
msa = AlignIO.read("msa.phy", "phylip")
calculator = DistanceCalculator("identity")       # 或 "blastn"/自定义矩阵
dm = calculator.get_distance(msa)
constructor = DistanceTreeConstructor(calculator, method="nj")   # "upgma" 亦可
tree2 = constructor.build_tree(msa)
tree2.format("newick")                            # 序列化为字符串
Phylo.write(tree2, "out.nwk", "newick")

# 绘图（matplotlib；无头环境用 Agg 后端 + savefig）：
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8, 5), dpi=100)
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree2, axes=ax)
fig.savefig("tree.png")
# 用 networkx 布局：
Phylo.to_networkx(tree)   # 返回 networkx 图（networkx 3.6.1 已装）
```

## 四、BLAST 与数据库检索

### 在线 BLAST（首选：浏览器 UA 的 NCBI URLAPI；本机 2026-07 实测修复）

⚠️ **必读**：`Bio.Blast.NCBIWWW.qblast`（UA=`BiopythonClient`）与默认 Python urllib 会被 NCBI BLAST 服务端**间歇性直接断连**（`RemoteDisconnected: Remote end closed connection without response`），对 nr 大库 + `entrez_query` 的重请求几乎必现。DNS/TCP/GET 均正常——是服务端按 User-Agent 拒连，不是防火墙问题；Entrez（eutils）不受影响。小库（swissprot）偶发能过时 `qblast` 仍可作快捷路径，失败即切回本方案。

**修复方案：带浏览器 UA 走 BLAST URLAPI（CMD=Put → 轮询 CMD=Get）+ 重试。** 现成脚本已随本技能分发（纯标准库，不依赖 Bio）：

- 脚本：`C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\biopython-analyses\resources\ncbi_blast.py`
- 命令行：
  `& "C:\Program Files\Python313\python.exe" "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\biopython-analyses\resources\ncbi_blast.py" blastp nr query.fasta --entrez-query "Viridiplantae[ORGN]" --expect 1e-10 --hits 10 -o result.xml`
- 编程调用：
```python
import sys; sys.path.insert(0, r"C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\biopython-analyses\resources")
from ncbi_blast import ncbi_blast
xml = ncbi_blast("blastp", "nr", ">q\nMSPQTETKAGAGF...", entrez_query="Viridiplantae[ORGN]", expect=1e-10, hitlist=10)
```
- 核心内联实现（逻辑同脚本；Put 失败放 `for attempt in range(6)` + `time.sleep` 重试）：
```python
import re, time
from urllib.request import Request, urlopen
from urllib.parse import urlencode
URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36")
put = {"CMD":"Put","PROGRAM":"blastp","DATABASE":"nr","QUERY":query_fasta,
       "EXPECT":"1e-10","HITLIST_SIZE":"10","ALIGNMENTS":"10","DESCRIPTIONS":"10",
       "FORMAT_TYPE":"XML","MATRIX_NAME":"BLOSUM62","GAPCOSTS":"11 1","FILTER":"F",
       "ENTREZ_QUERY":"Viridiplantae[ORGN]","email":"your@email.com","tool":"dsh"}
page = urlopen(Request(URL, urlencode(put).encode(), {"User-Agent": UA}), timeout=90).read().decode("utf-8", "replace")
rid = re.search(r"RID = (\S+)", page).group(1)
get = {"CMD":"Get","FORMAT_TYPE":"XML","RID":rid,"ALIGNMENTS":"10","DESCRIPTIONS":"10"}
while True:                       # 轮询：20 s 起，之后 60 s（NCBI 官方限频）
    body = urlopen(Request(URL + "?" + urlencode(get), headers={"User-Agent": UA}), timeout=90).read().decode("utf-8", "replace")
    if "Status=" not in body: break   # 就绪：最终 XML 无 Status 标签
    time.sleep(60)
# body 即 BLAST XML
```
- 限频纪律：两次 Put 间隔 ≥10 s、同一 RID 轮询 ≥60 s；**不要**并发轰炸。
- NCBI FTP（ftp.ncbi.nlm.nih.gov）大文件被带宽整形（~1 KB/s 级）：改走 Datasets API（`https://api.ncbi.nlm.nih.gov/datasets/v2alpha/...`）或 Entrez efetch 的 `seq_start`/`seq_stop` 分段下载。

### 解析 BLAST 结果（SearchIO；1.87 注意属性位置）

```python
from Bio import SearchIO
from io import StringIO
q = SearchIO.read(StringIO(xml), "blast-xml")   # 单 query；多 query 用 parse()
for hit in q[:10]:
    hsp = hit[0]                # evalue/bitscore/identity 在 hsp 上，不在 hit 上！
    ident = hsp.ident_num / max(hsp.aln_span, 1) * 100
    cov = hsp.hit_span / max(hit.seq_len, 1) * 100
    print(hit.id, "evalue=", hsp.evalue, f"ident={ident:.1f}% cov={cov:.0f}%", hit.description)
```

### 解析本地 BLAST 输出（SearchIO，替代 NCBIXML）
```python
from Bio import SearchIO
# 支持 blast-xml / blast-tab（-outfmt 6）/ blast-text 等
qres = SearchIO.read("out.xml", "blast-xml")          # 单条 query
qresults = list(SearchIO.parse("out.tsv", "blast-tab"))  # 多条 query（tab 分隔）
for q in qresults:
    print(q.id, "hits:", len(q))
    for hit in q[:5]:
        print("  ", hit.id, "evalue=", hit[0].evalue, "bitscore=", hit[0].bitscore)
# 注意：blast-tab 文件若含多个 query 必须用 parse()，用 read() 会报错
```

### Entrez 数据库（NCBI，需要网络；必须设置邮箱）
```python
from Bio import Entrez
Entrez.email = "your@email.com"      # NCBI 要求，否则可能被限流/封禁
handle = Entrez.esearch(db="nucleotide", term="Homo sapiens[Organism] AND COX1[Gene]", retmax=10)
record = Entrez.read(handle)
print(record["Count"], record["IdList"])
# 下载序列：
h2 = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta", retmode="text")
from Bio import SeqIO
rec = SeqIO.read(h2, "fasta")
# efetch genbank：rettype="gbwithparts"；蛋白质：db="protein"
# 注意 NCBI 限流：esearch/efetch 间 sleep 0.34 秒以上
import time; time.sleep(0.4)
```

## 五、蛋白质结构：Bio.PDB

```python
from Bio.PDB import PDBParser, MMCIFParser, PDBList, Selection, Superimposer
from Bio.PDB import PDBIO, Select

p = PDBParser(QUIET=True)                       # QUIET 抑制警告
struct = p.get_structure("1A8O", r"D:\biopython\code\Tests\PDB\1A8O.pdb")
# 或 mm = MMCIFParser(QUIET=True); struct = mm.get_structure("1A8O", "1A8O.cif")
model = struct[0]                               # 第一个模型
for chain in model:                             # 遍历链
    for res in chain:                           # 遍历残基
        res.id[1]; res.resname                  # 残基序号、残基名
for atom in struct.get_atoms():                 # 遍历原子
    atom.name; atom.coord; atom.get_bfactor()
# 选择/提取（原子级选择器）：
sel = Selection.unfold_entities(struct, "A")    # 取全部原子
# 以残基范围取子链（1.87 注意：Chain 不支持切片索引，用 get_list()）：
for res in model["A"].get_list()[10:20]:
    print(res.resname, res.id[1])

# 结构比对（Superimposer，计算 RMSD）：
sup = Superimposer()
atoms1 = [a for a in model["A"].get_atoms()][:50]
sup.set_atoms(atoms1, atoms1)     # 需要两套等价原子
print("RMSD:", sup.rms)           # 平移/旋转在 sup.rotran

# 溶剂可及表面积（SASA）：
from Bio.PDB.SASA import ShrakeRupley
sasa = ShrakeRupley()
sasa.compute(model, level="R")    # level: A(原子)/R(残基)/C(链)
for res in model["A"]:
    print(res.resname, res.id[1], res.sasa)

# 二级结构（DSSP，需要外部 dssp 可执行文件，本机未必有；无则跳过）：
# from Bio.PDB import DSSP; dssp = DSSP(model, pdb_file)

# 写回 PDB（配合 Select 过滤）：
class ChainA(Select):
    def accept_chain(self, chain): return chain.id == "A"
io = PDBIO(); io.set_structure(struct); io.save("chainA.pdb", ChainA())

# 从网上下载 PDB：
# pdbl = PDBList(); pdbl.retrieve_pdb_file("1A8O", pdir=".", file_format="pdb")
```

## 六、限制性内切酶分析：Bio.Restriction

```python
from Bio.Restriction import EcoRI, BamHI, HindIII, Analysis, RestrictionBatch, AllEnzymes
from Bio.Seq import Seq

seq = Seq("GAATTCAAGCTTGGATCC")
EcoRI.search(seq)                    # 识别位点位置列表，如 [2]
EcoRI.catalyse(seq)                  # 酶切后的片段
# 多酶批量分析：
rb = RestrictionBatch([EcoRI, BamHI, HindIII])
an = Analysis(rb, seq)
an.print_as("map")                   # 打印酶切图谱
# AllEnzymes：内置全部 1088 种酶（1.87 实测）
# 引物/序列设计辅助：Bio.SeqUtils.nt_search 找简并位点
```

## 七、模体（motif）分析：Bio.motifs

```python
from Bio import motifs
from Bio.Seq import Seq

# 由一组实例创建 PWM/PFM：
m = motifs.create([Seq("ACGT"), Seq("ACGA"), Seq("ACGT")])   # 1.87：传 Seq 列表
m.counts                      # 位置计数矩阵
m.pwm                         # 频率矩阵
m.consensus                   # 共识序列
# 在序列上搜索模体（1.87 注意：m.search 已移除，用 pssm.search）：
for pos, score in m.pssm.search(Seq("ACGTTTTACGT"), threshold=0.0, both=True):
    print(pos, score)         # 负 pos = 反向链命中；both=False 只看正向
# 从文件解析（MEME/JASPAR/TRANSFAC 格式）：
# m = motifs.parse(open("meme.txt"), "MEME")
# 读取 JASPAR 数据库条目：m = motifs.read(open("MA0004.1.jaspar"), "jaspar")
```

## 八、选择压力分析：Bio.Phylo.PAML（Bio.PAML 已移除！）

**1.87 关键：** `Bio.PAML` 模块已**移除**，全部使用 `Bio.Phylo.PAML`（codeml / yn00 / baseml 包装器）。运行需要外部程序 `codeml`/`yn00`/`baseml` 在本机 PATH 中（本机未必装有，写代码时先 `Get-Command codeml` 确认，没有则告知用户）。

```python
from Bio.Phylo.PAML import codeml, yn00, baseml

cml = codeml.Codeml(alignment="aligned.fasta", tree="tree.nwk", out_file="codeml.out")
cml.read_ctl_file("codeml.ctl")      # 读控制文件；或直接设选项（关键字参数，注意不是 dict）：
cml.set_options(NSsites=[0, 1, 2, 7, 8], model=0, seqtype=1)
cml.run(verbose=True)
results = codeml.read("codeml.out")  # 解析结果（参数估计、似然值等）
# yn00（核苷酸替代速率 dN/dS 两两比较）：
yn = yn00.Yn00(alignment="aln.fasta", tree="tree.nwk", out_file="yn00.out")
yn.set_options({"verbose": 1}); yn.run(verbose=True)
# 测试数据见 D:\biopython\code\Tests\PAML\（Control_files、Trees\lysin.trees）
```

## 九、图形输出

```python
# 1) 基因组圈图/线性图（GenomeDiagram，需 reportlab，已装）：
from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
rec = SeqIO.read("NC_005816.gb", "genbank")
gd = GenomeDiagram.Diagram("genome")
track = gd.new_track(1, name="CDS")
feat = track.new_set()
for f in rec.features:
    if f.type == "CDS":
        feat.add_feature(f, color=colors.blue, label=True)
gd.draw(format="linear", orientation="landscape", pagesize="A4")
gd.write("genome.pdf", "PDF")   # 本机实测：PDF/EPS 正常；PNG 会报 RenderPMError
# （原因：reportlab 缺 rlPyCairo 栅格后端，PNG 输出不可用，输出 PDF 即可）

# 2) 进化树绘图见第三节（Phylo.draw + matplotlib savefig，PNG 正常）
```

## 十、其他常用模块速览

| 模块 | 用途 | 1.87 实测要点 |
|------|------|--------------|
| `Bio.SwissProt` | 解析 SwissProt 文本 | `SwissProt.read(handle)` 返回记录对象 |
| `Bio.ExPASy` | UniProt/SwissProt 在线下载 | 需网络 |
| `Bio.SearchIO` | 检索结果解析（通用） | blast-xml/blast-tab 均支持；多 query 用 parse() |
| `Bio.AlignIO` | 比对文件读写 | 见第二节格式列表 |
| `Bio.SeqFeature` | 注释特征 | CompoundLocation 拼接位置 |
| `Bio.Data.CodonTable` | 遗传密码子表 | 见 biopython skill 第六节 |
| `Bio.Graphics.GenomeDiagram` | 基因组图 | 需 reportlab（已装） |
| `Bio.motifs` | 模体 | `motifs.create([Seq(...)])`；SeqInstances 不存在 |
| `Bio.Phylo` | 系统发育 | 支持 newick/nexus/phyloxml |
| `Bio.Blast.NCBIWWW` | 在线 BLAST（简易） | qblast 会被 NCBI 间歇性拒连（RemoteDisconnected）；首选第四节浏览器 UA URLAPI 方案（`resources/ncbi_blast.py`） |

## 十一、常见错误速查

| 报错 | 原因与解法 |
|------|-----------|
| `BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated` | 改用 `Bio.Align.PairwiseAligner` |
| `ModuleNotFoundError: No module named 'Bio.PAML'` | 1.87 已移除；用 `Bio.Phylo.PAML` |
| `AttributeError: module 'Bio.motifs' has no attribute 'SeqInstances'` | `motifs.create` 直接传 `Seq` 列表 |
| `ValueError: More than one query result found in handle` | 多 query 结果用 `SearchIO.parse()` 而非 `read()` |
| `BiopythonDeprecationWarning: SummaryInfo has been deprecated` | 用 `msa.alignment` 新式对象 + `[:, col]` 切片 |
| `ValueError: Keyword argument "branch_length=None" is not in the format...` | `Phylo.draw(tree, axes=ax)` 不要传 `branch_length=None` |
| `UserWarning: FigureCanvasAgg is non-interactive` | 正常，改用 `fig.savefig()` 而非 `plt.show()` |
| codeml/yn00 运行报找不到程序 | 需要外部 PAML 程序；先 `Get-Command codeml` 确认再写代码 |
| `RemoteDisconnected: Remote end closed connection without response`（qblast / urllib POST BLAST） | NCBI BLAST 服务端拒非浏览器 UA；改用第四节浏览器 UA URLAPI 方案（`resources/ncbi_blast.py`），带重试 |
| `AttributeError: 'Hit' object has no attribute 'evalue'` | SearchIO 1.87 中 evalue/bitscore 在 `hsp` 上（`hit[0].evalue`）；identity 用 `hsp.ident_num / hsp.aln_span` |
| ftp.ncbi.nlm.nih.gov 下载极慢/中断 | NCBI FTP 带宽整形（~1 KB/s 级）；改走 Datasets API 或 Entrez `seq_start`/`seq_stop` 分段 efetch |

## 十二、什么时候用这个 skill

- 双序列/多序列比对、比对文件格式转换
- 构建/解析/可视化系统发育树，dN/dS 等选择压力分析（PAML）
- 在线 BLAST、NCBI Entrez 检索下载、本地 BLAST 结果解析
- 蛋白质三维结构解析（PDB/mmCIF）、RMSD 比对、SASA 计算
- 限制酶切、模体发现、基因组图形绘制
