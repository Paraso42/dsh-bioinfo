---
name: biopython
description: "本机 Biopython 1.87 核心技能：序列处理与文件 I/O（环境引导、Seq/SeqIO/SeqRecord、格式转换、GC 与翻译、常见坑）"
---
# Biopython 核心技能（本机环境 · 序列处理与文件 I/O）

本技能教会你熟练调用本机标准生信 Python 代码库（Biopython 1.87）解决 DNA/RNA/蛋白质的**序列处理、文件读写、格式转换、序列统计**类问题。所有 API 均在本机实测验证（Biopython 1.87）。

## 一、环境引导（写任何代码前必读）

本机生信环境（固定路径，已在 DSH 中验证）：

| 项目 | 值 |
|------|-----|
| Python 解释器 | `C:\Program Files\Python313\python.exe`（**必须是 3.13**） |
| Biopython | 1.87，安装于 `D:\biopython`（含 `Bio/`、`BioSQL/`） |
| numpy | 2.4.6（随 `D:\biopython` 捆绑，cp313 编译版） |
| 可选依赖（已装） | matplotlib 3.10.3、reportlab 4.5.1、scipy 1.15.3、networkx 3.6.1、pandas 2.2.3 |

**⚠️ 最重要的坑：** 默认 `python`（3.14）导入 `Bio` 会失败（`numpy._core._multiarray_umath` 报错，捆绑 numpy 是 cp313 版）。**必须用 Python 3.13**。不要试图用默认 python 重装 numpy 去"修复"——直接用 3.13。

### 运行方式（通过 pwsh 工具）

```powershell
# 方式 A：写脚本文件再运行（推荐，多行代码）
# 脚本开头固定两行：
#   import sys; sys.path.insert(0, r"D:\biopython")
& "C:\Program Files\Python313\python.exe" "你的脚本.py"

# 方式 B：单行
& "C:\Program Files\Python313\python.exe" -c "import sys; sys.path.insert(0, r'D:\biopython'); from Bio.Seq import Seq; print(Seq('ATGC').complement())"

# 方式 C：设置 PYTHONPATH 后运行（等效）
$env:PYTHONPATH = "D:\biopython"
& "C:\Program Files\Python313\python.exe" script.py
```

### 测试数据与文档（验证脚本非常好用）

- 测试数据：`D:\biopython\code\Tests\` 下有大量真实格式文件：`GenBank\NC_005816.gb`（细菌基因组 GenBank）、`Quality\*.fastq`、`PDB\1A8O.pdb` / `1A8O.cif`、`SwissProt\P62258.txt`、`TreeConstruction\msa.phy`、`PAML\Trees\lysin.trees` 等。
- 源码/文档：`D:\biopython\code\Doc\`（Tutorial 源文件）、`D:\biopython\code\Scripts\`（实用脚本）。
- 写分析脚本后，建议先用 `code\Tests` 里的样例数据跑通，再上真实数据。

## 二、Bio.Seq：序列对象与操作

```python
import sys; sys.path.insert(0, r"D:\biopython")
from Bio.Seq import Seq

dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
dna.complement()            # 互补（保留大小写）
dna.reverse_complement()    # 反向互补
rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
rna.transcribe()            # DNA->RNA（仅 T->U）；Seq("ATGC").transcribe()
Seq("ATG").back_transcribe()  # RNA->DNA
dna.translate()             # 蛋白质翻译（默认标准密码子表，遇终止密码子停止）
dna.translate(to_stop=True) # 翻译到第一个终止密码子即停
dna.translate(table="Bacterial")  # 指定遗传密码子表（"Bacterial"=细菌/线粒体等）
# translate(cds=True) 会校验起始/终止密码子，若序列不是完整 CDS 会抛 TranslationError
seq.lower(); seq.upper()
seq.count("AT"); len(seq)
seq.replace("-", "")        # 去 gap（注意：1.87 无 ungap 方法，用 replace）
seq.split("G"); seq.join([Seq("AA"), Seq("TT")])
Seq("ACGT") == Seq("ACGT")  # 按字母比较
```

**注意 1.87 变更：** `Seq.ungap()` 已移除，用 `seq.replace("-", "")`；`Bio.SeqUtils.GC()` 已移除，用 `gc_fraction()`（见下文）。

## 三、SeqRecord / SeqFeature：带注释的序列记录

```python
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

rec = SeqRecord(Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"),
                id="gene1", name="gene1", description="example gene")
rec.annotations["molecule_type"] = "DNA"   # GenBank/EMBL 写盘必需！
rec.annotations["organism"] = "Homo sapiens"
rec.features.append(SeqFeature(FeatureLocation(0, 15), type="gene",
                               qualifiers={"gene": ["gene1"]}))
# 拼接型位置（join 型 CDS）：
join_loc = CompoundLocation([FeatureLocation(0, 3), FeatureLocation(5, 8)], operator="join")
```

`rec.id`、`rec.name`、`rec.description`、`rec.seq`、`rec.features`、`rec.annotations`、`rec.letter_annotations`（FASTQ 质量值放这里）、`rec.dbxrefs`、`rec.format("fasta")` 是最常用字段/方法。

## 四、Bio.SeqIO：序列文件读写（核心中的核心）

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 读取（parse=多记录迭代器；read=单记录）
for rec in SeqIO.parse("input.fasta", "fasta"):
    print(rec.id, len(rec.seq))
rec = SeqIO.read("single.gb", "genbank")       # 仅一条记录时用 read，多条会报错
records = list(SeqIO.parse("reads.fastq", "fastq"))

# 写入
SeqIO.write(records, "out.fasta", "fasta")
# FASTQ 写入需要每条记录有质量值：
rec.letter_annotations["phred_quality"] = [30, 30, 30, 30]  # 缺质量值写 fastq 会 ValueError

# 格式转换（一行搞定；注意转 GenBank/EMBL 需要 molecule_type 注释）
SeqIO.convert("in.fasta", "fasta", "out.gb", "genbank")

# 索引大文件（内存友好，惰性访问；参数必须是文件路径，不能是句柄）
idx = SeqIO.index("big.fasta", "fasta")
seq = idx["some_id"]            # 按 id 随机访问
SeqIO.index_db("big.idx", ["big.fasta"], "fasta")  # 持久化索引，跨进程可用

# 字典化（小文件）
d = SeqIO.to_dict(SeqIO.parse("in.fasta", "fasta"))

# 支持格式（1.87 实测）：fasta, fastq, fastq-sanger, fastq-solexa, fastq-illumina,
#   genbank/gb, embl, ig, swiss, uniprot-xml, phd, qual, sff, abi, ace, tab,
#   gff3(部分), clustal, emboss, fasta-2line, fasta-blast, fasta-m10, fasta-pearson, ...
```

**1.87 注意：** `SeqIO.format()` 已移除（用 `SeqIO.write` 到 `io.StringIO()`，或 `rec.format("fasta")`）；`SeqIO.convert` 转 GenBank/EMBL 前务必给记录补 `annotations["molecule_type"]`，否则报 `ValueError: missing molecule_type`。

## 五、序列统计与工具（Bio.SeqUtils）

```python
from Bio.SeqUtils import gc_fraction, molecular_weight, nt_search
from Bio.Seq import Seq

gc_fraction(Seq("ATGCGC"))              # 0.666...（1.87 的新名字；GC() 已移除）
molecular_weight(Seq("ATGC"))           # 核酸分子量
nt_search("ACGTACGTACG", "ACG")         # 查找模式，返回 [模式, 位置...]
# 蛋白质理化性质（需 Bio.SeqUtils.ProtParam 和蛋白质序列）：
from Bio.SeqUtils.ProtParam import ProteinAnalysis
pa = ProteinAnalysis("MAIVMGRKGAR")
pa.molecular_weight(); pa.gravy(); pa.aromaticity(); pa.isoelectric_point()
pa.count_amino_acids(); pa.get_amino_acids_percent()
```

## 六、翻译与遗传密码子表（Bio.Data）

```python
from Bio.Data import CodonTable, IUPACData
std = CodonTable.unambiguous_dna_by_name["Standard"]
std.forward_table["ATG"]      # 'M'
std.stop_codons               # ['TAA', 'TAG', 'TGA']
std.back_table["M"]           # 编码该氨基酸的密码子
CodonTable.unambiguous_rna_by_name["Vertebrate Mitochondrial"]  # 线粒体表
IUPACData.ambiguous_dna_values  # 简并碱基字典，如 'M' -> 'AC'
```

## 七、日常任务模板（直接可用）

### 1) 统计 FASTA 中每条序列的 GC 含量
```python
import sys; sys.path.insert(0, r"D:\biopython")
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
for rec in SeqIO.parse("sequences.fasta", "fasta"):
    print(rec.id, round(gc_fraction(rec.seq) * 100, 2))
```

### 2) FASTQ → FASTA（可配合质量过滤）
```python
from Bio import SeqIO
with open("clean.fasta", "w") as out:
    for rec in SeqIO.parse("reads.fastq", "fastq"):
        if min(rec.letter_annotations["phred_quality"]) >= 20:  # Q20 过滤
            SeqIO.write(rec, out, "fasta")
```

### 3) 从 GenBank 提取全部 CDS 翻译成蛋白
```python
from Bio import SeqIO
rec = SeqIO.read("NC_005816.gb", "genbank")
for f in rec.features:
    if f.type == "CDS" and f.location:
        cd = f.extract(rec.seq)
        prot = cd.translate(table="Bacterial", to_stop=True)
        print(f.qualifiers.get("gene", ["?"])[0], str(prot)[:40])
```

### 4) 多条序列批量翻译/反向互补后输出
```python
from Bio import SeqIO
from Bio.Seq import Seq
recs = []
for rec in SeqIO.parse("dna.fasta", "fasta"):
    rec.seq = rec.seq.reverse_complement()
    recs.append(rec)
SeqIO.write(recs, "rc.fasta", "fasta")
```

### 5) 查找限制性内切酶位点（详见 biopython-analyses skill）
```python
from Bio.Restriction import EcoRI
from Bio.Seq import Seq
EcoRI.search(Seq("GAATTCAAGCTT"))   # [2]  —— 切位点索引
```

## 八、常见错误速查

| 报错 | 原因与解法 |
|------|-----------|
| `ModuleNotFoundError: No module named 'numpy._core._multiarray_umath'` | 用了默认 python 3.14；改用 `C:\Program Files\Python313\python.exe` |
| `ValueError: missing molecule_type in annotations` | SeqIO 写 GenBank/EMBL 前补 `rec.annotations["molecule_type"]="DNA"` |
| `SeqIO.write` 写 fastq 报无质量值 | 给每条记录加 `letter_annotations["phred_quality"]` |
| `AttributeError: 'Seq' object has no attribute 'ungap'` | 1.87 用 `seq.replace("-", "")` |
| `ImportError: cannot import name 'GC' from 'Bio.SeqUtils'` | 1.87 用 `gc_fraction` |
| `AttributeError: module 'Bio.SeqIO' has no attribute 'format'` | 用 `SeqIO.write` 到 StringIO 或 `rec.format()` |
| `SeqIO.index` 报 TypeError | 参数必须是文件路径字符串，不能是文件句柄 |

## 九、什么时候用这个 skill

- 序列操作：互补/翻译/剪接/去 gap
- 任何 FASTA/FASTQ/GenBank/EMBL 文件的读写、转换、过滤、统计
- 引物/探针的 GC、Tm、反向互补计算
- 序列注释（features）提取

若任务是**多序列比对、进化树、BLAST/Entrez 数据库、蛋白质结构(PDB)、PAML 选择压力分析**，请加载 `biopython-analyses` 技能。
