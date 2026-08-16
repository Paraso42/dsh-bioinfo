---
name: bio-data-hub
description: "公共生物数据库 REST 客户端(纯标准库、内置重试与限速):UniProt(注释/FASTA/检索)、RCSB PDB(元数据/结构/序列/检索)、STRING(互作网络/ID 映射)、KEGG(find/get/link);入口 uniprot_fetch.py / pdb_fetch.py / string_fetch.py / kegg_fetch.py"
---

# 生物数据获取中心(bio-data-hub)

四个纯标准库 REST 客户端,覆盖结构/序列/互作/通路四大数据库。均内置浏览器 UA、退避重试、超时与输出落盘;可直接用系统 Python 跑(无第三方依赖)。

## 〇、环境

- 解释器:`C:\Program Files\Python313\python.exe`(纯标准库;venv 也可)
- 脚本目录:`C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\bio-data-hub\resources\`
- 已在线实测:UniProt get/search/fasta、PDB download/meta(1YPH/3PTB)、STRING map/network(TP53)、KEGG find/get/link

## 一、uniprot_fetch.py — UniProt

```powershell
$U = "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\bio-data-hub\resources\uniprot_fetch.py"
& 'C:\Program Files\Python313\python.exe' $U get P04637 --out p53.json          # 注释摘要 + 全量 JSON
& 'C:\Program Files\Python313\python.exe' $U fasta P00648 P00649 --out seqs.fasta
& 'C:\Program Files\Python313\python.exe' $U search "organism:9606 AND reviewed:true AND keyword:Kinase" `
    --size 10 --fields accession,id,gene_names,protein_name --format tsv --out hits.tsv
```

- 查询语法自动纠错:`organism:` → `organism_id:`(taxonomy id);其他 Solr 字段照写
- **按物种名检索用 `organism_name:`**(如 `organism_name:"Riccia fluitans"`);**经 PowerShell 传参时整段查询用单引号包裹**(内部双引号保留):`$U search 'organism_name:"Riccia fluitans" AND reviewed:true' ...` —— 用双引号嵌套会导致参数碎裂(`unrecognized arguments: fluitans`)
- `--format list` 只出 accession 列表(便于管道);`--format json` 全量
- 429/5xx 自动指数退避重试(最多 5 次)

## 二、pdb_fetch.py — RCSB PDB

```powershell
$P = "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\bio-data-hub\resources\pdb_fetch.py"
& 'C:\Program Files\Python313\python.exe' $P meta 3ptb --json 3ptb_meta.json   # 标题/分辨率/方法/链/配体
& 'C:\Program Files\Python313\python.exe' $P download 3ptb --format pdb --out 3ptb.pdb
& 'C:\Program Files\Python313\python.exe' $P download 3ptb --format cif --out 3ptb.cif
& 'C:\Program Files\Python313\python.exe' $P fasta 3ptb --out 3ptb.fasta
& 'C:\Program Files\Python313\python.exe' $P search '{"query":{"type":"terminal","service":"text","parameters":{"value":"trypsin","attribute":"struct.title.pdbx"}},"return_type":"entry"}' --out hits.json
```

- meta 的链/配体信息:core/entry 端点无链数据时自动从结构文件解析兜底(`source` 字段标注)
- meta 输出含 **`chain_residues_polymer`**:逐链聚合物残基数(结构文件 ATOM 计数、去结晶水,只取第一个 MODEL)——大/小亚基一眼分辨(如 8RUC:783 vs 207)
- search 接受 RCSB 查询 JSON 对象或 `@file.json`

## 三、string_fetch.py — STRING-DB(互作网络)

```powershell
$S = "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\bio-data-hub\resources\string_fetch.py"
& 'C:\Program Files\Python313\python.exe' $S network P04637 --species 9606 --score 900 --limit 50 --out tp53_net.json
& 'C:\Program Files\Python313\python.exe' $S map P04637 ENSG00000141510 --species 9606   # 外部 ID → STRING ID
& 'C:\Program Files\Python313\python.exe' $S image P04637 MDM2 --species 9606 --score 700 --out net.png
```

- network 输出:partner A/B、综合 score、experiments/database/textmining 分项
- **计数与显示名 caveat**:打印的 `partners: N` 是 `--limit`(默认 50)截断后的实际行数,**不是**该节点的全部互作数;`preferredName` 可能是 STRING 的怪异别名(如 P10896 显示 `MTI20.21` ≠ 蛋白名)——**以 map 的映射行为准,显示名仅供参考**
- map 逐 ID 请求(v12 端点不接受多行 ID)并自动限速 1s
- 注意:跨物种查询会返回空结果(如鸡 P00698 查 9606),先确认 taxonomy

## 四、kegg_fetch.py — KEGG

```powershell
$K = "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\bio-data-hub\resources\kegg_fetch.py"
& 'C:\Program Files\Python313\python.exe' $K find genes barnase --org bam --limit 10       # 基因检索
& 'C:\Program Files\Python313\python.exe' $K find pathway apoptosis --org hsa
& 'C:\Program Files\Python313\python.exe' $K get hsa:5594 hsa:207 --out kegg_entries.txt    # 条目(自动分批)
& 'C:\Program Files\Python313\python.exe' $K link pathway hsa:5594 --out mapk1_paths.json   # 基因→通路
& 'C:\Program Files\Python313\python.exe' $K link genes path:hsa04210                         # 通路→基因
```

- `--org` 用 3-4 字母生物代码(bsu/bam/eco/hsa…),find 时以 `+` 拼接
- get/link 单请求 ≤10 ID(自动分批 + 0.6s 限速)

## 五、常见错误

| 报错/现象 | 原因与解法 |
|---|---|
| UniProt `'organism' is not a valid search field` | 已自动转 `organism_id:`;其他字段查 UniProt 帮助页 |
| STRING map 空结果 | 物种不匹配(先查 taxonomy);或 ID 类型不在 STRING 覆盖内 |
| KEGG `HTTP 400 /find/...` | org 写法:`find genes query+org`(脚本已处理);org 代码查 KEGG 生物列表 |
| RCSB meta chains/ligands 为空 | 自动回退解析结构文件;仍空则该条目确无聚合物链 |
| PowerShell 下 uniprot search `unrecognized arguments`(如 fluitans) | 查询串内嵌双引号与外层双引号冲突导致参数碎裂;整段查询用**单引号**包裹、内部用双引号(见第一节) |
| STRING network 显示名是怪异别名(如 MTI20.21) | STRING `preferredName` ≠ 蛋白名;以 map 映射行/stringId 为准 |
| 429 反复 | 脚本已退避;降低并发(这些端点 1 req/s 量级足够) |
