---
name: protein-quality
description: "预测结构质量评估闭环:TM-score(与官方 TMalign 交叉验证)/CA-与全原子 RMSD/lDDT/GDT-TS/GDT-HA/DockQ(Fnat/iRMS/LRMS)复合物评分 + PRODIGY 结合亲和力(ΔG/Kd);入口 struct_eval.py 与 prodigy_affinity.py,见 D:\\bioai 部署"
---

# 蛋白质结构质量评估与结合亲和力(protein-quality)

预测模型好不好的**定量闭环**:结构级指标(struct_eval.py,与官方工具交叉验证)+ 结合强度(prodigy_affinity.py)。

## 〇、环境

- 解释器:`D:\bioai\venv\Scripts\python.exe`(venv 内含 numpy/scipy/biopython/prodigy-prot/prodigy-lig)
- 脚本目录:`C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\protein-quality\resources\`
- 已实测验收数据:`D:\bioai\jobs\acceptance\`(1brs 晶体与 AF2 预测模型、评估 JSON)

## 一、struct_eval.py — 结构质量评估

```powershell
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\struct_eval.py" `
    --model predicted.pdb --ref native.pdb `
    --complex --model-chains A B --ref-chains A D --out eval.json
```

### 指标与判定

| 指标 | 含义 | 参考阈值 |
|---|---|---|
| TM-score | 结构相似度 0~1(对长度/错配稳健) | >0.5 同折叠;>0.7 高相似;<0.2 无关 |
| CA-RMSD | 叠加后对齐核心的 CA 均方根偏差(Å) | <1 优秀;<2 良好;<5 尚可 |
| lDDT | 局部距离差检验(叠加无关,现代 CASP 主指标) | ≥0.7 高置信;0.5-0.7 中;<0.5 低 |
| GDT_TS / GDT_HA | 全局距离检验(CASP 经典) | 越高越好 |
| coverage / seq_identity | 映射对上的覆盖与一致率;**homology 模式(默认)= 全长同源比对统计**;identical 模式 = 相同残基对子集 | 诊断用 |
| DockQ(复合物) | 对接质量 0~1(Fnat+1/(1+(iRMS/1.5)²)+1/(1+(LRMS/8.5)²))/3 | ≥0.23 可接受;≥0.49 中;≥0.80 高(CAPRI 标准) |

输出 JSON 含逐链 lDDT、残基映射(ref→model 编号)、Fnat/iRMS/LRMS 与原生接触数。

### 关键行为

- **残基自动映射**:模型与参考编号不一致(colabfold 输出从 1 重编号)时自动按序列比对映射,无需预处理
- **映射模式 `--mapping`**(`auto` 默认 / `homology` / `identical`):
  - 序列相同时直接对角映射;序列不同时默认 **homology**——全长同源全局比对(match=2/mismatch=-1/gap=-2/-0.5),coverage/seq_identity 反映真实全长一致率(远缘同源物 ~55% 一致时覆盖≈100%、id≈55%)
  - **解读 caveat(重要)**:旧版(现 `--mapping identical`)只保留"相同残基对"——远缘同源物(如 RnRBCS1A vs 菠菜 SSU,全长一致 ~55%)覆盖率只有 ~68% 且 id=100%,**TM-score/lDDT 是在相同残基子集上算的,低 TM 未必等于折叠错误,可能只是映射稀疏**;遇到低覆盖先看 coverage/mapping_mode 字段,必要时用 homology 重跑或改链选择
- **TM-score 交叉验证**:本机实现忠实复刻官方 TMalign 20240303(get_initial 片段 20/100 + 跳步 15、全局 NW 仿射 gap -0.6/0、TMscore8 多种子迭代、d0_search 4.5-8、score_d8=1.5·L^0.3+3.5);WSL 内官方二进制对照:1brs 三组测试对误差 ≤0.02(典型 <0.01)
- 差模自检:AF2(single_sequence)对 1brs 晶体给出 TM=0.48/0.32、lDDT=0.36、DockQ=0.02 → grade 正确判"low confidence",与 pLDDT=36 一致
- 大蛋白(>200 残基)单链评估约 0.5-1 分钟;`--out` JSON 供下游使用

### 用法速查

- 单体:`--model m.pdb --ref r.pdb`(默认各自第一条链)
- 复合物:`--complex --model-chains A B --ref-chains A D [--rec-ref A --lig-ref D --rec-model A --lig-model B]`
- 映射:`--mapping auto|homology|identical`(远缘同源物解读见"关键行为")
- 自比校验(应 TM=1、lDDT=1):`--model x.pdb --ref x.pdb --complex`

## 二、prodigy_affinity.py — 结合亲和力(Bonvin 实验室 PRODIGY)

```powershell
# 蛋白-蛋白(ΔG + Kd)
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\prodigy_affinity.py" `
    --complex complex.pdb --chains A B --temperature 25 --out affinity.json
# 蛋白-小分子(配合 vina_dock / virtual_screen.py 使用)
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\prodigy_affinity.py" `
    --complex receptor.pdb --ligand pose.sdf --out lig_affinity.json
```

- 输出:ΔG(kcal/mol)、Kd(M)、温度、IC 接触分类、raw 文本
- 实测:1brs 晶体复合物 ΔG=-11.3 kcal/mol、Kd=5e-9 M(barnase-barstar 实测 ~10⁻¹³ M,IC 模型偏保守;用于**排序对比**而非绝对 ΔG)
- 依赖:venv 内 prodigy-prot 2.4.0 / prodigy-lig 1.1.4(PyPI 名**不是** `prodigy`,勿装错;freesasa 已随装)

## 三、典型组合拳(质量闭环)

1. AF2/ESMFold 预测 → 2. `struct_eval` 对照晶体/参考评估(不合格直接判废)→ 3. `pp_interact` 界面分析 → 4. `prodigy_affinity` 结合强度 → 5. 写进结果报告

## 四、常见错误

| 报错/现象 | 原因与解法 |
|---|---|
| `KeyError: 'A'`(链不存在) | 先 `--ref-chains/--model-chains` 指定正确链;用 `pdb_fetch.py meta` 或 Biopython 列链 |
| TM 与文献差很多 | 确认 `--model` 在前(移动)、`--ref` 在后(固定);编号不同没关系(自动映射) |
| prodigy `UnicodeEncodeError: gbk` | 脚本已强制子进程 UTF-8;勿改系统控制台编码 |
| `No module named prodigy` | PyPI 名为 prodigy-prot/prodigy-lig;`pip install prodigy-prot prodigy-lig` |
| DockQ=None | iRMS/LRMS 无有效对(受体或配体对 <3);检查链映射是否正确 |
| 远缘同源物 TM 很低但 coverage 也低(且 id≈100%) | 旧式"相同残基对"映射稀疏(identical 模式):低 TM 未必是折叠错误;用 `--mapping homology` 重跑,看全长一致率与覆盖 |
