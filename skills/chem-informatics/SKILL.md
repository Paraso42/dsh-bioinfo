---
name: chem-informatics
description: "RDKit 化学信息学与 AutoDock Vina 虚拟筛选:分子性质/标准化/相似性/子结构/构象/2D 图(mol_tools.py)+ 化合物库批量对接管道(virtual_screen.py,断点续跑、阳性对照已验收);环境 D:\\bioai\\venv(rdkit/meeko)+ D:\\bioai\\bin\\vina"
---

# 化学信息学与虚拟筛选(chem-informatics)

小分子层面的完整工具链:分子描述与相似性(mol_tools.py)+ 批量对接(virtual_screen.py)。

## 〇、环境

- 解释器:`D:\bioai\venv\Scripts\python.exe`(rdkit 2026.3 + meeko 0.7.1)
- Vina 二进制:`D:\bioai\bin\vina_1.2.7_win.exe`(脚本自动探测)
- 脚本目录:`C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\chem-informatics\resources\`

## 一、mol_tools.py — RDKit 分子套件

```powershell
$M = "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\chem-informatics\resources\mol_tools.py"
# 1) 性质表(Lipinski/QED/TPSA/可旋转键/芳香环数)
& 'D:\bioai\venv\Scripts\python.exe' $M describe --library lib.csv --out props.csv
& 'D:\bioai\venv\Scripts\python.exe' $M describe --smiles "CC(=O)Oc1ccccc1C(=O)O" --json
# 2) 标准化(去盐/规范互变异构;药物筛选前必做)
& 'D:\bioai\venv\Scripts\python.exe' $M canonical "CC(=O)O.[Na+]"        # -> CC(=O)O
# 3) 2D 图(可高亮 SMARTS 子结构)
& 'D:\bioai\venv\Scripts\python.exe' $M depict --smiles "CC(=O)Oc1ccccc1C(=O)O" --highlight "c1ccccc1" --out mol.png
# 4) 相似性检索(Morgan2/3、AtomPair;按 Tanimoto 排序)
& 'D:\bioai\venv\Scripts\python.exe' $M similarity --query "NC(=N)c1ccccc1" --library lib.csv --top 20 --out sim.csv
# 5) 子结构匹配(计数 + 高亮图)
& 'D:\bioai\venv\Scripts\python.exe' $M substructure --smiles "..." --smarts "c1ccccc1" --png sub.png
# 6) 构象生成(ETKDGv3 + MMFF 优化,能量排序 SDF)
& 'D:\bioai\venv\Scripts\python.exe' $M conformers --smiles "..." --n 10 --out confs.sdf
```

- 库文件约定:CSV,列名含 `smiles`(可用 `--smiles-col` 覆盖);坏行跳过并在结果中记录 error
- 已实测:7 化合物性质表、相似性检索(苯甲脒自比 1.0)、子结构高亮、5 构象生成

## 二、virtual_screen.py — Vina 批量虚拟筛选管道

流程:受体 PDB 自动刚性 PDBQT(排除水与指定残基)→ 逐配体 meeko 准备 → Vina 对接 → 增量写 results.csv(断点续跑)→ 排序汇总 + 前 N 姿势转 PDB + 可选 PRODIGY-Lig 亲和力。

```powershell
$V = "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\chem-informatics\resources\virtual_screen.py"
# 盒子自动取自共晶配体(推荐);--exclude-res 排除配体/离子/水
& 'D:\bioai\venv\Scripts\python.exe' $V `
    --receptor 3ptb.pdb --ligands lib.csv `
    --ref-ligand 3ptb_ben.pdb --exclude-res HOH,WAT,BEN,CA,SO4 `
    --exhaustiveness 8 --num-modes 3 --top 5 `
    --outdir D:\bioai\jobs\vscreen1 --out report.json
# 或显式盒子:--center x y z --size x y z
# 附加结合亲和力:--prodigy(对前 N 姿势跑 prodigy-lig)
```

- 输出:`results.csv`(rank/id/smiles/best_score/dg/kd/error)、`top_NN_*.pdb` 姿势、JSON 报告
- **断点续跑**:默认开启;已完成的配体自动跳过,中断后重跑同一 `--outdir` 即可;`--no-resume` 重头来
- **阳性对照已验收**(3PTB 胰蛋白酶):苯甲脒(共晶配体)+ 6 个诱饵(苯/甲苯/苯酚/扑热息痛/布洛芬/阿司匹林),exhaustiveness 8 → **苯甲脒以 -5.90 kcal/mol 排名第一**,诱饵 -3.7~-5.7
- 大库建议:`--exhaustiveness 8 --num-modes 3`(初筛),命中后再 16/9 精对接
- 受体准备等价于 `pdb_to_pdbqt.py`(元素映射,未知元素回退 C;见 protein-modeling 技能);脚本内联实现以保持独立运行

## 三、与结构预测的衔接

1. 蛋白靶标:`pdb_fetch.py download` 取晶体 → `--exclude-res` 去掉共晶配体/离子
2. 无晶体口袋时:AF2 预测结构 + `virtual_screen`(注意侧链精度)
3. 命中验证:`vina_dock.py` 精对接 + `prodigy_affinity.py --ligand` 亲和力 + `pp_interact`-风格接触分析

## 四、常见错误

| 报错/现象 | 原因与解法 |
|---|---|
| `cannot parse SMILES` | 去盐/标准化后再试;库文件坏行会自动跳过 |
| meeko `produced no setup` | 分子过大/罕见元素;检查 SMILES 合法性 |
| vina `Atom type XX is not valid` | 受体含罕见元素;`--exclude-res` 排除该残基,或换 `pdb_to_pdbqt.py`(Se 已自动回退) |
| `vina binary not found` | 确认 `D:\bioai\bin\vina_1.2.7_win.exe` 存在(pip 无 cp313 轮子,勿用 pip 装) |
| 所有配体打分接近 | 盒子不对齐口袋;用 `--ref-ligand` 自动取盒,或扩大 `--size` |
| `--prodigy` 报错 | prodigy-lig 未装:`pip install prodigy-lig`(PyPI 名非 prodigy) |
