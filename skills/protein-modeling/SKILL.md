---
name: protein-modeling
description: "本机科研级蛋白质结构预测与互作技能:LocalColabFold(AF2/AF2-Multimer,WSL2 GPU)预测、ESMFold 通道、ESM-2 嵌入特征、AutoDock Vina 对接、LightDock 蛋白-蛋白互作、Biopython 界面分析(pp_interact.py)、OpenMM MM-GBSA 结合自由能与显式溶剂 MD 协议(md_mmgbsa.py);部署契约见 D:\\bioai\\deploy-plan.md"
---

# 蛋白质结构预测与互作技能(AF2 · ESMFold · ESM-2 · Vina · LightDock · OpenMM)

本技能让 Agent 熟练调用本机部署的科研级结构预测与互作工具链(全部位于 **D 盘** `D:\bioai\`,C 盘零大文件;部署契约见 `D:\bioai\deploy-plan.md`)。

## 〇、环境与布局(速查)

- **GPU**:NVIDIA RTX 4060 Laptop **8GB**;AF2/AF2-Multimer 可行(reduced 设置),**AF3/Boltz-2 因显存硬性排除**
- **解释器铁律**:生信脚本用 `C:\Program Files\Python313\python.exe`(默认 python 是 3.14,勿用);`PYTHONPATH=D:\biopython` 提供 Biopython 1.87
- **布局**:

| 路径 | 内容 |
|---|---|
| `D:\bioai\miniforge3` | Windows Miniforge;colabfold 环境(WSL 内亦有一套) |
| `D:\bioai\venv` | Python 3.13 venv:vina + meeko + rdkit + openmm + numpy/scipy/biopython/pandas/matplotlib/seaborn/logomaker/pycirclize/mdtraj + prodigy-prot/prodigy-lig(化学信息学、可视化、MM-GBSA/MD、亲和力共用) |
| `D:\bioai\venv-esm` | 独立 venv:torch(CPU) + fair-esm 2.0 + numpy(ESM-2 嵌入) |
| `D:\bioai\models` | AF2 参数(COLABFOLDDIR)、torch hub 缓存(TORCH_HOME,ESM-2 权重) |
| `D:\bioai\jobs` | 预测/对接/模拟任务输出目录 |
| `D:\bioai\pip-cache` | PIP_CACHE_DIR 重定向(所有 pip 必须带) |
| `D:\bioai\preset-maintenance` | 技能备份与恢复(restore-bioinfo-skills.ps1) |
| WSL2(Ubuntu) | colabfold 1.5.5 + jax 0.4.22 GPU 栈 + TMalign(20240303,TM-score 交叉验证用) |

- **环境变量纪律**:`$env:PIP_CACHE_DIR='D:\bioai\pip-cache'`;`$env:TORCH_HOME='D:\bioai\models\torch-hub'`(ESM 权重缓存)

## 一、结构预测:LocalColabFold(AF2 / AF2-Multimer)—— 主力

**研究级定位**:AF2 是同源建模与复合物预测主力;AF2-Multimer 直接预测蛋白-蛋白复合物结构(互作预测首选)。**已完整部署在 WSL2 + GPU,验收通过**(1brs 复合物 28.9s/模型)。

```powershell
# 唯一入口(内置完整 WSL 环境:HOME=/root、cd /root、unset PYTHONPATH、COLABFOLDDIR、
# XDG_CACHE_HOME、XLA_PYTHON_CLIENT_PREALLOCATE=false、--disable-unified-memory)
& 'C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\protein-modeling\resources\run_colabfold.ps1' `
    -Fasta query.fasta -OutDir D:\bioai\jobs\af2_job1 `
    -ModelType alphafold2_multimer_v3 -NumModels 1 -NumRecycle 3 -MsaMode mmseqs2_uniref_env
```

- **复合物预测(互作,首选)**:fasta 中多条序列用 `:` 连接:
  ```
  >complex_A_B
  MSEQUENCEOFCHAINAAA...:MSEQUENCEOFCHAINBBB...
  ```
- `-ModelType` 取值(**1.5.5 实测**,`multimer` 是无效值):复合物 `alphafold2_multimer_v3`(默认)/ v1/v2;单体 `alphafold2_ptm` / `alphafold2`;或 `auto`
- `-MsaMode`:`mmseqs2_uniref_env`(MMseqs2 服务器,高质量,需外网)/ `single_sequence`(离线快速,精度低,验收时用此模式故 pLDDT≈36)
- **8GB 显存纪律**:`-NumModels 1`、`-NumRecycle 3`;OOM 加 `--disable-unified-memory`(脚本已内置)+ `XLA_PYTHON_CLIENT_PREALLOCATE=false`
- 结果:`*_unrelaxed_rank_001_*.pdb` + pLDDT/PAE PNG + scores JSON;质量评估用 `protein-quality` 技能的 `struct_eval.py`(TM-score/lDDT/DockQ,与官方 TMalign 交叉验证)

## 二、ESMFold 快速通道(云 API,零本地依赖)

- 客户端:`C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\protein-modeling\resources\esmfold_api.py`(纯标准库)
- 用法:
  ```powershell
  & 'C:\Program Files\Python313\python.exe' "...\resources\esmfold_api.py" --fasta query.fasta --out D:\bioai\jobs\esm_out.pdb
  ```
- 编程调用:`from esmfold_api import fold_sequence; pdb_text = fold_sequence(seq)`
- 特点:免费无 key、分钟级、无 MSA;适合快速粗模/断网兜底之外的"快速通道";速率受限时退避重试已内置

## 三、ESM-2 嵌入特征(残基级 / 序列级)

蛋白质序列 → 数值向量(下游 ML/聚类/突变效应/特征工程)。

```powershell
$env:TORCH_HOME='D:\bioai\models\torch-hub'
& 'D:\bioai\venv-esm\Scripts\python.exe' "...\resources\esm_embed.py" `
    --fasta query.fasta --model esm2_t6_8M_UR50D --outdir D:\bioai\jobs\esm_embed1
```

- 输出:`<id>_residue_L<层>.npy`(L×D 残基级)+ `<id>_sequence_L<层>.npy`(D 维 mean-pool)+ metadata.json
- 模型:t6_8M(D=320,~30MB,默认)/ t12_35M / t30_150M / **t33_650M(D=1280,科研级,2.5GB 下载需加速器)**
- 复合物 fasta(`:` 连链)自动按链拆分嵌入,ID 加 `_chainN` 后缀
- 权重缓存于 TORCH_HOME;下载中断用 `curl -L -C -` 断点续传(实测 dl.fbaipublicfiles.com 需加速器)

## 四、蛋白-蛋白互作

1. **复合物预测**:AF2-Multimer(第一节)——互作界面结构的最强手段
2. **对接(LightDock,WSL 内)**:
   ```powershell
   wsl.exe -d Ubuntu -- bash -lc "source /root/miniforge3/etc/profile.d/conda.sh && conda activate <env> && lightdock3_setup.py receptor.pdb ligand.pdb"
   ```
3. **HADDOCK3(云端,领域标杆)**:REST API 需注册账号/积分;适合 LightDock 不满足精度时
4. **界面分析(本技能分发脚本,必用)**:`resources/pp_interact.py`
   ```powershell
   & 'C:\Program Files\Python313\python.exe' "...\resources\pp_interact.py" --complex predicted_complex.pdb --chains A B --out interface.json
   ```
   输出:接触对(NeighborSearch 5Å)、界面残基、**埋藏表面积 BSA**(ShrakeRupley ΔSASA)、复合物 SASA;JSON + 人类可读报告。验收基线:1brs 晶体 55 接触 / BSA 1280.9 Å²
5. **质量与亲和力闭环**(protein-quality 技能):`struct_eval.py --complex`(TM-score/Fnat/iRMS/LRMS/DockQ 对照晶体)+ `prodigy_affinity.py`(ΔG/Kd)

## 五、蛋白-小分子对接:AutoDock Vina + meeko

```powershell
# 1) 配体准备(SMILES/SDF → PDBQT)
& 'D:\bioai\venv\Scripts\python.exe' -m meeko --input ligand.sdf --output ligand.pdbqt
# 2) 受体准备:推荐用 vina_dock.py 自动完成(PDB → 刚性 PDBQT,元素映射+0.00 电荷,
#    未知元素如 Se 自动回退 C;pdb_to_pdbqt.py 可单独调用)
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\pdb_to_pdbqt.py" receptor.pdb -o receptor.pdbqt
# 3) 对接(vina 官方 1.2.7 二进制;pip 无 cp313 轮子,勿用 pip 装)
& 'D:\bioai\bin\vina_1.2.7_win.exe' --receptor receptor.pdbqt --ligand ligand.pdbqt --center_x 10 --center_y 10 --center_z 10 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 16 --out docked.pdbqt
```
- 便捷包装(推荐):`resources/vina_dock.py` — SMILES/SDF → meeko 配体准备 → 对接 → 打分排序 → 姿势 PDB + JSON 报告:
```powershell
& 'D:\bioai\venv\Scripts\python.exe' "...\resources\vina_dock.py" --receptor-pdbqt receptor.pdbqt --smiles "CCO" --center 10 10 10 --size 20 20 20 --outdir D:\bioai\jobs\dock1
```
- **批量筛选**:`chem-informatics` 技能的 `virtual_screen.py`(化合物库 → 排序 CSV,断点续跑;3PTB 苯甲脒阳性对照排名第一)
- 8GB 卡跑 Vina 用 CPU 即可;筛选大库用 `--exhaustiveness 8`

## 六、MM-GBSA 结合自由能与 MD 协议(OpenMM)

```powershell
$M = "C:\deepseek-harness\.dsh\.agent-presets\bioinfo\skills\protein-modeling\resources\md_mmgbsa.py"
# 1) MM-GBSA(OBC2 隐式溶剂,分钟级;蛋白-蛋白)
& 'D:\bioai\venv\Scripts\python.exe' $M --mode gb --complex complex.pdb --rec-chains A --lig-chains D --out mmgbsa.json
# 2) 显式溶剂 MD(加热 → 平衡 → 采样 → RMSD/RMSF 分析)
& 'D:\bioai\venv\Scripts\python.exe' $M --mode md --complex complex.pdb --steps 100000 `
    --heating-steps 10000 --equil-steps 50000 --outdir D:\bioai\jobs\md1 [--platform OpenCL]
```

- gb 输出:ΔG_bind = E(complex) − E(receptor) − E(ligand),分解 internal(键/角/二面)与非键(vdW+GB+SA);**绝对值过稳定(1brs 实测 -36.7 vs 实验 -18.9 kcal/mol),用于排序/对比**
- md 输出:trajectory.dcd + state.csv(温度/势能/体积)+ rmsd.csv/rmsf.csv + rmsd_rmsf.png + system.prmtop
- 内置 PDB 清洗(全部实测配方):MSE→MET、去 HETATM/结晶水、**丢弃原子集不完整的残基**、断口自动 TER 分链 + OXT 封端(晶体编号断档/缺失侧链是 OpenMM 模板匹配最常见的坑)
- 平台:默认 CPU(确定性好);`--platform OpenCL` 提速(显式溶剂大体系建议)
- 时间尺度:采样 100k 步 = 200 ps;真实项目建议 ≥50 ns(2500 万步,需数小时,后台跑)

## 七、常见错误速查

| 报错/现象 | 原因与解法 |
|---|---|
| `CUDA out of memory` | 8GB 显存:减 `--num-recycle`/模型数,或 `--msa-mode single_sequence` |
| colabfold 参数下载到 C 盘 | 未设重定向变量;按第〇节设置后重跑(部署时已实测变量名) |
| `RemoteDisconnected`(ESM Atlas) | 已内置浏览器 UA + 退避重试;仍失败则稍后再试或切 AF2 |
| vina 报受体/配体 PDBQT 不合法 | 用 meeko 重做配体;受体去水/去配体后再 rigid 转换 |
| LightDock Windows 不可用 | 走 WSL 或改 AF2-Multimer 复合物预测 + pp_interact 分析 |
| `ModuleNotFoundError: Bio` | 生信脚本必须用 Python313 + `PYTHONPATH=D:\biopython`,不要用 venv 里的 python 跑 Bio 脚本 |
| OpenMM `No template found for residue (MSE)` | 力场无 Se 模板:脚本已自动 `MSE→MET`、` SE → SD ` |
| OpenMM `No template found for residue (HOH)` | 结晶水无模板:脚本已自动删除 |
| OpenMM `missing N H atoms` | PDB 无氢:`Modeller.addHydrohens`(脚本已内置) |
| OpenMM `No template ... matches CXXX, external C too many/few` | 编号断档/缺原子/中链 OXT:md_mmgbsa.py 的清洗层已处理(完整原子集过滤 + TER 分链 + OXT 封端),勿绕过 load_pdb |
| vina `Atom type SE is not a valid AutoDock type` | 受体含硒代甲硫氨酸;pdb_to_pdbqt.py 已自动回退 C,勿用其他转换器 |
| colabfold `CUDA_ERROR_OUT_OF_MEMORY`(34GB unified alloc) | 加 `--disable-unified-memory` + `XLA_PYTHON_CLIENT_PREALLOCATE=false`(run_colabfold.ps1 已内置) |
| colabfold `module 'jax.random' has no attribute 'default_prng_impl'` | jax/haiku 版本错配;**本机已验证组合:jax==0.4.22 + jaxlib==0.4.22 + jax-cuda12-plugin/pjrt==0.4.22 + dm-haiku==0.0.10 + pandas<2**(0.4.25 会触发 RNG 不兼容) |
| colabfold `module 'jax' has no attribute 'linear_util'` | haiku 0.0.10 补丁:dot.py 顶部加 `from jax.extend import linear_util`,并把 `jax.linear_util` 全部替换为 `linear_util`(本机已应用) |
| colabfold 启动即挂(开 socket 写 params_model_1.npz) | 参数必须放**两个目录**:`COLABFOLDDIR/params/` 与 `XDG_CACHE_HOME/colabfold/params/`,各含全部 npz + 对应 `download_*_finished.txt` 标记 |
| colabfold `argument --model-type: invalid choice: 'multimer'` | 1.5.5 用 `alphafold2_multimer_v3`(run_colabfold.ps1 已内置);MSA 用 `--msa-mode mmseqs2_uniref_env` 或 `single_sequence` |
| ESM-2 权重下载中断 | `curl -L -C -` 断点续传至 TORCH_HOME/hub/checkpoints(需加速器) |
| esm_embed `KeyError: ':'` | 复合物 fasta 的 `:` 连链格式;脚本已自动按链拆分 |
| GCS 直连时好时坏/0 字节 | 多前端 IP 部分被墙;加速器开启时下载加 `-Proxy http://127.0.0.1:7897`(实测 2-15MB/s);或 `D:\bioai\bin\parallel-download.ps1` 钉健康 IP + 持久分块续传 |
| pip 下载占满 C 盘 | 忘设 `PIP_CACHE_DIR`;pip 缓存必须重定向(第〇节) |

## 八、什么时候用这个 skill

- 序列 → 三维结构(单链或复合物)预测与质量评估(pLDDT/PAE + TM-score/lDDT/DockQ)
- 序列 → 嵌入特征(ESM-2,残基级/序列级)
- 蛋白-蛋白互作:复合物预测、对接、界面残基/接触/埋藏面积分析、结合亲和力(ΔG/Kd)
- 蛋白-小分子对接与虚拟筛选、姿势打分排序
- 结合自由能(MM-GBSA)与显式溶剂 MD(加热/平衡/采样 + RMSD/RMSF)
