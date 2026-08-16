# 从零部署:生信模式(dsh-bioinfo)完整复刻指南

目标:在**你自己的 DSH 上**,复刻出与本项目开发机**行为一致**的「生信模式」预设 agent。

## 0. 复刻基准(重要)

以下为**基准布局**——按本指南原样部署即可实现零改动、逐文件一致的复刻:

| 项 | 基准值 |
|---|---|
| DSH 预设根 | `C:\deepseek-harness\.dsh\.agent-presets\bioinfo\`(DSH_HOME=`C:\deepseek-harness\.dsh`) |
| Python 3.13 | `C:\Program Files\Python313\python.exe` |
| Biopython | `D:\biopython`(1.87 源码) |
| 工具栈 | `D:\bioai\(venv、venv-esm、bin、models、jobs、pip-cache、wsl)` |

偏离基准时的处理:

- **插件侧**:用 `BIO_TOOLS_*` 环境变量覆盖(见 §5),无需改代码。
- **技能文档/脚本内路径**:SKILL.md 示例与部分脚本默认值引用上述路径,非基准布局时全局替换预设目录字符串:
  ```powershell
  $old='C:\deepseek-harness\.dsh\.agent-presets\bioinfo'
  $new='<你的预设目录>'
  Get-ChildItem $new -Recurse -Include *.md,*.js,*.ps1 |
    ForEach-Object { (Get-Content $_.FullName -Raw) -replace [regex]::Escape($old), $new |
      Set-Content $_.FullName -NoNewline }
  ```

前置条件:Windows 10/11、管理员权限(WSL 启用用)、可选 NVIDIA GPU(仅 `af2_predict` 推荐,无 GPU 亦可跑但慢,可用 `esmfold_predict` 云端通道替代)、磁盘 ≥ 30 GB(AF2 参数 7.6 GB 在 D 盘)。

## 1. 安装 DSH

```powershell
npm install -g @deepseek-ai/dsh
```

具体启动/配置以 DSH 官方文档为准。

## 2. 挂载预设

把仓库中以下内容复制到 `<DSH_HOME>\.agent-presets\bioinfo\`:

```
preset.yml
agent.cordis.yml
plugins\
skills\
```

然后启动 DSH,新会话选择「生信模式」。

> **运维纪律(必读)**:预设组合(`agent.cordis.yml`、`plugins/`、`skills/`)在
> **DSH 宿主进程启动时只挂载一次**,之后修改任何预设文件都必须**重启宿主进程**
> (仅新开会话不会重载,内存中仍是旧版本)。典型症状:生信模式会话开场即报
> `Invalid schema for function 'af2_predict' ... got 'type: null'`——磁盘已是
> 修复版但进程未重启。

## 3. Python 3.13 + Biopython 1.87

```powershell
# 1) 安装 Python 3.13 到 C:\Program Files\Python313(官方安装包)
# 2) 取 Biopython 1.87 源码(基准布局):
git clone --branch biopython-187 https://github.com/biopython/biopython.git D:\biopython
# 3) 依赖(与开发机一致):
& 'C:\Program Files\Python313\python.exe' -m pip install numpy==2.4.6
# 4) 验证:
$env:PYTHONPATH='D:\biopython'
& 'C:\Program Files\Python313\python.exe' -c "import Bio; print(Bio.__version__)"   # -> 1.87
```

## 4. D:\bioai 工具栈

### 4.1 目录与 pip 缓存

```powershell
New-Item -ItemType Directory -Force D:\bioai\bin, D:\bioai\jobs, D:\bioai\models, D:\bioai\pip-cache, D:\bioai\wsl
$env:PIP_CACHE_DIR='D:\bioai\pip-cache'   # 所有 pip 安装带上(防 C 盘膨胀)
```

中国大陆网络:配置 pip TUNA 镜像(`pip.ini` index-url = `https://pypi.tuna.tsinghua.edu.cn/`)与 conda 镜像(见 `deploy\wsl-bootstrap.sh` 内置 `.condarc`);海外用户跳过。

### 4.2 主 venv(结构评估/对接/筛选/亲和力/MD/可视化)

开发机精确版本(`pip freeze` 实录):

```powershell
& 'C:\Program Files\Python313\python.exe' -m venv D:\bioai\venv
& 'D:\bioai\venv\Scripts\python.exe' -m pip install `
  biopython==1.88 rdkit==2026.3.5 meeko==0.7.1 openmm==8.5.2 numpy==2.4.6 `
  scipy==1.18.0 pandas==3.0.5 matplotlib==3.11.1 seaborn==0.13.2 logomaker==0.8.7 `
  pycirclize==1.10.1 parmed==4.3.1 mdtraj==1.11.1.post2 gemmi==0.7.5 freesasa==2.2.1 `
  prodigy-prot==2.4.0 prodigy-lig==1.1.4
```

注意:PyPI 包名是 `prodigy-prot` / `prodigy-lig`(不是 `prodigy`)。若镜像缺 `rdkit==2026.3.5`,可用最近的 2025.x 版本,验收数值以本清单版本为准。

### 4.3 venv-esm(ESM-2 嵌入)

```powershell
& 'C:\Program Files\Python313\python.exe' -m venv D:\bioai\venv-esm
& 'D:\bioai\venv-esm\Scripts\python.exe' -m pip install `
  torch==2.9.1 --index-url https://download.pytorch.org/whl/cpu
& 'D:\bioai\venv-esm\Scripts\python.exe' -m pip install fair-esm==2.0.0 numpy==2.5.2
```

### 4.4 AutoDock Vina 二进制

从官方 release 下载 Windows 版 1.2.7 → `D:\bioai\bin\vina_1.2.7_win.exe`(pip 无 cp313 轮子,**不要** pip 装)。

### 4.5 WSL2 + LocalColabFold(`af2_predict` 用)

```powershell
# 管理员 PowerShell(仅一次,随后重启;DISM 走本地组件库,无需联网):
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
# 重启后:
pwsh -File deploy\wsl-setup.ps1     # 自动 import Ubuntu 22.04 rootfs + Miniconda + colabfold + LightDock
```

- rootfs 可用 `deploy\parallel-download.ps1` 预下载(TUNA 镜像 URL 见 `wsl-setup.ps1` 报错提示)。
- `wsl-bootstrap.sh` 内置 TUNA/USTC 源,海外用户可改回官方源。
- **colabfold 1.5.5 实测版本组合**:`jax==0.4.22 + jaxlib==0.4.22 + dm-haiku==0.0.10 + pandas<2`(0.4.25 有 RNG 不兼容);haiku 需打 dot.py 补丁(`from jax.extend import linear_util`)。

### 4.6 AF2 参数(~7.6 GB)

```powershell
# 推荐(与开发机一致):WSL 内分块直连 GCS,落 D:\bioai\models\colabfold\params
wsl -d Ubuntu -- bash /mnt/d/bioai/bin/wsl-download-params.sh
# 或 Windows 侧走代理:
pwsh -File deploy\download-af2-params.ps1 -Proxy http://127.0.0.1:7897
```

参数必须同时放在 `COLABFOLDDIR/params/` **和** `XDG_CACHE_HOME/colabfold/params/` 两个目录,各含 npz 与对应 `download_*_finished.txt` 标记,否则 colabfold 启动即挂。

### 4.7 deploy 脚本落位

```powershell
Copy-Item deploy\*.ps1, deploy\*.sh D:\bioai\bin\ -Force
```

## 5. 环境变量覆盖(可选,偏离基准时用)

| 变量 | 默认值 | 含义 |
|---|---|---|
| `BIO_TOOLS_PYTHON` | `C:\Program Files\Python313\python.exe` | 带 Biopython 的 Python 3.13 |
| `BIO_TOOLS_VENV_PY` | `D:\bioai\venv\Scripts\python.exe` | 主 venv 解释器 |
| `BIO_TOOLS_RES_DIR` | `<预设>\skills\protein-modeling\resources` | protein-modeling 后端 |
| `BIO_TOOLS_RES_PQ_DIR` | `<预设>\skills\protein-quality\resources` | struct_eval 后端 |
| `BIO_TOOLS_RES_CI_DIR` | `<预设>\skills\chem-informatics\resources` | virtual_screen 后端 |
| `BIO_TOOLS_JOBS_DIR` | `D:\bioai\jobs` | 默认输出根 |
| `BIO_TOOLS_BIOPYTHON` | `D:\biopython` | Bio 导入 PYTHONPATH |

## 6. 验收:证明复刻一致

1. **自洽性**:`pwsh -File scripts\verify-layout.ps1` → `LAYOUT OK`
2. **插件 schema**:仓库内 `npm test` → `ALL SCHEMAS OK (7 tools)`
3. **分项阳性对照**(夹具在 `fixtures\acceptance\`,预期值详见 `fixtures/README.md`):
   - `struct_eval`(1brs 预测 vs 晶体):TM-score 对官方 TMalign 误差 ≤ 0.02
   - `prodigy_affinity`:1brs ΔG = **-11.3 kcal/mol**
   - `md_run --mode gb`:MM-GBSA dG_bind = **-36.7 kcal/mol**
   - `vscreen_run`(3PTB,7 化合物):苯甲脒排名 **#1(-5.90 kcal/mol)**
   - `esm_embed`:barnase/barstar 各链 110/90 × 320 嵌入
4. **端到端**:`pwsh -File deploy\run-acceptance.ps1`(AF2-Multimer 预测 + 界面分析;对照晶体基线 **55 contacts / BSA 1280.9 Å²**)

## 7. 常见坑(开发机实测记录)

- colabfold 参数目录/标记文件要求(§4.6);启动报 `argument --model-type: invalid choice` 时用 `alphafold2_multimer_v3`(`run_colabfold.ps1` 已内置)。
- GCS 在 Windows 侧被墙:走 WSL(镜像网络模式 + 禁 IPv6),见 `.wslconfig [wsl2] networkingMode=mirrored`。
- 显式溶剂 MD:`app.DCDFile` 只是写入器,轨迹读取/对齐用 **mdtraj**(parmed 4.3.1 读不了 DCD);分析必须用模拟时的同一拓扑。
- logomaker 0.8.7 用 `width=`(非 `stack_width`);pycirclize 1.10 的 `Track.text(text, x=, r=)`;pandas 3.0 标量 DataFrame 构造行为变化。
- prodigy 输出 GBK 编码错误:脚本已强制子进程 UTF-8,勿改控制台编码。
- Vina `Atom type SE is not a valid AutoDock type`:受体含硒代甲硫氨酸,`pdb_to_pdbqt.py` 已自动回退。

## 8. 合规提醒

受限组件(学术限用)见 `THIRD_PARTY_NOTICES.md`:TM-align/TMalign、PRODIGY、ESMFold 本地权重、KEGG。商用前须取得相应授权;仓库不分发这些组件本身。
