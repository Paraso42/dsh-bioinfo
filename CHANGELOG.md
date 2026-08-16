# 更新日志 (Changelog)

## v0.2.0 — 2026-08(用户反馈修复与行为变更版本)

**致全体 dsh-bioinfo(生信模式)用户:**

本版本基于一项真实科研项目(浮苔 Rubisco 全流程:序列获取 → 比对/树/热图 → 结构界面 → 互作网络 → AF2 预测 → 结构评估 → 引物设计)的完整使用反馈发布,共修复关键缺陷 3 项、改进评估语义 1 项、完善技能文档 4 份。全部变更已通过自动化校验与验收基线回归测试,并同步至 GitHub 仓库(提交 `c4fdb20`)。

### 一、修复内容

1. **`af2_predict` 原生工具整体不可用(WSL 发行版名解码缺陷,重大)。**
   现象:运行即报 `WSL_E_DISTRO_NOT_FOUND`。原因:`wsl -l -q` 输出被按 UTF-16 解码,发行版名内残留 NUL 字符。处理:`run_colabfold.ps1` 在解析前剥离 `\u0000`,并固化修复至仓库、现场 preset 与维护备份。实测验证:181 aa 单体 `single_sequence` 模式 25 秒、MSA 模式 42 秒。

2. **`pp_interact` 工具结果序列化失败(工具层)。**
   现象:工具返回 `value is not lossless JSON`,同一输入经 CLI 调用正常。原因:报告含 NaN/Infinity 等非有限浮点值。处理:双层清洗——脚本落盘启用 `allow_nan=False` 并清洗非有限值;插件读取层对全部 7 个工具的报告统一清洗(文本层 + 树层)。1brs 验收基线回归一致(55 接触 / 1280.9 Å² 逐位复现)。

3. **ESM Atlas 通道长期不可用。**
   现象:连续多轮请求全部返回 HTTP 504(3 轮共 12 次)。处理:失败信息内置本地 AF2 兜底指引(`run_colabfold.ps1 -MsaMode single_sequence`),`esmfold_predict` 工具失败时自动附提示;文档中将该通道降级标注为**备用通道**。

### 二、行为变更与兼容性说明(敬请阅读)

- **`struct_eval` 残基映射策略变更。** 旧版对序列不同的模型与参考仅映射"相同残基对"(覆盖率偏低、`seq_identity` 恒为 100%),TM-score/lDDT 在该子集上计算,对远缘同源物可能低估结构质量。新版默认采用**全长同源全局比对**(`--mapping auto|homology`),`coverage`/`seq_identity` 反映真实全长一致率;旧行为可通过 `--mapping identical` 显式启用。每链报告新增 `mapping_mode` 字段。
- **历史结果甄别方法:** 报告中无 `mapping_mode` 字段且 `seq_identity≈1.0`、`coverage<1.0` 者为旧版输出。涉及远缘同源对比且位于决策边界(如 TM-score 跨 0.5/0.7 阈值、模型判废、DockQ 分级)时,建议以新版重算(单链约 1 分钟)。
- **`pdb_fetch meta` 输出增强:** 新增 `chain_residues_polymer` 字段,提供逐链聚合物残基数(去结晶水、取首模型),大/小亚基可一眼分辨。
- **其余工具算法未做任何变更:** `pp_interact`、`vina_dock`、`vscreen_run`、`md_run` 的历史结果继续有效,无需重算。

### 三、质量验证

| 项目 | 结果 |
|---|---|
| `npm test`(7 工具模式校验) | 全部通过 |
| Python 语法校验(4 个修改脚本) | 全部通过 |
| pp_interact 验收基线回归(1brs 晶体) | 55 接触 / 1280.92 Å²,逐位一致 |
| 映射行为测试(合成 55% 一致序列对) | homology:覆盖 100%/一致率 59.3%;identical:覆盖 61%/一致率 100%(复现旧版症状,新默认已修复) |

### 四、用户行动项

1. **更新 preset 文件并重启 DSH 主机进程**:插件层变更须重启生效(仅新建会话不生效);CLI 脚本无需重启。
2. **同步刷新维护备份**:如保有 `preset-maintenance` 技能备份,请一并刷新,避免还原旧版时覆盖本次修复。
3. **记录所用版本**:建议在项目报告中标注 preset 版本(commit `c4fdb20`)与 `mapping_mode`,确保结果可追溯。

### 五、数据通道现状(实测)

| 通道 | 状态 |
|---|---|
| NCBI Entrez / BLAST(浏览器 UA)/ Datasets | 稳定 |
| UniProt / RCSB PDB / STRING | 稳定 |
| MMseqs2(colabfold 服务器,经 WSL) | 可用 |
| NCBI FTP | 带宽整形(~1 KB/s 级),请使用 Datasets API 或 efetch 分段 |
| ESM Atlas | 连续 504,已降级为备用通道 |

### 六、后续计划

- 蛋白家族 × 物种面板一键工作流(UniProt 优先 → BLAST-TSA 兜底 → 同源矩阵/热图/NJ 树);
- 克隆引物设计模块(酶切位点冲突检测 + 兼容酶回退 + Tm_NN 选长 + CSV 输出);
- UniProt→CDS 全自动核验脚本(EMBL 交叉引用 → GenBank 全长 mRNA → 转运肽检查)。

反馈与问题请通过 GitHub Issues 提交;安全漏洞请走私有通告通道(见 SECURITY.md)。

---

## v0.1.0 — 2026-08-15(初始开源版本)

完整预设套件首次发布:preset 身份与人格、`protein-tools` 插件(7 个模型工具)、7 技能库(14 个后端脚本)、部署脚本(WSL2 LocalColabFold + AF2 参数下载)、验收夹具与自证脚本,以及 MIT 许可与 THIRD_PARTY_NOTICES 合规文档(提交 `946444f` 至 `f6832e3`)。
