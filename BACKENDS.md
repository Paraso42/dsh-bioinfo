# Backend reference

The `protein-tools` plugin is a frontend: every tool shells out to a backend
script shipped under `skills/`. This file is the authoritative map of what
each tool runs and what it consumes.

## Layout keys (env-overridable in the plugin)

| Variable | Default | Meaning |
|---|---|---|
| `BIO_TOOLS_PYTHON` | `C:\Program Files\Python313\python.exe` | Python 3.13 with Biopython |
| `BIO_TOOLS_VENV_PY` | `D:\bioai\venv\Scripts\python.exe` | venv: RDKit/meeko/OpenMM/mdtraj/PyMOL 3.1.0/… |
| `BIO_TOOLS_RES_DIR` | `<preset>\skills\protein-modeling\resources` | protein-modeling backends |
| `BIO_TOOLS_RES_PQ_DIR` | `<preset>\skills\protein-quality\resources` | struct_eval / prodigy backends |
| `BIO_TOOLS_RES_CI_DIR` | `<preset>\skills\chem-informatics\resources` | virtual_screen / mol_tools backends |
| `BIO_TOOLS_JOBS_DIR` | `D:\bioai\jobs` | default output root |
| `BIO_TOOLS_BIOPYTHON` | `D:\biopython` | `PYTHONPATH` for Bio imports |

## Tool → backend table

| Tool | Backend script | Interpreter | Key third-party deps |
|---|---|---|---|
| `esmfold_predict` | `<RES>\esmfold_api.py` | `BIO_TOOLS_PYTHON` (stdlib only) | ESM Atlas API (network) |
| `pp_interact` | `<RES>\pp_interact.py` | `BIO_TOOLS_PYTHON` | Biopython 1.87 |
| `vina_dock` | `<RES>\vina_dock.py` (+ `pdb_to_pdbqt.py`) | `BIO_TOOLS_VENV_PY` | RDKit, meeko, AutoDock Vina binary |
| `af2_predict` | `<RES>\run_colabfold.ps1` | pwsh → WSL2 colabfold | LocalColabFold 1.5.5, jax 0.4.22, MMseqs2 server |
| `struct_eval` | `<RES_PQ>\struct_eval.py` | `BIO_TOOLS_VENV_PY` | numpy/scipy/Biopython; TM-score reimplementation (academic, see NOTICES) |
| `vscreen_run` | `<RES_CI>\virtual_screen.py` | `BIO_TOOLS_VENV_PY` | RDKit, meeko, Vina binary; optional prodigy-lig |
| `md_run` | `<RES>\md_mmgbsa.py` | `BIO_TOOLS_VENV_PY` | OpenMM, mdtraj, parmed |

## argv / output contract (what the plugin passes, what it reads back)

### esmfold_predict
- argv: `<sequence> --out <pdb> [--retries N]`
- back: exit 0 + PDB written to `--out`
- result: `{exitCode, stdout, stderr, pdbPath}`
- **availability**: the ESM Atlas endpoint has been repeatedly down (HTTP 504,
  2026-08). Treat as a fallback channel; on failure the tool appends a hint to
  use `af2_predict` with `msaMode: "single_sequence"` (offline AF2).

### pp_interact
- argv: `--complex <pdb> --chains <A> <B> [--cutoff F] --out <json>`
- back: JSON report (contacts / interface residues / BSA)
- result: `{exitCode, stdout, stderr, jsonPath, report}`

### vina_dock
- argv: `(--receptor-pdbqt <p> | --receptor <pdb>) (--smiles <s> | --ligand <f> | --ligand-pdbqt <f>) --center x y z --size x y z [--exhaustiveness N] --outdir <dir> --name dock --out <report.json>`
- back: poses + `<outdir>\report.json`
- result: `{exitCode, stdout, stderr, outdir, report}`

### af2_predict
- argv: `-Fasta <f> -OutDir <d> -ModelType <t> [-NumModels N] [-NumRecycle N] [-MsaMode m]`
- back: colabfold layout (ranked PDBs, PAE/pLDDT PNGs, scores JSON)
- result: `{exitCode, stdout, stderr, outdir}`

### struct_eval
- argv: `--model <pdb> --ref <pdb> [--complex] [--model-chains .. --ref-chains ..] [--rec-ref/--lig-ref/--rec-model/--lig-model ..] [--mapping auto|homology|identical] --out <json>`
- back: JSON (TM-score, CA/all-atom RMSD, lDDT, GDT-TS/HA, DockQ, grade); per-chain `mapping_mode`, `coverage`, `seq_identity`
- mapping: `auto` (direct when sequences identical, else full homology alignment) / `homology` / `identical` (legacy: identical-residue pairs only — sparse coverage on distant homologs; metrics are computed on the mapped subset)
- result: `{exitCode, stdout, stderr, jsonPath, report}`

### vscreen_run
- argv: `--receptor <pdb> --ligands <csv> (--ref-ligand <pdb> | --center x y z --size x y z) [--exclude-res "HOH,WAT"] [--exhaustiveness N] [--top N] --outdir <dir> --out <report.json>`
- back: resumable `results.csv`, ranked `report.json`, top poses as PDB
- result: `{exitCode, stdout, stderr, outdir, report}`

### md_run
- gb argv: `--mode gb --complex <pdb> [--rec-chains A] [--lig-chains B] [--platform P] --out <json>` → MM-GBSA dG report; result `{exitCode, stdout, stderr, jsonPath, report}`
- md argv: `--mode md --complex <pdb> [--steps N] [--platform P] --outdir <dir>` → `md_report.json`, DCD trajectory, RMSD/RMSF plots; result `{exitCode, stdout, stderr, outdir, report}`

## Common result envelope

`{exitCode, stdout, stderr}` plus the tool-specific fields above; when the
backend writes a JSON report the plugin embeds it verbatim as `report`.
Backend JSON is sanitized on read: non-finite floats (`NaN`/`Infinity`, e.g.
empty interfaces in `pp_interact`) become `null` so tool results stay lossless
JSON. Missing backends surface the shell error to the model — no crash.

## Data channel reliability (field-tested 2026-08)

| Channel | Status |
|---|---|
| NCBI Entrez / BLAST (browser UA) / Datasets | stable |
| UniProt / RCSB PDB / STRING | stable |
| MMseqs2 (colabfold server, via WSL) | **unreachable 2026-08** (official API blocked, backup site abandoned) — see local MSA DB installer below |
| GWDG colabfold DB mirror (Göttingen) | reachable at file level (2026-08) — source for `scripts/install-local-msa.ps1` |
| NCBI FTP | bandwidth-shaped (~1 KB/s); use Datasets API / efetch slices |
| ESM Atlas (`esmfold_predict`) | **frequently down (repeated 504)** — fallback channel; local AF2 `single_sequence` is the offline fallback |

## Local MSA database (offline MSA search)

`scripts/install-local-msa.ps1` installs mmseqs2 (apt → conda via TUNA
conda-forge mirror) plus UniRef30 + colabfold_envdb (~70 GB, resume-safe,
`-Proxy` supported). After installation the two-step local recipe is
`colabfold_search` (against `D:\bioai\msa-db`) → `colabfold_batch`. Auto-wiring
a local-MSA mode into `run_colabfold.ps1`/`af2_predict` is pending the DB
download + end-to-end smoke test.

## Visualization

- `pymol_render.py` (`<RES>\pymol_render.py`, `BIO_TOOLS_VENV_PY`): headless
  PyMOL ray-trace renders (publication/cartoon/rainbow/surface/line styles,
  hetatm coloring, semi-transparent surface overlay). PyMOL 3.1.0 open-source
  installed from cgohlke's cp313 wheel — the official pymol.org Windows
  bundles are discontinued and the PyPI Windows wheel is broken (both verified
  2026-08).
- matplotlib figures (`stat_plots.py`, `seq_logo.py`, `md_mmgbsa.py`): CJK
  font auto-selection (Microsoft YaHei → SimHei → Noto Sans CJK SC → SimSun →
  DejaVu fallback with warning).
