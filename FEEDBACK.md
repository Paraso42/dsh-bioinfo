# Feedback log

Real-world usage feedback and its disposition. Each entry records the round
(project, date), the findings, and what happened to them: **fixed** (code
change in this repo), **documented** (skill docs), or **backlog** (accepted,
not yet implemented).

## Round 1 — Riccia fluitans Rubisco project (2026-08)

Full pipeline exercised: sequence acquisition → alignment/tree/heatmap →
structure & interface → interaction network → AF2 prediction → structure
evaluation → primer design. Overall verdict: docs 95 / scripts 90 / tool
wrappers 75 — the three fixes below close the main gap.

### Fixed

1. **`af2_predict` unusable — WSL distro-name decode bug (critical).**
   `wsl -l -q` output can be re-encoded as UTF-16, leaking NUL bytes into the
   distro name (`U\0b\0u\0n\0t\0u\0`), so `wsl -d` failed with
   `WSL_E_DISTRO_NOT_FOUND`. Fixed in
   `skills/protein-modeling/resources/run_colabfold.ps1` by stripping `\u0000`
   before splitting (verified: 181-aa monomer, 25 s single_sequence / 42 s MSA).
   The same script backs the native `af2_predict` tool. Also recorded in the
   protein-modeling skill error table; note that restoring from an **old**
   `preset-maintenance` skills backup would reintroduce the bug.

2. **`pp_interact` tool result "value is not lossless JSON".**
   Backend reports can contain non-finite floats (`NaN`/`Infinity`, e.g. empty
   interfaces). Two-layer fix: `pp_interact.py` now dumps with
   `allow_nan=False` after sanitizing, and the plugin's JSON reader sanitizes
   NaN/Infinity tokens (text level) plus any non-finite numbers (tree level)
   for **all** tool reports. This protects `vina_dock`/`struct_eval`/`md_run`
   report reads too.

3. **ESM Atlas channel repeatedly 504 (12/12 attempts across 3 rounds).**
   The retry/backoff already worked as designed; the service itself was down.
   `esmfold_api.py` errors now include the local fallback command
   (`run_colabfold.ps1 -MsaMode single_sequence`), the `esmfold_predict` tool
   appends the same hint on failure, and the docs (protein-modeling skill,
   BACKENDS.md, README) downgrade ESMFold to a fallback channel.

### Improved in code

4. **`struct_eval` residue-mapping semantics (important caveat).** The old
   aligner weights (`match=1/mismatch=0/gap=0`) maximized identical-residue
   pairs only, so a distant homolog (RnRBCS1A vs spinach SSU, ~55% full-length
   identity) mapped just 84/123 residues with `seq_identity=100%` — TM-score
   was computed on that identical subset, which can understate quality for
   reasons unrelated to the fold. `map_residues` now uses standard homology
   weights by default (`--mapping auto|homology`): full-length coverage and
   honest identity. The legacy behavior remains available as
   `--mapping identical`; per-chain reports carry `mapping_mode`,
   `coverage`, `seq_identity`. Documented in the protein-quality skill and
   tool description.

5. **`pdb_fetch.py meta` now reports `chain_residues_polymer`** — per-chain
   polymer-only residue counts (ATOM records, first model, waters excluded),
   so large/small subunits are distinguishable at a glance (8RUC: 783 vs 207).

### Documented (skill docs)

- **biopython-analyses**: primer-design subsection — `Bio.SeqUtils.MeltingTemp.Tm_NN`
  usage + `Bio.Restriction` internal cut-site clash detection + annealing-region
  length selection; quick-check entry: `.translate()` works on `Seq` only, not
  `str`.
- **bio-data-hub**: UniProt `search` CLI quoting rules
  (`organism_name:"..."` inside single-quoted PowerShell strings);
  `organism_name:` field tip; STRING `network` "partners: N" is the
  `--limit`-truncated row count, not the full interaction count; STRING
  `preferredName` can be a weird alias (P10896 → `MTI20.21`) — the `map`
  mapping row is authoritative.
- **protein-modeling**: WSL NUL/UTF-16 error-table entry; ESMFold availability
  note + offline fallback; new section "overexpression donor-sequence checks"
  (see below).
- **protein-quality**: mapping-mode caveat, updated metrics table and
  quick-reference.
- **BACKENDS.md**: data-channel reliability table (Entrez/BLAST/Datasets,
  UniProt/PDB/STRING stable; MMseqs2 usable; NCBI FTP bandwidth-shaped; ESM
  Atlas frequently down).

### Backlog (accepted, not yet implemented)

- **Protein-family × species panel workflow** (item 8): family + species list →
  UniProt reviewed-first → TrEMBL → BLAST-TSA fallback for unannotated species
  (the Riccia fluitans RBCS/RCA recovery route) → identity matrix + heatmap +
  NJ tree. Candidate: new `analyze_panel.py` skill script.
- **Cloning primer design module** (item 9): cut-site conflict detection with
  compatible-enzyme fallback (NcoI→BspHI, SpeI→XbaI→BstEII), Tm_NN-based
  annealing-region length selection, CSV output. Candidate: extend
  biopython-analyses resources.
- **UniProt→CDS route** (item 10, partially documented): the checklist
  (EMBL cross-reference → GenBank mRNA preferred over `Entrez.elink`
  protein→nuccore; reject partial-CDS records missing the transit peptide,
  e.g. Q9AT38 vs full-length precursor AB034748.1) is documented in the
  protein-modeling skill; a helper script that automates the EMBL-cross-ref →
  GenBank → full-length CDS verification is backlog.

## Round 2 — day-long field test, MSA-outage recovery run (2026-08)

User-reported verdict: **结果可靠、过程透明、环境有小缺陷**. The day proved the
mode's resilience: with every MSA server down (official API blocked, backup
site abandoned, ESM Atlas 504), a manual MSA of 1557 UniProt homologs
recovered the structure — final TM-score 0.832, RMSD 1.76 Å vs crystal.

### Already fixed (verified this round, no code change)

1. **`esmfold_predict` tool result "value is not lossless JSON".** The wasted
   run predates v0.2.0. All 7 tools read backend JSON through the plugin's
   two-layer sanitizer (`protein-tools.js` `readJson`: NaN/Infinity tokens at
   text level + non-finite values at tree level), and `esmfold_predict` reads
   its report through `readJson` (L221) — verified, no gap. The backend writes
   PDB text only, so no `allow_nan` change applies there.

### Fixed this round (v0.2.1)

2. **matplotlib Chinese tofu blocks.** All three figure-producing scripts
   (`stat_plots.py`, `seq_logo.py`, `md_mmgbsa.py`) now auto-select a CJK font
   (Microsoft YaHei → SimHei → Noto Sans CJK SC → SimSun → DejaVu fallback
   with a stderr warning) and set `axes.unicode_minus=False`. Verified:
   matplotlib sees YaHei/SimHei/SimSun; heatmap with Chinese row labels
   (浮苔/菠菜/拟南芥) renders correctly. bio-visualization and
   biopython-analyses skill docs ship the same recipe for ad-hoc code.

3. **No molecular renderer on the machine.** PyMOL 3.1.0 open-source installed
   into `D:\bioai\venv` from cgohlke's cp313 wheel (GitHub release v2025.2.2,
   sha256 `62471028…`, via gh-proxy — github.com CDN is unreachable from this
   network). Route facts verified 2026-08: official pymol.org Windows bundles
   discontinued (installers directory emptied), PyPI carries only 3.2.0a0
   Windows wheels (user-verified broken). New headless renderer
   `skills/protein-modeling/resources/pymol_render.py`: 5 styles (cartoon /
   publication / rainbow / surface / line), hetatm element coloring, optional
   semi-transparent surface, 300 dpi ray-trace. Smoke-tested on the 1brs AD
   complex (1781 atoms, 3 styles). PyMOL 3.x API note: all `by*` color names
   are removed (`Unknown color`); ligand coloring uses `util.cnc` — recorded
   in the protein-modeling error table.

### Backlog (accepted; ready-to-run scripts shipped)

4. **Local MSA database (~70 GB) to end server dependency.** All MSA channels
   were unreachable today. `scripts/install-local-msa.ps1` installs mmseqs2 in
   WSL (apt → conda via TUNA conda-forge mirror) and downloads UniRef30 2302 +
   colabfold_envdb 202108 with resume support (reuses deploy/parallel-download.ps1,
   `-Proxy` option). Source status probed 2026-08 from this machine: GWDG
   Göttingen reachable at file level; steineggerlab workers and mmseqs.com
   time out; NJU/TUNA/USTC do not mirror the colabfold DBs (404). WSL has no
   mmseqs yet. Auto-wiring a local-MSA mode into `run_colabfold.ps1`/
   `af2_predict` is pending the DB download + an end-to-end smoke test; the
   two-step `colabfold_search` → `colabfold_batch` recipe is documented in the
   protein-modeling skill and printed by the installer.

- Renderer follow-ups (optional): predicted-vs-crystal alignment overlay,
  per-chain surface presets, batch PNG gallery.
