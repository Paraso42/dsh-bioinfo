# Acceptance fixtures

Small positive controls copied from the reference deployment
(`D:\bioai\jobs\acceptance`, ~3.3 MB) so a fresh replica can be compared
against known-good values. Inputs are for re-running the same checks; outputs
are the reference results from the machine this preset was developed and
acceptance-tested on.

| File(s) | Purpose | Reference result |
|---|---|---|
| `1brs_complex.fasta`, `1brs_ref.pdb`, `1brs_AF_complex.pdb`, `1brs_AD_complex.pdb` | AF2-Multimer input / crystal & prediction models | crystal baseline: 55 contacts, BSA 1280.9 Å² |
| `af2_out/` | reference AF2-Multimer run (ranked PDB + PAE/pLDDT + scores JSON) | layout produced by `deploy/run-acceptance.ps1` |
| `eval_af2_vs_crystal.json`, `eval_self.json` | `struct_eval` TM-score / RMSD / lDDT / GDT / DockQ | TM-score vs official TMalign: error ≤ 0.02 |
| `affinity_1brs.json` | `prodigy_affinity` (prodigy-prot) | 1brs ΔG = -11.3 kcal/mol |
| `mmgbsa_1brs.json` | `md_run --mode gb` (OpenMM MM-GBSA) | dG_bind = -36.7 kcal/mol |
| `3ptb.pdb`, `3ptb_ben.pdb`, `liglib.csv`, `vscreen_3ptb_report.json` | Vina virtual-screening positive control | benzamidine ranked #1 (-5.90 kcal/mol) |
| `props.csv`, `sim_benzamidine.csv`, `aspirin.png`, `aspirin_confs.sdf`, `aspirin_sub.png` | `mol_tools` (RDKit) | — |
| `deg_demo.csv`, `heat_demo.csv`, `circos_demo.csv`, `volcano.png`, `ma.png`, `heatmap.png`, `circos.png` | `stat_plots` | volcano 40 up / 20 down; heatmap 24×6 |
| `barnase_accs.txt`, `barnase_homologs.fasta`, `barnase_aligned.fasta`, `barnase_logo_bits.png`, `barnase_logo_prob.png` | biopython-analyses + `seq_logo` | barnase family logo |
| `1yph.pdb` | bio-data-hub PDB fetch sample | — |
| `string_tp53.json` | STRING-DB fetch sample | — |
| `esm_embed_out/` | `esm_embed` (ESM-2 t6) | per-chain 110 / 90 × 320 embeddings |
| `sanitized_1brs.pdb` | `md_mmgbsa` PDB sanitization | — |

Excluded intentionally (large / non-distributable / regenerable): MD
trajectories (`*.dcd`), TMalign source, the Vina binary, AF2 params.
