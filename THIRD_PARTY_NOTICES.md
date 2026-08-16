# Third-party notices

Code in this repository is MIT ([LICENSE](LICENSE)). The backends and the
deployment toolchain use third-party software itemized below; provisioning the
environment makes their terms yours to comply with. Items marked
**ACADEMIC-ONLY** must not be used commercially without the listed permission.

## ACADEMIC-ONLY components

### TM-align / TM-score (Zhang lab) — `struct_eval.py`
`skills/protein-quality/resources/struct_eval.py` contains a faithful
re-implementation of the TMalign algorithm (Zhang lab), cross-validated
against the official binary (max error 0.02). TMalign is **free for academic
use only**; commercial use or redistribution requires permission from the
Zhang lab. The official TMalign source is **not distributed** in this
repository; for the cross-validation binary install `tmalign` from bioconda
(academic use).
Cite: Zhang Y, Skolnick J. TM-align: a protein structure alignment algorithm
based on the TM-score. Nucleic Acids Res. 2005;33(7):2302-2309. PMID 15808260.

### PRODIGY (Bonvin lab, Utrecht University) — `prodigy_affinity.py`, `virtual_screen.py --prodigy`
prodigy-prot 2.4.0 / prodigy-lig 1.1.4 (PyPI) are **free for academic use**;
commercial use requires a license from Utrecht University.
Cite: Xue LC, Rodrigues JP, Kastritis PL, Bonvin AM, Vangone A. PRODIGY: a web
server for predicting the binding affinity of protein-protein complexes.
Bioinformatics. 2016;32(23):3676-3678. PMID 27153664.

### ESMFold / ESM Metagenomic Atlas (Meta)
`esmfold_predict` calls the free ESM Atlas API (no key; the service's terms
apply). Local ESMFold model weights are released by Meta under a
**non-commercial research license**: obtain them yourself (Hugging Face
`facebook/esmfold_v1`) and do not redistribute.
Cite: Lin Z, et al. Evolutionary-scale prediction of atomic-level protein
structure with a language model. Science. 2023;379(6637):1123-1130.
PMID 36927031.

### KEGG (Kanehisa labs) — `bio-data-hub/kegg_fetch.py`
KEGG REST access is free for academic use; commercial/redistribution use
requires the KEGG license: https://www.kegg.jp/kegg/legal.html

## Permissive components

| Component | Reference version | License | Citation |
|---|---|---|---|
| Biopython | 1.88 (venv) / 1.87 (`D:\biopython`) | MIT / BSD-3-Clause dual | Cock PJ et al. Bioinformatics 2009;25(11):1422-3. PMID 19304878 |
| AutoDock Vina | 1.2.7 binary | Apache-2.0 | Eberhardt J et al. J Chem Inf Model 2021;61(8):3891-8. PMID 34278794 |
| RDKit | 2026.3.5 | BSD-3-Clause | https://www.rdkit.org |
| meeko | 0.7.1 | LGPL-3.0 | https://github.com/forlilab/Meeko |
| OpenMM | 8.5.2 | MIT | Eastman P et al. PLoS Comput Biol 2017;13(7):e1005659. PMID 28746339 |
| MDTraj | 1.11.1 | LGPL-2.1 | McGibbon RT et al. Biophys J 2015;109(8):1528-32. PMID 26488642 |
| ParmEd | 4.3.1 | LGPL-3.0 | https://github.com/ParmEd/ParmEd |
| FreeSASA | 2.2.1 | MIT | Mitternacht S. F1000Res 2016;5:189. PMID 26973785 |
| gemmi | 0.7.5 | MPL-2.0 | https://gemmi.readthedocs.io |
| LocalColabFold | 1.5.5 | MIT | Mirdita M et al. Nat Methods 2022;19(6):679-82. PMID 35637307 |
| AlphaFold parameters | 2022-12-06 / 2022-03-02 / 2021-07-14 | CC-BY-4.0 (DeepMind) | Jumper J et al. Nature 2021;596(7873):583-9. PMID 34265844 |
| ESM-2 / fair-esm | 2.0.0, `esm2_t6_8M_UR50D` | MIT (code & weights) | Lin Z et al. Science 2023 (above) |
| LightDock (optional) | pip | GPL-3.0 | Jimenez-Garcia B et al. Bioinformatics 2018;34(20):3461-9. PMID 29718115 |
| pyCirclize | 1.10.1 | MIT | — |
| Logomaker | 0.8.7 | MIT | Tareen A, Kinney JB. Bioinformatics 2020;36(7):2272-4. PMID 31821414 |
| seaborn / matplotlib / pandas / numpy / scipy | reference venv | BSD-3 / PSF / BSD-3 / BSD-3 / BSD-3 | — |

## Data services used by backends

- ESM Metagenomic Atlas API (`esmfold_predict`) — Meta service terms.
- MMseqs2 server (`af2_predict` MSA mode) — free academic API; the server
  software itself is GPLv3 (API use unaffected).
- NCBI BLAST API (`biopython-analyses/ncbi_blast.py`) — NCBI usage policies.
- UniProt REST (CC-BY-4.0 data), RCSB PDB (CC0), STRING-DB (CC-BY-4.0),
  KEGG REST (see above).

## Not distributed

- `vina_1.2.7_win.exe` — download from the official AutoDock Vina release.
- TMalign source/binary — academic-only; install via bioconda for validation.
- AF2 params (~7.6 GB) and ESM/ESMFold weights — fetched by the deploy
  scripts at install time.
- WSL rootfs / Miniforge installers — fetched per `docs/INSTALL.md`.
