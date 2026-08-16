# Contributing to dsh-bioinfo

Thanks for your interest! This repository is the open-source release of the
生信模式 (Bioinformatics Mode) agent preset for DeepSeek Harness.

## Quick development loop

```powershell
# nothing to install for the plugin layer; Node >= 18 is enough
npm test                          # node --check + tool-schema validator

# syntax-check every Python backend (any Python 3)
python -m py_compile (Get-ChildItem skills -Recurse -Filter *.py).FullName
```

- `plugins/protein-tools.js` — the DSH plugin (7 model tools).
- `skills/*/resources/*.py` — the scientific backends the tools call.
- `deploy/` — PowerShell/bash deployment scripts (canonical layout: `D:\bioai`).

## Rules for adding or changing a tool

1. Register **through `defineToolDef()`** — never pass a flat parameter map to
   `ctx.tools.register` (upstream bug: the raw preset-file path ships
   `definition.parameters` verbatim to the model API, which rejects schemas
   without a root `type: "object"`).
2. Keep every machine-specific path behind the `BIO_TOOLS_*` environment
   variables with the canonical-layout default (see the config block at the
   top of `apply()`).
3. Update `BACKENDS.md` for the new backend (interpreter, argv, output
   contract) and `THIRD_PARTY_NOTICES.md` for any new third-party dependency.
4. `npm test` must pass; add your tool to the `EXPECTED` list in
   `validate-plugin-schemas.js`.
5. Windows/pwsh is the only supported shell today. If you add a Linux/bash
   backend, mirror it for both shells and note it in the README.

## Pull requests

- One logical change per PR; describe what and why.
- Use the pull-request template.
- If you touch deployment scripts, say how you verified them (the acceptance
  runner is `deploy/run-acceptance.ps1`; reference values live in
  `fixtures/README.md`).
- License: your contribution is accepted under the repository's MIT license.
  Do not commit code you cannot license this way; ACADEMIC-ONLY components
  (TM-align, PRODIGY, ESMFold weights) must never be bundled — they are
  documented in `THIRD_PARTY_NOTICES.md` and fetched at install time.

## Language

README/BACKENDS/notices are English; `docs/INSTALL.md` and the skill
documents are Chinese (the project's primary user base). Either is fine in
PRs — but keep each document internally consistent.
