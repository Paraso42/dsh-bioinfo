## Summary

<!-- What does this PR do and why? One or two sentences. -->

## Checklist

- [ ] `npm test` passes (plugin syntax + tool-schema validator)
- [ ] Python backends compile (`python -m py_compile` over `skills/**/*.py`)
- [ ] New tools register through `defineToolDef()` (never a flat parameter map)
- [ ] Machine-specific paths stay behind `BIO_TOOLS_*` env vars with canonical defaults
- [ ] `BACKENDS.md` updated (interpreter / argv / output contract)
- [ ] `THIRD_PARTY_NOTICES.md` updated for any new dependency (license + citation)
- [ ] Deployment/acceptance behavior verified (how?)

## Notes for reviewers

<!-- Anything unusual: schema pitfalls, license constraints, platform caveats. -->
