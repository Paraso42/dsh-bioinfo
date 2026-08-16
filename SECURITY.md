# Security Policy

## Supported versions

Only the default branch (`main`) of this repository is supported. The
repository is a **preset kit** — it contains no server-side services, no
credential handling, and no network listeners. The deployment scripts under
`deploy/` execute on your machine with your permissions; review them before
running, as you would any shell script.

## Reporting a vulnerability

Please report security problems privately via GitHub's advisory system:

1. Go to the **Security** tab of this repository.
2. Choose **Report a vulnerability** (draft security advisory).
3. Describe the issue, the affected files, and how to reproduce it.

Do **not** open a public issue for a suspected vulnerability. We will respond
in the advisory, and — with your agreement — publish a fixed release and a
public advisory after the fix is available.

## Notes

- No API keys, tokens, or passwords are stored in this repository (the
  backends use keyless public services; credentials, if ever needed, must be
  supplied by the operator via environment variables).
- Keep machine-specific secrets out of commits: never add files like
  `.credentials.yaml`, proxy configs, or personalized SSH/Python caches.
