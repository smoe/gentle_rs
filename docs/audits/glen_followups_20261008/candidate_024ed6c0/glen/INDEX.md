# Glen small evidence handoff for Codex

Prepared 2026-10-14 as a file-level handoff. This directory contains copied,
byte-identical text/JSON/CSV receipts only; it is not a replacement for either
retained evidence tree and does not retarget their source revisions.

## comparison-A

- `REPORT-A.md`: lazy-discovery comparison at exact before/after revisions
  `e723f46e` and `64baba4f`.
- `analysis-A/`: included-run tables and exclusion-aware statistics.
- `runs/`: the 12 included runs plus the two excluded cache-fill runs used by
  the final warm measurements. Only metadata, external timings and trace
  summaries are copied; images and raw traces remain in the retained tree.
- `SHA256SUMS`: original complete-tree manifest.

Original retained root:
`/home/clawbio/.openclaw/workspace/artifacts/gentle-help-startup-audit-20261013`

## final-024ed6c0

- `REPORT.md`: consolidated VKORC1/dbSNP and comparison summary.
- `REPORT-B.md`: dev versus package-opt1 comparison at exact `024ed6c0`.
- `analysis-B/` and `help-startup-B/runs/`: included-run tables, statistics,
  metadata, external timings and trace summaries, including excluded fills.
- `focused/`: exact focused commands, exit codes and stdout/stderr, preserving
  both original stack-overflow failures and labelled diagnostic reruns.
- `ci/`: normal-CI run/job receipts.
- `vkorc1/`: retained input hash manifest and small real-reference command,
  placement, identity and comparison receipts. The 2.6-MiB project state and
  raw refSNP payload are intentionally omitted.
- `SHA256SUMS`, `MANIFEST-IDENTITY` and `manifest-verification.txt`: original
  complete-tree manifest identity and verification record.

Original retained root:
`/home/clawbio/.openclaw/workspace/artifacts/gentle-consolidated-audit-024ed6c0-20261013`

## Integrity

`HANDOFF-SHA256SUMS` binds every copied file in this small handoff. Paths are
relative to this directory. The original `SHA256SUMS` files continue to bind
the complete retained trees, including omitted binaries, screenshots and raw
traces.
