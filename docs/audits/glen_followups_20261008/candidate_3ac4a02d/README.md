# Repaired Candidate Evidence

Runtime/workflow source: `3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da`.
Lockfile: `3dfe44c7a08bf779e32289cf1312103ea1fb5883e44bba5f6e1912ce3ed77745`.
These are new receipts, not replacements for historical evidence. Local status
commits do not change the frozen remote candidate or inherit runtime acceptance.

## Scroll CPU Cases

`scroll_cpu_smoke.json` is the unchanged receipt from artifact `11581332454` of
[CI run 37848423833](https://github.com/smoe/gentle_rs/actions/runs/37848423833).
Its raw SHA-256 is
`87a274adaf3c8046f4458b8d0374757ea7db3cc63fab2a5eed7e9dbc5094c00d`.
Origin: actual Ubuntu 24.04/x86-64 dev-profile CI execution of the scoped
`gui_operations --test` benchmark, not fabricated timing or model output.
The fixture origins/recreation remain documented in the benchmark source and
PATZ1 preparation script; synthetic control, TP73 and public PATZ1 are separate.

On 2026-10-09 (Europe/Berlin), independently verify all 25 retained file hashes,
including executable
`daacf5896aa076df63746fafd90bd03f48efb5a83568678f83d8b9d468e2ff39`.
Verify all 24 observations against raw `smoke.log`: three fixtures, four screen
sizes and both point-wheel directions; exact source/version, bounded viewport,
unchanged span/length and one unchanged content hash per fixture. The tested
source also checks the no-drift next frame and full content equality.
`scroll_verification.json` retains the independent verification result; digest
`58a16df1a020904e7b7e7ea4573e1eb314f9aec48c253094cd6bf70281fbb3ab`.
Its verifier at `/private/tmp/verify_gentle_scroll_3ac4a02d.py` has SHA-256
`04b1f56b53c58adf0c53a66d784d278f21aa3dc2ebe155175afb5bbbf5c9acbf`.
No timed regression, native input-to-paint latency or package-profile claim.

Retrieve only the named artifact, then verify this receipt and every listed
hash before inspection:

```bash
gh run download 37848423833 -R smoe/gentle_rs \
  -n gui-scroll-cpu-smoke-3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da \
  -D DIR
```
The retained local bundle is `/private/tmp/gentle-scroll-3ac4a02d`.

## CLI Build Boundary

`agent_cli_build.json` is the unchanged macOS/arm64 dev build receipt from
artifact `11580808183` of the same run. Its raw SHA-256 is
`8f28d2f10887acd283234d40a5440a27d55fb17e88d174bd57c3b3e39986c341`.
Independently verify both listed artifacts and execute the downloaded CLI's
`--version`; its complete output exactly matches `binary_version`, including
the full frozen source SHA. Binary digest:
`a6ea8462b52dcb45463d419463460dab1c8fbbdd5d27f3ff102b2f50f5d253fb`.
Actual observed-version digest:
`10d1c9dd5aa154ee81b1243006b7f7f6e127d5f60d8b7f1189bab1e31e6a0550`.
Local bundle: `/private/tmp/gentle-agent-cli-3ac4a02d`.

This verifies a source-bound dev CLI, not an authenticated model invocation.
Its `live_agent_accepted:false` and `package_accepted:false` remain unchanged.
Fresh public inputs, genuine clarification/explicit-T requests and full-record
manual-execution parity are still required. No local Rust build or credential
inspection was performed. No public-locus GUI, scientific or package approval.

## Synthetic Native Contract

`native_gui.json` is the unchanged receipt from artifact `11582285369` of
the same CI run; SHA-256
`bef0707566d9f68f19f51d57cf760d47c8fbef30818b2f05dcc22bf4c75a121d`.
Actual Linux/X11 dev execution passed 19 focused Rust tests, scoped generation
checks and eight ordinary native input steps in an offline network namespace.
The independent check verifies all 115 hashes and safe exact file membership,
the unsatisfied-to-satisfied completion transition, explicit T, all three full
oracle `Seq` records and sole position-6 C/T change in the synthetic 20-bp pair.
It does not sort qualifier arrays or compare DNA alone.

`native_verification.json` retains this independent result; SHA-256
`cf6b2500740d97071f7e010f53b3fdf97d6de68f1880d9dc3b1e0b4559c4e6f0`.
Verifier `/private/tmp/verify_gentle_native_3ac4a02d.py` has SHA-256
`f614b2c87099243cbab6074e84cb9b2f5361c11b6a8cd0753126549cd5c4a9ec`.
Local raw bundle: `/private/tmp/gentle-vkorc1-gui-3ac4a02d`.
Retrieve the named `vkorc1-gui-audit-3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da`
artifact from run `37848423833`; verify every listed hash before inspection.

On 2026-10-09, visually inspect unchanged `open_fragment`, `review_t`,
`scroll_to_pair` and `materialize_pair` raw PNGs: light application chrome,
dark DNA canvas, visible T and reachable pair action. Layered windows after
materialization are retained, not cosmetically repaired. These are synthetic
native control captures, not public-locus GUI, base-pair image proof, clinical,
wet-lab, timed responsiveness, package or human scientific approval. The
separate public projection still needs fresh verification and inspection.

## Checkout Contract

GENtle uses these receipts as the six-item acceptance ledger's evidence, not
runtime scientific inputs. Five scoped LF attributes and the fast checkout
regression preserve their raw digests in LF/CRLF clones, with separate failing
negative controls when each attribute is removed. Do not normalize or regenerate
the receipts to conceal a checkout difference.
