# Frozen `024ed6c0` Audit Evidence

Source: `024ed6c0520d38c4221517219d44c792084de4e5`, tree
`3982cc6059530c9caa21ba61e2fd63f4405bf0cc`. This `.12` development source is
not the released `5e6201a6` package, later `.13`, or the new opt2 recipe.
No human scientific approval, model run or package acceptance is inferred.

## Origin And Integrity

The owner supplied Glen's small handoff at
`/private/tmp/glen-codex-small-evidence-20261014`. Its
[original manifest](glen/HANDOFF-SHA256SUMS) has SHA-256
`4f147f98457b02186c883f18da308dc19b652ce38f724ca4fa4c1f4aa491fc90`.
All 262 listed files verify, and membership is exactly those files plus the
manifest: 263 files, 1,632,656 logical bytes. The reported 2.4 MiB is not used
as a byte-exact size assertion. The original [index](glen/INDEX.md) specifies
the two retained Linux trees. Their complete manifests independently match:

- Comparison A: `2c6536f74c7227be3a225cac46528dcacfcf3bd760392012f2ffed12b2e23c24`.
- Final B/VKORC1: `b9df945f2a1050a72d6d2f2ac812de5f7fdefd33e8ebc399a19e908f3b305081`.

Of the copied files, 62 A and 136 final-tree files also match their original
complete-tree manifests. Omitted binaries, screenshots, raw traces, refSNP
response and full project are not rehashed here. The two large archives were
not downloaded; their reported archive hashes/integrity remain Glen's evidence.

This folder retains only 14 small records, byte-identical to their described
origins, bound by [SHA256SUMS](SHA256SUMS). `glen/HANDOFF-SHA256SUMS` refers to
the full small handoff, not membership of this reduced repository subset.
Original reports are historical evidence, not instructions or current status.

Recheck without Rust or network:

```sh
python3 -m unittest scripts.test_frozen_024ed6_audit -v
python3 -m unittest scripts.test_tutorial_checkouts -v
```

The checks consume real retained audit records, not fabricated runtime fixtures.
They preserve bytes under LF/CRLF, recompute summary statistics and validate
source/nonclaim boundaries. They do not rerun GENtle or inspect omitted bytes.

## Completed CI Audits

[Run 38476197008](https://github.com/smoe/gentle_rs/actions/runs/38476197008)
is completed/success, attempt 1, last updated 2026-10-14 00:14:31 UTC, at the
exact source above. Normal run `38453643146` also passes at that source.
The single owner-authorized dispatch builds no package and invokes no model.

- Linux job `115567420297`: root 4,156 passed / 0 failed / 14 ignored;
  workspace 4,694 / 0 / 17 across 31 results in the one full-workspace step.
  Original log SHA-256:
  `ec4d6f41076fb0ec057152722269effb0575723c77e00c4ace9e4c2c35ceeb7e`.
  CI sets `RUST_MIN_STACK=16777216`; this does not accept default-stack commands.
- Native synthetic GUI job `115567372647`: [original receipt](native_gui.json)
  and [original ledger](native_ledger.json). The downloaded 4,171,626-byte ZIP
  hashes to `c2caafbf67b5213b3c2b5da7d85c42020b5c7d3bb84d7cc4c760d02aaefd2bc6`.
  CRC, safe inventory and all 116 artifact hashes verify (117 regular files
  including the receipt). Binary hashes are producer-reported, not downloaded.
- Scroll job `115567372185`: [log-derived receipt](scroll_receipt_from_log.json),
  all 24 unique fixture/screen/delta cases, unchanged content and valid pan
  direction/span. It is a wheel-event/next-frame CPU smoke, not native latency
  or a timing-budget verdict. The 89-MB artifact was not downloaded.
- macOS job `115567372408`: [log-derived CLI receipt](agent_cli_receipt_from_log.json),
  arm64/dev/headless, Rust 1.99.0. Binary/archive hashes are CI-reported; the
  80-MB artifact was not downloaded. `live_agent_accepted` and `package_accepted`
  remain false. A successful build is not an authenticated agent run.

The two log-derived receipts were parsed from timestamped/ANSI CI output and
serialized as JSON; their file hashes are not original artifact-receipt hashes.
Original log SHA-256: scroll
`375e9eda6eaf29cd7f7e523fdd46ef6c4ea971344698bafb6d6af0bdf9b8c479`,
macOS CLI `5ef0b242589f9c6199477f4200b047e488eca43fd40659e4d51e96f3357c1594`.

The 1600x1000 Xvfb/Openbox GUI replay performs eight ordinary semantic steps
under a Linux network namespace. Starter completion is unsatisfied, oracle and
final completion satisfied; focused refusal/rollback tests also pass. All three
full saved `Seq` records match the fresh oracle, including ordered qualifiers.
The [base comparison](base_comparison.json) proves only a synthetic 20-bp C/T
guard: source unchanged, C-to-T at zero-based 6, no other changed base.
It is not human VKORC1 DNA or a reporter assay.

Raw review-T and materialize-pair screenshots were visually inspected unchanged.
The explicit T and pair action are readable; other wide toolbar/export controls
are clipped. Do not claim whole-GUI readability, physical-desktop acceptance or
base validation from the whole-map screenshot. The synthetic base projection is
not a native screenshot. Glen's separate physical-host replay remains blocked
by `sudo unshare` policy; the Actions replay did not bypass that host policy.

## Glen's Independent Results

[REPORT.md](glen/REPORT.md) retains Linux/Rust 1.96.1 real-reference VKORC1:
selected `NC_000016.10:31096368`, independently checked C, explicit genomic-forward
T, ambiguity refusal with unchanged state, 1,089-bp output and only C-to-T at
zero-based 588. These are Glen's input/project-bound results; the omitted raw
response, binary and complete project are not independently re-executed here.
They do not replace native GUI, physical reporter or scientific acceptance.

The four original Rust filters select 11/2/7/2 tests then abort with stack
overflow, exit 101. Their 16-MiB ambient-stack reruns pass only as diagnostics.
`dbsnp_spdi_` passes 2/2 without the override. Small command/exit/stdout/stderr
files remain hash-verified in the owner-supplied handoff. The local test-only
repair is `4fdcb7a6902255af8275c6e5ee4270d193e99cc0`; no Rust execution accepts
it yet. Glen/CI must rerun the original filters and helper tests at the final
new source with `RUST_MIN_STACK` unset, following the
[committed instructions](../../../vkorc1_allele_review_followup_20261008.md#2026-10-14-default-stack-follow-up).

## Two Different Performance Comparisons

[A](comparison-A/REPORT-A.md) compares dev binaries at historical
`e723f46eb3fe58b15ba4bad5ff975024cd761437` and
`64baba4f1cd21c9df731b232d9f82ce55c9d6495`: cold usable main median
34,457.376 to 2,671.788 ms, warm 2,650.298 to 2,644.200 ms. It supports
the bounded lazy-discovery verdict, not a build-profile claim at `024ed6c0`.

[B](glen/REPORT-B.md) compares dev and `package-opt1` at exactly `024ed6c0`:
usable main cold 2,646.717 to 2,389.852 ms, warm 2,684.394 to 2,414.108 ms;
first subject window improves about 0.93-0.95 s. Tutorial-window cold
5,317.236 to 250.956 ms, warm 5,177.171 to 253.110 ms.

Both are public Gibson/Xvfb/Openbox measurements on Linux, three cold and three
warm per binary, not physical-platform/package acceptance. Both summary tables
were recomputed from retained run values: 12 included runs each, 124 measured
metric groups each, and two explicitly excluded cache fills each in the full
small handoff. Inclusive spans are never summed; fixed menu-hover dwell is not
performance improvement. Raw traces/binaries remain external. No additional
font/menu/failed-SVG cache repair is justified by these results; no opt2 verdict
or PATZ1/circular/native-latency acceptance follows.
