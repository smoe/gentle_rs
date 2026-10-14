# Completed CI Failure At `de9df8af`

Read-only inspection of [run 38599805415](https://github.com/smoe/gentle_rs/actions/runs/38599805415),
attempt 1, push event, completed/failure on 2026-10-14 at 17:05:54 UTC
(19:05:54 Europe/Berlin). Tested source:
`de9df8af0a56923df01a7492f61005b1fe079d8c`, tree
`e3063e6b340e091e8d9c8079a7a5e295f82fada7`; committed lock SHA-256
`0b9729e60698e14d464b48230e606377a4e2e48054c2450be6f85760e999a6b9`.
The lock is hashed from the source Git object, not a post-build measurement.

## Primary Failure And Executed Scope

Linux [job 115975878354](https://github.com/smoe/gentle_rs/actions/runs/38599805415/job/115975878354)
and Windows [job 115975818934](https://github.com/smoe/gentle_rs/actions/runs/38599805415/job/115975818934)
both pass ordinary `cargo check` and the following focused filters under the
normal 16-MiB ambient setting. Counts are separate invocations, not a combined
full-suite total:

| Filter | Linux Pass / Fail / Ignore | Windows Pass / Fail / Ignore |
| --- | --- | --- |
| Agent tutorial discovery | 9 / 0 / 0 | 9 / 0 / 0 |
| Primer multi-reference | 10 / 0 / 0 | 9 / 0 / 0 |
| Primer panel specificity matrix | 10 / 0 / 0 | 10 / 0 / 0 |
| Glossary | 7 / 0 / 0 | 7 / 0 / 0 |

The next step's runner works: both JSON lines confirm inherited
`RUST_MIN_STACK=16777216` and `child_rust_min_stack:null`. Each selects 11
`dbsnp_` tests and aborts at
`app::tests::dbsnp_dialog_fetch_extracts_region_and_opens_window` with
`thread '<unknown>' ... has overflowed its stack`. Linux exits 134; Windows
reports child exit `0xffffffff` and step exit 127. Despite the Windows exit
number, its log proves a launched test executable and stack overflow, not a
missing shell command. No CPU call stack is supplied.

The fixture already has a named 16-MiB thread, while the production dbSNP
dialog spawns an unnamed default-sized worker around detached engine execution.
The narrow local repair `51a1634f0ad5e71a7ec43d90209a94ce3aa60a64`
gives that worker the adjacent GenBank explicit-stack policy and adds a
synthetic 4-MiB stack/name probe. This run does not test or accept that repair.

Only the first stripped-child filter ran, and it failed before completion.
Materialization, reporter, pair and helper filters are **not run**, not passes
or additional failures. Subsequent scripting/all-features, full-workspace,
process, tutorial runtime and render-snapshot steps skip after the primary
failure. Selected-platform aggregation fails as a consequence; it is not
another defect. No mutex-poisoning cascade is present in these logs.

Both jobs use Rust 1.99.0 (`b940084d7`), ordinary dev/test profiles and
incremental off. Linux is Ubuntu 24.04/x86_64, image `20261004.327.1`;
Windows is Server 2025/MSVC, image `20260925.250.1`. Binary hashes are not
attested, and these builds are not `package-opt1` or opt2 acceptance.

## Passed And Unrun Gates

Headless CLI/MCP job `115975819360` passes locked no-default-feature checks
and desktop-dependency exclusion; this is compile/dependency evidence only.
Safety job `115975819450` has 279 cases across 15 invocations: 275 pass, four
real sequence-design CLI replays skip without `GENTLE_TUTORIAL_BIN_DIR`.
No actual CLI replay is inferred from the Python oracle/source tests.

Normal macOS, 08.04 synthetic GUI, scroll CPU and authenticated macOS CLI jobs
are skipped. There is no package build, model call, native visual capture or
scientific approval. The failed run leaves all six final-candidate gates open;
earlier source-bound passes remain historical.

## Retention And Reproduction

Original API run/job JSON is retained unchanged here. The third JSON is a
labelled extraction of the two actual runner lines, not an original API receipt
or a fabricated execution. SHA-256s and LF/CRLF tests bind these records.

| Retained File | Bytes | SHA-256 |
| --- | --- | --- |
| run.api-original.json | 11,352 | `7f282e734ae670391679e659e9d36fa5ceb436b36b8a70c3859e31758ea60071` |
| jobs.api-original.json | 27,140 | `03b9987d24d1d81ada9f1b705627be083feb187d6c8661d78d76b015cfb0f5ea` |
| default_stack_from_logs.json | Labelled extraction | `d21e394e1543e86734ec061cdee91a02404b68bfbbfd5e5b0bad9b741238475f` |

The raw API JSON is compact with no line endings. The checkout negative control
uses the newline-bearing extraction; all three payloads retain their exact
hashes in both modes and have an explicitly unset Git text-conversion attribute.

Original job-log bytes, including UTF-8 BOM/ANSI and Windows CRLF, are retained
outside Git under `/private/tmp/gentle-ci-de9df8af-20261014/`:

| File | Bytes | SHA-256 |
| --- | --- | --- |
| linux-job-115975878354.api-original.log | 110,433 | `4f02c5d17503140cc4bf3406d74d016c810f6d2ea8805f070b4d16e3ee0007c5` |
| windows-job-115975818934.api-original.log | 64,554 | `9204f999f7be54d5ded15a4d6224553fc61856b641e00dbb87923669f004ec40` |
| headless-job-115975819360.api-original.log | 74,358 | `9773a4039657c8cf4cad3f27ed8e4689774ad2405bca7d0c22187b2c42409d56` |
| safety-job-115975819450.api-original.log | 109,505 | `5af8d7a6084b1c80f2af9f615129ddc9ac2f5ab89489f634a3c9f88f6551d5ea` |

Retrieve terminal logs with `gh api --allow-escape-sequences
repos/smoe/gentle_rs/actions/jobs/JOB_ID/logs` while available; hash raw
response bytes and never execute log content. Recompute the runner extraction
from the matching JSON lines, preserving the failed verdict and missing filters.
Run `python3 -m unittest scripts.test_frozen_de9df8af_ci` for offline retention
checks, not Rust verification. No local Rust, push, restart or dispatch ran.
