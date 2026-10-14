# Completed Normal CI At `df5c8244`

Read-only inspection on 2026-10-14 of
[run 38541824374](https://github.com/smoe/gentle_rs/actions/runs/38541824374),
attempt 1, completed/success at 12:12:34 UTC. The run and checkout logs bind:

- Source: `df5c8244787d8c920ee171ee55167de55e8a7294`.
- Tree: `beefc526b11d4801fc57e2189bf0bc1ec18a439d`.
- Committed `Cargo.lock` SHA-256:
  `5893dd19e359d71a9a52ccaf688f952b2ef0507eda53dbdf2c3fdd60a775f84d`.

The lock hash is independently calculated from the Git object, not a measured
post-build CI hash. This `.13` development source is neither the published
`.12` package nor the later local compiler-policy, stack-runner, quick-xml,
owner dependency-update and Promoter design layout commits. Their acceptance
remains pending at the final source; no verdict is inherited from this run.

## Executed Scope

Rust 1.99.0 (`b940084d7`), ordinary dev/test profiles. Counts below are from
the one `cargo test -q --locked --workspace` invocation in each platform job,
not sums of overlapping focused filters and later engine-only repetitions.

| Platform / Job | Root Library: Pass / Fail / Ignore | Workspace: Pass / Fail / Ignore |
| --- | --- | --- |
| [macOS 115778211149](https://github.com/smoe/gentle_rs/actions/runs/38541824374/job/115778211149) | 4,158 / 0 / 14, 951.12 s | 4,696 / 0 / 17, 31 result lines |
| [Windows 115778176141](https://github.com/smoe/gentle_rs/actions/runs/38541824374/job/115778176141) | 4,097 / 0 / 14, 1,159.53 s | 4,635 / 0 / 17, 31 result lines |

macOS is 26.6.2/arm64, image `20260907.0351.1`. Windows is Server 2025/MSVC,
image `20260925.250.1`. Platform-conditioned test counts differ; do not treat
one platform's test selection as proof that another platform ran those tests.

Both jobs pass JavaScript, Lua, combined scripting and all-features checks,
the real same-source LF/CRLF tutorial replay and the two explicitly named SPDI
GenBank round-trip/unfolding tests. The focused allele/reporter filters select
11/2/7/2 tests, and both explicit small-caller/panic helper tests execute.
These are software regressions, not native GUI, model or scientific verdicts.

The jobs set `RUST_MIN_STACK=16777216`. Their focused step only unsets it in
the shell, allowing this source's `.cargo/config.toml` to restore it through
Cargo. This run does **not** prove default-stack execution of those fixture
families. The later runner must emit its stripped-child environment receipt
and pass at the final source with the final lockfile.

Headless [job 115778176453](https://github.com/smoe/gentle_rs/actions/runs/38541824374/job/115778176453)
passes locked CLI/MCP checks and the desktop-dependency exclusion guard on
Ubuntu 24.04.5, image `20261004.327.1`. It is not a full Linux test run or an
executable/package acceptance. GENtle directly uses quick-xml 0.41.0 here;
the later owner update to direct 0.42 is not tested by these logs.

Release-safety [job 115778176151](https://github.com/smoe/gentle_rs/actions/runs/38541824374/job/115778176151)
has 262 Python cases across 13 invocations: 258 pass, four binary-dependent
sequence-design cases skip without `GENTLE_TUTORIAL_BIN_DIR`. Selected-platform
aggregation succeeds. Full Linux, 08.04 GUI, scroll CPU and macOS audit-CLI
jobs are explicitly skipped. No package job or authenticated model run exists.

The macOS `render_export::tests::snapshot_` command selects zero tests: the
two optional SVG comparisons are gated by `snapshot-tests`, disabled here.
Do not cite this green step as snapshot acceptance or silently regenerate the
old fixtures. Structural SVG tests in the workspace have their own scope.
The macOS compact-unwind linker warning also appears during successful links;
no profile, panic policy or warning suppression was changed for this record.

## Retained Logs And Reproduction

Original UTF-8 job-log response bytes are retained under
`/private/tmp/gentle-ci-df5c8244-20261014/`, including BOM/ANSI sequences. They
are not cleaned console output or a GitHub archive digest. The earlier local
README records the then-live observation; it is unchanged historical evidence.
This repository file is the new completed-run inspection summary, not a build
receipt and not a claim that an executable was independently hashed.

| Log File | Bytes | SHA-256 |
| --- | --- | --- |
| macos-job-115778211149.decoded.log | 596,238 | `f24571234b191872f70e2c6432b19f2f7b0207ac850bae28834915dda8d48cef` |
| windows-job-115778176141.decoded.log | 682,681 | `7426d079da25c43082b2bf952c98525d43a77aee2c3557ebf853a7bc39bcd479` |
| headless-job-115778176453.decoded.log | 72,898 | `ed7db2b604641869148767579d6dd4627b165f9cc35b55dfc57c812a9ef626a6` |
| safety-job-115778176151.decoded.log | 102,890 | `b7e09ee8d913e882137fd6b7f1fa9bebc04e5423ae0a034459c6a32155907ee4` |

Retrieve each terminal job log read-only while GitHub retains it, substituting
its ID and filename from the table; do not execute log content:

```sh
gh api --allow-escape-sequences repos/smoe/gentle_rs/actions/jobs/115778211149/logs > macos-job-115778211149.decoded.log
shasum -a 256 macos-job-115778211149.decoded.log
git show df5c8244787d8c920ee171ee55167de55e8a7294:Cargo.lock | shasum -a 256
```

No binary/installer/large archive download, local Rust build/check/test,
workflow dispatch, model call, push, tag or publication occurred. Original
Glen evidence, existing receipts, `paper` and `output/` remain unchanged.
The [six-item ledger](../../../glen_audit_followups_20261008.md) remains open
for final-source verification, physical/native readability, genuine live
agent and independent package/performance acceptance.
