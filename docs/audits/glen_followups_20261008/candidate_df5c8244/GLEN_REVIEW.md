# Glen's Frozen Review And Native-X11 Supplement

Owner delivered both small trees on 2026-10-14. This record supplements the
[completed normal-CI inspection](README.md); it does not amend that earlier
record or Glen's original files. Tested source is exclusively
`df5c8244787d8c920ee171ee55167de55e8a7294`, tree
`beefc526b11d4801fc57e2189bf0bc1ec18a439d`, lock SHA-256
`5893dd19e359d71a9a52ccaf688f952b2ef0507eda53dbdf2c3fdd60a775f84d`.
Later source, the owner-updated lockfile and opt2 are not accepted by this audit.

## Retention And Integrity

Original small tree: `/tmp/glen-df5c8244-small-evidence-20261014/`.
All 94 handoff-listed files verified independently; the tree has 96 files,
including its manifest and producer's verification output. Retain a 40-file
subset under [glen/](glen/REPORT.md): original report, summary, environment,
cleanup, handoff manifest and all seven Rust command/exit/stdout/stderr/time
record sets. Its original handoff manifest SHA-256 is
`51bc128734d18f6ead063be80efffc7a24a43806a4dc83b7dab9a4b8be53a61a`.
The manifest lists omitted files too; do not interpret this subset as a full
archive. Original paths in records are historical, not runnable local paths.

Retain the entire 22-file [native-X11 supplement](native-x11-supplement/INDEX.md).
All 21 listed files verify; its handoff manifest SHA-256 is
`867e89d7e8ce5782707d45a572bd4d116d8a06e6427b2c2c2fe05880d42a6732`.
The producer reports archive hash
`861f43980ba230927bbe934a52c372e26fef23d222bc6d0b5e3f825e7d30a15a`;
no archive was supplied or rehashed here. No binaries, screenshots, states,
credentials or reruns are included. Binary/input hashes of omitted artifacts
remain producer-reported, not independent local verification.

Both full-tree manifest versions verify and are preserved separately:

| Version | Lines | SHA-256 |
| --- | --- | --- |
| Pre-cleanup | 194 | `0b6a9d95c689a7cb0f0ce2becf118b2a9c6192ed39738879f0b867cb39b16430` |
| Final | 195 | `713a8d034e2c7964d29e6350e3fa79d026726e7d28b15686d66597f5f9754e1a` |

Their only difference is the final `cleanup.txt` entry, hash
`0f7cfc444537b81e073a39068d393b80316016ad054f1500bd1a5da2c420aff8`.
The retained cleanup bytes match. Every earlier entry is unchanged; the
194-line companion still names its original path, now occupied by the final
version on Glen's host. Do not overwrite it or relabel either version. This
resolves the version discrepancy, not verification of the omitted full tree.

The scoped Git checkout policy disables conversion and whitespace linting only
for original records, including historical trailing spaces, mixed LF/CRLF logs
and empty streams. Source-code checks are unchanged. Offline integrity checks:

```sh
python3 -m unittest scripts.test_frozen_df5c8244_audit
python3 -m unittest scripts.test_tutorial_checkouts
```

These inspect original bytes and manifest membership, not GUI/Rust execution.

## Executed Scope And Stack Qualification

Glen used clean, isolated Debian x86_64 Linux, Rust/Cargo 1.96.1, dev/default
features, jobs=1, incremental off and dev debug=0. All six Rust filters report
28 passes, zero failures; locked `cargo check` exits zero. Both named SPDI
regressions execute. The reported 17:07 and 10,732,312-KiB peak belong to the
cold build, not test runtime or GUI rendering.

Glen's unchanged summary says `rust_min_stack:null`; the commands remove it
from the parent shell. At this frozen SHA `.cargo/config.toml` restores
`RUST_MIN_STACK=16777216` when absent. These are valid ordinary test passes,
but not proof that the families ran without an ambient stack override. Keep
the two explicit small-caller helper tests' bounded proof separate. The final
source still needs the corrected executable-boundary runner and its
`child_rust_min_stack:null` receipt. Do not rewrite Glen's original summary.

The producer reports 109 Python passes, 66 catalog entries, 30 executable
chapters, 53 workflows and no LF/CRLF/generated drift. Full projection streams
are omitted from this retained subset. Same-source CI completed successfully;
its independently inspected scope remains in the earlier README.

## Native Findings And Limits

| Area | Exact Frozen Verdict / Remaining Boundary |
| --- | --- |
| FLNA | Glen reports a normal-size virtual Linux/X11 presentation pass, genomic-minus/local-plus orientation, 36 transcripts/103 exons and both UniProt rows. Raw captures are omitted here and have not been inspected by Codex. Physical macOS/Windows, historical release provenance and scientific confirmation remain open; the draft stays inactive. |
| 07.04 | Public 83,686-bp locus overview opened. Opening its sequence window left X11 probes timing out; at 131 s GENtle was 36.0% CPU, 1,281,072 KiB RSS. Native navigation/results/reopen remain blocked in that environment. |
| 10.03 | Synthetic project and actual Patterns menu route were used. At 114 s GENtle was 34.1% CPU, 412,176 KiB RSS; X11 probe timed out. Native calibrated/null results/reopen remain blocked, not replaced by SVG. |
| 08.16 | Report-bound viewer adoption was not established. A positional GenBank argument is not a project import; no lane capture or renderer-hang conclusion follows. |
| 04.09 / 06.07 / 03.05 / 10.02 | Not run: no prepared exact-source GUI state; full primer import remains dependent on an evidenced performance repair. |

Both retained GUI stderr files contain only Mesa/EGL DRI3 acceleration warnings,
not Rust panics or CPU stacks. Openbox logs include a missing menu and an XIO
message; Xvfb logs label the keyboard warnings nonfatal. These do not establish
the causal blocker or justify an application/backend change. Process samples
are historical tool-call transcriptions, not new sampling or a CPU profile.
Keep graphics/backend limits, application work and harness/adoption errors
distinct until an exact-source observation identifies the stalled boundary.

No new native, default-stack, scroll, package, model or scientific acceptance
is claimed. No local Rust build/check/test, push, dispatch or large download
was performed. The six-item goal remains open at its documented boundaries.
