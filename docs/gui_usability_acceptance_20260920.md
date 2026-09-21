# GUI usability acceptance after the `.11` source update

Date: 2026-09-20

Upstream baseline: `885fac493a091313e7777a73de3736af4aae81c2`

Final tested code: `677c796ab0269b67b6eeece3d69e4112cd57e97b`

## Scope and verdicts

- **Updated source contracts:** accepted for the focused assembly-bound
  annotation, source-comparison, PATZ1, tutorial and release-version checks.
- **Synthetic TSS GUI lifecycle:** accepted after restoring the two previously
  reviewed GUI fixes. Two independent, network-isolated runs completed all 18
  checkpoints from the starter project.
- **General GUI usability/performance:** not accepted. The fixed workflow is
  operable, but this does not establish a responsive or generally usable GUI.
- **`.11` release tag:** not ready at this source revision. The package version
  remains `0.1.0-internal.10`; the release coordinator correctly rejects the
  proposed tag `v0.1.0-internal.11` until an explicit version/metadata commit is
  made.
- **Release artifacts:** no tag was created and no release workflow was
  dispatched. The artifact release therefore remains a separate exercise.

## Frozen environment and binaries

- `rustc 1.95.0 (59807616e 2026-04-14)`
- `cargo 1.95.0`
- Linux `7.1.3+deb14-amd64`, `x86_64`
- `Cargo.lock` SHA-256:
  `67a96b07e6e2d5beb5b0fd1eb2ba57e1c11b6d1a4a355dab5bce6ca830a25236`
- final debug GUI SHA-256:
  `d578f7120e5313cf6d311a1682a6320226087c6394494742a2b2cbad3a6f068f`
- final debug CLI SHA-256:
  `dbd97e1ed0cf280569a022685898fce5ba9614eb6ef883ecaec0021ae3a4653f`
- final debug documentation helper SHA-256:
  `5349545173b44796f7aad9758570f89cd431ec7bd1e9146ff3a2a7a90ea684c1`

The GUI runs used fresh HOME/XDG/cache/runtime/temp directories, Xvfb at
`1600x1000x24`, Openbox as an EWMH window manager, explicit X11 and screenshot
tools, and a fresh Linux network namespace. No private data or inner-agent path
was used.

## Reproduced defect and narrow fixes

At the upstream baseline, the isolated run failed at `open_tfbs_menu`: the
registered TFBS control was present in semantic state but not reachable by
ordinary X11 input. Baseline ledger SHA-256:
`e4570531faa1f39f11b38c837f6be8d25c7606919ecbd875540f05887c865acb`.

The final branch restores:

1. a usable initial/minimum native DNA viewport and bounded, resizable toolbar
   allocation;
2. a named 16 MiB worker for TSS snapshot/preview work, with bounded failure
   reporting;
3. an explicit root repaint after native child windows enqueue TSS collection
   open requests.

These presentation/scheduling changes do not alter sequence content, project
state, collection approval, stale-member detection or scientific verification.
Their four fully qualified Rust regressions pass.

## GUI acceptance evidence

| Run | Result | Ledger SHA-256 | Sum of recorded step time |
| --- | --- | --- | ---: |
| baseline debug | failed at step 2 | `e4570531faa1f39f11b38c837f6be8d25c7606919ecbd875540f05887c865acb` | 32.2 s |
| final debug A | 18/18 pass | `5656aa454da36bbbd1ed95263e99c39664e40b18796734c8d15377cf89585157` | 68.7 s |
| final debug B | 18/18 pass | `36daffd8c66d84edd52e37fc289e48ffc4c2ef9ef40ca39ee0ec6ee397450f4d` | 69.2 s |

Both final runs verified typed reports and oracle sequence content as well as
visible claims. They cover three 701 bp windows, shared transcript membership,
strand/local-TSS geometry, exactly four subject-bound viewers without
duplicates, cancel/Forget semantics, current persisted metadata and Undo.

The timings are acceptance-harness observations, not a formal performance
benchmark: screenshot/snapshot capture and fixed waits are included. However,
both final runs spent about 32 seconds in application startup plus opening the
first locus. The remaining operations were usually about 1--3 seconds each.
Passing the contract therefore does not answer the user's broader lag report.

## Authentic PATZ1 rendering and resize benchmark

The follow-up benchmark code was tested at
`2d6789d1f29a2157629e1215b271dde8f1947d2a`. It loads the public, offline PATZ1
project through the real `ProjectState` path and the portable source-comparison
report through the same binding-checked loader as the GUI. Setup fails unless
the fixture contains 20,802 bp, at least 75 loaded features, 13 Ensembl
transcripts and 17 joined source records (13 Ensembl and four RefSeq). All
source-comparison rows are expanded in the real DNA-window presentation; the
benchmark does not replace the locus with a feature-free surrogate.

Optimized `bench-audit` Criterion means on this host were:

| Operation | Viewport | Mean |
| --- | --- | ---: |
| eager DNA-window construction | n/a | 27.852 ms |
| deferred UI-thread hydration | n/a | 32.742 ms |
| first embedded frame | 820x520 | 29.469 ms |
| first embedded frame | 1200x800 | 29.496 ms |
| first embedded frame | 1600x1000 | 28.396 ms |
| first embedded frame | 1920x1080 | 28.098 ms |
| steady embedded frame | all four sizes | 0.786--0.793 ms |
| first frame after resize from 1200x800 | to 820x520 | 3.944 ms |
| first frame after resize from 1200x800 | to 1600x1000 | 3.839 ms |
| first frame after resize from 1200x800 | to 1920x1080 | 3.834 ms |

These are CPU-side egui preparation/painting costs. They exclude native window
creation, GPU upload, X11 event delivery and compositor scheduling. Viewport
size has little effect on PATZ1's steady paint cost, while initial
construction/hydration and first paint are well above a 16.7-ms frame budget.

Two fresh-profile Xvfb/Openbox runs then resized the real subject-bound PATZ1
DNA viewer through `820x520`, `1200x800`, `1600x1000`, `1920x1080`, and back to
`820x520`. Semantic snapshots confirmed every requested client size. The debug
GUI took respectively:

- run A: 409.6, 63.1, 229.4, 163.7 and 163.3 ms;
- run B: 425.7, 62.3, 187.3, 203.7 and 183.4 ms.

Those native timings include the `gui-test-support` semantic snapshot writer
and therefore are not release-binary frame times. They nevertheless reproduce
visibly laggy resize-to-confirmed-content behavior and show that the missing
time lies outside the optimized steady egui paint loop. A discarded harness
attempt wrote per-frame snapshots directly to CIFS and blocked in kernel I/O;
valid runs used local NVMe and copied evidence only after exit.

The optimized benchmark's first Thin-LTO build/link was itself expensive; the
measured matrix ran quickly once the binary existed. This supports a separate
`.12` task for a prebuilt or faster-linking GUI performance runner, without
misclassifying build time as runtime latency.

The hash-bound external benchmark bundle is retained under
`/mnt/storage-box-1-gentle/Glen/gui-benchmark-patz1-20260920-2d6789d1/`.
Its `SHA256SUMS` manifest has SHA-256
`69c89f6b802f52bd577522c77c6ce7383d47f720bd2d01776b3ca1498f76b2ec`.
After upstream advanced to `6aa71b113181d7c61181154175560067f97b570b`,
the two benchmark commits were rebased. The measured benchmark/window source
files remained byte-identical, and all 39 non-statistical benchmark cases
passed again at the post-rebase smoke revision
`15611fd1151b2aa1f4981df3f84622ebac4f422c`.

Raw public synthetic evidence is retained outside the repository under
`/mnt/storage-box-1-gentle/Glen/gui-usability-885fac49/`. Its 652-file
`SHA256SUMS` manifest verifies and has SHA-256
`ad9389224e7b85146700551b759f67e3bb80a0b6a4ba81e6823233f3905e331e`.

## Deterministic checks

- assembly-bound local projection: 3/3
- source-annotation presentation: 3/3
- authentic PATZ1 focused case: 1/1
- release-version consistency: 4/4
- Python GUI/screenshot/checkout checks: 33/33
- restored GUI regressions: 4/4
- component-crate suite: 267 passed, 2 explicit visual writers ignored
- release-candidate/package coordinator suite: 37/37
- catalog: 58/58 entries
- tutorial manifest/check: 29/29 chapters
- locked Cargo check, formatting and whitespace: pass

The tutorial checker retains two pre-existing human-review staleness warnings;
this work does not renew those reviews.

## Incomplete optimized-run investigation

Two optimized local attempts were deliberately stopped and are not counted as
passes or failures:

- `bench-audit`/Thin-LTO built dependencies and the optimized library, but an
  individual binary link was still active after about 30 minutes;
- reusing debug dependencies while compiling only GENtle at `opt-level=2`
  remained in the monolithic root crate after more than 35 minutes.

Both processes were CPU-active and emitted no compiler error. This establishes
a build-feedback problem, not a runtime verdict. A dedicated UX performance
profile must be evaluated as its own small build-system change rather than
being guessed into this acceptance fix.

## `.12` priorities

The [roadmap](roadmap.md#12-priorities) owns the release ordering.
The [DNA feature latency plan](dna_feature_rendering_latency_plan.md) owns
the implementation sequence, scope decisions and audit gates, including startup
and build-feedback work. This document retains historical acceptance evidence,
not a second mutable task list.
The PATZ1 results narrow the next investigation; they do not establish that
feature caching is the solution or permit subtracting debug X11 latency from
optimized embedded-frame timings to assign a runtime cause.

These are `.12` usability goals, not claims that the two narrow fixes resolve
the reported general lag.
