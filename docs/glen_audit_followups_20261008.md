# Glen Audit Follow-Ups

## Scope And Authority

Owner-authorized 2026-10-08 follow-up to
[Glen's retained report](glen_tutorial_parity_gui_audit_20261008.md), originally
supplied as `REPORT (2).md`. Its exact SHA-256 is
`186e75a987fb50cb50ed883fa2b6fa5434718fcb553ceec7aad90a8af7c1cfce`.
The reported candidate is `30f23aa4cb084694870a139714d0f0fef1726378`;
its committed `Cargo.lock` SHA-256, independently calculated in this checkout,
is `3dfe44c7a08bf779e32289cf1312103ea1fb5883e44bba5f6e1912ce3ed77745`.

Work on `codex/glen-audit-followups-20261008`. Separate narrowly scoped commits,
branch pushes and non-publishing CI are authorized. No local Rust build, test or
check, Claude consultation, merge into main, tag change or release publication
is authorized. Leave `paper`, `output/` and historical evidence unchanged.
Human scientific approval and an external auditor's performance verdict remain
independent of code, replay, screenshot and packaging checks.

## Six-Item Checklist

An authored contract is not an executed pass. Every completion needs an exact
source SHA, platform/profile, command or run URL and retained evidence. New
candidate results never overwrite or inherit the older candidate's verdict.

The owner-requested rebase is now on local main
`b1547544980cb423dbda9ee57dfce81f8b193f5b`. The isolated codex worktree retained
all 34 audit commits and the subsequent status commit. Changelog and roadmap
resolutions preserve both histories and pending gates. The overlapping pair
execution keeps main's existing guarded rollback helper unchanged, with both
main/audit test families; the mixed-SNV/indel guard is unchanged. Evidence at
`6d8b4db7` stays accepted only for that
historical source. The current candidate's checklist is:

The new frozen runtime/workflow candidate is
`3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da`, with the unchanged lockfile above.
It includes separate container-identity (`57cc62bc`) and eleven-line primer
initialization (`3ac4a02d`) repairs. The remote codex branch is deliberately
held at that SHA; subsequent local evidence/status commits are not new runtime
acceptance. Fresh non-publishing runs are
[Linux/Windows plus focused audits 37848423833](https://github.com/smoe/gentle_rs/actions/runs/37848423833)
and [macOS/Windows 37848434104](https://github.com/smoe/gentle_rs/actions/runs/37848434104).
Both run metadata independently report the exact candidate SHA. After the
fresh focused Rust and native replay steps passed, dispatch
[native packages 37850876088](https://github.com/smoe/gentle_rs/actions/runs/37850876088)
and [build/load-only container 37850885456](https://github.com/smoe/gentle_rs/actions/runs/37850885456)
with the same explicit candidate and `publish=false`; all builds remain
pending acceptance. Do not reuse the successful historical `6d8b4db7`
native-package run as current acceptance.

At this frozen SHA, all 157 local offline Python tests pass on macOS (23.617 s),
log SHA-256 `319f6f43384043334ae0383aea910176fdfa35e4bab3396487937c33fb47f3b1`.
The remote policy job passes all 180 Python tests, raw log SHA-256
`dd07dd8952e28576bb078a282dab049471b0db4e9789c202866fb4c8ab50dd52`.
The headless CLI/MCP Rust check also passes on Linux, raw log SHA-256
`6745dfd679f944370d2107c98bb10ef20e2434d97afc4d1426ccd6f2cea5414a`.
Neither result accepts native tests, binary execution, extracted packages,
container identities or scientific behavior. Scoped formatting passes and
session-close reports five OK, one plan-fidelity reminder and zero failures;
the necessary container/compile repairs are explicitly scoped above.

| Item | Current Status | Required Evidence |
| --- | --- | --- |
| 1. Audit preservation | Complete; rebase does not change original bytes | [Original archive receipt](audits/glen_followups_20261008/original_archive_30f23aa4/README.md), pinned report/archive/manifest and all 721 safe file hashes. |
| 2. Tutorial 08.04 GUI/agent contracts | Complete: bounded native contract at `3ac4a02d` | [Fresh retained receipts](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#synthetic-native-contract): all 115 hashes and safe exact membership, 19 focused Rust tests and eight ordinary Linux/X11 dev input steps. Explicit T, completion transition and all three full oracle `Seq` records verified. Synthetic 20-bp guard, not public-locus GUI or scientific approval. |
| 3. Live-agent acceptance | Fresh CLI verified; actual acceptance pending | [Unchanged macOS/arm64 dev build receipt](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#cli-build-boundary), two artifact hashes and actual full version verified at `3ac4a02d`. Require fresh public starter, genuine clarification and unexecuted explicit-T suggestions, unchanged requests/refusal and full-record manual-execution parity. |
| 4. Honest visual evidence | Post-rebase verification pending | New source-bound unedited synthetic native captures and labelled public base projection, independently hashed and visually inspected. No public-locus GUI or scientific claim. |
| 5. Scroll benchmarks | Complete: deterministic CPU cases at `3ac4a02d` | [Fresh scroll verification](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#scroll-cpu-cases), Ubuntu 24.04/x86-64 dev. All 25 retained hashes including executable, all 24 bidirectional fixture/size cases, bounded viewport and unchanged content verified against the raw log. No timed regression or native-latency claim. |
| 6. .12 package acceptance | Exact `3ac4a02d` build-only gates running | Native run `37850876088` and container run `37850885456`, both explicit `publish=false`. Require all three extracted native packages, five native/three container entrypoints and full Linux/macOS/Windows gates at the frozen SHA/lockfile/`package-opt1` recipe. Actual embedded container revisions/hashes required; no publication or old verdict inheritance. |

On 2026-10-09, independent native/scroll checks pass without local Rust
execution. Four unchanged synthetic PNGs are visually inspected: T and pair
action are visible; layered windows and light chrome/dark canvases are retained.
The public projection and genuine agent request remain separate pending gates.
The new five receipts have scoped LF attributes and individual unprotected-CRLF
negative controls. The initial three-receipt seven-suite local check passes
158 tests (24.539 s), log SHA-256
`8a64b289c669388a5d510118e5ebcb1039583993876e90e9a1d0addd8ca1a33d`;
the expanded five-receipt checkout suite passes 26/26 (21.640 s), log SHA-256
`fddda1416a55309164d5ad2e57370ff9ee798728a66273156e6edcbe79813b73`.
Session-close reports four OK, two warnings (this intentional pending commit
and manual plan fidelity), zero failures; roadmap remains 400 lines.
Local evidence commits stay unpushed while the package source is frozen.

The post-rebase `1e12fd0ff70c57413dffe91d92a7a85d92a2f7b7` policy job
([37846506619](https://github.com/smoe/gentle_rs/actions/runs/37846506619))
passed 172 Python tests. Its raw log SHA-256 is
`acab72653732a960fdb8f67e8af3ca4b145c61521238e184f80cc70f18ed55fa`.
It does not accept Rust/native/package behavior. Subsequent container review
found that `.git` is excluded and no source argument reached `build.rs`;
help-only smoke and image labels cannot verify embedded binary identity.
The scoped repair passes the frozen SHA and checks all three actual versions
and digests in the loaded image before writing its receipt. Eight new offline
regressions pass with the existing 15 container tests. A new exact-SHA freeze
and fresh gates are required; `1e12fd0f` and `6d8b4db7` verdicts are not inherited.
On `1e12fd0f` plus this scoped container/status diff, all 157 tests in the seven
offline policy/checkout/GUI-helper suites pass on macOS (23.661 s); log SHA-256
`0633fd497ba2b8442441066f79f55d0b24334b4304829dc19f5bc7fff818359e`.
This is not actual container or Rust acceptance. Fresh `1e12fd0f` CI failed
compilation in the unchanged main primer handoff family (missing filesystem
import and result-field initializers); address those separately before freeze.
Both fresh CI runs
([Linux 37846506619](https://github.com/smoe/gentle_rs/actions/runs/37846506619),
[macOS 37846515094](https://github.com/smoe/gentle_rs/actions/runs/37846515094))
failed at `1e12fd0f`; the retained Linux/Windows/headless failure log SHA-256 is
`90d460742a7418ab0237febf7ded1f90afc6ab2e0fb704576ba3839e61ab166b`.
The separate eleven-line repair adds the missing `std::fs` import and optional
`None` field in all ten production/test literals. Scoped formatting and diff
checks pass; existing deterministic regressions and remote Rust checks/tests
still need execution. Public preparation's missing CLI is a build cascade,
not an independent public-data or scientific failure.

### Historical `6d8b4db7` Snapshot

The following snapshot predates the requested rebase. Its immutable receipts,
verdicts and run links are historical evidence, not current-candidate passes or
live workflow status; they must not be overwritten or inherited.

| Item | Status | Evidence And Completion Boundary |
| --- | --- | --- |
| 1. Audit preservation | Complete: original bytes independently verified | [New integrity receipt](audits/glen_followups_20261008/original_archive_30f23aa4/README.md) pins the unchanged `30f23aa4` archive, manifest and report. Exact 319,263,072-byte archive digest, `zstd -t`, all 721 regular-file hashes and safe exact membership pass on macOS. Original project, raw images and both Criterion trees inspected; no historical binaries executed or performance/scientific verdict inferred. |
| 2. Tutorial 08.04 GUI/agent contracts | Complete: bounded native contract at `6d8b4db7` | [New native receipt](audits/glen_followups_20261008/candidate_6d8b4db7/README.md) independently verified: all 115 artifact hashes, 18 focused Rust tests, locked GUI-support/generation checks and eight ordinary Linux/X11 dev-profile input steps. Explicit T and vertically reachable pair action inspected; all three full `Seq` records match the independent oracle. This is an offline synthetic 20-base guard, not public-locus GUI or scientific acceptance. |
| 3. Live-agent acceptance | Complete: fresh public starter and real agent at `6d8b4db7` | [Executed receipt](audits/glen_followups_20261008/candidate_6d8b4db7/README.md#genuine-live-agent) `aa5a3d8e...`, verified macOS/arm64 dev CLI and existing authorized login. Genuine clarification followed by exactly two unexecuted `ask` suggestions; requests and bare-alternate refusal preserve raw project bytes. Reviewed manual execution preserves both sources and all four full `Seq` records equal both fresh CI routes, including ordered dbSNP provenance. The older-starter pilot fails full-record parity and remains historical, not accepted. No credentials inspected/copied or particular model ID claimed. |
| 4. Honest visual evidence | Complete: raw native inspection and labelled public projection at `6d8b4db7` | [Source-bound evidence](audits/glen_followups_20261008/candidate_6d8b4db7/README.md#public-reference-and-base-projection): 115 native and 44 public hashes verified. Four unedited synthetic native captures retain light chrome/dark canvases; visible T and pair action. Fresh public SVG/JSON independently bound to the 1,089-bp exports; visually inspected 33-base window highlights sole 588 C/T, with display-only uppercase and explicit projection/non-claims. Not human-locus GUI, reporter-map or scientific proof. Historical WIP previews stay untouched. |
| 5. Scroll benchmarks | Complete: deterministic CPU cases at `6d8b4db7` | [New same-candidate receipt](audits/glen_followups_20261008/candidate_6d8b4db7/README.md), Ubuntu 24.04/x86-64, Rust 1.99.0, dev profile. All 25 artifact hashes, including the downloaded executable, and all 24 bidirectional fixture/size/content checks independently verified. The [older `d1bf3fd3` receipt](audits/glen_followups_20261008/README.md) stays unchanged. No timed regression, native input-to-paint latency or package acceptance claim. |
| 6. .12 package acceptance | Frozen `6d8b4db7` builds and full platform gates running | After focused handoff/native tests pass, dispatch [native packages 37833274999](https://github.com/smoe/gentle_rs/actions/runs/37833274999) and [build/load-only container 37833283228](https://github.com/smoe/gentle_rs/actions/runs/37833283228) with `publish=false`. Both candidate receipts independently agree on source/workflow `6d8b4db7`, unchanged lockfile, `package-opt1` and `validate_only`; digest `b7afbf01b953a65a19230b7a89a72f5d06ef868c55569d08c8a63a817c04c50a`. Require all three extracted native packages, five native/three container entrypoints and full Linux/macOS/Windows gates. No older verdict inheritance or publication; native GUI/scientific approval remain separate. |

## Retained Baseline

Glen reports at `30f23aa4` on Linux/X11: checkout tests 18/18, CLI/MCP
walkthroughs 9/9, splicing protocol 10, renderer 9 plus one ignored manual
smoke, and root splicing regressions 108 serially. Catalog/manifest/workflow
counts are 61/29/51. The 08.04 online-cache-prepared workflow took 32.18 s and
produced 11 sequences, nine containers and no arrangements. These are reported
Linux results, not fresh Codex execution, Windows/macOS acceptance or validation
of the stripped binaries against the unstripped package-opt1 recipe.

Glen's GUI resize times measure window-manager acknowledgement, and wheel
dispatch is an input smoke. The TP73 first-frame CPU estimates of 50.001-59.002
ms do not establish native input-to-paint latency or a historical regression.
Do not change optimization, theme or global serialization ordering in this
follow-up merely to improve those observations.

## Execution Record

- 2026-10-08: created the persistent goal and branch at `30f23aa4`; existing
  `paper`/`output/` changes are untouched. GitHub account `smoe` is authenticated
  with repository/workflow access, and `codex login status` reports ChatGPT
  authentication. No credential contents were read.
- Requested the raw audit archive location and approved Linux/X11 runner from
  the owner; implementation and non-publishing CI work continue independently.
- Archive preparation on `30f23aa4` plus this scoped diff: original and retained
  report SHA-256 match exactly; `python3 -m unittest scripts.test_tutorial_checkouts
  -v` passed 21/21 on macOS. The new regression includes LF/CRLF checkouts and a
  failing-hash negative control without the report's LF attribute.
- Audit archive commit `5bea5ace`: the focused report LF/CRLF regression was
  rerun successfully on macOS. The scroll follow-up adds an opt-in Linux CI smoke
  with all three fixtures and 24 directional observations; Rust execution is
  still pending. Scoped `rustfmt`, YAML parsing and `git diff --check` pass.
- Scroll commit `ac18df5e5c91c04931103bfcd9907b89a98c66c6` is pushed. Non-publishing
  [CI run 37763578792](https://github.com/smoe/gentle_rs/actions/runs/37763578792)
  targets that exact SHA. The headless CLI/MCP job passed; Linux/Windows and
  the scroll smoke remain in progress. Release-policy tests rejected the new
  artifact upload's v4 action: it must use the repository-reviewed v6 action.
  The 08.04 GUI pair path exposed a reference-only partial mutation on alternate
  refusal; a scoped rollback repair and three synthetic regressions are authored,
  with execution still pending.
- No new Rust, live-model, native GUI, benchmark or package execution is yet
  claimed. Do not mark the goal complete while an item remains unverified.
- On `4f8194fc` plus the scoped semantic-control diff, Python GUI-runner tests
  pass 25/25 on macOS. Exact feature bindings and single-base policies are
  authored; they are not a native GUI pass. The `ac18df5e` scroll executable
  built on Linux, but smoke stopped at PATZ1 after 16 passing observations.
  Its retained binary SHA-256 is
  `3ed56a89a51f7998617ce51c23398221b14b9fddc95c4ae8accf67c97f84219c`.
  Investigate the visible wheel target and bind the observation's versioned
  source identity correctly before rerunning; do not weaken the movement check.
- Rerun [37766941695](https://github.com/smoe/gentle_rs/actions/runs/37766941695)
  is verified against `d1bf3fd321c651483125ee6f588f8129490e0cd5`. The headless
  and release-policy jobs passed; the clipped-map wheel repair still awaits its
  smoke result. The 08.04 guard is authored on this base with two independent
  synthetic workflows and an exact one-base persistence test. Python runner
  and focused LF/CRLF checks pass 26/26 on macOS; Rust/native/public-data and
  generated-document acceptance remain pending.
- On `9bb9e0b3` plus the scoped CI helper diff, two synthetic base-proof tests
  and all 37 release-policy tests pass on macOS. The native job compiles/tests
  only on Actions, restores historical baseline bytes in its disposable
  generation checkout, and retains source-bound projections and raw captures.
  The `d1bf3fd3` scroll job is green; its downloaded receipt still needs
  independent hash and 24-observation inspection before checklist acceptance.
- The `d1bf3fd3` Linux/Windows lib-test builds both reject the new sibling-module
  tests' access to the private pair handler (`E0624`), not a runtime allele or
  tutorial-discovery failure. Match the adjacent handlers' parent-only visibility
  before dispatching the focused native replay; do not broaden the public API.
- Verified the downloaded successful scroll artifact on macOS without executing
  or rebuilding Rust. Its receipt is retained unchanged with SHA-256
  `60add8b059510d660c36a25c1b7115138b3fbc8693dcecf2b42fd2316ddceab5`.
  All listed hashes and all 24 direction/content invariants match. Item 5's
  deterministic-case acceptance is complete at `d1bf3fd3`; timing and native
  responsiveness remain external questions, not inherited passes.
- [Native replay CI 37769190082](https://github.com/smoe/gentle_rs/actions/runs/37769190082)
  is verified against `1aee220c32f6f2d498edcc27488a423c29f7042e`. Its build and
  focused tests are in progress. The pre-existing local CLI identifies itself
  as `04debbe0eb229e9758090d1964663849621470b9`; it is not the new candidate and
  is not used to claim candidate acceptance. No local Rust build was run.
- On `93414fc89734fe70a70772a1c1c241a8661d4c4a`, all 23
  `scripts.test_tutorial_checkouts` regressions passed on macOS, including the
  scoped report, synthetic guard and scroll-receipt LF/CRLF negative controls.
  The `1aee220c` scroll job also passed; this does not substitute for inspecting
  that run's artifact or accept the still-running native replay.
- Real-agent authentication pilot on macOS used the existing
  `target/debug/gentle_cli`, self-reporting source
  `04debbe0eb229e9758090d1964663849621470b9`, binary SHA-256
  `441c608aa4325d81d5ccbcf7f7754629d63bdee61c8d3703b49bf74e3387e4e4`.
  Load the public `test_files/tp73.ncbi.gb` as `tp73_public`, then call
  `agents ask codex_local_stdio --catalog assets/agent_systems.json
  --timeout-secs 180 --max-retries 0` with an identifier/length-only prompt,
  `execution=ask`, and an explicit prohibition on mutation, fetching, nested
  agents and OS commands. `CODEX_BIN` identified the already-authenticated
  packaged Codex executable; no credentials were inspected or copied.
  The actual `external_json_stdio` response returned `tp73_public: 83,686 bp.`
  and one unexecuted suggestion:
  `introspect facts --domain project --seq-id tp73_public`. After review,
  that exact read-only GENtle command exited zero and returned the expected
  `sequence.length=83686` project fact. No command ran automatically.
  Raw files are retained in `/private/tmp/gentle-agent-auth-04debbe0`;
  `agent-result.json` SHA-256 is
  `79d26e92b9a1ce8b81860a46b07c9b94f8feb139e56d433e8d6839654bcc2d21`
  and `reviewed-command-result.json` SHA-256 is
  `a19747758f82844cd762c0cb0d4391998374e5f5388a675e3d2b5a1cad946f2d`.
  This older-binary pilot makes no build-profile, fresh 08.04, package,
  Linux-authentication or scientific-approval claim.
- The `1aee220c` native job compiled its three binaries and passed all three
  pair regressions. Its starter/oracle test rejected an incorrect uppercase
  fixture expectation: GenBank retains lowercase input, and materialization
  uppercases only the selected reference/alternate base. Correct the test and
  synthetic visual-proof helper while keeping exact source preservation,
  unchanged nonvariant bytes and the sole reference/alternate C/T difference.
  Generation and native input replay did not run after that failure.
- Repair run [37771959476](https://github.com/smoe/gentle_rs/actions/runs/37771959476)
  is verified against `9100235ce19f58f23872079e1de89ad8ec0194d5`. Native replay
  remains in progress. On that SHA, 107 Python release-candidate, desktop-package,
  container, GUI-runner and synthetic base-proof tests passed on macOS; these
  Python fixture checks do not run Rust or certify real installers.
- Add an opt-in `audit_public_vkorc1` CI preparation step to obtain the complete
  catalogued reference outside this disk-constrained checkout. The helper
  retains one raw public NCBI response and explicitly replays that exact file
  through GENtle's existing override. It compares direct/shared-shell reports,
  checks unchanged persisted state on ambiguous alternate refusal and retains
  a public starter plus the explicit C/T outputs. Complete reference input
  hashes/manifests are retained, but the multi-gigabyte cache is not uploaded.
  No new public replay, live-agent or native/public-GUI pass is yet claimed.
- On `9100235c`, all focused Rust filters and the locked GUI-support check
  passed on Linux. Tutorial generation then failed before native input replay;
  the original diagnostic artifact upload also failed because a Rust filter
  containing `::` became a nonportable filename. Sanitize only the log basename,
  keep the test filter unchanged and expose a bounded retained error excerpt.
  Do not count generation or native replay as accepted. Run
  [37774017701](https://github.com/smoe/gentle_rs/actions/runs/37774017701) is
  verified against `8d81f41714a5ce6774606d7753e4046a8ab55eb0` and includes
  public preparation; that still-running run predates this retention repair.
- Retention rerun
  [37775099176](https://github.com/smoe/gentle_rs/actions/runs/37775099176)
  is verified against `4de826d13047fb22144fe75de830d503654cae7a`.
  A public base-window SVG is authored on that base, guarded by exact raw-DNA
  hashes, one C/T difference and full source/binary identities. Seven public
  assertion/view and four synthetic/retention Python tests pass on macOS.
  Only successful public reference/parity validation may emit this labelled
  projection; it is not a native capture or an already-executed public pass.
- On `8ac39726` plus the scoped CI addition, 38 release-policy Python tests
  pass on macOS. The opt-in audit now builds a headless macOS CLI for later
  use with the host's existing authorized login. CI records full source,
  lockfile, binary, architecture and toolchain identities; no model invocation
  or credential transfer occurs there. Fresh agent execution remains pending.
- Run [37776460982](https://github.com/smoe/gentle_rs/actions/runs/37776460982)
  is verified against `f99c44a1a8a48d024630f00704ab98177dc2284d`; headless and
  release-policy jobs passed, with native/public preparation and macOS CLI
  still running. Retained diagnostics from `4de826d1` bind all four passing
  focused Rust filters and the generation failure to exact binary/lock hashes:
  `tutorial-generate` rejects the new companion's missing use-case context.
  Add that teaching metadata without weakening the validator; generation,
  projection import and native replay still need a new executed pass.
- The public preparation at `8d81f417` passed the complete-reference operation
  but failed before the slice at the helper's obsolete `/beta/` NCBI endpoint
  (HTTP 404). GENtle's existing non-beta endpoint was independently fetched on
  macOS: rs9923231, GRCh38 `NC_000016.10` C with A/G/T, raw SHA-256
  `080131663dad501faf0101fc1fd50b842b64d000d17c17a4222c654f4c16052a`.
  Raw response is `/private/tmp/gentle-public-refsnp-9923231-20261008.json`.
  Align the helper only, add an engine-default
  agreement test and fetch the retained response before costly preparation.
  The source-bound macOS dev CLI from `f99c44a1` independently matches its
  receipt and self-reported full SHA; binary SHA-256 is
  `b09b92a7c8353fb4c798922e1a35e2ea34767fcf22ea79e2fac62e36a9403da6`.
  It has not yet run a genuine public 08.04 model request or package acceptance.
- Public repair run
  [37778513095](https://github.com/smoe/gentle_rs/actions/runs/37778513095) is
  verified against `a1fb305a9d3c2dcdad1f75d732340fb0aa08ea9d`. The superseded
  metadata-only rerun `37778064938` was cancelled to avoid duplicate builds;
  cancellation is not acceptance. On `a1fb305a`, all 23 fast LF/CRLF checkout
  regressions pass on macOS without Rust execution.
- On `a1fb305a`, complete public reference preparation, direct/shared-shell
  reports, ambiguous refusal and the explicit C/T pair pass in run `37778513095`.
  Generation checks also pass, but GENtle exits 101 before `open_fragment` in
  native replay. Root-owned mode-private projects then prevent both receipt
  hashing and artifact upload, hiding the panic detail. Drop to the runner's
  UID/GID inside the isolated network namespace and upload native diagnostics
  before public preparation; do not relabel native startup as accepted.
  Run `37776460982` was cancelled after retaining its passing macOS CLI;
  that cancellation does not invalidate its recorded synthetic pilot or
  constitute overall CI acceptance.
- Genuine synthetic agent pilot on the independently verified `f99c44a1`
  macOS CLI used `codex_local_stdio`, the existing authorized host login and
  `execution=ask`, with no auto/execute-all flags or credential inspection.
  Input is the unchanged hand-crafted 20-base guard, SHA-256
  `4eed8b3c4977eb5369a591366d993753cc61b5c2cb00dd382dd4f40e4d100394`.
  The first real response asks which A/G/T alternate is intended and suggests
  no mutation. The explicit-T response requests command syntax instead of
  inventing it. After receiving this binary's `help variant materialize-allele`,
  a third response proposes exactly two parser-valid `ask` commands with the
  named outputs and explicit T; both remain unexecuted until Codex review.
  Bare-alternate engine refusal then leaves the saved state byte-identical.
  Reviewed execution creates two 20-base inserts with one C/T difference at
  zero-based position 6; all other raw bases and the source record are unchanged.
  Raw responses/projects are `/private/tmp/gentle-agent-pair-f99c44a1`;
  ambiguity response SHA-256 is
  `c151b9e71f4bc3b612a7dea8c423e0f56cd521b7e07fd84f9165caefcdc768d7`,
  initial explicit-T response is
  `4b63a1809662b20f06c76733353cd631dd707d08f865a86c64657e3e0403a9a6`,
  help-grounded response is
  `51661027985910dcfa9c1bf15e54f1a3c1104dd03a7ab13d8101382f9aa72f0a`
  and final project is
  `76acd4aecc4713a45afc3069193f6a39422959dae589c74951b62a7e8062e4c5`.
  This is real model-planning evidence on synthetic DNA, not public 08.04,
  final package, Linux authentication or human scientific approval.
- Native retention repair run
  [37782548046](https://github.com/smoe/gentle_rs/actions/runs/37782548046)
  is verified against `0cba2e9c79d5563bd40f7fb8b65fa282ec98708c`.
  It omits repeat public preparation, builds only on Actions, and keeps native
  acceptance pending. All 44 scoped helper/release-policy Python tests pass on
  macOS; no local Rust build/check/test was run.
- Retained native failure at `0cba2e9c`: independently verified all 50 artifact
  hashes and receipt SHA-256
  `2fe4448749634f5324f6f00f3f9682e676add4861e3e064b0a84962cfa22e2c7`.
  The GUI stderr identifies `xkbcommon-dl` failing to load
  `libxkbcommon-x11.so` before the first `open_fragment` input; this is an
  audit-runner dependency failure, not an allele result. The three binaries,
  all focused Rust filters and four generator checks pass on Linux/dev at
  that SHA. Add `libxkbcommon-x11-dev` only to the native audit dependencies;
  native captures and tutorial-input acceptance still require a fresh run.
- Import the five required projections retained by the passing `0cba2e9c`
  generator checks. The companion's `executed` flag refers to its engine
  workflow, not the failed native replay. Restore only the ledger's original
  Unicode/LF presentation; historical PATZ1 bytes and checksum values stay
  unchanged. Test the real generated chapter and hub hashes in both LF/CRLF
  checkouts, including an unprotected negative control.
- At `4ca26d7a` plus the scoped projection/serializer diff, 145 Python tests
  pass on macOS across audit helpers, public assertions, GUI runner, LF/CRLF
  checkouts, release policy, native packaging and container wiring. Independently
  verify both imported Markdown hashes and all three unchanged historical
  PATZ1 report bytes/checksums against `30f23aa4`. No local Rust build/check/test
  or real installer/container execution is represented by these fixture tests.
- Independently verified the `a1fb305a` public receipt SHA-256
  `d5865df2a9940075d6edd9c3729ab68487b347301fd6ae67d8e1ea2c30f47820`
  and all 44 listed artifact hashes from run `37778513095`. The raw C/T inserts
  are 1,089 bp and differ only at zero-based position 588. Direct/shared-shell
  full sequence records and scientific report content agree. The public SVG
  renders legibly with `rsvg-convert`; it explicitly remains an exported data
  view, not a native screenshot or a full reporter-construct map. Bundle:
  `/private/tmp/gentle-vkorc1-public-a1fb305a` (also the run's public artifact).
- Genuine public-input agent pilot used this same `a1fb305a` macOS/arm64 dev
  CLI, binary SHA-256
  `f678e85c9be3a08e85ccd8422fbe3392bd11b639a81d580049ab1ea5833f9619`,
  with the authorized host login and `codex_local_stdio`. The retained starter
  SHA-256 is `35ee53107482d13cb5291c8408e3c220294fb5d797cf6e52afb3df1b839a58f5`.
  The model asks which A/G/T alternate to use; after actual GENtle help and
  the reviewed genomic-forward T choice, it proposes exactly two parser-valid
  `ask` commands. Both are unexecuted until Codex review. Bare-alternate refusal
  leaves persisted bytes unchanged; reviewed execution preserves both source
  records and produces full insert records identical to the CI route, including
  the sole C/T difference at position 588. Response SHA-256 values are
  `13a1683010c27730f92b185326b0efae759a576fa4671897c4e04b8c562b3c35`
  (clarification) and
  `d215a897763a64e3d2a861c8247db771c15bde90d9bc39687a7703a45e5be8bc`
  (reviewed-T suggestions); final project SHA-256 is
  `a7a317131d14902400b79abbe27ae2a8916f578948a8a1c04fc236cb9e5cfda7`.
  Raw evidence is `/private/tmp/gentle-agent-public-a1fb305a`. This accepts the
  bounded public-input live-model boundary at this SHA, not the final package,
  native GUI, clinical effect, wet-lab suitability or human scientific approval.
- Candidate `62f7180b40471bd6d33e71b64f3fab2a98066efa` is frozen for
  [Linux/native/public CI 37787662254](https://github.com/smoe/gentle_rs/actions/runs/37787662254),
  [macOS CI 37787672955](https://github.com/smoe/gentle_rs/actions/runs/37787672955),
  [build-only installers 37787687283](https://github.com/smoe/gentle_rs/actions/runs/37787687283)
  and [build/load-only container 37787699458](https://github.com/smoe/gentle_rs/actions/runs/37787699458).
  Installer/container candidate receipts independently agree on this source,
  workflow revision, unchanged lockfile, `package-opt1`, `publish=false` and
  `mode=validate_only`; receipt SHA-256 is
  `3d436f866a4f83d7d13a351444bf7320c04ef316f882b1d25fd20b0e0f86dc9d`.
  Candidate resolution is not extracted-package acceptance.
- At `62f7180b`, six native steps pass through explicit T, then `scroll_to_pair`
  times out: `promoter.materialize_pair` remains enabled but horizontally
  clipped. Retained raw capture and semantic snapshot identify the off-screen
  action row; neither oracle mutation nor a relaxed verifier is an acceptable
  substitute. Receipt SHA-256 is
  `d9bcb400b9b52e3fc19da24492ab176544f55dac30b47f2584ef1c3f40aa51c6`.
  The upload has 91 independently verified hashes but omitted the hashed
  `native/vkorc1_warfarin_promoter_luciferase_gui/profile/home/.gentle_gui_settings.json`.
  Include hidden files only in the isolated audit evidence directory and rerun;
  preserve this incomplete upload as failed evidence.
- On `6f5f2365` plus the scoped layout diff, move only the pair action out of
  the wide toolbar into its own left-aligned row. A synthetic egui regression
  sends only vertical wheel events and requires the real semantic button to
  become visible/enabled without mutation. The native contract and independent
  starter/oracle remain unchanged. New source-bound CI/package execution is
  required; no local Rust build/check/test is run.
- Owner identified two additional Downloads documents. Read
  `GENTLE_DOCS_AUDIT_2026-10-07.md` (an older documentation audit at `04debbe0`)
  and `GENTLE_SMOE_BOT_BRANCH_HANDOFF_2026-10-08.md` (branch inventory at
  `a02ef777`), as well as the unchanged `REPORT (2).md`. Neither additional
  document supplies the original `30f23aa4` raw bundle or a remote-runner access
  route. Their cleanup/merge suggestions are out of this six-item scope and
  are not executed. Fresh Actions evidence does not verify Glen's raw archive.
- Linux and macOS workflow-example runtime filters at `62f7180b` each pass
  94/95 tests, failing `tss_tutorial_agent_drafts_parse_and_keep_mutations_reviewed`
  because its four-tutorial expectation predates the new 08.04 contract. The
  retained `tutorial-check.log` independently reports five tutorials, 27 cases,
  eight declared and 15 parser-classified mutations, with no parity findings.
  Update only these exact counts; parser validation, mutation review and the
  08.04 ambiguity refusal remain unchanged. New Rust execution is pending.
- On `1f5610fd` plus this count-only repair, 65 Python checkout, GUI-runner and
  public-audit helper regressions pass on macOS; scoped `rustfmt --check` and
  `git diff --check` pass. These checks do not execute the repaired Rust test.
- The `1f5610fd` focused GUI build fails with `E0063`: the new visibility test
  omits egui's required wheel `phase`. Its native step is skipped, not passed.
  Set `TouchPhase::Move`, matching the maintained benchmark and adjacent GUI
  tests; retain vertical-only input and all visibility/state assertions.
  On `91038fe0` plus this one-field repair, 34 Python helper/GUI-runner tests,
  scoped `rustfmt --check` and `git diff --check` pass on macOS; new Rust
  execution remains pending.
- Glen supplies the immutable original archive location:
  `/mnt/storage-box-1-gentle/Glen/handoffs/gentle/30f23aa4/tutorial-parity-30f23aa4-20261008.tar.zst`,
  SHA-256 `d257b98b95359c91d4639eb18d8a11008575d306fe7c8e890189721086437fe1`,
  with a sibling `.sha256` file. He reports 305 MB compressed,
  1,406,382,080 bytes unpacked and successful `zstd -t`. The original host-local
  tree is `/home/clawbio/.openclaw/workspace/artifacts/tutorial-parity-30f23aa4-20261008`.
  These are supplied archive facts, not independent Codex verification.
  Request an approved read-only SSH/SFTP alias or access-controlled URL and a
  small file/hash inventory; the owner need not transfer the archive on a train.
- Owner supplies the approved [audit index](https://gentle.functional.domains/openclaw/canvas/gentle-30f23aa4-audit/)
  and asks to avoid large transfers. Stop the obsolete `62f7180b` installer
  download; retrieve only three metadata files, totalling 140,532 bytes.
  The [complete per-file manifest](https://gentle.functional.domains/openclaw/canvas/gentle-30f23aa4-audit/tutorial-parity-30f23aa4-20261008.files.sha256)
  is 139,612 bytes and independently hashes to
  `0f4f723115453cb1ed491654ac7671cd392116cb1f3b8c601decd31e8d128361`.
  All 721 paths are relative, traversal-free and unique. Its 11 small-inventory
  hashes agree with the owner; `REPORT.md` agrees with the unchanged committed
  report. It lists 208 `criterion-bench-audit` and 312 `criterion-patz1` files,
  not an independent verification of those file contents.
  The 756-byte bundle manifest hashes to
  `dd0998fafaae4fc86dac0e6edbf16b98cc2d5e1e29d42754502875c3bc3903fa`;
  its source, metadata sizes and metadata hashes agree with the retrieved bytes.
  It reports the unchanged archive as 319,263,072 bytes with the supplied
  `d257b98b...` digest. The 164-byte checksum file independently hashes to
  `afcf0917ef17ee8696d19d834e328c37d79b6c03fbdea65949f9dfcbce29edd7`
  and records that same archive digest and Storage Box path. Retain the small
  files in `/private/tmp/gentle-original-audit-30f23aa4`; do not execute anything
  from untrusted metadata. No archive, binary, screenshot or Criterion payload
  was fetched, and no independent `zstd -t` or raw-content pass is claimed.
- Reverify retained `62f7180b` evidence without further transfers: all 44 public
  artifact hashes and 24 directional/content scroll observations agree with
  their receipts. Genuine agent responses remain review-first, and reviewed
  execution preserves both public source records; all four full sequence
  records agree with that SHA's fresh public CI outputs, with one C/T difference
  at zero-based position 588 in the 1,089-bp inserts. Older-candidate results
  never certify either a nearby source SHA or the rebased main integration.
- Pre-rebase candidate `b067ae3be024d2c2b7e4ea1ee841db8461211806` runs
  [Linux/native/public CI 37812324779](https://github.com/smoe/gentle_rs/actions/runs/37812324779)
  and [macOS CI 37812338675](https://github.com/smoe/gentle_rs/actions/runs/37812338675).
  Its focused native filters, locked GUI-support check, generation and ordinary
  native input replay pass, as do headless, policy, scroll and macOS CLI-build
  jobs. Raw evidence is not independently inspected. Sizes were checked without
  downloading: native 4,655,065 bytes, macOS CLI 79,759,483 bytes, scroll bundle
  88,700,925 bytes. Exact-SHA [installers 37814554493](https://github.com/smoe/gentle_rs/actions/runs/37814554493)
  and [container 37814563508](https://github.com/smoe/gentle_rs/actions/runs/37814563508)
  were dispatched with `.12` as a label and `publish=false`, without creating
  a tag or downloading packages. The subsequent rebase supersedes them.
- Owner requests finishing the interrupted rebase onto local main
  `84cc8e1654967edad0b1095672dfd4ff549b3990`. Replay all 29 audit commits;
  resolve changelog/roadmap by retaining both work streams, and GUI docs/code
  by retaining main's pure allele preflight plus the existing rollback pair.
  Main's dbSNP evidence checks, local/source cDNA orientation and boundary cache
  regressions remain intact. Range-diff confirms no unrelated audit patch
  changes; the already-upstream handler visibility needs no second code change.
  Rebased HEAD before this CI-coverage/status addition is
  `89a3e821df839e7169db570cedb9fec83c413a1e`. Its 103 Python checkout, GUI-runner,
  helper, public-audit and release-policy tests pass on macOS, together with
  scoped `rustfmt --check`; lock/report hashes and historical generated bytes
  remain unchanged. No local Rust build/check/test is run.
- Add main's `variant_followup_allele_pair_` retry/preflight regression to the
  native audit's existing `promoter_pair_` filter family and guard both in the
  fast Python policy suite. Freeze and execute a fresh rebased SHA before any
  new GUI, live-agent or package acceptance; a rewritten commit does not inherit
  its original execution verdict. Push only the codex branch with an exact
  remote-SHA lease; leave main, `paper`, `output/` and historical evidence alone.
- Candidate `b13f2b652c7b966ba83a52b34605ddecdf3f1252` runs
  [Linux/native/public CI 37816644255](https://github.com/smoe/gentle_rs/actions/runs/37816644255)
  and [macOS/Windows CI 37816649864](https://github.com/smoe/gentle_rs/actions/runs/37816649864).
  The focused GUI/native/public, scroll, headless, macOS CLI-build and policy
  jobs pass; raw native/public evidence is not independently inspected.
  Full Linux and Windows suites fail at
  `reporter_construct_handoff_binds_unique_loaded_alternates_or_requires_review`,
  with `key must be a string` at `reporter_ops.rs:2156`: raw JSON in generated
  `op` commands loses its quotes in the unchanged shared tokenizer.
  Linux reports 4,112 passes/one failure; the companion Windows job reports
  4,052 passes/one failure. Retain these verdicts at their exact SHA rather than
  dispatching known-red packages. Quote all four handoff operation payloads
  with `quote_shell_arg`, add a synthetic actual-parser/executor regression
  with quoted identifiers and a temporary FASTA path, and include the handoff
  family in focused native CI. No local Rust execution or parser relaxation.
- On `b13f2b65` plus the scoped handoff repair, all 104 focused Python checkout,
  GUI-runner, audit-helper and release-policy tests pass on macOS. Scoped
  `rustfmt --check`, YAML parsing and `git diff --check` pass. The new Rust
  parser/executor regression is authored, not executed locally; acceptance
  requires the newly committed source on Actions.
- After the owner allowlisted the hotel's IP and restored normal connectivity,
  retrieve the unchanged approved archive. Verify its pinned size/digest,
  `zstd -t`, all 721 regular-file hashes, 313 directories and exact safe manifest
  membership without extraction or execution. Retain the new receipt with
  SHA-256 `1a41c227b6e91587e6d310114d6672d875008f51170bf75428da303b3b50dbb5`.
  Separately inspect 82 bounded data files, original native/export images and
  Criterion slope records. Item 1 is complete for preservation only; original
  Linux execution, source/profile authenticity and performance remain reported.
  This evidence-only commit does not change the frozen runtime/package
  candidate `6d8b4db772a5ed08f5fd7926ec4f1635a14a236a` or inherit its verdicts.
  On that candidate plus this evidence-only diff, all 24 Python checkout tests
  pass on macOS, including the new LF/CRLF missing-attribute negative control;
  `cmp` confirms receipt byte identity and `git diff --check` passes. No local
  Rust build/check/test was run.
- Frozen `6d8b4db7` [native/public CI 37830325283](https://github.com/smoe/gentle_rs/actions/runs/37830325283)
  passes all 18 focused Rust tests, including the new quoted-handoff actual
  parser/executor case, both main/audit pair families and semantic/starter
  checks. Independently verify all 115 native artifact hashes, eight ordinary
  input steps, completion transition and three full oracle `Seq` records.
  Inspect four raw captures without edits; explicit T and the pair action are
  visible. Retain receipt `c0b0c89a...` unchanged. Item 2's bounded synthetic
  offline GUI contract is complete, not a public-locus or scientific GUI pass.
- At that same SHA, independently verify the fresh scroll executable and all
  25 listed artifacts, plus all 24 direction, viewport-bound and unchanged-DNA
  checks. Receipt `420e2c28...` remains dev-profile CPU evidence, not timed
  performance or package-opt1 acceptance. Retain both new receipts with scoped
  LF attributes and missing-attribute LF/CRLF negative controls; older receipts
  and historical evidence remain untouched.
- Independently download/hash the `6d8b4db7` macOS/arm64 headless dev CLI and
  check its actual `--version`: binary SHA-256
  `e6129b249aa061837829df7b7f1f591f9a47e57b2fb91d5c2c18c314b736e120`.
  Genuine `codex_local_stdio` requests retain clarification response
  `2125c304e31eb8e27aeae5c55ea780384e0d972e5db06e88e81dd1840135a9c3`
  and explicit-T/help response
  `04bcfafd13b620da359a6f1484630add7de2d9c809cea5a049b8dd2b74d663f7`.
  Both leave the public starter's bytes unchanged; no automatic execution,
  nested agent, fetch, OS command or credential inspection. A bare alternate
  is refused with exit 1 and unchanged bytes. Manually reviewed execution
  creates the 1,089-bp pair with only 588 C/T and preserves both source records.
  Final project hash is
  `eb83a93c434766b0719bc0cb102e2bd8981bdb085d41330a22dfad9dcc5b134d`;
  raw pilot is `/private/tmp/gentle-agent-public-6d8b4db7`. Its input is the
  hash-verified retained `62f7180b` public starter, not a new download. Full
  sequence-record comparison with fresh same-SHA public CI remains required.
- Freeze runtime/package SHA `6d8b4db772a5ed08f5fd7926ec4f1635a14a236a`
  for full [Linux/native/public CI 37830325283](https://github.com/smoe/gentle_rs/actions/runs/37830325283)
  and [macOS/Windows CI 37830335706](https://github.com/smoe/gentle_rs/actions/runs/37830335706).
  Dispatch installer/container runs only with `publish=false` after focused
  native/handoff success; their identical `b7afbf01...` candidate receipts are
  identity resolution, not extracted-package acceptance. All full platform
  and package verdicts remain pending; no tag or publication is performed.
- On archive-evidence commit `d17d830d` plus this native/scroll retention diff,
  all 25 Python checkout tests pass on macOS. Both new receipts exactly match
  their downloaded raw bytes; their missing-LF-attribute negative controls and
  LF/CRLF replays pass. `git diff --check` is clean; no local Rust execution.
- On committed evidence HEAD `6fcbffc769f55165978800448961acf9e88c8526`,
  all 149 focused Python checkout, GUI-runner, audit-helper, public-reference,
  release-candidate, native-package and container tests pass on macOS.
  Retained log SHA-256 is
  `dab18e48959a8e6d22a197c5554ece452da50054086ef710393b0b4cf55eb624`.
  The owner switched the shared checkout to separately advanced main
  `e04f6044`; continue the codex branch in the isolated managed worktree
  `/Users/u005069/.codex/worktrees/glen-followup-evidence/gentle_rs` without
  modifying main, its `paper`/`output/` changes or the frozen runtime SHA.
- At frozen `6d8b4db7`, public CI passes complete reference preparation,
  direct/shared reports and exact insert records. Independently verify all 44
  downloaded artifacts and retain unchanged receipt `d65c02e6...`, JSON and SVG.
  Render and inspect the unedited 920-by-430 SVG preview: three 33-base rows,
  source/binary hashes, sole 588 C/T, display-only uppercase and explicit
  exported-base projection/non-claims. The synthetic raw GUI captures and
  human-locus base projection remain distinct evidence classes.
- Comparing the older-starter live pilot against fresh public CI finds missing
  dbSNP assembly/reference qualifiers in all four full records. Preserve that
  failed comparison and its original files, rather than sorting annotations,
  removing provenance or comparing DNA alone. Repeat genuine requests from the
  fresh starter `7eade3a6...` with the independently verified 6d macOS CLI.
  Clarification response `b890be0a...` and explicit-T response `1d6c9122...`
  leave raw project bytes unchanged; no automatic execution. Bare alternate
  returns exit 1 with unchanged bytes. Two reviewed manual commands exit zero;
  all four full ordered `Seq` records exactly equal both fresh CI routes and
  both originals are preserved. Retain executed receipt `aa5a3d8e...` separately
  from the unchanged build receipt's `live_agent_accepted:false`. Items 3 and 4
  are complete within these bounded contracts, not package or scientific
  approval. Final project SHA-256 is
  `75d9158a14ef10bf84d75e585a47b2eecb7ef8d5fb85d436bb6752125ed27534`.
- On `6fcbffc7` plus this evidence-only diff, all 25 Python checkout tests
  pass with seven individually exercised missing-LF-attribute controls and
  both LF/CRLF replays. The first test attempt exposed a temporary-directory
  name collision between the new JSON and SVG; using the full filename instead
  of its stem repairs the harness without weakening any digest check.
  Log SHA-256 is
  `336aec55899f6ac63391c4ea008b6730b20736a826dccce142d4db2ef07bfc22`.
  `git diff --check` passes; no local Rust build/check/test is run.
- On committed evidence HEAD `24bf5f4ff4d6981a4c973abb04ec5ade6227419c`,
  all 149 focused Python checkout, GUI-runner, public-reference, release-candidate,
  native-package and container tests pass on macOS; log SHA-256 is
  `84431e0fd61a1982a2551c0637b9269c9a75eb807fb8263712a264e48a24c9f0`.
- At the owner's request, rebase the isolated codex branch onto local main
  `51f31afacef0fc942f523aea16d56f32790f3dee`, including its rack and primer
  updates. Rebased HEAD before this status-only addition is
  `1d8774e42fc63b3d6741e240a06c1d7af98db5e9`. Verify all 34 commits with
  `git range-diff`: 32 unchanged patches and two changelog-only context/placement
  resolutions preserving both histories. No runtime conflict resolution or
  primary-checkout edits. Keep the safety branch
  `codex/glen-audit-before-main-e04f6044` and all original receipt bytes.
  Reset current-candidate items 2-6 to pending, preserving historical `6d8b4db7`
  verdicts. Freeze a new post-rebase SHA for fresh CI, live-agent and package
  acceptance; do not mark rewritten commits accepted from their old runs.
- On rebased `1d8774e4` plus this status-only diff, all 149 focused Python tests
  pass on macOS (50.758 s); retained log SHA-256 is
  `b6d82f10d5552a3abd345254050c510548963b51c8f864da3c1c3324685d5d86`.
  Scoped `rustfmt --edition 2024 --config skip_children=true --check` and
  `git diff --check` pass. Incoming main's rack and primer implementation files
  exactly match main, and original report, historical receipts and lockfile
  remain byte-identical. The roadmap remains 400 lines. No local Rust build,
  check or test was run; fresh Actions source-bound acceptance is still pending.
- Main advanced again to `b1547544980cb423dbda9ee57dfce81f8b193f5b` while
  the initial rebase was being verified. Include that completed merge before
  pushing or dispatching fresh CI. Rebased HEAD before this final status-only
  addition is `2783d80b029fe6778c3f5532e0583e9887c7ecbe`; all 35 picks remain.
  Preserve main's reusable `apply_variant_followup_allele_pair` helper, its
  non-panicking lock/output checks, preflight, active-view preservation and
  late-refusal/history regression. Retain all four `promoter_pair_` and both
  `variant_followup_allele_pair_` tests; native CI already selects both families.
  The mixed-SNV/indel and incoming primer/rack implementations stay unchanged.
  Main's [M1-M3/L1 verification handoff](vkorc1_allele_review_followup_20261008.md)
  remains applicable; this is conflict reconciliation, not executed acceptance.
- On rebased `2783d80b` plus this final status-only diff, all 149 fast Python
  tests pass on macOS (24.807 s); retained log SHA-256 is
  `5013e2f90b07cb55bf458aa2d569299a0dcfdf382eac185f51577040fe5fc9fc`.
  Scoped Rust formatting, whitespace and conflict-marker checks pass. Original
  report, archive/acceptance receipts and lockfile are unchanged; incoming
  rack, primer, mixed-allele and engine-test source matches `b1547544` exactly.
  The roadmap remains 400 lines; primary main still has only the owner's
  untouched `paper`/`output/` changes. No local Rust build/check/test was run.
