# Glen Audit Follow-Ups

## Scope And Authority

Owner-authorized 2026-10-08 follow-up to
[Glen's retained report](glen_tutorial_parity_gui_audit_20261008.md), originally
supplied as `REPORT (2).md`. Its exact SHA-256 is
`186e75a987fb50cb50ed883fa2b6fa5434718fcb553ceec7aad90a8af7c1cfce`.
The reported candidate is `30f23aa4cb084694870a139714d0f0fef1726378`;
its committed `Cargo.lock` SHA-256, independently calculated in this checkout,
is `3dfe44c7a08bf779e32289cf1312103ea1fb5883e44bba5f6e1912ce3ed77745`.

The original goal used `codex/glen-audit-followups-20261008`, with separate
narrowly scoped commits, branch pushes and non-publishing CI authorized. On
2026-10-09 the owner explicitly requested checkout of `main` and integration of
this work, then requested integration of `gentle_rs_2_main` and `gentle_rs_3`
and a normal push of main to GitHub. These requests supersede the original
prohibitions on merging and pushing main for those requested integrations.
On 2026-10-11 the owner requires a new explicit request before any GitHub push
and requests continued development on `main`. These latest instructions
supersede the original branch and automatic-push authority: make narrowly
scoped local commits on `main`, without pushing unpublished work or treating
an earlier push request as continuing permission.
No local Rust build, test or check, Claude consultation, tag change or release
publication is authorized. Leave `paper`, `output/` and
historical evidence unchanged.
Human scientific approval and an external auditor's performance verdict remain
independent of code, replay, screenshot and packaging checks.

## Owner-Requested Main Integration

Integrate the 45 commits from codex tip
`4c394efefddd51c83bea93525e5024383d2cf372` into local main
`b1547544980cb423dbda9ee57dfce81f8b193f5b` with an explicit merge commit.
The codex branch already includes that main base; the merge has no conflicts
and preserves both histories. The only integration-specific edits record the
changed authority and pending acceptance in this ledger, roadmap and changelog.
Existing `paper` and untracked `output/` changes are not staged or altered.
Main was not pushed during that integration; the owner subsequently authorized
the development-branch merges and push described below. Neither integration
completes the six-item goal.
The codex remote remains frozen at `e42daa02` for its running verification;
those and older source-bound results do not accept the new main merge SHA.
Fresh exact-source native/full-platform and package gates remain required.
Before committing the integration, all 181 scoped offline Python tests pass
on macOS (44.214 s), log
`/private/tmp/gentle-main-audit-merge-python-unsandboxed.log`, SHA-256
`0dc847be996a1301ddbb91e1bd1d18acf47a605cd324885582a5035cd1eb9168`.
The initial sandboxed run has four release-monitor assertion failures because
`ps` raises `PermissionError`; a synthetic-process reproduction identifies that
boundary and the full unsandboxed rerun passes without code changes. Retain its
failing log SHA-256
`78862d452d8878bfce36abb8183d9a04d3c77e051226acda22d1480152d248b4`.
Scoped rustfmt, CI/container YAML parsing and both staged/unstaged whitespace
checks pass. The merged code tree equals the codex tip; only the three
integration-status documents differ. No local Rust build/check/test was run.

The subsequent owner-requested integration starts from main `1ff8aa0d`.
`gentle_rs_2_main` tip `b1547544` is already an ancestor; merge `gentle_rs_3`
tip `4c24ce70` with both histories preserved. Resolve the changelog by retaining
both histories and the primer module by retaining both required imports and
main's stronger portable raw-path regression. A normal push of this combined
main is authorized, not a force push, tag or release publication. New
multi-reference behavior and generated catalog/parity projections need fresh
Rust-backed checks at the exact merged SHA; the frozen codex runs do not certify
this integration or close `.12` acceptance.
The incoming commit adds two `OpResult` fields without initializing its ten
existing production/test constructors. Complete only those empty defaults and
add an offline source guard before pushing. Run the 04.09/source checks in the
existing fast CI policy job; compiler verification remains CI-only.
Pre-commit offline checks on macOS pass 184/185 in 107.902 s; raw log
`/private/tmp/gentle-main-development-merge-python.log` has SHA-256
`0f63be7f161eeb83088ebb6f8d85baf5ed01bb6397aeecb694703d836a723a0d`.
The unchanged `test_runner_forwards_termination_and_fails` returns `-15`, not
`143`, both in the full run and isolated rerun. Its fixed 0.1-s startup wait
remains a separate follow-up, not a silently waived pass or a change to this
merge. All changed Rust files pass scoped edition-2024 rustfmt; YAML parsing
and whitespace pass. No Rust build/check/test runs locally. An unrelated
`Cargo.lock` dependency update appeared during verification; leave it unstaged
and out of the push. The committed lockfile and older evidence remain unchanged.

## Six-Item Checklist

An authored contract is not an executed pass. Every completion needs an exact
source SHA, platform/profile, command or run URL and retained evidence. New
candidate results never overwrite or inherit the older candidate's verdict.

### Final-Candidate Review At `f35a58e6`

Glen executed the consolidated request at exact
`f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`, tree
`f678fa505010119b0b546db116a82438e07f6405`, excluding the exercise tag.
Codex independently verifies the unchanged 103,172,261-byte archive, all 213
file hashes and exact safe membership. Preserve the
[source-bound receipt, original report and manifest](audits/glen_followups_20261008/candidate_f35a58e6/README.md).
No archived executable is run locally. New ledger commits do not become the
tested source; all preceding snapshots remain historical.

| Item | Status At `f35a58e6` | Evidence / Remaining Gate |
| --- | --- | --- |
| 1. Audit preservation | Complete; old and new bytes preserved | Original `30f23aa4` evidence stays unchanged; new report/archive/213-file manifest independently verified. |
| 2. Tutorial 08.04 GUI/agent contracts | Focused software checks pass; native replay blocked | 19 Rust tests and projections pass. X11 receipt fails at interactive `sudo unshare` policy, not a product assertion; do not bypass isolation or claim native acceptance. |
| 3. Live-agent acceptance | Not run at this source | No authorized authenticated session; earlier `3ac4a02d` pass remains historical. |
| 4. Honest visual evidence | Source/placement/export contracts pass; native captures open | Tutorial/visual/provenance, TSS export and FLNA 2400-px render evidence retained. No ordinary native readability or human scientific approval. |
| 5. Scroll benchmarks | Deterministic CPU contract complete at this source | 24 exact-SHA fixture/size/direction observations and executable hash independently checked; native latency explicitly not measured. |
| 6. .12 package acceptance | Open | Exact main CI is now fully green, including macOS, but full Linux/optional GUI jobs are skipped; no same-SHA `package-opt1` native/container/extracted-package run. |

08.16 real-DuckDB, 08.18 LocalProject, 07.04/10.03 learner and focused
04.09/parity checks pass within their documented software scope. Native views
remain untested. FLNA still lacks historical release evidence and human
teaching-exon/interpretation approval; keep its draft inactive. Real BLAST
searches succeed, but imports of 5,999/12,202 HSP rows remain performance-blocked
without a completed specificity verdict. Preserve raw outputs and collect the
bounded import-only profile described in the receipt before changing rules.
Frozen `e981469e`/`dfadf2a5` audits remain separate. No local Rust, push,
dispatch, package publication or overall six-gate completion is claimed.

### Human FLNA Decision Supplement: 2026-10-13

The [retained human decision and inspection receipt](audits/glen_followups_20261008/flna_teaching_exon_f35a58e6/README.md)
close only teaching-exon selection at the historical `f35a58e6` source.
Codex independently verifies the original 992-byte record against the owner's
SHA-256 `5202f0605fca8cdd24bbd00e5066da855ffe41ff9aa6695ad09c8f75da4972be`
and preserves the small review manifest. Steffen selects local one-based closed
`19370-19393` (`[19369,19393)`, 24 bp). Retained CDS geometry differs by that
interval alone, consistent with the annotation-derived eight-amino-acid lesson.
The original execution report and all old evidence remain unchanged.

The loaded comparison annotation is `ENST00000360319.10`, while the retained
UniProt cross-reference and decision name `.9`. Do not assert version identity
or silently rewrite the human record. Glen's separate 773-byte label addendum
is now independently verified against owner-supplied SHA-256
`7ce94e64dbb6629dc0dbe3a6baeef6c3a567e5ca75c1a5766f145fa597f730a5`;
its v2 decision is byte-identical to the original. Retain the addendum and its
small manifest separately; the version-label request is closed, not a claim
that `.9` and `.10` are interchangeable. The choice is not reopened.
Native normal-size readability, historical Ensembl
release and scientific/experimental confirmation remain open; keep the draft
inactive. Browser review/export evidence does not close gate 4 or any other
native/agent/package gate. No local Rust or new candidate acceptance follows.

### Bounded Post-BLAST Measurements At `f35a58e6`

Glen returns the diagnostic requested by coordination commit `05285a5c`, while
testing only the earlier exact `f35a58e6`. Codex independently verifies the
29,227-byte archive, all 75 file hashes and exact membership; retain the
[original report/manifest and measured boundary](audits/glen_followups_20261008/post_blast_profile_f35a58e6/README.md).
Both ACTB/TP53 imports time out at 60.08 seconds with 99% CPU, approximately
915 MiB peak RSS and no saved state/stdout/surviving process. CPU sampling is
blocked, not executed: `perf` is absent and kernel policy is unchanged.

The ten-second trace independently recounts 511 opens: 22 cDNA BLAST-index
files, four annotation opens, two TSV opens and **zero prepared FASTA/FAI
opens**. Do not conflate BLAST index files with the earlier per-hit FAI
hypothesis. Repeated resource inspection/annotation reads are observed but
their CPU share and the later import phases remain unknown. Original input
hash lists and reference inspection JSONs are identical before/after; no real
specificity verdict follows. Next is kernel-independent phase timing, not
unmeasured policy changes or permission bypass. No runtime repair, local Rust,
push, dispatch, publication or new six-gate completion is claimed.

### Current Main Diagnostic Follow-Up: Awaiting Execution

Glen's 2026-10-13 clarification confirms the completed software/learner scope
at `f35a58e6`; the earlier `5e6201a6` worktree request is not a new full task.
The [active consolidated tutorial residuals](glen_final_candidate_execution_handoff.md#remaining-native-and-human-checks)
are native 08.16 lanes, native microarray/gel navigation/captures/reopen, normal
native FLNA/missing historical release (human exon choice now closed above), and native 04.09
presentation with completed real import deferred until performance repair.
No 08.18 CLI repeat or cohort GUI is requested; closed tutorial-CI/link repairs
stay closed. New-code diagnostics remain separate from these residuals and
from historical candidate acceptance.

Both requested branches are now integrated: `gentle_rs_2_main` at `dee6033d`
through merge `a76b8418`, then `gentle_rs_3` at `9e717bb4`. The reported matrix
MCP test argument is repaired separately at `90e634f1`; its execution remains
pending. These runtime changes do not inherit `f35a58e6`'s acceptance.

The [new phase-trace procedure](glen_post_blast_import_phase_trace.md) measures
the same unchanged real inputs with opt-in, bounded stderr wall-clock spans,
without a kernel profiler or weakened validation. Diagnostic source and inline
tests are authored, not Rust-verified locally. Glen needs the exact new commit
through an owner push or bundle, plus focused tests and bounded real imports.
All six gate verdicts above remain bound to their recorded source; no additional
native, live-agent, package, scientific or performance closure follows.

### Next Package Recipe After The Internal.12 Publication

On 2026-10-13 the owner requests a bounded increase after the `.12` artifact
cycle. Main merges `gentle_rs_1_bbe0` at `e40b5f77` and prepares `package-opt2`:
only opt-level changes from 1 to 2, with `package-opt1` unchanged for comparison.
This is implementation permission, not fresh acceptance or dispatch authority.
Item 6 retains `.12`'s actual published receipts and open independent native
boundary; the new `.13` recipe needs its own native/container/extracted receipts
and [same-SHA opt1/opt2 resource/runtime evidence](dna_feature_rendering_latency_plan.md#measurement-profiles).
The other five boundaries and exact-source audit authorization remain unchanged.
No local Rust, push, package build, model call or publication occurs here.

### New-Source CI Regression At `6c87f8d1`

Codex follows existing [run 38414971911](https://github.com/smoe/gentle_rs/actions/runs/38414971911),
attempt 1, without dispatch/restart. The completed
[Linux job 115356926933](https://github.com/smoe/gentle_rs/actions/runs/38414971911/job/115356926933)
checks out exact `6c87f8d1f160b4c3856323048c1cb334f7a4f63d`, tree
`314953cd7ed44cddf80ff847581f5be4978792db`. Runner is Ubuntu 24.04.5,
image `20261004.327.1`, Rust 1.99.0; ordinary unoptimized/debuginfo test
profile, root defaults plus workspace feature unification, `CARGO_INCREMENTAL=0`
and `RUST_MIN_STACK=16777216`. The committed source lock
hash remains `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`,
not a separately CI-measured lock hash or a package receipt.

Ten focused multi-reference tests pass. The full
`cargo test --locked --workspace -j1` stops at the root library: 4,147 pass,
eight fail, 14 ignored, 1,069.56 seconds. All five trace-recorder tests and
the traced refusal/synthetic import pass within that failed suite. Do not
call the workspace green or accept later test binaries that were not reached.
Every matrix failure is `NotFound` for `synthetic_panel___inspection`: the
fixture stores its raw spaces/apostrophe ID, unlike production construction.
Preserve the original decoded Linux log outside Git at
`/private/tmp/gentle-ci-6c87f8d1-20261013/linux-job-115356926933.decoded.log`,
UTF-8 representation SHA-256
`df1ade2fa646465ae66ab048a56ac28a2a27dda7251424bcb7f01aee8511c523`.
This is not a raw HTTP/archive digest. Windows remains live at observation.

Repair only fixture construction, retain raw unusual adapter requests, and
add canonical/raw lookup plus read-only assertions. An early matrix filter on
the existing platform jobs catches this shared primary failure before lengthy
adapter/full-suite work. Rust verification of the repair remains external.
No production identifier rule, biological policy or historical verdict changes.

The owner explicitly authorizes 08.04 GUI replay, scroll CPU and macOS audit-CLI
checks at `6c87f8d1` only after normal CI passes, excluding packages, public-data
preparation, model calls, pushes, tags and publication. The CI prerequisite is
unmet; do not dispatch or retry the unchanged failing fixture. A new source
needs owner push/access, successful exact-SHA CI and explicit audit retargeting.
Keep all six final-candidate gates open to their actual evidence boundaries.

For the local repair, all 61 offline tutorial/checkout checks pass in 33.818 s;
retain `/private/tmp/gentle-matrix-fixture-repair-offline.log`, SHA-256
`0067bbe86619ec419e4c6e58962d7be1b3ecceaf85439fea2dab25e7d3f239ae`.
Scoped rustfmt, workflow-structure and whitespace checks pass. This is source/
documentation validation, not execution of the new Rust fixture or regression.

### Repaired-Source Audit Authorization At `8b6f8186`

The owner pushes exact `8b6f8186c81636491dadb446e94ca68429ba0860` and explicitly
authorizes retargeting only the 08.04 synthetic native replay, scroll CPU audit
and macOS headless audit-CLI build after normal CI at that SHA passes. Follow
[push run 38432492850](https://github.com/smoe/gentle_rs/actions/runs/38432492850),
attempt 1. At 2026-10-13 16:54 UTC it remains live on macOS and Windows;
headless and safety jobs pass, full Linux is skipped by the existing sampler.
Do not dispatch based on partial results, restart live jobs, or claim repaired
Rust acceptance before the focused/full results are observed.

**Dated update, 2026-10-13, after the final 19:05:33 UTC run update:** run
`38432492850` is now `completed/success` at the exact `8b6f8186` source above.
Windows job `115417977228` passed at 18:47:34 UTC; macOS `115418026971` passed
at 19:05:26 UTC; `CI / selected platform` (`115485159934`) passed at
19:05:32 UTC. Headless CLI/MCP and release-safety jobs also passed. Full Linux
and the three optional audit jobs were skipped, not accepted. This closes
the normal-CI prerequisite at that source only; it does not validate later
main `024ed6c0`, the `.13` rollover, independent packages or scientific results.
No fresh failure diagnosis or Rust repair is needed for this successful run;
no audit dispatch was made in this session. Preserve the earlier observation
and its source-bound checklist below as history.

The old [Windows job 115356869784](https://github.com/smoe/gentle_rs/actions/runs/38414971911/job/115356869784)
is now terminal at `6c87f8d1`: 4,086 root tests pass, the same eight matrix
fixture lookups fail, 14 are ignored, 1,340.73 seconds. Windows Server 2025,
runner image `20260925.250.1`, Rust 1.99.0; no distinct Windows-only repair is
inferred from the same `synthetic_panel___inspection` primary failure.
Retain decoded log
`/private/tmp/gentle-ci-6c87f8d1-20261013/windows-job-115356869784.decoded.log`,
UTF-8 representation SHA-256
`617fbfab037ddb06743190540f4306c9e90438620584c64f1d1ef9d822875176`.
The old whole run concludes failure; neither its partial passes nor new
documentation commits accept `8b6f8186`.

| Item | Fresh `8b6f8186` Boundary At This Observation |
| --- | --- |
| 1. Audit preservation / Linux validation | Historical bytes remain verified and separately attributed. Fresh headless/safety checks pass; no full Linux execution at this SHA yet. |
| 2. 08.04 GUI/agent contracts | Authored contracts remain; optional native replay awaits successful normal CI and authorized dispatch. No online/scientific acceptance follows. |
| 3. Genuine live agent | Not run; authorization explicitly excludes model calls. A headless CLI build is not authenticated agent acceptance. |
| 4. Honest visual evidence | Historical teaching-exon decision and separate version addendum are preserved; ordinary source-bound native readability/captures remain open. |
| 5. Scroll benchmarks | Historical 24-case deterministic CPU pass stays at `f35a58e6`; fresh exact-SHA CPU audit awaits dispatch. No native latency claim. |
| 6. .12 package acceptance | Not run; package builds stay on the owner's normal release path, outside this audit authorization. Native/container/extracted acceptance remains open. |

No package build, public-genome preparation, model call, push, tag or
publication is included. Immediately before dispatch, verify a published ref
still points to this authorized SHA, and verify the created run's actual SHA.
The workflow's existing ordinary baseline jobs may run with the audit jobs;
request Linux to obtain fresh full Linux evidence, without widening audit scope.
If the ref drifts, request a new source choice rather than silently retargeting
or creating/pushing a ref. Local FLNA evidence commit `5bd1cd37` is separate
and unpushed; it does not replace the authorized target. No dispatch or local
Rust execution has occurred at this observation.

### Completed Normal CI At `8b6f8186` And Owner Retarget

Run `38432492850`, attempt 1, is now `completed/success`, independently read
on 2026-10-13 19:16 UTC. Its exact source remains
`8b6f8186c81636491dadb446e94ca68429ba0860`, not later main. Both native
jobs record Rust 1.99.0, exact checkout/binary source stamps and successful
LF/CRLF tutorial parity. Locked full-workspace tests finish successfully,
not merely the formerly failing matrix filter:

| Platform / Job | Root Library | Full-Workspace Result Lines |
| --- | --- | --- |
| [Windows 115417977228](https://github.com/smoe/gentle_rs/actions/runs/38432492850/job/115417977228) | 4,095 pass, zero fail, 14 ignored; 1,199.58 s | 4,633 pass, zero fail, 17 ignored across 31 result lines |
| [macOS 115418026971](https://github.com/smoe/gentle_rs/actions/runs/38432492850/job/115418026971) | 4,156 pass, zero fail, 14 ignored; 988.47 s | 4,694 pass, zero fail, 17 ignored across 31 result lines |

Workspace totals exclude duplicate focused/engine steps. Quiet aggregates do
not enumerate named tests, and zero-test result lines are not coverage. The
early matrix filter executes ten tests on both platforms; multi-reference
executes nine on Windows and ten on macOS. Unix-only mocks do not certify
native Windows or real-reference BLAST/import acceptance. macOS runner is
`macos-26-arm64`, image `20260907.0351.1`, host `aarch64-apple-darwin`;
Windows uses `x86_64-pc-windows-msvc`. Ordinary dev/test profile and workspace
feature unification are not `package-opt1` or native readability evidence.
The committed lock hash remains
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`;
these normal jobs do not independently measure retained binary/lock hashes.

Retain decoded logs outside Git in `/private/tmp/gentle-ci-8b6f8186-20261013/`:

| Log / Observation | UTF-8 Representation SHA-256 |
| --- | --- |
| `windows-job-115417977228.decoded.log` | `4b12d5c9e0b175bd530932e086c5ac9e72116545b5e444213b5004a9167b495d` |
| `macos-job-115418026971.decoded.log` | `d8e9fc3bb8a98b42f24a0ccce256f3a9425d938974933354b570d8e043684b59` |
| `headless-job-115417977116.decoded.log` | `6447d691fd5d7bd9eb31c78c4206580fc812b46ed56d3fd49a79216a8e9829f9` |
| `safety-job-115417977132.decoded.log` | `bec3a45fca851fb57d33b1a91eacca3ec07a68fbf9a7004b8561f82e6deb6762` |
| `completion-observation.json` | `2ef0c0b6e8c6fd6184f12dd5789d0a10fc8bca6deecb1367f69a3b83db86df43` |

These are decoded log/API representations, not raw HTTP/archive digests.
Text retention normalized Windows CRLF to LF; the first Windows hash above
identifies that derived text, not the API-original bytes. Preserve the latter
separately as `windows-job-115417977228.api-original-complete.log`:
687,832 bytes, SHA-256
`e0785ae196ac33d042a66936f7f25746ba24f8220e3a77d5c02d28b893208c14`.
Removing only its 9,639 CR bytes reproduces the derived hash exactly. The
CLI's ANSI-output guard initially refuses an empty retrieval; the successful
GET explicitly retains escapes to a file, never to the terminal. No product
failure or required isolation-policy change follows from that observer guard.
Lossless gzip copy
`/private/tmp/gentle-windows-ci-8b6f8186-20261013.api-original.log.gz` is
95,782 bytes, SHA-256
`9804daca2290576431cbb65185786ecfefced5a0064728ffdf59532f4cd490bd`; decompression reproduces the original hash.
The separate evidence-only archive
`/private/tmp/gentle-normal-ci-8b6f8186-20261013.tar.gz` is 212,380 bytes,
SHA-256 `4c59c51c15713696ed4f809133bdf8eb0661a829363be97c133a26e2f3872beb`.
Its six regular members are exactly these five files plus `SHA256SUMS`, whose
SHA-256 is `9ee0c1eb08ee2582ad9fdb50d8e0462bb3b197e9dda351478bdb6b87a9b78b1c`.
Gzip integrity, safe membership and all five hashes pass after isolated
extraction. This is not an executable, installer or new Glen audit bundle;
the old archives are untouched.
Linux headless check and 228 safety/helper tests pass, but full Linux and all
three optional audits are skipped. No native GUI/scroll run, package or live
agent follows from this successful normal run.

Immediately before optional dispatch, remote main has advanced to
`024ed6c0520d38c4221517219d44c792084de4e5`; no published branch still pins
`8b6f8186`. The preflight refuses before any dispatch command executes.
The owner then explicitly retargets only the same three audit checks to
`024ed6c0` after its exact-SHA normal CI passes. Follow already-live
[run 38453643146](https://github.com/smoe/gentle_rs/actions/runs/38453643146);
at 19:25 UTC macOS/Windows are live, headless/safety pass and full Linux is
skipped. Old passes do not fulfill this new prerequisite. Main includes the
owner-pushed FLNA evidence and experimental sequence-design audit tooling;
their source-only tests do not become runtime/native acceptance or a `.12` gate.

| Item | Current `024ed6c0` Acceptance Boundary |
| --- | --- |
| 1. Audit preservation / Linux validation | Original audit bytes and completed `8b6f8186` logs stay separately attributed. Fresh headless/safety jobs pass; full Linux at this source remains unexecuted. |
| 2. 08.04 GUI/agent contracts | Optional exact-source replay awaits whole normal CI success and dispatch. No public-data, model or scientific verdict is implied. |
| 3. Genuine live agent | Not run; the audit authorization excludes model calls. |
| 4. Honest visual evidence | Human teaching-exon/version records remain preserved; ordinary same-source native presentation/capture gates remain open. |
| 5. Scroll benchmarks | Fresh CPU audit awaits dispatch; historical 24-case `f35a58e6` pass is not this source or native latency. |
| 6. .12 package acceptance | Open; packages stay on the owner's normal release path and outside this dispatch authority. |

Recheck the exact published ref immediately before dispatch and the created
run's source immediately after it. Use existing Linux selection for a fresh
full-Linux baseline alongside the three requested audits; the workflow also
runs its established Windows/headless/safety jobs. No further source choice,
push, tag, package build, public-genome preparation, model call or publication
is implicit. No optional dispatch or local Rust execution has occurred.

### Verified Small Handoff And Completed `024ed6c0` Audits (2026-10-14)

This update supersedes the live status below without rewriting it. The single
authorized [audit run 38476197008](https://github.com/smoe/gentle_rs/actions/runs/38476197008)
is completed/success, attempt 1, at exactly
`024ed6c0520d38c4221517219d44c792084de4e5`; final update 00:14:31 UTC.
Linux/Windows, headless/safety and all three optional audit jobs pass. No package,
public-reference preparation, model call, push or publication was dispatched.
The [new frozen evidence subset](audits/glen_followups_20261008/candidate_024ed6c0/README.md)
records verified versus producer-reported bytes, original versus log-derived
receipts, Linux counts, artifact identities and remaining boundaries.

The owner-supplied small Glen bundle has exactly 263 files, all 262 manifest
entries verified, manifest SHA-256
`4f147f98457b02186c883f18da308dc19b652ce38f724ca4fa4c1f4aa491fc90`.
No large archive is downloaded. Recompute both included-run summary tables:
12 runs and 124 measured metric groups each; excluded cache fills remain
excluded. Preserve A's two historical dev SHAs separately from B's same-`024ed6c0`
dev/opt1 comparison. B supports roughly 257/270-ms cold/warm usable-window
improvement and 0.93-0.95-s first-locus improvement on Linux virtual X11, not
native package, PATZ1/circular, opt2 or scientific acceptance.

| Item | Tested Scope / Remaining Gate |
| --- | --- |
| 1. Audit preservation / Linux validation | Original evidence stays unchanged. Small handoff membership and hashes verify; exact-source Linux workspace passes 4,694 / 0 / 17, root 4,156 / 0 / 14. Omitted Glen binaries/raw traces/projects and large archives are not independently rehashed here. A later final candidate needs fresh validation. |
| 2. 08.04 GUI/agent contracts | Synthetic X11 replay passes with network-namespace isolation, explicit T, eight steps and full three-`Seq` oracle parity. Refusal/rollback regressions pass. Glen separately reports hash-bound real-reference C/T/refusal at this SHA, but his physical-host GUI replay remains policy-blocked. New default-stack test repair is authored, not Rust-accepted. |
| 3. Genuine live agent | macOS/arm64 exact-source headless dev CLI build passes; its receipt explicitly says `live_agent_accepted:false`. No model was called. Genuine current/final-source acceptance remains open and outside this audit authorization. |
| 4. Honest visual evidence | All 116 downloaded GUI artifact hashes verify; raw T/pair screenshots inspected, synthetic base projection explicitly not a native capture. Some wide controls clip; physical/native lesson readability and scientific review remain separate. Glen's newly merged route labels historical images as orientation only. |
| 5. Scroll benchmarks | All 24 unique wheel/next-frame CPU-smoke observations are content-bound with valid direction/span. The receipt is log-derived and the 89-MB archive/binary is not rehashed here. Native compositor/input latency and timing-budget acceptance remain open; A/B startup comparisons do not replace them. |
| 6. .12 package acceptance | Published `5e6201a6` opt1 receipts remain separately attributed. The independently verified [validate-only container receipt](audits/glen_followups_20261008/container_validate_5e6201a6/README.md) establishes the non-publishing headless component at that source, not native/extracted-package acceptance. B is a same-source opt1 performance comparison, not installer acceptance. These `024ed6c0` CI jobs build no package; `.13` opt2 needs new exact-source native/container/extracted smokes and opt1/opt2 comparison. |

On 2026-10-14, independently retrieve the original 2,425-byte `.12` container
receipt from artifact `11781064566` of tag-push run `38322028567`. Its 1,154-byte
ZIP matches GitHub's artifact digest; extracted bytes, source/workflow, original
lockfile and opt1 recipe verify. CI executed all three binary identities and
network-disabled help/helper smokes; `mode=validate_only`, `pushed=false` and
skipped GHCR publication remain explicit. Image/binary hashes are CI-observed,
not independently rehashed here. Preserve this supplement separately from the
later publication run and old audit records. No native installer, new-source,
live-agent, scientific or overall six-gate completion follows.

Glen's small failure logs reproduce the original Linux/Rust 1.96.1 default-stack
abort for 11/2/7/2 selected allele/reporter tests (exit 101); 16-MiB ambient
reruns pass diagnostically. This is not contradicted by CI's 16/64-MiB override.
Local `4fdcb7a6902255af8275c6e5ee4270d193e99cc0` moves only those fixture bodies
onto an explicit test stack and adds panic/small-caller regressions plus an
intended override-free CI step. A shell `unset` alone still allows Cargo to
restore its configured override; that initial step is not default-stack proof.
The corrected runner removes `RUST_MIN_STACK` just before the test executable
starts and emits a JSON environment receipt. Normal full-suite configuration
and historical records remain untouched. Offline subprocess/wiring tests pass,
but the final source and owner-updated lockfile need fresh external Rust evidence,
including the quick-xml SVG-prefix regressions. No production/parser broadening
or local Rust verification is implied.
[Glen must rerun the original commands at the new source](vkorc1_allele_review_followup_20261008.md#2026-10-14-default-stack-follow-up).

The owner confirms the dependency update is finished and directs Glen to the
final committed tip after remaining development, not an earlier checkpoint.
The unchanged updated lockfile is frozen with SHA-256
`0b9729e60698e14d464b48230e606377a4e2e48054c2450be6f85760e999a6b9`.
Its static structure checks do not establish Rust, native, model or package
acceptance. Do not dispatch a new run, publish, or transfer old verdicts merely
because dependencies are now committed; the owner's no-push instruction remains.

The retained `024ed6c0` raw `review_t` screenshot shows the analysis/export
controls clipping to the right after the parameter grid widens scrollable
content; the T field and separate pair action remain readable. Bound only
those two action rows to the visible viewport, preserving labels, enablement,
handlers and parameter horizontal scrolling. Two authored egui regressions
cover minimum/default geometry and all 23 real-window action labels with
vertical-only input and unchanged project/history. Include their filter in the
existing opt-in audit and [final-source native check](glen_final_candidate_execution_handoff.md#remaining-native-and-human-checks).
This is a scoped implementation repair, not an executed Rust/native verdict;
item 4 remains open and original screenshots/evidence are unchanged.

### Completed Normal CI At `df5c8244` (2026-10-14)

[Run 38541824374](https://github.com/smoe/gentle_rs/actions/runs/38541824374),
attempt 1, completes successfully at 12:12:34 UTC, exact source
`df5c8244787d8c920ee171ee55167de55e8a7294`. Preserve the new
[completed-run inspection and original-log hashes](audits/glen_followups_20261008/candidate_df5c8244/README.md)
separately from earlier observations and `024ed6c0` records. macOS has
4,696 workspace passes / 0 failures / 17 ignored, Windows 4,635 / 0 / 17,
each across 31 summaries in one full-workspace command. Root counts are
4,158 / 0 / 14 and 4,097 / 0 / 14 respectively. Headless/safety also pass.

| Item | What This New Run Establishes / Remaining Gate |
| --- | --- |
| 1. Audit preservation / Linux validation | Source/lock-bound terminal logs retained; Linux headless passes, full Linux is skipped. Old audit bytes remain unchanged. Later final-source validation is still required. |
| 2. 08.04 GUI/agent contracts | Named SPDI round-trip/unfolding and allele/reporter/helper tests pass on macOS/Windows. Synthetic GUI audit is skipped. Shell-only stack unset does not prove default-stack families; the corrected runner still needs execution. |
| 3. Genuine live agent | Not run; macOS audit-CLI job is skipped. No model acceptance or authorization follows from normal CI. |
| 4. Honest visual evidence | Real LF/CRLF tutorial replay passes, but native captures/readability are not run. The optional render-snapshot step selects zero tests; it is not snapshot acceptance. Local viewport repair needs fresh Rust/native checks. |
| 5. Scroll benchmarks | Job skipped; retain the 24-case `024ed6c0` CPU verdict separately. New-source wheel/next-frame and native timing boundaries remain open. |
| 6. .12 package acceptance | No package job. Released and validate-only `5e6201a6` evidence remains separately attributed; independent extracted/native acceptance and newer opt2 checks remain open. |

The six newer local commits through `5909a5dd`, including the completed owner
lockfile update, are not tested by this run. Glen targets the final available
committed tip after remaining development. No local Rust, push, dispatch,
model call or publication was performed, and no item is newly marked complete.

Glen's documentation commit `b6e54e0a` is integrated as local `99ac1715`, preserving
the `.13`/opt2 work and unchanged historical evidence. The new manual/reference
03.09 route has separate design/Gibson fixtures, four provenance-labelled visuals
and no continuous experimental/function claim. Its boundary tests now run in CI.
The catalog has 66 entries; executable chapters remain 30. No Golden Gate/RNA/guide
implementation is inferred from the new backlog; keep future work separate.

Glen's reported no-finding scan of 15 commits / 61 files at `024ed6c0` remains
source-bound: 10 matrix and five trace tests pass, 22 Python pass / four binary
checks skip, and both exact-source CI runs succeed. It does not close real
post-BLAST import performance, native screenshots, historical FLNA release
provenance or the other independent gates. Do not re-open accepted `f35a58e6`
08.18 CLI, tutorial registration or link repairs merely to repeat old requests.

### Completed Normal CI And Authorized Audit Dispatch At `024ed6c0`

On 2026-10-13 at 22:10:56 UTC, [normal run 38453643146](https://github.com/smoe/gentle_rs/actions/runs/38453643146)
attempt 1 completes successfully at exact
`024ed6c0520d38c4221517219d44c792084de4e5`. The earlier live observation above
is retained unchanged. Rust 1.99.0 (`b940084d7`), ordinary dev/test builds,
locked workspace tests and native LF/CRLF tutorial parity pass:

| Platform / Job | Root Library | Full-Workspace Result Lines |
| --- | --- | --- |
| [Windows 115490639067](https://github.com/smoe/gentle_rs/actions/runs/38453643146/job/115490639067) | 4,095 pass, zero fail, 14 ignored; 669.15 s | 4,633 pass, zero fail, 17 ignored across 31 result lines |
| [macOS 115490688166](https://github.com/smoe/gentle_rs/actions/runs/38453643146/job/115490688166) | 4,156 pass, zero fail, 14 ignored; 1,014.39 s | 4,694 pass, zero fail, 17 ignored across 31 result lines |

Totals cover only each full-workspace command, excluding later/focused filters.
Both early matrix filters execute ten tests; multi-reference filters execute
nine on Windows and ten on macOS. Quiet counts do not attest named tests or
real-reference/native GUI acceptance. Windows image is `windows-2025-vs2026`
`20260925.250.1`, target `x86_64-pc-windows-msvc`; macOS image is
`macos-26-arm64` `20260907.0351.1`, target `aarch64-apple-darwin`. Full Linux
and all optional audits are skipped in this normal run. Linux headless checks
the locked CLI/MCP no-default-feature dependency boundary, not runtime.
Safety executes 248 helper tests: 244 pass and four real sequence-design CLI
tests skip. The committed source lock hash is
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`,
not a separately measured CI binary/lock receipt. No import speedup is inferred.

Preserve API-original log bytes, including Windows CRLF and ANSI, outside Git
in `/private/tmp/gentle-normal-ci-024ed6c0-20261014/`:

| File | SHA-256 |
| --- | --- |
| `windows-job-115490639067.api-original.log` | `9e7da94659f3d23540a671393e71b4695d972bd32429708f45249ab3e2b7fec9` |
| `macos-job-115490688166.api-original-complete.log` | `66ab6872ffb0fdfe47044ccaa6a95a240b93f03ba3ca2bee7deb2e627bbbea46` |
| `headless-job-115490639102.api-original.log` | `4e63e8bdf8fa199e787110451aeb73c7a7d9392931c829cd9cc5e52fd6fec02c` |
| `safety-job-115490639076.api-original.log` | `9e105b95f7a0443eb75af6810ab15f5f540716426bdcc5f393ac86a9b764c43b` |
| `completion-observation.json` | `96f8944bf2a84167add3a90ebe999c51cf3124c7e3782e0ccba17ce744536bc9` |

The evidence-only archive
`/private/tmp/gentle-normal-ci-024ed6c0-20261014.tar.gz` is 204,624 bytes,
SHA-256 `d6c45af141bf47304106b158401a022c7133071fb6e96944b22785638dc39fbb`.
Its exact six regular members are the five files above plus `SHA256SUMS`,
digest `fa1db7b61497e157bf9e0e47498825e2e2682c8709af23debeb3c10b54b33ae1`.
Gzip integrity, safe membership and all five raw-file hashes pass after isolated
extraction. The earlier failed live-log observation is excluded; older archives
and evidence remain unchanged. This archive contains no installer or executable.

After confirming whole normal success, no existing dispatch and remote `main`
still at the approved SHA, Codex dispatches **once** at 22:12:36 UTC:
[audit run 38476197008](https://github.com/smoe/gentle_rs/actions/runs/38476197008).
Its observed `head_sha` matches the full source above. Requested inputs are
`platform=linux`, `audit_followups=true`, `audit_public_vkorc1=false`.
The established Windows/headless/safety jobs also run. No package workflow,
public VKORC1 preparation, model call, push, tag or publication is included.

| Item | Current Evidence / Remaining Boundary |
| --- | --- |
| 1. Audit preservation / Linux validation | Old bytes and new original logs verified separately. Normal headless succeeds; full Linux job `115567420297` is running in the audit run, not yet accepted. |
| 2. 08.04 GUI/agent contracts | Synthetic native replay job `115567372647` is running; inspect receipt, exact operation semantics, T/refusal, source preservation and output parity before acceptance. |
| 3. Genuine live agent | macOS CLI-build job `115567372408` is running; building it is not authenticated/model acceptance. Model calls remain outside this authorization. |
| 4. Honest visual evidence | Historical human/version records remain intact. New synthetic replay captures are pending; no normal native lesson readability or scientific claim follows from their presence alone. |
| 5. Scroll benchmarks | CPU audit job `115567372185` is running; require all 24 unique wheel/next-frame observations and hashes. Native compositor/input latency remains separate. |
| 6. .12 package acceptance | Released `5e6201a6` opt1 evidence stays source-bound, with independent native acceptance open. This run builds no package; later local `.13`/`package-opt2` is not accepted by these older-source passes. |

These are live handles, not passing verdicts. Follow this run rather than
restarting it or creating duplicate dispatches. Local `edc3ed5c` and the `.13`
rollover remain unpushed and separate from the authorized `024ed6c0` source.
No local Rust or overall six-gate completion is claimed.

### Completed Container Validation At `5e6201a6`

After the main integration, Codex follows the already-live
[container run 38322028567](https://github.com/smoe/gentle_rs/actions/runs/38322028567)
without dispatching, restarting or publishing. On 2026-10-13 00:43 UTC, its
first attempt is `completed/success` (last updated 00:27:55 UTC). Source and workflow
revision are exactly `5e6201a69a599cba87491c86cfd7311cb6fee173`, tree
`2379f8bdfd115e6c93c0614a2bfe6188bb93c020`; this predates merged candidate
`2fffd4c00727cd6c829da96a620fe3679214c570`. The tag-push event does not grant
publication: the downloaded candidate record says `publish=false` and
`mode=validate_only`, the build receipt says `pushed=false`, and the GHCR
publication job is explicitly skipped. No owner-confirmed release follows.

**Dated correction, 2026-10-13:** that observation preceded the owner's
01:02:04 UTC publication of the `.12` prerelease at the same `5e6201a6` tag.
The separate release-event [container run 38330177090](https://github.com/smoe/gentle_rs/actions/runs/38330177090)
published GHCR successfully; [installer run 38330176948](https://github.com/smoe/gentle_rs/actions/runs/38330176948)
attempt 2 completed package validation/publication by 10:05 UTC. See the
[published ledger and seven-asset inventory](release_notes/release_notes_v0.1.0-internal.12.md#publication-evidence).
The earlier validate-only record below remains unchanged. CI extracted-package
smokes passed, but this does not close Glen's independent native/package,
live-agent, scientific or performance gates, nor certify newer source.

[Build job 115052172732](https://github.com/smoe/gentle_rs/actions/runs/38322028567/job/115052172732)
runs on Ubuntu 24.04.5 (`ubuntu-24.04` image `20261004.327.1`) with Docker
Buildx 0.37.2, targeting `linux/amd64` / Debian forky. Its log installs Debian
`rustc` package `1.96.1+dfsg1-2`; no separate `rustc -vV` capture is supplied.
The actual locked three-binary, no-default-features build uses `package-opt1`
and `-j1`, finishing in 56m 15s. The receipt matches the explicit Cargo recipe:
opt-level 1, LTO off, 256 codegen units, incremental/debug off, assertions and
overflow checks on, unwind and no stripping. Lockfile SHA-256 remains
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.

The completed smoke step runs default/CLI/MCP/docs entrypoints as the
unprivileged runtime user with `--network none`, checks dependencies and absence
of GUI/JS/Lua binaries, imports the Python wrapper, and renders nonempty RNAPKIN
hairpin SVG/PNG. RNAPKIN 0.3.9
is built separately under unoptimized `dev`; its corrected helper lock digest
is `f8c593a086f1bf6277a1886aa879b90b39f57ce163964de7cbf851b927702dfb`.
This is packaging/rendering smoke, not RNA-folding or scientific validation.

| CI-Executed Binary | CI-Reported SHA-256 |
| --- | --- |
| `gentle_cli` | `c4dd32476b2c3b99097bfd8839cf83b66843e640a109961098fc835da632c5a1` |
| `gentle_mcp` | `1a28058d88a09e353d908c90ab8e0c2e6a5ef3905ac39d5524bcdb09e2ae1b9b` |
| `gentle_examples_docs` | `a2b7657c641db15f24a5141966542f981f82ebb900461a7a180cc4d9d669d2de` |

All three self-identities name `.12`, build `1791848560452` and exact `5e6201a6`.
The receipt binds image ID
`sha256:57351ae9cd1c53c3a4308680057cc4c807853b52a8990c1bc77144ad31cd58d1`
and image digest
`sha256:1328e7335a58e1f7c2ef92d4b06c95ed0371c4594a71b3b3f7db9ee413d2d3bd`.
No image or binaries are downloaded; these are CI-measured binary hashes, not
an independent local rehash of executables.

Retain the following small original ZIPs and lossless decoded log under
`/private/tmp/gentle-container-38322028567/`. Independently check ZIP size/digest,
single-member inventory and candidate/build/recipe/helper-lock/self-identity
bindings before recording the pass:

| Retained Evidence | Size / SHA-256 |
| --- | --- |
| `release-candidate.zip`, artifact `11778863187` | 401 bytes; `81421874f3ad7101c4adf2157be90d0a44ebaf6d456f898a3cfd1b78380c8f65`; raw JSON member `cf58c3650d9d5332fb03cd2267c4c8c983b2465d1f0484428794422b6a2bb087`. |
| `container-evidence.zip`, artifact `11781064566` | 1,154 bytes; `0128405dd13a7e5fa85a0c8262f0ac7fd22364e466dcce9e5a5ede22c394be9d`; raw `container-build.json` member `682d7d10df1b907dd80121c7055c62b3335952147fc7a3bf41988e7f13ffc42d` (2,425 bytes). |
| `job-115052172732.decoded-log.json` | Saved JSON `b30404cbbf4f735f83657bc5509dc4cd2b6a66097b9082e6972773b8bbbd554d`; decoded-log UTF-8 `9b2f13565a42e42957f7afd6815fc74fbbe6a7e4f650a05fd9d446f1cb941653`. These are representation hashes, not raw HTTP/archive digests. |

The extra Docker build-record artifact is not downloaded. This accepts only the
observed headless container build/smoke at `5e6201a6`. Three native packages,
independent extracted-package smokes, full Linux and fresh final-source
GUI/live-agent/scroll gates remain open. Remote main still names `5e6201a6`
and the full merged `2fffd4c0` has zero Actions runs at observation time; retain
its 69 offline Python checks separately. The new ledger-only commit is not a
runtime change or an accepted build source. Use the full final pushed SHA in
the [consolidated Glen request](glen_final_candidate_execution_handoff.md),
without inheriting this baseline result. No local Rust or publication occurs.

### Source-Bound CI Update: 2026-10-13

Codex independently reads the public Actions run, latest-attempt jobs, artifact
inventory and decoded job logs at 2026-10-12 23:52 UTC (2026-10-13 locally).
[Run 38295036951](https://github.com/smoe/gentle_rs/actions/runs/38295036951)
is `completed/success`, attempt 1, exact source
`5e6201a69a599cba87491c86cfd7311cb6fee173`, tree
`2379f8bdfd115e6c93c0614a2bfe6188bb93c020`. Its checkout and executable
self-identities in the native job logs name that same SHA. The committed
lockfile SHA-256 is
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
These are ordinary default dev/test checks, not `package-opt1` packages.

| Executed Job | Directly Observed Log Evidence |
| --- | --- |
| [macOS 114964003981](https://github.com/smoe/gentle_rs/actions/runs/38295036951/job/114964003981) | macOS 26.6.2/arm64, Rust 1.99.0 (`b940084d7`). Focused primer multi-reference 9/9 and glossary 7/7; root lib under locked workspace suite 4,140 passed/14 ignored; workflow examples 95 passed. Default/JS/Lua/script/all-feature checks and schema/parity/process steps succeed. |
| [Windows 114963959579](https://github.com/smoe/gentle_rs/actions/runs/38295036951/job/114963959579) | Windows Server 2025/x86-64 MSVC, same Rust. Focused primer multi-reference 8/8 and glossary 7/7; root lib under locked workspace suite 4,079 passed/14 ignored; workflow examples 94 passed. The Unix-only external-import fixture is not compiled on Windows; this does not execute its real-reference or native presentation gate. Feature/schema/parity/process steps succeed. |
| [Linux headless 114963960058](https://github.com/smoe/gentle_rs/actions/runs/38295036951/job/114963960058) | Ubuntu 24.04/x86-64, same Rust; `cargo check --locked --no-default-features --bin gentle_cli --bin gentle_mcp` and dependency-boundary guard succeed. This is compile/dependency evidence, not Linux runtime/native acceptance. |
| [Safety 114963960095](https://github.com/smoe/gentle_rs/actions/runs/38295036951/job/114963960095) | Nine offline helper suites execute 223 tests plus the separate one-test tutorial-discovery guard; all succeed. Includes synthetic public-audit/package/GUI-contract helpers and LF/CRLF checkout checks, not actual GUI/model/package execution. |

The green `Render snapshot tests` step executes **zero tests** on both native
jobs; its exit status supplies no visual evidence. Linux full-suite job
`114964006033`, scroll `114963961155`, authenticated-audit preparation
`114963961298` and 08.04 GUI `114963961543` are all explicitly skipped.
The artifact inventory is empty. No binary hashes, native images, genuine model
responses, 24 new scroll observations or package/extraction receipts are supplied.
The explicitly ignored real-DuckDB companion is not accepted by a default suite.

Retain the connector's complete decoded logs losslessly as JSON strings in
`/private/tmp/gentle-ci-5e6201a6-20261013/job-JOB_ID.decoded-log.json`.
Original LF/CRLF characters remain encoded in the strings. These are hashes of
the saved JSON and its decoded-log UTF-8 re-encoding, **not** raw HTTP/archive
digests or binary hashes:

| Job ID | Saved JSON SHA-256 | Decoded-Log UTF-8 SHA-256 |
| --- | --- | --- |
| `114964003981` | `c870f941c6c45f6a6b1e3b335766cb71dcf1c08310a4c4c824c8b5e243dead28` | `bb7b939fc23bc15f71517b2a424075875f4a66dfc92af5d6be072576ac6a0687` |
| `114963959579` | `bd3ba12b8f910e7e823fb79c084780470a1b02a27ce6d657e8fd02370f0520a8` | `d35d46fa80d512f9723abebce2d5b119bf4a2d8703aae6cc9292b45b25a6bec3` |
| `114963960058` | `f5bdc9c04b848c1479912eccd9c2aeb4b3db007dd7a85314c2b35f453e10de0e` | `bdf8ae21c061ce89113ded4c92520985c09b831d31892b56383772c97743ea7d` |
| `114963960095` | `eae07fd1186736e6bf9fe3c4170d0964f346a841a41607891cdcbb978662a3aa` | `c18bc2a85c9dfe84c2e0232b2d93617078cd728ac3fe58464ab8ced426057a8f` |

At this observation, public remote main still resolves to `5e6201a6`. The
Actions query for each full SHA `e981469e`, `dfadf2a5` and local coordination
commit `9e5c8a05719ca6d81c68c88c47927b1a97113c30` returns zero runs.
No live verification process is observed, so this is not a verified wait.
The [consolidated final-candidate request](glen_final_candidate_execution_handoff.md)
is committed locally; its target still needs the owner's full final pushed SHA.
Older frozen requests and all remaining native/agent/scroll/package gates stay
separate and open. No push, dispatch, local Rust or scientific approval occurs.

### Current Main Snapshot: 2026-10-12

#### Completed Tutorial CI At `5383f7ec`

The owner confirms completion, and Codex independently reads
[CI run 38223352242](https://github.com/smoe/gentle_rs/actions/runs/38223352242)
and its latest-attempt jobs on 2026-10-12: exact head
`5383f7eceeb761d778b7d26efc7d0f5c4486d941`, attempt 1, `completed/success`.
Release-candidate safety, Linux headless CLI/MCP, Windows, macOS and the
selected-platform aggregate all succeed. The safety job's tutorial discovery
metadata and LF/CRLF checkout regressions succeed; both native platform jobs'
tutorial schema/drift, ClawBio parity, core CLI runtime smoke and full-suite
steps succeed. The review-manifest/LF-CRLF tutorial-CI repair is complete at
this SHA, alongside the integrated link repair. Remove those two tasks from
the active worklist, not their canonical repository records.

This completed result supersedes only the earlier in-progress observation
below. Linux full-suite, scroll, macOS authenticated-audit preparation and
08.04 synthetic GUI jobs remain skipped. It is a selected-platform CI pass,
not full-platform, native GUI, live-agent or `package-opt1` acceptance, and it
does not execute the separate 08.16/08.18 DuckDB or FLNA acceptance gates.

The two active workstreams below remain: current tutorial/native/FLNA follow-ups
and the [frozen execution handoff](glen_e981469e_execution_handoff.md) for
`e981469e` plus separate `dfadf2a5` audits. No SHA-exact execution receipts for
those frozen targets are supplied by this confirmation. The six-item checklist
does not advance; historical evidence, attachments and prompt copies remain
intact. No local Rust, push, dispatch, tag or publication occurs.

#### Glen Interim Feedback At `5383f7ec`

The owner relays Glen's interim review of
`5383f7eceeb761d778b7d26efc7d0f5c4486d941`. Glen reports that his fetch
updated remote/local main without checking it out, changing his dirty Clariom
worktree or pushing; a final re-fetch found no further upstream drift.

| Handoff | Reported Disposition | Remaining Acceptance |
| --- | --- | --- |
| `release_handoff_link_repair_prompt.md` | Five repairs from `1ee02cba` integrated; Glen reports `release-gate` 4 OK / 0 errors. | Link repair only, not package/release acceptance. |
| [Tutorial input review](glen_tutorial_inputs_review_20261012.md) | Report complete; Glen reports the attachment byte-identical to the repository copy. | Current 08.16/08.18 execution, C/D native checks and FLNA provenance/exon/native/scientific gates. |
| [Tutorial gap preparation](glen_tutorial_gap_preparation_20261012.md) | Partially complete; source/fixtures and three explicit unreviewed entries integrated. Glen reports 66 Python tests and LF/CRLF checkout checks passing. | GENtle Rust/DuckDB execution, native GUI and FLNA gates. |
| [Frozen execution handoff](glen_e981469e_execution_handoff.md) | Active; Glen reports the attachment byte-identical and no Actions runs at `e981469e` or `dfadf2a5`. | Six-item audit and separately pinned GenBank/variant and sequence-design rechecks. |

These are attributed interim results, not a new Codex execution receipt. The
message supplies no complete command list, environment, binary hashes or raw
log manifest for the 66 tests; do not promote them to runtime acceptance.
The three new lessons remain unreviewed. The separately refreshed 08.16 Codex
source/readability review is not native or scientific approval.

Keep two active workstreams, using the existing handoffs rather than duplicating
their commands:

1. Current tutorials: freeze one final main SHA and use only matching binaries
   for 08.16/08.18 Rust/DuckDB and catalog/parity replay, C/D and enriched-TSS
   native displays, and FLNA release provenance, human teaching-exon selection,
   native normal-size/2400-px and independent scientific acceptance. Keep FLNA
   outside the catalog until its activation gates pass.
2. Frozen audits: retain the exact `e981469e` six-item target and separate
   `dfadf2a5` GenBank/variant and experimental sequence-design target. A newer
   main pass cannot accept either; sequence design remains outside `.12`.

On 2026-10-12 Codex independently reads
[CI run 38223352242](https://github.com/smoe/gentle_rs/actions/runs/38223352242)
and its latest-attempt jobs: exact head `5383f7ec`, overall `in_progress`,
release-candidate safety and Linux headless jobs successful, macOS and Windows
in progress. Linux full-suite, scroll, macOS authenticated-audit preparation and
08.04 synthetic GUI jobs are skipped. This is a read-time snapshot, not a final
or full-platform verdict; no dispatch or rerun is requested or performed.

Retain canonical repository documents and historical evidence. Glen's suggested
inbox cleanup is not an owner deletion request: Codex removes no attachment,
untracked prompt copy or audit artifact. The six-item checklist below is
unchanged; none of its open gates closes on this interim feedback.

The subsequent owner-requested integration of `gentle_rs_2_main` commit
`ae491a0c` adds imported-TSS/repeat tutorial fixtures and Rust regression source.
It merges without conflicts and passes 93 scoped offline tests, including real
DuckDB I/O; details and raw log hash are in the changelog. The `f56fdf84`
checkpoint below remains source-bound historical evidence, not acceptance of
this new merge. The frozen execution-handoff source is now
`e981469e4beb637d2def12617ff974c76661ace5`, tree
`05246928299ebe32a6008cf38e391c1817b4ee2a`, with unchanged `Cargo.lock`
SHA-256 `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
No matching new binaries, full-platform or package verdict is inferred.

All 263 scoped offline tests pass at this exact SHA on macOS/arm64, Python
3.14.7 (39.148 s), including two real DuckDB 1.5.4 fixture I/O checks:

```bash
env -u GENTLE_TUTORIAL_BIN_DIR python3 -m unittest \
  scripts.test_release_candidate scripts.test_package_desktop \
  scripts.test_release_build_monitor scripts.test_gui_startup_summary \
  scripts.test_container scripts.test_tutorial_checkouts \
  scripts.test_primer_specificity_multi_tutorial scripts.test_ci_vkorc1_gui_audit \
  scripts.test_tutorial_gui_acceptance scripts.test_prepare_public_vkorc1_audit \
  scripts.test_tutorial_acceptance \
  scripts.test_tutorial_walkthroughs.WalkthroughSourceTests \
  scripts.test_tutorial_gap_lessons scripts.test_regulatory_motif_fixture -q
```

Retained raw log: `/private/tmp/gentle-glen-e981469e-offline.log`, SHA-256
`d1357b08c27dba92e748b5450aa6bb19270155abb5f4069317eac3d1ef8b4738`.
These helper/source/fixture checks execute no Rust build, native GENtle or live
model, and supply no package or scientific acceptance.

The no-push source bundle `/private/tmp/gentle-glen-e981469e.bundle` is 46,038
bytes, SHA-256
`8bd274c6ee4d811401100dcd584ec5bdafd5572e4e5e0e5875bbcea3771402bd`.
It contains the committed delta from prerequisite
`f56fdf849957d03d6f232c041cc74347d5ab471b` to the frozen merge, not untracked
worktree data or an audit archive. `git bundle verify` passes. An independent
temporary repository fetched only that prerequisite and bundle, then checked
out the exact new HEAD/tree above with empty `git status --porcelain` and the
same raw lockfile hash. The bundle is prepared locally, not uploaded or pushed;
the owner relays it to Glen with the handoff below. Later status-only commits
do not change the selected execution source or inherit acceptance.

### Preceding Main Snapshot: f56fdf84

The preceding execution-handoff source is
`f56fdf849957d03d6f232c041cc74347d5ab471b`, the merge of `gentle_rs_2_main`
after the `gentle_rs_3` integration. This exact source is already available on
GitHub; Codex neither pushed it nor dispatched the observed run. Its unchanged
`Cargo.lock` SHA-256 is
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
It includes the primer binding/manifest repairs, monitor readiness regression,
public-audit identity guards and merged readiness/tutorial changes.

All 258 scoped offline tests pass at this exact SHA on macOS/arm64, Python
3.14.7 (36.605 s). The command deliberately excludes installed GENtle binaries:

```bash
env -u GENTLE_TUTORIAL_BIN_DIR python3 -m unittest \
  scripts.test_release_candidate scripts.test_package_desktop \
  scripts.test_release_build_monitor scripts.test_gui_startup_summary \
  scripts.test_container scripts.test_tutorial_checkouts \
  scripts.test_primer_specificity_multi_tutorial scripts.test_ci_vkorc1_gui_audit \
  scripts.test_tutorial_gui_acceptance scripts.test_prepare_public_vkorc1_audit \
  scripts.test_tutorial_acceptance \
  scripts.test_tutorial_walkthroughs.WalkthroughSourceTests \
  scripts.test_tutorial_gap_lessons -q
```

Retained raw log: `/private/tmp/gentle-glen-f56fdf84-offline.log`, SHA-256
`ec429ce020d2bf27c4d8c7cf2371643416563f459e6abfd6484d9ef9baf651ba`.
These are synthetic process/receipt, source, mocked-fetch and disposable
LF/CRLF checkout checks, not real GENtle, native GUI or live-model acceptance.

The already-running push-triggered
[run 38181932226](https://github.com/smoe/gentle_rs/actions/runs/38181932226)
checks out this same full SHA. Independently inspected completed job logs prove:

- [Release-safety job 114596149687](https://github.com/smoe/gentle_rs/actions/runs/38181932226/job/114596149687)
  passes nine offline suites (221 tests) and the separate one-test discovery
  guard on Ubuntu 24.04/x86-64, Python 3.12.15. Raw log
  `/private/tmp/gentle-glen-f56fdf84-policy-job.log`, SHA-256
  `426cfc3d5cdd208e7c11d1afbd97ef6a113d27969dc4b252d2420d65fbf99ff8`.
- [Linux headless job 114596149703](https://github.com/smoe/gentle_rs/actions/runs/38181932226/job/114596149703)
  passes `cargo check --locked --no-default-features --bin gentle_cli --bin
  gentle_mcp` and the normal-dependency desktop exclusion guard. Ubuntu
  24.04/x86-64, Rust 1.99.0 (`b940084d7`), dev/check only; no runtime executable
  or package acceptance. Raw log
  `/private/tmp/gentle-glen-f56fdf84-headless-job.log`, SHA-256
  `a062f59914023041f13078406b7b6271abe269f3c1d11a5b6df335c0778318c2`.

At this inspection, Windows `114596149585` and macOS `114596168821` have passed
their initial Cargo checks and are running agent tutorial discovery regressions;
neither full-platform verdict is accepted yet. The full Linux job and opt-in
native GUI, authenticated-audit CLI build and scroll jobs are skipped, not
passed. No fresh native/public/scroll or package artifacts are accepted.
Later local evidence-documentation commits do not change this frozen handoff
source or inherit its CI results. Historical evidence and unrelated worktree
files remain untouched; no local Rust execution or new push/dispatch occurs.

The existing local dev `gentle_cli`, `gentle_mcp` and `gentle_examples_docs`
each exit zero for `--version` but self-report source
`150a339074e7a2d282567d791a8f50834b70ff9c`, not the frozen handoff SHA.
They are not used for fresh runtime or live-agent acceptance. No local
`target/package-opt1` native payload is present. Checking existing binary
identity executes no Cargo build, check or test.

### Preceding Offline Snapshot: 1737b01e

The exact tested script/contract revision is
`1737b01e7b3f7083ce8f225bea5eb3b54b4d41e4` on local `main`. It includes
`dd9ba573`'s corrected synthetic primer bindings, the test-only monitor readiness
repair `4b6a2569`, and the clean-source public-audit repair `1737b01e`.
Its unchanged lockfile SHA-256 is
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
All 241 scoped offline tests pass on macOS/arm64, Python 3.14.7 (34.178 s).
The tested scripts, documentation inputs and checkout policies match that
commit; unrelated `paper`, `output/` and the untracked release-handoff draft
are not used, staged or altered. Retained raw log:
`/private/tmp/gentle-glen-followups-1737b01e-offline.log`, SHA-256
`b952cb7913458851f2401031d33027fc1c231f13c6e82ed31b366feb1b4a5002`.

```bash
python3 -m unittest scripts.test_release_candidate scripts.test_package_desktop \
  scripts.test_release_build_monitor scripts.test_gui_startup_summary \
  scripts.test_container scripts.test_tutorial_checkouts \
  scripts.test_primer_specificity_multi_tutorial scripts.test_ci_vkorc1_gui_audit \
  scripts.test_tutorial_gui_acceptance scripts.test_prepare_public_vkorc1_audit \
  scripts.test_tutorial_acceptance -q
```

These tests use synthetic executable/receipt data, disposable Git repositories
and mocked public fetches. They are neither current Rust/native acceptance nor
an execution of public biological data or a live model. The monitor's
fixed-startup-wait failure is reproduced independently with a 300-ms Python
startup delay; production monitoring is unchanged. The public audit now checks
the shared clean-checkout and exact CLI-source identity before evidence creation
or fetch, then rechecks inputs before a pass. Both newly byte-bound public
source/workflow inputs have scoped LF policies and LF/CRLF negative controls.
No historical receipt, image, project, Criterion tree or verdict is rewritten.
Syntax, help and whitespace checks pass; session-close reports four OK, two
warnings (unrelated worktree dirt and manual scope confirmation), zero failures.
No push, workflow dispatch, Rust build/check/test, tag or publication occurs.
Later documentation commits do not inherit runtime/package acceptance.

The preceding 2026-10-11 checkpoint and CI observations retain their original
scope below; none supplies a native/package pass for the new script revision.

The preceding locally tested runtime/workflow gate source was
`6b39e5d80e880f738e4bc0aace1a7c19d56aee55`, including the earlier fixture
repair `b2307cf4` after the owner's tutorial-visual-contract merge `7c952e98`.
Its lockfile SHA-256 is
`ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
This differs from the historical audit lockfile; older package receipts cannot
accept this source/lock pair. No repaired-source binaries or packages exist
as acceptance evidence yet. This snapshot and its commits remain local.

The preceding CI source is `f39bed3ce79211cdc5cbde710137e581e2adf129`, with
the same current lockfile. In
[macOS/Windows run 38159995322](https://github.com/smoe/gentle_rs/actions/runs/38159995322),
[Windows job 114530710419](https://github.com/smoe/gentle_rs/actions/runs/38159995322/job/114530710419)
completes successfully, including default/JS/Lua/all-feature checks, binary and
process smokes, full tests and catalog/parity checks. This is ordinary CI,
not `package-opt1` acceptance. macOS job `114530736672` and
[Linux job 114530773148](https://github.com/smoe/gentle_rs/actions/runs/38160010120/job/114530773148)
both stop at the Unix-only multi-reference fixture's manifest read:
`primer_specificity_multi.rs:1566:63`, `NotFound`, eight passes/one failure.
Their later full suites do not run. Both runs' Linux headless and policy jobs
pass; opt-in native/public/scroll audit jobs are skipped, not newly accepted.
The second run's Windows job `114530754013` is still running at this inspection;
the already completed Windows job above supplies the bounded Windows evidence.

`b2307cf4` repairs only the fixture's cache-root/install-directory confusion:
use the fresh cDNA preparation report's sequence parent for the manifest,
retaining actual-kind mismatch, atomicity, import/replay and tamper assertions.
On the macOS working tree committed as this SHA, scoped edition-2024 rustfmt,
whitespace and `python3 -m unittest scripts.test_primer_specificity_multi_tutorial
scripts.test_release_candidate scripts.test_tutorial_checkouts -q` pass
(93 tests, 22.873 s). No local Rust build/check/test runs. Corrected-source
Rust verification requires a newly authorized push before CI can see this
commit; it is not waived or substituted with these source-only tests.

`6b39e5d8` closes an authored package-contract gap, not package acceptance:
the old native smoke accepted 75 synthetic stale/unversioned/ambiguous cases
across all three layouts. Require all five binaries' exact package/source
observations and hashes, including the publication helper's new shared version
probe, in extracted-package smoke, build receipts and aggregate collection.
Retain receipts outside the payload without overwriting evidence. The new
publication process regression is authored but unrun under the Rust restriction.
Existing historical build receipts lack these observations and are not accepted
by the new collector; their bytes and original verdicts remain untouched.

On exact `6b39e5d8`, macOS offline package/policy/container/checkout/native-audit/
public-audit checks pass all 162 tests (26.485 s), using synthetic stand-ins and
no local Rust or real native GENtle execution. Retained log:
`/private/tmp/gentle-native-identities-6b39e5d8-offline.log`, SHA-256
`8041658649b0b2a6af05b6cebae6e803ed9e56c97eb999deee7e70da355b8508`.
Scoped Rust formatting, release YAML/embedded Python parsing and whitespace pass.
Session-close reports four OK, two warnings (intentional `paper`/`output/` dirt
and manual scope confirmation), zero failures. Fresh final-SHA CI/native audits
and package-opt1 builds/independent package checks remain pending; no new push
or workflow dispatch occurs. Later documentation commits do not inherit a
native/package verdict from these offline tests.

| Item | Current Status | Required Evidence |
| --- | --- | --- |
| 1. Audit preservation | Complete; original bytes unchanged | [Original archive receipt](audits/glen_followups_20261008/original_archive_30f23aa4/README.md), pinned report/archive/manifest and all 721 safe file hashes at `30f23aa4`. |
| 2. Tutorial 08.04 GUI/agent contracts | Historical bounded pass at `3ac4a02d`; current source pending | [Retained native receipt](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#synthetic-native-contract): all 115 hashes, 19 focused Rust tests and eight Linux/X11 dev input steps. Explicit T, completion transition and three full oracle `Seq` records verified. Synthetic 20-bp guard, not public-locus GUI or scientific approval; retain fresh final-source evidence. |
| 3. Live-agent acceptance | Historical genuine public-starter pass at `3ac4a02d`; current source pending | [Executed receipt](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#genuine-live-agent), macOS/arm64 dev CLI and authorized login. Clarification then two unexecuted `ask` suggestions; reviewed manual execution preserves both sources and all four full `Seq` records match fresh CI. No credential/model/scientific claim; require final-source acceptance. |
| 4. Honest visual evidence | Historical bounded pass at `3ac4a02d`; clean-source guard passes offline at `e981469e` and in earlier Linux CI at `f56fdf84`, fresh execution pending | [Retained evidence](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#public-reference-and-base-projection): 115 native/44 public hashes, four unchanged synthetic captures and exactly reproduced/inspected public SVG. Sole 588 C/T in 1,089-bp exports, display-only uppercase and explicit non-claims. New guards refuse dirty inputs and stale/changed binaries; they do not execute the public workflow or supply a public-locus GUI verdict. |
| 5. Scroll benchmarks | Historical deterministic CPU pass at `3ac4a02d`; current source pending | [Retained scroll receipt](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#scroll-cpu-cases), Ubuntu 24.04/x86-64 dev. All 25 hashes and 24 bidirectional fixture/size cases, bounded viewport and unchanged content. No timed regression/native-latency claim; rerun the final-source cases without replacing historical evidence. |
| 6. .12 package acceptance | Open; selected suites and independently inspected build-only `package-opt1` container pass at `5e6201a6`, not final-source acceptance | [Container receipt/logs](#completed-container-validation-at-5e6201a6) bind all three CI binary identities, recipe and network-disabled smoke; publication is skipped. Main CI's macOS/Windows primer suites pass 9/8 (Unix-only fixture absent); full Linux/native-agent-scroll are skipped and that CI run has no artifacts. Require fresh full-platform, three native packages, independent extracted-package smokes and container evidence at the owner's final pushed SHA. Older frozen results remain historical; no inherited verdict. |

### Glen Execution Handoff: 2026-10-12

The owner now requests a forwardable message and committed supporting reading,
which the owner will push. Use the focused
[e981469e execution handoff](glen_e981469e_execution_handoff.md) for the bundle,
commands and returned-evidence requirements. The subsequent `gentle_rs_3`
merge `5b599808` includes a production GenBank/SPDI fix absent from `e981469e`;
the requested older target cannot accept this newer main or a release.
The focused handoff also consolidates the owner's supplemental audits at exact
`dfadf2a554dd8a23afae923592baa65d37162bc9`: GenBank/variant safety and a separate
experimental sequence-design recheck. Return their evidence independently;
neither replaces these six gates, and sequence design is outside `.12` acceptance.

The owner assigns Glen assistance with all remaining execution gates. This
handoff supersedes the earlier `f56fdf84` target: use exact merged source
`e981469e4beb637d2def12617ff974c76661ace5`, not moving main or an older binary.
Relay the small checksum-bound bundle above; Codex performs no GitHub push or
credential transfer. Make a new isolated checkout on Glen's approved
external runner, never reset his primary tree. With `BUNDLE` pointing to the
received file and `CHECKOUT` naming a new directory, reconstruct it as follows:

```bash
git init "$CHECKOUT"
git -C "$CHECKOUT" fetch --depth=1 https://github.com/smoe/gentle_rs.git \
  f56fdf849957d03d6f232c041cc74347d5ab471b
git -C "$CHECKOUT" bundle verify "$BUNDLE"
git -C "$CHECKOUT" fetch "$BUNDLE" HEAD:refs/heads/glen-audit-source
git -C "$CHECKOUT" checkout --detach e981469e4beb637d2def12617ff974c76661ace5
git -C "$CHECKOUT" rev-parse HEAD HEAD^{tree}
git -C "$CHECKOUT" status --porcelain
```

Before fetching from the bundle, require its exact byte count and SHA-256 above.
Before building, require the exact HEAD/tree, empty status and raw lockfile hash
above. Use an absolute bundle path so Git resolves it outside the checkout.
Keep caches, credentials and new evidence outside the source checkout. If any
required platform/tool/authentication is unavailable, report that gate as unrun
rather than substituting older binaries or a nearby SHA. This is a prepared
handoff, not a claim that Glen has received or executed it. The owner selects
relay through this chat; leave the complete handoff here for the owner rather
than inventing a direct transport or asking again for a delivery channel.

Glen's requested work, all without publication:

1. Confirm the exact checkout, lockfile, platform/toolchain and binary hashes.
   Rerun the 263 offline gates above and the targeted multi-reference import/glossary
   regressions before full Rust suites. Include default, headless CLI/MCP,
   JS/Lua/all-features and fresh Linux/macOS/Windows results; report primary
   failures separately from cascades and native-Windows omissions.
2. Execute the authored 08.04 semantic GUI contract in an approved isolated
   Linux/X11 runner, using the starter/oracle and shared audit helpers in
   `.github/workflows/ci.yml` and `scripts/ci_vkorc1_gui_audit.py`. Retain ordinary
   input, explicit T, completion transition and full-record oracle comparison.
   Synthetic 20-bp evidence never establishes the public-locus GUI or an assay.
3. Prepare fresh public inputs with `scripts/prepare_public_vkorc1_audit.py`
   from a clean checkout and matching CLI, then conduct the genuine live-agent
   acceptance with existing authorized login. Retain clarification, unexecuted
   review-first suggestions, ambiguous-alternate refusal and full-record parity.
   Never copy credentials, bypass recursion guards or silently execute commands.
4. Retain raw native captures and source-bound public base projections, labelled
   distinctly; verify the sole C/T difference from persisted sequence data.
   Do not recolour screenshots, reinterpret identical whole-map images, or
   publish captures. Keep scientific approval explicitly unclaimed.
5. Run the maintained scroll-event/next-frame smoke with fresh PATZ1 input:
   three fixtures, four sizes, both directions, all 24 observations, unchanged
   sequence/annotations and stable next-frame viewport. Preserve dev/profile
   metadata. A deterministic CPU pass is not native latency or a timed verdict.
6. Build native and headless-container `package-opt1` artifacts at the same
   frozen source/lock/recipe with publishing disabled. Verify actual extracted
   archive membership/checksums, all five native and three container binary
   identities/entrypoints, and independent package receipts per `docs/release.md`.
   Do not tag, publish, push our local branch or accept a mixed-source package.

Return commands and exit codes, exact SHAs/hashes, environment/profile, failures,
and a retained manifest/receipt for each gate. Only matching executed evidence
can advance this checklist; pending platform, live-agent or package checks stay
pending. No human scientific or release approval is requested by this handoff.

Run the following only on Glen's approved external build host, with
`EVIDENCE_DIR` set to a new directory outside the clean checkout. The existing
workflow defines dependency installation and isolation; the local no-Rust
restriction remains in force for Codex. Start with the failed boundaries:

```bash
cargo test -q --locked --lib -j1 engine::primer_specificity_multi::tests:: -- --test-threads=1
cargo test -q --locked --lib -j1 engine_shell::tests::glossary_ -- --test-threads=1
cargo test -q --locked --lib -j1 execute_introspect -- --test-threads=1
```

For the native contract, reuse `.github/workflows/ci.yml`'s six focused filters
and the `desktop-gui,gui-test-support` build, then execute:

```bash
python3 scripts/ci_vkorc1_gui_audit.py --repo-root . --binary-dir target/debug \
  --evidence-dir "$EVIDENCE_DIR/native"
python3 scripts/prepare_public_vkorc1_audit.py --repo-root . \
  --gentle-cli target/debug/gentle_cli --evidence-dir "$EVIDENCE_DIR/public"
```

The public preparation command invokes no agent. Use the retained
[genuine-agent procedure](audits/glen_followups_20261008/candidate_3ac4a02d/README.md#genuine-live-agent)
on the newly prepared starter with a matching binary and authorized login;
repeat its two requests, raw-state checks, refusal and reviewed manual execution.
Retain new responses and four full ordered `Seq` comparisons; old responses are
not reusable acceptance. For scroll, reuse `gui_operations --test` and the
workflow's fresh PATZ1 setup, executable discovery and 24-observation receipt,
not a guessed executable path or a Criterion timing surrogate.

For packages, follow [build-only verification](release.md#build-only-candidate-verification)
with this exact candidate and `publish=false`. Do not dispatch release/container
workflows without explicit owner approval, or select moving main after a later
documentation commit. Require actual native archive and container validation
receipts; the currently green policy/headless jobs provide neither.

### Historical Candidate Gates

The following records retain their original source and inspection-time scope;
references to a then-current candidate or running job are historical, not the
current main snapshot above.

The owner-requested rebase is now on local main
`b1547544980cb423dbda9ee57dfce81f8b193f5b`. The isolated codex worktree retained
all 34 audit commits and the subsequent status commit. Changelog and roadmap
resolutions preserve both histories and pending gates. The overlapping pair
execution keeps main's existing guarded rollback helper unchanged, with both
main/audit test families; the mixed-SNV/indel guard is unchanged. Evidence at
`6d8b4db7` stays accepted only for that
historical source. The current candidate's checklist is:

The current frozen runtime/workflow candidate is
`e42daa0216a3bc78c2c84fb05d5b972781e9e459`, pushed only to the codex branch.
It adds the separately committed test-only macOS fixture repair and early
exact-test CI step; no production path validation or biological behavior changes.
Both fresh run records report this exact SHA:
[Linux/Windows with focused/native/public/scroll audits 37915942296](https://github.com/smoe/gentle_rs/actions/runs/37915942296)
and [macOS/Windows 37915951526](https://github.com/smoe/gentle_rs/actions/runs/37915951526).
They are running, not accepted. Keep the remote branch frozen through these
gates and later explicit `publish=false` package builds. Before the later
owner-requested main integration, no local Rust build, main update, tag or
publication occurred. Session-close after the repair commit
reports five OK, one manual plan-fidelity reminder and zero failures.
The lockfile remains byte-identical. Until fresh acceptance is retained, items
2-5 below refer only to the historical `3ac4a02d` source, not this new candidate.
At `e42daa02`, the remote policy job `113772019301` passes all 181 Python tests;
retained raw log `/private/tmp/gentle-e42daa02-policy.log` SHA-256
`232a34e12fbfd811377ee8b1f8bebba3e6122afc09b89e556eb3767ca8ca625d`.
Linux headless job `113772019137` passes
`cargo check --locked --no-default-features --bin gentle_cli --bin gentle_mcp`;
raw log `/private/tmp/gentle-e42daa02-headless-check.log` SHA-256
`03493f7debd2c9b1b8d9d342c3fd3ee21009038c7cb5e8809948db579081fad2`.
These are partial gates only: the repaired native test, full platform suites,
fresh audit artifacts and non-publishing package acceptance remain pending.

The previous frozen runtime/workflow candidate is
`3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da`, with the unchanged lockfile above.
It includes separate container-identity (`57cc62bc`) and eleven-line primer
initialization (`3ac4a02d`) repairs. The remote codex branch was deliberately
held at that SHA until both package builds finished; subsequent evidence/status
commits and the test-only macOS repair are not new runtime acceptance.
Fresh non-publishing runs are
[Linux/Windows plus focused audits 37848423833](https://github.com/smoe/gentle_rs/actions/runs/37848423833)
and [macOS/Windows 37848434104](https://github.com/smoe/gentle_rs/actions/runs/37848434104).
Both run metadata independently report the exact candidate SHA. After the
fresh focused Rust and native replay steps passed, dispatch
[native packages 37850876088](https://github.com/smoe/gentle_rs/actions/runs/37850876088)
and [build/load-only container 37850885456](https://github.com/smoe/gentle_rs/actions/runs/37850885456)
with the same explicit candidate and `publish=false`; both builds succeed,
with both publication jobs skipped. Package acceptance remains open because
macOS full CI failed and independent extracted-bundle verification is pending.
Do not reuse the successful historical `6d8b4db7`
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

On 2026-10-09, independent native/scroll checks pass without local Rust
execution. Four unchanged synthetic PNGs are visually inspected: T and pair
action are visible; layered windows and light chrome/dark canvases are retained.
The public projection and genuine agent requests now pass their separate gates.
All nine retained receipts have scoped LF attributes and individual unprotected-CRLF
negative controls. The initial three-receipt seven-suite local check passes
158 tests (24.539 s), log SHA-256
`8a64b289c669388a5d510118e5ebcb1039583993876e90e9a1d0addd8ca1a33d`;
the expanded five-receipt checkout suite passes 26/26 (21.640 s), log SHA-256
`fddda1416a55309164d5ad2e57370ff9ee798728a66273156e6edcbe79813b73`.
Session-close reports four OK, two warnings (this intentional pending commit
and manual plan fidelity), zero failures; roadmap remains 400 lines.
Local evidence commits stayed unpushed while those package builds were running.
The expanded nine-receipt seven-suite check passes 158/158 on macOS (41.865 s),
log SHA-256 `364e4b411bb3953c2354118babe531b9044357bf5af41c4da841ffbdaa660f7f`.

Full Linux/Windows run `37848423833` succeeds: Linux lib 4,122 passed/13 ignored,
Windows lib 4,062 passed/13 ignored. Retained complete run log SHA-256:
`5ebdceba689533466582d6a7880e1540e85c4110f05c4341d5c69624983ad46f`.
Run `37848434104` has successful Windows/headless/policy jobs but failed macOS:
4,121 passed, one failed, 13 ignored. The primary panic is
`specificity_multi_output_paths_refuse_traversal_and_existing_evidence` at
`primer_specificity_multi.rs:525`, rejecting a supposedly valid new temporary
bundle. Raw failing log SHA-256:
`36d585d071df1aa47e5504a15df852d3a46fecc66b404bdf919405f389f2d6d2`.
macOS `/var` is a symlink to `/private/var`; the production validator explicitly
rejects symlink ancestors. Repair the fixture root, not that safety boundary,
and retain raw traversal and deliberate alias rejection before rerunning.
The separate repair canonicalizes only the trusted temporary fixture root,
constructs raw traversal as `OsString` and asserts that the parent component
reaches the validator. Existing-evidence and deliberate Unix alias refusals
now assert `InvalidInput`; production validation is unchanged. This authored
test repair is not yet a Rust pass: scoped formatting/Python checks precede
fresh remote targeted/full-platform verification at the next frozen source.
All three native jobs now run this exact test before long suites; a focused
offline policy test guards the command, order and absence of skip/failure hiding.
The repair's seven-suite offline check passes 159/159 on macOS (39.621 s),
log `/private/tmp/gentle-primer-path-repair-python.log`, SHA-256
`281b05762db8e5703920bb922c73b407cb5e81e67339a828e1c736555b5f9ca6`.
Scoped `rustfmt --check`, Ruby CI-YAML parsing and `git diff --check` pass;
no local Rust build/check/test was run. Freeze the next exact source after
this separate repair commit and rerun all three native gates before accepting
new packages. The successful `3ac4a02d` runs retain their original limited scope.

The successful build-only container receipt is retained locally at
`/private/tmp/gentle-container-3ac4a02d/container-build.json`, SHA-256
`3a4efcc3f0e3f2e8e763ded55db38c287f5426af96cff693a75b9f42456e9315`;
all three observed versions embed `3ac4a02d` and each has a SHA-256. The native
aggregate attributes receipt is retained at
`/private/tmp/gentle-native-metadata-3ac4a02d/gentle-v0.1.0-internal.12-release-attributes.json`,
SHA-256 `d4588ed3331bf22c0849fe4b01582fe6aa196dbf18de2ace1690b2454e96bf7b`.
Neither a green build nor these metadata receipts closes the failed full-platform
gate or substitutes for independently checking the extracted native packages.

### Package Byte-Verifier Preparation

A task-local read-only verifier at `/private/tmp/verify_gentle_native_package.py`
(SHA-256 `cd1746d5a2e6983efe7619b7207243be918649613157a84e1226c656d2926966`)
is prepared for the new packages without changing their source or running Rust.
It checks explicit source/workflow/lock/profile/recipe and `validate_only`, the
outer archive digest, safe exact payload membership, every `SHA256SUMS` entry,
five nonempty native binary formats/hashes, packaged resources and revision.
Linux/Windows checks stream the extracted bytes; macOS uses a read-only mounted
app root. This byte check does not execute entrypoints or supply native/full-suite,
GUI, latency or scientific acceptance. Those require their separate exact-run
evidence and, for local execution, a matching native platform.

Seven synthetic test groups in `/private/tmp/test_verify_gentle_native_package.py`
(SHA-256 `8db314cf1d061092fe0c4f9106985fb1cadaa14a40d3d12d713cc2aef803c003`)
pass, including stale sources/recipes, membership/hash failures and traversal,
symlink and special-member refusal. Test log SHA-256:
`6027ca906dcab6eaf0e89e604173ee50f09abcb0325955d3c8b1ca7fa375e282`.
A read-only verifier-development trial on the already retained historical
`6d8b4db7` Linux archive verifies 1,099 payload files and five binary hashes;
receipt `/private/tmp/gentle-package-bytes-6d8b4db7-historical.json` SHA-256
`bed7e7d13966854d78d6373feb6bf991fdf327fd0459682d655ffd09c39f85ce`.
An explicit `e42daa02` invocation refuses that historical receipt before reading
the payload. No original/historical bytes were altered or new package downloaded;
this preparation closes none of item 6's pending candidate gates.

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
