# Glen: Next Verification And Return

## Native RNA Compile Repairs, 2026-10-25

The owner's build exposed two primary integration defects: missing explicit
paths for the native processing/evidence modules, and a numeric journal count
shadowing the optional command-admission baseline. The missing UI-method errors
were cascades of the module failure. Both are repaired without disabling UI
handlers or weakening project-commit/stale/cancel safeguards. Two new LF/CRLF
source regressions fail before repair and pass afterward; they are not Rust
verification. Include this repair at the next ordinary integrated source.

Before the native rehearsal, run the default `cargo check -q --locked --tests`
and headless check below, then the existing `command_execution`, `rna_screen_`
and `rna_native_` filters. Retain replaced-project refusal, retained-operation
counts, artifact-only no-commit behavior and the existing native handlers.
Report actual source/toolchain/platform, primary failures and cascades through
Steffen. Do not repeat completed campaigns or impose a new source freeze.

## Native RNA Processing And Tutorial 10.04, 2026-10-24

Integrate `4603d029` (execution/bindings), `a87fb9c6` (setup/monitor),
`6d7d6ef3` (inspector/export) and the subsequent rehearsal/parity commit at
ordinary current main; do not freeze another candidate. The earlier unrelated
goal remains paused. No local Codex Rust/native run, private sequencing access,
capture, provider call or push is claimed.

Run focused checks, then normal full CI on that exact integrated SHA:

```sh
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 rna_screen_ -- --test-threads=1
cargo test --locked --lib -j1 rna_native_ -- --test-threads=1
cargo test --locked --lib -j1 command_execution -- --test-threads=1
cargo test --locked --lib -j1 rna_artifact_writes -- --test-threads=1
cargo test --locked --test capability_registry_parity --test parity_matrix_freshness
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_examples_docs
cargo test --locked --lib --features gui-test-support -j1 rna_screen_native_semantic_targets -- --test-threads=1
cargo fmt --all -- --check
python3 -m unittest scripts.test_rna_panel_tutorial scripts.test_tutorial_checkouts scripts.test_capability_maturity
```

The `rna_screen_` filter includes draft-save agent confirmation. Preserve
existing RNA MCP, GUI intent/palette, approval and command-service contracts.
The authored app test holds the live project write lock across native screen,
monitor, bound reopening and summary export; require two synthetic samples,
four reads each, retained class/member sets and no live journal append. This
is not native GUI evidence. Keep the channel-controlled late-cancel/concurrent-
job tests below and ordinary mutator stale rejection.

With a verified source-matched docs binary, check catalog/manifest/parity
generation and tutorial replay using the repo's normal commands. The five
new native command cells and GUI totals were projected from the changed
glossary; exact Rust generator freshness still needs execution. Maturity
projections were refreshed by their Python generator. The synthetic demo now
contains stable toy `gene_id` qualifiers and a stored resolved gene set so
native preparation uses the real binding path. New demo/panel/report hashes
legitimately change; legacy omitted-field hashes must not.

Follow [the GUI-only rehearsal](tutorial/10-04_targeted_nanopore_fastq_panel.md#gui-only-rehearsal)
and its [typed checkpoints](tutorial/10-04_targeted_nanopore_fastq_panel.md#typed-checkpoints-and-capture-targets).
Retain permitted synthetic source/binary/snapshot-bound screenshots of setup,
progress/outcome, nomination prefix, summary ambiguity and table verification.
Test native and hosted narrow-window scrolling. Missing dialog/control automation
coverage is `harness_gap`, not guessed clicks or a biology verdict. This chapter
remains manual/hybrid, not a new generic acceptance-runner contract.

Use an observable source-bound synthetic job for scheduling checks: other windows
and agent conversation stay usable; close/reopen does not cancel; four occupied
execution slots do not block marker control/plain inspection; table verification
uses bounded execution capacity. Verify cancel before/during/after publication,
marker-write failure and post-start control_error retention. Reopen saved bindings
after restart without claiming worker recovery. On Windows test handles, quoted/
Unicode paths, raw traversal, fresh atomic draft save/no overwrite and LF/CRLF
request/cancel identity. Do not use private downloads to slow down the rehearsal.

Inner/outer-agent rehearsal covers `ui open rna-read-panel`, discovery/build,
metadata-only `panel-screen --request ... [--save ...]`, approved screening,
`introspect runtime`, request-bound cancel, bound `panel-show`, summary export
and table verification. Writes stay `ask`; draft saving via MCP uses
`agent_execute_plan`, scientific operations remain `op`/workflow/script accessible.
Compare canonical scientific payloads, not host-local job IDs or timestamps.
No model call was made here. Return SHA/toolchain/platform, exact commands,
ledger/hash locations and classified failures via the main-maintaining thread.
Queueing and short result-graph attachment stay separate.

Exercise existing size budgets explicitly: interactive result metadata is capped
at 16 MiB. An oversized artifact run retains its bound compact receipt/full
disk report; an oversized inspection can refuse rather than silently truncate.
Measure actual metadata size for the intended panel/file counts before claiming
large-library GUI coverage. Paging beyond that budget is not implemented here.

## File-Only RNA Execution And Result Bindings, 2026-10-24

At the ordinary integrated source, run the focused `rna_panel_` lib filter,
existing `command_execution::tests`, `agent_feedback`/approval/MCP parity tests,
glossary/capability/parity freshness and `cargo check -q --locked --tests`.
No new freeze, local Codex Rust build, private FASTQ run or native acceptance is
claimed. Regenerate/replay the synthetic tutorial with a source-matched binary;
new submitted-request fields legitimately change new report hashes, while
legacy omitted-field canonical bytes and cancel identities must remain stable.

New engine tests cover exact frozen legacy request bytes (LF/CRLF), reviewed
panel replacement before streaming, submitted/effective manifest normalization,
changed/deleted source files after publication, binding tampering, bound shared-
shell inspection and summary/legacy metadata. New service tests use channels to
hold the publication-to-retention window, request late cancellation and edit
live structure; require retained Completed evidence, zero committed operations
and no overwritten project edits. A separate held mutator must still commit
after artifact retention, proving there was no intervening live journal append.

Exercise callback cancellation and foreign-marker control errors after the
first verified file. Require Cancelled/Failed lifecycle states respectively,
retained bound partial reports and unchanged scientific outcomes. Fill all four
execution slots and hold the live project lock: marker control and metadata
inspection must still finish in their independent lanes; explicit expensive
verification must report execution-capacity refusal. Registered panel builds
remain normal project commits with correct mutation feedback. Oversized file-
only output must retain a small bound envelope and the full published report.

Preserve existing mutator stale/cancel rejection and bounded history tests;
do not infer file-only policy from `state_changed=false` or caller JSON. Review
`submit_parsed` freezing versus raw-text worker parsing, agent project-change
handling and compact-result feedback. Cross-process cancellation remains the
request-bound marker; these process-local job IDs are not portable handles.
Native setup/monitor/inspector/export and the rehearsal are implemented;
integrated runtime/native acceptance is described above and in the
[processing plan](rna_panel_gui_processing_plan.md).

## Gene-Set RNA Panel Preparation And Graph, 2026-10-24

At the ordinary integrated source, without another freeze, run:

```sh
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_mcp -- --test-threads=1
cargo test --locked -p gentle-shell ui_intent
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo check -q --locked --tests
```

Run the normal tutorial/catalog/maturity generators and freshness checks; the
source/metadata Python tests do not replace Rust execution. New authored tests
cover loaded-group discovery, immutable 10-to-20-locus extension, source-hash
drift, gene-set binding, registration/save/reload/undo/redo, GUI/SVG graph edges,
UI-intent lifecycle, stale workers and MCP read-only discovery. Require legacy
omitted-field digests and the synthetic screen/summary outputs to stay unchanged.

Follow tutorial 10.04 in a disposable native and hosted project: resolve an
explicit set in Gene Set Inspector; import its annotated loci; refresh, bind exact
stable IDs, review provenance/transcripts and build. Verify the set/source/panel
graph, reopen by right-click and double-click, save/restart, then extend with an
expanded set. Undo must remove only registration, not the external file; old
panels/screens must remain unchanged. A moved panel file must not prevent
snapshot inspection, but must prevent extension until explicitly resolved.
Missing/ambiguous IDs, unresolved members and stale source hashes are errors,
not zero transcripts. Other windows must remain usable during preparation.

Compare CLI `panel-build --list-sources [GENE_SET_ID]`, the copied build JSON,
MCP confirmed `op`, and inner-agent approved builds using the same scientific
payload. Discovery is read-only; file/project writes stay ask/confirmation-bound.
Live-model invocation and native semantics remain unverified locally. This is
panel preparation, not a new native FASTQ-screening form, isoform-discrimination
acceptance or a private sequencing run. Return exact source/platform/toolchain.

## PCR-Strip Test-Contract Repairs, 2026-10-23

Glen's two isolated failures from `fc75114a` are corrected in tests only:
16 successful exports append 16 external-operation receipts without changing
project state; invalid count types are refused during operation execution, not
generic shell parsing. Refusals must preserve journal length and existing model
bytes. Valid exports, quoted paths, legacy serialization and capability checks
remain covered. At the next ordinary integrated source, rerun only:

```sh
cargo test --locked --lib -j1 labware_pcr_strip_ -- --test-threads=1
```

Return actual source/toolchain/platform and results through Steffen. No repeated
mesh render, completed acceptance campaign or new source freeze is requested.
Codex ran no Rust verification; production journaling/parser behavior is unchanged.

## Mid-Run FASTQ Control-Error Retention, 2026-10-23

At the ordinary integrated SHA, run the focused `rna_panel_` and
`rna_panel_mcp` lib filters, agent-feedback/agent-execution and shared-shell
approval tests, existing glossary/parity/catalog freshness, and
`cargo check -q --locked --tests`. No new freeze or private chunk is requested.

New `rna_panel_mid_run_marker_errors_preserve_verified_units_with_explicit_outcomes`
injects faults through progress callbacks without sleeps. Require preserved
first-unit bytes/hashes and `status=control_error` at file boundary, during
intake and at finalization; malformed/nonregular markers follow the same policy.
Check zero verified units remain unassessed, matching cancellation remains
cancelled, and a later callback cancel cannot overwrite the control error.
The hard checksum/resource-error regression must still publish no result, even
after an earlier unit verified. Preflight mismatch remains a hard refusal.
The binding regression rejects missing diagnostics, mismatched paths and
relabelled success/cancel/legacy outcomes. MCP rehearsal now inspects the
retained control-error report and exports a partial summary with the same cause.

Replay tutorial 10.04 in a disposable project through CLI, MCP and the live
inner agent. Require a visible warning, failed screening feedback and retained
named missing inputs, while named inspection can complete. A summary keeps
`source_runs[].control_error` and remains partial even if finalization failed
after every file verified. The agent should explain fresh control/output paths
and explicit remaining-file retry, not silently restart or edit the marker.
Native Windows must verify closed handles before partial-unit cleanup and
publication; retain Unix nonregular/symlink coverage. Optional diagnostics are
omitted when absent so historical v1/v2 semantic digests remain unchanged.
No Rust/native checks were run locally; report exact source/toolchain/platform
and executed results, not a new biological or throughput verdict. The added
selection/cumulative-budget tutorial cautions do not implement the separate
index rewrite, resource-limited continuation or legacy-path follow-ups.

## Stored-Plan And RNA MCP Test Expectations, 2026-10-23

The owner's two latest library failures were obsolete approval expectations,
not evidence that either production gate should be relaxed. The shell test now
requires denial for a parameter write (including falsely read-only metadata),
unchanged denied state/history and one confirmed operation. A separate read-only
shell plan still executes without confirmation. The RNA MCP rehearsal first
denies raw `op` verification without `confirm=true`, checks untouched project/
report bytes, then explicitly approves it and retains scientific-payload parity
and unchanged semantic state. Shared-shell read-only exemptions are unchanged.

At the next ordinary integrated run, include these focused checks; do not
restart completed acceptance campaigns or impose an interim source freeze:

```sh
cargo test --locked --lib -j1 engine_shell::tests::execute_agents_execute_plan_ -- --test-threads=1
cargo test --locked --lib -j1 mcp_server::tests::rna_panel_mcp_rehearsal_matches_direct_scientific_payload -- --exact --test-threads=1
```

Use the normal run-owned target and return actual source/toolchain/platform,
failures/skips and retained logs through Steffen. Codex ran no Rust verification.

## RNA Table-Verifier Import Repair, 2026-10-23

The owner's compiler found E0432 in `rna_panel/verification.rs`: `sha2` was
not declared by the root crate. The verifier now uses the same already-declared
`ring::digest::Context` / SHA-256 as FASTQ intake, keeping bounded streaming and
all exact hash/size/path checks. No dependency or lockfile update is needed.
The previous repair's two unused `Engine` imports are removed. Include the
existing `rna_panel_summary_verification_is_explicit_bounded_and_hash_bound`
Rust regression in the ordinary `rna_panel_` filter below; source/format checks
are not compiler acceptance. No new campaign or repeated historical checks.

## RNA Confirmation, Starter And Tutorial Freshness Repairs, 2026-10-23

The reported five failures have three causes: the stored-plan executor trusted
shell candidate flags instead of shared command mutation policy; the demo saved
unreconciled lineage/containers; and 10.04 lacked a review entry while its catalog
projection had noncanonical indentation. The repair checks parsed commands for
every candidate kind, saves an engine-reconciled synthetic starter and registers
10.04 explicitly as unreviewed. It preserves V2 co-nomination, cancellation and
table-verification contracts from `874ff657`.

Do not relax confirmation, unchanged-state assertions or freshness checks.
New Rust tests exercise mislabeled shell/raw-op writes denied before artifacts,
successful explicit approval, unchanged blocked state/history and read-only
inspection without approval. The MCP rehearsal now also checks that loading
the saved starter is idempotent before its unchanged-state comparison.

At the next ordinary owner-pushed source, run the existing RNA/MCP checks below
plus `cargo test --locked --lib -j1 agent_execution::tests:: -- --test-threads=1`.
Re-run the three originally failing `workflow_examples::tests` cases:
`tutorial_source_units_generate_committed_catalog`,
`tutorial_catalog_check_passes_on_committed_tree` and
`tutorial_check_passes_on_committed_tree`. Keep canonical generator checks;
do not repeat completed gel or earlier glossary/protocol acceptance campaigns.
Report failures/skips and actual source/toolchain/platform through Steffen.

Codex used the existing `f37f7039` documentation binary for catalog regeneration
and freshness diagnostics; its generator source is unchanged at the integrated
revision. This is not new-source Rust compilation or acceptance of the revised
agent/demo code. Both existing-binary catalog and full tutorial freshness checks
pass after the repair; the two historical stale-review warnings remain. All 109
Python/source/LF-CRLF checkout checks pass without skips, alongside formatting
and documentation/maturity checks. No local Rust build/test harness, push,
dispatch or private-read run is performed. Existing scientific artifacts and
tutorial report stay intact.

## Narrowed FASTQ Repairs, 2026-10-23

The follow-up keeps the existing seed/chaining algorithm and fail-closed
budget policy. V2 records retain class plus members; the renamed counts are
partial co-nomination, never whole-read compatibility. The demo_1 regression
requires A's junction `[40]` and B's `[]` while both remain nominated. V1
screen aggregation uses `legacy_class_unresolved`; v1 summary JSON must keep
its original semantic digest. Named sample coverage distinguishes incomplete
and wholly unassessed libraries. Progress counts processed reads and compressed
read-ahead bytes, not verified read fractions. Matching pre-start cancellation
works; stale/mismatched markers fail, and polls are throttled. Explicit table
verification reads three allowlisted tables only, capped at 256 MB; ordinary
display and bounded nomination prefixes never silently rehash large sidecars.

At the integrated SHA, using the usual run-owned target, Glen/CI should run:

```sh
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_mcp -- --test-threads=1
cargo test --locked --test glossary_cli_usage
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo check -q --locked --tests
```

Replay all of 10.04 synthetically, through CLI and MCP/shared-agent commands,
including `panel-show summary/report.json --verify-artifacts`. Native Windows
must cover quoted paths, closed handles, TSV LF/CRLF tampering and publication;
Unix adds symlink refusal. Check named missing sample B after cancellation
following sample A, and wholly unassessed samples with blank counts. Keep the
original `artifact_content_verified=false` for bounded nomination prefixes.
The live inner agent must discover the new flag and explain co-nomination,
unknown coverage and pending integrity separately; no OS helper is required.
Retain exact SHA/toolchain/platform/test counts and distinguish failures from
cascades. No local Rust or native execution is claimed here; no private reads,
performance/biological verdict, new freeze, Glen contact, push or dispatch.

## FASTQ Supporting-Read Summary Follow-Up, 2026-10-23

The new shared `ExportRnaReadPanelSummary` / `rna-reads panel-export-summary`
consumes only stored screen reports and retained panels; no original FASTQ or
large nomination TSV is reopened. It exports sample/transcript support and
exact ambiguous sets, preserves cancelled/zero-verified-file uncertainty,
refuses duplicate input hashes and incompatible panel/policy/sample bindings,
and publishes a content-bound JSON marker after closed, hashed LF TSVs.
`panel-show` reopens summaries without journaling. Tutorial 10.04 now supplies
an engine-created summary request and teaches non-additive ambiguous support.
The demo's two library headers are deliberately distinct; identical byte
inputs would otherwise correctly be refused as duplicate evidence.

Run the existing `rna_panel_` filter, glossary, capability/registry parity,
canonical catalog/manifest freshness and `cargo check -q --locked --tests`
listed below at the ordinary integrated source. Seven new `rna_panel_summary_`
tests cover actual aggregation, source deletion/no reread, LF/CRLF report JSON,
stable ordering/table hashes, chunk/sample pooling, duplicates, cancellation,
tampering/count overflow, traversal/no-overwrite and shared-shell inspection.
The existing MCP/agent-approval rehearsal now covers export and reopening.
These Rust tests are authored, not run locally; neither source checks nor
rustfmt establish their compilation or execution.

Native Windows must exercise spaced/quoted paths and handle closure before
publication. In the live inner agent, discover the exact typed request,
require approval for export, reopen it via `panel-show` and explain that the
numbers are supporting reads, not molecule copies or confirmed isoforms.
Do not silently redistribute or add ambiguous per-member support. Keep prior
native/background/real-library acceptance open; do not process private chunks
or contact Glen from this worktree without explicit authorization.

## Provider Preview Fixture Follow-Up, 2026-10-22

After integrating the owner's Windows LF commit `1cfea9e3`, three failures in
the supplied macOS log are already addressed in source: the numeric glossary
fixture (`254091c7`) and optional-source discovery in two provider tests
(`bfe685ed`). Keep their focused regressions; no new parser broadening or
missing-explicit-source tolerance is intended.

The remaining preview-reopen test mixed its synthetic preview/apply Context
with host-wide discovery in the writer. It now injects the same fixture source
list into the shared writer implementation. The public entry point still
discovers and checks every configured source; raw traversal, symlink, source-
destination and no-overwrite refusals are unchanged, not bypassed for macOS.

On the next ordinary verification, use a clean checkout and run-owned target:

```sh
cargo test --locked --lib -j1 external_service_providers::updates::tests:: -- --test-threads=1
```

Include macOS for the original reopen failure and the synthetic symlink-source
ancestor guard; retain native Windows coverage for the existing drive/verbatim
tests. Report executed test counts, failures/skips, actual source/lockfile/
toolchain/platform and logs through Steffen. The completed glossary/all-features
return below is historical evidence, not a second full campaign. No interim
freeze, bundle, push, dispatch, native capture, provider call or package build
is requested here. Local source checks do not certify the Rust paths.

## Targeted FASTQ Panel, 2026-10-22

Main integrated `9a2594fc` on `2ea67b1b`. Its recursive formatting check exposed
unresolved bare child-module paths under the path-loaded `rna_panel.rs`;
explicit `rna_panel/{input,tools,tests}.rs` paths now resolve. A new LF/CRLF
source guard fails before that repair and passes afterward. Recursive rustfmt
is not compiler/runtime acceptance; use the integrated source for the checks
below, retaining the original feature handoff's limits.

The owner's subsequent local compiler exposed E0425/E0433 for `Value` in
`tools.rs` and `tests.rs`. Repair `f37f7039` imports `serde_json::Value` in both
modules. Six source tests and recursive formatting pass, not a Rust check;
include that repair in the next ordinary RNA verification.

Changed areas: `src/engine/analysis/rna_panel.rs` and its `input.rs`, `tools.rs`,
`tests.rs`; `src/engine/protocol/rna_panel.rs`; thin enum/result/parser/MCP/
feedback/introspection wiring (including required `OpResult` literal fields),
legacy FASTA format refusal, tutorial 10.04/source/catalog/index, glossary/parity
and derived maturity source hashes. The only dependency change makes the
already-locked `md-5` version directly available. No external binary added.

Authored, not Rust-verified: `rna_panel_` covers independent source/member
identities, 20 loci / 200 ambiguous members without a 64-posting cap, both
orientations, junction flanks, strict FASTQ/gzip including truncation, quality
`>`, changed files, same-pass MD5/SHA-256, overlapping sources, cancellation,
file-boundary continuation, no overwrite, report binding and shared-shell/MCP
scientific payloads. Run normal checks at the integrated source, not a new
intermediate acceptance freeze:

```bash
cargo test --locked --lib rna_panel_
cargo test --locked --lib glossary_
cargo test --locked --lib multi_gene_sparse
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo run --locked --bin gentle_examples_docs -- tutorial-catalog-check
cargo run --locked --bin gentle_examples_docs -- tutorial-manifest-check
cargo check -q --locked --tests
```

Confirm the catalog/parity projections with the canonical Rust generators;
their Python/source checks are not generator execution. Follow [10.04](tutorial/10-04_targeted_nanopore_fastq_panel.md)
on synthetic data before any explicitly authorized real chunk. Native Windows
must exercise spaced paths, gzip failure and handle closure before publication.
For live inner-agent acceptance: each tutorial stage must use the shared route,
approve external writes/project loading, preserve compatibility sets and
support runtime inspection/cancellation without inventing completion. Native
responsiveness is unestablished; no new persistent job scheduler was added.
Any throughput/background run must identify exact source/policy/input hashes;
twenty-target capacity is a synthetic functional fixture, not a measured claim.
Do not process ongoing downloads, POD5 or private libraries without approval.

## Glossary And All-Features Follow-Up, 2026-10-22

The two repairs are independently verified at
`5ab0b02aef2c05edc15e8434be2cf1429ebf1ce5`, tree
`ac970706ac5ee4804b5022050f33479f15cb4871`; do not repeat the completed campaign
or transfer this scoped pass to later source. The previous failures were the
glossary fixture and E0521 at `70969635` / CI `39739876557`.

Retained owner-supplied archive:
`/tmp/gentle-two-repair-5ab0b02a-20261022.tar.gz` (128,815 bytes), SHA-256
`ed724dbe04adfe480872e2e072ba76a5f792f98fdc1979bdbd3ef7d7093dbed0`.
Codex validated all 26 manifested payloads in place; `SHA256SUMS` SHA-256 is
`5e50695983ef47ac6d1db6635306421fa35c7693b2f820951439d77a2be72ead`.
`REPORT.md`, `COMMANDS.md`, four command streams and CI JSON receipts preserve
the scope; binaries are identified by hash but not included.

- Linux: glossary 8/8; protocol UI 2/2, including seven widgets with one static
  step ID, distinct expected pseudonymous scopes and unchanged engine state.
- Linux: `cargo check --tests --features gui-test-support` and `--all-features`
  both returned 0. Checks compile targets; they do not execute tests.
- Toolchain/profile: Rust/Cargo 1.97.1, Debian forky x86_64, locked dev/test,
  one build job/test thread and repository-configured 16 MiB `RUST_MIN_STACK`.
  This is not default-stack evidence or native GUI acceptance.
- Same-source [CI 39774820033](https://github.com/smoe/gentle_rs/actions/runs/39774820033)
  passed both repaired gates on Windows/macOS. Windows completed all later
  gates; Linux headless and release safety passed. Full Linux was skipped.
- Overall CI remains **FAIL**: macOS later failed `Full test suite` (exit 101).
  Retained public logs returned HTTP 403 and disclose no failing test. Do not
  attribute it to either repair or claim the later provider-fixture fix
  resolves it without the named failure/results.

Next action is normal current-source CI/RNA verification plus the separate
provider-fixture follow-up above, not another glossary/protocol survey. The
RNA glossary filter covers newly added routes. Obtain the named macOS failure
from authorized logs if it remains reproducible. Disk-full warnings and a
failed final fetch are retained environment limits; confirm adequate free space
before Glen's next run, without cleaning a shared cache. A write-free remote
comparison established that the tested source was an ancestor of `1cfea9e3`,
not evidence for that descendant. No local Rust, push, dispatch or rerun was
performed when reconciling this return.

Preserve Glen's reported `70969635` passes separately: root test-target check,
protocol labware 11/11, renderer labware 2/2, root labware 16/16, rack mutation
2/2 and planning 20/20. Do not repeat those surveys solely to relabel them.
Legacy raw baselines, source-matched meshes and native carrier controls remain
open; missing OpenSCAD skips are not acceptance. The proposed commercial
product/offer/protocol slice and further bench models are separate work, not
part of these CI repairs. Historical procedures below remain source-scoped.

## Catalog-Parity Repair And Gel Return, 2026-10-20

Preserve Glen's [scoped gel pass and partial native return at bb1686bd](gel_workbench.md#integrated-linux-return-at-bb1686bd-2026-10-20).
The integration and normal [CI 39346038531](https://github.com/smoe/gentle_rs/actions/runs/39346038531)
remain failed. Both Linux/Windows primer-matrix logs show the same invalid
`catalog doctor` override assertion; do not invent a second primer defect.

After an ordinary owner-approved push, use a clean checkout and run-owned target
to verify the narrow canonical-reachability repair. Keep default GUI features
for the root test check, record the actual SHA/toolchain/platform/profile/stack
and retain complete unfiltered results:

```sh
cargo check -q --locked --tests
cargo test --locked -p gentle-protocol -j1 -- --test-threads=1
cargo test --locked --lib -j1 specificity_matrix_ -- --test-threads=1
cargo test --locked --lib -j1 introspection_distinguishes_drafts_from_validated_reports -- --test-threads=1
cargo test --locked --lib -j1 catalog_mcp_read_only_parity_and_import_confirmation -- --test-threads=1
cargo test --locked -j1 --test parity_matrix_freshness
```

Require normal Linux/Windows CI beyond the formerly failing matrix step;
later skips and a suite with the six poisoned-registry tests filtered are not
acceptance. No new interim source freeze or workflow dispatch is requested.
The repaired registry keeps all four catalog commands shell-reachable on
MCP/JS/Lua without fake engine operations or n/a exceptions; provider update
confirmation remains unchanged. Do not re-pin matrix counts merely to pass.

Review reconciliation, 2026-10-24: the repair is already in `91070ae4`, not a
new gel defect to implement. The [density rebase return](gel_workbench.md#comparison-density-verification-2026-10-20)
records all six protocol `gel_workbench` tests passing on `09639034` plus the
density slice, including the canonical parity projection. The earlier failure
above remains evidence for its older source, not today's gate result. That
filtered pass does not certify a current unfiltered protocol/CI run; keep the
commands above in normal consolidated verification, without repeating a full
campaign merely to relabel the SHA.

For gel native follow-up, actual project replacement, explicit selection
switch and separate Cancel-button behavior remain from Glen's `bb1686bd` return.
At the next consolidated native run, also check the later density UI explicitly:

- Select a saved comparison: summary and QQ must wait for **Verify and display
  saved comparison**, rather than appearing immediately. Successful replay
  reveals the shared profile panel, summary and QQ together. A deliberately
  mismatching copy must show the failure reason without a trusted curve,
  summary or QQ; raw CLI `inspect_comparison` remains diagnostic, not validation.
- Select **Signal per size width**, verify, then explicitly select **Legacy bin
  integrals** and verify again. Mode/selection/revision changes must remove the
  previous view; exercise project replacement, Cancel and close during work so
  no late result returns. Confirm the selector also drives PDF export. Use
  explicit density in the CLI for GUI/PDF comparison; omission remains legacy.
- Preserve saved integrals/statistics, zero versus gaps and the common-axis
  acquisition caveat. The separate unscaled custom-only reference is not an
  algebra coefficient or normalization. Real-data subtraction remains undecided.

Retain honest platform/X11/Xvfb versus physical-host scope and original captures;
green automated guards do not close these paths. Do not repeat completed gel
arithmetic, six-command 10.03, JSON/PDF or five-page legacy comparisons just
to relabel them; recheck affected paths if their source changes. Any native
capture requires the existing owner-approved scope, not private desktop data.
Return pass/fail/not-run, command logs, source/binary/input/artifact hashes and
locations through Steffen. No live model, provider contact, package or scientific
approval is requested. The carrier checks below remain separately outstanding.

## Carrier Follow-Up, 2026-10-20

After the next ordinary owner push, verify the new carrier changes without
repeating completed sequence-design/expression/gel campaigns. Record the source
SHA, toolchain/profile, lock/binary/input hashes, commands, pass/fail/not-run
and artifact locations. No new interim freeze, publication or scientific approval.
Use a clean checkout and explicit run-owned `CARGO_TARGET_DIR`.

```sh
cargo check -q --locked --tests
cargo test --locked -p gentle-protocol labware
cargo test --locked -p gentle-render labware
cargo test --locked --lib labware_ -- --test-threads=1
cargo test --locked --lib planning_wetlab_ -- --test-threads=1
cargo test --locked --lib rack_openscad_all_templates_preserve_bottom_floor_and_open_top -- --test-threads=1
cargo check -q --locked
python3 -m unittest scripts.test_labware_models -v
```

Keep default GUI features enabled for the test-target check: library-only
`cargo check` does not compile `src/app/tests.rs`. The missing carrier-template
import there is repaired, but the new source still needs external Rust checks.
Claude's reported 20/20 ranking run used a temporary import fix; it is not a
pass for the unpatched `28d964e6` or an acceptance claim for this branch.

Check both Safe-Lock support modes, exact 4x9 dimensions/A1/D9, refusal reasons,
generic-op/shared-shell/GUI labels, unverified disclosures and unchanged project
state. Verify the canonical parity projection: only template arguments change,
no operation or command rows are added. Retain native GUI tube/carrier export
control acceptance separately; authored GUI tests are not live acceptance.

The `labware_` filter also covers these new, unrun export regressions:

Claude reports three root `labware_` failures at `09639034`: the two tests
below and `labware_carriers_preserve_legacy_json_and_export_both_support_modes`
used no placements, but carrier labels require placement scope. Each fixture
now has one unoccupied A1 guard placement. Claude's temporary macOS probe passed
12/12 after this change and was reverted; this committed repair has not been
run here. Keep the engine's empty-scope refusal intact and rerun the root
`labware_` filter at the integrated SHA; no new acceptance is claimed.

- `labware_pcr96_carrier_exports_preserve_floor_and_disclose_nonclaims`: exactly
  96 floor-preserving/open-top cuts, fixed height/opening/floor and SCAD
  assumptions/nonclaims; each of the four SVG kinds has exactly one description
  and one visible unverified disclosure. Project state remains unchanged.
- `labware_legacy_carrier_exports_omit_drawing_disclosures`: all three legacy
  templates retain the original SCAD first line, none of the new SCAD comments
  and neither disclosure in any of the four SVG kinds. This guards disclosure
  leakage, not full byte identity; the historical baselines remain required.

The earlier `rack_openscad_all_templates_preserve_bottom_floor_and_open_top`
already resizes its rack to 8x12 and counts 96 PCR96 floor-preserving cuts.
Retain that coverage; the new test adds total-cut, header and SVG-count guards.
All new Rust checks remain Glen/CI work, not a locally executed pass.

Capture and commit the nine original carrier baselines from `1c5be317` source-
matched binaries: SCAD, simulation JSON and the embedded physical spec for each
legacy template on one fixed synthetic rack. Document exact recreation and raw
hashes; add scoped LF rules, LF/CRLF coverage in `scripts.test_tutorial_checkouts`
and Rust byte-comparison tests. Compare against new-source raw exports, without
normalizing or regenerating away differences. Glen's earlier 9/9 comparison is
historical evidence, not a committed fixture or this source's pass.

Render both new Safe-Lock carriers in OpenSCAD from real source-matched GENtle
exports. Require finite nonempty meshes, expected bounds, 2-mm floors and open
bores. Retain original SCAD/spec/mesh hashes and logs; optional tool skips are
not acceptance. Physical fit, cap clearance and thermal/material suitability
remain unverified even with successful mesh renders.

Include the conditional PCR slice, whose source gate is now met: Eppendorf
`0012 595 670-00/2026-05`, page 4, raw SHA-256
`67b84dea22bcb2cbea26907b8f1116c506210f1917744feb64b2fafb6745b40d`.
The [source record](../data/labware/README.md#pcr-02-ml-drawing-and-96-position-carrier)
explains maximum OD 6.0 versus main-cylinder 5.8 mm, unknown seating collar,
6.6-mm opening, 13-mm floor bore and checked nominal standard coordinates.
The same filters cover new footprint/A1/H12, 8x11/12x8/16x24 refusals,
no-overwrite, shell/generic-op and GUI labels. Run the renderer's complete
tests too (`cargo test --locked -p gentle-render`): its older source-binding
test name does not include `labware`. Verify clean `InvalidInput` when exporting
the PCR tube visualization (unknown collar), not a crash or invented dimension.
Render a third, PCR96 carrier mesh using an exact 8x12 grid. Keep all legacy
and Safe-Lock checks; inspect new-template SVG disclosures without inferring
full ANSI/SLAS compliance. No dense spacing, boxes, 0.1 ml profile or thermal
contact was implemented in the original carrier tranche. The offline mesh
harness uses literal parameters and
actual SCAD format strings, not a freshly built Rust exporter.

Owner requirements clarified 2026-10-20: offer **adjustable** compact spacing
below 20 mm in a separate extension without changing existing templates.
Storage prioritizes dense transport/freezer organization rather than treatment
access. Commercial paper/cardboard boxes should be modeled as possible storage
destinations, not presumed fabrication targets. Box dimensions/dividers/lids,
closed-cap clearance, fit and material/temperature qualification still need
their own evidence/contract. Box formats must be variable, not fixed to 9x9.
Compact spacing is now authored separately; no physical box model or new release
gate follows.

### Compact Carrier Follow-Up

Run the existing protocol/root `labware` filters and `cargo check --tests` on the
integrated commit, then canonical parity freshness and glossary/interface checks.
New authored tests (not run here):

- `labware_compact_pitch_bounds_preserve_conservative_edges`
- `labware_compact_pitch_decimal_and_legacy_rack_bytes_are_exact`
- `labware_compact_pitch_exports_geometry_persists_and_clears_without_drift`
- `labware_compact_pitch_refuses_invalid_modes_and_other_templates_without_writes`
- `labware_compact_pitch_shell_and_generic_operation_share_mutation_contract`
- `labware_compact_pitch_keeps_variable_grid_formats_independent_of_spacing`

Existing GUI-label coverage also checks pitch-control visibility. Review the
numeric field, exact mm/micrometre conversion, Clear, seating-mode refusal and
unsupported-template message natively. Source-matched OpenSCAD renders should
include one compact layout. Retain conservative/legacy byte baselines; none was
regenerated here. The parity Markdown projection was updated without running
its Rust generator: Glen/CI must confirm byte-exact canonical freshness.
The new operation is shared by GUI, CLI, generic MCP/scripting/workflow routes;
no agent auto-execution/approval exemption is added.

### Rack-Mutation Confirmation Follow-Up

Claude reports protocol 11/11 and render 2/2, but root `labware_` 15/16 at
`17f846eb`: `labware_compact_pitch_shell_and_generic_operation_share_mutation_contract`
correctly rejects the missing mutation classification. All eleven rack mutators,
not just compact pitch, were absent from `ShellCommand::is_state_mutating`, so
`--allow-auto-exec` could bypass the assistant's mutation-confirmation gate.
The repair adds exactly those eleven variants; no geometry, export strings,
non-rack classification or MCP tool permission changes. The shared assistant,
GUI and planner admission now use the corrected flag. Claude's reverted
267-test probe is not execution of this committed source or current acceptance.

New Rust tests, authored but not run here:

- `rack_mutators_are_state_mutating_and_block_agent_auto_exec`: parse all eleven
  commands; retain read-only `racks show`; block an actual auto-intent pitch
  suggestion with `AGENT_MUTATION_CONFIRMATION_REQUIRED` and unchanged rack,
  state, journal and revision; explicit index selection can still execute it.
- `rack_mutation_admission_cases_stay_in_sync_with_executor_dispatch`: bind the
  per-variant cases and classifier to the actual rack mutation handler. This
  source-set guard plus adjacent comments catches a future dispatch-list drift;
  it does not substitute for the execution regression above.

After ordinary integration, use the same clean checkout/run-owned target rules
above, default GUI features and the exact tested SHA. In addition to the
protocol/render and carrier checks, run the **full root library suite**, not
only the `labware_` filter:

```sh
cargo check -q --locked --tests
cargo test --locked --lib rack_mutat -- --test-threads=1
cargo test --locked --lib labware_ -- --test-threads=1
cargo test --locked --lib -- --test-threads=1
cargo check -q --locked
```

Keep agent-parity/mutation/GUI results and primary failures separate from any
mutex-poisoning cascade. The source-only Python LF/CRLF membership/negative
control checks are not Rust, live-agent or native acceptance. OpenSCAD renders,
nine legacy baseline files and native GUI review remain outstanding. Include
11.75..12.049 mm: setter admission uses the collar limit; floor-seated export
must clearly refuse it. Optional positional disclosure cleanup was deferred.

Export classification inventory (unchanged): all seven named rack SVG/SCAD/JSON
exports return `is_state_mutating=false`, as do named `SaveProject`, `ExportPool`
and `LaddersExport`. Named guide oligo/protocol exports return `true`, as does
the conservative generic `op`/`workflow` route. Thus rack exports match some
file-writing routes, not a uniform export policy. `false` means no declared
project mutation here, not absence of filesystem side effects; no broader
export-policy repair is part of this commit.

The clarified future box requirement is floor-supported 1.5/2.0 ml tubes with
snug divider contact in openings that are square in top view. Cells are
separated by wall thickness; literal tube-to-tube contact is not required.
The [Semadeni example](../data/labware/README.md#different-storage-box-formats)
supplies a source lead, not interior/divider/lid dimensions or Safe-Lock fit.
Variable formats and all physical qualification remain separate future work.

## Historical Covering Request, 2026-10-19

Updated 2026-10-19: after his push, Steffen should forward the
**[concise parity/expression follow-up](glen_parity_expression_followup_20261019.md)**.
Only that message is active for this request. The execution checklists and
earlier capture requests below are retained background, not renewed authority
or a request to repeat completed checks. No release freeze is selected.

The original covering request was prepared 2026-10-18. Its source-bound return
and historical procedures follow; continue development independently.

## Returned Results And Remaining Request

Glen executed the pre-inspector checks below at
`394dd5d0523ea91c19a8eb94d5ed8ad7bd5299b2`; see the
[verified evidence and six-item status](glen_audit_followups_20261008.md#focused-verification-at-394dd5d-2026-10-18).
README 4/4, sequence-design Rust 67/67 plus actual CLI suite 22/22, starter
suite 6/6, workflows 95/95 and rack/labware 127/127 pass. The nine historical
carrier outputs are byteidentical. Keep those checks closed at that source;
do not reissue this entire original checklist as unfinished work.

One reported defect is repaired locally: the parity matrix's two missing
tube-operation rows/counts. Codex applies Glen's canonical generated delta and
preserves newer inspector surfacing without running local Rust. On the next
ordinary owner-pushed revision, run `parity-matrix-check` with its source-matched
docs binary, `cargo test --locked -j1 --test parity_matrix_freshness`, and inspect
normal CI. Run the seven inspector regressions and `gentle-shell` UI-intent
checks in original item 3: `5893ebb0` postdates the tested source. Its updated
09.03 tutorial/catalog projection also needs source-matched freshness checking.
No new bundle or intermediate freeze is requested; final-candidate affected
checks remain on the normal release path. OpenSCAD mesh execution remains
`not_run` on Glen's host, not a failed Rust exporter. Native, agent, scroll and
package boundaries below receive no new acceptance from this return.

## Source And Isolation

Fetch published `smoe/gentle_rs main` once and record the resolved full SHA and
tree. Test that source in a clean detached checkout with an explicitly set,
run-owned absolute `CARGO_TARGET_DIR` outside the checkout. Record toolchain,
platform, features/profile and raw lockfile SHA-256. Use only binaries built
from this checkout; preserve the primary checkout, shared caches and historical
evidence. A later main commit does not invalidate this scoped report or require
another synchronization request, but it is not automatically accepted by it.

No push, workflow dispatch, tag, release publication, provider contact, live model
or package/container build is requested. The owner's follow-up additionally
requests only the four synthetic inspector tutorial screenshots described below;
no unrelated GUI/desktop/private-data capture is authorized. Existing normal CI
can be inspected without restarting it. Keep `package-opt1-fallback` separate;
the selected .13 packaging recipe remains opt2. Codex runs no local Rust.

## Original Focused Checks And Inspector Follow-Up

1. **README metadata repair:** verify the `57a2403a` repair is present. Run the
   release-version regression below and inspect the new normal push CI when
   available. Historical run `39066178748` at `1c5be317` failed the README
   label guard on Linux and Windows; it did not test the repair. Both root
   libraries and all five stripped-child families passed at that older SHA.
   Preserve those bounded results rather than calling the failed run green.

   ```sh
   cargo test --locked -j1 --test release_version_consistency
   ```

2. **Sequence-design follow-up:** verify `8a0e2670` and `958b2238` are included.
   The first retains a freshly validated best after edited-candidate match-cap
   overflow, marking search incomplete without optimality/infeasibility proofs
   instead of returning Invalid and discarding it. The second adds independent
   180 single-symbol and 10,800 paired-symbol IUPAC cases without changing the
   production matcher. Eighteen offline oracle/harness passes are not Rust
   acceptance; execution subsequently passed at `394dd5d` as recorded above.
   Execute the core/protocol tests that the failed older workspace run did not
   reach, plus the shared engine and CLI boundaries:

   ```sh
   cargo test --locked -p gentle-sequence-design -j1
   cargo test --locked -p gentle-protocol -j1 sequence_design
   cargo test --locked --lib --no-default-features -j1 sequence_design -- --test-threads=1
   cargo test --locked --bin gentle_cli --no-default-features -j1 sequence_design -- --test-threads=1
   cargo check -q --locked -j1
   cargo build --locked --no-default-features -j1 --bin gentle_cli --bin gentle_examples_docs
   ```

   Replay tutorial 06.07 and the four actual CLI acceptance tests using the
   [functional procedure](sequence_design_acceptance.md#exact-candidate-functional-checks).
   For these functional checks the dev binaries just built are sufficient;
   record that profile, without claiming package or performance acceptance.
   Set `GENTLE_TUTORIAL_BIN_DIR` to the absolute run-owned target's `debug`
   directory, `GENTLE_SEQUENCE_DESIGN_CANDIDATE` to the recorded full HEAD, and
   `GENTLE_SEQUENCE_DESIGN_EVIDENCE_DIR` to a fresh absolute evidence directory:

   ```sh
   python3 -m unittest scripts.test_sequence_design_audit_oracle scripts.test_sequence_design_tutorial -v
   ```

   All four actual CLI tests must execute; skips and the four Rust CLI unit
   tests are not substitutes. Check retained-best behavior after match-cap
   overflow, incomplete status without false proof flags, cancellation clearing
   output, lexical ties and independent IUPAC coordinates/strands/palindromes.
   Preserve both strategies' joint motif/global/window-GC, tampering,
   materialization, save/reopen and same-session undo/redo checks. Historical
   `de9df8af` prototype acceptance stays separately attributed.

3. **Integrated .13 expression/storage/tutorial paths:** follow the focused
   commands, generator freshness/parity checks and two real-CLI starter replays
   in [the preparation handoff](internal_13_preparation_handoff.md#external-completion-before-candidate-selection).
   Its owner dispatch/publication section is not part of this request. Use the
   same freshly built binaries with `scripts.test_release_workflow_tutorials`;
   no optional-binary skip is acceptance. Retain LF/CRLF and quoted-path results,
   deterministic container order, conservative expression refusals and unchanged
   review state/files. Report drift; do not regenerate historical evidence.

   The merged `gentle_rs_3` inspector at `5893ebb03b4915ba1af2f77f1d7681904266b6f9`
   supersedes the older preparation note deferring all inspector work. Include:

   ```sh
   cargo test --locked --lib -j1 protein_expression_inspector_ -- --test-threads=1
   cargo test --locked -p gentle-shell -j1 ui_intent
   cargo test --locked --lib -j1 mcp_server::protein_expression::tests:: -- --test-threads=1
   cargo test --locked --lib -j1 protein_expression_inner_agent_suggestion_ -- --test-threads=1
   cargo test --locked --lib -j1 protein_expression_agent_card_ -- --test-threads=1
   cargo test --locked --lib --features gui-test-support -j1 protein_expression_capture_targets_ -- --test-threads=1
   ```

   Require all seven inspector regressions to execute. Verify shared-report
   parity, explicit/absent inputs, typed withholding, copy-only suggestions,
   stale/closed-worker rejection and unchanged project/history/files. Check
   canonical tutorial/catalog/parity freshness at this source; static metadata
   edits are not generator acceptance. Native window review stays at release
   preparation. The later owner request additionally asks for the four
   [09.03 screenshot checkpoints](tutorial/09-03_protein_expression_handoff_cli.md#required-tutorial-screenshots)
   from a verified binary and clean synthetic profile. Four manual macOS PNGs
   at `9f01098b` are now retained with binary/input/image hashes and limits in
   the linked archive. They do not contain the uncompiled parity follow-up.
   Review them; retain typed GUI reports and semantic snapshots at the integrated
   candidate separately, without treating these images as automated acceptance.
   No provider/model contact, unrelated capture or automatic capture is allowed.
   Also execute the named MCP `protein_expression_handoff` discovery/report/
   refusal regressions, the deterministic inner-agent command comparison and
   the `gui-test-support` status/locator test. The new control card is guidance,
   not evidence of a model's behavior; compare actual reports on saved state.

4. **Merged rack safety and tube reference models:** verify the integrated
   `gentle_rs_2_main` commit `5b1da7d3981e7a2a4280cadc4b86b9345d82896c`.
   Run its protocol/renderer, headless engine/shell and default GUI regressions,
   then existing rack export and capability/parity checks:

   ```sh
   cargo test --locked -p gentle-protocol -j1 labware::tests::
   cargo test --locked -p gentle-render -j1 labware::tests::
   cargo test --locked --lib --no-default-features -j1 labware_ -- --test-threads=1
   cargo test --locked --lib -j1 labware_ -- --test-threads=1
   cargo test --locked --lib -j1 rack_ -- --test-threads=1
   python3 -m unittest scripts.test_labware_models -v
   ```

   Require positive, overflow-checked admission through 1,536 slots, refusal
   above the limit, unchanged unsupported saved racks on reopen, and refusal
   before overwriting any existing export. Exercise generic `op` inspection
   and all three tube exports with quoted paths: deterministic SCAD, legacy
   result decoding, explicit artifact-write confirmation and unchanged project
   state/history. Generator/parity checks must include the new operation IDs.
   OpenSCAD mesh checks alone do not execute the Rust exporter; missing
   OpenSCAD is a recorded skip, not mesh acceptance. Keep source URL/digest/
   revision/page and exact/approximate/unknown dimensions intact.

   Follow [the provenance and legacy-baseline boundary](../data/labware/README.md#verification-boundary).
   In separate isolated checkouts, compare raw SCAD, simulation JSON and spec
   JSON for the same fixed rack and all three existing carrier templates,
   before this feature at `1c5be3175323eec41359b980b5a41e22c7edaec1` and at
   the tested integrated source. Retain input and before/after hashes outside
   Git; no source-only or template-only result proves byte identity. Native
   Rack controls remain a release-stage check, not new capture authority here.
   Tube models are body/collar visualizations with cap/hinge/cavity omitted,
   not printable vessels or fit evidence. New carriers, support/pitch,
   96-position interpretation and storage boxes remain deferred owner choices.

Use `--locked` throughout. If dependencies are already cached, `--offline` may
also be used; otherwise record retrieval failures as environmental. Do not
change the lockfile or widen resource limits silently to obtain a pass. Missing
tools, zero selected tests, timeouts and unrun platforms remain explicit.

## Remaining Six-Gate Boundaries

Use the [six-item ledger](glen_audit_followups_20261008.md) as the status owner,
and the [release-stage request](glen_final_candidate_execution_handoff.md) only
when release preparation is authorized. Audit preservation, 08.04 contracts,
live-agent acceptance, honest visual evidence, scroll benchmarks and package
acceptance retain their individual evidence boundaries. These focused checks
do not close all six. Existing tutorial/link repairs and historical audits
must not be reissued as unfinished work. Native/physical-host, model, scroll,
package and scientific verdicts stay separate; this request authorizes none
of their additional execution.

## How To Return Results

Send **one concise response to Steffen for forwarding into the main Codex
chat**. Do not rely on direct Glen-to-Codex messaging or assume a path on one
host exists on the other. Include:

- Tested full SHA/tree, platform, Rust/Cargo versions, features/profile,
  lockfile and binary hashes, and any same-source CI URL/conclusion.
- Per-check `pass`, `fail`, `not_run` or `blocked`, actual selected test counts,
  skips, exact commands/exit codes and primary failure versus cascades.
- Remaining six-gate boundaries, distinguishing new results from historical
  passes and software execution from native/scientific/release acceptance.
- Evidence location and SHA-256 of its manifest, plus any concrete defect and
  minimal proposed repair. Do not modify or push source as part of this check.

Supply a small folder/archive containing `REPORT.md`, command/exit-code logs,
raw stdout/stderr, CI/CLI receipts and a per-file SHA-256 manifest. Steffen can
copy it to `/tmp/glen-verification-<short-sha>-<date>/` on the Codex host and
relay the path and manifest hash here. Keep bulky binaries, raw captures and
reference inputs on Glen's host with hashes/locations; transfer them only on
request. Never include credentials or private projects. Record content hashes
as identity checks, not authenticated provenance or human scientific approval.

## Preserved Prompt Records: Not Active Instructions

The following received files are committed verbatim for history:

- [Release-link repair prompt](release_handoff_link_repair_prompt.md): repair
  already closed; its preparation-time failure and old runner/version status
  are historical, not a request to restore the deleted .12 planning ledger.
- [Sequence-design closure prompt](sequence_design_closure_codex_prompt_20261014.md):
  closure was reconciled in the plan/acceptance records. The current
  [compiler decision](decisions.md#workspace-compiler-compatibility) supersedes
  its broader proposed workspace minimum/publication gate; select no MSRV.
  The [issued handoff's SHA correction](glen_sequence_design_acceptance_handoff_20261014.md)
  remains intact.
- [Original Claude VKORC1 prompt](vkorc1_acceptance_codex_prompt_20261014.original-claude.md):
  retain the original separately from the [reconciled task](vkorc1_acceptance_codex_prompt_20261014.md).
  Its in-progress CI, undecided .13 aim and other snapshot claims are not
  current status. `--offline` is not a substitute for `--locked`, nor proof
  that a lockfile was rewritten. No old embedded dispatch is authorized here.
