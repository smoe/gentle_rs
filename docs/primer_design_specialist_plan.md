# Primer Design Specialist: Current Status And Remaining Scope

Slice B was authorized on 2026-11-05 for implementation and synthetic evaluation.
The owner authorized only tail-free selected-pair preparation into simple PCR on
2026-11-19; see the implementation ledger below. The other specialist proposals
are not implementation authorization or scientific acceptance.

## Verified Starting Point

Engine and shared-shell support for primer-pair design is already available via
`DesignPrimerPairs` plus report inspection/export (`primers design`,
`primers list-reports`, `primers show-report`, `primers export-report`).

The shared UI destination is `ui open|focus|close pcr-design`; a headless
acknowledgement does not establish that a native window opened.

The dedicated PCR Designer already provides pair-design forms, alternative
pair/report inspection, backend provenance, Primer3 preflight and explicit PCR
handoffs. The transcript-panel and gene-assay-study workspaces are additional,
narrower entry points. See [the GUI guide](gui.md) and
[PATZ1 tutorial 04.08](tutorial/04-08_gene_assay_study_gui.md). Their existence
does not establish native acceptance at the current candidate.

The remaining specialist scope is a general report-centric entry point,
persisted filtered views and the separately proposed single-primer pane. Do
not implement another pair-design engine or duplicate the PCR Designer just
to satisfy this historical plan. First reconcile each proposed control with
its existing owner.

Report navigation is not absent: `GENtleApp::open_lineage_analysis_artifact`
already routes primer, specificity, qPCR, RNA-read and other saved artifacts
to their existing viewers. Lineage offers dedicated open actions, and PCR
Designer lists/shows/exports its reports. A new all-domain report dashboard or
an untyped global `REPORT_ID` is not the starting requirement.

Whole-reference specificity is a separate evidence workflow. Source-bound
intended-target geometry and companion per-dimension reference-selection
receipts now prevent a newer report on another database from silently replacing
the active evidence. Legacy unbound evidence is inspectable, not passing.
The standalone multi-reference handoff/import/show/list core and
[tutorial 04.09](tutorial/04-09_reference_bound_primer_specificity.md) are
verified for the bounded Linux pair cases in Glen's report below. The standalone
pair contract does not attach anything to order readiness. Panel orchestration
is a separate operation; automatic tail selection and engine-native single-primer
design remain separate.

Read-only `primers transcript-assay-specificity-matrix PANEL_REPORT_ID` now
projects retained panel receipts with separate design/reference/policy scopes,
binding diagnostics, active/superseded history and unlinked same-oligo citations.
PCR Designer now offers the same read-only review in place; Rust and native checks remain with Glen/CI.
This is inspection only, not panel fan-out,
new handoff/import machinery, reference selection or readiness attachment.

## Conditional Slice B: Bound Panel Fan-Out

Status: implemented for synthetic evaluation, with historical Linux execution
and canonical parity-repair receipts recorded in the Glen handoff. This is not
blanket current-candidate, native or real-reference panel acceptance. No dedicated
GUI, new BLAST execution, changes to the matrix or readiness attachment are included.

### Motivation And Entry Conditions

A retained matrix answers what evidence exists; a single-reference panel handoff
cannot bind an inventory across multiple references for one atomic cross-reference import.
The owner authorized this bounded composition after [Glen's phase attribution](glen_next_verification_handoff.md#post-blast-phase-trace-return-2026-11-05)
at `dfca3261` and the [Linux reader verification](audits/glen_followups_20261008/reader_verify_e4bc7386/REPORT.md)
at `e4bc7386`. All four single-pair imports completed and reopened. TP53
untraced took 58.26 s of a 60 s budget; peak RSS reached about 1.08 GiB.
That is neither a measured panel throughput nor a guaranteed speedup ratio.
Do not repeat that campaign or infer native Windows or scientific acceptance.

### Implemented Contract

Prepare one existing child per available reference covering all selected assays.
Bind current design, source, full pairs/tails, derived intended geometry, common
policy and exact command inventory. A scientific ID is independent of the bundle
location; a separate transport digest locks paths and commands. Preflight required
indexes before reserving a fresh directory. Optional unavailable references remain
explicit. No inferred liftover, caller geometry override or external execution.

Import checks all command rows (including pending/failed/cancelled), validates
successful retained outputs before interpretation, and uses the existing pair-level
interpreter, reference freshness checks and per-assay dimension rules. One store
commit retains standalone diagnostics, never the panel or its selection receipts.
Do not call `record_primer_specificity_reference_selection` or finalize children.
Readiness attachment needs a later policy decision.

### Admission And Verification

Owner-selected draft ceilings **for synthetic use only**: 1-24 assays, 1-8
references (at least one required), at most 384 commands (two primers per assay
per prepared reference), 128 MiB per output and **128 MiB total across the entire parent**,
25,000 nonblank TSV rows total. Effective limits travel with the bundle; tighter
recorded ceilings apply. Check declared byte totals before touching outputs:
unknown sizes are not zero. The importer works serially; input ceilings do not
bound RSS, reference decoding, alignment/join work or run time. Real-reference
panel use requires separate resource measurements and approval. Do not add
cross-assay reference caching on the strength of the single-pair trace.

See [the protocol](protocol.md#panel-multi-reference-specificity), [tutorial
04.09](tutorial/04-09_reference_bound_primer_specificity.md) and the
[verification handoff](glen_next_verification_handoff.md#slice-b-panel-multi-reference-specificity).
Glen supplied canonical parity generation/checking and scoped Rust results at
`27a8de33`; the generated repair was imported. Do not request that historical
Linux campaign again as missing. Normal CI checks affected integrated source;
unproved native scenarios and real-reference panel runtime/RSS remain separate.
The old single-pair phase-trace/reader blocker is resolved, not a reason to defer
all Slice B verification. No production-capacity or biological pass is implied.

## Remaining Proposal: Persisted Filtered Views

Session selection and display controls are not a portable, named filtered-view
contract. Such a contract is still proposed, not registered as `primers view-*`
commands. Do not advertise those commands or infer that exporting an existing
report implements them.

Before implementation, define an engine-owned view referencing the exact report
ID/content hash, filter criteria, stable ordering/tie-breaks and selected pair
ranks. A view must not alter scores or scientific gate outcomes. Retain the
unfiltered report and make stale report bindings explicit. GUI, shared shell,
MCP and other intended adapters must consume the same view and export contract.
No second design algorithm or adapter-local ranking is warranted.

Required deterministic tests for that future slice:

- View create/inspect/export/delete and project save/reopen preserve the same
  criteria, report binding and selected ranks.
- Invalid criteria and stale report hashes fail consistently across adapters.
- Equal-score ties have a stable ordering; a filtered export states its subset
  without claiming completeness of the original search.
- Explicit selected-pair handoff preserves the original backend and evidence;
  filtering never turns generated assays into order-ready ones.

Fresh native report browsing and handoff checks remain acceptance work, separate
from these proposed commands. Engine-native single-primer design and broader
specificity orchestration are outside this residual plan.

## Other Remaining Scope (Requires Separate Approval)

- Reuse the PCR Designer and shared report paths for a general specialist entry
  point, only where a verified presentation gap remains.
- Support alternative-result filtering with persisted filtered views.
- Add a UI-derived single-primer inspection pane (derived from pair reports).
- Add explicit handoff actions into PCR, PCR Advanced, and PCR Mutagenesis
  workflows.
- Preserve UI-intent and shared-shell parity expectations for agent/MCP/CLI
  pathways.

## Proposed Next Review Slice, 2026-11-16

The owner forwarded Mistral's seven-component proposal for roadmap discussion.
Keep its UX ideas, not its stale phase-trace blocker, blanket pending Slice B
status, unverified absence claims or day estimates. Existing source and the
historical receipts above are the starting evidence. This is planning, not
authorization for a new implementation or release gate.

**First: audit selected-pair preparation for the three existing PCR forms.**
Inventory basic PCR, PCR Advanced and PCR Mutagenesis separately: existing
control, source report/pair/template binding, destination fields, missing inputs,
overwrite/approval behavior and tests. Reuse `render_primer_design_report_preview`,
the existing report stores and PCR operation forms. Record unknowns as unknown,
not as implemented or absent. If an existing control already provides the
behavior, prefer navigation/documentation repair over another command or button.

Do not conflate three different actions: opening a saved report, building
`BuildExperimentalAssayHandoff` readiness cards from a transcript panel, and
preparing a selected pair in a PCR execution form. The first two already have
GUI paths; their existence does not prove the third for every mode.

If the audit proves a gap, propose one bounded report-to-existing-form flow:

- Preserve exact report identity/content binding, selected pair rank, template,
  backend/evidence, complete oligos/tails and annealing geometry. Do not silently
  strip tails or guess advanced-PCR options to fit a simpler form.
- Preparation is not PCR execution, readiness approval, project/history mutation
  or output publication. Explain missing/stale evidence without rescoring it;
  require an explicit decision before replacing edited destination fields.
- Mutagenesis still requires explicit position/reference/alternate inputs.
  Selecting a pair must not invent a mutation or silently change PCR defaults.
- Route any shared preparation contract through existing UI-intent/shell
  boundaries; adapters must not acquire their own biological transformations.
  New shell spellings are proposals until implemented and registered.
- Require deterministic synthetic cases for exact selection and provenance,
  stale/wrong-template/missing-pair refusal, tails and missing options, edited-form
  protection, unchanged project/history and no implicit operation execution.
  Shared-route/quoted-path and native GUI checks follow the actual change;
  source inspection alone is not their acceptance.

**Second: assess a primer-scoped report locator**, reusing Lineage, report stores
and existing viewers before considering an all-domain hub. Define typed report
ownership, ambiguous IDs, missing sources and stale bindings before proposing
`ui open report ...`. Do not promise arbitrary file-backed RNA reports are in
the project store or silently import/open unrelated files.

Persisted filtered views and the pair-derived single-primer inspector remain
separate proposals. They are not independent primer design, score changes or
order-readiness promotion. Priorities here are provisional; no effort estimate
is justified before the path audit.

A bounded Mistral Vibe critique was attempted with all tools disabled and only
selected public source facts; startup failed before a review result. This plan
is Codex-authored, not Mistral-reviewed or newly Claude-reviewed. Optional Claude
read-only consultation was offered; it has not been invoked.

## Simple PCR Pair Preparation Ledger, 2026-11-19

Implemented the owner's first bounded slice, not the global report hub or the
three-mode Mistral draft. A button in the existing saved primer-report preview
uses one read-only engine helper to validate the exact saved/displayed report,
explicit rank, tail-free oligos and current template binding sites. Reverse
oligos are copied unchanged. Source report/backend and preparation-time digests
remain visible; no report schema, Operation or shell/MCP command was added.

The controller fills only simple-PCR forward/reverse fields, opens/scrolls the
existing form and keeps execution separate. Different occupied fields require
an explicit replacement decision. Confirmation is bound to the engine instance,
source snapshot, whole preparation-time template and previous destination values.
Invalidation/cancellation does not change fields. The action does not call the
UI metadata writer, apply an operation, create a product or publish a file.
Advanced/Mutagenesis inputs and the existing Unique policy are untouched.

The read-only audit found a shared Annealing length in the Advanced/Mutagenesis
forms and independent per-primer lengths in the engine. General tailed-pair
handoff therefore remains a separate slice, as do explicit mutation selection
and the later primer-scoped report locator. Do not silently populate historical
mutation defaults as if the user chose them. No new scientific logic is needed.

Legacy primer reports contain a template ID, not its original sequence digest.
Current coordinate/oligo validation and the new confirmation snapshot do not
retroactively authenticate an old whole-template or whole-reference assessment.
This control makes no readiness or laboratory acceptance claim.

Sixteen synthetic Rust regressions (six engine, ten GUI/controller) cover
selection, orientation, tails/geometry, source/destination/project drift,
confirmation/cancellation, idempotence, unchanged project/History and paint-only
behavior. They are authored, **not run locally**. Fast source guards and parsing
are not type checking. Use the
[remaining verification](glen_next_verification_handoff.md#simple-pcr-selected-pair-preparation-2026-11-19)
for exact-source Rust and native GUI/Windows/macOS checks. No new Claude review,
local Rust build, push or publication is claimed here.

Local recording, 2026-11-20: 12 Python source/tutorial cases, documentation
status, new-module formatting, parent-module syntax parsing and whitespace pass.
The initial source guard's whitespace-sensitive call match was corrected for
rustfmt line breaks; no production check was relaxed. Rust/native results remain
unrun. Session-close has four OK, two intentional warnings and no failures.

### Compile Repair After Glen's Return, 2026-11-20

Glen's two source-bound Linux returns at `ccd839dc` establish default GUI/test
compile failures, not runtime failures or acceptance. The outer and pending UI
state lacked required `Clone, Debug`, and the paint test used egui's removed
`Context::run`; the engine-test wildcard import was unused. The headless library
check passed, but zero default Rust tests and no native GUI checks executed.

Add the derives, remove only the redundant import and use the current `run_ui`
root-Ui callback with explicit disposal of unused texture deltas. A new Rust
regression checks cloned pending-state independence, trait bounds and unchanged
project content. There are now 17 authored tests (six engine, eleven GUI), plus
four fast source guards. Thirteen scoped Python source/tutorial cases and scoped
format/documentation checks pass; none establishes repaired-source compilation.
Keep the [focused compile/native handoff](glen_next_verification_handoff.md#simple-pcr-compile-repair-2026-11-20)
open for Glen/CI without replaying completed unrelated campaigns. No new biological
policy, command, report schema, termination change, local Rust build or push.
