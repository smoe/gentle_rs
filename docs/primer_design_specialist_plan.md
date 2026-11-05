# Primer Design Specialist: Current Status And Remaining Scope

Slice B was authorized on 2026-11-05 for implementation and synthetic evaluation.
The remaining specialist proposals are not implementation authorization or
scientific acceptance.

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

Status: implemented for synthetic evaluation, pending Glen/CI verification of the
new source. No dedicated GUI, new BLAST execution, changes to the matrix or
readiness attachment are included.

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
Canonical parity regeneration and checking, Rust execution and Windows behavior
remain with Glen/CI. No production-capacity or biological pass is implied.

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
