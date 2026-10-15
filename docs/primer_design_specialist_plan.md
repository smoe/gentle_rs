# Primer Design Specialist: Current Status And Remaining Scope

Documentation reconciled 2026-10-15; preserve the 2026-10-13 source and conditional
fan-out plan below. This residual plan is not implementation authorization or
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
implemented pending Glen/CI Rust verification. This core does not add panel fan-out or
automatically attach its aggregate to order readiness. Those integrations,
automatic tail selection and engine-native single-primer design remain separate.

Read-only `primers transcript-assay-specificity-matrix PANEL_REPORT_ID` now
projects retained panel receipts with separate design/reference/policy scopes,
binding diagnostics, active/superseded history and unlinked same-oligo citations.
Rust execution remains with Glen/CI. This is inspection only, not panel fan-out,
new handoff/import machinery, reference selection or readiness attachment.

## Conditional Slice B: Bound Panel Fan-Out

Status: retained proposal, not implementation authorization, a `.12` gate or
an acceptance claim. Slice A landed at `dee6033d`; its Rust/native/real-reference
verification remains with Glen/CI. The owner must separately approve B after
reconciling this outline with current source; offer a read-only Claude review.

### Motivation And Entry Conditions

Slice A answers "what evidence have we retained?", but neither requests nor
executes the searches still needed for a panel. With several isoform assays and
genomic/cDNA references, repeated manual preparation/import can omit one pair,
mix references or import results against a changed design. B would give the
external runner one bounded, reproducible inventory of required work and the
reviewer one assay-by-reference result with inspectable failures and gaps.
It improves orchestration, not primer design, biological specificity or ordering.

1. Have Glen/CI verify A and the existing pair/panel specificity contracts at
   one frozen candidate, without borrowing nearby revisions' acceptance.
2. Complete phase attribution through the main-maintaining thread, following
   the [bounded real-reference import profile](glen_post_blast_import_profiling_f35a58e6.md)
   and its [retained measurements](audits/glen_followups_20261008/post_blast_profile_f35a58e6/README.md).
   Both imports timed out; the sampled trace had zero prepared FASTA/FAI opens,
   so repeated FASTA-index loading is not established as the cause. Repair a
   demonstrated bottleneck before multiplying the workload; these measurements
   do not authorize a speculative optimization or identify a dominant phase.
3. Confirm that retained-evidence inspection and existing single-reference panel
   handoffs do not meet the actual batching need. If they do, do not build B.
   Otherwise approve its final contract and numerical admission limits separately.

### Proposed Minimal Composition

- Reuse `TranscriptAssayPanelSpecificityHandoff` children: one saved panel,
  one explicit reference and every selected assay per child. Add only the parent
  binding/manifest needed for atomic cross-reference import and shared limits.
  Do not recreate existing pair-level verdicts, applicability or interpretation.
- Bind the normalized panel design, source snapshot, assay identities, full
  oligos, annealing/tail boundaries, exact references, common policy and commands.
  Preserve derived intended-product geometry for its own reference/target space;
  additional caller mappings cannot override it or transfer coordinates by alias.
  Reject a hidden single-target selector in the common policy; record explicit
  per-child target normalization rather than silently changing that policy.
- Preflight all required references before publishing a fresh bundle. Optional
  unavailable references stay explicit. No download, indexing, BLAST execution,
  scheduler, macro loop or generic collection abstraction is introduced.
- Extract non-persisting validation/result construction from existing import
  seams. Repeated calls to persisting pair import or panel finalization are not
  atomic composition. Validate all retained bytes/bindings first, then perform
  one guarded store commit; tampering or late source conflicts leave no writes.
- Keep partial/cancelled batches structurally separate from panel acceptance:
  completed evidence may be retained as immutable standalone diagnostics, never
  through `record_primer_specificity_reference_selection` or the state-writing
  `finalize_transcript_assay_panel_specificity_handoff`. Active receipts and
  readiness remain unchanged. Readiness attachment needs a later policy decision.
- Reuse required/optional and genomic/transcriptome verdict semantics. Required
  failure dominates, missing/stale/incomplete required evidence cannot pass, and
  optional biological failures remain visible without becoming required failures.
  No universal specificity, isoform discrimination, oligo-QC or order claim.

### Limits And Acceptance To Set Before Implementation

Reuse the existing 1-8 reference admission and per-output byte ceiling through
shared constants rather than duplicated literals. Proposed initial cumulative
TSV budget: **128 MiB total across the entire parent**, derived from the current
pair-level `MULTI_MAX_OUTPUT_BYTES`, not 128 MiB multiplied by every command.
Check declared and actual retained-byte totals before allocation/interpretation;
unknown sizes are not zero and exceeding a limit is not accepted truncation.
The earlier 32-assay / 512-command proposal is not an approved limit. Choose
assay, command, parsed-HSP and memory/work caps from bounded profile evidence,
record exact limits in the effective request, and refuse excess before work.

Author synthetic tests for omitted/duplicate cells, common-policy conflicts,
source/pair/tail/reference replacement, intended geometry on both target spaces,
required/optional verdicts, partial execution, cumulative budgets, manifest/hash
tampering, LF/CRLF and quoted paths, stable IDs, save/reload and atomic rollback.
Explicitly assert unchanged active selections/readiness and no partial panel
attachment. All adapters use the shared engine; no dedicated GUI is in this slice.
Extend tutorial 04.09 only with implemented routes after approval. Glen/CI owns
Rust execution and exact-candidate acceptance; no private report regeneration.

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
