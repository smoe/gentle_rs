# DNA-Centred Regulatory Evidence And Reporter Planning

Status: discussion draft, not approved for implementation.
Initial source review: `d978b0c4368071f839ceecf8d494212ad9b992ec`, 2026-09-13.
Review follow-up: `e7ca9ad133326aae68fa08481e9e845a9cce62af`, 2026-09-13.
Claude feedback: supplied by the user and checked against the source. Its
substantive findings and Codex's revised scope are recorded below. The earlier
local CLI attempt failed authentication; it did not produce this review.
The revised plan has not been approved by the user or re-reviewed by Claude.
This draft does not add requirements to the unreleased `.10` acceptance gate.

## Recommendation

The DNA viewer should be the GUI's starting point and coordinate reference.
Expose the existing locus evidence inspector directly from it, then connect
that inspector to the existing Promoter design planners. Keep the same entry
points discoverable and invokable through the inner agent's shared UI intents.

The user journey is:

```text
DNA locus + exact reference
    -> inspect supplied evidence and annotated TSS/transcript choices
    -> select a saved genomic region OR a transcript-aware architecture
    -> compare reporter hypotheses in Promoter design
    -> review an exact product proposal
    -> separately approve creation of designed sequences
```

The core work integrates existing capabilities, not a new CUT&RUN engine,
region store, reporter designer, or all-purpose specialist window. Native
TSS-profile inspection is a separate new GUI surface, not part of that first
integration. CUT&RUN is an
optional evidence layer: absent data must not prevent viewing an existing
report or planning a clearly labelled hypothesis. Evidence requirements apply
only when the selected scientific policy explicitly requires them.

## Verified Current State

| Area | Already present | Remaining integration problem |
| --- | --- | --- |
| DNA entry | Feature/map context menus open Promoter design; `Regions...` opens the persistent region manager. | The combined locus inspector is presented under Splicing Expert rather than as a direct DNA evidence view. |
| Locus graphics | Interactive exon/CDS, translation, occupancy/chromatin, regulatory-site, TF-score, saved-region and reporter lanes; hover, zoom, selection and DNA navigation. | Composition is configured through numerous IDs/paths; fresh composition requires an imported isoform panel. |
| Reporter documents | `LocusDocument` loads both locus and reporter-comparison envelopes, preserving reporter rows and exports with live-sequence checks. | Splicing selection changes clear the displayed report; the document, presentation and lifecycle need one sequence-scoped owner independent of that selection. |
| CUT&RUN | Shared catalog/status/prepare/project/read-report operations; GUI support tables and window capture; background BED/BigWig import; an existing TP73 proof-oriented Evidence Preparation dialog. | Generalize the existing dialog into a locus-scoped chooser rather than create a second preparation workflow; do not reuse its synchronous execution for long jobs. |
| Regions | Persistent assembly-bound region sets, import/export, colour, selection/CUT&RUN/Ensembl capture, conservation, and reporter A/B selectors. | Region choices, evidence visibility and planner inputs are spread across views; overlay selection currently feeds Locus figure settings. |
| Reporter planning | Architecture comparison, promoter panels, regulatory-fragment panels, exact product proposals and separate approval-gated creation. | Feature opening replaces all planner UI state and defaults score tracks to the full sequence. Exact-region handoff is not implemented; some panels remain JSON-first, and evidence preparation uses Splicing settings. |
| TSS profiles | Shared strict accession-pinned scoring, context, renderer and receipt-bound SVG/PNG/PDF/GenBank export; reachable through GUI Shell. | No native `TssProfileReport` GUI consumer was found. A first TSS surface and full composite-document parity require separate design, not just exposing an existing widget. |
| Inner agent | Shared biological commands and some GUI destinations already exist. | The shared target catalog lacks first-class locus-evidence and promoter/reporter destinations. This is a discoverability gap, not missing engine reachability. |

Source anchors:

- [DNA entry points](../src/main_area_dna.rs#L6788) and
  [feature-to-Promoter design](../src/main_area_dna/variant_followup.rs#L104).
- [Existing inspector and regressions](../src/main_area_dna/locus_inspector.rs#L88),
  [portable document envelopes](../src/locus_report.rs#L25), and
  [sequence-bound import](../src/main_area_dna/locus_inspector.rs#L628).
- [Isoform-panel prerequisite](../src/main_area_dna/auxiliary_workspaces.rs#L7981),
  [composition request](../src/main_area_dna/auxiliary_workspaces.rs#L8081),
  [composition controls](../src/main_area_dna/auxiliary_workspaces.rs#L9486).
- [CUT&RUN support](../src/main_area_dna/cutrun_support.rs#L74),
  [background track import](../src/app.rs#L10275), and
  [persisted track subscriptions](../src/engine.rs#L7358).
- [Region capture](../src/main_area_dna/genomic_regions_ui.rs#L278),
  [region overlay settings](../src/main_area_dna/genomic_regions_ui.rs#L411),
  [existing A/B selectors](../src/main_area_dna/variant_followup.rs#L6008).
- [Reporter preparation](../src/main_area_dna/variant_followup.rs#L1993),
  [JSON-first panel](../src/main_area_dna/variant_followup.rs#L5649), and
  [architecture controls](../src/main_area_dna/variant_followup.rs#L6542).
- [Shared GUI targets](../crates/gentle-shell/src/ui_intent.rs#L123),
  [current synchronous operation helper](../src/main_area_dna.rs#L21055),
  [report composition boundary](integrated_locus_tss_profiles.md).
- [Whole-state feature seeding](../src/main_area_dna/variant_followup.rs#L243),
  [exact saved-region binding](../src/main_area_dna/variant_followup.rs#L2233),
  [region-to-specialist pattern](../src/main_area_dna/genomic_regions_ui.rs#L106).
- [Report UI wrapper](../src/main_area_dna/auxiliary_workspaces.rs#L9063),
  [Splicing-driven invalidation](../src/main_area_dna/auxiliary_workspaces.rs#L9955),
  [existing preparation dialog](../src/app/evidence_preparation_ui.rs#L227),
  [its synchronous runner](../src/app/evidence_preparation_ui.rs#L682).
- [Manual digest gate](../src/main_area_dna/variant_followup.rs#L2169),
  [strict TSS panel adapter](../src/tfbs_track_panel.rs#L1), and
  [shared UI-intent rules](decisions.md#dec-011-ui-intent-shared-catalog).

This was a source review, not a new live-GUI usability or latency acceptance.
The inspected reporter paths concern designed constructs. A dedicated workflow
for plate-reader measurements, luciferase normalization and measured assay
results was not identified and is not part of this integration proposal.

## Original Codex Plan Sent To Claude

1. Expose/reuse Locus figure from the DNA viewer, with one sequence-scoped owner
   and explicit reference, selected transcript/TSS, region and report identity.
2. Add a guided evidence chooser over existing local tracks, CUT&RUN catalog
   lifecycle and saved read reports; make acquisition explicit.
3. Connect existing region IDs/digests to planner choices, keeping contiguous
   genomic selections distinct from spliced or fusion architectures.
4. Replace JSON-first common paths with guided forms over existing planners;
   retain advanced JSON and separate exact-product approval.
5. Add native read-only TSS detail navigation; decide the full composite report
   contract separately instead of duplicating the Python compositor in GUI code.
6. Include inner-agent/UI-intent parity, background execution, cancellation,
   stale-result rejection and semantic acceptance tests in every increment.

## Claude Feedback And Codex Assessment

- Confirmed: native TSS inspection is new GUI work. There is no source-grounded
  estimate that it exceeds all other work combined, but it must be separated
  from the connective milestones and given its own acceptance contract.
- Confirmed: feature opening replaces `variant_followup_ui`, including score
  bounds; arbitrary saved regions cannot use the feature-only entry unchanged.
  The revised handoff below distinguishes exact-region inputs from newly
  derived transcript-promoter windows instead of promising they are equivalent.
- Confirmed: Splicing refresh clears the report. Move the document, associated
  caches and invalidation together; merely adding a DNA button is insufficient.
- Qualified: a report wrapper and interactive painter are already separate
  functions. The required refactor is to remove their Splicing context/state
  dependency and expose them without the composition form, not build another
  renderer. Splicing-specific assay continuation can remain in its adapter.
- Confirmed: the shared target catalog lacks these destinations. The original
  draft did mention shared UI intents and approval separation, but buried them
  too late. Make intent contracts the first dependency and approval invariants
  explicit; do not describe these existing project rules as newly invented.
- Confirmed: generalize the existing Evidence Preparation dialog and retain its
  TP73 proof as an explicit preset. Its current synchronous executor is not a
  suitable model for responsive long-running preparation.
- Confirmed: guided forms must never fill the approval digest. Architecture
  kinds and geometry already belong to engine types; the GUI passes those
  records through rather than recreating their biological rules.

The following is Codex's revised proposal after that review, not a claim of
Claude approval. The initial implementation boundary is milestones 0 and 1
only; exact-region handoff, evidence preparation and guided design follow in
separate increments.

## Revised Milestones

### 0. Define Shared Navigation Before GUI Wiring

- Follow DEC-011/012 and AGENTS.md: define locus-inspector and Promoter design
  targets in the shared `ui ...` catalog/parser, with structured subject
  arguments and readiness. Menus, inner agent and MCP use the same resolution;
  native handlers only open/focus/close or apply a resolved UI context.
- Design sequence/report identity and feature/region arguments together, but
  advertise executable routes only as their handlers land. The first slice
  delivers the locus route; feature and exact-region planner routes accompany
  milestone 2. Do not add a general Splicing destination just for this task.
- Reuse existing report/reference and saved-region records. Validate sequence,
  exact anchor, report digest and, when supplied, region-set/region digests.
  Feature indices require matching live sequence and annotation identity;
  ambiguous mappings yield a choice or typed unavailability, never a guessed
  nearest transcript.
- Headless invocation returns an inspectable UI intent, not a false claim that
  a window opened. Navigation cannot prepare evidence, save regions, run a
  comparison or approve materialization. Keep those separate operations.

Acceptance: catalog, shared-shell parsing/resolution, MCP/agent discovery and
GUI handlers agree on the delivered destination, subject and readiness. Invalid
or stale subjects cause no writes; UI close/focus does not reseed a planner.

### 1. Separate Inspection From Splicing Composition

- Add the DNA action `Regulatory evidence...` over the existing inspector and
  `LocusDocument` loader. Remove the need to open Splicing Expert or traverse
  its composition inputs just to inspect an existing document. Do not copy
  `render_splicing_locus_evidence_tab` or its large form into a new window.
- Give the document, full reporter envelope, presentation, selection, preview
  and pending load one sequence-scoped owner, independent of Splicing target
  refresh. Keep Splicing's composition form and prerequisites; its inspector
  entry opens/focuses the same owner. A newly composed report replaces the
  inspected document only through an explicit action, not selection sync.
- Reuse the existing painter and portable export paths. Detach the common
  report wrapper from mandatory `SplicingExpertView`; retain Splicing-only
  qPCR actions in their existing adapter when their context is required.
- Loading a compatible stored report needs no isoform panel or CUT&RUN setup.
  Keep exact live-binding checks and historical-report limitations. Changed DNA
  or anchors invalidate live actions; changing a Splicing selection alone does
  not erase the independently inspected report or its reporter overlays.
- Read/parse/preview work runs outside the paint callback, with cancellation,
  progress and stale-completion rejection. Use one semantic owner so opening
  via DNA and Splicing cannot duplicate test identifiers.
- Bound the first slice to inspector close/reopen while the DNA viewer remains
  alive. Saved regions retain their existing project persistence; imported
  documents may be explicitly reloaded after the DNA viewer or project closes.
  Durable document-reference restoration is later work,
  not a new project-store requirement hidden in this navigation change.

Acceptance: a synthetic bound reporter document opens directly from DNA or a
shared intent, retains its rows/exports/hover, stages and explicitly saves a
region, and survives Splicing retargeting and inspector close/reopen. Test
wrong sequence/anchor, historical input, DNA mutation, duplicate entry points
and obsolete load completion. This slice does not yet promise planner handoff,
fresh evidence composition, native TSS pages or automatic project restoration.

### 2. Make Planner Handoff Explicit And Exact

Decision: preserve an exact saved region only in a destination that already
accepts it. Do not substitute a feature-derived promoter window under the same
label. Use the existing conservation/promoter-similarity region handoffs as
the identity, digest and changed-context template.

- Provide distinct actions: open a selected annotated feature for
  transcript-derived promoter analysis, or use a saved region as the candidate
  in the existing regulatory-fragment panel. The former explicitly derives
  new promoter bounds; it does not claim to transfer the original ROI.
- For the exact-ROI route, initially require an explicitly resolved supported
  feature as the Promoter design window context. Resolve by report/provenance
  identity or ask the user to choose; no automatic nearest-feature assignment.
  The feature is navigation context, not evidence that the ROI regulates it.
- Use a narrow `seed -> apply validated destination inputs` sequence on an
  explicit new/retarget action. For the regulatory-fragment destination, apply
  the saved candidate region binding after seeding. Its existing request
  builder already consumes exact regions; do not make every downstream planner
  accept an unowned span or silently retarget all score-track requests.
- Focus/reopen of the same context preserves the draft. Explicit retargeting
  must warn before replacing edited inputs, clear stale proposals/approvals,
  and show the exact candidate coordinates, reference and identity afterwards.
  Missing/ambiguous feature context keeps the ROI inspectable/saveable and
  explains the unavailable planner route. Feature-free planner opening is a
  separate extension, not an inferred biological association.
- Keep report identity for traceability. Merely carrying a source report does
  not attach all its lanes to the planner; evidence inputs remain explicit.

Acceptance: plus/minus exact ROI bindings survive initial seed and focus/reopen;
shared intent and native action yield the same request. Test multiple/no
feature candidates, changed region digest, edited-draft retarget and missing
reference release. A derived-promoter action visibly reports different bounds.

### 3. Generalize The Existing Evidence Preparation Dialog

- Reuse `src/app/evidence_preparation_ui.rs` as the preparation home. Separate
  its TP73 proof defaults into an explicit example preset; ordinary opening
  binds the requested sequence/reference and never injects TP73 paths or data.
  Keep its existing proof/tutorial behavior available without rebuilding
  unrelated array/repeat workflows in this increment.
- Offer prepared CUT&RUN sources, saved read reports, projected/subscribed
  tracks and explicitly chosen BED/BigWig files through existing shared routes.
  List/status is read-only; prepare/download and project/import are separate
  explicit actions. Retain optional paths/JSON as advanced configuration.
- Display exact reference compatibility, source identity/hash, assay/factor or
  chromatin mark, available sample/control/replicate metadata and scale/units.
  Unknown metadata stays unknown; controls are not inferred from a similar name.
- Reuse typed availability states and existing calibration policies. H3K4me3
  remains chromatin context, not factor occupancy; different sources do not gain
  quantitative comparability by sharing a colour, label or scale.
- Keep display visibility separate from evidence included in a planning request.
  A hidden lane is not deleted evidence. A requested but unavailable lane remains
  visible in readiness/summary, not represented as zero.
- Provide nearby setup/import help at the point of use. Do not require a full
  raw-read alignment/acquisition workflow for the first processed-evidence GUI.
- Replace long synchronous calls on these new paths with the existing task and
  import patterns; generalizing a form must not reproduce its blocking runner.
  Return explicit evidence selection to the composition request, not global
  Splicing fields. No new biological operation is expected for the chooser.

Acceptance: local signal plus an unavailable control gives the same typed
report as shell execution; incompatible assembly cannot be attached silently;
cancelled preparation leaves no falsely ready source or stale report adoption.
Opening a second non-TP73 locus cannot inherit the TP73 preset. Keep the current
public proof runnable. Persistent report-reference restoration, if included,
must use inspectable project metadata with reload/relocate on changed files;
never persist preview caches, workers or approval tokens.

### 4. Guided Reporter Planning, Not A New Designer

- Reuse the existing saved-region selectors and exact transcript/TSS choices.
  Overlapping TSS windows remain separately selectable; shared exon structure
  or occupancy must not silently choose the biologically active promoter.
- Keep two entry types: contiguous source region and transcript-aware
  architecture. A spliced leader or fusion must retain ordered source segments,
  transcript identity and explicit ATG/vector policy; never flatten it to a ROI.
  Reuse `PromoterReporterArchitectureKind` and the engine's geometry records;
  these are existing type/validation guarantees, not new GUI calculations.
- Add guided controls for common architecture and promoter-panel requests,
  exact vector identity, controls, boundaries and effective mutation policy.
  Route comparison, readiness, planning and materialization through the existing
  engine operations. Do not add GUI-side coordinate or mutation calculations.
- Audit the architecture-to-product handoff before promising a one-click route.
  If a chosen architecture lacks a supported product-planner input, expose that
  limit; add a bounded shared conversion only after its scientific contract is
  agreed, rather than synthesizing geometry in an event handler.
- Preserve immutable proposals and separate digest-bound creation approval.
  Input/evidence drift invalidates readiness and approval. Imported reports are
  inspectable evidence, not permission to recreate or modify their experiments.
  Guided controls may fill requests only: never prefill, copy into, restore or
  auto-submit the human approval field from `proposal.proposal_digest`. Retain
  the current explicit exact-match gate and add a regression for it.
- Use the existing gene-set study composer for future multi-gene entry; no
  hidden GUI loops with different per-gene policies.

Acceptance: one native-sequence hypothesis reaches an exact product proposal
through GUI and shell with matching scientific fields. Test stale inputs,
missing vector/control requirements, noncontiguous architecture handling,
cancelled work and explicit approval before creation. No wet-lab validation claim.

## Separate Follow-Up: First Native TSS Profile Surface

This is a new GUI consumer, outside milestones 0-4 and requiring its own plan
and approval. Existing locus score lanes do not constitute a native
`TssProfileReport` inspector. No delivery-size comparison is established yet.

- Add a native view of stored `TssProfileReport` windows, selecting by promoter
  ID and checked reference/gene/strand/TSS identity. Different-size overview,
  evidence and detail windows must not be joined by equal sequence hashes.
- Reuse the checked geometry and stored scores, full hover provenance, printed
  coordinate summaries, distinct TSS bands and shared sequence exports.
  Redraw, zoom, selection and export never trigger hidden rescoring.
- Preserve the fail-closed `TfbsTrackPanel` adapter contract for score grammar,
  clipping, strand and calibration policies. A new GUI cannot silently convert
  heterogeneous tracks into a panel the shared TSS renderer would reject.
- Separate native imported-report inspection from the larger composite-document
  contract. The existing Python compositor is not already an engine-native GUI
  document; do not launch it silently or reimplement its receipt checks in a
  frontend. Decide that promotion as a bounded follow-up if needed.

Acceptance: two overlapping TSSs, including a negative-strand case, retain
distinct identities and exact motif/occupancy positions in native inspection
and shared SVG/export. GenBank bases and feature locations match the bound
report. A PDF/SVG alone cannot reconstruct the missing scientific report.

## Cross-Cutting Requirements

- Inner agent and direct GUI actions resolve the same subject and readiness,
  using shared catalog/operation paths. New destination names are proposals,
  not commands that can already be executed. No recursive agent calls or
  automatic acquisition, materialization or screenshot sharing.
- Reuse cancellation/progress machinery. New expensive paths must not use the
  synchronous engine-write-lock helper from a paint callback. Do not solve this
  by copying an entire large project per frame. Cache immutable presentation,
  virtualize dense lanes and reject results bound to changed inputs.
- Preserve source annotations, motif predictions, occupancy, chromatin context,
  response/expression and measured reporter outcomes as separate evidence.
  "Confirmed" legacy occupancy status must be explained as overlap support,
  not validated TF binding or promoter activity.
- Test shared operations and GUI semantic navigation on public synthetic inputs
  first. Reuse existing fixtures with provenance. Glen independently checks
  live GUI responsiveness, plus/minus real-data reports and exact-revision
  acceptance; screenshot appearance alone is not a scientific gate.

## Decisions For The Discussion

1. Recommended home: DNA viewer for context and selection; existing Promoter
   design for hypotheses and construction. No new top-level reporter specialist.
2. Recommended first approval: milestones 0 and 1 only, delivering direct
   report inspection and stable in-session ownership through GUI/shared intent.
   No planner handoff or persistent document store is smuggled into that slice.
3. Recommended handoff policy: exact region into an existing region-aware
   planner after explicit feature-context resolution and seeding; separately
   labelled derivation for transcript-promoter analysis. Do not lose selected
   bases for the sake of a cheaper navigation shortcut.
4. Generalize the existing Evidence Preparation dialog. Keep first native TSS
   inspection, durable document restoration and measured assay analysis
   separately bounded. None becomes an extra `.10` release gate by this draft.
5. The supplied Claude review is incorporated, with the qualifications above.
   A further read-only review can focus on this narrower first slice and the
   exact-region handoff policy before implementation is authorized.
