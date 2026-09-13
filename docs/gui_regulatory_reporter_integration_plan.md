# DNA-Centred Regulatory Evidence And Reporter Planning

Status: discussion draft, not approved for implementation.
Source review: `d978b0c4368071f839ceecf8d494212ad9b992ec`, 2026-09-13.
Claude review: **pending**. The user reports a separate Claude review is in
progress. This session's non-interactive read-only CLI attempt returned
`401: OAuth access token expired`; no feedback or approval has been incorporated.
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

This is integration of existing capabilities, not a new CUT&RUN engine, region
store, reporter designer, or all-purpose specialist window. CUT&RUN is an
optional evidence layer: absent data must not prevent viewing an existing
report or planning a clearly labelled hypothesis. Evidence requirements apply
only when the selected scientific policy explicitly requires them.

## Verified Current State

| Area | Already present | Remaining integration problem |
| --- | --- | --- |
| DNA entry | Feature/map context menus open Promoter design; `Regions...` opens the persistent region manager. | The combined locus inspector is presented under Splicing Expert rather than as a direct DNA evidence view. |
| Locus graphics | Interactive exon/CDS, translation, occupancy/chromatin, regulatory-site, TF-score, saved-region and reporter lanes; hover, zoom, selection and DNA navigation. | Composition is configured through numerous IDs/paths; fresh composition requires an imported isoform panel. |
| Reporter documents | `LocusDocument` loads both locus and reporter-comparison envelopes, preserving reporter rows and exports with live-sequence checks. | Loaded reports live in window caches; viewing a report and configuring a new comparison are not one seamless handoff. |
| CUT&RUN | Shared catalog/status/prepare/project/read-report operations; GUI support tables and window capture; background BED/BigWig import. | No integrated locus-scoped chooser for prepared datasets, projected tracks, controls and read reports. |
| Regions | Persistent assembly-bound region sets, import/export, colour, selection/CUT&RUN/Ensembl capture, conservation, and reporter A/B selectors. | Region choices, evidence visibility and planner inputs are spread across views; overlay selection currently feeds Locus figure settings. |
| Reporter planning | Architecture comparison, promoter panels, regulatory-fragment panels, exact product proposals and separate approval-gated creation. | Some panels remain JSON-first; architecture preparation depends on Splicing settings to attach quantitative locus evidence. |
| TSS profiles | Shared strict accession-pinned scoring, context, renderer and receipt-bound SVG/PNG/PDF/GenBank export; reachable through GUI Shell. | The native locus loader does not accept `TssProfileReport`; native TSS detail navigation and full composite-document parity are separate work. |
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

No Claude critique has yet been incorporated. The following ordering is Codex's provisional
minimal implementation sequence, not a Claude-approved revision.

## Proposed Milestones

### 1. Connect The Existing Views

- Add one DNA-viewer action, provisionally `Regulatory evidence...`, that
  exposes the existing locus inspector without first navigating Splicing tabs.
  Keep the old Splicing entry as an alias, not a second state owner or painter.
- Start with loading an existing exact-compatible locus/reporter document.
  This read-only path needs no newly configured isoform panel or CUT&RUN data.
  Fresh composition still uses its existing prerequisites and explicit setup.
- From the inspector, stage/save a region using the existing manager, then
  explicitly open the appropriate Promoter design section with that selection.
  Do not silently save on drag or reset unrelated planner choices on focus.
- Resolve the small handoff state explicitly: sequence/reference binding,
  source report identity/digest, selected region ID/set digest, and optional
  transcript/TSS/architecture ID. Reuse current records; do not invent a new
  biological study schema simply to move between views.
- Define reopening: persist bounded inputs and report references/bindings using
  project conventions, never PNG caches, worker handles or approval tokens.
  Missing/changed referenced files show a reload/relocate state, not an automatic
  scientific recomputation. Any stored-reference extension must remain
  inspectable through the shared project/agent contract.
- Add one shared UI destination and semantic control IDs for this first path.
  Read/parse/render work must not block the frame; results return only to the
  still-matching sequence/report selection. Reuse existing worker patterns.

Acceptance: an existing synthetic bound reporter report opens from DNA, retains
its reporter rows, stages and saves a region, then opens the existing planner
on that exact region. No new design algorithm, download or construct creation.
Test wrong sequence/anchor, historical unbound report, duplicate entry points,
reopening, and obsolete background completion. One UI owner must also prevent
duplicate semantic registrations when both entry paths are used.

### 2. Attach Evidence Without Knowing IDs Or JSON

- Offer prepared CUT&RUN sources, saved read reports, projected/subscribed tracks
  and explicitly chosen local BED/BigWig files through existing shared routes.
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

Acceptance: local signal plus an unavailable control gives the same typed
report as shell execution; incompatible assembly cannot be attached silently;
cancelled preparation leaves no falsely ready source or stale report adoption.

### 3. Guided Reporter Planning, Not A New Designer

- Reuse the existing saved-region selectors and exact transcript/TSS choices.
  Overlapping TSS windows remain separately selectable; shared exon structure
  or occupancy must not silently choose the biologically active promoter.
- Keep two entry types: contiguous source region and transcript-aware
  architecture. A spliced leader or fusion must retain ordered source segments,
  transcript identity and explicit ATG/vector policy; never flatten it to a ROI.
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
- Use the existing gene-set study composer for future multi-gene entry; no
  hidden GUI loops with different per-gene policies.

Acceptance: one native-sequence hypothesis reaches an exact product proposal
through GUI and shell with matching scientific fields. Test stale inputs,
missing vector/control requirements, noncontiguous architecture handling,
cancelled work and explicit approval before creation. No wet-lab validation claim.

### 4. TSS Detail And Export Parity

- Add a native view of stored `TssProfileReport` windows, selecting by promoter
  ID and checked reference/gene/strand/TSS identity. Different-size overview,
  evidence and detail windows must not be joined by equal sequence hashes.
- Reuse the checked geometry and stored scores, full hover provenance, printed
  coordinate summaries, distinct TSS bands and shared sequence exports.
  Redraw, zoom, selection and export never trigger hidden rescoring.
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
2. Recommended first approval: milestone 1 only. Evidence-source ergonomics
   follow immediately, but are not a prerequisite for testing the connection.
3. Keep native TSS supplements and assay measurement/normalization separate from
   that first slice. Do not quietly enlarge the `.10` release gate.
4. Ask Claude to challenge the first slice's size, persistence boundary,
   Splicing prerequisite removal from navigation, and discontinuous-architecture
   handoff before implementation. Incorporate the pending review before proceeding.
