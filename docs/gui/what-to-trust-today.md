## What To Trust Today

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

If you are opening GENtle as a biologist rather than as a contributor, use
this as a confidence map for the current GUI surface.

### Recommended now

- Single-insert Gibson specialist work:
  `Patterns -> Gibson...` for destination-first planning, preview, apply,
  reopen-from-lineage, and SVG/cartoon export.
- Core PCR / primer-pair / qPCR routes when you already know the intended task
  and want deterministic execution plus inspectable outputs.
- Prepared-genome retrieval/extraction flows once the reference has already
  been prepared locally.
- Visual explanation/export paths such as lineage SVG, protocol cartoons,
  dotplots, and isoform architecture.
- Bench-facing cloning handoffs through `File -> Export Lab Assistant
  Report...`, with ODT as the open editable default plus DOCX/Markdown export
  choices from the same engine report.
- External-service quote/handoff inspection through `Services -> External
  Services...`, using the same provider catalog, preflight, and quote payload
  contracts as CLI/agent/ClawBio routes. Metabion and GeneArt appear only via
  the shared provider catalog; the GUI does not duplicate vendor logic.
  Discovery-only suppliers, including Promega, remain `catalog_only` and cannot
  earn quote/submission eligibility. The shared catalog JSON exposes country
  and dated HTML/text-PDF protocol references for filtering; no new graphical
  country filter or protocol executor is implied. [Tutorial 09.04](../tutorial/09-04_cloning_time_cost_and_protocol_suppliers.md)
  explains these boundaries alongside advisory cloning time and cost.
- TP73 evidence-viewer material preparation through the command palette action
  `Evidence Preparation`, which guides local repeat, array, BED, TFBS, and
  proof-export steps while showing copyable CLI/R handoff commands for
  external CEL/probeset preparation.

### Works with caveats

- Multi-insert Gibson preview is useful, but execution currently requires a
  defined destination opening; `existing_termini` remains the single-fragment
  handoff path.
- Sequencing confirmation now has a dedicated specialist that accepts called
  reads, direct ABI/AB1/SCF trace import, optional baseline context, and
  both variant-focused chromatogram review and trace-base browsing, but it is
  still not a full chromatogram editor.
- Primer3-backed workflows are available, but the internal backend is still
  the more predictable default while parity hardening continues.
- PrimerBank records are reachable from the GUI Shell through
  `primers primerbank search|show|test-cdna`. The same shared route can parse a
  saved PrimerBank result, export `gentle.primerbank_search.v1`, and test one
  selected pair against project transcript models. A dedicated catalog browser
  is not required for parity; GENtle does not mirror the external database or
  turn catalog presence into a validation claim. Search results expose the
  requested/observed species cross-check, and cDNA testing requires a confirmed
  human or mouse PrimerBank match. The same report separately checks any
  `organism` annotation on the selected project sequence and refuses a known
  mismatch; missing target-species metadata stays visibly unresolved.
- Manual GUI walkthroughs in Help/Tutorial are useful for orientation, but the
  generated executable tutorials remain the stronger reproducibility baseline.
- The command palette and feature/map context menus open **Feature Editor**.
  Its Location tab edits exact simple annotations and exact segments in flat
  `Join`/`Order` compounds. Its Create tab appends one exact forward/reverse
  feature with ordered, duplicate-capable, optionally valueless qualifiers;
  its Delete tab previews and removes a complete feature of any existing
  location shape. Its Split tab divides an exact simple feature at an internal
  base boundary while preserving its record fields on both explicit outputs;
  its Merge tab combines only exactly touching simple records whose kind,
  strand, and ordered qualifiers already match. Every mode requires Preview
  before Apply, refuses stale previews, and uses shared engine operations.
  Curation previews list overlap and recognized shared identifiers as review
  evidence without changing related annotations. Split/Merge reject
  nested/fuzzy locations, gaps, overlaps, and metadata conflicts rather than
  guessing a biological repair.
  The [offline Feature Editor tutorial](../tutorial/02-05_feature_editor_gui_cli.md)
  walks through every tab on a 120 bp synthetic starter, including preview
  invalidation, unchanged DNA and save/reopen checks.

### Exploratory / not yet first choice

- Broader routine-family coverage outside the strongest current Gibson and
  restriction-centered paths.
- Propagated transcript-boundary curation, topology-changing compound edits,
  and nested/fuzzy feature-location editing.
- guideRNA workflows, deeper assay families, and richer virtual-PCR /
  off-target analysis paths.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Detailed Reference](detailed-reference.md) | [Next: Start the GUI](start-the-gui.md)
<!-- manual-page-nav:end -->
