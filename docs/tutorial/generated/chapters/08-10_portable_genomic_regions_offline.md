---
chapter_id: "portable_genomic_regions_offline"
title: "Save and Share Genomic Regions (Offline SERPINE1 Example)"
tier: "core"
example_id: "portable_genomic_regions_offline"
source_example: "docs/examples/workflows/portable_genomic_regions_offline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-09-04"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/portable_genomic_regions_offline"
---

# Save and Share Genomic Regions (Offline SERPINE1 Example)

Preserve a selected genomic span as a portable, assembly-bound region; capture source-specific Ensembl and CUT&RUN evidence without flattening its provenance; export lossless JSON and BED6 plus a binding manifest; and reuse the same region set on a locus figure.

A coordinate copied from a figure is easy to misread. BED uses zero-based, half-open coordinates, while a colleague will usually expect a one-based, inclusive interval. A chromosome name alone also does not identify an assembly. GENtle therefore stores one canonical genomic interval with its species, assembly, contig identity, strand, and explicit coordinate convention. A sequence-local projection is optional and remains bound to the exact sequence digest and genome anchor that produced it. Human text and BED are projections of that record, not independent sources of truth.

The reason for selecting a region is equally important, but it is not the geometry. In this example one region comes from a manual selection, one from a synthetic Ensembl Regulation row, and one from a synthetic CUT&RUN support window. Each retains a typed selection method and source evidence. The Ensembl association remains a provider annotation; CUT&RUN enrichment remains occupancy evidence rather than affinity or causal activity. JSON is the lossless exchange format. BED6 is accompanied by a content-bound manifest because BED alone cannot preserve assembly accessions, evidence, derivation, or local projection identity.

**Prerequisites:** Read [Chapter 29: Promoter-Reporter Panel Planning (Offline Approval-Gated Demo)](./08-09_promoter_reporter_panel_planning_offline.md), [Chapter 27: PATZ1 Gene-Locus Evidence Composition (Offline Synthetic Demo)](./06-06_patz1_gene_locus_evidence_offline.md) first.

## Parameters That Matter

- `start_0based / end_0based_exclusive` (where used: gentle.genomic_region_of_interest.v1 and BED6)
  - Why it matters: This is the canonical machine geometry. The start base is included and the end base is excluded, so length is always end minus start.
  - How to derive it: For a human interval A-B (1-based inclusive), use BED start A-1 and BED end B. Do not subtract one from the end.
- `assembly_name=GRCh38 / assembly_accession=GCA_000001405.15` (where used: every saved interval and the BED manifest)
  - Why it matters: Chromosome 7 is not a portable coordinate without a reference assembly; accessions make similarly named builds distinguishable.
  - How to derive it: Use the exact verified genome anchor or an explicit caller-supplied reference. Never infer it from a filename or gene symbol.
- `selection_method` (where used: manual, Ensembl-regulatory, and CUT&RUN captures)
  - Why it matters: It records how the interval entered the set without promoting that source to evidence of affinity or causal regulation.
  - How to derive it: Choose the typed method matching the source operation; derived unions/intersections/hulls use a separate explicit operation.
- `manifest_path` (where used: BED export and import)
  - Why it matters: BED6 cannot carry complete assembly, evidence, projection, and digest semantics. The manifest binds those records to the exact BED bytes.
  - How to derive it: Keep the generated `.bed.manifest.json` beside the BED and supply both paths when importing.

## When This Routine Is Useful

- You want to send a promoter or enhancer candidate to a colleague without an off-by-one ambiguity.
- You want to retain why a region was selected without turning an annotation or occupancy signal into a verdict.
- You need BED for another tool but also need a lossless round trip back into GENtle.
- You want selected regions to remain reusable after the original evidence file is unavailable.
- You want to compare saved candidates on the same strand-aware locus axis as transcripts, CDS, occupancy, motif scores, and reporter architectures.

## What You Learn

- Distinguish canonical genomic coordinates, human display coordinates, BED coordinates, and sequence-local coordinates.
- Keep genomic strand distinct from a locus viewer's local/display orientation.
- Explain why region geometry, selection purpose, and evidence provenance are separate fields.
- Capture manual, provider-annotation, and occupancy-derived regions through the same engine-owned contract.
- Use canonical JSON for lossless sharing and BED6 plus its manifest for external-tool interoperability.
- Recognize assembly mismatch as an error that requires an explicit liftover workflow rather than an automatic coordinate change.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.

## At a Glance

1. Open the pinned SERPINE1 tutorial locus and select a short genomic span in th...
2. In Saved genomic regions, choose a set ID, label the selection, and save it. ...
3. Use Copy human coordinates (1-based inclusive), then Copy BED row (0-based ha...
4. Copy the canonical ROI JSON and confirm it contains zero_based_half_open, the...
5. Capture the synthetic Ensembl regulatory row from the Locus figure table and ...
6. Export JSON... for lossless exchange, then BED + manifest... for a BED6 consu...
7. Select serpine1_regions in the Locus figure Saved region set IDs field and co...
8. Import the BED with its manifest into a fresh project and inspect the region-...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open the pinned SERPINE1 tutorial locus and select a short genomic span in th...

GUI: Open the pinned SERPINE1 tutorial locus and select a short genomic span in the DNA view. Open its context menu and choose `Save/share selected genomic region...`.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/portable_genomic_regions_offline.json
```

> Expected: The selected local span becomes one canonical GRCh38 genomic interval; local and genomic coordinates are both retained and explicitly labelled.

### Step 2: In Saved genomic regions, choose a set ID, label the selection, and save it. ...

GUI: In `Saved genomic regions`, choose a set ID, label the selection, and save it. Confirm the table shows assembly, 1-based inclusive coordinates, strand, `manual_span`, and evidence availability.

CLI:

```bash
jq '.regions[] | {region_id,label,interval,selection_method,evidence,identity_sha256,content_sha256}' artifacts/portable_genomic_regions.region_set.json
```

> Expected: Saving the manual span uses the shared `CaptureGenomicRegion` operation rather than a GUI-only feature annotation.

### Step 3: Use Copy human coordinates (1-based inclusive), then Copy BED row (0-based ha...

GUI: Use `Copy human coordinates (1-based inclusive)`, then `Copy BED row (0-based half-open)`. Compare the boundaries: BED start is one less; the inclusive human end equals the BED-exclusive end.

CLI:

```bash
cat artifacts/portable_genomic_regions.bed
```

> Expected: The human representation is 1-based inclusive and the BED representation is 0-based half-open; a one-base region remains one base in both forms.

### Step 4: Copy the canonical ROI JSON and confirm it contains zero_based_half_open, the...

GUI: Copy the canonical ROI JSON and confirm it contains `zero_based_half_open`, the GRCh38 assembly/accession, contig identity, local projection, selection method, and both identity and content digests.

CLI:

```bash
jq '{bed_sha256,coordinate_convention,columns,region_set_sha256:.region_set.content_sha256}' artifacts/portable_genomic_regions.bed.manifest.json
```

> Expected: Canonical JSON carries the complete portable record. Label and notes can change the content digest without changing the immutable identity digest.

### Step 5: Capture the synthetic Ensembl regulatory row from the Locus figure table and ...

GUI: Capture the synthetic Ensembl regulatory row from the Locus figure table and the synthetic CUT&RUN support window from its support-window table. Confirm the three rows retain different methods and evidence statements.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-region-import.json shell 'regions import {"path":"artifacts/portable_genomic_regions.bed","format":"bed","manifest_path":"artifacts/portable_genomic_regions.bed.manifest.json"}'
```

> Expected: The Ensembl and CUT&RUN captures retain source IDs and conservative non-claims rather than becoming generic prose or causal assertions.

### Step 6: Export JSON... for lossless exchange, then BED + manifest... for a BED6 consu...

GUI: Export `JSON...` for lossless exchange, then `BED + manifest...` for a BED6 consumer. Keep the manifest beside the BED when the region set may return to GENtle.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-region-import.json shell 'regions list'
```

> Expected: The BED manifest binds the exact BED bytes, BED6 column contract, complete region set, assembly identity, evidence, and region-set digest.

### Step 7: Select serpine1_regions in the Locus figure Saved region set IDs field and co...

GUI: Select `serpine1_regions` in the Locus figure `Saved region set IDs` field and compose the view. Confirm the three compact saved-region overlays align with the transcript/CDS, occupancy, TP73 score, and reporter rows.

CLI:

```bash
cargo run --bin gentle_cli -- --state /tmp/gentle-region-import.json shell 'regions inspect {"set_id":"serpine1_regions","region_id":"serpine1_manual_candidate"}'
```

> Expected: The saved-region overlay uses the locus report's explicit local-axis mapping and appears once on the plus-strand SERPINE1 locus.

![A pinned public Ensembl-116 SERPINE1 locus with three saved regions: a manual candidate, a synthetic Ensembl-regulatory annotation, and a synthetic CUT&RUN support window. Their shared position is for inspection and does not establish binding or causal regulation.](../artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.locus.svg)

*Figure: A pinned public Ensembl-116 SERPINE1 locus with three saved regions: a manual candidate, a synthetic Ensembl-regulatory annotation, and a synthetic CUT&RUN support window. Their shared position is for inspection and does not establish binding or causal regulation. Regenerate with `cargo run --bin gentle_examples_docs -- tutorial-generate`.*

> SVG text labels: `SERPINE1 promoter-reporter architectures | 5' -> 3' display, genomic 101121158 -> 101139263 | strand + | assembly GRCh38 | annotation Ensembl 116 pinned public fixture | upstrea...`. If the embedded preview omits text in the GUI, open the linked SVG or use these labels as the figure legend.

### Step 8: Import the BED with its manifest into a fresh project and inspect the region-...

GUI: Import the BED with its manifest into a fresh project and inspect the region-set digest. Try the BED alone only with an explicit assembly/reference request; GENtle does not infer GRCh38 from chromosome 7.

> Expected: A manifest-bound import round-trips losslessly; bare BED requires an explicit genome reference and does not trigger liftover or network access.


## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/portable_genomic_regions_offline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/portable_genomic_regions_offline.json'
```

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/portable_genomic_regions_offline.json
cargo run --bin gentle_examples_docs -- tutorial-generate
cargo run --bin gentle_examples_docs -- tutorial-check
```

## Checkpoints

- The workflow executes offline against the pinned public SERPINE1 sequence and annotation; no provider is contacted.
- The region-set JSON contains exactly three regions with `manual_span`, `ensembl_regulatory_feature`, and `cutrun_support_window` methods.
- Every interval is bound to Homo sapiens, GRCh38/GCA_000001405.15, chromosome 7/NC_000007.14, and `zero_based_half_open` coordinates.
- The BED has exactly three BED6 rows and its SHA-256 matches the manifest.
- The Ensembl evidence retains the feature ID and canonical provider URL; CUT&RUN evidence states occupancy/enrichment but not affinity or causality.
- The locus SVG contains stable saved-region markers for all three methods alongside transcript/CDS, occupancy, TP73 score, and reporter-architecture rows.
- The retained receipt and report omit workstation-local paths by default while retaining source hashes.

## What This Chapter Produces

- [`artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.locus.svg`](../artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.locus.svg)

  - Embedded above near Step 7; kept here as an audit link.

> SVG text labels: `SERPINE1 promoter-reporter architectures | 5' -> 3' display, genomic 101121158 -> 101139263 | strand + | assembly GRCh38 | annotation Ensembl 116 pinned public fixture | upstrea...`. If this embedded preview omits text in the GUI, open the linked SVG or use these labels as the figure legend.

- [`artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.bed.manifest.json`](../artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.bed.manifest.json) - schema: `gentle.genomic_region_bed_manifest.v1`
- [`artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.region_set.json`](../artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.region_set.json) - schema: `gentle.genomic_region_set.v1`
- [`artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.bed`](../artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.bed) - `7 101126000 101126120 Synthetic_TP73_CUT_RUN_support 0 .`

## Tutorial Provenance

- Chapter id: `portable_genomic_regions_offline`
- Tier: `core`
- Example id: `portable_genomic_regions_offline`
- Tutorial source JSON: `docs/tutorial/sources/08-10_portable_genomic_regions_offline.json`
- Workflow file: `docs/examples/workflows/portable_genomic_regions_offline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/portable_genomic_regions_offline`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-09-04`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Save and Share Genomic Regions (Offline SERPINE1 Example)`
- Tutorial/chapter id: `portable_genomic_regions_offline`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
