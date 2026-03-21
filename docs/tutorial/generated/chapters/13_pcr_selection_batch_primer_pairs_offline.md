# Selection-first PCR batch primer design (offline)

- Chapter id: `pcr_selection_batch_primer_pairs_offline`
- Tier: `core`
- Example id: `pcr_selection_batch_primer_pairs_offline`
- Source example: `docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Queue TP73-AS2 promoter-proximal ROIs from local TP73 sequence context, then run deterministic batch primer-pair reports with optional extracted-region copies.

This chapter is PCR-focused and biology-anchored without requiring genome mapping mechanics. You work directly on local TP73 locus sequence (`test_files/tp73.ncbi.gb`), queue promoter-proximal regions around TP73-AS2, run one shared primer-constraint set across queued regions, and review/export deterministic batch results (`..._r01`, `..._r02`). This prepares a practical assay baseline before entering the separate luciferase cloning tutorial.

## When This Routine Is Useful

- You want to combine one ad hoc mouse/text selection with one or more annotated features in one PCR batch run.
- You want to amplify TP73-AS2 promoter-proximal windows from local TP73 sequence context.
- You want to design primer pairs for several selected regions without retyping ROI coordinates.
- You want PCR-first preparation that cleanly feeds later promoter-luciferase planning.
- You need copied sequence artifacts for each PCR target before forwarding to other workflows.
- You want deterministic report IDs for review and export.

## What You Learn

- Capture PCR queue regions from both ad hoc selection and selected features.
- Run deterministic multi-region primer design in one batch action with stable `_rNN` report suffixes.
- Use batch-results actions (`Show`, `Export`, `Open`) for validation and downstream handoff.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md).
  - Reoccurs in: [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md).
  - Reoccurs in: [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md).
  - Reoccurs in: [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).

## GUI First

1. Load `test_files/tp73.ncbi.gb` and open `Engine Ops -> Primer and qPCR design reports`.
2. Locate `TP73-AS2` annotation and make one promoter-proximal ROI selection near the antisense transcript 5' boundary (for example around genomic coordinates `61720..62120`).
3. Create one map/text ROI selection, then use DNA-window `PCR ROI -> Add current selection to queue`.
4. Select the `TP73-AS2` feature (or another nearby regulatory feature) and use DNA-window `PCR ROI -> Add selected feature(s) to queue`.
5. In `Design primer pairs`, set shared constraints and a `report_id` base (for example `tp73_as2_promoter_batch`).
6. Optionally enable `Also create extracted region copies`, then run `Design Primer Pairs for queued regions`.
7. In `PCR batch results`, use `Show` for report summary, `Export` for JSON output, and `Open` to jump to the extracted copy (or template when copy mode is off).

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json'
```

## Parameters That Matter

- `TP73-AS2 promoter-proximal ROI window (for example `61720..62120`)` (where used: Design primer pairs ROI fields and queue sourcing around TP73-AS2)
  - Why it matters: Keeps PCR targets tied to a biologically meaningful promoter-proximal context that can later feed luciferase planning.
  - How to derive it: Anchor on TP73-AS2 annotation in `test_files/tp73.ncbi.gb` (gene/ncRNA at lines ~1257-1270), then choose a local promoter-proximal window around the antisense 5' boundary.
- `PCR ROI source action (`Add current selection` vs `Add selected feature(s)`)` (where used: DNA-window `PCR ROI` menu before queue execution)
  - Why it matters: Defines whether region bounds come from manual selection or feature coordinates and keeps queue-source labels interpretable.
  - How to derive it: Use current selection for ad hoc boundaries; use selected features for annotation-driven boundaries.
- `DesignPrimerPairs.report_id base` (where used: batch execution in Design primer pairs form)
  - Why it matters: Batch runs derive deterministic report IDs from this base with `_rNN` suffixes.
  - How to derive it: Pick one short base tied to experiment context (for example `tp73_as2_promoter_batch`).
- `PCR queue copy toggle (`Also create extracted region copies`)` (where used: batch execution in Design primer pairs form)
  - Why it matters: Enables per-region `ExtractRegion` artifacts before/alongside primer report generation.
  - How to derive it: Enable when downstream steps need sequence artifacts; leave off for report-only planning (default off).

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
cargo run --bin gentle_cli -- shell 'primers list-reports'
cargo run --bin gentle_cli -- shell 'primers show-report tp73_as2_promoter_batch_r01'
cargo run --bin gentle_cli -- shell 'primers export-report tp73_as2_promoter_batch_r01 tp73_as2_promoter_batch_r01.json'
```

## Checkpoints

- PCR queue rows include both current-selection and selected-feature source labels.
- Two primer reports are created with deterministic `_r01`/`_r02` suffixes and appear in batch results.
- If your runtime falls back to the internal primer backend, reports may validly contain zero accepted pairs under strict constraints; report creation/export behavior remains deterministic.
- If copy toggle is enabled, each batch-results row includes a copy ID and `Open` jumps to that extracted sequence.
- At least one report can be shown and exported from the batch-results row actions without manual JSON editing.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json`
- Inspect this JSON file directly when you need full option-level detail.
