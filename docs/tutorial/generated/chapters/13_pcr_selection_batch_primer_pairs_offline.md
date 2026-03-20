# Selection-first PCR batch primer design (offline)

- Chapter id: `pcr_selection_batch_primer_pairs_offline`
- Tier: `core`
- Example id: `pcr_selection_batch_primer_pairs_offline`
- Source example: `docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Queue multiple ROIs from one template and run deterministic primer-pair reports with optional extracted-region artifacts.

Regular PCR planning often starts from visual selections in the DNA window rather than manual coordinate typing. This chapter focuses that selection-first workflow: choose two regions, optionally materialize extracted-region copies for downstream handoff, and run one primer-design report per region with deterministic naming (`..._r01`, `..._r02`).

## When This Routine Is Useful

- You want to design primer pairs for several selected regions without retyping ROI coordinates.
- You need copied sequence artifacts for each PCR target before forwarding to other workflows.
- You want deterministic report IDs for review and export.

## What You Learn

- Use DNA-window selections as first-class PCR queue inputs.
- Run deterministic multi-region primer design in one batch action.
- Understand optional extraction-copy artifacts as downstream handoff objects.

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

1. Load `test_files/pGEX_3X.fa` and open `Engine Ops -> Primer and qPCR design reports`.
2. Use DNA-window `PCR ROI` actions to add current selection and selected features to the PCR region queue.
3. In `Design primer pairs`, set shared constraints once and run `Design Primer Pairs for queued regions`.
4. Optionally enable `Also create extracted region copies`, then export one generated report via `Export report_id...`.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json'
```

## Parameters That Matter

- `DesignPrimerPairs.report_id base` (where used: batch execution in Design primer pairs form)
  - Why it matters: Batch runs derive deterministic report IDs from this base with `_rNN` suffixes.
  - How to derive it: Pick one short base tied to experiment context (for example `tp73_promoter_batch`).
- `PCR queue copy toggle (`Also create extracted region copies`)` (where used: batch execution in Design primer pairs form)
  - Why it matters: Enables per-region `ExtractRegion` artifacts before/alongside primer report generation.
  - How to derive it: Enable when downstream steps need sequence artifacts; leave off for report-only planning.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json
cargo run --bin gentle_cli -- shell 'primers list-reports'
```

## Checkpoints

- Two primer reports are created with deterministic `_r01`/`_r02` suffixes.
- Optional extracted-region copy artifacts are created when the queue toggle is enabled.
- At least one report can be exported from GUI report helpers without manual JSON editing.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json`
- Inspect this JSON file directly when you need full option-level detail.
