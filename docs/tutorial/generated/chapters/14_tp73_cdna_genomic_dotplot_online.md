# Compare TP73 cDNA against TP73 genomic context via dotplot (online)

- Chapter id: `tp73_cdna_genomic_dotplot_online`
- Tier: `online`
- Example id: `tp73_cdna_genomic_dotplot_online`
- Source example: `docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.

Fetch TP73 cDNA, extract TP73 from GRCh38, and compute a pair-forward dotplot with high-sensitivity defaults (`word=7`, `step=1`, `mismatches=0`).

This chapter captures a practical cDNA-vs-genomic verification route for transcript structure interpretation. The focus is a reproducible first-pass map that reveals exon-aligned block patterns while preserving one shared operation path across GUI and CLI interfaces.

## When This Routine Is Useful

- You want a deterministic cDNA-vs-genomic control for exon/intron-aware interpretation.
- You need a reproducible TP73 baseline for demonstrating dotplot settings to collaborators.
- You want one canonical workflow file that mirrors your GUI tutorial run.

## What You Learn

- Use one deterministic route from accession retrieval to pairwise dotplot artifact generation.
- Understand why high-sensitivity seed sampling (`step=1`, low word size) improves cDNA-vs-genomic block visibility.
- Trace parity between GUI actions and canonical workflow JSON execution.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md).
  - Reoccurs in: [Chapter 15: Prepare a Gibson specialist testing baseline (offline)](./15_gibson_specialist_testing_baseline.md).
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md).
  - Reoccurs in: [Chapter 15: Prepare a Gibson specialist testing baseline (offline)](./15_gibson_specialist_testing_baseline.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md).
  - Reoccurs in: [Chapter 15: Prepare a Gibson specialist testing baseline (offline)](./15_gibson_specialist_testing_baseline.md).
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
  - Status: reinforced from [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md).
  - Reoccurs in: no later chapter.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
  - Status: reinforced from [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Fetch GenBank accession `NM_001126241.3` as `tp73_cdna`.
2. Prepare `Human GRCh38 Ensembl 116`, retrieve gene `TP73` as `tp73_genomic`.
3. Open `tp73_cdna`, switch to `Dotplot map`, set pair mode against `tp73_genomic`, and compute with `word<=7`, `step=1`, `mismatches=0`.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json'
```

## Parameters That Matter

- `ComputeDotplot.word_size / step_bp / max_mismatches` (where used: operation 4)
  - Why it matters: These settings control sensitivity for dispersed exon-aligned signal in long genomic spans.
  - How to derive it: Start with `word=7`, `step=1`, `max_mismatches=0`; reduce to `word=6` only if additional sensitivity is needed.
- `ComputeDotplot.reference_seq_id` (where used: operation 4)
  - Why it matters: Pairwise mode requires explicit reference identity for y-axis mapping.
  - How to derive it: Use the exact ID emitted by the gene-extraction step (`tp73_genomic` in this chapter).

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'dotplot show tp73_cdna_vs_genomic_dotplot'
```

## Checkpoints

- Dotplot artifact is created and can be listed/shown by id (`tp73_cdna_vs_genomic_dotplot`).
- Pair-forward map shows exon-block style structure rather than one continuous full-length diagonal.
- GUI tutorial and workflow JSON describe the same parameter baseline.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json`
- Inspect this JSON file directly when you need full option-level detail.
