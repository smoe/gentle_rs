# Prepare a reference genome cache (online)

- Chapter id: `prepare_reference_genome_online`
- Tier: `online`
- Example id: `prepare_reference_genome_online`
- Source example: `docs/examples/workflows/prepare_reference_genome_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.

Document network-dependent preparation without destabilizing offline default checks.

Reference-genome preparation is crucial for genome-anchored cloning interpretation, but it depends on online resources and local tool setup. This chapter keeps that path explicit and opt-in so default tutorial checks remain robust.

## When This Routine Is Useful

- You need genome-anchored extraction around gene/promoter context.
- You want to prepare local cache/index assets for repeated anchor operations.
- You need to understand which routines are intentionally online-only in CI defaults.

## What You Learn

- Recognize which tutorial flows require network access and external tools.
- Use `GENTLE_TEST_ONLINE` to opt into online chapter execution.
- Preserve offline CI reliability while still documenting online capabilities.

## Concepts and Recurrence

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md).
  - Reoccurs in: [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./13_tp73_cdna_genomic_dotplot_online.md).
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
  - Status: introduced in this chapter.
  - Reoccurs in: [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./13_tp73_cdna_genomic_dotplot_online.md).
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
  - Status: reinforced from [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md).
  - Reoccurs in: [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./13_tp73_cdna_genomic_dotplot_online.md).

## GUI First

1. Open prepared-reference controls from the GUI menus.
2. Select the target genome and start preparation with explicit cache settings.
3. Confirm prepared status in the GUI before attempting extraction workflows.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/prepare_reference_genome_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/prepare_reference_genome_online.json'
```

## Parameters That Matter

- `PrepareGenome.genome_id` (where used: operation 1)
  - Why it matters: Selects the exact reference build and annotation set.
  - How to derive it: Choose the genome build matching your experimental system and downstream coordinate system.
- `PrepareGenome.catalog_path / cache_dir` (where used: operation 1)
  - Why it matters: Controls source catalog and local cache destination.
  - How to derive it: Use repository defaults unless your environment requires custom catalogs or cache locations.

## Follow-up Commands

```bash
GENTLE_TEST_ONLINE=1 cargo run --bin gentle_examples_docs -- tutorial-generate
```

## Checkpoints

- Genome preparation runs only when GENTLE_TEST_ONLINE is enabled.
- Offline generation still emits the chapter with execution status noted.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/prepare_reference_genome_online.json`
- Inspect this JSON file directly when you need full option-level detail.
