# Load pGEX and digest with BamHI/EcoRI

- Chapter id: `load_and_digest_pgex`
- Tier: `core`
- Example id: `load_and_digest_pgex`
- Source example: `docs/examples/workflows/load_and_digest_pgex.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Introduce restriction digest planning and deterministic fragment products.

Restriction digest is a core molecular cloning routine used for vector linearization, insert release, and diagnostic fragment checks. This chapter focuses on how digest parameters map to reproducible fragment sets that can later feed ligation or analysis steps.

## When This Routine Is Useful

- You want to verify expected restriction fragments before ordering primers or designing ligations.
- You need a reproducible digest baseline to compare against wet-lab gel expectations.
- You plan to reuse fragment IDs in later operations.

## What You Learn

- Execute digest operations and inspect created fragment IDs.
- Reason about multi-product lineage from one parent sequence.
- Identify why stable IDs matter for follow-up ligation/extraction steps.

## Concepts and Recurrence

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md).
  - Reoccurs in: [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md).
  - Reoccurs in: [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md).

## GUI First

1. Load `test_files/pGEX-3X.gb` in the GUI and inspect annotated features.
2. Run Digest from the DNA window using enzymes `BamHI` and `EcoRI`.
3. Review created fragment entries and confirm they are stored as independent sequence products.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_and_digest_pgex.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/load_and_digest_pgex.json'
```

## Parameters That Matter

- `Digest.enzymes` (where used: operation 2)
  - Why it matters: The enzyme set determines cut positions and resulting fragment repertoire.
  - How to derive it: Choose enzymes based on cloning strategy (diagnostic digest vs insert release vs vector opening).
- `Digest.output_prefix` (where used: operation 2)
  - Why it matters: Controls deterministic fragment ID namespace.
  - How to derive it: Use a short routine-specific prefix (e.g., `frag`, `d`, `eco_bam`).
  - Omit when: Omit only if auto-generated IDs are acceptable for ad hoc inspection.

## Checkpoints

- Digest operation completes and creates fragment sequence IDs.
- Fragment IDs are deterministic across repeated runs.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/load_and_digest_pgex.json`
- Inspect this JSON file directly when you need full option-level detail.
