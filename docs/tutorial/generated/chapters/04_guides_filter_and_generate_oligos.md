# Guide practical filtering and oligo generation

- Chapter id: `guides_filter_and_generate_oligos`
- Tier: `core`
- Example id: `guides_filter_and_generate_oligos`
- Source example: `docs/examples/workflows/guides_filter_and_generate_oligos.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Apply practical guide constraints and produce cloning-ready oligo candidates.

For CRISPR-style cloning, guide quality control is where many downstream failures are prevented. This routine demonstrates how to encode practical constraints directly in operations, then generate oligos from the passed candidates.

## When This Routine Is Useful

- You need a transparent guide filtering step before oligo ordering.
- You want to document why specific guides were excluded.
- You need deterministic oligo generation from a reusable guide set.

## What You Learn

- Create and filter guide sets through explicit engine operations.
- Generate oligo records from filtered guide candidates.
- Connect guide workflows to the same deterministic operation model used for cloning steps.

## Concepts and Recurrence

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md).
  - Reoccurs in: [Chapter 5: Digest -> Ligation -> ExtractRegion minimal slice](./05_digest_ligation_extract_region_minimal.md), [Chapter 8: Prepare a reference genome cache (online)](./08_prepare_reference_genome_online.md).
- **Guide Design Pipeline** (`guide_design_pipeline`): Guide sets can be created, filtered, expanded to oligos, and exported with protocol context.
  - Status: introduced in this chapter.
  - Reoccurs in: [Chapter 6: Guide oligo export (CSV + protocol)](./06_guides_export_csv_and_protocol.md).

## GUI First

1. Open the guides workflow controls in GENtle and create/import a guide set for a target region.
2. Apply practical filters (GC range, homopolymer limits, U6 terminator avoidance).
3. Generate oligos from passed guides and inspect the resulting oligo set IDs.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/guides_filter_and_generate_oligos.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/guides_filter_and_generate_oligos.json'
```

## Parameters That Matter

- `FilterGuidesPractical.config.gc_min / gc_max` (where used: operation 2)
  - Why it matters: GC bounds trade off stability and synthesis/efficiency behavior.
  - How to derive it: Start with literature-typical bounds (e.g., 0.3-0.7) and tighten per assay constraints.
- `GenerateGuideOligos.template_id` (where used: operation 3)
  - Why it matters: Template controls overhang/adaptor context for your cloning backbone.
  - How to derive it: Select the template matching your vector and cloning strategy.

## Checkpoints

- Guide set and passed guide set are present in metadata.
- Oligo set generation succeeds for passed guides.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/guides_filter_and_generate_oligos.json`
- Inspect this JSON file directly when you need full option-level detail.
