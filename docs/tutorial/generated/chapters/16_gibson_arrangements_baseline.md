# Gibson Arrangements Starter Project (offline)

- Chapter id: `gibson_arrangements_baseline`
- Tier: `core`
- Example id: `gibson_specialist_testing_baseline`
- Source example: `docs/examples/workflows/gibson_specialist_testing_baseline.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Build the same stable Gibson starter project, but open directly into the arrangement-focused guide for container, arrangement, and gel-export review.

This chapter reuses the deterministic pGEX + insert starter baseline from the Gibson specialist tutorial, but changes the guide handoff so `File -> Open Tutorial Project...` can land directly on arrangement inspection and gel-export workflow instead of overlap/primer setup.

## When This Routine Is Useful

- You want to inspect the arrangement that Gibson apply creates without first navigating through the earlier specialist tutorial.
- You want a reproducible starter state for checking singleton output containers and arrangement-level gel export.
- You want one offline tutorial-project entry that opens directly on the arrangement guide in Help/Tutorial.

## What You Learn

- Use one executable starter project to reach the arrangement walkthrough directly.
- Recognize that the same deterministic Gibson baseline can support more than one tutorial guide.
- Replay the same arrangement-focused setup from GUI and CLI without changing the biological inputs.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md).
  - Reoccurs in: no later chapter.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md).
  - Reoccurs in: no later chapter.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Open `File -> Open Tutorial Project...` and choose `Gibson Arrangements Starter Project (offline)`.
2. Confirm the opened project contains `gibson_destination_pgex` (circular) and `gibson_insert_demo` (linear).
3. The Help window should open automatically on `Gibson Arrangements Tutorial`; continue there from `Step 2` onward.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json'
```

## Parameters That Matter

- `destination seq_id = gibson_destination_pgex` (where used: Gibson destination lane in the arrangement walkthrough)
  - Why it matters: Keeps the original vector lane stable across the arrangement guide and CLI replay.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `insert seq_id = gibson_insert_demo` (where used: Insert lane in the arrangement walkthrough)
  - Why it matters: Provides the same deterministic insert lane that the Gibson arrangement is expected to reference.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `workflow example = gibson_specialist_testing_baseline` (where used: `Open Tutorial Project...` and CLI workflow replay)
  - Why it matters: The arrangement starter project intentionally reuses the same stable baseline as the Gibson specialist starter project.
  - How to derive it: Select the chapter from `Open Tutorial Project...` or run the canonical workflow JSON directly.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json
```

## Checkpoints

- The tutorial project opens with exactly the expected destination and insert IDs already loaded.
- The Help window lands on the arrangement tutorial rather than the specialist testing tutorial.
- The arrangement guide can start from this baseline without any additional sequence imports.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/gibson_specialist_testing_baseline.json`
- Inspect this JSON file directly when you need full option-level detail.
