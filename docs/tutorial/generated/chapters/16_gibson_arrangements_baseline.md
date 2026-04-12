# Gibson Arrangements Starter Project (offline)

- Chapter id: `gibson_arrangements_baseline`
- Tier: `core`
- Example id: `gibson_arrangements_baseline`
- Source example: `docs/examples/workflows/gibson_arrangements_baseline.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Open directly into an arrangement-ready Gibson result with vector, insert, assembled product, and the stored three-lane lane set already present.

This chapter now uses its own deterministic workflow example instead of reusing the Gibson specialist starter. `File -> Open Tutorial Project...` builds the pGEX + insert starter, applies the canonical single-insert Gibson plan, and then opens the arrangement guide so the user can inspect singleton outputs, stored arrangements, and gel export without first repeating the cloning step.

## When This Routine Is Useful

- You want to inspect the arrangement that Gibson apply creates without first navigating through the earlier Gibson-specialist apply walkthrough.
- You want a reproducible starter state for checking singleton output containers, the assembled product, and arrangement-level gel export.
- You want one offline tutorial-project entry that opens directly on the arrangement guide in Help/Tutorial.

## What You Learn

- Use one executable starter project to reach an arrangement-ready walkthrough directly.
- Recognize the separation between the Gibson specialist apply tutorial and the downstream arrangement/gel-inspection tutorial.
- Replay the same arrangement-focused setup from GUI and CLI without requiring a second manual Gibson apply.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md).
  - Reoccurs in: [Chapter 17: TP53 UniProt domain mapping and feature-coding DNA query (online)](./17_tp53_uniprot_projection_online.md).
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md).
  - Reoccurs in: no later chapter.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Open `File -> Open Tutorial Project...` and choose `Gibson Arrangements Starter Project (offline)`.
2. Confirm the opened project contains `gibson_destination_pgex` (circular), `gibson_insert_demo` (linear), and the assembled product `gibson_destination_pgex_with_gibson_insert_demo`.
3. The Help window should open automatically on `Gibson Arrangements Tutorial`; continue there from `Step 1` onward.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_arrangements_baseline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/gibson_arrangements_baseline.json'
```

## Parameters That Matter

- `destination seq_id = gibson_destination_pgex` (where used: Stored vector lane in the arrangement walkthrough)
  - Why it matters: Keeps the original vector lane stable across the arrangement guide and CLI replay.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `insert seq_id = gibson_insert_demo` (where used: Insert lane in the arrangement walkthrough)
  - Why it matters: Provides the same deterministic insert lane that the Gibson arrangement is expected to reference.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `assembled product seq_id = gibson_destination_pgex_with_gibson_insert_demo` (where used: Stored Gibson result already present when the arrangement tutorial opens)
  - Why it matters: Makes the arrangement tutorial disjunct from the separate Gibson-specialist apply walkthrough while keeping one deterministic product identity.
  - How to derive it: The chapter workflow applies the canonical single-insert Gibson plan before the guide opens.
- `workflow example = gibson_arrangements_baseline` (where used: `Open Tutorial Project...` and CLI workflow replay)
  - Why it matters: The arrangement starter now has its own canonical example so the tutorial can begin with the cloned result and stored lane arrangement already available.
  - How to derive it: Select the chapter from `Open Tutorial Project...` or run the canonical workflow JSON directly.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_arrangements_baseline.json
```

## Checkpoints

- The tutorial project opens with the expected destination, insert, and assembled-product IDs already loaded.
- The Help window lands on the arrangement tutorial rather than the specialist testing tutorial.
- The arrangement guide can start from this baseline without an additional manual Gibson apply.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/gibson_arrangements_baseline.json`
- Inspect this JSON file directly when you need full option-level detail.
