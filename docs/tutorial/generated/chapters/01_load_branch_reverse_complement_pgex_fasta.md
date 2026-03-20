# Load FASTA, branch, and reverse-complement

- Chapter id: `load_branch_reverse_complement_pgex_fasta`
- Tier: `core`
- Example id: `load_branch_reverse_complement_pgex_fasta`
- Source example: `docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Start with a simple cloning-prep routine to learn product IDs and sequence lineage.

In real cloning projects, the first question is often whether a sequence in hand is represented correctly and whether you can derive alternative views without losing provenance. This routine is the smallest possible example of that practice: load one sequence, create an explicit branch, and generate the reverse complement as a separate product.

## When This Routine Is Useful

- You received a FASTA plasmid or amplicon sequence and want a reproducible project starting point.
- You need a reverse-complemented working copy but want to preserve the original sequence identity.
- You are onboarding to GENtle and need to understand how derived IDs are created and reused.

## What You Learn

- Run one canonical workflow from source JSON and interpret the result IDs.
- Understand that branch and reverse-complement operations create explicit derivative sequences.
- Recognize the shared-engine contract across GUI, CLI, and shell entry points.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: introduced in this chapter.
  - Reoccurs in: [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: introduced in this chapter.
  - Reoccurs in: [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: introduced in this chapter.
  - Reoccurs in: [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).

## GUI First

1. Open GENtle and load `test_files/pGEX_3X.fa` via `File -> Open`.
2. In the DNA window, create a branch copy from the loaded sequence (Branch action).
3. Apply reverse-complement to the branch and confirm a new sequence entry appears in lineage/table views.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json'
```

## Parameters That Matter

- `LoadFile.path` (where used: operation 1)
  - Why it matters: Defines the exact source sequence. Wrong file path means a different biological starting point.
  - How to derive it: Use the file you just loaded in the GUI for this routine (`test_files/pGEX_3X.fa`).
- `Branch.output_id / ReverseComplement.output_id` (where used: operations 2 and 3)
  - Why it matters: Stable IDs make downstream commands unambiguous.
  - How to derive it: Pick short descriptive IDs tied to intent (e.g., branch copy vs reverse-complement product).
  - Omit when: You can omit IDs only when you accept auto-generated names.

## Checkpoints

- Workflow executes without warnings or errors.
- Derived sequence IDs include a branch and reverse-complement product.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`
- Inspect this JSON file directly when you need full option-level detail.
