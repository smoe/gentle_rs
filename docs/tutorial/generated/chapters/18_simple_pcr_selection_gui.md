# Simple PCR From a Selected Core Region

- Chapter id: `simple_pcr_selection_gui`
- Tier: `core`
- Example id: `simple_pcr_selection_gui`
- Source example: `docs/examples/workflows/simple_pcr_selection_gui.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Open a ready local TP73 sequence and walk through the smallest useful PCR story: select the core ROI, limit primer distance from the core, and cap the amplicon length.

This chapter is intentionally minimal. The tutorial project opens one stable local TP73 locus so you can stay in the GUI, paint or select one region of interest, launch `Simple PCR from selection`, and reason about primer placement in beginner terms before moving on to the richer batch-PCR chapter.

## When This Routine Is Useful

- You want the shortest possible path from one selected region to one PCR design attempt.
- You want to learn the meaning of core ROI, maximum primer distance from the core, and maximum amplicon length without queueing or painting multiple regions.
- You want one offline tutorial project that opens directly into local TP73 context before using the PCR Designer.

## What You Learn

- Start a simple PCR directly from one sequence selection.
- Interpret the beginner PCR controls as deterministic flank-window constraints around the selected core ROI.
- Review saved primer pairs in terms of amplicon length and left/right distance from the core ROI.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md), [Chapter 17: TP53 UniProt domain mapping and feature-coding DNA query (online)](./17_tp53_uniprot_projection_online.md).
  - Reoccurs in: no later chapter.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Open `File -> Open Tutorial Project... -> Core -> 18. Simple PCR From a Selected Core Region`.
2. In the opened TP73 project, keep the DNA map in linear mode and drag over one short core region that must be included in the amplicon.
3. Right-click the selection and choose `Simple PCR from selection`.
4. In `PCR Designer`, adjust `max primer distance from core` and `max amplicon`, then click `Apply simple flank windows` if you changed the distance.
5. Run `Design Primer Pairs` and inspect the in-panel primer report preview for left/right distance from the core ROI and whether the pair cleanly flanks the core.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/simple_pcr_selection_gui.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/simple_pcr_selection_gui.json'
```

## Parameters That Matter

- `core ROI` (where used: map selection and PCR Designer starter block)
  - Why it matters: This is the part of the biology that must definitely be included in the PCR product.
  - How to derive it: Select only the indispensable sequence interval, not the full desired amplicon.
- `max primer distance from core` (where used: PCR Designer starter block)
  - Why it matters: This sets how far upstream and downstream primers may move away from the required core region.
  - How to derive it: Choose the largest flank width you are comfortable allowing for primer search on each side of the core ROI.
- `max amplicon` (where used: Design primer pairs form)
  - Why it matters: This caps product length so the simple-PCR result stays experimentally practical.
  - How to derive it: Set the longest acceptable product length for the assay before running primer design.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/simple_pcr_selection_gui.json
cargo run --bin gentle_cli -- shell 'ui open pcr-design'
cargo run --bin gentle_cli -- shell 'primers list-reports'
```

## Checkpoints

- The tutorial project opens with local TP73 sequence `tp73_locus` already loaded.
- A non-empty selection exposes `Simple PCR from selection` in the DNA-window context menu.
- The PCR Designer shows the `Simple PCR starter` block.
- After primer design, the report preview shows left/right distance from the core ROI and whether the pair cleanly flanks the core.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/simple_pcr_selection_gui.json`
- Inspect this JSON file directly when you need full option-level detail.
