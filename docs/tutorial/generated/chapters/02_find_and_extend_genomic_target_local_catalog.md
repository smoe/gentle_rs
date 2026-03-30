# Find and extend the right genomic target (local catalog)

- Chapter id: `find_and_extend_genomic_target_local_catalog`
- Tier: `core`
- Example id: `prepare_extract_extend_localproject_gene`
- Source example: `docs/examples/workflows/prepare_extract_extend_localproject_gene.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Use assets/genomes.json + annotation-driven gene targeting, then extend the anchored sequence directly in the DNA window.

Many investigations start with "which exact genomic interval should I look at" rather than "which cloning operation should I run". This chapter introduces the practical sequence-selection loop: choose a prepared genome from `assets/genomes.json`, narrow candidates with annotation-backed gene filtering, extract the target interval, and then widen context with 5'/3' anchor extension as needed for promoter/splice/regulatory interpretation. The same workflow also clarifies where GenBank/EMBL imports fit: exported GenBank/EMBL records are loadable directly, and GenBank entries carrying NCBI `ACCESSION ... REGION:` metadata can be mapped to genome coordinates and extended by the same anchor controls.

## When This Routine Is Useful

- You need to start from a biologically relevant gene locus instead of an arbitrary coordinate slice.
- You imported a GenBank/EMBL record and want to recover surrounding genomic context around that entry.
- You need to add upstream/downstream bases iteratively while preserving deterministic provenance.

## What You Learn

- Understand that `assets/genomes.json` defines selectable genome sources for preparation/retrieval.
- Use annotation-backed gene filtering to find target regions reproducibly.
- Extend anchored sequences in GUI/CLI without losing provenance or sequence lineage.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md).
  - Reoccurs in: [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md).
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md).
  - Reoccurs in: [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md).
  - Reoccurs in: [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md).
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
  - Status: introduced in this chapter.
  - Reoccurs in: [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).

## GUI First

1. Open `File -> Prepare Reference Genome...`, set catalog to `assets/genomes.json`, choose `LocalProject`, and prepare it.
2. Open `File -> Retrieve Genome Sequence...`, use `Gene filter` (regex) to narrow annotation hits, pick one gene match, and extract it.
3. In the resulting DNA window, use the `Extend 5'` / `Extend 3'` anchor controls (next to `Genome anchor`) to add flanking context.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/prepare_extract_extend_localproject_gene.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/prepare_extract_extend_localproject_gene.json'
```

## Parameters That Matter

- `PrepareGenome.genome_id / catalog_path / cache_dir` (where used: operation 1)
  - Why it matters: These fields define which reference source is prepared and where deterministic cache/index files are written.
  - How to derive it: Start with `assets/genomes.json`; use `LocalProject` for offline walkthroughs, then switch to your organism/build.
- `Retrieve dialog Gene filter (regex) and Top matches` (where used: GUI extraction step)
  - Why it matters: Filtering keeps candidate lists interpretable on large annotations and avoids selecting the wrong locus.
  - How to derive it: Use exact regex (`^GENE$`) when known; broaden (`^GENE_PREFIX`) when exploring families.
  - Omit when: Only omit filtering for tiny local annotations where the candidate list is already small.
- `ExtendGenomeAnchor.side / length_bp / output_id` (where used: operation 3 and DNA-window anchor controls)
  - Why it matters: Defines which biological flank is added and by how much; output IDs keep downstream comparisons clear.
  - How to derive it: Choose `5p` for upstream-context questions and `3p` for downstream-context questions; start with 200-2000 bp and iterate.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- genomes list --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes genes LocalProject --catalog assets/genomes.json --cache-dir data/genomes --filter '^etp' --limit 20
cargo run --bin gentle_cli -- genomes extend-anchor local_etpc 5p 500 --output-id local_etpc_ext5_more --catalog assets/genomes.json --cache-dir data/genomes
```

## Checkpoints

- PrepareGenome succeeds for `LocalProject` without network access.
- ExtractGenomeGene produces `local_etpc` from annotation-backed lookup.
- ExtendGenomeAnchor produces `local_etpc_ext5` with widened genomic interval provenance.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/prepare_extract_extend_localproject_gene.json`
- Inspect this JSON file directly when you need full option-level detail.
