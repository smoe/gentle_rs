# TP53 UniProt protein mapping expert view (online)

- Chapter id: `tp53_uniprot_projection_online`
- Tier: `online`
- Example id: `tp53_uniprot_projection_online`
- Source example: `docs/examples/workflows/tp53_uniprot_projection_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.

Build a TP53 locus project, fetch UniProt P04637, and open/export the shared protein-mapping expert view from one persisted projection route.

This chapter demonstrates the new UniProt projection UX. Instead of importing a curated panel JSON, you fetch the reviewed UniProt TP53 entry, project its reference-protein intervals onto one extracted TP53 locus, and inspect the result through the same shared isoform/protein expert canvas used by GUI and CLI. The emphasis is on one persisted projection record that can be reopened, rendered, and audited later.

## When This Routine Is Useful

- You want to compare one gene locus against the reviewed UniProt reference protein without opening a first-class protein sequence window.
- You want a deterministic TP53 example that maps transcript/CDS geometry onto UniProt protein coordinates.
- You want a replayable expert-view SVG export that comes from stored projection state rather than manual figure editing.

## What You Learn

- Use one persisted UniProt projection as the canonical bridge from reviewed protein annotation to locus-level transcript/CDS geometry.
- Understand that the GUI protein expert is a thin presentation layer over `gentle.uniprot_genome_projections.v1`, not a separate GUI-only mapping model.
- Export the exact same UniProt protein-mapping view through the same shared engine route from either the GUI specialist or CLI.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md).
  - Reoccurs in: no later chapter.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
  - Status: reinforced from [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).
  - Reoccurs in: no later chapter.
- **UniProt Projection Mapping** (`uniprot_projection_mapping`): Reviewed UniProt protein annotations can be projected onto transcript/CDS geometry and reopened as one persisted expert-view artifact.
  - Status: introduced in this chapter.
  - Reoccurs in: no later chapter.
- **Expert View Parity** (`expert_view_parity`): The same expert-view payloads should be inspectable and renderable from GUI, CLI, and other adapters without frontend-only projection logic.
  - Status: introduced in this chapter.
  - Reoccurs in: no later chapter.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
  - Status: reinforced from [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md).
  - Reoccurs in: no later chapter.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
  - Status: reinforced from [Chapter 7: Guide oligo export (CSV + protocol)](./07_guides_export_csv_and_protocol.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Prepare `Human GRCh38 Ensembl 116` and extract gene `TP53` into `grch38_tp53`.
2. Open `File -> UniProt Mapping...`, fetch `P04637`, keep `entry_id=P04637`, choose sequence `grch38_tp53`, and run `Project To Sequence`.
3. Use `Render Protein Mapping SVG...` to export the stored projection directly, or press `Open Protein Expert` first if you want to inspect the transcript/protein geometry before exporting.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp53_uniprot_projection_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp53_uniprot_projection_online.json'
```

## Parameters That Matter

- `FetchUniprotSwissProt.query / entry_id` (where used: operation 3)
  - Why it matters: Pins the reviewed UniProt entry that will contribute the protein coordinate system and interval annotations.
  - How to derive it: Use a stable accession or reviewed UniProt id. For the canonical TP53 example, `P04637` is the expected reviewed entry.
- `ProjectUniprotToGenome.projection_id / transcript_id` (where used: operation 4)
  - Why it matters: The projection id becomes the durable handle reused by the GUI recent-projection list and the shared expert-render route. Leaving `transcript_id=null` lets the UniProt Ensembl transcript xrefs drive the full TP53 projection set.
  - How to derive it: Use a stable project-specific id such as `tp53_uniprot_p04637`. Supply `transcript_id` later only when you intentionally want a single-transcript projection.
- `RenderFeatureExpertSvg.target / path` (where used: operation 5)
  - Why it matters: Confirms that the projection can be reopened through the shared feature-expert route and exported reproducibly.
  - How to derive it: Use `UniprotProjection` with the stored projection id and a stable project-relative output path such as `exports/tp53_uniprot_projection.svg`.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'uniprot projection-show tp53_uniprot_p04637'
cargo run --bin gentle_cli -- inspect-feature-expert grch38_tp53 uniprot-projection tp53_uniprot_p04637
cargo run --bin gentle_cli -- render-feature-expert-svg grch38_tp53 uniprot-projection tp53_uniprot_p04637 exports/tp53_uniprot_projection.svg
```

## Checkpoints

- UniProt fetch/import reports the reviewed TP53 entry and the projection operation stores `tp53_uniprot_p04637` in project metadata.
- Open Protein Expert shows one protein reference axis with projected transcript/CDS lanes from the TP53 locus.
- Protein-mapping SVG export succeeds directly from the UniProt specialist without requiring a separate protein sequence window.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/tp53_uniprot_projection_online.json`
- Inspect this JSON file directly when you need full option-level detail.
