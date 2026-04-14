# TP53 UniProt domain mapping and feature-coding DNA query (online)

- Chapter id: `tp53_uniprot_projection_online`
- Tier: `online`
- Example id: `tp53_uniprot_projection_online`
- Source example: `docs/examples/workflows/tp53_uniprot_projection_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.

Build a TP53 locus project, fetch UniProt P04637, map its domains onto the locus, and query which genomic DNA plus exon or exon pair encode one mapped feature.

This chapter demonstrates the UniProt projection workflow as an inspectable bridge from protein annotation back to genomic DNA. Instead of importing a curated protein panel JSON, you fetch the reviewed UniProt TP53 entry, project its reference-protein intervals onto one extracted TP53 locus, inspect the mapped domains through the shared expert canvas, and then query one mapped feature such as `DNA-binding` to recover the exact spliced genomic coding DNA, exon attribution, and an optional translation-speed-oriented coding alternative.

## When This Routine Is Useful

- You want to compare one gene locus against reviewed UniProt domains or regions without opening a first-class protein sequence window.
- You want to know which spliced genomic DNA and which exon or exon pair encode one mapped feature such as TP53 DNA-binding.
- You want a deterministic TP53 example that keeps one persisted projection record reusable for expert rendering, audit, and follow-up DNA queries.

## What You Learn

- Use one persisted UniProt projection as the canonical bridge from reviewed protein annotation to locus-level transcript/CDS geometry and back to coding DNA.
- Understand that both the protein expert and the feature-coding DNA query are thin views over stored engine state, not separate GUI-only mapping models.
- Recover exact genomic coding DNA plus exon attribution for one mapped protein feature, and compare it with an optional translation-speed-oriented codon choice.

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
- **Feature Coding-DNA Attribution** (`feature_coding_dna_attribution`): A persisted UniProt projection can be queried for the exact coding DNA, exon attribution, and splice-junction exon pairs that encode one mapped protein feature.
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
2. Open `File -> Protein Evidence...`, fetch `P04637`, keep `entry_id=P04637`, choose sequence `grch38_tp53`, and run `Project To Sequence`.
3. Inspect the stored projection with `Open Protein Expert` or export it with `Render Protein Mapping SVG...` so you can verify how UniProt domains/regions landed on the TP53 transcripts.
4. In `Feature coding DNA query`, enter `DNA-binding`, leave `feature transcript` empty unless you want to pin one isoform, choose `mode=both`, and press `Query Coding DNA`.
5. Read the result panel to see the amino-acid span, genomic coding DNA, optional translation-speed optimized DNA, and the reported exon or exon pair for each matching transcript feature span.

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
- `query_uniprot_feature_coding_dna.feature_query / transcript_id` (where used: GUI follow-up or shell follow-up after projection)
  - Why it matters: The feature query chooses which mapped UniProt interval you want to trace back to coding DNA, while `transcript_id` narrows the report to one isoform when the projection stored multiple TP53 transcripts.
  - How to derive it: Use a case-insensitive mapped feature substring such as `DNA-binding`, `activation`, or `DOMAIN`. Leave `transcript_id` empty until you need to pin one transcript.
- `query_uniprot_feature_coding_dna.mode / translation_speed_profile` (where used: GUI follow-up or shell follow-up after projection)
  - Why it matters: Controls whether you inspect only the exact genomic coding DNA or also a preferred-codon translation-speed-oriented alternative for the same amino-acid interval.
  - How to derive it: Use `mode=both` for this tutorial so you can compare the genomic sequence with the optimized alternative. Keep the speed profile on `Auto` in the GUI or choose `human` explicitly in shell/CLI for TP53.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'uniprot projection-show tp53_uniprot_p04637'
cargo run --bin gentle_cli -- inspect-feature-expert grch38_tp53 uniprot-projection tp53_uniprot_p04637
cargo run --bin gentle_cli -- render-feature-expert-svg grch38_tp53 uniprot-projection tp53_uniprot_p04637 exports/tp53_uniprot_projection.svg
cargo run --bin gentle_cli -- shell 'uniprot feature-coding-dna tp53_uniprot_p04637 DNA-binding --mode both --speed-profile human'
```

## Checkpoints

- UniProt fetch/import reports the reviewed TP53 entry and the projection operation stores `tp53_uniprot_p04637` in project metadata.
- Open Protein Expert shows one protein reference axis with projected transcript/CDS lanes from the TP53 locus.
- Protein-mapping SVG export succeeds directly from the UniProt specialist without requiring a separate protein sequence window.
- The feature-coding DNA query reports the coding-strand DNA for `DNA-binding` together with transcript-specific exon attribution, including an exon pair when the feature spans a splice junction.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/tp53_uniprot_projection_online.json`
- Inspect this JSON file directly when you need full option-level detail.
