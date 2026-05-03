# Promoter Design Artifact Slice (Offline Synthetic TP73 Locus)

- Chapter id: `promoter_design_artifact_slice_offline`
- Tier: `core`
- Example id: `promoter_design_artifact_slice_offline`
- Source example: `docs/examples/workflows/promoter_design_artifact_slice_offline.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Run a tiny TP73-like promoter locus through the Promoter design artifact set: alternative-promoter grouping, promoter evidence matrix, isoform promoter comparison, promoter expression evidence, TFBS score-track SVG, TFBS similarity JSON, and a component manifest.

This chapter is a release-signoff slice rather than a biological claim. The input is a synthetic 249 bp TP73-labeled locus with two transcripts sharing one 5' boundary, one alternative-start transcript, and local annotation evidence. The point is to verify that GENtle produces the portable artifacts a GUI user, CLI user, or ClawBio-style consumer can inspect without requiring an online genome fetch or a hand-built project.

## When This Routine Is Useful

- You want a fast offline sanity check for the Promoter design window before an internal release.
- You want to verify that transcript-derived promoter windows collapse by DNA span instead of stacking duplicate promoter symbols.
- You want a compact promoter evidence matrix with transcript, TFBS, variant, repeat, and CUT&RUN-style interval evidence.
- You want to compare common and differential promoter-region evidence between isoform starts of the same gene.
- You want to attach downstream transcript/gene expression evidence to promoter candidates without forcing GENtle to write the final story.

## What You Learn

- Validate the GUI-facing Promoter design controls against shared engine operations.
- Inspect the boundary between GENtle-owned component artifacts and downstream ClawBio/OpenClaw narrative assembly.
- Recognize the difference between transcript-level promoter interpretation and collapsed DNA-level promoter candidates.
- Compare common and unique promoter-region evidence between isoforms of the same gene.
- Attach expression-level evidence to promoter groups through a structured report instead of free-text reasoning.
- Use one synthetic input to exercise JSON and SVG exports deterministically.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md), [Chapter 17: TP53 UniProt domain mapping and feature-coding DNA query (online)](./17_tp53_uniprot_projection_online.md), [Chapter 18: Simple PCR From a Selected Core Region](./18_simple_pcr_selection_gui.md), [Chapter 19: Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)](./19_tp73_uniprot_projection_audit_cli.md).
  - Reoccurs in: no later chapter.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./03_load_and_digest_pgex.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 5: Guide practical filtering and oligo generation](./05_guides_filter_and_generate_oligos.md), [Chapter 6: Digest -> Ligation -> ExtractRegion minimal slice](./06_digest_ligation_extract_region_minimal.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Selection-first PCR batch primer design (offline)](./13_pcr_selection_batch_primer_pairs_offline.md), [Chapter 14: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./14_tp73_cdna_genomic_dotplot_online.md), [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md), [Chapter 16: Gibson Arrangements Starter Project (offline)](./16_gibson_arrangements_baseline.md), [Chapter 18: Simple PCR From a Selected Core Region](./18_simple_pcr_selection_gui.md).
  - Reoccurs in: no later chapter.
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
  - Status: reinforced from [Chapter 7: Guide oligo export (CSV + protocol)](./07_guides_export_csv_and_protocol.md), [Chapter 10: TP53 isoform architecture expert panel (online)](./10_tp53_isoform_architecture_online.md), [Chapter 17: TP53 UniProt domain mapping and feature-coding DNA query (online)](./17_tp53_uniprot_projection_online.md), [Chapter 19: Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence (CLI Tutorial)](./19_tp73_uniprot_projection_audit_cli.md).
  - Reoccurs in: no later chapter.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.
  - Status: reinforced from [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Open `docs/examples/assets/tp73_promoter_artifact_demo.gb` via `File -> Open Sequence...`.
2. Open `Promoter design` from the `TP73` gene or one of the `TP73-demo-*` mRNA features.
3. Set `Gene label` to `TP73`, `promoter upstream bp` to `40`, and `promoter downstream bp` to `15`.
4. Click `Annotate promoter windows`, then `Compare alternative promoters`; confirm that three transcript-level interpretations collapse into two DNA-level promoter windows.
5. Click `Build evidence matrix`; confirm the shared promoter row reports `2 tx` and that evidence kinds include promoter geometry, transcript support, promoter annotation, TFBS, variant, repeat, and CUT&RUN-style overlap evidence.
6. Run the isoform promoter comparison; confirm that the shared TSS transcripts and alternative-start transcript are compared as separate promoter groups with differential evidence signatures.
7. Load or paste expression rows for the TP73 demo transcripts; confirm that expression evidence attaches to the matching promoter groups rather than becoming a GUI-only note.
8. Set TF motifs to `SP1,TP53,TP63,TP73`, run `Show TF score tracks`, then export `TF score tracks SVG...` for the visual artifact.
9. Set TFBS similarity anchor to `SP1`, compare against `TP53,TP63,TP73,CTCF`, run `Show TFBS similarity ranking`, then export the JSON ranking.
10. Inspect the component manifest to see which JSON/SVG artifacts were produced; downstream tools can choose their own presentation order.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/promoter_design_artifact_slice_offline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/promoter_design_artifact_slice_offline.json'
```

## Parameters That Matter

- `promoter_upstream_bp=40 / promoter_downstream_bp=15` (where used: AnnotatePromoterWindows, SummarizeAlternativePromoterComparison, SummarizePromoterEvidenceMatrix, SummarizeIsoformPromoterComparison, SummarizePromoterExpressionEvidence)
  - Why it matters: The synthetic locus is deliberately tiny; these values create readable promoter windows around TSS positions 101 and 141.
  - How to derive it: Use the fixed tutorial values so the shared TSS pair collapses to local span `60..116` and the alternative start produces `100..156`.
- `expression rows for `ENSTTP73DEMO1`, `ENSTTP73DEMO2`, and `ENSTTP73DEMO3`` (where used: SummarizePromoterExpressionEvidence)
  - Why it matters: Expression evidence is deliberately externalized as rows/artifact references so GENtle can attach it to promoter groups without inventing a biological conclusion.
  - How to derive it: Use transcript IDs from the synthetic mRNA features; real workflows would supply RNA-seq, qPCR, or ClawBio-provided abundance rows.
- `target span `60..158`` (where used: SummarizeTfbsScoreTracks, RenderTfbsScoreTracksSvg, SummarizeTfbsTrackSimilarity)
  - Why it matters: This span covers both promoter windows plus the synthetic SP1 and TP73-like TFBS sites.
  - How to derive it: Use the coordinate range printed by the workflow or the Promoter design score-track range seeded from the evidence rows.
- `motifs `SP1,TP53,TP63,TP73` and similarity candidates `TP53,TP63,TP73,CTCF`` (where used: TF score tracks and TFBS similarity ranking)
  - Why it matters: The motif set keeps the demo close to TP73/p53-family promoter reasoning while still producing a compact ranking table.
  - How to derive it: Use the exact tokens in this chapter; they resolve through GENtle's shared local JASPAR query layer.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/promoter_design_artifact_slice_offline.json
cargo run --bin gentle_cli -- shell 'features promoter-evidence-matrix tp73_promoter_artifact_demo --gene-label TP73 --promoter-upstream-bp 40 --promoter-downstream-bp 15 --path artifacts/tp73_promoter_artifact_demo.evidence_matrix.json'
cargo run --bin gentle_cli -- shell 'features promoter-isoform-comparison tp73_promoter_artifact_demo --gene-label TP73 --promoter-upstream-bp 40 --promoter-downstream-bp 15 --path artifacts/tp73_promoter_artifact_demo.isoform_promoter_comparison.json'
cargo run --bin gentle_cli -- shell 'features promoter-expression-evidence tp73_promoter_artifact_demo --gene-label TP73 --promoter-upstream-bp 40 --promoter-downstream-bp 15 --source-label synthetic_demo --expression-json {"transcript_id":"ENSTTP73DEMO1","value":18.0,"unit":"a.u."} --path artifacts/tp73_promoter_artifact_demo.promoter_expression_evidence.json'
cargo run --bin gentle_cli -- shell 'features tfbs-track-similarity tp73_promoter_artifact_demo --anchor SP1 --candidate TP53 --candidate TP63 --candidate TP73 --candidate CTCF --range 60..158 --score-kind llr_background_tail_log10 --path artifacts/tp73_promoter_artifact_demo.tfbs_similarity.json'
```

## Checkpoints

- The workflow executes offline without warnings beyond the expected alternative-promoter collapse note.
- `alternative_promoters.json` reports `transcript_window_count=3` and `collapsed_window_count=2`.
- `evidence_matrix.json` reports two promoter candidates and includes `cutrun_peak_overlap`, `repeat_context`, `tfbs_annotation`, and `variant_overlap` among observed evidence kinds.
- `isoform_promoter_comparison.json` reports two promoter groups and surfaces differential evidence signatures for the shared versus alternative-start promoter.
- `promoter_expression_evidence.json` assigns three synthetic expression rows to the two promoter groups.
- `promoter_artifact_manifest.json` reports all required promoter component artifacts as present.
- `tfbs_score_tracks.svg` is written and opens as a compact promoter score-track figure.
- `tfbs_similarity.json` ranks four candidates against SP1 using `smoothed_spearman`.

## Retained Outputs

- [`artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.alternative_promoters.json`](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.alternative_promoters.json)
- [`artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.evidence_matrix.json`](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.evidence_matrix.json)
- [`artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.isoform_promoter_comparison.json`](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.isoform_promoter_comparison.json)
- [`artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.promoter_artifact_manifest.json`](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.promoter_artifact_manifest.json)
- [`artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.promoter_expression_evidence.json`](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.promoter_expression_evidence.json)
- [`artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.tfbs_score_tracks.svg`](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.tfbs_score_tracks.svg)
- [`artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.tfbs_similarity.json`](../artifacts/promoter_design_artifact_slice_offline/artifacts/tp73_promoter_artifact_demo.tfbs_similarity.json)

## Canonical Source

- Workflow file: `docs/examples/workflows/promoter_design_artifact_slice_offline.json`
- Inspect this JSON file directly when you need full option-level detail.
