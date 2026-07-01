# Tutorial Coverage Audit

Phase 0 audit for the tutorial presentation overhaul. This page is a
review gate, not a schema migration. Do not assign decimal ids, rename files,
or regroup menus until the decisions at the end of this page have been signed
off.

Generated from repository inspection on 2026-05-22.

## Snapshot

| Surface | Count | Notes |
| --- | ---: | --- |
| Catalog entries | 38 | `docs/tutorial/catalog.json` |
| Review manifest entries | 38 | Matches current catalog ids |
| Tutorial source JSON files | 38 | Filename prefixes run `01` through `33`, with duplicate prefixes |
| Executable manifest chapters | 20 | `docs/tutorial/manifest.json` |
| Generated chapter Markdown files | 20 | Current filenames use `01` through `19`, then `24` |
| Top-level Markdown under `docs/tutorial/` | 24 | Includes `README.md`, five uncataloged tutorial-like pages, and cataloged guided pages |
| Tutorial Markdown outside `docs/tutorial/` | 1 | `docs/tutorial/01-01_agent_interfaces.md`, cataloged as `agent_interfaces` |

Current checks pass before the overhaul:

- `cargo run --bin gentle_examples_docs -- tutorial-catalog-check`
- `cargo run --bin gentle_examples_docs -- tutorial-manifest-check`
- `cargo run --bin gentle_examples_docs -- tutorial-check`

## Existing Numbering Conflicts

The source JSON filename prefixes are not a trustworthy ordering key:

| Prefix | Files |
| --- | --- |
| `02` | `02_tp73_promoter_luciferase_gui.json`, `02_vkorc1_variant_followup_expert_gui.json` |
| `04` | `04_agent_interfaces.json`, `04_gibson_specialist_testing_gui.json` |
| `06` | `06_gibson_physical_rack_gui.json`, `06_sequencing_confirmation_trace_cli.json`, `06_simple_pcr_selection_gui.json` |
| `10` | `10_geneart_external_service_handoff_gui_cli.json`, `10_metabion_external_service_handoff_gui_cli.json` |

One filename is also semantically stale:

- `docs/tutorial/sources/08-04_vkorc1_warfarin_promoter_luciferase_gui.json` now describes
  `vkorc1_warfarin_promoter_luciferase_gui`, the maintained VKORC1 reporter
  tutorial.

These conflicts are exactly why the decimal-id migration should number source
units rather than preserving current filename prefixes.

## Coverage Findings

There are no catalog entries without source JSON files, and no source JSON ids
without catalog entries. The current gap is Markdown coverage: five
top-level Markdown pages are linked or tutorial-like but not registered in the
catalog/review manifest.

| Markdown path | Current coverage | Proposed Phase 1 classification | Notes |
| --- | --- | --- | --- |
| `docs/tutorial/02-02_stateless_sequence_inspection_gui_cli.md` | Not cataloged | Add as active guided source unit | Already linked from the tutorial README; covers state-optional sequence inspection across GUI, CLI, and ClawBio. |
| `docs/tutorial/08-01_tfbs_similarity_ranking_gui.md` | Not cataloged | Add as active guided source unit | Already linked from the tutorial README; covers GUI sign-off for TFBS similarity ranking. |
| `docs/tutorial/tp73_promoter_luciferase_gui.md` | Not cataloged | Treat as deprecated pointer, not active unit | Resolved after sign-off: moved to `docs/tutorial/archive/tp73_promoter_luciferase_gui.md`. |
| `docs/tutorial/tp73_promoter_luciferase_showcase.md` | Not cataloged | Treat as deprecated pointer, not active unit | Resolved after sign-off: moved to `docs/tutorial/archive/tp73_promoter_luciferase_showcase.md`. |
| `docs/tutorial/vkorc1_warfarin_factorx_followup_planning.md` | Not cataloged | Keep as planning note until implemented | Resolved after sign-off: folded into `docs/roadmap.md` parking lot. |

`docs/tutorial/README.md` is the hub and should remain a meta page, not a
tutorial source unit.

## Cataloged Source Units

This table is the current source-of-truth alignment before grouping. `Generated`
means the same id appears in `docs/tutorial/manifest.json` and therefore has an
executable reference chapter.

| Catalog order | Tutorial id | Type | Status | Catalog path | Source JSON | Generated |
| ---: | --- | --- | --- | --- | --- | --- |
| 1 | `generated_hub` | executable_collection | generated+checked | `docs/tutorial/generated/README.md` | `docs/tutorial/sources/ref_generated_hub.json` | - |
| 2 | `vkorc1_warfarin_promoter_luciferase_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/08-04_vkorc1_warfarin_promoter_luciferase_gui.md` | `docs/tutorial/sources/08-04_vkorc1_warfarin_promoter_luciferase_gui.json` | - |
| 3 | `tp73_cdna_genomic_dotplot_gui` | gui_walkthrough | manual | `docs/tutorial/02-03_tp73_cdna_genomic_dotplot_gui.md` | `docs/tutorial/sources/02-03_tp73_cdna_genomic_dotplot_gui.json` | - |
| 4 | `gibson_specialist_testing_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/03-05_gibson_specialist_testing_gui.md` | `docs/tutorial/sources/03-05_gibson_specialist_testing_gui.json` | - |
| 5 | `gibson_arrangements_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/03-07_gibson_arrangements_gui.md` | `docs/tutorial/sources/03-07_gibson_arrangements_gui.json` | - |
| 6 | `gibson_physical_rack_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/03-08_gibson_physical_rack_gui.md` | `docs/tutorial/sources/03-08_gibson_physical_rack_gui.json` | - |
| 7 | `sequencing_confirmation_trace_cli` | cli_walkthrough | manual/hybrid | `docs/tutorial/10-01_sequencing_confirmation_trace_cli.md` | `docs/tutorial/sources/10-01_sequencing_confirmation_trace_cli.json` | - |
| 8 | `sequencing_confirmation_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/10-02_sequencing_confirmation_gui.md` | `docs/tutorial/sources/10-02_sequencing_confirmation_gui.json` | - |
| 9 | `rna_read_batch_gene_support_cli` | cli_walkthrough | manual/hybrid | `docs/tutorial/07-01_rna_read_batch_gene_support_cli.md` | `docs/tutorial/sources/07-01_rna_read_batch_gene_support_cli.json` | - |
| 10 | `tutorial_landscape_overview` | operational_reference | manual/reference | `docs/tutorial/landscape_overview.md` | `docs/tutorial/sources/ref_tutorial_landscape_overview.json` | - |
| 11 | `agent_interfaces` | operational_reference | manual/reference | `docs/tutorial/01-01_agent_interfaces.md` | `docs/tutorial/sources/01-01_agent_interfaces.json` | - |
| 12 | `metabion_external_service_handoff_gui_cli` | gui_cli_walkthrough | manual/hybrid | `docs/tutorial/09-01_metabion_external_service_handoff_gui_cli.md` | `docs/tutorial/sources/09-01_metabion_external_service_handoff_gui_cli.json` | - |
| 13 | `geneart_external_service_handoff_gui_cli` | gui_cli_walkthrough | manual/hybrid | `docs/tutorial/09-02_geneart_external_service_handoff_gui_cli.md` | `docs/tutorial/sources/09-02_geneart_external_service_handoff_gui_cli.json` | - |
| 14 | `simple_pcr_selection_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/04-01_simple_pcr_selection_gui.md` | `docs/tutorial/sources/04-01_simple_pcr_selection_gui.json` | `18` |
| 15 | `protein_transcript_native_expert_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/06-02_protein_transcript_native_expert_gui.md` | `docs/tutorial/sources/06-02_protein_transcript_native_expert_gui.json` | - |
| 16 | `protein_reverse_translation_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/06-01_protein_reverse_translation_gui.md` | `docs/tutorial/sources/06-01_protein_reverse_translation_gui.json` | - |
| 17 | `vkorc1_variant_followup_expert_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/08-02_vkorc1_variant_followup_expert_gui.md` | `docs/tutorial/sources/08-02_vkorc1_variant_followup_expert_gui.json` | - |
| 18 | `tp73_uniprot_projection_audit_cli` | cli_walkthrough | manual/hybrid | `docs/tutorial/06-05_tp73_uniprot_projection_audit_cli.md` | `docs/tutorial/sources/06-05_tp73_uniprot_projection_audit_cli.json` | `19` |
| 19 | `qpcr_exon_junctions_gui` | gui_walkthrough | manual/hybrid | `docs/tutorial/04-03_qpcr_exon_junctions_gui.md` | `docs/tutorial/sources/04-03_qpcr_exon_junctions_gui.json` | - |
| 20 | `load_branch_reverse_complement_pgex_fasta` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/02-01_load_branch_reverse_complement_pgex_fasta.md` | `docs/tutorial/sources/02-01_load_branch_reverse_complement_pgex_fasta.json` | `01` |
| 21 | `find_and_extend_genomic_target_local_catalog` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/05-01_find_and_extend_genomic_target_local_catalog.md` | `docs/tutorial/sources/05-01_find_and_extend_genomic_target_local_catalog.json` | `02` |
| 22 | `load_and_digest_pgex` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/03-01_load_and_digest_pgex.md` | `docs/tutorial/sources/03-01_load_and_digest_pgex.json` | `03` |
| 23 | `gibson_two_fragment_overlap_preview` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/03-03_gibson_two_fragment_overlap_preview.md` | `docs/tutorial/sources/03-03_gibson_two_fragment_overlap_preview.json` | `04` |
| 24 | `guides_filter_and_generate_oligos` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/04-04_guides_filter_and_generate_oligos.md` | `docs/tutorial/sources/04-04_guides_filter_and_generate_oligos.json` | `05` |
| 25 | `digest_ligation_extract_region_minimal` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/03-02_digest_ligation_extract_region_minimal.md` | `docs/tutorial/sources/03-02_digest_ligation_extract_region_minimal.json` | `06` |
| 26 | `guides_export_csv_and_protocol` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/04-05_guides_export_csv_and_protocol.md` | `docs/tutorial/sources/04-05_guides_export_csv_and_protocol.json` | `07` |
| 27 | `contribute_to_gentle_development` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/01-02_contribute_to_gentle_development.md` | `docs/tutorial/sources/01-02_contribute_to_gentle_development.json` | `08` |
| 28 | `prepare_reference_genome_online` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/05-02_prepare_reference_genome_online.md` | `docs/tutorial/sources/05-02_prepare_reference_genome_online.json` | `09` |
| 29 | `tp53_isoform_architecture_online` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/06-03_tp53_isoform_architecture_online.md` | `docs/tutorial/sources/06-03_tp53_isoform_architecture_online.json` | `10` |
| 30 | `tp63_anchor_extension_online` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/05-03_tp63_anchor_extension_online.md` | `docs/tutorial/sources/05-03_tp63_anchor_extension_online.json` | `11` |
| 31 | `tp53_multi_gene_sparse_mapping_online` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/07-02_tp53_multi_gene_sparse_mapping_online.md` | `docs/tutorial/sources/07-02_tp53_multi_gene_sparse_mapping_online.json` | `12` |
| 32 | `pcr_selection_batch_primer_pairs_offline` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/04-02_pcr_selection_batch_primer_pairs_offline.md` | `docs/tutorial/sources/04-02_pcr_selection_batch_primer_pairs_offline.json` | `13` |
| 33 | `tp73_cdna_genomic_dotplot_online` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/02-04_tp73_cdna_genomic_dotplot_online.md` | `docs/tutorial/sources/02-04_tp73_cdna_genomic_dotplot_online.json` | `14` |
| 34 | `gibson_specialist_testing_baseline` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/03-04_gibson_specialist_testing_baseline.md` | `docs/tutorial/sources/03-04_gibson_specialist_testing_baseline.json` | `15` |
| 35 | `gibson_arrangements_baseline` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/03-06_gibson_arrangements_baseline.md` | `docs/tutorial/sources/03-06_gibson_arrangements_baseline.json` | `16` |
| 36 | `tp53_uniprot_projection_online` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/06-04_tp53_uniprot_projection_online.md` | `docs/tutorial/sources/06-04_tp53_uniprot_projection_online.json` | `17` |
| 37 | `promoter_design_artifact_slice_offline` | executable_chapter | generated+checked/human-pending | `docs/tutorial/generated/chapters/08-03_promoter_design_artifact_slice_offline.md` | `docs/tutorial/sources/08-03_promoter_design_artifact_slice_offline.json` | `24` |
| 38 | `reporter_construct_handoff_cli` | cli_walkthrough | manual/hybrid | `docs/tutorial/08-05_reporter_construct_handoff_cli.md` | `docs/tutorial/sources/08-05_reporter_construct_handoff_cli.json` | - |

## Phase 0 Gate Decisions

Please sign off or adjust these before Phase 1:

1. Add `stateless_sequence_inspection_gui_cli` as an active source unit.
2. Add `tfbs_similarity_ranking_gui` as an active source unit.
3. Keep the two TP73 promoter/luciferase pages out of active decimal numbering,
   unless we deliberately introduce a deprecated/replaced-by tutorial state in
   the catalog view.
4. Keep `vkorc1_warfarin_factorx_followup_planning` out of active decimal
   numbering until it becomes an executable or manually validated tutorial.
5. Treat `generated_hub` and `tutorial_landscape_overview` as meta/reference
   units, not biological workflow tutorials.
6. Keep `agent_interfaces` as a tutorial unit even though its Markdown path is
   outside `docs/tutorial/`.

Once these are approved, Phase 1 can add group metadata and v4 schema support
without guessing what should be numbered.

## Phase 0 Sign-Off

Resolved on 2026-05-22 after external read-only review.

Accepted decisions:

1. Add `stateless_sequence_inspection_gui_cli` as an active guided source unit.
2. Add `tfbs_similarity_ranking_gui` as an active guided source unit.
3. Move `tp73_promoter_luciferase_gui.md` and
   `tp73_promoter_luciferase_showcase.md` out of active tutorial numbering;
   the preferred destination is `docs/tutorial/archive/`, with no catalog
   entry and no decimal id.
4. Move the planned
   `vkorc1_warfarin_factorx_followup_planning.md` content into the roadmap
   as a planned-tutorial note rather than cataloging it as a tutorial.
5. Keep `generated_hub` and `tutorial_landscape_overview` catalogued as
   reference/navigation units with no decimal id. This means catalog/manifest
   v2 must allow `decimal_id` to be absent/null.
6. Keep `agent_interfaces` as a tutorial unit, and move
   `docs/tutorial/01-01_agent_interfaces.md` into `docs/tutorial/` during the
   dedicated file-rename phase.
7. When the stale source file
   `docs/tutorial/sources/08-04_vkorc1_warfarin_promoter_luciferase_gui.json` is renamed,
   use its content-correct slug `vkorc1_warfarin_promoter_luciferase_gui`.

The audit also verified that `simple_pcr_selection_gui` and
`tp73_uniprot_projection_audit_cli` are both genuinely one source unit with two
surfaces: each has one source JSON containing both a catalog entry and a
generated chapter.

## Phase 1 Group Assignment Gate


Assigned on 2026-05-22. This is the Phase 1 review gate: all active tutorial source units are grouped, 38 tutorial units receive decimal ids, and the two reference/navigation units remain grouped but unnumbered.

| Unit id | Group | Position | decimal_id | Rationale |
| --- | --- | ---: | --- | --- |
| `generated_hub` | `01` | null | `null` | Reference/navigation unit; grouped with interfaces but intentionally unnumbered. |
| `tutorial_landscape_overview` | `01` | null | `null` | Reference/navigation unit; grouped with interfaces but intentionally unnumbered. |
| `agent_interfaces` | `01` | 1 | `01.01` | Interface tutorial rather than a biological workflow, so it starts the interfaces group. |
| `contribute_to_gentle_development` | `01` | 2 | `01.02` | Could be contributor/process material; grouped with interfaces because it teaches how humans and agents approach GENtle. |
| `load_branch_reverse_complement_pgex_fasta` | `02` | 1 | `02.01` | - |
| `stateless_sequence_inspection_gui_cli` | `02` | 2 | `02.02` | Could fit restriction/TFBS topics; grouped here because the primary lesson is state-optional sequence inspection. |
| `tp73_cdna_genomic_dotplot_gui` | `02` | 3 | `02.03` | Could fit genome or transcript topics; grouped here because the primary lesson is comparing two sequence representations. |
| `tp73_cdna_genomic_dotplot_online` | `02` | 4 | `02.04` | Could fit genome or transcript topics; grouped here because the primary lesson is comparing cDNA and genomic sequence representations. |
| `load_and_digest_pgex` | `03` | 1 | `03.01` | - |
| `digest_ligation_extract_region_minimal` | `03` | 2 | `03.02` | - |
| `gibson_two_fragment_overlap_preview` | `03` | 3 | `03.03` | - |
| `gibson_specialist_testing_baseline` | `03` | 4 | `03.04` | - |
| `gibson_specialist_testing_gui` | `03` | 5 | `03.05` | - |
| `gibson_arrangements_baseline` | `03` | 6 | `03.06` | - |
| `gibson_arrangements_gui` | `03` | 7 | `03.07` | - |
| `gibson_physical_rack_gui` | `03` | 8 | `03.08` | - |
| `simple_pcr_selection_gui` | `04` | 1 | `04.01` | - |
| `pcr_selection_batch_primer_pairs_offline` | `04` | 2 | `04.02` | - |
| `qpcr_exon_junctions_gui` | `04` | 3 | `04.03` | Could fit splicing/expression; grouped here because qPCR assay design is the primary learning intent. |
| `guides_filter_and_generate_oligos` | `04` | 4 | `04.04` | Could fit cloning guides; grouped here because the user-facing product is filtered oligos. |
| `guides_export_csv_and_protocol` | `04` | 5 | `04.05` | Could fit external handoff; grouped here because it extends the oligo-generation workflow. |
| `find_and_extend_genomic_target_local_catalog` | `05` | 1 | `05.01` | - |
| `prepare_reference_genome_online` | `05` | 2 | `05.02` | Could fit external services; grouped here because it prepares genome-coordinate context for later work. |
| `tp63_anchor_extension_online` | `05` | 3 | `05.03` | - |
| `protein_reverse_translation_gui` | `06` | 1 | `06.01` | - |
| `protein_transcript_native_expert_gui` | `06` | 2 | `06.02` | - |
| `tp53_isoform_architecture_online` | `06` | 3 | `06.03` | - |
| `tp53_uniprot_projection_online` | `06` | 4 | `06.04` | Could fit external services; grouped here because UniProt is used to teach protein-to-transcript projection. |
| `tp73_uniprot_projection_audit_cli` | `06` | 5 | `06.05` | Could fit external services; grouped here because the audit teaches protein feature projection and coding-sequence attribution. |
| `rna_read_batch_gene_support_cli` | `07` | 1 | `07.01` | - |
| `tp53_multi_gene_sparse_mapping_online` | `07` | 2 | `07.02` | Could fit genome context; grouped here because the primary output is RNA-read support across gene cohorts. |
| `tfbs_similarity_ranking_gui` | `08` | 1 | `08.01` | - |
| `vkorc1_variant_followup_expert_gui` | `08` | 2 | `08.02` | Could fit genome/variant context; grouped here because variant interpretation feeds regulatory reporter follow-up. |
| `promoter_design_artifact_slice_offline` | `08` | 3 | `08.03` | - |
| `vkorc1_warfarin_promoter_luciferase_gui` | `08` | 4 | `08.04` | Could fit cloning/assembly; grouped here because the biological endpoint is promoter-reporter design. |
| `reporter_construct_handoff_cli` | `08` | 5 | `08.05` | Could fit external handoff or cloning; grouped here because it closes a regulatory reporter-design workflow. |
| `metabion_external_service_handoff_gui_cli` | `09` | 1 | `09.01` | Could fit oligo/PCR workflows; grouped here because the primary lesson is service handoff. |
| `geneart_external_service_handoff_gui_cli` | `09` | 2 | `09.02` | Could fit cloning/protein workflows; grouped here because the primary lesson is external service handoff. |
| `sequencing_confirmation_trace_cli` | `10` | 1 | `10.01` | - |
| `sequencing_confirmation_gui` | `10` | 2 | `10.02` | - |
