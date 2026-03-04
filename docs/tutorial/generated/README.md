# GENtle Tutorial (Generated)

This folder is generated from:
- `docs/tutorial/manifest.json`
- `docs/examples/workflows/*.json`

Regenerate with:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
```

Validate committed generated output:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
```

Recommended progression:
- Start chapters in the GUI to understand biological intent and visual context.
- Then use the command equivalents for repeatable runs and automation.

Online execution was disabled (`GENTLE_TEST_ONLINE=0` during generation).

## Chapters

- 1. [Load FASTA, branch, and reverse-complement](./chapters/01_load_branch_reverse_complement_pgex_fasta.md) - `core` - example `load_branch_reverse_complement_pgex_fasta` - executed `yes`
- 2. [Find and extend the right genomic target (local catalog)](./chapters/02_find_and_extend_genomic_target_local_catalog.md) - `core` - example `prepare_extract_extend_localproject_gene` - executed `yes`
- 3. [Load pGEX and digest with BamHI/EcoRI](./chapters/03_load_and_digest_pgex.md) - `core` - example `load_and_digest_pgex` - executed `yes`
- 4. [Guide practical filtering and oligo generation](./chapters/04_guides_filter_and_generate_oligos.md) - `core` - example `guides_filter_and_generate_oligos` - executed `yes`
- 5. [Digest -> Ligation -> ExtractRegion minimal slice](./chapters/05_digest_ligation_extract_region_minimal.md) - `core` - example `digest_ligation_extract_region_minimal` - executed `yes`
- 6. [Guide oligo export (CSV + protocol)](./chapters/06_guides_export_csv_and_protocol.md) - `advanced` - example `guides_export_csv_and_protocol` - executed `yes`
- 7. [Contribute to GENtle development](./chapters/07_contribute_to_gentle_development.md) - `advanced` - example `contribute_gentle_development_baseline` - executed `yes`
- 8. [Prepare a reference genome cache (online)](./chapters/08_prepare_reference_genome_online.md) - `online` - example `prepare_reference_genome_online` - executed `no`
- 9. [TP53 isoform architecture expert panel (online)](./chapters/09_tp53_isoform_architecture_online.md) - `online` - example `tp53_isoform_architecture_online` - executed `no`
- 10. [Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./chapters/10_tp63_anchor_extension_online.md) - `online` - example `tp63_extend_anchor_online` - executed `no`

## Concepts and Where They Recur

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Appears in: [Chapter 1: Load FASTA, branch, and reverse-complement](./chapters/01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./chapters/02_find_and_extend_genomic_target_local_catalog.md), [Chapter 7: Contribute to GENtle development](./chapters/07_contribute_to_gentle_development.md), [Chapter 9: TP53 isoform architecture expert panel (online)](./chapters/09_tp53_isoform_architecture_online.md), [Chapter 10: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./chapters/10_tp63_anchor_extension_online.md).
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
  - Appears in: [Chapter 1: Load FASTA, branch, and reverse-complement](./chapters/01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./chapters/02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./chapters/03_load_and_digest_pgex.md), [Chapter 4: Guide practical filtering and oligo generation](./chapters/04_guides_filter_and_generate_oligos.md), [Chapter 5: Digest -> Ligation -> ExtractRegion minimal slice](./chapters/05_digest_ligation_extract_region_minimal.md), [Chapter 8: Prepare a reference genome cache (online)](./chapters/08_prepare_reference_genome_online.md).
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
  - Appears in: [Chapter 1: Load FASTA, branch, and reverse-complement](./chapters/01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./chapters/02_find_and_extend_genomic_target_local_catalog.md), [Chapter 3: Load pGEX and digest with BamHI/EcoRI](./chapters/03_load_and_digest_pgex.md), [Chapter 5: Digest -> Ligation -> ExtractRegion minimal slice](./chapters/05_digest_ligation_extract_region_minimal.md), [Chapter 10: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./chapters/10_tp63_anchor_extension_online.md).
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
  - Appears in: [Chapter 2: Find and extend the right genomic target (local catalog)](./chapters/02_find_and_extend_genomic_target_local_catalog.md), [Chapter 8: Prepare a reference genome cache (online)](./chapters/08_prepare_reference_genome_online.md), [Chapter 9: TP53 isoform architecture expert panel (online)](./chapters/09_tp53_isoform_architecture_online.md), [Chapter 10: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./chapters/10_tp63_anchor_extension_online.md).
- **Guide Design Pipeline** (`guide_design_pipeline`): Guide sets can be created, filtered, expanded to oligos, and exported with protocol context.
  - Appears in: [Chapter 4: Guide practical filtering and oligo generation](./chapters/04_guides_filter_and_generate_oligos.md), [Chapter 6: Guide oligo export (CSV + protocol)](./chapters/06_guides_export_csv_and_protocol.md).
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
  - Appears in: [Chapter 6: Guide oligo export (CSV + protocol)](./chapters/06_guides_export_csv_and_protocol.md), [Chapter 9: TP53 isoform architecture expert panel (online)](./chapters/09_tp53_isoform_architecture_online.md).
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.
  - Appears in: [Chapter 7: Contribute to GENtle development](./chapters/07_contribute_to_gentle_development.md).
- **Contribution Loop** (`contribution_loop`): Contributions should couple code edits with docs updates and deterministic tests.
  - Appears in: [Chapter 7: Contribute to GENtle development](./chapters/07_contribute_to_gentle_development.md).
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
  - Appears in: [Chapter 8: Prepare a reference genome cache (online)](./chapters/08_prepare_reference_genome_online.md), [Chapter 9: TP53 isoform architecture expert panel (online)](./chapters/09_tp53_isoform_architecture_online.md), [Chapter 10: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./chapters/10_tp63_anchor_extension_online.md).
- **Isoform Architecture Panels** (`isoform_architecture_panels`): Curated transcript/protein architecture overlays can be imported and rendered as deterministic expert-view SVG outputs.
  - Appears in: [Chapter 9: TP53 isoform architecture expert panel (online)](./chapters/09_tp53_isoform_architecture_online.md).

## Source Summary

- Tutorial schema: `gentle.tutorial_manifest.v1`
- Chapter count: `10`
- Generation report: [`report.json`](./report.json)
