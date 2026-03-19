# TP53 isoform architecture expert panel (online)

- Chapter id: `tp53_isoform_architecture_online`
- Tier: `online`
- Example id: `tp53_isoform_architecture_online`
- Source example: `docs/examples/workflows/tp53_isoform_architecture_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.

Build a TP53 locus project and export Figure-1-style transcript/protein isoform architecture from one deterministic engine route.

This chapter demonstrates a publication-oriented use case: derive TP53 from a prepared GRCh38 reference, import a curated isoform panel, and render transcript/protein architecture from the same expert-view payload used by GUI and shell interfaces. The focus is parity and provenance, not manual figure drawing.

## When This Routine Is Useful

- You want transcript and protein isoform architecture from one sequence context with explicit panel provenance.
- You need a deterministic SVG baseline for Figure-1-style TP53 isoform presentation.
- You want to preserve the exact panel import and rendering operations in sequence history.

## What You Learn

- Execute an end-to-end isoform-panel workflow using shared engine operations.
- Understand how curated panel JSON joins genome-derived transcripts without GUI-only logic.
- Export deterministic architecture SVG suitable for downstream figure styling.

## Concepts and Recurrence

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
  - Status: reinforced from [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md), [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md), [Chapter 8: Contribute to GENtle development](./08_contribute_to_gentle_development.md).
  - Reoccurs in: [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./13_tp73_cdna_genomic_dotplot_online.md).
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
  - Status: reinforced from [Chapter 2: Find and extend the right genomic target (local catalog)](./02_find_and_extend_genomic_target_local_catalog.md), [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md).
  - Reoccurs in: [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./13_tp73_cdna_genomic_dotplot_online.md).
- **Isoform Architecture Panels** (`isoform_architecture_panels`): Curated transcript/protein architecture overlays can be imported and rendered as deterministic expert-view SVG outputs.
  - Status: introduced in this chapter.
  - Reoccurs in: no later chapter.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
  - Status: reinforced from [Chapter 9: Prepare a reference genome cache (online)](./09_prepare_reference_genome_online.md).
  - Reoccurs in: [Chapter 11: Retrieve TP63 and extend the displayed region by +/-2 kb (online)](./11_tp63_anchor_extension_online.md), [Chapter 12: Map TP53 locus reads with multi-gene sparse indexing (online)](./12_tp53_multi_gene_sparse_mapping_online.md), [Chapter 13: Compare TP73 cDNA against TP73 genomic context via dotplot (online)](./13_tp73_cdna_genomic_dotplot_online.md).
- **Artifact Exports** (`artifact_exports`): Representative outputs (CSV/protocol/SVG/text) are retained for auditability and sharing.
  - Status: reinforced from [Chapter 7: Guide oligo export (CSV + protocol)](./07_guides_export_csv_and_protocol.md).
  - Reoccurs in: no later chapter.

## GUI First

1. Prepare `Human GRCh38 Ensembl 116` and extract gene `TP53` into `grch38_tp53`.
2. Open DNA window `Engine Ops -> Isoform architecture panels`, import `assets/panels/tp53_isoforms_v1.json`.
3. Open `Isoform Expert` and export SVG from the same panel context.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp53_isoform_architecture_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp53_isoform_architecture_online.json'
```

## Parameters That Matter

- `ImportIsoformPanel.panel_path / panel_id / strict` (where used: operation 3)
  - Why it matters: Defines which curated panel is loaded and whether transcript mapping mismatches should fail hard.
  - How to derive it: Use `assets/panels/tp53_isoforms_v1.json` and keep `strict=false` for exploratory mapping; switch to `strict=true` when curation is locked.
- `RenderIsoformArchitectureSvg.path` (where used: operation 4)
  - Why it matters: Controls deterministic export location used for tutorial artifact retention and figure review.
  - How to derive it: Use a stable project-relative path such as `exports/tp53_isoform_architecture.svg`.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'panels inspect-isoform grch38_tp53 tp53_isoforms_v1'
cargo run --bin gentle_cli -- save-project tp53_isoform_architecture.project.gentle.json
```

## Checkpoints

- Panel import reports mapped transcript lanes and any unresolved transcript warnings.
- Isoform architecture SVG export succeeds with deterministic lane ordering.
- Saved project contains TP53 sequence plus imported panel metadata for replay.

## Retained Outputs

- None for this chapter.

## Canonical Source

- Workflow file: `docs/examples/workflows/tp53_isoform_architecture_online.json`
- Inspect this JSON file directly when you need full option-level detail.
