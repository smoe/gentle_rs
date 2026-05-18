---
chapter_id: "tp53_multi_gene_sparse_mapping_online"
title: "Map TP53 locus reads with multi-gene sparse indexing (online)"
tier: "online"
example_id: "tp53_multi_gene_sparse_mapping_online"
source_example: "docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json"
example_test_mode: "online"
executed_during_generation: false
---

# Map TP53 locus reads with multi-gene sparse indexing (online)

Prepare GRCh38, extract TP53, then map RNA reads with origin_mode=multi_gene_sparse against TP53-family targets.

This chapter extends the TP53 genome-targeting path toward read-origin mapping in one deterministic route. After extracting the TP53 locus from GRCh38, you run `InterpretRnaReads` with `multi_gene_sparse` so the index can include additional target genes (for example TP53/TP63/TP73) in one run. The goal is to keep genome anchoring and multi-gene interpretation in the same reproducible workflow contract.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. The workflow prepares/extracts GRCh38 Ensembl 116 from Ensembl FTP, then runs the multi-gene sparse RNA-read interpretation against locally derived TP53-family transcript templates.

## Parameters That Matter

- `InterpretRnaReads.seed_feature_id` (where used: operation 3)
  - Why it matters: Defines which feature seeds the splicing view and baseline transcript scope for sparse expansion.
  - How to derive it: In GUI this is implicit from the selected transcript. For direct workflow/CLI editing, inspect TP53 feature indices first and set the correct mRNA feature id.
- `InterpretRnaReads.origin_mode / target_gene_ids` (where used: operation 3)
  - Why it matters: Controls whether indexing stays baseline (`single_gene`) or is expanded with matched target-gene transcripts (`multi_gene_sparse`).
  - How to derive it: Use TP53-family IDs (`TP53, TP63, TP73`) for contrast runs, then narrow/expand based on biological question.
- `InterpretRnaReads.roi_seed_capture_enabled` (where used: operation 3)
  - Why it matters: Tracks request intent for future annotation-independent ROI capture layer.
  - How to derive it: Keep `false` for current runtime behavior; set `true` only when you want the deterministic pending-feature warning in provenance.

## When This Routine Is Useful

- You want one TP53-based run that can contrast seed support across TP53-family targets.
- You need a reproducible baseline for comparing single-gene vs multi-gene sparse indexing behavior.
- You want GUI and CLI routes to produce the same multi-gene request payload.

## What You Learn

- Map one gene locus from reference genome preparation through read-interpretation in a single operation chain.
- Understand what `multi_gene_sparse` changes at runtime (expanded local transcript-template indexing).
- Interpret deterministic report provenance for `origin_mode`, target-gene set, and planned ROI capture flag.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.

## GUI First

1. Prepare `Human GRCh38 Ensembl 116` and extract `TP53` into `grch38_tp53`.
2. Open Splicing Expert for a TP53 transcript, set `Origin mode` to `multi_gene_sparse`, and set `Target genes` to `TP53, TP63, TP73`.
3. Run Nanopore interpretation from the same panel and inspect sparse-origin warnings/summary fields in the report section.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json'
```

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json
cargo run --bin gentle_cli -- shell 'rna-reads list-reports grch38_tp53'
cargo run --bin gentle_cli -- shell 'rna-reads show-report tp53_family_sparse_template'
```

## Checkpoints

- TP53 locus extraction and RNA-read interpretation operations are both present in one workflow run history.
- Report summary shows `origin_mode=multi_gene_sparse` and target-gene provenance fields.
- Warnings remain explicit for any missing target genes and for ROI-capture requests when enabled.

## Canonical Source

- Chapter id: `tp53_multi_gene_sparse_mapping_online`
- Tier: `online`
- Example id: `tp53_multi_gene_sparse_mapping_online`
- Workflow file: `docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect this JSON file directly when you need full option-level detail.
