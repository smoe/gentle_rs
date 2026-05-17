---
chapter_id: "tp73_cdna_genomic_dotplot_online"
title: "Compare TP73 cDNA against TP73 genomic context via dotplot (online)"
tier: "online"
example_id: "tp73_cdna_genomic_dotplot_online"
source_example: "docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json"
example_test_mode: "online"
executed_during_generation: false
---

# Compare TP73 cDNA against TP73 genomic context via dotplot (online)

Fetch TP73 cDNA, extract TP73 from GRCh38, and compute a pair-forward dotplot with high-sensitivity defaults (`word=7`, `step=1`, `mismatches=0`).

This chapter captures a practical cDNA-vs-genomic verification route for transcript structure interpretation. The focus is a reproducible first-pass map that reveals exon-aligned block patterns while preserving one shared operation path across GUI and CLI interfaces.

> **How to Run This Locally**
> Set `GENTLE_TEST_ONLINE=1` and run from the repository root. This chapter fetches TP73 cDNA `NM_001126241.3` from NCBI/GenBank and prepares GRCh38 Ensembl 116 from Ensembl FTP before computing the local dotplot.

## When This Routine Is Useful

- You want a deterministic cDNA-vs-genomic control for exon/intron-aware interpretation.
- You need a reproducible TP73 baseline for demonstrating dotplot settings to collaborators.
- You want one canonical workflow file that mirrors your GUI tutorial run.

## What You Learn

- Use one deterministic route from accession retrieval to pairwise dotplot artifact generation.
- Understand why high-sensitivity seed sampling (`step=1`, low word size) improves cDNA-vs-genomic block visibility.
- Trace parity between GUI actions and canonical workflow JSON execution.

## Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.

## GUI First

1. Fetch GenBank accession `NM_001126241.3` as `tp73_cdna`.
2. Prepare `Human GRCh38 Ensembl 116`, retrieve gene `TP73` as `tp73_genomic`.
3. Open `tp73_cdna`, switch to `Dotplot map`, set pair mode against `tp73_genomic`, and compute with `word<=7`, `step=1`, `mismatches=0`.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json'
```

## Parameters That Matter

- `ComputeDotplot.word_size / step_bp / max_mismatches` (where used: operation 4)
  - Why it matters: These settings control sensitivity for dispersed exon-aligned signal in long genomic spans.
  - How to derive it: Start with `word=7`, `step=1`, `max_mismatches=0`; reduce to `word=6` only if additional sensitivity is needed.
- `ComputeDotplot.reference_seq_id` (where used: operation 4)
  - Why it matters: Pairwise mode requires explicit reference identity for y-axis mapping.
  - How to derive it: Use the exact ID emitted by the gene-extraction step (`tp73_genomic` in this chapter).

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'dotplot show tp73_cdna_vs_genomic_dotplot'
```

## Checkpoints

- Dotplot artifact is created and can be listed/shown by id (`tp73_cdna_vs_genomic_dotplot`).
- Pair-forward map shows exon-block style structure rather than one continuous full-length diagonal.
- GUI tutorial and workflow JSON describe the same parameter baseline.

## Canonical Source

- Chapter id: `tp73_cdna_genomic_dotplot_online`
- Tier: `online`
- Example id: `tp73_cdna_genomic_dotplot_online`
- Workflow file: `docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect this JSON file directly when you need full option-level detail.
