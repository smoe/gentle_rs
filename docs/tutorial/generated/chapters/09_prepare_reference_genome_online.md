---
chapter_id: "prepare_reference_genome_online"
title: "Prepare a reference genome cache (online)"
tier: "online"
example_id: "prepare_reference_genome_online"
source_example: "docs/examples/workflows/prepare_reference_genome_online.json"
example_test_mode: "online"
executed_during_generation: false
---

# Prepare a reference genome cache (online)

Document network-dependent preparation without destabilizing offline default checks.

Reference-genome preparation is crucial for genome-anchored cloning interpretation, but it depends on online resources and local tool setup. This chapter keeps that path explicit and opt-in so default tutorial checks remain robust.

## When This Routine Is Useful

- You need genome-anchored extraction around gene/promoter context.
- You want to prepare local cache/index assets for repeated anchor operations.
- You need to understand which routines are intentionally online-only in CI defaults.

## What You Learn

- Recognize which tutorial flows require network access and external tools.
- Use `GENTLE_TEST_ONLINE` to opt into online chapter execution.
- Preserve offline CI reliability while still documenting online capabilities.

## Concepts

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.

## GUI First

1. Open prepared-reference controls from the GUI menus.
2. Select the target genome and start preparation with explicit cache settings.
3. Confirm prepared status in the GUI before attempting extraction workflows.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/prepare_reference_genome_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/prepare_reference_genome_online.json'
```

## Parameters That Matter

- `PrepareGenome.genome_id` (where used: operation 1)
  - Why it matters: Selects the exact reference build and annotation set.
  - How to derive it: Choose the genome build matching your experimental system and downstream coordinate system.
- `PrepareGenome.catalog_path / cache_dir` (where used: operation 1)
  - Why it matters: Controls source catalog and local cache destination.
  - How to derive it: Use repository defaults unless your environment requires custom catalogs or cache locations.

## Follow-up Commands

```bash
GENTLE_TEST_ONLINE=1 cargo run --bin gentle_examples_docs -- tutorial-generate
```

## Checkpoints

- Genome preparation runs only when GENTLE_TEST_ONLINE is enabled.
- Offline generation still emits the chapter with execution status noted.

## Canonical Source

- Chapter id: `prepare_reference_genome_online`
- Tier: `online`
- Example id: `prepare_reference_genome_online`
- Workflow file: `docs/examples/workflows/prepare_reference_genome_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect this JSON file directly when you need full option-level detail.
