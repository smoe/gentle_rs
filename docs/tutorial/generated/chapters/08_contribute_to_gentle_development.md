---
chapter_id: "contribute_to_gentle_development"
title: "Contribute to GENtle development"
tier: "advanced"
example_id: "contribute_gentle_development_baseline"
source_example: "docs/examples/workflows/contribute_gentle_development_baseline.json"
example_test_mode: "always"
executed_during_generation: true
---

# Contribute to GENtle development

Executable contributor onboarding from baseline run to test/documentation checks.

Contributing effectively requires keeping biological behavior, command contracts, and docs in sync. This chapter provides a concrete contributor routine: run a baseline operation chain, then validate examples, tutorial generation, and tests before opening a change request.

## When This Routine Is Useful

- You want to submit a feature or fix without breaking cross-interface behavior.
- You need to verify that tutorial/docs stay synchronized with executable workflows.
- You want a deterministic pre-PR checklist for local validation.

## What You Learn

- Understand the expected contributor loop from change to verification.
- Use tutorial and example checks to prevent documentation drift.
- Identify the minimum validation commands expected before handoff or PR.

## Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.
- **Contribution Loop** (`contribution_loop`): Contributions should couple code edits with docs updates and deterministic tests.

## GUI First

1. Run the baseline sequence routine in the GUI and inspect resulting lineage entries.
2. Locate the same routine as a canonical workflow JSON example in `docs/examples/workflows`.
3. Use this mapping to validate that your planned code change affects shared engine behavior, not only one interface.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/contribute_gentle_development_baseline.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/contribute_gentle_development_baseline.json'
```

## Parameters That Matter

- `Workflow file path @docs/examples/workflows/...` (where used: command equivalent)
  - Why it matters: Ensures you are validating the same canonical routine used by docs/tests.
  - How to derive it: Select the chapter's `example_id` and open the matching JSON file under `docs/examples/workflows`.

## Follow-up Commands

```bash
cargo check -q
cargo test -q workflow_examples -- --test-threads=1
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_cli -- capabilities
```

## Checkpoints

- Baseline contribution workflow executes without failures.
- Contributor validation commands pass in a clean environment.
- Tutorial-generated content remains synchronized with source manifests/examples.

## Canonical Source

- Chapter id: `contribute_to_gentle_development`
- Tier: `advanced`
- Example id: `contribute_gentle_development_baseline`
- Workflow file: `docs/examples/workflows/contribute_gentle_development_baseline.json`
- Example test_mode: `always`
- Executed during generation: `yes`
- Inspect this JSON file directly when you need full option-level detail.
