---
chapter_id: "contribute_to_gentle_development"
title: "Contribute to GENtle development"
tier: "advanced"
example_id: "contribute_gentle_development_baseline"
source_example: "docs/examples/workflows/contribute_gentle_development_baseline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "unreviewed"
review_stale: false
codex_reviewed_at: null
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: "Tutorial confusion"
review_issue_template_path: ".github/ISSUE_TEMPLATE/tutorial-confusion.md"
generated_artifact_dir: "docs/tutorial/generated/artifacts/contribute_to_gentle_development"
---

# Contribute to GENtle development

Executable contributor onboarding from baseline run to test/documentation checks.

Contributing effectively requires keeping biological behavior, command contracts, and docs in sync. This chapter provides a concrete contributor routine: run a baseline operation chain, then validate examples, tutorial generation, and tests before opening a change request.

## What You Will Accomplish

- Understand the expected contributor loop from change to verification.
- Use tutorial and example checks to prevent documentation drift.
- Identify the minimum validation commands expected before handoff or PR.

## Before You Start

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

**Useful when:**

- You want to submit a feature or fix without breaking cross-interface behavior.
- You need to verify that tutorial/docs stay synchronized with executable workflows.
- You want a deterministic pre-PR checklist for local validation.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Run the baseline sequence routine in the GUI and inspect resulting lineage entries

**GUI**

Run the baseline sequence routine in the GUI and inspect resulting lineage entries.

**CLI (terminal)**

```bash
gentle_cli workflow @docs/examples/workflows/contribute_gentle_development_baseline.json
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Run the baseline sequence routine in the GUI and inspect resulting lineage entries. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The baseline workflow creates `contrib_seed`, `contrib_branch`, and `contrib_branch_rc` through the shared engine.

### Step 2: Locate the same routine as a canonical workflow JSON example in docs/examples/workflows

**GUI**

Locate the same routine as a canonical workflow JSON example in `docs/examples/workflows`.

**CLI (terminal)**

```bash
cargo run --bin gentle_examples_docs -- --check
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Locate the same routine as a canonical workflow JSON example in `docs/examples/workflows`. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Example generation/checking confirms the canonical workflow JSON remains parseable and adapter snippets are current.

### Step 3: Use this mapping to validate that your planned code change affects shared engine behavior, not only one interface

**GUI**

Use this mapping to validate that your planned code change affects shared engine behavior, not only one interface.

**CLI (terminal)**

```bash
cargo check -q
cargo run --bin gentle_examples_docs -- tutorial-check
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Use this mapping to validate that your planned code change affects shared engine behavior, not only one interface. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The Rust check and tutorial drift check catch interface or documentation changes that escaped the GUI path.


## Interpretation and Reference

## Parameters That Matter

- `Workflow file path @docs/examples/workflows/...` (where used: command equivalent)
  - Why it matters: Ensures you are validating the same canonical routine used by docs/tests.
  - How to derive it: Select the chapter's `example_id` and open the matching JSON file under `docs/examples/workflows`.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Tutorial Drift Checks** (`tutorial_drift_checks`): Tutorial content is generated from executable examples and verified in automated checks.
- **Contribution Loop** (`contribution_loop`): Contributions should couple code edits with docs updates and deterministic tests.

## Follow-up Commands

```bash
cargo check -q
cargo test -q workflow_examples -- --test-threads=1
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
gentle_cli capabilities
```

## Checkpoints

- Baseline contribution workflow executes without failures.
- Contributor validation commands pass in a clean environment.
- Tutorial-generated content remains synchronized with source manifests/examples.

## Tutorial Provenance

- Chapter id: `contribute_to_gentle_development`
- Tier: `advanced`
- Example id: `contribute_gentle_development_baseline`
- Tutorial source JSON: `docs/tutorial/sources/01-02_contribute_to_gentle_development.json`
- Workflow file: `docs/examples/workflows/contribute_gentle_development_baseline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/contribute_to_gentle_development`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Contribute to GENtle development`
- Tutorial/chapter id: `contribute_to_gentle_development`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
