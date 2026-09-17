---
chapter_id: "gibson_two_fragment_overlap_preview"
title: "Gibson two-fragment overlap planning baseline"
tier: "core"
example_id: "gibson_two_fragment_overlap_preview"
source_example: "docs/examples/workflows/gibson_two_fragment_overlap_preview.json"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/gibson_two_fragment_overlap_preview"
---

# Gibson two-fragment overlap planning baseline

Use the built-in Gibson routine baseline to validate overlap assumptions and generate deterministic preview assemblies.

Gibson assembly in practice depends on overlap design quality. This chapter anchors that design logic in GENtle through one repeatable workflow: build two overlapping fragments, run Gibson-specific preflight checks on overlap compatibility, and produce forward-order preview sequences for review and communication.

## What You Will Accomplish

- Understand where Gibson overlap assumptions are checked in shared-shell preflight.
- Interpret overlap mismatch diagnostics before state mutation.
- Use deterministic output IDs for clear cloning communication.

## Before You Start

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

**Useful when:**

- You want reproducible documentation of intended Gibson fragment order and overlap assumptions.
- You need fast preflight feedback before primer ordering or wet-lab execution.
- You want one routine that can be explained to collaborators through GUI + shell parity.

## Walkthrough: GUI, CLI and Inner Agent

The three routes below describe the same operation. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Create two overlapping pGEX fragments (gibson_left

**GUI**

Create two overlapping pGEX fragments (`gibson_left`, `gibson_right`) from `test_files/pGEX_3X.fa` using region extraction.

**CLI / GUI Shell**

```bash
gentle_cli workflow @docs/examples/workflows/gibson_two_fragment_overlap_preview.json
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Create two overlapping pGEX fragments (`gibson_left`, `gibson_right`) from `test_files/pGEX_3X.fa` using region extraction. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The workflow loads `pgex_fasta`, extracts `gibson_left` and `gibson_right`, and leaves deterministic preview IDs for inspection.

### Step 2: Import the Gibson two-fragment overlap preview routine from Patterns ->

**GUI**

Import the Gibson two-fragment overlap preview routine from `Patterns -> Routine catalog`.

**CLI / GUI Shell**

```bash
gentle_cli shell 'macros template-import assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json'
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Import the Gibson two-fragment overlap preview routine from `Patterns -> Routine catalog`. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The routine catalog registers `gibson_two_fragment_overlap_preview`, so the next shell command can validate that specific fragment order.

### Step 3: Validate the Gibson overlap preview from Shell

**GUI**

Validate the Gibson overlap preview from `Shell`, then rerun the same template without `--validate-only` to create outputs.

**CLI / GUI Shell**

```bash
gentle_cli shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Validate the Gibson overlap preview from `Shell`, then rerun the same template without `--validate-only` to create outputs. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The validate-only run reports executable bindings for the chosen overlap; rerun with `--transactional` to create the named preview assembly.


## Interpretation and Reference

## Parameters That Matter

- `left_seq_id / right_seq_id` (where used: template preflight + run)
  - Why it matters: Order defines the expected assembled product direction and adjacent overlap checks.
  - How to derive it: Choose fragment IDs in intended 5'->3' assembly order.
- `overlap_bp` (where used: Gibson family preflight check)
  - Why it matters: Controls suffix/prefix overlap length validation between adjacent fragments.
  - How to derive it: Use planned primer-homology length (commonly 15-40 bp for many Gibson workflows).
- `assembly_prefix / output_id` (where used: template run output naming)
  - Why it matters: Stable IDs make review/export/communication unambiguous.
  - How to derive it: Use project-specific names (for example `tp73_gibson_round1`).
  - Omit when: Omit only when default IDs are acceptable for exploratory runs.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## Follow-up Commands

```bash
gentle_cli shell 'macros template-import assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json'
gentle_cli shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
gentle_cli shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --transactional'
```

## Checkpoints

- Validate-only run returns `can_execute=true` for matching overlap fragments.
- Mismatch overlap inputs return explicit Gibson overlap diagnostics.
- Mutating run creates deterministic preview outputs (`${assembly_prefix}_*` and `${output_id}`).

## Tutorial Provenance

- Chapter id: `gibson_two_fragment_overlap_preview`
- Tier: `core`
- Example id: `gibson_two_fragment_overlap_preview`
- Tutorial source JSON: `docs/tutorial/sources/03-03_gibson_two_fragment_overlap_preview.json`
- Workflow file: `docs/examples/workflows/gibson_two_fragment_overlap_preview.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/gibson_two_fragment_overlap_preview`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Gibson two-fragment overlap planning baseline`
- Tutorial/chapter id: `gibson_two_fragment_overlap_preview`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
