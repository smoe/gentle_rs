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
codex_reviewed_at: null
human_reviewed_at: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/gibson_two_fragment_overlap_preview"
---

# Gibson two-fragment overlap planning baseline

Use the built-in Gibson routine baseline to validate overlap assumptions and generate deterministic preview assemblies.

Gibson assembly in practice depends on overlap design quality. This chapter anchors that design logic in GENtle through one repeatable workflow: build two overlapping fragments, run Gibson-specific preflight checks on overlap compatibility, and produce forward-order preview sequences for review and communication.

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

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

## When This Routine Is Useful

- You want reproducible documentation of intended Gibson fragment order and overlap assumptions.
- You need fast preflight feedback before primer ordering or wet-lab execution.
- You want one routine that can be explained to collaborators through GUI + shell parity.

## What You Learn

- Understand where Gibson overlap assumptions are checked in shared-shell preflight.
- Interpret overlap mismatch diagnostics before state mutation.
- Use deterministic output IDs for clear cloning communication.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Create two overlapping pGEX fragments (gibson_left, gibson_right) from test_f...

GUI: Create two overlapping pGEX fragments (`gibson_left`, `gibson_right`) from `test_files/pGEX_3X.fa` using region extraction.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_two_fragment_overlap_preview.json
```

> Expected: The workflow loads `pgex_fasta`, extracts `gibson_left` and `gibson_right`, and leaves deterministic preview IDs for inspection.

### Step 2: Import the Gibson two-fragment overlap preview routine from Patterns -> Routi...

GUI: Import the Gibson two-fragment overlap preview routine from `Patterns -> Routine catalog`.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json'
```

> Expected: The routine catalog registers `gibson_two_fragment_overlap_preview`, so the next shell command can validate that specific fragment order.

### Step 3: Validate the Gibson overlap preview from Shell, then rerun the same template ...

GUI: Validate the Gibson overlap preview from `Shell`, then rerun the same template without `--validate-only` to create outputs.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
```

> Expected: The validate-only run reports executable bindings for the chosen overlap; rerun with `--transactional` to create the named preview assembly.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --transactional'
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
