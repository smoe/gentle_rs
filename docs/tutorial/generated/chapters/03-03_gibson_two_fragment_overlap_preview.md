---
chapter_id: "gibson_two_fragment_overlap_preview"
title: "Gibson two-fragment overlap preflight and concatenation preview"
tier: "core"
example_id: "gibson_two_fragment_overlap_preview"
source_example: "docs/examples/workflows/gibson_two_fragment_overlap_preview.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-10-01"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/gibson_two_fragment_overlap_preview"
---

# Gibson two-fragment overlap preflight and concatenation preview

Validate a 20 bp two-fragment overlap and generate deterministic, explicitly non-final concatenation previews.

Gibson assembly in practice depends on overlap design quality. This chapter anchors that design logic in GENtle through one repeatable workflow: build two overlapping fragments, run Gibson-specific preflight checks on overlap compatibility, and produce forward-order concatenation previews for review and communication. GENtle does not yet provide a dedicated overlap-eliding Gibson assembly operation: two 320 bp inputs with a shared 20 bp overlap produce a 640 bp preview, not a 620 bp Gibson product. Use the preview to review order and diagnostics, never as an order-ready construct sequence.

## What You Will Accomplish

- Understand where Gibson overlap assumptions are checked in shared-shell preflight.
- Interpret overlap mismatch diagnostics before state mutation.
- Use deterministic output IDs for clear cloning communication without confusing a concatenation preview with an overlap-elided product sequence.

## Before You Start

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

**Useful when:**

- You want reproducible documentation of intended Gibson fragment order and overlap assumptions.
- You need fast overlap and fragment-order feedback before primer ordering or wet-lab execution, while retaining the distinction between a preview and a final construct.
- You want one routine that can be explained to collaborators through GUI + shell parity.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Create two overlapping pGEX fragments (gibson_left, gibson_right) from test_files/pGEX_3X.fa using region extraction

**GUI**

Create two overlapping pGEX fragments (`gibson_left`, `gibson_right`) from `test_files/pGEX_3X.fa` using region extraction.

**CLI (terminal)**

```bash
gentle_cli workflow @docs/examples/workflows/gibson_two_fragment_overlap_preview.json
```

**Ask the inner agent**

> Inspect `test_files/pGEX_3X.fa` and propose the exact load plus `ExtractRegion` operations for `gibson_left` at `[0,320)` and `gibson_right` at `[300,620)`. Confirm that each input is 320 bp and their intended shared span is 20 bp. Explain that the workflow's 640 bp concatenation previews retain both overlap copies; do not execute until I approve.

**Expected**

> The workflow loads `pgex_fasta`, extracts 320 bp `gibson_left` and `gibson_right` inputs sharing the intended 20 bp source span, and leaves deterministic 640 bp concatenation-preview IDs for inspection.

### Step 2: Import the Gibson two-fragment overlap preview routine from Patterns -> Routine catalog

**GUI**

Import the Gibson two-fragment overlap preview routine from `Patterns -> Routine catalog`.

**CLI (terminal)**

```bash
gentle_cli shell 'macros template-import assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json'
```

**Ask the inner agent**

> Check whether `gibson_two_fragment_overlap_preview` is registered. If it is absent, propose the exact template import from `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`; treat registry mutation separately and do not import until I approve.

**Expected**

> The routine catalog registers `gibson_two_fragment_overlap_preview`, so the next shell command can validate that specific fragment order.

### Step 3: Validate the Gibson overlap preview from Shell, then rerun the same template without --validate-only to create outputs

**GUI**

Validate the Gibson overlap preview from `Shell`, then rerun the same template without `--validate-only` to create outputs.

**CLI (terminal)**

```bash
gentle_cli shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
```

**Ask the inner agent**

> Propose a validate-only run for `gibson_left`, `gibson_right`, `overlap_bp=20`, `assembly_prefix=gibson_demo` and `output_id=gibson_demo_forward`. Return the complete preflight report and stop. If I later request execution, propose a separate transactional run and state that its 640 bp output is a concatenation preview, not an overlap-elided Gibson product.

**Expected**

> The validate-only run reports executable bindings for the chosen overlap; a separately approved `--transactional` run creates the named 640 bp concatenation preview, not a final Gibson product.


## Ask an Outer Agent (MCP or ClawBio/OpenClaw)

An outer agent does not inherit the unsaved GUI project. Give it this chapter's canonical workflow and an explicit disposable state path; ask it to retain the structured result, artifacts and reproducibility receipt instead of replacing them with prose.

> Use GENtle's `gentle-cloning` skill to replay `docs/examples/workflows/gibson_two_fragment_overlap_preview.json` against a new disposable state. First report the exact workflow, inputs, state path, outputs and whether the selected route needs confirmation. Do not infer state from an open GUI. Return the structured result, produced artifacts and reproducibility receipt, and state any unmet prerequisite.

Equivalent direct structured request for the generic wrapper:

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "workflow",
  "state_path": "/tmp/gentle-gibson-two-fragment-overlap-preview.state.json",
  "workflow_path": "docs/examples/workflows/gibson_two_fragment_overlap_preview.json",
  "timeout_secs": 300
}
```

Submitting a direct structured request is an explicit wrapper invocation. If natural language selects a narrower delegated skill and that route mutates state, selects biological material or writes artifacts, the caller must preserve that skill's proposal/approval boundary and approve only the exact bound digest. This tutorial generation step does not invoke an agent or grant approval.

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
  - Omit when: Omit only when default IDs are acceptable for exploratory runs; regardless of naming, these outputs remain concatenation previews.

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

- Validate-only run returns `can_execute=true` for the two 320 bp fragments with a matching 20 bp overlap.
- Mismatch overlap inputs return explicit Gibson overlap diagnostics.
- A separately approved mutating run creates deterministic 640 bp preview outputs (`${assembly_prefix}_*` and `${output_id}`); an overlap-elided 620 bp Gibson product is intentionally not claimed.
- The canonical outer-agent workflow receipt covers input preparation and concatenation previews only; retain the separate validate-only macro preflight report as the overlap evidence.

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
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-10-01`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Gibson two-fragment overlap preflight and concatenation preview`
- Tutorial/chapter id: `gibson_two_fragment_overlap_preview`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
