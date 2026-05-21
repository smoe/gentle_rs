---
chapter_id: "gibson_specialist_testing_baseline"
title: "Gibson Specialist Starter Project (offline)"
tier: "core"
example_id: "gibson_specialist_testing_baseline"
source_example: "docs/examples/workflows/gibson_specialist_testing_baseline.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "unreviewed"
codex_reviewed_at: null
human_reviewed_at: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/gibson_specialist_testing_baseline"
---

# Gibson Specialist Starter Project (offline)

Build a ready-to-open tutorial project with one circular destination and one linear insert so the Gibson specialist can be tested without setup drift.

This chapter exists to remove setup friction from Gibson specialist testing. The executable workflow loads a known circular pGEX destination and one small synthetic insert under stable IDs, so `File -> Open Tutorial Project...` can hand you a reproducible starting state before you continue in the dedicated GUI specialist.

See also: guided walkthrough [docs/tutorial/gibson_specialist_testing_gui.md](../../gibson_specialist_testing_gui.md). Use that page first when you want a human-led path; this chapter is the executable reference.

**Prerequisites:** Read [Chapter 4: Gibson two-fragment overlap planning baseline](./04_gibson_two_fragment_overlap_preview.md) first.

## Parameters That Matter

- `destination seq_id = gibson_destination_pgex` (where used: Gibson specialist destination selection after tutorial project open)
  - Why it matters: Keeps the circular vector choice stable across tutorial-project, GUI, and CLI replay paths.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `insert seq_id = gibson_insert_demo` (where used: Gibson specialist insert selection after tutorial project open)
  - Why it matters: Provides one deterministic linear insert for overlap and primer derivation checks.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `workflow example = gibson_specialist_testing_baseline` (where used: `Open Tutorial Project...` and CLI workflow replay)
  - Why it matters: This is the canonical executable setup layer behind the hand-written Gibson testing walkthrough.
  - How to derive it: Select the chapter from `Open Tutorial Project...` or run the canonical workflow JSON directly.

## When This Routine Is Useful

- You want to test `Patterns -> Gibson...` without manually importing files first.
- You want a stable tutorial-project baseline that can be opened from the GUI and replayed from the CLI.
- You want one deterministic starting state for checking overlaps, primer suggestions, cartoon rendering, and export behavior.

## What You Learn

- Use executable tutorial projects as deterministic GUI test baselines.
- Recognize the stable sequence IDs that the Gibson testing guide expects.
- Replay the same baseline from GUI and CLI without changing biological setup.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open the Gibson specialist starter project from File -> Open Tutorial Project...

GUI: Open the Gibson specialist starter project from `File -> Open Tutorial Project...`.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json
```

> Expected: The workflow loads the same starter inputs that the GUI tutorial-project menu opens, without requiring manual sequence import.

### Step 2: Confirm the starter contains gibson_destination_pgex (circular) and gibson_in...

GUI: Confirm the starter contains `gibson_destination_pgex` (circular) and `gibson_insert_demo` (linear).

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json
```

> Expected: The workflow result contains `gibson_destination_pgex` from `test_files/pGEX-3X.gb` and `gibson_insert_demo` from the tutorial insert FASTA.

### Step 3: Continue in the automatically opened Gibson Specialist Testing Tutorial from ...

GUI: Continue in the automatically opened `Gibson Specialist Testing Tutorial` from `Step 3` onward.

CLI:

```bash
cargo run --bin gentle_examples_docs -- tutorial-catalog-check
```

> Expected: The tutorial catalog remains linked so the generated starter chapter points readers to the hand-written Gibson specialist walkthrough.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json
```

## Checkpoints

- The tutorial project opens with exactly the expected destination and insert IDs already loaded.
- The destination remains circular and the insert remains linear before entering `Patterns -> Gibson...`.
- The hand-written Gibson testing guide can start from this baseline without any additional sequence imports.

## Tutorial Provenance

- Chapter id: `gibson_specialist_testing_baseline`
- Tier: `core`
- Example id: `gibson_specialist_testing_baseline`
- Tutorial source JSON: `docs/tutorial/sources/25_gibson_specialist_testing_baseline.json`
- Workflow file: `docs/examples/workflows/gibson_specialist_testing_baseline.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/gibson_specialist_testing_baseline`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Gibson Specialist Starter Project (offline)`
- Tutorial/chapter id: `gibson_specialist_testing_baseline`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
