---
chapter_id: "gibson_arrangements_baseline"
title: "Gibson Arrangements Starter Project (offline)"
tier: "core"
example_id: "gibson_arrangements_baseline"
source_example: "docs/examples/workflows/gibson_arrangements_baseline.json"
example_test_mode: "always"
executed_during_generation: true
---

# Gibson Arrangements Starter Project (offline)

Open directly into an arrangement-ready Gibson result with vector, insert, assembled product, a stored three-lane lane set, and an assistant-facing bench handoff already present.

This chapter now uses its own deterministic workflow example instead of reusing the Gibson specialist starter. `File -> Open Tutorial Project...` builds the pGEX + insert starter, applies the canonical single-insert Gibson plan, exports a lab-assistant handoff Markdown file, and then opens the arrangement guide so the user can inspect singleton outputs, stored arrangements, gel export, and bench-facing instructions without first repeating the cloning step.

See also: guided walkthrough [docs/tutorial/gibson_arrangements_gui.md](../../gibson_arrangements_gui.md). Use that page first when you want a human-led path; this chapter is the executable reference.

**Prerequisites:** Read [Chapter 15: Gibson Specialist Starter Project (offline)](./15_gibson_specialist_testing_baseline.md) first.

## Parameters That Matter

- `destination seq_id = gibson_destination_pgex` (where used: Stored vector lane in the arrangement walkthrough)
  - Why it matters: Keeps the original vector lane stable across the arrangement guide and CLI replay.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `insert seq_id = gibson_insert_demo` (where used: Insert lane in the arrangement walkthrough)
  - Why it matters: Provides the same deterministic insert lane that the Gibson arrangement is expected to reference.
  - How to derive it: This workflow assigns the ID directly during `LoadFile`.
- `assembled product seq_id = gibson_destination_pgex_with_gibson_insert_demo` (where used: Stored Gibson result already present when the arrangement tutorial opens)
  - Why it matters: Makes the arrangement tutorial disjunct from the separate Gibson-specialist apply walkthrough while keeping one deterministic product identity.
  - How to derive it: The chapter workflow applies the canonical single-insert Gibson plan before the guide opens.
- `workflow example = gibson_arrangements_baseline` (where used: `Open Tutorial Project...` and CLI workflow replay)
  - Why it matters: The arrangement starter now has its own canonical example so the tutorial can begin with the cloned result and stored lane arrangement already available.
  - How to derive it: Select the chapter from `Open Tutorial Project...` or run the canonical workflow JSON directly.
- `lab assistant handoff = artifacts/gibson_lab_assistant_handoff.md` (where used: ClawBio/Telegram demo artifact and CLI replay)
  - Why it matters: Provides non-IT bench-facing instructions tied to the exact deterministic design outputs.
  - How to derive it: The workflow calls `ExportLabAssistantInstructions` after applying the Gibson plan.

## When This Routine Is Useful

- You want to inspect the arrangement that Gibson apply creates without first navigating through the earlier Gibson-specialist apply walkthrough.
- You want a reproducible starter state for checking singleton output containers, the assembled product, and arrangement-level gel export.
- You want a deterministic example of GENtle handing a designed cloning experiment to a non-IT lab assistant.
- You want one offline tutorial-project entry that opens directly on the arrangement guide in Help/Tutorial.

## What You Learn

- Use one executable starter project to reach an arrangement-ready walkthrough directly.
- Recognize the separation between the Gibson specialist apply tutorial and the downstream arrangement/gel-inspection tutorial.
- Inspect the generated lab-assistant handoff as the bench-facing counterpart to the design state.
- Replay the same arrangement-focused setup from GUI and CLI without requiring a second manual Gibson apply.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open the Gibson arrangements starter project from File -> Open Tutorial Proje...

GUI: Open the Gibson arrangements starter project from `File -> Open Tutorial Project...`.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_arrangements_baseline.json
```

> Expected: The workflow builds the same arrangement-ready starter state that the GUI tutorial-project menu opens.

### Step 2: Confirm the starter contains gibson_destination_pgex (circular), gibson_inser...

GUI: Confirm the starter contains `gibson_destination_pgex` (circular), `gibson_insert_demo` (linear), and the assembled product `gibson_destination_pgex_with_gibson_insert_demo`.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_arrangements_baseline.json
```

> Expected: The workflow result includes the vector, insert, assembled Gibson product, stored arrangement context, and the retained handoff artifact.

### Step 3: Continue in the automatically opened Gibson Arrangements Tutorial from Step 1...

GUI: Continue in the automatically opened `Gibson Arrangements Tutorial` from `Step 1` onward.

CLI:

```bash
cargo run --bin gentle_examples_docs -- tutorial-catalog-check
```

> Expected: The tutorial catalog remains linked so the generated arrangement starter points readers to the hand-written Gibson Arrangements walkthrough.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/gibson_arrangements_baseline.json
```

## Checkpoints

- The tutorial project opens with the expected destination, insert, and assembled-product IDs already loaded.
- The workflow writes `artifacts/gibson_lab_assistant_handoff.md` with material IDs, design-derived bench steps, checkpoints, and safety scope.
- The Help window lands on the arrangement tutorial rather than the specialist testing tutorial.
- The arrangement guide can start from this baseline without an additional manual Gibson apply.

## What This Chapter Produces

- [`artifacts/gibson_arrangements_baseline/artifacts/gibson_lab_assistant_handoff.md`](../artifacts/gibson_arrangements_baseline/artifacts/gibson_lab_assistant_handoff.md) - `# Gibson assembly lab assistant handoff`

## Canonical Source

- Chapter id: `gibson_arrangements_baseline`
- Tier: `core`
- Example id: `gibson_arrangements_baseline`
- Workflow file: `docs/examples/workflows/gibson_arrangements_baseline.json`
- Example test_mode: `always`
- Executed during generation: `yes`
- Inspect this JSON file directly when you need full option-level detail.
