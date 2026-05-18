---
chapter_id: "simple_pcr_selection_gui"
title: "Simple PCR From a Selected Core Region"
tier: "core"
example_id: "simple_pcr_selection_gui"
source_example: "docs/examples/workflows/simple_pcr_selection_gui.json"
example_test_mode: "always"
executed_during_generation: true
---

# Simple PCR From a Selected Core Region

Open a ready local TP73 sequence and walk through the smallest useful PCR story: select the core ROI, limit primer distance from the core, and cap the amplicon length.

This chapter is intentionally minimal. The tutorial project opens one stable local TP73 locus so you can stay in the GUI, paint or select one region of interest, launch `Simple PCR from selection`, and reason about primer placement in beginner terms before moving on to the richer batch-PCR chapter.

See also: guided walkthrough [docs/tutorial/simple_pcr_selection_gui.md](../../simple_pcr_selection_gui.md). Use that page first when you want a human-led path; this chapter is the executable reference.

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./01_load_branch_reverse_complement_pgex_fasta.md) first.

## Parameters That Matter

- `core ROI` (where used: map selection and PCR Designer starter block)
  - Why it matters: This is the part of the biology that must definitely be included in the PCR product.
  - How to derive it: Select only the indispensable sequence interval, not the full desired amplicon.
- `max primer distance from core` (where used: PCR Designer starter block)
  - Why it matters: This sets how far upstream and downstream primers may move away from the required core region.
  - How to derive it: Choose the largest flank width you are comfortable allowing for primer search on each side of the core ROI.
- `max amplicon` (where used: Design primer pairs form)
  - Why it matters: This caps product length so the simple-PCR result stays experimentally practical.
  - How to derive it: Set the longest acceptable product length for the assay before running primer design.

## When This Routine Is Useful

- You want the shortest possible path from one selected region to one PCR design attempt.
- You want to learn the meaning of core ROI, maximum primer distance from the core, and maximum amplicon length without queueing or painting multiple regions.
- You want one offline tutorial project that opens directly into local TP73 context before using the PCR Designer.

## What You Learn

- Start a simple PCR directly from one sequence selection.
- Interpret the beginner PCR controls as deterministic flank-window constraints around the selected core ROI.
- Review saved primer pairs in terms of amplicon length and left/right distance from the core ROI.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.

## At a Glance

1. Open File -> Open Tutorial Project... -> Core -> 18. Simple PCR From a Select...
2. In the opened TP73 project, keep the DNA map in linear mode and drag over one...
3. Right-click the selection and choose Simple PCR from selection.
4. In PCR Designer, adjust max primer distance from core and max amplicon, then ...
5. Run Design Primer Pairs and inspect the in-panel primer report preview for le...

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open File -> Open Tutorial Project... -> Core -> 18. Simple PCR From a Select...

GUI: Open `File -> Open Tutorial Project... -> Core -> 18. Simple PCR From a Selected Core Region`.

CLI:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/simple_pcr_selection_gui.json
```

> Expected: The starter project loads local TP73 context as `tp73_locus` and pins the primer backend to the deterministic internal implementation.

### Step 2: In the opened TP73 project, keep the DNA map in linear mode and drag over one...

GUI: In the opened TP73 project, keep the DNA map in linear mode and drag over one short core region that must be included in the amplicon.

CLI:

```bash
# GUI-only selection gesture. For a fully scripted primer-design payload, use Chapter 13's batch PCR route.
```

> Expected: The selected core ROI is the biological interval that must be included; scripted users should encode the same interval explicitly in a primer-design request.

### Step 3: Right-click the selection and choose Simple PCR from selection

GUI: Right-click the selection and choose `Simple PCR from selection`.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'ui open pcr-design'
```

> Expected: The PCR Designer opens from the same shared UI-intent route used by agent and shell callers.

### Step 4: In PCR Designer, adjust max primer distance from core and max amplicon, then ...

GUI: In `PCR Designer`, adjust `max primer distance from core` and `max amplicon`, then click `Apply simple flank windows` if you changed the distance.

CLI:

```bash
# The scripted equivalent is a DesignPrimerPairs request with explicit core_start/core_end and flank limits.
```

> Expected: The primer-design request should carry flank-window and maximum-amplicon constraints derived from the selected core ROI.

### Step 5: Run Design Primer Pairs and inspect the in-panel primer report preview for le...

GUI: Run `Design Primer Pairs` and inspect the in-panel primer report preview for left/right distance from the core ROI and whether the pair cleanly flanks the core.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'primers list-reports'
```

> Expected: After primer design, report rows expose amplicon length plus left/right distance from the core ROI.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/simple_pcr_selection_gui.json
cargo run --bin gentle_cli -- shell 'ui open pcr-design'
cargo run --bin gentle_cli -- shell 'primers list-reports'
```

## Checkpoints

- The tutorial project opens with local TP73 sequence `tp73_locus` already loaded.
- A non-empty selection exposes `Simple PCR from selection` in the DNA-window context menu.
- The PCR Designer shows the `Simple PCR starter` block.
- After primer design, the report preview shows left/right distance from the core ROI and whether the pair cleanly flanks the core.

## Canonical Source

- Chapter id: `simple_pcr_selection_gui`
- Tier: `core`
- Example id: `simple_pcr_selection_gui`
- Workflow file: `docs/examples/workflows/simple_pcr_selection_gui.json`
- Example test_mode: `always`
- Executed during generation: `yes`
- Inspect this JSON file directly when you need full option-level detail.
