---
chapter_id: "simple_pcr_selection_gui"
title: "Simple PCR From a Selected Core Region"
tier: "core"
example_id: "simple_pcr_selection_gui"
source_example: "docs/examples/workflows/simple_pcr_selection_gui.json"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/simple_pcr_selection_gui"
---

# Simple PCR From a Selected Core Region

Open an 800-base local TP73 extract and walk through the smallest useful PCR story: select the core ROI, limit primer distance from the core, and cap the amplicon length.

This chapter uses the 800-base interval [61520, 62320) of the committed TP73 locus, shared with the scripted PCR oracle. Open `tp73_locus`, the compact extract; `simple_pcr_source_locus` retains the full input for provenance. The fixed smoke core is extract bases 201..600 (1-based inclusive), or [200, 600) in engine coordinates. This bounds template and flank search sizes without relaxing primer rules or precomputing a GUI result.

See also: guided walkthrough [docs/tutorial/04-01_simple_pcr_selection_gui.md](../../04-01_simple_pcr_selection_gui.md). Use that page first when you want a human-led path; this chapter is the executable reference.

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

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
2. Open tp73_locus (800 bases) and select =201 .. 600 in linear mode.
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

> Expected: The starter project contains the 800-base TP73 extract `tp73_locus` and its full source locus; open the extract. The primer backend is pinned to the deterministic internal implementation.

### Step 2: Open tp73_locus (800 bases) and select =201 .. 600 in linear mode

GUI: Open `tp73_locus` (800 bases) and select `=201 .. 600` in linear mode.

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

![Focused selection context menu with the Simple PCR action and enough sequence-map context for orientation.](../../../screenshots/tutorial_gui_acceptance/simple_pcr_selection_gui/open_selection_context.context.svg)

*Figure: Focused selection context menu with the Simple PCR action and enough sequence-map context for orientation. Screenshot captured 2026-09-08.*

![Whole-screen orientation for opening the selected core region's context menu.](../../../screenshots/tutorial_gui_acceptance/simple_pcr_selection_gui/open_selection_context.orientation.svg)

*Figure: Whole-screen orientation for opening the selected core region's context menu. Screenshot captured 2026-09-08.*

### Step 4: In PCR Designer, adjust max primer distance from core and max amplicon, then ...

GUI: In `PCR Designer`, adjust `max primer distance from core` and `max amplicon`, then click `Apply simple flank windows` if you changed the distance.

CLI:

```bash
# The scripted equivalent is a DesignPrimerPairs request with explicit core_start/core_end and flank limits.
```

> Expected: The primer-design request should carry flank-window and maximum-amplicon constraints derived from the selected core ROI.

![PCR Designer starter controls after seeding the selected core region.](../../../screenshots/tutorial_gui_acceptance/simple_pcr_selection_gui/seed_simple_pcr.context.svg)

*Figure: PCR Designer starter controls after seeding the selected core region. Screenshot captured 2026-09-08.*

### Step 5: Run Design Primer Pairs and inspect the in-panel primer report preview for le...

GUI: Run `Design Primer Pairs` and inspect the in-panel primer report preview for left/right distance from the core ROI and whether the pair cleanly flanks the core.

CLI:

```bash
cargo run --bin gentle_cli -- shell 'primers list-reports'
```

> Expected: After primer design, report rows expose amplicon length plus left/right distance from the core ROI.

![Whole-screen orientation after primer design, with the PCR Designer report and project lineage visible together.](../../../screenshots/tutorial_gui_acceptance/simple_pcr_selection_gui/design_primers.orientation.svg)

*Figure: Whole-screen orientation after primer design, with the PCR Designer report and project lineage visible together. Screenshot captured 2026-09-08.*


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/simple_pcr_selection_gui.json
cargo run --bin gentle_cli -- shell 'ui open pcr-design'
cargo run --bin gentle_cli -- shell 'primers list-reports'
```

## Checkpoints

- The tutorial project contains the 800-base TP73 extract `tp73_locus`; the full source remains available for provenance.
- A non-empty selection exposes `Simple PCR from selection` in the DNA-window context menu.
- The PCR Designer shows the `Simple PCR starter` block.
- After primer design, the report preview shows left/right distance from the core ROI and whether the pair cleanly flanks the core.

## Tutorial Provenance

- Chapter id: `simple_pcr_selection_gui`
- Tier: `core`
- Example id: `simple_pcr_selection_gui`
- Tutorial source JSON: `docs/tutorial/sources/04-01_simple_pcr_selection_gui.json`
- Workflow file: `docs/examples/workflows/simple_pcr_selection_gui.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/simple_pcr_selection_gui`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Simple PCR From a Selected Core Region`
- Tutorial/chapter id: `simple_pcr_selection_gui`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
