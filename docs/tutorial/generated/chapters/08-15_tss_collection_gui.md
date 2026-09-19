---
chapter_id: "tss_collection_gui"
title: "From annotated transcript starts to reusable TSS windows (offline)"
tier: "core"
example_id: "tss_collection_gui_oracle"
source_example: "docs/examples/workflows/tss_collection_gui_oracle.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "codex_reviewed"
review_stale: false
codex_reviewed_at: "2026-09-19"
human_reviewed_at: null
human_reviewer: null
review_stale_reason: null
review_issue_template: null
review_issue_template_path: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/tss_collection_gui"
---

# From annotated transcript starts to reusable TSS windows (offline)

Preview and approve annotated starts, validate and reopen their windows, then cancel or confirm metadata-only forgetting and undo it.

Two transcripts can share a first base without being the same transcript. In this artificial 2 kb locus, plus_a and plus_b share start 601; plus_c starts at 901; minus_a starts at 1500 on the opposite strand. Their common TOY label selects two separate gene IDs without merging them. Three 701 bp windows result, oriented in transcript direction with the annotated start at local base 501. These annotations are not experimental evidence for active promoters. The repetitive toy DNA is not suitable for primer or biological TFBS conclusions. The generated passing status refers to engine workflow replay only, not to a completed Linux/Xvfb acceptance run.

See also: guided walkthrough [docs/tutorial/08-15_tss_collection_gui.md](../../08-15_tss_collection_gui.md). Use that page first when you want a human-led path; this chapter is the executable reference.

## What You Will Accomplish

- Interpret transcript-oriented starts without inferring a preferred biological TSS.
- Use preview-bound approval and identity-preserving collection operations.
- Retain typed evidence independently from screenshots.

## Before You Start

**Useful when:**

- Learn why a transcript count differs from a distinct-start count.
- Separate an unchecked stored record from a currently validated collection.
- Practice without network access, private genes, model services or a prepared human genome.

## At a Glance

1. Prepare the synthetic starter with the command in the companion walkthrough; open its project and double-click tss_locus. Never open the oracle as the starter.
2. Choose TFBS scan > Transcript starts / TSS windows.
3. Set Gene to TOY; keep Collection ID tss_windows and upstream/downstream 500/200. Inspect starts without changing project state.
4. Check starts 601 (+), 901 (+), 1500 (-), with plus_a/plus_b sharing the first row. Select available and explicitly approve creation of all three windows.
5. Refresh collections: readable / not checked is expected, not pass. Inspect stored collection; compare the validated report with the independent oracle.
6. Open TSS collection twice. Expect exactly three member viewers plus the original locus, with no duplicate windows. Compare persisted members with the independent oracle.
7. Request Forget registry entry, then Cancel. Verify the collection still validates. Request again and confirm the named entry: the registry disappears, but all sequences and open viewers remain.
8. Choose Edit > Undo in the main window. Inspect stored collection again: its report and member content must match the original oracle.
9. Continue the separate manual checklist for save/close/reopen and deliberately editing a member. Rejection of a stale member after undoing forget is not yet certified by the automated GUI subset.

## Walkthrough: GUI, CLI and Inner Agent

### Step 1: Prepare the synthetic starter with the command in the companion walkthrough; open its project and double-click tss_locus. Never open the oracle as the starter

**GUI**

Prepare the synthetic starter with the command in the companion walkthrough; open its project and double-click tss_locus. Never open the oracle as the starter.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Prepare the synthetic starter with the command in the companion walkthrough; open its project and double-click tss_locus. Never open the oracle as the starter. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The starter has one locus and no TSS collection.

### Step 2: Choose TFBS scan > Transcript starts / TSS windows

**GUI**

Choose TFBS scan > Transcript starts / TSS windows.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Choose TFBS scan > Transcript starts / TSS windows. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Opening the form does not approve or derive windows.

### Step 3: Set Gene to TOY; keep Collection ID tss_windows and upstream/downstream 500/200. Inspect starts without changing project state

**GUI**

Set Gene to TOY; keep Collection ID tss_windows and upstream/downstream 500/200. Inspect starts without changing project state.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Set Gene to TOY; keep Collection ID tss_windows and upstream/downstream 500/200. Inspect starts without changing project state. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Preview exposes three exact starts and its approval digest.

### Step 4: Check starts 601 (+), 901 (+), 1500 (-), with plus_a/plus_b sharing the first row. Select available and explicitly approve creation of all three windows

**GUI**

Check starts 601 (+), 901 (+), 1500 (-), with plus_a/plus_b sharing the first row. Select available and explicitly approve creation of all three windows.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Check starts 601 (+), 901 (+), 1500 (-), with plus_a/plus_b sharing the first row. Select available and explicitly approve creation of all three windows. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Only the reviewed preview and explicit selection authorize derivation.

### Step 5: Refresh collections: readable / not checked is expected, not pass. Inspect stored collection; compare the validated report with the independent oracle

**GUI**

Refresh collections: readable / not checked is expected, not pass. Inspect stored collection; compare the validated report with the independent oracle.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Refresh collections: readable / not checked is expected, not pass. Inspect stored collection; compare the validated report with the independent oracle. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The validating engine report contains three 701 bp windows; not checked and validated are distinct.

### Step 6: Open TSS collection twice. Expect exactly three member viewers plus the original locus, with no duplicate windows. Compare persisted members with the independent oracle

**GUI**

Open TSS collection twice. Expect exactly three member viewers plus the original locus, with no duplicate windows. Compare persisted members with the independent oracle.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Open TSS collection twice. Expect exactly three member viewers plus the original locus, with no duplicate windows. Compare persisted members with the independent oracle. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Repeated opening reuses the same four subject-bound DNA viewers; queued is not yet opened.

### Step 7: Request Forget registry entry, then Cancel. Verify the collection still validates. Request again and confirm the named entry: the registry disappears, but all sequences and open viewers remain

**GUI**

Request Forget registry entry, then Cancel. Verify the collection still validates. Request again and confirm the named entry: the registry disappears, but all sequences and open viewers remain.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Request Forget registry entry, then Cancel. Verify the collection still validates. Request again and confirm the named entry: the registry disappears, but all sequences and open viewers remain. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Cancellation changes nothing; confirmed forgetting removes only collection registry metadata, not the member sequences.

### Step 8: Choose Edit > Undo in the main window. Inspect stored collection again: its report and member content must match the original oracle

**GUI**

Choose Edit > Undo in the main window. Inspect stored collection again: its report and member content must match the original oracle.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Choose Edit > Undo in the main window. Inspect stored collection again: its report and member content must match the original oracle. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Undo restores the collection; a fresh inspection validates the restored members.

### Step 9: Continue the separate manual checklist for save/close/reopen and deliberately editing a member. Rejection of a stale member after undoing forget is not yet certified by the automated GUI subset

**GUI**

Continue the separate manual checklist for save/close/reopen and deliberately editing a member. Rejection of a stale member after undoing forget is not yet certified by the automated GUI subset.

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Continue the separate manual checklist for save/close/reopen and deliberately editing a member. Rejection of a stale member after undoing forget is not yet certified by the automated GUI subset. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Manual restart and stale-member GUI acceptance remain separately recorded.


## Complete Workflow Replay

When an individual GUI gesture has no standalone shell command, replay the complete canonical workflow:

```bash
gentle_cli workflow @docs/examples/workflows/tss_collection_gui_oracle.json
gentle_cli shell 'workflow @docs/examples/workflows/tss_collection_gui_oracle.json'
```

## Interpretation and Reference

## Parameters That Matter

- `upstream_bp / downstream_bp` (where used: TSS inventory and materialization)
  - Why it matters: The annotated start itself adds one base: 500 + 1 + 200 = 701. Minus-strand windows reverse-complement their genomic interval.
  - How to derive it: Use the reviewed biological question and loaded flanks; this tutorial fixes 500/200 for an entirely synthetic locus.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## Checkpoints

- Starter completion is unsatisfied; the independent oracle satisfies it.
- Typed collection report and persisted sequences agree, including both strands.
- No claim of live Linux or TP73 acceptance is made until Glen runs and reviews the retained evidence.

## Tutorial Provenance

- Chapter id: `tss_collection_gui`
- Tier: `core`
- Example id: `tss_collection_gui_oracle`
- Tutorial source JSON: `docs/tutorial/sources/08-15_tss_collection_gui.json`
- Workflow file: `docs/examples/workflows/tss_collection_gui_oracle.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/tss_collection_gui`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `codex_reviewed`
- Codex reviewed at: `2026-09-19`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `From annotated transcript starts to reusable TSS windows (offline)`
- Tutorial/chapter id: `tss_collection_gui`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
