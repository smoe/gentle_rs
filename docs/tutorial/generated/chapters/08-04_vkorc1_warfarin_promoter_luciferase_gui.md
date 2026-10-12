---
chapter_id: "vkorc1_warfarin_promoter_luciferase_gui"
title: "VKORC1 / rs9923231 PGx Alert -> Mammalian Luciferase Reporter (GUI Tutorial with Matching CLI Commands)"
tier: "core"
example_id: "vkorc1_pair_gui_starter"
source_example: "docs/examples/workflows/vkorc1_pair_gui_starter.json"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/vkorc1_warfarin_promoter_luciferase_gui"
---

# VKORC1 / rs9923231 PGx Alert -> Mammalian Luciferase Reporter (GUI Tutorial with Matching CLI Commands)

Offline synthetic guard for the online tutorial's explicit multiallelic insert-pair action. The 20-base input is not human VKORC1 DNA.

This companion deliberately isolates step 4 from genome retrieval and assay interpretation. Its one variation is C at zero-based position 6 with alternate candidates A,G,T. Nothing selects T automatically. The GUI must start without outputs, take an explicitly reviewed T, and save both inserts with bases and annotations matching the separate engine oracle. Refusal atomicity is also covered by adapter regressions; a native positive pass is not live-agent, online-genome, reporter or scientific acceptance.

See also: guided walkthrough [docs/tutorial/08-04_vkorc1_warfarin_promoter_luciferase_gui.md](../../08-04_vkorc1_warfarin_promoter_luciferase_gui.md). Use that page first when you want a human-led path; this chapter is the executable reference.

## What You Will Accomplish

- Distinguish genomic-forward allele choice from transcript strand or clinical nomenclature.
- Verify one-base pair identity through persisted engine data instead of whole-map screenshot similarity.

## Before You Start

**Useful when:**

- Practice an explicit genomic-forward T choice in a synthetic multiallelic example before working with the public VKORC1 locus.
- Check that GUI-created reference and alternate inserts preserve the source and differ at one reviewed base, without interpreting the toy sequence as an assay.

## At a Glance

1. Open the synthetic offline starter project, then open vkorc1_rs9923231_promoter_fragment (20 bp, not human DNA).
2. Expand the variation feature group and select its exact rs9923231 feature (index 0).
3. Use Open Promoter Design from the selected feature's description panel.
4. Scroll to Alternate base and explicitly enter genomic-forward T; never treat a bare multiallelic alternate as an approved choice.
5. Scroll to Make reference/alternate inserts, create the pair, save, and compare the source and both outputs with the independent oracle.

## Walkthrough: GUI, CLI and Inner Agent

1. Open the synthetic offline starter project, then open vkorc1_rs9923231_promoter_fragment (20 bp, not human DNA).
2. Expand the variation feature group and select its exact rs9923231 feature (index 0).
3. Use Open Promoter Design from the selected feature's description panel.
4. Scroll to Alternate base and explicitly enter genomic-forward T; never treat a bare multiallelic alternate as an approved choice.
5. Scroll to Make reference/alternate inserts, create the pair, save, and compare the source and both outputs with the independent oracle.

## Complete Workflow Replay

When an individual GUI gesture has no standalone shell command, replay the complete canonical workflow:

```bash
gentle_cli workflow @docs/examples/workflows/vkorc1_pair_gui_starter.json
gentle_cli shell 'workflow @docs/examples/workflows/vkorc1_pair_gui_starter.json'
```

## Interpretation and Reference

## Parameters That Matter

- This chapter intentionally avoids additional command options; run the canonical workflow unchanged.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.

## Checkpoints

- Starter has no pair outputs.
- Explicit T creates one C/T difference at position 6 without editing the source.
- Synthetic/native evidence does not close the online or human-scientific gates.

## Tutorial Provenance

- Chapter id: `vkorc1_warfarin_promoter_luciferase_gui`
- Tier: `core`
- Example id: `vkorc1_pair_gui_starter`
- Tutorial source JSON: `docs/tutorial/sources/08-04_vkorc1_warfarin_promoter_luciferase_gui.json`
- Workflow file: `docs/examples/workflows/vkorc1_pair_gui_starter.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/vkorc1_warfarin_promoter_luciferase_gui`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `VKORC1 / rs9923231 PGx Alert -> Mammalian Luciferase Reporter (GUI Tutorial with Matching CLI Commands)`
- Tutorial/chapter id: `vkorc1_warfarin_promoter_luciferase_gui`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
