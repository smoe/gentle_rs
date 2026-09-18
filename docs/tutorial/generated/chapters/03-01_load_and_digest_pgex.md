---
chapter_id: "load_and_digest_pgex"
title: "Load pGEX and digest with BamHI/EcoRI"
tier: "core"
example_id: "load_and_digest_pgex"
source_example: "docs/examples/workflows/load_and_digest_pgex.json"
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
generated_artifact_dir: "docs/tutorial/generated/artifacts/load_and_digest_pgex"
---

# Load pGEX and digest with BamHI/EcoRI

Introduce restriction digest planning and deterministic fragment products.

Restriction digest is a core molecular cloning routine used for vector linearization, insert release, and diagnostic fragment checks. This chapter focuses on how digest parameters map to reproducible fragment sets that can later feed ligation or analysis steps.

## What You Will Accomplish

- Execute digest operations and inspect created fragment IDs.
- Reason about multi-product lineage from one parent sequence.
- Identify why stable IDs matter for follow-up ligation/extraction steps.

## Before You Start

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

**Useful when:**

- You want to verify expected restriction fragments before ordering primers or designing ligations.
- You need a reproducible digest baseline to compare against wet-lab gel expectations.
- You plan to reuse fragment IDs in later operations.

## Walkthrough: GUI, CLI and Inner Agent

Each step pairs GUI instructions with related terminal commands or guidance and a review-only inner-agent prompt. Some steps require GUI interaction; a listing command only inspects results, it does not perform the design. CLI snippets assume an installed `gentle_cli` and use GENtle's default `.gentle_state.json` unless stated otherwise. From a source checkout, replace `gentle_cli` with `cargo run --bin gentle_cli --`. Add `--state PATH` or `--project PATH` for an explicit sandbox; a separate CLI process does not inherit the open GUI project's unsaved state.

In the **GUI Shell**, enter only the shared command inside `gentle_cli shell '...'`, without the executable prefix or outer quotes. Run UI-opening commands there to open windows: a headless CLI returns the UI intent but does not open a GUI. Other terminal commands are not automatically GUI Shell commands. Inner-agent examples request a proposal for review; they are not executed during tutorial generation.

### Step 1: Load test_files/pGEX-3X.gb in the GUI and inspect annotated features

**GUI**

Load `test_files/pGEX-3X.gb` in the GUI and inspect annotated features.

**CLI (terminal)**

```bash
gentle_cli op '{"LoadFile":{"path":"test_files/pGEX-3X.gb","as_id":"pgex"}}'
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Load `test_files/pGEX-3X.gb` in the GUI and inspect annotated features. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The state contains the annotated pGEX sequence as `pgex`.

### Step 2: Open Sequence Tools from the DNA window, expand Core cloning operations, keep enzymes BamHI,EcoRI, set prefix frag, and run Digest

**GUI**

Open Sequence Tools from the DNA window, expand Core cloning operations, keep enzymes `BamHI,EcoRI`, set prefix `frag`, and run Digest.

**CLI (terminal)**

```bash
gentle_cli op '{"Digest":{"input":"pgex","enzymes":["BamHI","EcoRI"],"output_prefix":"frag"}}'
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Open Sequence Tools from the DNA window, expand Core cloning operations, keep enzymes `BamHI,EcoRI`, set prefix `frag`, and run Digest. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> The digest operation creates deterministic fragment sequence IDs using the `frag` prefix.

![Whole-screen orientation after opening Sequence Tools from the DNA viewer.](../../../screenshots/tutorial_gui_acceptance/load_and_digest_pgex/open_tools.orientation.svg)

*Figure: Whole-screen orientation after opening Sequence Tools from the DNA viewer. Screenshot captured 2026-09-08.*

![Digest controls with the enzyme list and output prefix in interaction context.](../../../screenshots/tutorial_gui_acceptance/load_and_digest_pgex/set_prefix.context.svg)

*Figure: Digest controls with the enzyme list and output prefix in interaction context. Screenshot captured 2026-09-08.*

### Step 3: Review created fragment entries and confirm they are stored as independent sequence products

**GUI**

Review created fragment entries and confirm they are stored as independent sequence products.

**CLI (terminal)**

```bash
gentle_cli workflow @docs/examples/workflows/load_and_digest_pgex.json
```

**Ask the inner agent**

> In the current GENtle project, help me perform this tutorial step: Review created fragment entries and confirm they are stored as independent sequence products. Show the exact GENtle operation or command and its expected result for my review. State any missing input. Do not execute it until I approve.

**Expected**

> Replaying the workflow reproduces the same loaded sequence and fragment-product lineage.

![Whole-screen orientation after running the digest and publishing its fragment products.](../../../screenshots/tutorial_gui_acceptance/load_and_digest_pgex/digest.orientation.svg)

*Figure: Whole-screen orientation after running the digest and publishing its fragment products. Screenshot captured 2026-09-08.*


## Interpretation and Reference

## Parameters That Matter

- `Digest.enzymes` (where used: operation 2)
  - Why it matters: The enzyme set determines cut positions and resulting fragment repertoire.
  - How to derive it: Choose enzymes based on cloning strategy (diagnostic digest vs insert release vs vector opening).
- `Digest.output_prefix` (where used: operation 2)
  - Why it matters: Controls deterministic fragment ID namespace.
  - How to derive it: Use a short routine-specific prefix (e.g., `frag`, `d`, `eco_bam`).
  - Omit when: Omit only if auto-generated IDs are acceptable for ad hoc inspection.

## Applied Concepts

- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.

## Checkpoints

- Digest operation completes and creates fragment sequence IDs.
- Fragment IDs are deterministic across repeated runs.

## Tutorial Provenance

- Chapter id: `load_and_digest_pgex`
- Tier: `core`
- Example id: `load_and_digest_pgex`
- Tutorial source JSON: `docs/tutorial/sources/03-01_load_and_digest_pgex.json`
- Workflow file: `docs/examples/workflows/load_and_digest_pgex.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/load_and_digest_pgex`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Load pGEX and digest with BamHI/EcoRI`
- Tutorial/chapter id: `load_and_digest_pgex`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
