# Gibson Specialist Testing Tutorial

> Type: `GUI walkthrough + CLI parity`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> shared `gibson preview` route, exported plan JSON, and deterministic local
> tutorial inputs.

This tutorial is a practical test script for the dedicated `Patterns ->
Gibson...` specialist window.

It is not meant to produce a bench-ready cloning design. Instead, it gives you
one stable local path for checking that the Gibson specialist does what it
should do:

1. accepts a circular destination plus one linear insert,
2. derives left and right junction overlaps from an opening,
3. suggests Gibson PCR primers from high-level overlap/Tm targets,
4. prepares the factual Gibson cartoon payload for export,
5. exports machine-readable plan/preview outputs, and
6. matches the shared-shell/CLI `gibson preview` result.

Current limitation to keep in mind while testing: multi-insert Gibson execution
currently requires a defined destination opening; `existing_termini` remains
the single-fragment handoff path.

## Fastest Setup

Recommended route:

1. `File -> Open Tutorial Project...`
2. choose `Gibson Specialist Starter Project (offline)`
3. the Help window should open automatically on this guide
4. continue below from `Step 3`

That tutorial-project baseline is backed by the executable workflow example:

- [`docs/examples/workflows/gibson_specialist_testing_baseline.json`](../examples/workflows/gibson_specialist_testing_baseline.json)

It loads these stable sequence IDs for you:

- `gibson_destination_pgex`
- `gibson_insert_demo`

## What You Will Test

By the end of this tutorial, you should have verified all of these:

- the `Patterns -> Gibson...` window opens and stays Gibson-specific
- destination and insert selection work from local sequences
- the specialist help text calls out the current multi-insert `defined opening`
  guardrail before apply
- defined-site opening coordinates are accepted
- `Preview Gibson Plan` yields:
  - two resolved junctions
  - two primer suggestions
  - a resolved Gibson cartoon payload that can be exported as SVG
- `Export Plan JSON...`, `Export Preview JSON...`, `Export Primer Summary...`,
  and `Export Cartoon SVG...` all work
- the exported plan can be replayed through:
  - `gentle_cli ... gibson preview @plan.json --output preview.json`
- GUI and CLI preview payloads are identical

## Inputs

Use only local committed inputs:

- destination vector:
  [`test_files/pGEX-3X.gb`](../../test_files/pGEX-3X.gb)
- synthetic insert:
  [`docs/tutorial/inputs/gibson_insert_demo.fa`](./inputs/gibson_insert_demo.fa)

Optional synthetic sanity-check pair for primer/overlap reasoning:

- circular destination:
  [`docs/tutorial/inputs/gibson_synthetic_c100_a_g100_vector.gb`](./inputs/gibson_synthetic_c100_a_g100_vector.gb)
- linear insert:
  [`docs/tutorial/inputs/gibson_synthetic_poly_a_insert.fa`](./inputs/gibson_synthetic_poly_a_insert.fa)
- companion note:
  [`docs/tutorial/inputs/README.md`](./inputs/README.md)

Why these inputs:

- `pGEX-3X.gb` gives you a known circular destination sequence.
- `gibson_insert_demo.fa` is a small synthetic linear insert built only for
  functional testing of the specialist UI and preview/export flow.
- the synthetic `C^100 A G^100` + poly-A pair gives you a deliberately simple
  sequence context for reasoning about overlap derivation and insert-primer
  structure without richer plasmid annotations getting in the way

## Step 1: Load the Two Input Sequences

Skip this step if you started from `Open Tutorial Project...`.

GUI:

1. start GENtle
2. `File -> Open Sequence...`
3. load
   [`test_files/pGEX-3X.gb`](../../test_files/pGEX-3X.gb)
4. `File -> Open Sequence...`
5. load
   [`docs/tutorial/inputs/gibson_insert_demo.fa`](./inputs/gibson_insert_demo.fa)

What to verify:

- the destination sequence reports as circular
- the insert sequence reports as linear

## Step 2: Save a Test Project

Do this before the actual Gibson preview so the later CLI parity check can read
the same sequence IDs and state.

GUI:

1. `File -> Save Project As...`
2. save to a temporary location, for example:
   `~/Desktop/gibson_ui_test.project.gentle.json`

## Step 3: Open the Gibson Specialist

GUI:

1. `Patterns -> Gibson...`

What to verify:

- the specialist window opens as a dedicated Gibson window
- the window contains these sections:
  - `Destination`
  - `Insert`
  - `Design Targets`
  - `Tm Model`
  - `Review`
  - `Outputs`

## Step 4: Fill the Destination Section

Use the circular `pGEX` sequence as the destination.

If you started from the tutorial-project baseline, the expected destination
ID is `gibson_destination_pgex`.

GUI:

1. in `Destination`, choose the circular `pGEX` sequence you loaded in Step 1
   or the tutorial-project destination `gibson_destination_pgex`
2. keep `opening = defined cut/opening`
3. preferred route: in `Suggested openings`, click `Use` on the
   `SmaI (MCS) | blunt | ...` row if it is offered from the MCS annotation
4. fallback route if you want a fixed coordinate-only test:
   - `left cut edge = 941`
   - `right cut edge = 941`

Optional prefill check:

1. in the `pGEX` DNA window, make a short selection
2. return to `Patterns -> Gibson...`
3. click `Use Active Selection`
4. confirm that the destination and opening coordinates update from the focused
   DNA window

Notes:

- if the destination already carries MCS/unique-cutter knowledge, the specialist
  should now start with unique cutters named by the MCS and only reveal other
  unique cutters when you ask for them explicitly
- the `Feature` column should show where the cutter sits, for example the MCS
  location and/or an overlapping gene label
- equal `left cut edge` / `right cut edge` means a cleavage point rather
  than a removed span
- for sticky cutters, different left/right cut edges correspond to the two
  primer-relevant termini on the destination arms
- `941..941` is a stable coordinate-only fallback for the `SmaI` cutpoint on
  the bundled pGEX tutorial destination when you do not want to rely on the
  suggestion table.
- This tutorial checks the mechanics of the specialist, not whether this exact
  insertion site is biologically optimal.

## Step 5: Fill the Insert Section

GUI:

1. in `Insert`, choose the `gibson_insert_demo` sequence
2. keep `orientation = forward`

## Step 6: Set Design Targets

For a tighter, easier-to-review preview, use fixed tutorial values instead of
the broader defaults.

GUI:

1. in `Design Targets`, set:
   - `overlap min/max bp = 30 / 30`
   - `minimum overlap Tm C = 60`
   - `priming Tm window C = 58 / 68`
   - `priming length bp = 18 / 35`
   - `output id hint = gibson_ui_test_product`

What this means:

- both terminal junctions should resolve to `30 bp` overlaps
- primer suggestions still choose their exact `3'` priming segment from the
  insert template within the allowed Tm/length window

## Step 7: Preview the Gibson Plan

GUI:

1. click `Preview Gibson Plan`

What to verify in `Review`:

- `preview schema` is `gentle.gibson_assembly_preview.v1`
- `can_execute=true`
- no blocking errors are present
- `Resolved junctions` contains exactly `2` rows
- each resolved junction shows:
  - one terminal destination/insert member pair
  - `30 bp` overlap length
  - an explicit overlap sequence
- `Primer suggestions` contains exactly `2` rows
- each primer row shows:
  - one full primer sequence
  - one `5' overlap`
  - one `3' priming` segment
- `Cartoon preview` may or may not render inline; for this tutorial, the
  textual review blocks and the exported SVG are the canonical checks

Expected interpretation:

- the destination opening defines the left and right vector ends
- the overlap sequences are derived from the destination flanks
- the insert does not need to already contain those overlaps because the primer
  suggestions add them at the `5'` end

## Step 8: Apply and Export the Outputs

GUI:

1. click `Apply Gibson Cloning`
2. verify that new sequence nodes appear for:
   - the left insert primer
   - the right insert primer
   - the assembled product
3. in lineage, click the Gibson operation itself and confirm that the
   specialist reopens with the saved plan loaded again
4. click `Export Plan JSON...`
5. save as:
   `~/Desktop/gibson_ui_test.plan.json`
6. click `Export Preview JSON...`
7. save as:
   `~/Desktop/gibson_ui_test.preview.json`
8. click `Export Primer Summary...`
9. save as:
   `~/Desktop/gibson_ui_test.primers.txt`
10. click `Export Cartoon SVG...`
11. save as:
   `~/Desktop/gibson_ui_test.cartoon.svg`

What to verify:

- applying creates the expected primer/product sequence nodes
- clicking the Gibson operation reopens the specialist without silently
  executing again
- all four files are written successfully
- the primer summary is human-readable
- the SVG opens and matches the resolved cartoon logic described by the
  textual review blocks

## Step 9: Check Shared CLI Parity

Run the exported plan through the shared CLI path against the same saved
project:

```bash
cargo run --quiet --bin gentle_cli -- \
  --project "${HOME}/Desktop/gibson_ui_test.project.gentle.json" \
  gibson preview \
  @"${HOME}/Desktop/gibson_ui_test.plan.json" \
  --output "${HOME}/Desktop/gibson_ui_test.preview.from_cli.json"
```

The two input files are different on purpose:

- `gibson_ui_test.project.gentle.json` is the project you saved via
  `Save Project As...`
- `gibson_ui_test.plan.json` is the Gibson plan you exported earlier via
  `Export Plan JSON...` in Step 8

Then compare the GUI-exported preview with the CLI-generated preview:

```bash
diff -u \
  "${HOME}/Desktop/gibson_ui_test.preview.json" \
  "${HOME}/Desktop/gibson_ui_test.preview.from_cli.json"
```

Expected result:

- `diff` prints nothing
- the shared GUI and CLI routes are therefore using the same deterministic
  Gibson preview contract

## Step 10: Check the Current Handoff Limitation

GUI:

1. return to `Patterns -> Gibson...`
2. click `Open in Routine Assistant`

Expected result in this tutorial:

- the status line should explain that this preview is still preview-only for
  the current opening mode

Why this is expected:

- this tutorial uses a circular destination with a `defined site` opening
- current Gibson specialist v1 only offers direct routine handoff for the
  already-linear `existing termini` path

## Pass/Fail Checklist

Mark the tutorial successful if all of these are true:

- [ ] `Patterns -> Gibson...` opens and shows all five sections
- [ ] `Tm Model` is visible as a dedicated box
- [ ] circular destination + linear insert can be chosen without ambiguity
- [ ] the specialist warning text says multi-insert execution currently
      requires a defined destination opening and that `existing_termini`
      remains the single-fragment handoff path
- [ ] defined-site opening `941..941` is accepted
- [ ] `Opening sketch` shows the exact cut sequence plus the destination after opening
- [ ] preview returns `2` resolved junctions
- [ ] preview returns `2` primer suggestions
- [ ] `Opening sketch` updates to show the two resolved 5' destination ends with overlap sequences
- [ ] cartoon SVG export succeeds and matches the resolved review text
- [ ] `Cancel` closes the specialist without applying anything
- [ ] Gibson apply creates primer/product sequence nodes
- [ ] clicking the Gibson operation reopens the specialist
- [ ] plan, preview, primer summary, and cartoon export all succeed
- [ ] CLI `gibson preview` reproduces the GUI preview byte-for-byte
- [ ] routine handoff reports the current v1 limitation clearly instead of
      failing silently

## Optional Extension

If you also want to test the active-context helpers more directly:

1. focus the `pGEX` DNA window
2. click `Use Active Sequence`
3. make a short selection in the DNA window
4. click `Use Active Selection`

That confirms the specialist can prefill from the currently focused DNA window
instead of forcing you to type everything manually.
