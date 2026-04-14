# Reverse Translate an Imported Protein and Audit the Result

> Type: `GUI walkthrough`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> current `Protein Evidence...` specialist, the persisted reverse-translation
> report, and the lineage reopen flow.

This tutorial is meant as the second **manual check** for the newer
protein-design path.

It verifies that GENtle can:

- import one first-class protein sequence
- reverse translate it into coding DNA
- show resolved translation-table and speed-profile provenance
- record the result in lineage so the created coding sequence remains easy to
  reopen and audit

## What You Will Check

By the end of this walkthrough, you should be able to confirm that:

- Ensembl protein fetch/import still works in `Protein Evidence...`
- reverse translation creates one coding DNA sequence
- the result panel states the resolved translation-table and speed-profile
  provenance clearly
- the created coding sequence is linked back into lineage as a reverse-
  translation analysis artifact

## Inputs

Start from the same local project fixture as the first tutorial:

- [`test_files/tp73.project.gentle.json`](../../test_files/tp73.project.gentle.json)

This walkthrough also requires online access for the Ensembl fetch step.

Known live smoke-check identifier already used in this codebase:

- `ENSP00000288602`

## Fastest Path

1. open [`test_files/tp73.project.gentle.json`](../../test_files/tp73.project.gentle.json)
2. open `File -> Protein Evidence...`
3. in the Ensembl section, fetch `ENSP00000288602`
4. import it as a first-class protein sequence
5. in `Reverse translate protein`, select that protein
6. choose a speed profile / speed mark and run reverse translation
7. inspect the provenance panel
8. confirm the created coding DNA opens and is represented in lineage

## Step-by-Step

### Step 1: Open the Tutorial Project

GUI:

1. `File -> Open Project...`
2. choose [`test_files/tp73.project.gentle.json`](../../test_files/tp73.project.gentle.json)

This gives us a stable project container for the imported protein and the
reverse-translated coding DNA.

### Step 2: Fetch One Ensembl Protein Entry

GUI:

1. `File -> Protein Evidence...`
2. in the Ensembl section, enter `ENSP00000288602`
3. click `Fetch Ensembl`

Expected result:

- the recent Ensembl table gains one entry
- the selected-entry panel shows transcript/gene/species context

### Step 3: Import the Protein Sequence

GUI:

1. keep `entry_id = ENSP00000288602`
2. optionally set `output_id = ensp00000288602_protein`
3. click `Import Sequence`

Expected result:

- one first-class protein sequence is created in the project
- that sequence becomes available in the `Reverse translate protein` dropdown

### Step 4: Reverse Translate the Protein

In the same `Protein Evidence...` window:

1. in `Reverse translate protein`, choose the imported protein sequence
2. optional settings for a meaningful manual check:
   - `output_id = ensp00000288602_coding`
   - `speed profile = Human`
   - `speed mark = Slow`
   - leave `translation table` empty unless you want to test an explicit override
   - `target anneal Tm = 58.0`
   - `window bp = 9`
3. click `Reverse Translate`

Expected result:

- one coding DNA sequence is created and opened
- the result panel below the controls refreshes immediately

### Step 5: Inspect the Provenance Panel

Stay in `Protein Evidence...` and inspect the reverse-translation result block.

Manual checks:

- `output` shows the created coding-sequence id
- `length` shows the protein-to-coding length relationship
- `translation table` shows:
  - resolved table
  - label
  - source
  - organism/organelle context when available
- `speed profile` shows:
  - requested profile
  - resolved profile
  - source
  - reference species
- `speed mark` and `anneal heuristic` reflect the options you chose
- inline `Coding DNA` text is present

### Step 6: Inspect the Created Coding Sequence

Expected result after reverse translation:

- a new sequence window opens on the created coding DNA
- this is an ordinary sequence window, not only transient dialog state

Manual checks:

- the product is DNA, not protein
- the created sequence id matches the reverse-translation result panel
- the sequence remains available after closing and reopening the specialist

### Step 7: Verify Lineage Recording

GUI:

1. return to the main project window
2. open the lineage `Table` view if needed
3. find the new reverse-translation analysis row
4. use `Open Coding Sequence`

Expected result:

- the reverse-translation analysis appears in lineage as its own artifact row
- `Open Coding Sequence` reopens the same created coding DNA product

This is the key audit check: the reverse translation should not live only as
one ephemeral dialog result.

## What This Tutorial Is Actually Checking

This walkthrough is the practical sanity pass for:

- imported first-class proteins
- reverse translation into coding DNA
- visible provenance for translation-table and translation-speed resolution
- lineage persistence of the reverse-translation artifact

## Command Equivalent (After GUI)

The GUI path is the main point here, but the same engine route is shared.

If you want a follow-up parity check later, use the same operation through
`gentle_cli op ...` once dedicated shell sugar lands.

## Checkpoints

- Ensembl fetch succeeds for `ENSP00000288602`.
- Importing that entry creates one first-class protein sequence.
- Reverse translation creates one coding DNA sequence and opens it.
- The provenance panel states the resolved translation-table and
  speed-profile/source/reference-species story clearly.
- Lineage shows a reverse-translation analysis row that can reopen the created
  coding sequence.

## If Something Looks Wrong

Report these separately because they indicate different failure classes:

1. `Ensembl fetch/import fails`
   - likely online/provider or import-path problem
2. `Reverse translation creates DNA but provenance panel is incomplete`
   - likely GUI/report presentation drift
3. `Reverse translation looks fine in the dialog but not in lineage`
   - likely persisted-report or lineage-materialization regression
4. `Open Coding Sequence from lineage points to the wrong product`
   - likely analysis-artifact linkage regression

## Related Tutorial

If you want to validate the transcript-native expert path without depending on
online protein fetch/import, use:

- [`docs/tutorial/protein_transcript_native_expert_gui.md`](./protein_transcript_native_expert_gui.md)
