# Transcript-Native Protein Expert Sanity Check

> Type: `GUI walkthrough`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> current `Protein Evidence...` specialist and the shared transcript-first
> Protein Expert route.

This tutorial is meant as a **manual check** after the recent transcript-native
protein work.

It verifies the path that does **not** depend on UniProt or Ensembl as the
source of truth:

- open one local project sequence
- derive transcript/product geometry on demand
- inspect the shared Protein Expert
- export the same view as SVG through the shared engine route

## What You Will Check

By the end of this walkthrough, you should be able to confirm that:

- `Protein Evidence...` can open the derived Protein Expert without any stored
  external protein evidence
- transcript-native product rows appear and stay readable
- the details grid exposes translation table, translation-speed provenance, and
  derivation mode
- derived-only SVG export succeeds from the same expert payload

## Local Input

Use the bundled local project fixture:

- [`test_files/tp73.project.gentle.json`](../../test_files/tp73.project.gentle.json)

Why this file:

- it is already a project, not only a bare sequence file
- it contains one stable sequence id, `tp73.ncbi`
- it avoids online dependencies for the transcript-native check

Fixture note:

- provenance for this historical TP73 fixture is tracked in
  [`test_files/README.md`](../../test_files/README.md)

## Fastest Path

1. open [`test_files/tp73.project.gentle.json`](../../test_files/tp73.project.gentle.json)
2. open `File -> Protein Evidence...`
3. in `Project entry to sequence`, choose `seq_id = tp73.ncbi`
4. leave `transcript` empty for the first pass
5. click `Open Derived Protein Expert`
6. confirm the expert opens with transcript-native product rows
7. return to `Protein Evidence...`
8. click `Render Derived Protein SVG...`
9. save the SVG and confirm export succeeds

## Step-by-Step

### Step 1: Open the Tutorial Project

GUI:

1. `File -> Open Project...`
2. choose [`test_files/tp73.project.gentle.json`](../../test_files/tp73.project.gentle.json)

Expected result:

- the project opens with sequence `tp73.ncbi`
- the project/lineage window is available for later inspection if needed

### Step 2: Open Protein Evidence

GUI:

1. `File -> Protein Evidence...`
2. in `Project entry to sequence`, choose `tp73.ncbi`
3. leave `transcript` empty for the first pass

Why this matters:

- this path checks transcript-native derivation directly from the current
  sequence state
- it does not require a stored UniProt projection or Ensembl entry

### Step 3: Open the Derived Protein Expert

GUI:

1. click `Open Derived Protein Expert`

Expected result:

- one Protein Expert window opens
- the lower product section is populated from transcript-native derivation
- no external-protein fetch/import step is required

Manual checks:

- the view should not be empty
- transcript/product rows should appear without needing a UniProt projection id
- the details grid should expose fields such as:
  - derived protein length
  - translation table / source
  - translation-speed profile / source / reference species
  - derivation mode

### Step 4: Narrow to One Transcript (Optional but Useful)

Back in `Protein Evidence...`:

1. fill the `transcript` field with one transcript id visible in the expert
   details or transcript list
2. click `Open Derived Protein Expert` again

Expected result:

- the same expert route opens again, now narrowed to one transcript
- this helps verify that transcript filtering still affects the transcript-first
  expert deterministically

### Step 5: Export the Derived Protein SVG

Back in `Protein Evidence...`:

1. click `Render Derived Protein SVG...`
2. save the SVG somewhere convenient, for example `/tmp/tp73_derived_protein.svg`

Expected result:

- export completes without needing UniProt or Ensembl evidence
- the SVG reflects the same transcript-native expert view rather than a
  separate rendering model

## What This Tutorial Is Actually Checking

This tutorial checks the **derived-only** protein path that should remain valid
when external evidence is absent or intentionally ignored.

It is the right manual check today for:

- transcript-native CDS/protein derivation in the GUI
- transcript filtering into the Protein Expert
- derived-only Protein Expert SVG export

It does **not** yet check first-class protein-sequence materialization from a
GUI button, because that direct GUI action is not the main user-facing path at
this moment.

## Command Equivalent (After GUI)

The same expert/export surface is shared with CLI/shell:

```bash
cargo run --bin gentle_cli -- inspect-feature-expert tp73.ncbi protein-comparison
cargo run --bin gentle_cli -- render-feature-expert-svg tp73.ncbi protein-comparison /tmp/tp73_derived_protein.svg
```

Use these after the GUI pass if you want to confirm adapter parity.

## Checkpoints

- `Protein Evidence...` opens and accepts `tp73.ncbi` as the target sequence.
- `Open Derived Protein Expert` opens a non-empty transcript-native protein view.
- The expert details expose translation/provenance fields instead of hiding
  them inside adapter-only status text.
- `Render Derived Protein SVG...` succeeds and produces one file from the same
  shared expert route.

## If Something Looks Wrong

Report these separately because they point at different layers:

1. `Derived Protein Expert opens but is empty`
   - likely transcript/admission or expert-payload regression
2. `Details grid misses translation/speed fields`
   - likely GUI presentation drift
3. `SVG export fails but the expert opens`
   - likely export-route mismatch rather than derivation failure

## Related Next Step

To check first-class protein import plus reverse translation and lineage audit,
continue with:

- [`docs/tutorial/protein_reverse_translation_gui.md`](./protein_reverse_translation_gui.md)
