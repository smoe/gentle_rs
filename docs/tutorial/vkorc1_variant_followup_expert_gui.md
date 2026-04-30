# Promoter Design for VKORC1 / rs9923231

> Type: `GUI walkthrough`
> Status: `manual/hybrid`
> Goal: use the dedicated `Promoter design` window to go from one dbSNP
> marker to one exportable promoter-reporter handoff bundle for later work in
> human cells.

This is the shortest GUI path for the current ClawBio-facing story.

It starts from one pharmacogenomic alert around warfarin and `VKORC1`, then
uses the dedicated `Promoter design` window to:

- classify `rs9923231` as promoter-proximal
- propose one transcript-aware promoter fragment
- materialize matched reference and alternate inserts
- preview a mammalian luciferase reporter pair
- export one portable handoff bundle

If you want the longer explanation with matching CLI commands and more biology
background, use:

- [docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md](./vkorc1_warfarin_promoter_luciferase_gui.md)

## What This Tutorial Produces

The endpoint is not a wet-lab claim. It is one reproducible design package.

By the end you should have:

- one promoter-context summary for `rs9923231`
- one recommended `VKORC1` promoter fragment
- one reference insert
- one alternate insert
- one matched mammalian luciferase reporter pair
- one export bundle containing:
  - promoter-context JSON
  - promoter-candidate JSON
  - promoter-context SVG
  - reference reporter SVG
  - alternate reporter SVG
  - `report.md`
  - `result.json`
  - `commands.sh`

## Biological Framing

Keep the story narrow and honest:

- ClawBio motivates the follow-up:
  `rs9923231` is associated with warfarin sensitivity and sits upstream of
  `VKORC1`, so a regulatory follow-up is plausible.
- GENtle builds the deterministic follow-up:
  fetch the locus, derive promoter windows, suggest the fragment, materialize
  allele-matched inserts, preview the reporter pair, and export artifacts.

This is a human-cell reporter-design story, not a bacterial-expression story.

## Prerequisites

1. GENtle desktop app is running.
2. `Human GRCh38 Ensembl 116` is already prepared, or can be prepared from
   `File -> Prepare Reference Genome...`
3. dbSNP fetch is reachable for `rs9923231`.
4. The pinned local backbone is present:
   [data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb](/Users/u005069/.codex/worktrees/47dd/gentle_rs/data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb)

## Step 1: Fetch the Local SNP Context

1. Open `File -> Fetch GenBank / dbSNP...`
2. In the dbSNP part of the dialog:
   - rsID = `rs9923231`
   - genome = `Human GRCh38 Ensembl 116`
   - flank = `3000`
   - output id = `vkorc1_rs9923231_context`
3. Click `Fetch Region`

Expected result:

- one DNA window opens for `vkorc1_rs9923231_context`
- the view includes a visible `variation` marker for `rs9923231`
- the feature tree shows at least `variation`, `gene`, and `mRNA`

## Step 2: Open the Dedicated Promoter Design Window

1. In the context sequence window, click the `variation` marker for
   `rs9923231`
2. Open `Promoter design`
   You can reach it from:
   - the description pane
   - the feature-tree context menu
   - the map context menu

Expected default seed values:

- `source_seq_id = vkorc1_rs9923231_context`
- `variant_label_or_id = rs9923231`
- `gene_label = VKORC1`
- promoter defaults `1000 / 200`

If those are already filled, the window is seeded correctly.

## Step 3: Annotate Promoter Windows

1. In `Promoter design`, click `Annotate promoter windows`

What GENtle is doing here:

- it derives promoter windows from transcript TSS geometry
- it writes them as generated `promoter` features
- it keeps them distinct from source-imported annotations

What to sanity-check in the source view afterward:

- `VKORC1` is reverse-strand
- the promoter windows appear as generated annotations, not as imported truth

## Step 4: Summarize Promoter Context

1. Click `Summarize promoter context`

What you want to see:

- the SNP is classified as promoter-overlapping or promoter-candidate context
- the chosen gene is `VKORC1`
- the signed TSS distance is plausible for reverse-strand geometry
- suggested assay ids include
  `allele_paired_promoter_luciferase_reporter`

This step is the real handoff pivot:
GENtle is now saying the SNP is worth a promoter-reporter follow-up rather than
just showing a nearby locus.

## Step 5: Propose the Reporter Fragment

1. Keep the current defaults:
   - `retain_downstream_from_tss_bp = 200`
   - `retain_upstream_beyond_variant_bp = 500`
   - `max_candidates = 5`
2. Click `Propose reporter fragment`

What the recommended candidate should achieve:

- keep a little sequence past the selected TSS
- keep the SNP well inside the insert
- extend beyond the SNP on the biologically upstream side

For the current baseline, the top candidate should usually map to:

- local `2412..3501` in the `vkorc1_rs9923231_context` parent sequence

## Step 6: Extract the Recommended Fragment

1. Click `Extract recommended fragment`

Expected result:

- one extracted promoter fragment appears under the suggested output id
- the source context sequence remains the anchor for the expert window

This is important: the expert should stay attached to the original context,
even after derived fragments and reporter products appear.

## Step 7: Make the Reference and Alternate Inserts

1. Click `Make reference/alternate inserts`

Expected products:

- reference insert
- alternate insert

Both inserts should:

- share identical boundaries
- differ only at the SNP position

If either insert is missing, stop here and check that the selected fragment
still carries the variant feature.

## Step 8: Preview the Mammalian Reporter Pair

1. Click `Preview luciferase pair`

What GENtle does here:

- loads the pinned local promoterless mammalian luciferase backbone
- combines the reference insert with the same backbone
- combines the alternate insert with the same backbone
- opens or records the paired preview ids

Why this backbone is appropriate:

- it is a mammalian reporter baseline
- luciferase output is the readout for promoter activity
- the comparison stays focused on allele-dependent promoter output in human
  cells rather than on bacterial protein expression

## Step 9: Export the Handoff Bundle

1. Click `Export handoff bundle`
2. Choose a parent folder

GENtle creates one bundle directory named from the current gene + variant.

The bundle should contain:

- `variant_promoter_context.json`
- `promoter_reporter_candidates.json`
- one promoter-context SVG
- one reference reporter SVG
- one alternate reporter SVG
- `report.md`
- `result.json`
- `commands.sh`

## Step 10: Review the Bundle

Use these files first:

- `report.md`
  - bench-facing summary
  - explicit assumptions
  - next actions
- `result.json`
  - machine-readable design summary
  - stable ids and paths
- promoter-context SVG
  - quick promoter geometry review
- paired reporter SVGs
  - quick construct review

## What a Bench-facing Colleague Should Understand

By the end of this tutorial, the handoff should communicate:

- which SNP is being followed up
- which gene/transcript context was chosen
- why the fragment is treated as promoter-relevant
- what the reference and alternate constructs are
- which mammalian backbone was used
- which files reproduce the current design decision

That is enough to continue into practical construct review and later assay
planning without pretending the wet-lab outcome is already known.

## Troubleshooting

If `Promoter design` does not appear:

- make sure you selected a `variation` feature, not a gene or transcript row

If `Summarize promoter context` looks empty:

- confirm the context sequence still has `variation`, `gene`, and `mRNA`
  features visible

If `Preview luciferase pair` fails:

- confirm the local backbone file exists
- confirm both allele inserts were materialized successfully

If `Export handoff bundle` writes files but the contents look incomplete:

- rerun `Summarize promoter context`
- rerun `Propose reporter fragment`
- rerun `Preview luciferase pair`
- then export again
