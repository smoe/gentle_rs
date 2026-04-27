# qPCR Across Exon Junctions

> Type: `GUI walkthrough + shared-shell parity`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it follows the current
> Splicing Expert -> PCR Designer qPCR route and the shared
> `DesignQpcrAssays` engine operation.

This tutorial is for designing qPCR assays from transcript annotation rather
than from a manually typed genomic interval.

The important distinction is:

- **shared-transcript qPCR** asks for an assay that reports the selected
  gene/group across its transcript set
- **transcript-specific qPCR** asks for an assay that distinguishes one
  selected transcript, optionally requiring exon-junction evidence

For strict exon-junction work, use transcript-specific qPCR with
`Specificity evidence = Junction only`. GENtle currently enforces this as a
junction-spanning primer requirement and reports whether the retained forward,
reverse, probe, and amplicon context span or cover annotated junctions.

## What You Will Do

By the end of this tutorial, you should be able to:

- open a transcript-rich local TP73 example
- open Splicing Expert for the TP73 gene/group
- launch shared-transcript qPCR directly from Splicing Expert
- launch transcript-specific qPCR from one selected transcript lane
- require junction-only evidence for transcript-specific assays
- inspect the saved qPCR report for transcript support and junction context
- export the report and the qPCR protocol cartoon

## Local Input

Use the bundled TP73 locus:

- [`test_files/tp73.ncbi.gb`](../../test_files/tp73.ncbi.gb)

This file is useful here because it has one large TP73 gene feature, multiple
TP73 mRNA/CDS features, and antisense transcript features nearby. The tutorial
does not require internet access.

## Mental Model

GENtle splits this workflow into two windows.

- `Splicing Expert`
  - owns the transcript group and the selected transcript lane
  - knows which exon chains and exon-exon junctions belong to each transcript
- `PCR Designer`
  - owns qPCR constraints, primer/probe limits, report ids, and execution
  - receives transcript-aware context from Splicing Expert

Do not start by manually typing exon coordinates unless you are deliberately
testing a low-level coordinate case. For biological exon-junction assays, seed
from Splicing Expert so the qPCR request carries transcript identity.

## Fastest Path

1. Open `File -> Open Tutorial Project... -> Core -> 18. Simple PCR From a Selected Core Region`.
2. Open the loaded `tp73_locus` sequence window.
3. In the feature tree, select the TP73 gene or one TP73 transcript feature.
4. Open `Splicing Expert`.
5. For total/group-level TP73 qPCR, click `Design shared-transcript qPCR`.
6. For one isoform/transcript, select a transcript lane in Splicing Expert,
   then click `Design transcript-specific qPCR`.
7. In `PCR Designer`, confirm the mode selector is on `qPCR`.
8. For transcript-specific exon-junction work, set
   `Specificity evidence` to `Junction only`.
9. Set practical qPCR bounds, for example:
   - `min amplicon = 70`
   - `max amplicon = 180`
   - `max assays = 50`
10. Click `Design qPCR Assays`.
11. Inspect the qPCR report preview for transcript targeting, supported
    transcripts, covered junctions, and whether any oligo spans a junction.

If you prefer to start from a file rather than the tutorial menu, use
`File -> Open Sequence...` and choose
[`test_files/tp73.ncbi.gb`](../../test_files/tp73.ncbi.gb).

## Step-by-Step

### Step 1: Open the TP73 Locus

GUI:

1. `File -> Open Sequence...`
2. choose [`test_files/tp73.ncbi.gb`](../../test_files/tp73.ncbi.gb)
3. keep the map in linear mode
4. make sure `Gene` and `mRNA` feature layers are visible

You should see TP73 as the main gene/group on the loaded locus.

### Step 2: Open Splicing Expert

GUI:

1. select the TP73 gene or a TP73 mRNA/transcript feature in the feature tree
2. use the feature detail action `Open Splicing Window`

What to check:

- the window title and locus summary mention TP73
- the transcript count is non-zero
- exon/junction rows are visible
- transcript lanes can be selected

If you selected a neighboring antisense feature by accident, choose the TP73
gene/group again before continuing.

### Step 3: Choose the qPCR Intent

Use one of these two paths.

#### A. Shared-Transcript qPCR

Use this when you want a gene/group-level qPCR assay that should be supported
by as much of the selected TP73 transcript group as possible.

GUI:

1. in Splicing Expert, click `Design shared-transcript qPCR`
2. let GENtle open or focus `PCR Designer`
3. in `PCR Designer`, confirm:
   - mode selector: `qPCR`
   - transcript targeting: `Shared across transcripts`
   - transcript context mentions TP73 and the current sequence

Interpretation:

- GENtle prefers exon or exon-chain context shared across the transcript group.
- If no fully shared context is assayable, the report records a deterministic
  fallback to the largest supported transcript subset.
- This is good for "measure this gene/group" questions, but it is not the
  strictest isoform-specific assay.

#### B. Transcript-Specific Junction qPCR

Use this when you want one selected transcript to be distinguishable from
other transcripts in the same TP73 group.

GUI:

1. in Splicing Expert, select the transcript lane you care about
2. click `Design transcript-specific qPCR`
3. let GENtle open or focus `PCR Designer`
4. in `PCR Designer`, confirm:
   - mode selector: `qPCR`
   - transcript targeting: `Specific transcript`
   - selected transcript is not `<select in Splicing Expert>`
5. set `Specificity evidence` to `Junction only`

Interpretation:

- `Junction only` asks GENtle to retain assays with a primer spanning a
  transcript-unique exon-exon junction.
- `Either (prefer junction)` is a fallback option when a useful transcript
  assay exists but strict junction-only evidence cannot be satisfied.
- `Unique exon / exon chain only` is useful for transcript discrimination, but
  it is not the strict exon-junction tutorial path.

### Step 4: Tighten qPCR Constraints

In `PCR Designer -> Design qPCR assays`, start with practical small-amplicon
settings:

- `min amplicon`: `70`
- `max amplicon`: `180`
- `max primer Tm delta`: keep the default unless the report is too sparse
- `max probe Tm delta`: keep the default unless the report is too sparse
- `max assays`: `50`
- `report_id`: use a memorable name, for example:
  - `tp73_shared_junction_qpcr`
  - `tp73_nm005427_junction_qpcr`

Leave the forward, reverse, and probe side constraints at their defaults for a
first pass. Only tighten primer/probe length, GC, Tm, motifs, or fixed
positions after you have seen whether the transcript-aware search can find a
basic assay.

### Step 5: Run qPCR Design

GUI:

1. click `Design qPCR Assays`
2. wait for progress messages to finish
3. click `Show report_id` if the preview did not update automatically

Expected result:

- a saved qPCR report appears under the chosen `report_id`
- the preview shows ranked assays with forward primer, reverse primer, probe,
  amplicon, and transcript context

### Step 6: Review Junction and Transcript Evidence

For a transcript-specific junction assay, inspect the saved qPCR report for:

- `transcript targeting`
  - should say `Specific transcript`
  - should name the selected transcript id
- `requested evidence`
  - should say `Junction only`
- `realized specificity evidence`
  - should also say `Junction only`
- per-assay transcript context
  - supported transcript ids
  - covered junction labels
  - whether the forward primer spans a junction
  - whether the reverse primer spans a junction
  - whether the probe crosses a junction
  - whether the assay satisfies the requested targeting intent

For a shared-transcript assay, inspect:

- support count and support fraction
- whether a shared-mode fallback was used
- whether the top assay is supported by the transcript set you intended to
  measure

The most important question is:

> Does the retained assay satisfy the biological intent, not just the primer
> length and Tm filters?

### Step 7: Export the Handoff Artifacts

GUI:

1. in `PCR Designer`, click `Export report_id...`
2. choose an output path such as:
   `exports/tp73_nm005427_junction_qpcr.report.json`
3. click `Export qPCR Cartoon SVG...`
4. choose an output path such as:
   `exports/tp73_nm005427_junction_qpcr.protocol.svg`

The JSON report is the important design handoff. The SVG cartoon is the
readable protocol strip for a notebook, README, or planning discussion.

## Shared-Shell Parity

The GUI path above is the recommended way to choose the transcript because it
lets you see the transcript lanes first. The same engine route is available
from the shared shell.

Start from the local tutorial state:

```bash
STATE=/tmp/gentle_qpcr_exon_junctions.state.json

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  workflow @docs/examples/workflows/simple_pcr_selection_gui.json
```

Then ask the shell to build a qPCR seed request from a splicing feature:

```bash
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  shell 'primers seed-qpcr-from-splicing tp73_locus FEATURE_ID --mode distinguish_transcript --transcript-id TRANSCRIPT_ID --specificity-evidence junction_only' \
  > /tmp/tp73_junction_qpcr_seed.json
```

Replace:

- `FEATURE_ID` with the GENtle feature index for the TP73 splicing group
- `TRANSCRIPT_ID` with the transcript id selected in Splicing Expert, for
  example an `NM_...` id from the TP73 transcript table

Run the design:

```bash
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  primers design-qpcr @/tmp/tp73_junction_qpcr_seed.json --backend internal
```

Inspect and export:

```bash
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  primers list-qpcr-reports

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  primers show-qpcr-report REPORT_ID

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  primers export-qpcr-report REPORT_ID /tmp/tp73_junction_qpcr.report.json

cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg pcr.assay.qpcr /tmp/tp73_junction_qpcr.protocol.svg
```

The shared-shell command emits the same transcript-aware
`DesignQpcrAssays` request that the GUI uses. That keeps GUI, CLI, ClawBio,
JS, and Lua on the same deterministic engine contract.

## Troubleshooting

### `Specific transcript` says no transcript is selected

Go back to Splicing Expert, select one transcript lane, then click
`Design transcript-specific qPCR` again. The selection lives in Splicing
Expert, not in the PCR Designer table.

### No assays are accepted in `Junction only` mode

Try this order:

1. increase `max amplicon` from `180` to `250`
2. increase `max assays`
3. relax primer/probe Tm or GC only after checking the first two
4. switch to `Either (prefer junction)` if a strict junction-only assay is not
   biologically available for that transcript

Do not silently reinterpret an empty strict result as a usable assay. Empty can
mean the requested junction-specific evidence does not exist under the current
constraints.

### Shared-transcript mode falls back to a subset

That means the entire transcript group did not share an assayable exon/exon
chain under the current settings. Either accept the reported subset, change
the biological question, or move to transcript-specific qPCR.

### You want to avoid genomic DNA amplification

Exon-junction placement helps, but still record standard assay controls in the
wet-lab plan:

- no-template control
- no-RT control for RNA-derived cDNA assays
- DNase treatment decision
- melt curve or probe-specific validation plan, depending on assay chemistry

GENtle designs and explains the candidate assay; it does not replace those
experimental controls.

## What To Save In Your Notes

For each qPCR assay you keep, record:

- sequence/project id
- selected gene/group
- selected transcript id if transcript-specific
- qPCR report id
- assay rank
- forward primer, reverse primer, and probe sequences
- amplicon coordinates and length
- covered junction labels
- realized specificity evidence
- exported JSON report path
- exported protocol-cartoon SVG path

That makes the design review reproducible before primers/probes are ordered.

## Related Next Steps

- For the beginner pair-PCR version of ROI seeding:
  [`docs/tutorial/simple_pcr_selection_gui.md`](./simple_pcr_selection_gui.md)
- For RNA-read mapping evidence around the same transcript/exon context:
  [`docs/tutorial/rna_read_batch_gene_support_cli.md`](./rna_read_batch_gene_support_cli.md)
- For the generated executable PCR batch chapter:
  [`docs/tutorial/generated/chapters/13_pcr_selection_batch_primer_pairs_offline.md`](./generated/chapters/13_pcr_selection_batch_primer_pairs_offline.md)
