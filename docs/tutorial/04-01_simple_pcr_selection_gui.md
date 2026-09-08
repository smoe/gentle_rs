# Simple PCR From a Selected Core Region

> Type: `GUI walkthrough`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is tied to the current
> selection-first PCR Designer and the shared primer-design engine route.

See also: executable reference chapter
[`18 Simple PCR From a Selected Core Region`](./generated/chapters/04-01_simple_pcr_selection_gui.md).

This tutorial is intentionally small.

It is for the most common first PCR question:

1. what is the core region I need to include,
2. how far away may the primers sit from that core, and
3. how long may the product become?

GENtle now supports that flow directly from a selection context menu, so you do
not have to translate the selection into PCR form fields by hand before you
start.

## What You Will Do

By the end of this tutorial, you should be able to:

- select one core region on a sequence map
- start a simple PCR directly from that selection
- set a maximum primer distance from the core region
- set a maximum amplicon length
- run primer-pair design and inspect the returned candidates

## Local Input

Use one local committed sequence so the tutorial works offline:

- [`test_files/tp73.ncbi.gb`](../../test_files/tp73.ncbi.gb)

The tutorial workflow extracts bases [61520, 62320) (0-based, end-exclusive)
from that file into the 800-base sequence `tp73_locus`. It retains the full
sequence as `simple_pcr_source_locus` for provenance. Open the compact extract
for this walkthrough, not the full source. The scripted PCR oracle uses the
same bases under the name `simple_pcr_template`.

The fixed smoke selection is `=201 .. 600` in the selection formula field
(1-based inclusive), which becomes engine ROI [200, 600). This leaves 200 bases
on each side for primer search. It is a real TP73 interval used for teaching,
not a validated assay or evidence of whole-genome specificity.

## Fastest Path

1. open `File -> Open Tutorial Project... -> Core -> 18. Simple PCR From a Selected Core Region`
2. open the 800-base `tp73_locus` extract and keep its map in linear mode
3. enter `=201 .. 600` in the selection formula field and apply it, or drag-select that core
4. right-click that selection
5. choose `Simple PCR from selection`
6. in `PCR Designer`, adjust:
   - `max primer distance from core`
   - `max amplicon`
7. click `Design Primer Pairs`

For the scripted starter, run the workflow linked by the executable reference
chapter. Opening the full GenBank locus directly is a separate, larger-input
exercise; its runtime is not covered by the bounded beginner smoke.

## Step-by-Step

### Step 1: Open One Sequence

GUI:

1. open the tutorial project through the menu above
2. open the compact `tp73_locus` sequence (800 bases)

### Step 2: Select the Core Region

GUI:

1. switch the DNA view to `Linear` if needed
2. apply `=201 .. 600` in the selection formula field, or drag over those bases

This selection is your **core ROI**.

Keep it simple:

- do not try to select the whole desired amplicon
- select only the region that must definitely be covered

### Step 3: Start Simple PCR From the Selection

GUI:

1. right-click on the map while the selection is active
2. choose `Simple PCR from selection`

What GENtle does for you:

- copies the selection into the PCR ROI fields
- enables `require ROI flanking`
- sets `min amplicon` to at least the core length
- applies default forward/reverse flank windows using the current
  `max primer distance from core`
- opens or focuses the dedicated `PCR Designer`

### Step 4: Interpret the Three Main Inputs

Inside `PCR Designer`, the beginner path is now the `Simple PCR starter` block.

Use it like this:

- `core ROI`
  - already seeded from your map selection
- `max primer distance from core`
  - this defines how wide the primer-search windows are on the left and right
    side of the core ROI
  - after you change it, click `Apply simple flank windows`
- `max amplicon`
  - this is the longest allowed product
  - set it in the normal primer-pair form just below the starter block

Important translation:

- GENtle turns `max primer distance from core = D` into:
  - forward primer search window: `[ROI start - D .. ROI start]`
  - reverse primer search window: `[ROI end .. ROI end + D]`

So the simple controls are still deterministic and inspectable: they just write
the existing forward/reverse side-window fields for you.

### Step 5: Run Primer Design

GUI:

1. keep or adjust `max amplicon`
2. set `max pairs` to `5` for this walkthrough
3. click `Design Primer Pairs`

`max pairs` limits the returned report, not the search effort. The 800-base
template and its 200-base flanks bound the search here. The Linux acceptance
runner still requires a non-empty primer report within its ten-minute compute
budget; a timeout is a failure, not a skipped or successful tutorial.

### Step 6: Review the Result

Inspect the returned primer-pair report for:

- forward primer
- reverse primer
- amplicon length
- left/right distance from the core ROI
- whether the ROI is flanked
- Tm / GC values
- any advisory notes

GENtle now keeps that beginner wording visible in two places:

- the in-panel `Primer report preview` inside `PCR Designer`
- shared-shell `primers show-report REPORT_ID` as `simple_pcr_pairs`

For this beginner flow, the most important question is simply:

> Does this primer pair cleanly flank my selected core region without exceeding
> my maximum product length?

## What To Adjust First If No Good Pair Appears

Try these in order:

1. increase `max primer distance from core`
2. increase `max amplicon`
3. widen primer length/Tm ranges only after the first two do not help

That keeps the beginner workflow intuitive before moving into advanced primer
constraints.

## Why This Tutorial Exists

The full PCR Designer also supports:

- painted ROI/upstream/downstream intervals
- PCR batch queues
- richer side constraints
- qPCR assay design

Those are powerful, but they are not the first story most people want.

This tutorial keeps the initial mental model deliberately narrow:

- select the biology you care about
- define how far primers may move away from it
- define the longest acceptable product

## Related Next Steps

- If you want the richer queue-and-paint workflow, continue with:
  [`docs/tutorial/generated/chapters/04-02_pcr_selection_batch_primer_pairs_offline.md`](./generated/chapters/04-02_pcr_selection_batch_primer_pairs_offline.md)
- If you want GUI navigation help for the wider application:
  [`docs/tutorial/README.md`](./README.md)

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context copied from GENtle Help -> Tutorial -> Copy Feedback Context.

- Tutorial title:
- Tutorial/chapter id:
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
