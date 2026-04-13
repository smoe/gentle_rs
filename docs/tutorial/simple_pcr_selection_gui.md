# Simple PCR From a Selected Core Region

> Type: `GUI walkthrough`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is tied to the current
> selection-first PCR Designer and the shared primer-design engine route.

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

Any local sequence works, but `tp73.ncbi.gb` is a stable bundled example with
enough annotated context to make the selection easy.

## Fastest Path

1. open [`test_files/tp73.ncbi.gb`](../../test_files/tp73.ncbi.gb)
2. keep the map in linear mode
3. drag-select one short core region of interest
4. right-click that selection
5. choose `Simple PCR from selection`
6. in `PCR Designer`, adjust:
   - `max primer distance from core`
   - `max amplicon`
7. click `Design Primer Pairs`

## Step-by-Step

### Step 1: Open One Sequence

GUI:

1. `File -> Open Sequence...`
2. choose [`test_files/tp73.ncbi.gb`](../../test_files/tp73.ncbi.gb)

### Step 2: Select the Core Region

GUI:

1. switch the DNA view to `Linear` if needed
2. drag over the exact region that must be present in the PCR product

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
2. optionally reduce `max pairs` if you want a shorter result list
3. click `Design Primer Pairs`

### Step 6: Review the Result

Inspect the returned primer-pair report for:

- forward primer
- reverse primer
- amplicon length
- whether the ROI is flanked
- Tm / GC values
- any advisory notes

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
  [`docs/tutorial/generated/chapters/13_pcr_selection_batch_primer_pairs_offline.md`](./generated/chapters/13_pcr_selection_batch_primer_pairs_offline.md)
- If you want GUI navigation help for the wider application:
  [`docs/tutorial/README.md`](./README.md)
