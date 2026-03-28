# Gibson Arrangements Tutorial

> Type: `GUI walkthrough + CLI parity`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> shared `gibson apply`, arrangement persistence, and `render-pool-gel-svg`
> paths so the same setup can be replayed deterministically.

This tutorial shows how GENtle turns one successful Gibson apply step into a
reusable arrangement.

In current GENtle terms, an arrangement is a stored serial lane setup. It is
not yet the future physical rack/plate placement layer. The arrangement is the
semantic lane order that can be reopened, inspected, and exported as one gel.

For a single-insert Gibson run, the intended story is:

1. original vector,
2. insert,
3. assembled product,
4. with DNA ladders flanking the whole run on gel export.

## Fastest Setup

Recommended route:

1. `File -> Open Tutorial Project...`
2. choose `Gibson Arrangements Starter Project (offline)`
3. the Help window should open automatically on this arrangement guide
4. continue below from `Step 2`

The same starter baseline is described here:

- [`docs/tutorial/generated/chapters/15_gibson_specialist_testing_baseline.md`](./generated/chapters/15_gibson_specialist_testing_baseline.md)
- [`docs/tutorial/generated/chapters/16_gibson_arrangements_baseline.md`](./generated/chapters/16_gibson_arrangements_baseline.md)

## What You Will Test

By the end of this tutorial, you should have verified all of these:

- one Gibson apply step creates three new sequence outputs:
  - left insert primer
  - right insert primer
  - assembled product
- those three outputs land as separate singleton containers, not one pooled
  container
- the same Gibson apply step also records one reusable serial arrangement
- the arrangement lane order is:
  - original vector
  - ordered insert lane(s)
  - assembled product
- arrangement-level gel export adds DNA ladders to either side of the sample
  lanes
- the same result can be replayed from the CLI through `gibson apply` and
  `render-pool-gel-svg --arrangement`

## Inputs

Use the same deterministic local inputs as the Gibson specialist tutorial:

- destination vector:
  [`test_files/pGEX-3X.gb`](../../test_files/pGEX-3X.gb)
- synthetic insert:
  [`docs/tutorial/inputs/gibson_insert_demo.fa`](./inputs/gibson_insert_demo.fa)
- deterministic single-insert Gibson plan for CLI parity:
  [`docs/figures/gibson_single_insert_readme.plan.json`](../figures/gibson_single_insert_readme.plan.json)

## Step 1: Open the Starter Project

If you already completed the Gibson specialist tutorial and still have the
project open, you can skip to `Step 2`.

GUI:

1. `File -> Open Tutorial Project...`
2. choose `Gibson Arrangements Starter Project (offline)`

What to verify:

- `gibson_destination_pgex` is present and circular
- `gibson_insert_demo` is present and linear
- the Help window opened on a Gibson tutorial page

## Step 2: Apply One Gibson Cloning Run

Use the same stable single-insert route that the Gibson specialist tutorial
already exercises.

GUI:

1. `Patterns -> Gibson...`
2. in `Destination`, choose `gibson_destination_pgex`
3. keep `opening = defined cut/opening`
4. click `Use` on `SmaI (MCS)` if it is offered
   - fallback: set `left cut edge = 941` and `right cut edge = 941`
5. in `Insert`, choose `gibson_insert_demo`
6. keep `orientation = forward`
7. in `Design Targets`, use the deterministic tutorial values:
   - `overlap min/max bp = 30 / 30`
   - `minimum overlap Tm C = 60`
   - `priming Tm window C = 58 / 68`
   - `priming length bp = 18 / 35`
   - `output id hint = gibson_ui_test_product`
8. click `Preview Gibson Plan`
9. confirm `can_execute=true`
10. click `Apply Gibson Cloning`

What to verify immediately:

- the status area says that Gibson cloning was applied
- the status area also mentions a created Gibson serial arrangement
- closing the Gibson window is fine; the arrangement is now stored in project state

## Step 3: Inspect the New Containers

Return to the main project window.

GUI:

1. stay on the lineage page
2. scroll to `Containers`

What to verify:

- you now have separate singleton containers for:
  - the left insert primer
  - the right insert primer
  - the assembled product
- these are separate vials conceptually; they are not collapsed into one pool

Why this matters:

- the two primers are synthesized separately
- the assembled product is a different physical artifact again
- the arrangement will reference those exact container objects when it builds
  the lane order

## Step 4: Inspect the New Arrangement

GUI:

1. stay on the main project window
2. scroll to `Arrangements`

What to verify in the new row:

- `Mode = serial`
- `Lanes = 3` for this single-insert case
- `Ladders` is either `auto` or one/two named DNA ladders
- the arrangement name should be similar to:
  - `Gibson lane set: gibson_ui_test_product`
  - or a deterministic product-derived fallback if you used a different output hint

How to read that row:

- `Lanes = 3` counts only the sample lanes:
  - original vector
  - insert
  - assembled product
- the ladders are not counted as sample lanes; they are added at gel-render time
  as flanking reference lanes

## Step 5: Open the Lane Containers

GUI:

1. in the arrangement row, click `Open Lanes`

What to verify:

- GENtle opens the lane containers referenced by the arrangement
- the set should correspond to the Gibson story:
  - original vector
  - insert
  - assembled product

This is a quick sanity check that the arrangement is pointing at the expected
containers instead of some unrelated latest pool.

## Step 6: Export the Arrangement Gel

GUI:

1. in the same arrangement row, click `Export Gel`
2. save the SVG somewhere easy to inspect

What to verify in the SVG:

- sample lanes appear in the arrangement order:
  - vector
  - insert
  - assembled product
- DNA ladders appear on both sides of the sample set
- the export is one arrangement-level gel, not three unrelated single-lane
  exports

This is the practical payoff of storing the arrangement: you no longer have to
reconstruct the lane order by hand each time.

## Step 7: Inspect the Graph View

GUI:

1. switch the main lineage page to `Graph`
2. inspect the Gibson operation hub and the arrangement node

What to verify:

- one `Gibson cloning` hub sits between the vector/insert inputs and the three
  created outputs
- the arrangement is represented separately as an arrangement node
- clicking the `Gibson cloning` node reopens the Gibson specialist

Interpretation:

- the Gibson hub is the cloning operation
- the arrangement is a reusable downstream plan for lane ordering and gel export
- keeping them separate is useful because the same output sequences could later
  participate in more than one arrangement

## CLI Parity

You can replay the same arrangement-creation path from the command line with a
fresh temporary state file.

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_arrangement_tutorial.state.json \
  workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json

cargo run --bin gentle_cli -- \
  --state /tmp/gibson_arrangement_tutorial.state.json \
  gibson apply @docs/figures/gibson_single_insert_readme.plan.json
```

At that point, the fresh tutorial state should contain one new arrangement,
which is typically `arrangement-1` in a clean state file. Export its gel with:

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_arrangement_tutorial.state.json \
  render-pool-gel-svg - /tmp/gibson_arrangement_tutorial.svg --arrangement arrangement-1
```

Optional audit step:

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_arrangement_tutorial.state.json \
  export-state /tmp/gibson_arrangement_tutorial.export.json
```

What to verify in the exported state:

- `container_state.arrangements` contains one serial arrangement
- the arrangement lane container IDs correspond to:
  - destination vector container
  - insert container
  - assembled product container
- the arrangement also carries ladder recommendations

## Advanced: Manual Arrangement Creation

You do not have to rely on Gibson auto-arrangement creation.

If you already know which containers belong in one lane set, you can create a
serial arrangement explicitly:

```bash
cargo run --bin gentle_cli -- \
  arrange-serial container-1,container-2,container-5 \
  --id arrangement-demo \
  --name "Gibson run A" \
  --ladders "NEB 100bp DNA Ladder,NEB 1kb DNA Ladder"
```

Then export the arrangement later with:

```bash
cargo run --bin gentle_cli -- \
  render-pool-gel-svg - /tmp/gibson_arrangement_manual.svg --arrangement arrangement-demo
```

This is the more general route for future mixed-experiment racks, reused lane
plans, or explicitly curated gel batches.

## Current Limits and Meaning

Important current meaning of `arrangement` in GENtle:

- it is a reusable serial lane setup
- it is not yet the future rack/plate placement UI
- it already gives us the right semantic layer for:
  - gel export
  - lane reuse
  - later rack mapping
  - later label printing

So if you are thinking ahead to racks, plates, and bench-day planning, this is
already the right conceptual starting point. It just is not yet the full
physical placement layer.

## Expected Outcome

If everything worked, you now have:

- one Gibson operation with two inputs and three outputs
- three concrete singleton output containers
- one reusable serial arrangement referencing the vector, insert, and assembled
  product
- one arrangement-level gel export with ladders flanking the sample lanes

That is the intended current workflow for “getting the arrangement done” in
GENtle.
