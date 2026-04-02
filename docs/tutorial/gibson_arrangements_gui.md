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
4. continue below from `Step 1`

The same starter baseline is described here:

- [`docs/tutorial/generated/chapters/15_gibson_specialist_testing_baseline.md`](./generated/chapters/15_gibson_specialist_testing_baseline.md)
- [`docs/tutorial/generated/chapters/16_gibson_arrangements_baseline.md`](./generated/chapters/16_gibson_arrangements_baseline.md)

Difference from the Gibson specialist tutorial:

- chapter 15 opens with only the destination and insert so you can test the
  Gibson specialist itself
- chapter 16 opens after one deterministic Gibson apply has already been
  executed, so arrangement review and gel export are disjunct from the cloning
  tutorial

## What You Will Test

By the end of this tutorial, you should have verified all of these:

- the arrangement starter already contains the result of one deterministic
  Gibson apply step, including three new sequence outputs:
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
- the same result can be replayed from the CLI through the dedicated starter
  workflow and `render-pool-gel-svg --arrangement`

## Inputs

The arrangement starter is built from the same deterministic local inputs as
the Gibson specialist tutorial, plus the canonical single-insert Gibson plan:

- destination vector:
  [`test_files/pGEX-3X.gb`](../../test_files/pGEX-3X.gb)
- synthetic insert:
  [`docs/tutorial/inputs/gibson_insert_demo.fa`](./inputs/gibson_insert_demo.fa)
- deterministic single-insert Gibson plan for CLI parity:
  [`docs/figures/gibson_single_insert_readme.plan.json`](../figures/gibson_single_insert_readme.plan.json)

This arrangement tutorial starts after that plan has already been applied once.
So the product and arrangement should be present as soon as the tutorial
project opens.

## Step 1: Open the Starter Project

GUI:

1. `File -> Open Tutorial Project...`
2. choose `Gibson Arrangements Starter Project (offline)`

What to verify:

- `gibson_destination_pgex` is present and circular
- `gibson_insert_demo` is present and linear
- `gibson_destination_pgex_with_gibson_insert_demo` is already present as the
  assembled product
- the Help window opened on `Gibson Arrangements Tutorial`

If you want to reproduce the cloning step itself rather than inspect its
result, switch to the Gibson specialist tutorial instead:

- [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md)

## Step 2: Inspect the New Containers

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

## Step 3: Inspect the New Arrangement

GUI:

1. stay on the main project window
2. scroll to `Arrangements`

What to verify in the new row:

- `Mode = serial`
- `Lanes = 3` for this single-insert case
- `Ladders` is either `auto` or one/two named DNA ladders
- the arrangement name should be similar to:
  - `Gibson lane set: gibson_destination_pgex_with_gibson_insert_demo`
  - or a deterministic product-derived fallback if the internal naming path
    changes later

How to read that row:

- `Lanes = 3` counts only the sample lanes:
  - original vector
  - insert
  - assembled product
- the ladders are not counted as sample lanes; they are added at gel-render time
  as flanking reference lanes

## Step 4: Open the Lane Containers

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

## Step 5: Export the Arrangement Gel

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

## Step 6: Inspect the Graph View

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

You can recreate the same arrangement-ready project from the command line with
a fresh temporary state file in one workflow replay.

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_arrangement_tutorial.state.json \
  workflow @docs/examples/workflows/gibson_arrangements_baseline.json
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
