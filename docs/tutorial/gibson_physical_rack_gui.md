# Gibson Physical Rack Tutorial

> Type: `GUI walkthrough + CLI parity`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> shared rack-placement, `racks isometric-svg`, `racks fabrication-svg`,
> carrier-label, and OpenSCAD export paths so the same result can be replayed
> deterministically.

This tutorial shows how GENtle turns an arrangement-ready Gibson result into a
physical rack view and one README-worthy exported figure.

The biology setup is the same deterministic single-insert Gibson example used
by the arrangement tutorial. What changes here is the downstream projection:

1. the arrangement already tells GENtle the semantic sample order,
2. the linked rack draft turns that order into physical A1-style positions,
3. the physical template turns that placement into:
   - a pseudo-3D/isometric SVG,
   - a fabrication/planning SVG,
   - carrier labels,
   - and a parameterized OpenSCAD file.

The figure result this tutorial aims for is:

- [`docs/figures/gibson_single_insert_rack_isometric.svg`](../figures/gibson_single_insert_rack_isometric.svg)

That SVG is suitable as a bench-facing communication artifact and as a README
hero figure because it is exported from the same shared project state rather
than drawn by hand.

## Fastest Setup

Recommended route:

1. `File -> Open Tutorial Project...`
2. choose `Gibson Arrangements Starter Project (offline)`
3. the Help window should open automatically on the arrangement guide
4. switch to this physical-rack guide and continue below from `Step 1`

Why reuse that starter:

- it already contains the deterministic Gibson result,
- it already contains the stored arrangement,
- and that arrangement already has a linked default rack draft.

So this tutorial can focus entirely on the physical-rack layer.

## What You Will Test

By the end of this tutorial, you should have verified all of these:

- the arrangement starter already contains one linked rack draft
- the rack draft preserves the arrangement order instead of inventing a new
  sample order
- GENtle can export one pseudo-3D/isometric SVG for that rack
- the same rack can also export:
  - fabrication/planning SVG
  - carrier labels SVG
  - OpenSCAD
- the isometric SVG can be regenerated deterministically from the CLI and used
  as a hero figure in the top-level README

## Inputs

This tutorial starts from the same deterministic local setup as the arrangement
tutorial:

- destination vector:
  [`test_files/pGEX-3X.gb`](../../test_files/pGEX-3X.gb)
- insert:
  [`docs/tutorial/inputs/gibson_insert_demo.fa`](./inputs/gibson_insert_demo.fa)
- arrangement-ready starter workflow:
  [`docs/examples/workflows/gibson_arrangements_baseline.json`](../examples/workflows/gibson_arrangements_baseline.json)

## Step 1: Open the Arrangement-Ready Gibson Starter

GUI:

1. `File -> Open Tutorial Project...`
2. choose `Gibson Arrangements Starter Project (offline)`

What to verify:

- `gibson_destination_pgex` is present
- `gibson_insert_demo` is present
- `gibson_destination_pgex_with_gibson_insert_demo` is already present
- the project already contains one arrangement and one linked rack draft

## Step 2: Inspect the Arrangement Row

GUI:

1. return to the main lineage page
2. scroll to `Arrangements`

What to verify:

- there is one Gibson-derived serial arrangement
- the row shows a linked physical rack status
- the row still tells the semantic story:
  - vector
  - insert
  - assembled product

Interpretation:

- the arrangement remains the semantic experiment order
- the rack is the physical carrier projection layered on top of it

## Step 3: Open the Rack

GUI:

1. in the arrangement row, click `Open Rack`

What to verify:

- a dedicated `Rack` window opens
- the rack uses A1-style coordinates
- occupied positions correspond to the Gibson story
- the rack is still driven by arrangement order, not by a separate manually
  typed sample list

## Step 4: Choose the Physical Template

In the `Rack` window:

1. find `Physical template`
2. choose `Pipetting PCR tube rack`

Why this template is the best README-facing export:

- it reads as a sturdier wet-lab carrier
- it makes the pseudo-3D rack geometry visually clearer
- it is the most natural bridge from GENtle’s digital arrangement model toward
  real bench handling and later 3D printing

You can also compare it against `Storage PCR tube rack` if you want the thinner
storage-oriented carrier variant.

## Step 5: Export the Hero Figure

In the same `Rack` window:

1. click `Isometric SVG...`
2. save it as:
   - `docs/figures/gibson_single_insert_rack_isometric.svg`

What to verify in the exported SVG:

- it is clearly not a screenshot; it is a deterministic vector export
- the rack body is rendered as a pseudo-3D/isometric carrier
- occupied positions are visible as colored sample caps
- ladders are visually distinct from the sequence-backed sample positions
- the export still carries the physical template identity

This is the result intended for README reuse.

## Step 6: Export the Other Physical Projections

Still in the `Rack` window, optionally export:

1. `Fabrication SVG...`
2. `Carrier labels SVG...`
3. `OpenSCAD...`

What each one is for:

- `Fabrication SVG`
  - planning / bench communication
  - top-view carrier geometry
- `Carrier labels SVG`
  - front strip + module labels
- `OpenSCAD`
  - parameterized printable geometry for later refinement or printing

Practical macOS note:

- before installing anything just to inspect or render the exported `.scad`
  file, first check whether this app binary already exists:
  `/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD`
- if it does, prefer using that direct app binary
- do not assume the Homebrew `openscad` cask is the right path

This is the important conceptual point:

- one arrangement
- one linked rack draft
- several downstream physical exports

No second sample-description system is needed.

## CLI Parity

Build the same arrangement-ready state in a fresh temporary project:

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_rack_tutorial.state.json \
  workflow @docs/examples/workflows/gibson_arrangements_baseline.json
```

Inspect the rack first:

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_rack_tutorial.state.json \
  racks show rack-1
```

Then export the README-facing figure:

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_rack_tutorial.state.json \
  racks isometric-svg \
  rack-1 \
  docs/figures/gibson_single_insert_rack_isometric.svg \
  --template pipetting_pcr_tube_rack
```

Optional physical follow-up exports:

```bash
cargo run --bin gentle_cli -- \
  --state /tmp/gibson_rack_tutorial.state.json \
  racks fabrication-svg \
  rack-1 \
  /tmp/gibson_single_insert_rack.fabrication.svg \
  --template pipetting_pcr_tube_rack

cargo run --bin gentle_cli -- \
  --state /tmp/gibson_rack_tutorial.state.json \
  racks carrier-labels-svg \
  rack-1 \
  /tmp/gibson_single_insert_rack.carrier.svg \
  --template pipetting_pcr_tube_rack \
  --preset front_strip_and_cards

cargo run --bin gentle_cli -- \
  --state /tmp/gibson_rack_tutorial.state.json \
  racks openscad \
  rack-1 \
  /tmp/gibson_single_insert_rack.scad \
  --template pipetting_pcr_tube_rack
```

## Why This Figure Matters

This is a small but important shift in what GENtle can communicate.

The same shared project state can now answer not only:

- what construct was planned,
- what primers were created,
- what product resulted,

but also:

- how the physical carrier is laid out,
- what labels belong on it,
- and what manufacturable geometry corresponds to that layout.

That is why this figure is a good README hero candidate: it ties GENtle’s
process abstraction back to the real world without leaving the deterministic
engine path.
