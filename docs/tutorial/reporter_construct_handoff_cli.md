# Reporter Construct Handoff Tutorial

This tutorial shows how to use GENtle's read-only reporter construct handoff
planner after a promoter-reporter candidate set has already been saved.

The planner connects three existing pieces:

- a `gentle.promoter_reporter_candidates.v1` JSON report
- the offline reporter recommender
- the built-in macro template
  `allele_paired_promoter_luciferase_reporter`

It does not create constructs, fetch live registries, optimize codons, or make
wet-lab claims. Its job is to produce one auditable JSON plan that says which
inputs are ready, which are derivable, which are missing, and which commands a
user or agent may run after review.

## When To Use This

Use this tutorial when you already have a saved promoter-fragment candidate
report and want to answer:

- which promoter candidate is selected
- which luciferase reporter GENtle recommends from the local catalog
- whether the reference insert, alternate insert, and reporter backbone are
  ready in the current project state
- which macro command would validate or build the paired reporter constructs

For the VKORC1 tutorial data, the candidate report is already committed here:

- [`docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json`](./reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json)

## Step 1: Generate The Handoff Plan

Run the direct CLI route from the repository root:

```bash
mkdir -p /tmp/gentle-reporter-handoff

cargo run --quiet --bin gentle_cli -- \
  reporters plan-handoff docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json \
  --backbone-seq-id gentle_mammalian_luciferase_backbone_v1 \
  --backbone-path data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb \
  --reference-fragment-seq-id vkorc1_rs9923231_promoter_reference \
  --alternate-fragment-seq-id vkorc1_rs9923231_promoter_alternate \
  --output-prefix vkorc1_rs9923231_reporter_pair \
  --output /tmp/gentle-reporter-handoff/reporter_construct_handoff.json
```

The command also prints the normal `OpResult` JSON to stdout. The file named by
`--output` contains only the portable handoff plan:

```bash
python3 -m json.tool /tmp/gentle-reporter-handoff/reporter_construct_handoff.json
```

## Step 2: Read The Plan

Start with these fields:

- `schema`
  - should be `gentle.reporter_construct_handoff.v1`
- `status`
  - `ready`: all required macro ports are already available in the inspected
    state
  - `needs_inputs`: the planner found a compatible route, but at least one
    input still needs to be derived or loaded
  - `no_compatible_reporter`: no local catalog reporter matched the constraints
  - `no_compatible_macro_route`: the requested reporter class is outside V1's
    luciferase macro route
- `provenance`
  - records the candidate-set path, source report ids when present, reporter
    catalog path, and macro template id
- `selected_fragment`
  - names the chosen promoter fragment and the reference/alternate fragment ids
    expected by the macro
- `selected_reporter`
  - records the offline luciferase recommendation and score
- `backbone`
  - says whether the reporter backbone is already in state, unresolved, or must
    be loaded manually
- `port_bindings`
  - lists each macro port with one typed status:
    `ready`, `provided_missing_from_state`, `derivable`, or `missing`
- `commands`
  - gives the reviewed follow-up commands in execution order

For a fresh direct CLI run, `needs_inputs` is expected: the command reads the
candidate-set JSON and local reporter catalog, but it does not silently load
the original sequence, materialize allele fragments, or load the reporter
backbone into a persistent project.

## Step 3: Use A Persistent State When You Want Readiness

If you want the plan to distinguish `ready` from `derivable`, run it through
the shared operation route against the same state file you used for promoter
candidate generation:

```bash
STATE=path/to/vkorc1_promoter_project.gentle.json

cargo run --quiet --bin gentle_cli -- --state "$STATE" \
  op '{"PlanReporterConstructHandoff":{"candidate_set_path":"docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json","reporter_backbone_seq_id":"gentle_mammalian_luciferase_backbone_v1","reporter_backbone_load_path":"data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb","reference_fragment_seq_id":"vkorc1_rs9923231_promoter_reference","alternate_fragment_seq_id":"vkorc1_rs9923231_promoter_alternate","output_prefix":"vkorc1_rs9923231_reporter_pair","path":"/tmp/gentle-reporter-handoff/reporter_construct_handoff_stateful.json"}}'
```

This route is still read-only. It inspects state and writes the plan, but it
does not create or modify sequences.

## Step 4: Review Before Running Follow-Up Commands

The plan's `commands` array is intentionally explicit. Review the commands
before running them, especially any row with `"mutating": true`.

For the VKORC1 flow, the command list should include:

- extract the selected promoter fragment
- materialize the reference-allele promoter insert
- materialize the alternate-allele promoter insert
- load the local luciferase backbone if needed
- import the built-in macro templates
- validate the macro binding with
  `allele_paired_promoter_luciferase_reporter`
- run the macro only after validation

The validation command is the useful pause point:

```bash
cargo run --quiet --bin gentle_cli -- --state "$STATE" \
  shell "macros template-run allele_paired_promoter_luciferase_reporter --bind reference_fragment_seq_id=vkorc1_rs9923231_promoter_reference --bind alternate_fragment_seq_id=vkorc1_rs9923231_promoter_alternate --bind reporter_backbone_seq_id=gentle_mammalian_luciferase_backbone_v1 --bind overlap_bp=20 --bind output_prefix=vkorc1_rs9923231_reporter_pair --transactional --validate-only"
```

Only after that succeeds should you decide whether to run the non-validate
macro command from the same plan.

## Step 5: Change The Candidate Or Catalog Deliberately

To inspect a non-default promoter candidate, pass `--candidate-id`:

```bash
cargo run --quiet --bin gentle_cli -- \
  reporters plan-handoff docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json \
  --candidate-id vkorc1_rs9923231_context_enst00000420057_2091_3501_promoter_fragment \
  --output /tmp/gentle-reporter-handoff/reporter_construct_handoff_candidate_2.json
```

To use a curated local reporter catalog rather than the built-in catalog, pass
`--catalog PATH`. Normal recommendations remain offline; live source checks
belong to explicit catalog refresh or provenance-review work, not this planner.

## Expected Outcome

At the end of this tutorial you should have:

- a `gentle.reporter_construct_handoff.v1` JSON plan
- a selected promoter fragment and selected local luciferase reporter
- typed readiness for the macro's fragment and backbone inputs
- an explicit validate-only macro command
- a reviewed next step for construct creation, without GENtle silently creating
  constructs during planning

Use the longer VKORC1 walkthrough for the surrounding GUI story:

- [`docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`](./vkorc1_warfarin_promoter_luciferase_gui.md)

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
