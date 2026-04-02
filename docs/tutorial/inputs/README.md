# Tutorial Input Fixtures

This folder contains small committed local inputs used by GUI tutorials and
manual contributor checks.

## Gibson Synthetic Sanity-Check Pair

Files:

- `gibson_synthetic_c100_a_g100_vector.gb`
- `gibson_synthetic_poly_a_insert.fa`

Purpose:

- Minimal destination/insert pair for inspecting Gibson overlap and primer
  logic without the richer pGEX annotation context.
- The circular destination is exactly `(C)^100 A (G)^100`.
- The insert is a short linear poly-A stretch (`23 bp`).

Recommended teaching opening:

- Use a defined opening with left/right cut edges `100` and `101` (0-based).
- That removes the single central `A` and exposes a pure `C` left flank and a
  pure `G` right flank.

Provenance:

- Both files are hand-crafted synthetic fixtures created directly in this
  repository for deterministic Gibson specialist testing and discussion.
- They are not derived from any biological plasmid or public accession.

Current intended use in GENtle:

- Load them via `File -> Open Sequence...` or `LoadFile` operations.
- Use them to inspect the `Patterns -> Gibson...` specialist window, especially
  overlap derivation, primer suggestions, and cartoon behavior on an easily
  understood sequence context.

## Sequencing Confirmation Trace Demo Construct

Files:

- `sequencing_confirmation_trace_demo_construct.fa`

Purpose:

- Minimal expected construct for the CLI trace-aware sequencing-confirmation
  tutorial.
- Lets one bundled ABI/AB1 trace confirm a known local junction target without
  adding another binary tutorial artifact.

Provenance:

- Hand-crafted local tutorial input derived from the first `48` called bases of
  the committed ABI fixture
  `test_files/fixtures/sequencing_confirmation/3100.ab1`.
- Exact sequence:
  `CAAGATTGCATTCATGATCTACGATTACTAGCGATTCCAGCTTCATAT`

Current intended use in GENtle:

- Load it as `trace_demo_construct` in
  `docs/tutorial/sequencing_confirmation_trace_cli.md`.
- Pair it with the bundled `3100.ab1` trace via `seq-trace import` and
  `seq-confirm run --trace-id ...`.
