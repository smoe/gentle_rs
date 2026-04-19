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
- `sequencing_confirmation_trace_demo_baseline.fa`

Purpose:

- Minimal expected-construct plus baseline/reference pair for the sequencing-
  confirmation tutorials.
- Lets one bundled ABI/AB1 trace confirm a known local junction target and one
  intended tutorial SNP without adding another binary tutorial artifact.

Provenance:

- Both files are hand-crafted local tutorial inputs derived from the first
  `48` called bases of the committed ABI fixture
  `test_files/fixtures/sequencing_confirmation/3100.ab1`.
- Exact expected sequence (`trace_demo_construct`):
  `CAAGATTGCATTCATGATCTACGATTACTAGCGATTCCAGCTTCATAT`
- Exact baseline/reference sequence (`trace_demo_baseline`):
  `CAAGATTGCATTCATGATCTACGGTTACTAGCGATTCCAGCTTCATAT`
- The two differ by one tutorial SNP at position `24` (1-based):
  expected=`A`, baseline=`G`.

Current intended use in GENtle:

- Load `trace_demo_construct` in
  [`docs/tutorial/sequencing_confirmation_trace_cli.md`](../sequencing_confirmation_trace_cli.md)
  for the shell/CLI path.
- Load both `trace_demo_construct` and `trace_demo_baseline` in
  [`docs/tutorial/sequencing_confirmation_gui.md`](../sequencing_confirmation_gui.md)
  for the GUI-first path.
- Pair them with the bundled `3100.ab1` trace via either:
  - `seq-trace import` + `seq-confirm run --trace-id ...`, or
  - GUI `Patterns -> Sequencing Confirmation...` raw trace import.

## Stateless Sequence-Inspection Demo

File:

- `inline_sequence_inspection_demo.fa`

Purpose:

- Minimal synthetic local sequence for the direct restriction-site / TFBS hit /
  TFBS score-track tutorial.
- Lets the GUI path use one tiny FASTA while the CLI/ClawBio parity route uses
  the same bases through the shared `inline_sequence` operand.

Provenance:

- Hand-crafted synthetic tutorial input created directly in this repository.
- Exact sequence:
  `GAATTCCCGGGATCCGGGCGGGGCGCATGTGTAACAGGGGCGGGGC`
- The sequence intentionally contains:
  - `EcoRI` (`GAATTC`)
  - blunt `SmaI` (`CCCGGG`)
  - `BamHI` (`GGATCC`)
  - one GC-rich `SP1` teaching region plus one `TP73` teaching region

Current intended use in GENtle:

- Load in
  [`docs/tutorial/stateless_sequence_inspection_gui_cli.md`](../stateless_sequence_inspection_gui_cli.md)
  for the GUI toolbar path (`RE scan`, `TFBS scan`, `TFBS score tracks`).
- Compare against the matching workflow example
  [`docs/examples/workflows/inline_sequence_inspection_stateless_offline.json`](../../examples/workflows/inline_sequence_inspection_stateless_offline.json),
  which uses the same bases through state-optional inline-sequence operands.
