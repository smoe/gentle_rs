# Synthetic UTR Core Contract Cases

`cases.json` is hand-crafted, synthetic test data authored for the first shared
UTR core slice on 2026-11-08. No provider annotation, actual gene, read, barcode,
private sequencing data or biological measurement is included. Names, hashes,
coordinates, exon/CDS/stop intervals and query blocks are literal fixtures.

The corpus binds V2 models and gives compact observation/path inputs plus
independently stated expected metric triples `[value_or_null, status, reason]`.
Numeric censored values are lower bounds, never ordinary comparator values.
The fixture schema is test-only, not a product wire or observed-SAM importer.

Recreate the inputs from the literal JSON; no external download or generator
is needed. Two consumers independently expand its compact inputs:

- Rust `observed_transcript_structures::read_ends::tests` constructs existing
  typed `ObservedSamAnalysis` assessments and exercises the new pure core.
- `tests/test_rna_read_end_core_contract.py` uses the unchanged standalone Python
  projector and checks the same literal outcomes under LF and CRLF encodings.

Run the Python half with
`python3 -m unittest discover -s tests -p test_rna_read_end_core_contract.py -v`.
Glen/CI runs the Rust half with
`cargo test --locked --lib -j1 read_ends::tests -- --test-threads=1`.
Neither passing Python alone nor source inspection proves that Rust compiles or
matches it. No native, file-publication or biological acceptance is implied.
