# Synthetic PATZ1 transcript-assay fixture

## Origin

`patz1_assay_minus_strand.gb` is hand-crafted synthetic test data. It is not a
reference sequence and must not be used to order laboratory primers. Its
GRCh38.p14-like local anchor (`chr22:31325800..31326039`) and three minus-strand
transcript models intentionally mirror the geometry in
`test_files/fixtures/isoform_evidence/patz1/`, while the sequence is diverse
enough to exercise primer selection.

## Deterministic recreation

The file is recreated by concatenating the six 40 nt blocks printed in its
`ORIGIN` section. The transcript locations are fixed GenBank joins:

- `PATZ1-201`: `complement(join(21..60,101..140,181..220))`
- `PATZ1-202`: `complement(join(21..60,181..220))`
- `PATZ1-203`: `complement(join(21..60,121..140,181..220))`

Thus the annotated mature cDNAs are 120, 80, and 100 nt. Their common outer
exons permit endpoint primer pairs to expose three explicit band sizes, while
the skipped-middle `PATZ1-202` junction matches the synthetic
Clariom JUC geometry in `patz1_probe_evidence.json`.

## GENtle use

The fixture is used by transcript-assay engine tests and the offline PATZ1
endpoint/SYBR workflow example. It validates minus-strand transcript ordering,
the annotated first x terminal reaction matrix, JUC-required short SYBR design,
and oligo-dT interpretation warnings.
