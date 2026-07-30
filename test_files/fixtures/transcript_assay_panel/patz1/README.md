# Synthetic PATZ1 transcript-assay fixture

## Origin

`patz1_assay_minus_strand.gb` is hand-crafted synthetic test data. It is not a
reference sequence and must not be used to order laboratory primers. Its
GRCh38.p14-like local anchor (`chr22:31325800..31326039`) and three minus-strand
transcript models intentionally mirror the geometry in
`test_files/fixtures/isoform_evidence/patz1/`, while the sequence is diverse
enough to exercise primer selection.

`patz1_assay_probe_evidence.json` is a mechanically adapted copy of the
synthetic `patz1_probe_evidence.json` fixture in that directory. Only `seq_id`
is changed so the same artificial PSR/JUC geometry can be composed with this
primer-design sequence. `patz1_assay_isoform_evidence.json` is the deterministic
`FeatureExpertView::IsoformEvidence` output generated from those two fixtures,
the shared synthetic isoform panel, and the shared synthetic expression table.
All values remain artificial and array overlap remains a design constraint, not
isoform validation.

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

Recreate the evidence inputs from the repository root:

```bash
jq '.seq_id = "patz1_transcript_assay_demo"' \
  test_files/fixtures/isoform_evidence/patz1/patz1_probe_evidence.json \
  > test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_probe_evidence.json

cargo run -q --bin gentle_cli -- --state /tmp/gentle-patz1-routine-fixture.json \
  op '{"LoadFile":{"path":"test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_minus_strand.gb","as_id":"patz1_transcript_assay_demo"}}'
cargo run -q --bin gentle_cli -- --state /tmp/gentle-patz1-routine-fixture.json \
  shell 'panels import-isoform patz1_transcript_assay_demo test_files/fixtures/isoform_evidence/patz1/patz1_isoform_panel.json --panel-id patz1_synthetic_v1'
cargo run -q --bin gentle_cli -- --state /tmp/gentle-patz1-routine-fixture.json \
  shell 'inspect-feature-expert patz1_transcript_assay_demo isoform-evidence patz1_synthetic_v1 --annotation-release synthetic_Ensembl116-like --probe-evidence test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_probe_evidence.json --expression-tsv test_files/fixtures/isoform_evidence/patz1/patz1_expression.tsv' \
  > test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_isoform_evidence.json
```

The generated isoform-evidence fixture has SHA-256
`267a8fc00da8fea653dd440259987b20fc68422554cfed8bbff61809d4eeaaf5`.
The workflow checks this digest before composing the routine.

## GENtle use

The fixture is used by transcript-assay engine tests and the offline PATZ1
endpoint/SYBR workflow example. It validates minus-strand transcript ordering,
the annotated first x terminal reaction matrix, JUC-required short SYBR design,
oligo-dT interpretation warnings, and pure-read composition of common-control,
junction-validation, and endpoint panels into
`gentle.gene_transcript_assay_routine.v1`.
