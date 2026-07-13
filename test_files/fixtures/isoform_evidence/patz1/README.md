# Synthetic PATZ1 isoform-evidence fixture

## Origin

All files in this directory are hand-crafted synthetic test data. They are not
experimental PATZ1 measurements and must not be cited as biological evidence.
The coordinate frame is a 240 bp artificial slice labelled as
`GRCh38.p14:chr22:31325800..31326039`; transcript geometry is deliberately on
the minus strand so tests exercise the distinction between genomic-ascending
and transcript 5'-to-3' order.

The probe rows use Clariom-style identifiers only to exercise the existing
probe-evidence adapter. The separately maintained vendor-derived minimal table
is documented in
`test_files/fixtures/affymetrix_clariom_d_human_na36_hg38_subset/README.md`.
No Thermo Fisher row was copied into this synthetic bundle.

## Deterministic recreation

1. Create a 240 nt all-`A` GenBank region with accession-region metadata
   `31325800..31326039` and the three transcript joins recorded in
   `patz1_minus_strand.gb`.
2. Map local 1-based coordinates to the labelled genome slice with
   `genomic = 31325800 + local - 1`.
3. Mint geometry ids as
   `EXF:GRCh38.p14:start-end:-` and
   `JCT:GRCh38.p14:low-high:-`.
4. Keep expression, cDNA/EST, and probe rows exactly as recorded in the local
   TSV/JSON files. No timestamps or network results enter the fixture.

## GENtle use

- Engine tests load the GenBank slice, import `patz1_isoform_panel.json`, and
  inspect `FeatureExpertTarget::IsoformEvidence`.
- The cDNA/EST resource tests observed-vs-unknown junction handling.
- Its explicit GRCh37 mismatch row is a negative control: matching text in a
  geometry id must not override incompatible coordinate provenance.
- The probe report tests that ambiguous array geometry remains
  `constraint_only`, never direct validation.
- The expression TSV tests dataset-relative abundance without inferring a
  cross-platform absolute scale.
