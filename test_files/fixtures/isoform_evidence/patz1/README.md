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
The declared `Clariom_D_Human` platform, `synthetic_case_vs_control` contrast,
and log2-fold-change values are hand-crafted metadata used to verify that
GENtle preserves selection evidence without claiming exact probe reuse.

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

Keep `patz1_probe_evidence.json` as UTF-8 with LF line endings, including on
Windows (`.gitattributes` enforces this). Its exact-byte SHA-256 is
`11e2a7b4cbb83af5868902b1834ecbfe6c5e7ddb80e521c1060abf7a4d8418a9`.
The PATZ1 endpoint/SYBR tutorial's junction and common-region reports retain
this digest; converting LF to CRLF changes the evidence identity even though
the parsed JSON is equivalent. Engine hashing deliberately does not normalize
user evidence. `patz1_probe_fixture_digest_matches_retained_tutorial_provenance`
checks the checkout against the retained report, and the focused PATZ1 workflow
test compares all three regenerated assay reports with the tutorial snapshots.

## GENtle use

- Engine tests load the GenBank slice, import `patz1_isoform_panel.json`, and
  inspect both `FeatureExpertTarget::IsoformEvidence` and the composed
  `FeatureExpertTarget::GeneLocusEvidence` report.
- The cDNA/EST resource tests observed-vs-unknown junction handling.
- Its explicit GRCh37 mismatch row is a negative control: matching text in a
  geometry id must not override incompatible coordinate provenance.
- The probe report tests that ambiguous array geometry remains
  `constraint_only`, never direct validation.
- The expression TSV tests dataset-relative abundance without inferring a
  cross-platform absolute scale.
- Synthetic projected occupancy intervals, motif scoring, and a persisted
  qPCR candidate are added by the engine test at runtime; they are not encoded
  in these fixture files or presented as PATZ1 measurements.
