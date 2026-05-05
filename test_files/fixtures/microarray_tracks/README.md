# Microarray Track Fixtures

## `clariomd.synthetic.manifest.json`

- Origin: hand-crafted synthetic Clariom D-like manifest for GENtle microarray
  track projection tests.
- Deterministic recreation:
  - create a 100 bp synthetic sequence anchored to `hg38 chr1:1001..1100`
  - define two probeset-level contrasts in manifest order
    `AdTAp73alpha-AdGFP`, then `AdTAp73beta-AdGFP`
  - add rows on `chr1` that overlap the anchor, one row outside the anchor,
    and one row on `chr2` to exercise strict streaming filters
  - use small fixed `logFC` and `adj.P.Val` values for heat-color and tooltip
    assertions
- Primary usage:
  - engine manifest parsing and contrast ordering tests
  - genome-anchor coordinate compatibility rejection tests
  - forward/reverse microarray projection tests
  - GUI helper tests for array tooltip/detail text and feature-tree filters
- Runtime relevance: mirrors the compact TSV contract produced by the Rostock
  Clariom D analysis script without committing raw CEL data.
