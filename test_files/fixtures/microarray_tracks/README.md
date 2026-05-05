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

## `clariomd.synthetic.hg19_projected.manifest.json`

- Origin: hand-crafted companion manifest using the same synthetic probeset TSVs
  as native `hg19` / `GRCh37` coordinates.
- Deterministic recreation:
  - copy the direct synthetic manifest
  - set `coordinate_system` to `hg19`
  - add `clariomd.synthetic.hg19-to-hg38.tsv`, a one-block interval projection
    from `hg19 chr1:1001..1100` to `hg38 chr1:2001..2100`
  - anchor the target test sequence to `hg38 chr1:2001..2100`
- Primary usage:
  - engine tests for build-mismatched microarray projection through an explicit
    coordinate map
  - tooltip/detail tests for native interval versus displayed interval
- Runtime relevance: models the intended GRCh37-to-GRCh38 array-track path
  without committing large chain files or raw array data.
