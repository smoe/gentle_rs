# TP73 Evidence Viewer Fixtures

These fixtures support the TP73 genome-anchored evidence-viewer release proof.
They are intentionally tiny and deterministic so CI can exercise repeat,
microarray, CUT&RUN-style, TFBS, and SVG display paths without downloading a
full genome, full UCSC `rmsk`, raw CEL files, or SRA reads.

## Source Anchor

- `test_files/tp73.ncbi.gb`
  - public NCBI GenBank TP73 locus already committed in the repository
  - sequence anchor: chromosome 1, GRCh38.p14 Primary Assembly,
    3652516..3736201
  - used as the real genome-anchored TP73 locus for the proof workflow

## Fixtures

- `tp73_evidence_viewer.rmsk.hg38.txt`
  - hand-crafted UCSC `rmsk`-style rows on the TP73 anchor interval
  - coordinates are chosen inside or just beyond the committed TP73 locus to
    exercise overlap projection and clipping
  - values are synthetic and must not be interpreted as biological repeat calls
- `tp73_evidence_viewer.rmsk.hg38.json`
  - generated resource snapshot from the text fixture
- `tp73_evidence_viewer.rmsk.hg38.interval-index.json`
  - generated interval index used by `features materialize-repeats`
- `clariomd.tp73_evidence_viewer.manifest.json`
  - synthetic Clariom D-style track manifest using public dataset label
    `E-MTAB-14704`
  - coordinates are GRCh38.p14-compatible with the committed TP73 locus
- `clariomd.tp73_evidence_viewer.AdTAp73alpha-AdGFP.tsv`
  and `clariomd.tp73_evidence_viewer.AdTAp73beta-AdGFP.tsv`
  - synthetic per-contrast probeset rows for deterministic projection tests
  - values are proof data only, not derived from CEL files
- `tp73_cutrun_demo.bed`
  - synthetic CUT&RUN-style BED6 intervals on the TP73 locus
  - names mirror TP73-GFP proof-track use but are not derived from raw reads

## Regeneration

Regenerate the UCSC repeat resource sidecars from the repository root:

```bash
cargo run --quiet --bin gentle_cli -- resources sync-ucsc-rmsk \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.txt \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.json \
  --assembly hg38

cargo run --quiet --bin gentle_cli -- resources prepare-ucsc-rmsk-index \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.json \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.interval-index.json
```

The microarray TSVs and BED file are hand-authored fixtures. Update this README
whenever a fixture is replaced by a public-data-derived subset. The committed
`tp73_evidence_viewer.rmsk.hg38.json` snapshot pins `fetched_at_unix_ms` to
`123` so fixture diffs stay deterministic; the CLI regeneration command records
the current wall-clock time, so normalize that field before committing a
regenerated snapshot.
