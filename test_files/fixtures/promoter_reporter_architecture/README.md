# Ensembl 116 SERPINE1 promoter-reporter fixture

## Origin

This acceptance fixture is a pinned, read-only derivative of the official
Ensembl REST records for human `ENSG00000106366` on GRCh38. It was retrieved
on 2026-08-30, when `/info/data` reported release 116.

The fixture contains:

- `ensembl_116_serpine1_geometry.json`: the gene, all 15 protein-coding
  transcript models, exon spans, and translation boundaries needed by the
  transcript-aware architecture test.
- `ensembl_116_serpine1_region.fa`: the exact forward-strand sequence for
  chromosome 7:101121158-101139263 (1-based inclusive).

Upstream response hashes:

- `/info/data`: `86e87c46ae1a424702676493d024d0508ce60f6c951b0b2803ce10269e09c20c`
- expanded lookup response: `f1013a8cb853de4011114a13adb812564742db320259fe1ddd35137e94568f7c`
- region sequence response: `440d6b10a5affeda66358d2d4c102ebc0b0a8a68a79eedfcb94748a02b6297e0`

## Deterministic retrieval

```sh
curl -s 'https://rest.ensembl.org/info/data?content-type=application/json' \
  -o ensembl116_info_data.json
curl -s 'https://rest.ensembl.org/lookup/id/ENSG00000106366?expand=1;utr=1;content-type=application/json' \
  -o ensembl116_serpine1_lookup.json
curl -s 'https://rest.ensembl.org/sequence/region/human/7:101121158..101139263:1?content-type=text/plain' \
  -o ensembl116_serpine1_region.txt
shasum -a 256 ensembl116_info_data.json ensembl116_serpine1_lookup.json \
  ensembl116_serpine1_region.txt
```

The committed geometry JSON is a deterministic projection of the expanded
lookup, not a second biological annotation. Its `provenance` object retains
the upstream URLs, release, assembly, and raw-response hashes.

## GENtle use

`src/engine/analysis/promoter_reporter_architecture.rs` uses the fixture for
an offline acceptance test. The test rebuilds an annotated locus, binds the
GRCh38 interval and source hashes as extraction provenance, and verifies the
observed TSS classes, common CDS start, and representative 5-prime-UTR/intron
geometry. It does not claim cell-context promoter use or contain CUT&RUN
evidence.

To emit a deterministic, inspectable report/SVG pair while running only that
acceptance test:

```sh
GENTLE_SERPINE1_REPORTER_ARTIFACT_DIR=/tmp/gentle-serpine1-reporter \
  cargo test ensembl_116_serpine1_fixture_matches_observed_representative_geometry \
  --lib -- --nocapture
```

The ordinary test run does not write artifacts. The portable request example
for a project containing this annotated locus is
`docs/examples/assets/serpine1_promoter_reporter_architecture_request.json`.
