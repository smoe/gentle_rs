# Synthetic Regulatory/TSS Motif Package

All coordinates, scores, ownership and regulatory annotations are hand-crafted.
The exact JASPAR accessions label 16/11-base test spans, not scores obtained from
the actual PFMs. The package is not evidence about TP73, PATZ1 or any organism.
The two source-taxon strings are deliberately synthetic. The reference/catalog
ID, assembly name and accession deliberately differ.

Layout source: `jaspar-mapping` commit
`77acb3846631195f154793ad449d612eeaa6fc84`,
`scripts/manage_regulatory_tfbs.py` (`prepare`, `finalize`) and
`scripts/export_regulatory_tfbs.py`. GENtle's generator reproduces their final
manifest, inventories, materialized catalog and Parquet column contracts. It
does not invoke that repository, download data or require the original atlas.

Recreate from GENtle's root with Python 3 and the DuckDB CLI:

```sh
python3 test_files/fixtures/regulatory_motif_subset/make_fixture.py '/tmp/regulatory toy'
```

The destination must not exist. Regeneration computes hashes over exact newly
written bytes; DuckDB versions may produce different Parquet/catalog bytes.
The contents and expected queries, not cross-version binary hashes, are stable.
There is no `complete.json`. Original-atlas inventory paths are deliberately
unavailable and their counts deliberately differ from delivered counts.

Used by `src/genomic_motif_evidence/regulatory/tests.rs` and the
[offline walkthrough](../../../docs/regulatory_motif_subset.md). Chromosome 1
has seven orientation records, two TSS windows and shared TSS ownership.
Chromosome 2 has valid zero-row Parquets. MT has a known-empty regulatory
intersection, null payload paths and null digests, not an absent genome scan.
TP73 retains scores down to -5; PATZ1 to -1; both orientations, zero and
fractional negative scores are represented. Windows include the TSS base and
use 700 upstream/300 downstream bases (1,001 bp before chromosome clipping).

Tests create temporary packages and deliberately corrupt their copies. The
fixture is never loaded during application startup.
