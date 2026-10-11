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

## Compatible 08.16 companion

The optional `--example tss-view` preset reuses the same writer/reader contract.
The default invocation retains the original reference, seven rows, matrices,
empty-contig cases and schema; it does not become a TSS-view fixture implicitly.
A named, bounded preset avoids unsupported arbitrary combinations of reference
IDs, matrix widths and coordinates.

```sh
python3 test_files/fixtures/regulatory_motif_subset/make_fixture.py /tmp/tss-import-lesson/package --example tss-view --duckdb duckdb
```

This preset deliberately matches `test_files/fixtures/tss_profiles/manifest.json`:
`synthetic-genome-v1`, `synthetic-assembly-v1`, `synthetic_chr`, plus BED window
`[96,107)` and minus `[292,303)`. It supplies five full-width six-base MA0004.1
footprints: plus-window hits at `[97,103)` on both genomic strands, minus-window
hits `[294,300)` genomic minus and `[296,302)` genomic plus, and genomic-plus
`[290,296)` clipped at the minus window's lower boundary. Scores 1.25, 0 and
-0.5 are **hand-crafted log2_relative_risk**, not calculated Arnt PFM scores,
affinities or occupancy. The synthetic source-audit placeholders describe no
real Ensembl audit or executed atlas producer.

Partial coverage is produced by a real query limited to `[294,301)`, not by
altering the package's completeness or hand-editing a query report. Keep full
and partial reports separate. The real export path hashes and embeds their
source bytes and validates geometry. Tutorial 08.16 gives exact query/export
commands into fresh output directories; existing teaching reports remain intact.

Consumers: `scripts.test_tutorial_gap_lessons` checks source geometry and refuses
invalid/clobbering preset requests without invoking DuckDB.
`scripts.test_regulatory_motif_fixture` optionally exercises actual DuckDB and
Parquet generation, including the unchanged default's logical tables; it is not
GENtle reader/export acceptance. `regulatory_tss_view_companion_query_export_preserves_geometry`
is an explicitly ignored real-DuckDB/shared-export Rust integration test. Glen
must run it and supply actual export receipts/native screenshots. No package
or execution receipt is manufactured by the Python source-plan tests.
