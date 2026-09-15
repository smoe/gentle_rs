# Synthetic collaborator package

`fixture.sql` is hand-crafted test data, not a JASPAR scan, experimental
measurement, or proprietary study. Coordinates, factors, samples and statistics
are artificial. It models two anchors, shared promoter ownership, one detailed
matrix and one overview-only matrix. Both use fabricated accessions.

Recreation: `src/promoter_cofactors_tests.rs` executes this repository-authored
SQL in a fresh in-memory DuckDB, exports each table to a temporary Parquet,
and writes fresh SHA-256 manifest/completion inventories. No binary data or
production package is committed. The tests never execute package-supplied SQL.
Set `GENTLE_TEST_DUCKDB=/absolute/path/to/duckdb` to require the real-Parquet
reader tests; absent that variable they skip explicitly. Pure request,
manifest, assembly, relocation, sparse/coordinate and integrity tests do not
require DuckDB. The fixture also supports the GUI browser teaching workflow in
`docs/promoter_cofactor_browser.md`.

`manifest.template.json` is hand-crafted scope/score metadata shared by the Rust
tests and `scripts/promoter_cofactor_tutorial.py`. Its inventories start empty;
both producers calculate hashes over freshly written fixture Parquet files.
The tutorial script writes only a new directory and exercises the real shared
CLI reader and `regions capture`/export. Recreate with the command in the
registered [08-14 walkthrough](../../../docs/tutorial/08-14_promoter_cofactor_browser.md).
`scripts/test_promoter_cofactor_tutorial.py` verifies metadata and, with explicit
CLI/DuckDB paths, the complete offline replay. DuckDB versions may encode
different Parquet bytes; logical fixture values, not fixed binary hashes, are
the reproducibility contract.

`handoff_report()` in `src/promoter_cofactors_tests.rs` separately constructs
small synthetic transport records with two gene links, a hit extending outside
a promoter and distinct TP73/control summaries. No experiment is represented.
It exercises saved-region provenance, forward/reverse sequence projection,
assembly rejection and GUI stale-result/copy binding without a DuckDB runtime.
