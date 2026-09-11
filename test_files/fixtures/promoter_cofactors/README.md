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
