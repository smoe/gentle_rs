# Probe-Region Adapter Input Fixtures

## `clariom_e_mtab_14704_tp73_glen_style.csv`

- Schema: `gentle.glen_probe_region_adapter_input.v1`
- Origin: compact, hand-curated rows from Glen's exploratory E-MTAB-14704 TP73
  probe-region analysis shape, reduced to the columns needed to regenerate
  GENtle's canonical helper-output contract.
- Raw CEL files, vendor annotation archives, R outputs, and Glen's full
  exploratory analysis snapshot are not committed.
- Coordinates are TP73 GRCh38.p14, anchored to `test_files/tp73.ncbi.gb`
  (`chr1:3652516..3736201`).
- The `P_SKMel29_*.CEL` values are log2-scale representative PM intensities.
  The adapter preserves these per-probe rows and aggregates parent
  probeset-region rows deterministically by arithmetic mean.
- Two compact junction/constraint rows are included to ensure the downstream
  interpretation report exercises exon-overlap, junction-spanning, and
  non-exonic transcript-span geometry for later figure-renderer tests.

Required columns:

- `probeset_id`, `probe_id`, `chromosome`, `parent_start`, `parent_stop`
- `start`, `stop`, `strand`
- `transcript_cluster_id`, `gene_symbol`, `intensity_source`
- `x`, `y`
- one or more sample columns ending in `.CEL`
- optional `mean_log2_*`, `sd_log2_*`, `log2FC_*`, and `adj_p_*` columns

Only sample, `mean_log2_*`, `sd_log2_*`, and `log2FC_*` columns are emitted in
the canonical GENtle CSV outputs. `adj_p_*` columns document the compact source
shape for future adapter extension but are intentionally not emitted yet.
