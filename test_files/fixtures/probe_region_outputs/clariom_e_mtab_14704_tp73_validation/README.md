# Clariom E-MTAB-14704 TP73 Validation Fixture

## Origin

- Source accession: `E-MTAB-14704`
- Source resource: `data/publication_resources/rostock_p73_clariomd_e_mtab_14704/`
- Raw CEL files are fetched explicitly by the committed `download.sh` and are not committed.
- This fixture is derived/synthetic-on-real-coordinates: TP73-compatible Clariom-style
  coordinates and public sample names are used with fixed synthetic intensity values.
- No CEL files, vendor annotation ZIPs, R, Bioconductor, or APT binaries were used to create
  this committed fixture.

## Anchor And Coordinate Assumptions

- TP73 anchor: `test_files/tp73.ncbi.gb`
- Genome build: GRCh38.p14
- Chromosome: `chr1`
- Anchor interval: 3652516..3736201
- Fixture coordinate system: `hg38`
- Rows are intentionally placed near the 5' TP73 anchor edge so deterministic tests can use a
  short synthetic sequence carrying the same genome-anchor metadata.

## Expected Contents

- `provenance.json`
  - schema: `gentle.probe_region_backend_provenance.v1`
  - `probe_intensity_source: probe_level_input`
- `normalized_feature_matrix_manifest.json`
  - schema: `gentle.probe_region_normalized_matrix_manifest.v1`
  - targets: `probeset`, `pm_probe`
- `region_intensity_chrom_order.csv`
  - 4 data rows
- `probe_intensity_chrom_order.csv`
  - 8 data rows
  - all rows are true PM probe input rows (`intensity_source=probe_level_input`)

Representative rows:

- `PSR0100145779.hg.1` / `PM_TP73_0001`: shared first-exon overlap,
  chr1:3652520..3652537
- `PSR0100145780.hg.1` / `PM_TP73_0003`: transcript-geometry-constraining overlap for one
  synthetic TP73 transcript model, chr1:3652572..3652586
- `PSR0100145781.hg.1` / `PM_TP73_0005`: complementary transcript-geometry-constraining overlap,
  chr1:3652594..3652610
- `JUC010000001.hg.1` / `PM_TP73_0007`: junction-spanning geometry review row,
  chr1:3652550..3652574

## Runtime Use

This fixture is a release-grade validation aid for presentation and routing:

- projection with `arrays project-probe-region-output --level pm_probe`
- interpretation with `arrays interpret-probe-region-evidence --level pm_probe`
- GUI action-to-capability parity tests

It is not biological evidence and must not be used to claim probe specificity, multi-hit status,
or isoform support.
