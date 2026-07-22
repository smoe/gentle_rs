# TP73 CUT&RUN Offline Release-Smoke Assets

These files are a tiny deterministic integration fixture. They are not
experimental CUT&RUN measurements and must not be used to make a biological
claim about TP73 occupancy.

## Origin

- `tp73_release_smoke.peaks.bed` is hand-crafted synthetic BED6 evidence on
  the public GRCh38.p14 TP73 anchor `NC_000001.11:3652516..3736201`.
- The paired FASTA records are exact 40 bp substrings of the committed public
  reference record `test_files/tp73.ncbi.gb`. R1 records use local 1-based
  starts 200, 215, 800, and 815. R2 records are reverse complements of the
  40 bp reference substrings beginning at local starts 261, 276, 861, and 876.
- The intervals and pair placements were selected solely to produce two small
  deterministic support windows for software testing.

## Deterministic Recreation

From the repository root, extract the GenBank `ORIGIN` sequence, remove digits
and whitespace, uppercase it, then emit the R1 substrings and reverse
complements of the R2 substrings at the positions listed above. The BED file is
hand-authored and should remain byte-stable unless the smoke geometry changes.

## GENtle Usage

`assets/cutrun.json` exposes these files as
`tp73_release_smoke_synthetic`. The TP73 release-proof workflow and CUT&RUN
engine tests use the entry to cover dataset discovery/preparation, BED
projection, paired-read ROI interpretation, support-window reasoning, and JSON
or TSV export without a genome or SRA download.
