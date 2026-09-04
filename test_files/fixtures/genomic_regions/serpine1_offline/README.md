# Portable genomic-region SERPINE1 fixture

## Origin

This public offline fixture combines the pinned Ensembl release 116 SERPINE1
annotation and sequence already documented in
`test_files/fixtures/promoter_reporter_architecture/README.md` with one
hand-crafted synthetic CUT&RUN-style interval. It contains no private Rostock
data and no measured regulatory result.

- `ensembl_116_serpine1_entry.json` is a deterministic adapter projection of
  `ensembl_116_serpine1_geometry.json` and
  `ensembl_116_serpine1_region.fa` into
  `gentle.ensembl_gene_entry.v1`. The source URLs, source hashes, assembly,
  release, gene/transcript/exon geometry, and exact sequence remain those of
  the pinned public fixture. Translation stable IDs were not present in the
  source geometry snapshot and remain empty; translation boundaries are
  retained.
- `synthetic_serpine1_cutrun.bed` is hand-crafted BED6 data on GRCh38
  chromosome 7. It is an artificial occupancy interval for data-flow and
  display testing only.

The tutorial workflow supplies a second hand-crafted synthetic Ensembl
Regulation feature and a manual sequence selection inline. These three inputs
exercise distinct selection methods without making a binding, affinity, or
causality claim.

The inline synthetic Ensembl-shaped overlap report uses deterministic stand-in
resource identities: its index digest is SHA-256 of the ASCII text
`synthetic_tutorial_index_v1`, and its interval digest is SHA-256 of
`synthetic_tutorial_intervals_v1`. They identify tutorial inputs only; they do
not claim that a provider index was downloaded.

## Deterministic recreation

1. Verify the source files and hashes described in
   `test_files/fixtures/promoter_reporter_architecture/README.md`.
2. Strip the FASTA header and newlines to recover the 18,106 bp forward-strand
   sequence for `GRCh38:7:101121158-101139263`.
3. Map each versioned transcript/exon ID in the geometry JSON to its stable ID
   plus numeric version, preserve every 1-based inclusive genomic range, and
   map `translation_start_1based` / `translation_end_1based` to the optional
   translation geometry in `gentle.ensembl_gene_entry.v1`.
4. Set import time to zero and retain the literal source geometry JSON in
   `raw_lookup_json`; no network retrieval is performed during tutorial runs.
5. Run `shasum -a 256` over the generated entry and BED file and compare the
   values recorded below after tutorial generation.

The generated entry SHA-256 is
`cfb894e8e1e878f36de0128f1e69e791671d6241661526fd015b2bad6a5ce46e`.
The BED SHA-256 is
`a65bcaf0048634a723896b0ce5c894634d95ad8671b388f305956301a54d9529`.
The BED file is recreated exactly as:

```text
7\t101126000\t101126120\tsynthetic_TP73_support\t950\t.
```

## GENtle use

`docs/examples/workflows/portable_genomic_regions_offline.json` loads this
entry without network access, saves manual, Ensembl-regulatory, and CUT&RUN
regions through the shared operations, exports lossless JSON plus BED6 and its
content-bound manifest, and renders the three regions on a SERPINE1 locus
figure. The engine and workflow tests verify coordinate conventions, source
provenance, explicit non-claims, and stable SVG markers.

The fixture validates contracts and rendering only. Synthetic signal is not
evidence for SERPINE1 regulation, TP73 occupancy, biochemical affinity, or
causal activity.
