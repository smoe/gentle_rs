# GENtle Tutorial (Generated)

This folder is generated from:
- `docs/tutorial/manifest.json`
- `docs/examples/workflows/*.json`

Regenerate with:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
```

Validate committed generated output:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
```

Online execution was disabled (`GENTLE_TEST_ONLINE=0` during generation).

## Chapters

- 1. [Load FASTA, branch, and reverse-complement](./chapters/01_load_branch_reverse_complement_pgex_fasta.md) - `core` - example `load_branch_reverse_complement_pgex_fasta` - executed `yes`
- 2. [Load pGEX and digest with BamHI/EcoRI](./chapters/02_load_and_digest_pgex.md) - `core` - example `load_and_digest_pgex` - executed `yes`
- 3. [Guide practical filtering and oligo generation](./chapters/03_guides_filter_and_generate_oligos.md) - `core` - example `guides_filter_and_generate_oligos` - executed `yes`
- 4. [Digest -> Ligation -> ExtractRegion minimal slice](./chapters/04_digest_ligation_extract_region_minimal.md) - `core` - example `digest_ligation_extract_region_minimal` - executed `yes`
- 5. [Guide oligo export (CSV + protocol)](./chapters/05_guides_export_csv_and_protocol.md) - `advanced` - example `guides_export_csv_and_protocol` - executed `yes`
- 6. [Prepare a reference genome cache (online)](./chapters/06_prepare_reference_genome_online.md) - `online` - example `prepare_reference_genome_online` - executed `no`

## Source Summary

- Tutorial schema: `gentle.tutorial_manifest.v1`
- Chapter count: `6`
- Generation report: [`report.json`](./report.json)
