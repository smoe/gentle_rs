# `tests/data`

Committed fixtures for `snapgene-reader`.

## `toy.small.dna`

- Origin: hand-crafted synthetic SnapGene `.dna` record that mirrors the
  existing GENtle `toy.small` import-parity sequence and annotations.
- Biological content:
  - 120 bp linear DNA sequence
  - `source` feature covering the full span
  - forward `gene` feature `toyA`
  - reverse `gene` feature `toyB`
  - one `misc_feature` note `barcode`
- Deterministic recreation:
  1. Start from the sequence and feature content in:
     - `test_files/fixtures/import_parity/toy.small.gb`
     - `test_files/fixtures/import_parity/toy.small.gbseq.xml`
  2. Encode a SnapGene packet stream consisting of:
     - cookie packet `0x09`
     - DNA packet `0x00`
     - features packet `0x0A`
     - notes packet `0x06`
  3. Use the same one-based feature ranges and synthetic description/accession
     metadata preserved in the committed fixture.
- Primary usage:
  - `packages/snapgene-reader/tests/reader.rs`
  - GENtle cross-format import tests in `src/dna_sequence.rs`
- Purpose: stable offline regression coverage for SnapGene `.dna` parsing and
  adapter parity without committing vendor-distributed example files.
