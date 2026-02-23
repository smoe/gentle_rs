# `test_files/fixtures`

Structured, committed fixtures used by parser/runtime tests and catalog smoke
checks.

## Layout

- `genomes/`
  - `AB011549.2.fa`
  - `AB011549.2.gb`
- `import_parity/`
  - `toy.small.fa`
  - `toy.small.gb`
  - `toy.small.gbseq.xml`
- `resources/`
  - `jaspar.edge.pfm`
  - `rebase.edge.withrefm`

## Provenance and usage

### `genomes/AB011549.2.fa` + `genomes/AB011549.2.gb`

- Origin: E. coli reference sequence/annotation pair used as a local catalog
  fixture.
- Primary usage:
  - `assets/genomes.json` `LocalProject` entry.
  - Genome catalog validation/smoke checks in CI.
- Purpose: deterministic local test genome that avoids remote network
  dependencies.

### `import_parity/toy.small.*`

- Origin: hand-crafted synthetic 120 bp sequence represented in FASTA, GenBank,
  and NCBI GenBank XML (`GBSet/GBSeq`) forms.
- Primary usage:
  - Cross-format import/parity tests.
  - Planned XML import normalization tests against current GenBank baseline.
- Purpose: tiny, reviewable fixture set for format-difference debugging.

### `resources/jaspar.edge.pfm`

- Origin: hand-crafted edge fixture from commit
  `2154eb966409b91f714c2139113c907028c42634`.
- Primary usage:
  - `src/resource_sync.rs` parser tests.
  - Runtime paths `resources sync-jaspar` / GUI import dialog.

### `resources/rebase.edge.withrefm`

- Origin: hand-crafted edge fixture from commit
  `2154eb966409b91f714c2139113c907028c42634`.
- Primary usage:
  - `src/resource_sync.rs` parser tests.
  - Runtime paths `resources sync-rebase` / GUI import dialog.

## Large exploratory XML samples

- Large one-off XML records are intentionally excluded from default committed
  fixtures to keep repository size controlled.
- If needed for local parser experiments, fetch/regenerate into a local
  scratch path and do not commit by default.
