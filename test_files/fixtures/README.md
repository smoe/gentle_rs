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
  - `toy.small.embl`
  - `toy.multi.embl`
  - `toy.small.gbseq.xml`
- `resources/`
  - `jaspar.edge.pfm`
  - `rebase.edge.withrefm`
- `primer3/`
  - `pairs.location_5_60.kv`
- `mapping/`
  - `ensembl_chimp_tp73_all.fasta`
  - `ensembl_human_tp73_all.fasta`
  - `ensembl_human_tp53_all.fasta`
  - `ensembl_mouse_trp73_all.fasta`

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
  EMBL, and NCBI GenBank XML (`GBSet/GBSeq`) forms.
- Primary usage:
  - Cross-format import/parity tests.
  - EMBL parser location/qualifier regression coverage.
  - XML import normalization tests against current GenBank baseline.
- Purpose: tiny, reviewable fixture set for format-difference debugging.

### `import_parity/toy.multi.embl`

- Origin: hand-crafted two-record EMBL fixture.
- Primary usage:
  - Multi-record EMBL ingestion tests.
  - Wrapped feature-location and wrapped qualifier regression coverage.
- Purpose: ensure parser behavior stays stable for realistic EMBL line wrapping
  and record boundaries.

### `resources/jaspar.edge.pfm`

- Origin: hand-crafted edge fixture from commit
  `2154eb966409b91f714c2139113c907028c42634`.
- Primary usage:
  - `src/resource_sync.rs` parser tests.
  - `src/engine_shell.rs` shared-shell execution tests for
    `resources sync-jaspar` + motif-registry reload side effects.
  - Runtime paths `resources sync-jaspar` / GUI import dialog.

### `resources/rebase.edge.withrefm`

- Origin: hand-crafted edge fixture from commit
  `2154eb966409b91f714c2139113c907028c42634`.
- Primary usage:
  - `src/resource_sync.rs` parser tests.
  - `src/engine_shell.rs` shared-shell execution tests for
    `resources sync-rebase`.
  - Runtime paths `resources sync-rebase` / GUI import dialog.

### `primer3/pairs.location_5_60.kv`

- Origin: hand-crafted synthetic Primer3 key/value response fixture.
- Deterministic self-creation:
  - create a plain-text file with:
    - `PRIMER_PAIR_NUM_RETURNED=1`
    - `PRIMER_LEFT_0=5,20`
    - `PRIMER_RIGHT_0=79,20`
    - `PRIMER_PAIR_0_PRODUCT_SIZE=75`
    - terminal `=`
- Primary usage:
  - `src/engine_shell.rs` adapter-equivalence tests for
    internal-vs-Primer3 primer report normalization using a fake local
    `primer3_core` script.
- Purpose: deterministic offline coverage for Primer3 normalization/provenance
  behavior without requiring a system Primer3 installation in CI.

### `mapping/*.fasta` (TP73/TP53 benchmark set)

- Origin: copied from the legacy TP73 mapping corpus in
  `test_files/mapping/True_TP73/` to provide a small, deterministic benchmark pack
  under committed fixtures.
- Current provenance status:
  - biological source is Ensembl transcript FASTA exports (human/chimp TP73,
    human TP53, mouse Trp73), as indicated by headers/filenames.
  - exact original export URLs/queries are not yet frozen in-repo; this is a
    known provenance follow-up item in `docs/roadmap.md`.
- Primary usage:
  - `src/engine.rs` seed-filter regression:
    `test_tp73_seed_filter_cross_species_and_tp53_specificity_sets`.
- Purpose:
  - keep CI/stable tests on compact mapping fixtures,
  - include TP53 as a close-family negative benchmark (as requested for better
    practical signal than TP73-only stress sets).

## Large exploratory XML samples

- Large one-off XML records are intentionally excluded from default committed
  fixtures to keep repository size controlled.
- If needed for local parser experiments, fetch/regenerate into a local
  scratch path and do not commit by default.
