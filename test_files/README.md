# `test_files` overview

This directory contains test/demo assets. Committed deterministic fixtures now
live under `test_files/fixtures/`.

## Key paths

- `fixtures/`: canonical committed fixtures used by tests and parser/runtime
  checks.
  - See `test_files/fixtures/README.md` for provenance and per-file usage.
  - Includes `fixtures/mapping/` for compact TP73/TP53 RNA-mapping
    benchmarks used in deterministic seed-filter tests.
- `pGEX-3X.gb`, `pGEX_3X.fa`, `tp73.ncbi.gb`:
  - historical sequence fixtures still referenced by existing tests/examples.
- `pGEX-3X.embl`:
  - EMBL export of ENA accession `U13852` (`pGEX-3X cloning vector, complete
    sequence`) from [https://www.ebi.ac.uk/ena/browser/view/U13852](https://www.ebi.ac.uk/ena/browser/view/U13852).
  - used for EMBL parser parity tests against `pGEX-3X.gb` in
    `src/dna_sequence.rs`.
- `MA1234.1.jaspar`:
  - minimal motif fixture for focused parser behavior.
- `demo_nonsensical.state.json`, `project.gentle.json`:
  - large demo/state snapshots for manual exploration and regression scenarios.
- `cloning_digest_ligation_extract.gsh`:
  - shell workflow example.
- `bioseq.rs`:
  - standalone XML parsing experiment source (not part of main runtime import
    path).
- `mapping/True_TP73/`, `mapping/False_TP73/`:
  - legacy TP73 mapping corpus and decoy sets used for exploratory runs.
  - these directories may contain larger local-only payloads and are not the
  preferred location for stable CI fixtures.
  - stable regression coverage should use `test_files/fixtures/mapping/`.

## Policy

- Prefer adding new deterministic fixtures to `test_files/fixtures/` with clear
  provenance notes.
- Avoid committing large one-off downloads unless they are required for a
  stable CI/test contract.
- For mapping benchmarks, prefer compact curated sets in
  `test_files/fixtures/mapping/` over full exploratory transcriptome dumps.
