# `test_files` overview

This directory contains test/demo assets. Committed deterministic fixtures now
live under `test_files/fixtures/`.

## Key paths

- `fixtures/`: canonical committed fixtures used by tests and parser/runtime
  checks.
  - See `test_files/fixtures/README.md` for provenance and per-file usage.
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

## Policy

- Prefer adding new deterministic fixtures to `test_files/fixtures/` with clear
  provenance notes.
- Avoid committing large one-off downloads unless they are required for a
  stable CI/test contract.
