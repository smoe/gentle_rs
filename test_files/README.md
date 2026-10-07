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
  - `tp73.ncbi.gb` is the public RefSeq `NC_000001.11` forward-strand excerpt
    at chr1:3652516..3736201, GRCh38.p14 (`GCF_000001405.40`); its record date is
    26-AUG-2024. The deterministic retained fixture is Git blob
    `8a17623518e659cfc2264351af982a6ae0e1971f`, SHA-256
    `9acebbd259f535a9f780006632e2fde3d796d4674783b313d5b28930ff898ab4`.
    Reproduce its bytes with `git show 8a17623518e659cfc2264351af982a6ae0e1971f`.
    The original reference interval is retrievable from NCBI EFetch with
    `db=nuccore&id=NC_000001.11&seq_start=3652516&seq_stop=3736201&strand=1&rettype=gbwithparts&retmode=text`;
    annotation updates may differ, so a fresh retrieval is not a replacement
    for the pinned test input.
  - `test_tp73_dnp73beta_terminal_exon_skip_preserves_annotated_cds_through_shell`
    in `src/engine/tests.rs` uses the runtime GenBank loader and shared shell
    planner/materializer. It resolves `NM_001126241.3` to `NP_001119713.1`
    without an explicit CDS transcript ID, checks the 1,353-bp full CDS and
    16-bp terminal coding contribution, and retains 1,337 coding bp within
    1,571-bp exon-skipped cDNA. It checks partial-codon/stop-loss warnings and
    saved annotation re-derivation, not wet-lab suitability or genome authority.
  - Simple-PCR GUI acceptance uses `ExtractRegion` on the committed
    `tp73.ncbi.gb` interval [61520, 62320) (0-based, end-exclusive), retaining
    the original locus and its projected annotations. Recreate it with
    `docs/examples/workflows/simple_pcr_selection_gui.json`; the scripted
    `simple_pcr_primer_design_offline.json` oracle extracts identical bases.
    The fixed ROI is [200, 600) of that 800-base template. This teaching
    fixture is not a validated assay or a whole-genome specificity proof.
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
