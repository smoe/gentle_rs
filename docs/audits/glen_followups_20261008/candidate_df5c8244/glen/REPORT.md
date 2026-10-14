# GENtle frozen `.13` verification at `df5c8244`

Date: 2026-10-14 (Europe/Berlin)

## Identity and isolation

- Source: `df5c8244787d8c920ee171ee55167de55e8a7294`
- Tree: `beefc526b11d4801fc57e2189bf0bc1ec18a439d`
- `Cargo.lock` SHA-256: `5893dd19e359d71a9a52ccaf688f952b2ef0507eda53dbdf2c3fdd60a775f84d`
- Platform: Debian Linux, `x86_64-unknown-linux-gnu`, kernel `7.2.7+deb14-amd64`
- Rust/Cargo: 1.96.1; `RUST_MIN_STACK` unset
- Effective local build settings: dev/default features, `CARGO_BUILD_JOBS=1`,
  `CARGO_INCREMENTAL=0`, `CARGO_PROFILE_DEV_DEBUG=0`, run-owned
  `CARGO_TARGET_DIR=/tmp/gentle-df5-audit-JuxMtC/target`
- Checkout: detached, clean, isolated worktree. The primary checkout and all
  historical evidence were left untouched.

## Repaired Rust boundary: PASS

The requested commands ran without `RUST_MIN_STACK` and without a diagnostic
replacement harness:

- `dbsnp_`: 11 passed, 0 failed
- `materialize_variant_allele`: 2 passed, 0 failed
- `reporter_construct_handoff_`: 7 passed, 0 failed
- `variant_followup_allele_pair_`: 2 passed, 0 failed
- `test_support::heavy_stack_tests`: 2 passed, 0 failed
- `release_version_consistency`: 4 passed, 0 failed
- `cargo check -q --locked`: passed

Total: 28 requested Rust tests passed. Both `dbsnp_spdi_*` regressions are in
the 11-test dbSNP filter. The cold default test harness required 17:07 wall and
10,732,312 KiB peak RSS; this is build cost, not test execution time.

## `.13` version and tutorial projections: PASS

- `gentle_examples_docs` reports `0.1.0-internal.13` and exact source
  `df5c8244787d8c920ee171ee55167de55e8a7294`.
- `gentle_examples_docs` SHA-256:
  `f5c5f99c3445ebfd700d2b28be527170d994b82ba9154097eb72f54cf407c3df`
- GUI binary SHA-256:
  `5575196c75653f58429007ccad8a430195d45af4715658bde2f021f334f342ae`
- `scripts.test_release_candidate`: 44 passed
- `scripts.test_synthetic_biology_tutorial`: 4 passed
- `scripts.test_frozen_024ed6_audit`: 5 passed
- `scripts.test_tutorial_checkouts`: 56 passed
- Helper `--check`, catalog, manifest and tutorial checks: passed
- Catalog: 66 entries; executable manifest: 30 chapters; workflows: 53
- LF and CRLF checkout replay: passed
- Generated catalog/manifest/tutorial bytes and `Cargo.lock`: no drift
- `git diff --check`: passed

Two expected historical review-hygiene warnings remain for 02.01 and 08.03.
They are not `.13` projection failures and were not rewritten.

## Remaining screenshot queue

All attempts used only retained public/synthetic inputs, a disposable copy,
the exact-source GUI binary, Xvfb/Openbox at 1600x1200, and no disabled capture
bridge or permission bypass.

### FLNA: PASS on virtual Linux/X11; physical platforms NOT RUN

The retained public project SHA-256
`88b40cf12d73d7ad1ca88a8d0fc59f9931abcf812bc0b869f81c47993d057452`
opened in the exact `.13` GUI. Raw captures retain:

- project/sequence overview;
- GRCh38 genomic-minus anchor versus local-plus splicing orientation;
- 36 transcripts and 103 exons;
- canonical `ENST00000369850` focus;
- both loaded UniProt-referenced isoform rows;
- the dedicated splicing-detail window at normal and 1500x1100 size.

This closes the requested Linux virtual-X11 presentation check only. It does
not provide historical annotation-release provenance, scientific confirmation,
or physical Windows/macOS acceptance.

### 07.04 microarray: BLOCKED in this environment

The retained 83,686-bp public locus project loaded and its lineage overview was
captured. Opening the sequence window stayed CPU-active for more than two
minutes (about 1.3 GiB RSS) while the X server ceased responding to bounded
window queries. The process was terminated; no source state was overwritten.
Native navigation/results/save-reopen and readable result capture remain open.

### 10.03 measured gel: BLOCKED in this environment

The retained synthetic project loaded, and the genuine
`Patterns -> Gel Image Analysis...` menu route was exercised. Opening the
editor stayed CPU-active for more than two minutes while X11 did not respond.
The process was terminated. Calibrated/null result display and GUI save/reopen
remain open; the existing headless SVG is not substituted as native evidence.

### 08.16 imported TSS: BLOCKED

The retained 11-bp GenBank export and complete/partial report bytes were used
without any DuckDB replay. The project GUI opened, but this bounded Xvfb route
did not establish a report-bound TSS viewer: the positional file argument did
not import the sequence, and the subsequent manual file-menu attempt opened an
unrelated specialist window rather than a validated sequence/report adoption.
No lane screenshot is claimed. Complete/partial imported-lane readability
remains open.

### 04.09, 06.07, 03.05 and 10.02: NOT RUN

No exact-source, pre-reviewed disposable GUI state exists in the retained
evidence for these interactions. Creating/applying new sequence-design,
Gibson, sequencing-confirmation or full primer-import state would broaden this
frozen verification and, for 04.09, violate the explicit performance-fix
prerequisite. They remain in `docs/tutorial/coverage_backlog.md`.

## CI and final identity

Public GitHub run 38541824374 completed successfully at the exact source SHA.
Windows, macOS, Linux/headless, Release-Safety and the selected-platform gate
passed. Optional full Linux, synthetic GUI, GUI-scroll and authenticated macOS
jobs were skipped and are not silently accepted.

A final fetch left `origin/main` at the same exact SHA. No push, workflow
dispatch, tag, publication, model call or package build was performed. No
`package-opt1` evidence is relabelled as `package-opt2` acceptance.

