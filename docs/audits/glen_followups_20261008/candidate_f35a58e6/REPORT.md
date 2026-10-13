# GENtle final-candidate execution report

## Scope and identity

- Authorized candidate: `origin/main` at `f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`.
- The `v0.1.0-internal.12` tag was explicitly excluded by the owner as a GitHub exercise.
- Final fetch at 2026-10-13T04:59:26+02:00 still resolved `origin/main` to the candidate with no diff.
- Git tree: `f678fa505010119b0b546db116a82438e07f6405`.
- Raw `Cargo.lock` SHA-256: `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
- Clean detached worktree; the dirty primary checkout was not changed.
- Host: Debian Linux x86_64, Rust/Cargo 1.96.1, Python 3.14.8, DuckDB 1.5.5.
- Profiles: ordinary candidate binaries used Cargo `dev`/default features; the focused 08.04 and scroll audits used CI's `dev`, `CARGO_BUILD_JOBS=1`, `CARGO_INCREMENTAL=0`, `CARGO_PROFILE_DEV_DEBUG=0` settings.
- Original default-feature binaries and hashes are in `logs/binary-identities.txt`. Focused 08.04 binaries are independently recorded in `artifacts/vkorc1-gui-audit/receipt.json`.

## Candidate checks

### 08.16 imported TSS evidence — PASS for source/shared commands; native NOT RUN

- Python fixture and the ignored real-DuckDB Rust regression passed, with exactly one Rust test executed.
- Learner queries returned 3 full-minus, 3 narrower-minus and 2 plus hits.
- Imported lanes preserve both strands, signed values and a genuine zero; the clipped footprint is retained in the export.
- Fresh export contains 12 hash-bound outputs. `report.json` SHA-256 is `104f5ad4ea7c99b74c5d97de0aef1dbc6389e96a7d61dc554293f242f222a0c9`.
- No ordinary native learner session or human readability review was available on this headless host.

### 08.18 repeat cohorts — PASS for software contract; native/biology NOT CLAIMED

- Rust regression passed with one actual test.
- `LINE/L1` selected two distinct AB011549.2 loci, both containing genomic-forward `AGCCAG` and opposite repeat strands.
- Two midpoint windows have nonempty score tracks and genomic peak positions.
- All three transcript-relative modes keep both loci and explicitly report `no transcript projection available`; RNA ranks remain unset and support lists empty.
- The bacterial sequence is genuine; repeat labels are synthetic. No dedicated cohort GUI, human repeat biology or RNA interpretation is claimed.

### 07.04 microarray and 10.03 gel — PASS for source/shared commands; native NOT RUN

- Microarray regression passed. Learner import projected five TP73 features with no skipped rows while retaining the unverified-local-genome and synthetic-statistics warnings.
- Gel fixture check and 26 Rust gel tests passed. Learner DNA result is 707.10678 bp; protein result is 35.3553 kDa. Both outside-calibration results remain null with `outside_calibrated_range`; TSV/SVG exports were retained.
- State files reopened across separate CLI invocations. Native configuration/results, navigation and screenshot readability were not run.

### FLNA — PASS for retained-state integrity and 2400-px rendering; BLOCKED for activation

- Stable selector `ENST00000369850.10` was rediscovered at feature index 300 from the state, not assumed.
- Current candidate reproduces 36 transcripts, 103 unique exons, 73 junctions, six event summaries and presentation fingerprint `sha256:742aa8462004ec391e38c739c4e1fe0addd272255ec142c8123dc0336f72b5e2`.
- Current SVG is byte-identical to the retained replay (`1828fdd8e07d21e79feab21519eccaa6aa4bbf43610784acf72db3705855e8de`). The inspected 2400x14808 PNG has SHA-256 `067693d140868689b3295493f2143e1e8873a7b49a2792aaa0c9f00e56cc9ec0`.
- Retained acquisition metadata, provider responses and notes contain neither a historical Ensembl annotation-release record nor a contemporaneous `/info/data` payload. No modern release was back-attributed.
- Human teaching-exon selection and scientific consequence interpretation were not supplied. Normal-size native GUI review was not run. The draft must remain inactive.

### 04.09 primer specificity — PASS for authored gates and external BLAST; PERFORMANCE-BLOCKED at real import

- Every command in the chapter's verification scope passed: 4 reference tests, 11 multi-reference tests, three focused single tests, 20 capability-registry tests and one parity-freshness test, plus candidate docs and Python checks.
- Real BLAST 2.17 ran all four argument-vector commands against prepared Ensembl 116 cDNA and GRCh38 references with exit code 0 and hash-bound manifests.
- Reference identities: cDNA fingerprint `sha256:3159519891bd8c0d4690e9c2f3de8eeb9621152b58c749e3d7d4920f4913263f` (465,769 sequences, 1,276,790,359 bp; subject-annotation SHA-1 `5f571d5775368aeb4e8a88c3ce4b1c8c138e3acf`); GRCh38 fingerprint `sha256:90a64ae8ff24fdbab67286bab563c95a6b2a04be37f50ad811474106c1ba3e7e` (706 contigs, 3,291,585,349 bp).
- Default ACTB pair produced 5,999 TSV HSP lines / 686,525 bytes. A stricter-policy replay produced the same raw BLAST volume. A 30-mer TP53 pair produced 12,202 HSP lines / 1,404,773 bytes.
- In all three cases, import remained CPU-bound; the first bounded attempt was stopped after about ten minutes and the two diagnostic variants after shorter confirmation runs. No state/report was committed by the interrupted imports. This is a concrete candidate performance limitation in post-BLAST import/evaluation, not a BLAST, reference, provenance or process-contract failure.
- Required/missing/stale/incomplete/history/save-reopen behavior is covered by the passing deterministic regressions. Native specialist presentation, native Windows filesystem behavior, scientific order acceptance and a completed real-reference verdict remain open.

## Visual, provenance and general regression evidence

- All 19 visual-placement decisions and all referenced assets were found. The exact-source checkout/freshness tests passed: 12 accepted placements, five retained in 02.03, two deliberately omitted.
- Public E-MTAB-15709 SDRF SHA-256 is `60a3cdd3b5ba8c3e2b7f3fa823f3911dbf5091c61b90a31dd453bc2efad1377d`; all six rows for the three named `pos_*_SKMel29` libraries are H3K4me3.
- Source hygiene, formatting and `cargo check --locked` passed.
- Exact-source tutorial catalog (65 entries), executable manifest (30 chapters / six GUI contracts), tutorial freshness and LF/CRLF checkout replay passed. `tutorial-check` retained only the two pre-existing stale-human-review warnings.
- Python: 60 scoped source/checkout tests plus 268 six-gate helper tests passed.
- Rust: 68 current tutorial/primer/parity tests and 19 focused 08.04 tests passed with nonzero counts.

## Six-gate ledger

1. Audit/source preservation: **PASS** for source, helper and freshness checks.
2. 08.04 synthetic contract: **PASS** for 19 focused Rust tests and regenerated projections; **BLOCKED** for the isolated X11 input replay because local `sudo unshare` required an interactive password. No network-isolation bypass was attempted.
3. Genuine live-agent acceptance: **NOT RUN**; no authorized authenticated session was supplied.
4. Honest visual evidence: **PASS** for placement/source contracts; **NOT RUN** for fresh ordinary native captures.
5. Scroll event/next-frame smoke: **PASS**. Exactly 24 observations cover three fixtures, four sizes and both directions; sequence/annotation content stayed unchanged. This is a CPU determinism contract, not native latency.
6. Non-publishing `.12` packages: **NOT RUN**. No workflow dispatch/package authorization was supplied; no package or container acceptance is inferred.

GitHub Actions run `38329717079` is bound to this exact SHA. At report assembly Linux/headless, Windows and release-safety jobs were successful; macOS was still running. Optional follow-up jobs were skipped by the push workflow and were not counted as passes.

## Environment limits and interrupted diagnostic

- No ordinary `DISPLAY`, Wayland session or consented screenshot bridge was available.
- A redundant default-feature six-gate harness build approached the temporary filesystem limit before running; it was stopped and only this run's 7.5-GiB incremental cache was removed. Its relevant primer tests had already passed in the focused 11-test run. The event is retained in the handoff rather than represented as a test pass.
- Final run-owned Cargo cleanup removed 8,740 files / 16.1 GiB from the external target directory.
- macOS native interaction, native Windows GUI/filesystem behavior, genuine live-agent use, scientific review and `package-opt1` packaging remain untested unless separately evidenced above by GitHub Actions.

## Mutation boundary

No source/document edits, commits, pushes, tags, workflow dispatches, releases, image publication or changes to the primary checkout were made. All new files are disposable evidence outside the tested source tree.
