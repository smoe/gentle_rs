# Consolidated Glen audit — exact candidate `024ed6c0`

## Scope and identity

- Tested source: `024ed6c0520d38c4221517219d44c792084de4e5`
- Final fetch after all measurements still resolved `origin/main` to this exact
  SHA; no tip drift occurred.
- Source tree: `3982cc6059530c9caa21ba61e2fd63f4405bf0cc`
- `Cargo.lock` SHA-256:
  `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`
- Host: Debian Linux x86_64, kernel `7.2.7+deb14-amd64`; Rust/Cargo
  `1.96.1`.
- Matching dev CLI SHA-256:
  `bd27733eced2b811a9835f56d9d59b05a28556d1d9cfdfde626e942209b04cb6`.
- Matching GUI binaries: dev
  `635620123bb64a492869eab82070031bf96a53fb77197a22cd9d253f96242a7a`,
  package-opt1
  `ab71e0c1c96e56f81276de1e37ec8fee6f9928d4d71da9265fad64240b0bd541`.
- Normal GitHub Actions run
  `38453643146` completed successfully at this exact SHA. Windows, macOS,
  Linux/headless and release-candidate safety jobs passed. Optional Linux
  aggregate, scroll, authenticated macOS and synthetic-GUI jobs were skipped
  according to workflow selection and are not silently counted as executions.
- The primary checkout, private projects and historical evidence were not
  changed. No push, dispatch, tag, publication or source repair was performed.
  All three clean detached audit worktrees were removed after the final fetch;
  the run-owned dev target (3.1 GiB) and package target (3.1 GiB) were cleaned.

## Findings, risk ordered

### 1. Medium — the requested Linux focused commands abort on default test-thread stack

The four requested original commands all selected nonzero tests but exited 101
after a test-thread stack overflow:

- `cargo test --locked --lib dbsnp_ -- --test-threads=1`: 11 selected; abort in
  `app::tests::dbsnp_dialog_fetch_extracts_region_and_opens_window`.
- `cargo test --locked --lib materialize_variant_allele -- --test-threads=1`:
  2 selected; abort in the first test.
- `cargo test --locked --lib reporter_construct_handoff_ -- --test-threads=1`:
  7 selected; five passed, then abort.
- `cargo test --locked --lib variant_followup_allele_pair_ -- --test-threads=1`:
  2 selected; abort in the first test.

This is a real failure of the requested Linux execution boundary. It is not
reported as a VKORC1 logic regression: diagnostic reruns with
`RUST_MIN_STACK=16777216` passed 11/11, 2/2, 7/7 and 2/2 respectively, and the
normal cross-platform CI is green. Those diagnostic runs do not replace the
failed original commands. A separate focused `dbsnp_spdi_` run passed 2/2, and
`cargo check -q --locked` passed. A future repair should be bounded to the test
harness/large stack use and rerun the original commands; no parser broadening
is justified by this evidence.

### 2. No real-reference VKORC1/dbSNP defect found

The retained NCBI refSNP response is 75,896 bytes with SHA-256
`080131663dad501faf0101fc1fd50b842b64d000d17c17a4222c654f4c16052a`.
It selects `NC_000016.10`, GRCh38.p14, zero-based position 31,096,367,
reference C and alternates A/G/T. Independent `samtools faidx` lookup of
`16:31096368` in the prepared Ensembl-116 GRCh38 FASTA returned C.

The exact-source CLI imported the response through a fixed local `file://`
provider into a fresh state in 24.54 s (484,580 KiB peak RSS). The 6,001-bp
context retained a matching reference/assembly, no fallback, 113 projected
annotations and the marker. Prepared payload SHA-1 values were
`294c432b69a88521437ffe2a9080e08d3f84902c` (sequence) and
`1cc995089649e7d6bc291b12947f2fe98989975f` (annotation), unchanged before and
after the run.

The shared-shell path selected a 1,089-bp minus-strand promoter fragment.
Ambiguous materialization was refused with A/G/T listed and byte-identical
state. Explicit genomic-forward T produced exactly one difference at local
zero-based position 588, C→T; fragment/reference/alternate lengths and feature
counts remained 1,089/20. The minus-strand gene did not silently complement T.
Assembly mismatch, unavailable reference and reference-base mismatch blockers
are covered by the focused contracts; their logic passed under the explicit
16-MiB diagnostic stack and in normal CI, while the original Linux command
boundary remains failed as described above.

### 3. Native acceptance remains partly blocked or unrun

- **08.04 physical/native Linux tutorial: BLOCKED.** The approved procedure
  requires the repository's isolated `sudo unshare` path; `sudo -n` requires a
  password on this host. No permission bypass was attempted. Consequently
  native explicit-selection/refusal/retry/rollback/save-reopen/reporter UI,
  raw captures and the deliberate two-Undo GUI observation remain unaccepted.
  The CLI sequence and full records are retained, but a new CLI process has no
  GUI-session undo history and is not a substitute.
- **Windows GenBank contract: PASS (automated native CI).** The exact-SHA
  Windows full suite passed. Both unguarded `dbsnp_spdi_*` regressions therefore
  ran on Windows and cover actual GenBank export/re-import, LF/CRLF, wrapped
  JSON, allele order and empty/long values, unchanged reference/assembly and
  absence of fabricated VCF alleles. Manual Windows GUI and extracted-package
  acceptance were not run and remain in the central release gate.
- No supported marker-import route was found demonstrably bypassing SPDI
  restoration. No speculative normalization change is proposed.

## Per-item status

1. Focused commands: **FAIL** for the four original test filters (stack
   overflow); **PASS** for `cargo check`; **PASS diagnostic only** with enlarged
   stack; `dbsnp_spdi_` **PASS 2/2**.
2. Actual rs9923231/GRCh38 replay: **PASS** for retained bytes, placement,
   independent reference base, reference binding, explicit T, ambiguity refusal
   and unchanged source/reference data.
3. Native tutorial 08.04: **BLOCKED** by approved-host isolation permission.
4. Native Windows GenBank round trip: **PASS** for exact-SHA automated native
   Windows tests; manual GUI/package layer **NOT RUN**.
5. Parser-scope review: **PASS / no demonstrated additional route defect**;
   no broadening proposed.

## Help/startup comparisons

### A. Lazy discovery — PASS as a separate comparison

Exact dev builds of `e723f46eb3fe58b15ba4bad5ff975024cd761437`
(binary SHA-256 `ec9e2886…bb74a3`) and
`64baba4f1cd21c9df731b232d9f82ce55c9d6495`
(`346e7b03…21f1d`) were measured with three cold and three warm included runs
each. Cold usable-main median improved from 34.457 s to 2.672 s; warm remained
2.650 s versus 2.644 s. The subject and first Help/Topic/menu measurements show
no startup saving transferred to first use. Full report:
`gentle-help-startup-audit-20261013/REPORT-A.md` (SHA-256
`6d7e27f0838661c880c9b081a2c8a807a1e7739ae286d72326292e6cf17ba831`).

### B. Build profile — PASS as a separate comparison

At exact `024ed6c0`, package-opt1 improves usable-main median versus dev by
256.865 ms cold (9.7%) and 270.286 ms warm (10.1%); first subject-correct locus
improves by about 0.93–0.95 s (25%). Tutorial-window latency fell from about
5.2 s to about 0.25 s in both cache conditions. This benefit is not deferred
to first Help use. All 12 included runs completed with 52 events and no losses,
failed spans or image failures; the two cache-filling runs are retained and
programmatically excluded. Full report: `REPORT-B.md`.

Neither A nor B establishes physical-desktop or release acceptance. The data
do not justify shared-font, failed-SVG or menu caching changes: no relevant
failure was observed and the measured menu-discovery span was already tiny.

## Evidence map

- Exact commands, exits and logs: `logs/focused/`
- Real refSNP and reference replay: `vkorc1/raw/`, `vkorc1/real-reference/`
- CI run/job receipts: `ci/run-final.json`, `ci/jobs-final.json`
- B raw traces, screenshots and statistics: `help-startup-B/`
- Exact binaries and build settings: `builds/`, `profile-definitions.txt`,
  `cargo-metadata-no-deps.json`
- A evidence (kept separate):
  `/home/clawbio/.openclaw/workspace/artifacts/gentle-help-startup-audit-20261013`
- B/VKORC1 evidence root:
  `/home/clawbio/.openclaw/workspace/artifacts/gentle-consolidated-audit-024ed6c0-20261013`
