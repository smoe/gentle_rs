# Glen: Next Verification And Return

## C2-C4 Conditional cDNA Interpretation, 2026-11-12

**New active scope after the owner's push:** the integrated source containing
`InterpretCdnaAlleles`, not a repeat of prior reader/BLAST, cDNA/locus, barcode,
packaging or real-data campaigns. C1 remains a separate artifact; comparison and
annotation-v1 wires/counts/policies must remain unchanged. See the
[contract](cdna_allele_interpretation.md) and [C1-C4 ledger](cdna_allele_interpretation_plan.md).
Record exact tested SHA, lock/binary hashes, toolchain, platform/features/profile,
commands, nonzero counts, failures and explicit skips. No private-data run,
provider call, push, dispatch or publication is requested.

### Glen Return At 2f1a3062

Exact source `2f1a30620ae04cba2c4def4f29fb46279b0b8f84` is **not accepted for
C2-C4**. Glen's Debian Linux/Rust-Cargo 1.97.1 run passed both compile gates,
canonical generation (zero diff bytes), the C1 downstream work-ceiling test
(1/1), glossary (16), registry/freshness (32 + 1), checkout/maturity (116) and
three source-level companion checks. Preserve those scoped results; the work-
ceiling test is not a comparison-byte or C2-C4 runtime pass.

Two primary boundaries failed independently:

- Six headless tests and both actual CLI cases stopped during C1 preparation
  with `coding_models.annotation_binding`. The script-interface build selected
  JS/Lua tests but failed at the same fixture before their VM route; MCP also
  never reached its intended assertions. Isolated reruns reproduced the refusal,
  so this is not mutex-poisoning fallout and not a VM pass or VM-engine defect.
- The `gentle-engine` test target failed E0282 at its unconstrained tuple key.
  Root `cargo check --tests` did not prove that separate crate's test target
  compiled. Its kernel tests and the historical old/new byte comparison remain
  blocked/not run, not passed.

Repair source now uses C1's existing annotation-reference selector namespace and
keeps the saved-analysis source distinct. C2 no longer equates those unrelated
IDs; shared raw annotation, exact gene, reference, model and C3 input-hash guards
remain. New guards exercise wrong selectors and mismatched bindings. The tuple
has an explicit `(String, String, String)` type. **These repairs are not yet
Rust/VM/process-verified.** No C1 schema or producer files were rewritten.

The supplied archive was checked without extraction/execution: all **22/22**
manifest entries match and its report is byte-identical to the standalone file.
Archive SHA-256: `b6a8677c9d87611684b2c90f04c31b98c887f9af3436d77c4d3034b241baf087`;
report SHA-256: `258b681580ae6d54eb24b20671b98b93e5dcabb8d3c928d908873709577b09af`.
The failed-source CLI hash was
`da95aebbc2173be7b28373c27f4b168f3c0a34a2357516678b2b0e32bcf6150d`.
CI `43223393601` was still running at Glen's closing observation; no terminal
or native-host verdict is inferred here. Glen's separate cleanup/archived unused
Nanopore import is not part of this repair.

### Existing Evidence, Not New-Source Acceptance

Two independent Glen repositories/runs must remain distinct. The owner-supplied
`GENtle-4966a3a5-verification-report.md` records a successful clean Linux run at
**`4966a3a514c1765b2e931c4b51cc078037e65634`**. The later routine scan's ENOSPC and
failed status update do not invalidate that pass or make completed checks unrun.

Glen's clean run used Debian Linux x86_64, Rust/Cargo 1.97.1, a clean detached
checkout and an explicitly run-owned target. Both locked compile gates passed;
canonical matrix generation produced no tracked diff. C1 headless: 14 passes;
GTF: 16; C1 with `script-interfaces`: 16, including one actual JS VM and one actual
Lua VM case. MCP confirmation-before-I/O passed in both profiles; glossary 16,
registry 31 and freshness one passed. Counts overlap across profiles: do not
sum them as independent tests.

The source-matched C1 CLI companion passed **4/4, zero skips**, including both
actual-process cases. A separate repeated direct-route probe was byte-identical,
kept inputs unchanged and created no project. All eight Read-End filter tests
passed, including explicitly configured `/usr/bin/samtools` 1.24 for indexed BAM
and reference-backed CRAM (no `NOT RUN`). C1 VM/CLI and real-samtools Linux
acceptance are therefore recorded, not pending historical campaigns.

- Report SHA-256 (the supplied Markdown only):
  `d7b592b24b2b9cf941bf5354d6c6f8bac4e648e33a6f8a83b54aedebbda7988c`.
- CLI SHA-256:
  `de8d5c29e1031363192198be800d9c9fcae962d91d72713abf5d2661f9723879`.
- Lockfile SHA-256:
  `d880f70bcb3bb8ebd8db53e0ef7a8456a060c841a891f97f8d4a2eaf159569fe`.

This is Glen-attributed evidence: the supplied file describes sealed raw logs,
but those logs were not independently inspected here. Its closing CI observation
was still in progress; retain the separate later terminal CI receipt below.

At the same SHA, normal
[CI 43122283124](https://github.com/smoe/gentle_rs/actions/runs/43122283124)
finished green for Linux, Windows, Linux/headless and release safety; macOS
skipped. The Linux log includes C1/core/GTF/shell/MCP and the seven nonoptional
Read-End fixes, 5,115 library passes (17 ignored), 31 registry and one freshness
pass. Glen's separate ENOSPC run reached no focused tests and failed to advance
its status file: record it as incomplete, not a product failure or extra pass.
Do not rebuild those old gates merely to repair the routine scan's missing local
receipt. These successful clean-run and CI receipts do not accept the newer
C2-C4 source or its shared-lowering/work-ceiling changes. Native GUI and unproved
Windows/macOS path/handle scenarios remain separate.

### Compile First, Then New Rust/VM/Parity Boundaries

Use a clean checkout, adequate scratch space and **explicit run-owned**
`CARGO_TARGET_DIR`; verify it before the first Cargo call. Never clean the shared
cache. A matching same-source/platform CI receipt can replace the same command.

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
cargo test --locked -p gentle-engine -j1 cdna_interpretation --no-run
cargo test --locked -p gentle-engine -j1 cdna_interpretation -- --test-threads=1
cargo test --locked --no-default-features --lib -j1 cdna_interpretation -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 cdna_interpretation -- --test-threads=1
cargo test --locked --no-default-features --lib -j1 coding_models_selector_rejects_annotation_reference_id_mismatch -- --test-threads=1
```

The first two root checks must not substitute for compiling the separate engine
crate's test target. Then execute the original failing kernel/service tests,
isolated JS/Lua/MCP cases and both process cases, requiring nonzero counts and
actual route entry. Check distinct namespace preservation and fail-closed wrong
selectors/shared annotation/gene/reference bindings, not renamed-equal fixtures.

Canonical generator (zero diff), registry/freshness, glossary, work-ceiling and
checkout/maturity are already passed at `2f1a3062`. This identity/test repair
does not change their registry/projection inputs; do not repeat those campaigns
as missing evidence. Normal newer-source CI remains distinct; recheck affected
projections if their inputs change. Require actual JS and Lua VM tests,
MCP denied/approved calls and project-free command policy; version smokes and
filtered-out tests do not prove those routes.

Verify plus/minus coding/UTR/noncoding/not-present coordinates, missing/inconsistent
accessions and incomplete CDS, split codons and all four TCC/TAC/TCA/TAA patterns.
Check below-threshold evidence, incomplete groups, sample isolation, inherited
counts/statuses, disjoint local blocks and explicit negative nonclaims. HGVS
follows the transcript reference, including cross-intron delins and the same-codon
one-base-gap exception; stop-loss extension counts the former stop as position 1.
Multi-codon protein composites stay withheld. Check strict frozen-label assembly,
alleles/source/UTC/license, tampering, cancellation and cumulative work/input/output.

C3 must match the original normalized evidence-request hash. No QNAME join,
cross-hypothesis allele/junction pairing or transitive splice support. Distinguish
compatible, competitor-excluding and selected-model-unique units. Retained-intron,
off-locus and foreign-contig geometry must not disappear under variant filtering.
For comparison byte preservation, additionally compare **old `4966a3a5` and new
source-matched binaries on the same prepared synthetic files/paths**. The new
same-process rerun regression alone is not a historical-byte proof. No old
private-data campaign or archive regeneration is needed.

### Real Processes And Remaining Host Boundaries

Build once per required feature profile, record binary SHA-256 and run:

```sh
GENTLE_CDNA_INTERPRETATION_CLI=/absolute/path/gentle_cli python3 -m unittest scripts.test_cdna_allele_interpretation -v
```

Both process cases must execute: real saved-evidence producer, C1, comparison,
C2/C3/C4, direct/shell JSON/TSV/summary, LF/CRLF raw bindings, quoted/apostrophe/
Unicode paths, deterministic repeats, refusal and unchanged files/project. Check
pure saved-JSON projection regeneration and confirmation before I/O. C1's old
CLI/VM campaign already passed at `4966a3a5`; do not rerun it as missing evidence.
These new process/VM tests exercise the affected newer source and reuse C1 within
the C2-C4 flow, rather than repeating the completed historical campaign.
The positive process case also invokes the actual `GentleClient.shell` wrapper
with an explicit same-source binary (no Cargo fallback), comparing JSON and both
text formats to the direct route. The separate mocked subprocess check proves
argument/format delegation only, not runtime or scientific acceptance. Owner
confirmed that Glen executes the remaining gates after the owner's push;
Codex does not push or substitute local Rust execution.

Native Windows requires actual path/handle checks; GUI Shell needs an available
native host, not inferred headless acceptance. No new GUI form is requested.
Real TP53 conditional HGVS receipt remains pending separately: only already
committed/authorized bound artifacts, no private source regeneration/other genes
or hard-coded observed counts. Do not let that optional data receipt imply a
biological verdict. OD-I1-I5 remain conservative defaults in the plan.

Repair-only Codex local evidence: combined companion/checkout/maturity suite
reports **120 passed, two same-source CLI skips** (122 cases), with the new
distinct-ID request regression and expanded LF/CRLF checkout coverage. Prior
Python-wrapper/C1 source checks remain historical, not runtime proof for this fix.
Documentation status, scoped formatting/syntax and whitespace pass. Synthetic
evidence only; repaired-source Rust, VM, historical-byte comparison, actual
process and native checks are **not_run**. The generator/freshness pass at
`2f1a3062` above remains source-bound; no new generator run is claimed.
Return primary diagnostics before cascade failures.

## Four Read-End Fixes At 18488691, 2026-11-12

Historical authoring checklist: all requested Linux filters, including configured
real samtools, passed in Glen's clean run at `4966a3a5` above. Do not repeat these
completed checks. Unproved native-platform scenarios remain separate.

**Original request, completed on Linux at `4966a3a5`:** test `184886914e9995f9af08a4ee83ae104e918e6121`
or record the exact merged SHA tested. This checklist is separate from C1.
The 208 checks at `d8656d46` predate all four fixes; no private cDNA data or
repeat of that campaign is needed.

Use a clean checkout and a run-owned `CARGO_TARGET_DIR` on Unix and native
Windows. Start with the locked compile gates; a same-SHA, same-platform CI
receipt for the same gate can be cited instead of rerunning it.

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
```

1. **samtools environment (`53558ec5`).** The native stub must see the
   parent-only variable in its control call, but not in GENtle's version probe,
   header read or region selection. All three must see `LC_ALL=C` and the
   private `REF_PATH` and `REF_CACHE`. Glen: also set `GENTLE_TEST_SAMTOOLS` to
   the absolute path of a real executable and record its version. Confirm the
   indexed BAM and reference-backed CRAM cases both execute on the synthetic
   fixture. If the variable is unset, the test prints `NOT RUN` and returns
   normally: record `not_run` for real samtools, regardless of Cargo's `ok`.
2. **Pileup sorting (`c5a91881`).** For 5'-UTR, 3'-UTR and read length,
   exact observations must precede right-censored lower bounds in both sort
   directions. An observed 20 bp must not tie with a lower bound of 20 bp.
   Check stable ties across pages and missing values first ascending, last descending.
3. **Selected SAM retention (`abc60048`).** Omitted or false
   `retain_selected_sam` must leave no SAM artifact. True must publish
   `artifacts/selection.sam` byte-for-byte with the receipt's size and SHA-256;
   verification must reject same-length tampering. Cancellation after staging
   but before report commit must leave no retained SAM on disk or in the
   cancelled report. Check cleanup on native Windows as well as Unix.
4. **Input budget (`18488691`).** Bound JSON must admit 8,387,584 raw bytes
   and refuse one more. Shell preflight must report the submitted and allowed
   byte counts and reject an overflowing declaration. A two-sample build must
   admit 8,386,560 charged raw bytes in total and refuse one more without
   publishing output. Discovery must distinguish the single-charge ceiling
   from the 128 MiB file cap. See [input accounting](rna_read_end_metrics.md#read-end-input-budget).

Run these headless library filters in checklist order:

```sh
cargo test --locked --no-default-features --lib -j1 indexed_locus_isolates_all_three_samtools_calls
cargo test --locked --no-default-features --lib -j1 indexed_locus_optional_real_samtools_requires_explicit_configuration -- --nocapture
cargo test --locked --no-default-features --lib -j1 read_end_pileup_value_sorts_separate_observations_from_bounds
cargo test --locked --no-default-features --lib -j1 indexed_locus_retains_selected_sam_only_on_opt_in
cargo test --locked --no-default-features --lib -j1 read_end_service_bound_json_input_byte_boundary
cargo test --locked --no-default-features --lib -j1 read_end_service_build_input_
cargo test --locked --no-default-features --lib -j1 read_end_oriented_discovery_names_bound_originals_and_preserves_typed_inspection
```

The `read_end_service_build_input_` filter must run both the shell preflight
and cumulative two-sample tests. Return the source SHA, OS, Rust/Cargo and
samtools versions, commands, nonzero test counts, failures and explicit
`not_run` cases per platform. This is a verification request, not a pass receipt.

## C1 Bound Coding Models, 2026-11-11

Historical authoring checklist: compile/core/GTF/shell/MCP/VM/real-process/parity
passed in Glen's clean Linux run at `4966a3a5` above. Only unproved native-host
scenarios and affected newer-source regressions remain;
the newer C2-C4 checklist supersedes this section's former C1-only scope limit.

**Original request, completed on Linux at `4966a3a5`:** verify `DeriveCodingTranscriptModels`
and `rna-reads coding-models-derive`, not another cDNA/exon-locus campaign.
The [C1 contract](coding_transcript_models.md) is separate from unchanged
`CompareCdnaAlleles` and normalized annotation v1. C2-C4/HGVS/consequences,
private reads and regeneration of the real fixture/analysis archives are not
authorized by this checklist. No fixed gene/barcode/provider rules.

### Compile, Grammar And Shared Routes

The owner's local compile of `d927ffa8` exposed two primary C1 errors before
runtime acceptance: E0277 (missing `EngineError` to shell `String` conversion)
and E0505 (a named input iterator retained request borrows until drop).
The narrow repair converts errors only at the existing shell boundary and
consumes the temporary iterator before moving the owned request into its report.
It adds no input clones or changes to bindings, budgets, schemas or biology.
Verify the repaired integrated source, not the original `d927ffa8` binary;
start with the compilation gates below before generator/VM/CLI execution.
New Rust regressions cover exact shell error text/state preservation and the
exact cumulative input-byte boundary with and without optional protein FASTA.
The actual-CLI tamper case now also requires the original reader diagnostic.
These Rust/process checks remain unrun by Codex; no new green compile is claimed.
Repair-only local checks: three Python/source/LF-CRLF cases pass; the two real
CLI cases skip without a supplied binary. Scoped formatting, documentation
status and whitespace checks pass; these do not replace a compiler run.

Use one clean source checkout and an explicitly run-owned target (never an
inherited shared `CARGO_TARGET_DIR`). Record the actual integrated SHA, lock
identity, Rust/Cargo, platform, features/profile, commands/counts and skips.

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
cargo run --locked --quiet --bin gentle_examples_docs -- parity-matrix-generate
cargo test --locked --no-default-features --lib -j1 coding_models -- --test-threads=1
cargo test --locked --no-default-features --lib -j1 gtf_preparation -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 coding_models -- --test-threads=1
cargo test --locked --lib -j1 glossary_ -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
```

The new matrix delta was derived from registry rules **without running its Rust
generator locally**. Confirm the canonical output/freshness, return any differing
patch, and distinguish primary compiler failures from cascades. Maturity catalogs
were regenerated with the Python projector; no evidence maturity was raised.
Require nonzero test counts, including actual JS/Lua VM and generic MCP cases,
the service observation/capacity test and typed discovery/parser regressions.
Existing GTF normalization tests are relevant because its shared grammar now
has opt-in coding capture; old annotation output must stay byte-identical.

### Synthetic Boundary And Actual CLI

Check plus/minus split codons without per-segment trimming, included/excluded
stop equivalence, explicit nuclear/table-1 refusal, malformed geometry and raw
tampering. Preserve separate mRNA/CDS completeness: `mRNA_end_NF` with complete
CDS must still translate, CDS-specific/missing-start/stop evidence must not.
Retain every biological reason, missing/inconsistent versions, optional protein
comparison, noncoding/undefined CDS, budgets, cooperative cancellation and
determinism. No fallback provider convention, repaired reference or HGVS output.

With a freshly built **same-source** CLI (record its SHA-256):

```sh
GENTLE_CODING_MODELS_CLI=/absolute/path/gentle_cli python3 -m unittest scripts.test_coding_transcript_models -v
python3 -m unittest scripts.test_tutorial_checkouts scripts.test_capability_maturity -v
```

The two actual-process cases must execute, not skip. They cover LF/CRLF raw
hashes versus equal model semantics, spaced/apostrophe paths, repeat output,
tamper refusal and absence of project/output files. Also use the direct
`gentle_cli rna-reads coding-models-derive @FILE` route. MCP's generic `op`
confirmation must remain before I/O; readonly discovery is not a bypass.
Windows path/handle and native GUI Shell observations need actual available
hosts; report `not_run` rather than inheriting Linux/macOS results. There is no
new dedicated C1 inspector or biological acceptance gate implied.

Local evidence is Python/source/checkout and formatting only. Rust, VM/MCP and
actual C1 CLI execution are **not_run locally**. Do not repeat the completed
208-check d8656d46 return, real Q20/Q30 analyses, locus exports, reader/BLAST,
barcode or packaging campaigns merely for this new read-only contract.
Initial C1 combined command: `python3 -m unittest scripts.test_coding_transcript_models
scripts.test_tutorial_checkouts scripts.test_capability_maturity -q` reports
119 cases: **117 passed, two same-source CLI skips**. Documentation-status and
whitespace checks pass. Session-close: four OK, two intentional warnings
(dirty/untracked work and plan-label/manual-fidelity reminder), no failures.

## Glen Integrated cDNA And Exon-Locus Return At d8656d46, 2026-11-11

Owner-forwarded result for **`d8656d46adb065d212231eef5f4ded071538a8d3`**,
not just `df47b05e`. It contains fixture/provenance commit `73fd19b1` and the
cDNA comparison/transitive-phasing commit `df47b05e`; the intervening final
commit changes only the release-note ledger. Glen reports technical satisfaction
and no source fix. This closes the reported focused checks at this source,
not native GUI, scientific acceptance or the selected-release ledger. Do not
repeat the completed campaigns or execute the historical authoring lists below
as new requests.

### Focused Return

Glen reports **208 focused test executions passed, zero skips**, plus
`cargo check -q --locked` passing:

| Reported scope | Passed |
| --- | ---: |
| `gentle-engine` library, including cDNA | 91 |
| Headless comparison/partitioning | 20 |
| SAM / frozen S2A / read-end service preservation | 12 / 4 / 12 |
| New fixture provenance/checkout boundaries, not the full 97-test suite | 10 |
| Actual cDNA CLI | 2 |
| Pysam fixture preparer | 10 |
| MCP parity | 1 |
| Exon-locus | 33 |
| Renderer | 7 |
| Actual LF/CRLF CLI rehearsal | 3 |
| Actual JS / Lua / MCP runtimes | 1 / 1 / 1 |

Partitioning retains original IDs, orientation, bases and qualities, combines
units before SNP/phasing and refuses duplicates. Supported transitive local
links do not require a full-marker read; weak/ambiguous/disconnected components
stay unresolved and K is not forced to two. The actual anchor regression keeps
query `[65,77)`, reference `[165,177)`, anchor **171** before window start 200,
complete JSON/TSV and no invalid glyph. All five export formats, repeat
determinism, pixel/no-clobber guards and 18/18 fixture hashes are reported green.

Glen retains sealed evidence at `/home/clawbio/work/evidence/gentle-d8656d46/`.
Only abbreviated CLI (`9ddf98b0...149e6`) and lock (`d880f70b...569fe`) hashes
were forwarded here. The complete manifests, raw commands, environment/profile
and full hashes have not been independently read/verified by Codex. Keep this
attribution; do not invent full hashes or transfer older binary identities.
No private archives were regenerated and the eight real-data analyses were
not rerun. Native GUI/display acceptance is `not_run`.

### Independent CI Observation

Codex read GitHub's run/job summaries on 2026-11-11 for the same exact SHA:
[run 42945255867](https://github.com/smoe/gentle_rs/actions/runs/42945255867),
attempt 1, **completed/success**, final update `2026-11-11T06:43:07Z`.
[Windows](https://github.com/smoe/gentle_rs/actions/runs/42945255867/job/130276401043)
and [macOS](https://github.com/smoe/gentle_rs/actions/runs/42945255867/job/130276447264)
pass, including their `Full test suite` workspace steps, adapter/all-features
checks and normal tutorial/parity gates. Linux/headless, release safety and the
selected-platform aggregate pass. **Full Linux is skipped**, as are optional
08.04 GUI, scroll CPU and authenticated-macOS audit jobs. These are independently
read status/step summaries, not downloaded full logs, package receipts or a
native GUI verdict. Glen's earlier `in_progress` observation remains historical;
no rerun/dispatch was needed or requested.
Passing configured Windows CI steps does not certify native GUI or additional
CLI/pysam/path/handle scenarios that were not selected by that workflow.

### Separate Routine Scan And Scientific Limits

The separate routine scan covers 17 commits after `7fa6b15a` through this same
SHA and reports no concrete regression: Rust 26 engine-cDNA, 18 ingress, six
locus-service and 26 locus-core passes; Python 11 passes and four expected
CLI/pysam skips. Do not add these overlapping scopes to the focused 208 or turn
the routine's skips into passes. Its first default-feature attempt was stopped
after discovering the shared `/home/clawbio/GENtle/target`; that attempt has no
test result and the shared cache was not cleaned. The isolated target/worktree
were removed (5,873 files / 7.2 GiB). No patch was recommended.

The [retained TP53/TP73 analysis](tp53_tp73_cdna_fixture_analysis_20261110.md)
already records the numerical limits. SK-MEL-28's local TP73 TAT/TCT classes
at Q30 have 21/8 assigned units but only 1/2 complete-vector units; transitive
support is not a complete chromosomal allele or independent molecule count.
WM793/TP73 stays **unassessed**, not a negative second-allele result. Shared
1:3722110 T evidence (73/73 versus 2/2 at Q20), sparse 1:3733753 callability and
three WM793 candidates outside all 20 selected mature models do not establish
cell-line differences. Absent switch/K>2 findings establish neither biological
absence nor normal copy number. More informative independent WM793 evidence, complete
mapping views/supplementaries and orthogonal DNA/library evidence remain
scientific follow-ups, not an authorized new data run or a product-code fix.

Remaining work is native GUI when available, full-Linux and candidate-bound
release/package/container/readiness/owner gates, plus the explicitly separate
scientific questions. No source changes, local Rust, private-data analysis,
push, dispatch, tag or publication were performed by Codex for this receipt.

## Real-cDNA Fixture Provenance And Checkout Repair, 2026-11-10

Glen's bot commits `a7d16efddae245955dccd4a8c174faeb853665a1` and
`42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649` are integrated from base
`249173ee`, with all 18 data hashes preserved. The corpus remains an explicitly
documented, owner-authorized asset, not an implicit private-data run or evidence
of phasing/biology. The archived Codex verification prompt does not grant local
Rust/build, provider, data-run or publication authority.

The new read-only verifier authenticates frozen binding metadata, then checks
actual size/full streamed SHA-256 of both CRAMs/receipts, reference, full GTF and
sample sheet before output creation. Header plus per-gene extracts must match
the frozen combined-GTF hash too. Regeneration writes computed verified facts;
relocated paths do not authorize new source bytes. Inputs must stay immutable;
no source lock or independently authenticated auxiliary index is claimed.

Fast regression scope: `scripts.test_tutorial_checkouts.CdnaFixtureSourceTests`
and `TutorialCheckoutTests.test_cdna_fixture_manifest_preserves_text_and_binary_bytes_in_both_checkouts`.
Only synthetic source bytes/tool stubs and retained-archive hash checks are used,
not CRAM analysis or private archival reads. The Bash-wrapper case is POSIX-only;
verifier and LF/CRLF Git checks are portable, not a claimed native Windows pass.
Local Darwin arm64 / Python 3.15.2 result:
`python3 -m unittest scripts.test_tutorial_checkouts -q` **97/97**,
zero skips (the ten new provenance/checkout tests are included, not additional).
`bash -n` passes; `shasum -a 256 -c SHA256SUMS` verifies **18/18** data files.
The data subtree is byte-identical to bot commit `42e819c9`; documentation and
whitespace checks pass. No Rust build/check/test or actual extractor ran.
Full regeneration against archival inputs was not run and needs separate owner
authorization if desired; do not repeat it or the completed earlier campaigns
as part of the remaining exon-locus/CLI/VM checks below. Return per-item results
and retain the actual source/lock/binary/input identities for those open checks.

## External CF01 SiLA Review Boundary, 2026-11-10

Context only, not a new GENtle test campaign. The owner supplied Claude's
read-only review of an external CF01 worktree at abbreviated revision
`7378a04`. That object and the named `cf01_sila2_runtime.py`,
`cf01_sila2_store.py`, gallery scripts and `docs/sila2-network.md` are not
present in this GENtle checkout. No external repository was fetched or changed;
the reported logs, hashes and historical-byte comparisons were not independently
verified here. Its LockController v1 result must not replace GENtle's explicitly
scoped v2-draft architecture review or count as GENtle/physical acceptance.

Relevant future evidence/journal safeguards are now recorded in the SiLA
architecture section. The external owner still needs revision-bound sealed
receipts, an actual old-profile runtime/reconciliation negative control with
pending-intent preservation, documented old-store inspection/recovery and a
semantic review of the changed batch source. Gallery near-match/negative checks
remain external test follow-ups, not GENtle regressions. A local laboratory demo
does not test the SiLA transport. Do not rerun those external campaigns or claim
their tests passed on GENtle; no network/device call or implementation is added.

## Claude Exon-Locus Return At b95aae07, 2026-11-10

Owner-forwarded Claude review/test report, not a new Codex execution or Glen
return. Tested base: `b95aae0774c105efcda808fa07e9e49f0b8bde8c`, macOS arm64,
headless lib-test build, committed Cargo.lock; the owner's then-uncommitted
manifest/lock updates were excluded. The scratch checkout/target was removed
(reported 6 GB freed). All nine repairs/defaults were source-reviewed positively.
No raw command streams, binary/lock hashes, toolchain version or manifest receipt
accompany this message; retain its attribution rather than inventing them.

Reported renderer `read_end_`: **7/7 pass**. Original `read_end_locus_`:
**31/32**, failing only the off-window internal-run test's expected anchor 165.
The unchanged production rule attaches aligned interval `[165,177)` at its
midpoint **171**, still below window start 200. A scratch overlay changing the
two expectations to 171 produced **32/32**. This is an overlay result, not a
pristine `b95aae07` pass, and an initially masked pipeline exit is not evidence.
Use the actual Cargo exit or explicit pipeline failure propagation in future
runs. The supplied two-line correction is retained on main; added
query/reference-interval and outside-window assertions have not run locally.
No production anchor/scoring/limit change is needed.

The scale fixture reports **5,543,670 work units**, 30 models/12 exons/20 reads
with deletion-split blocks and no originals. All models share one exon chain.
About 277,184 work/read and roughly 360 reads under 100M are arithmetic
extrapolations, not a measured refusal threshold, profile or throughput claim.
Opposite-strand comparison cost is suspected from source, not CPU-sampled.
OD-L4 remains open for distinct-isoform/real capacity and any separately
authorized optimization; do not raise caps or alter charging now.

The later owner dependency commit `e3b93bf3` changes Cargo.toml/Cargo.lock, not
these runtime sources; this return does not accept that lock or newer source.
Remaining scope: current-source locked/default checks, actual JS/Lua/MCP and
CLI LF/CRLF exports below, native GUI when available and normal Windows CI.
Keep these separate from the reported scoped passes; do not repeat the completed
7fa6b15a/reader/BLAST/barcode campaigns. `Membership` keeps its documented
kebab-case wire values; a cosmetic snake_case migration is a separate choice,
not part of correcting a test expectation. Codex follow-up checks: Python
read-end discovery **53 passed / 1 actual-CLI skip**, scoped rustfmt parsing,
documentation status and whitespace pass. No local Rust build/check/test, push
or dispatch; the added Rust assertions remain unexecuted here.

## Vector-PDF Writer Dependency Alignment, 2026-11-10

The owner's dependency update selected direct pdf-writer 0.15 while svg2pdf
0.13 still returns 0.12 `Chunk`/`Ref` values. The reported E0308 errors are a
typed dependency mismatch, not an `Extend` dispatch problem. Restore the direct
requirement to `0.12.1` (compatible patches allowed) and resolve the shared
writer. Preserve the other owner dependency changes and existing PDF semantics.
The offline `scripts.test_release_candidate.DependencyPinTests` now guards this
shared lock identity with LF/CRLF, split-version and missing/ambiguous controls;
it does not prohibit a future coordinated converter/writer upgrade.

Local no-build evidence: offline release-policy suite **57/57**, including the
eight dependency guards; `cargo tree --locked --offline -i pdf-writer` confirms
one 0.12.1 writer shared by GENtle and svg2pdf. Documentation and whitespace
checks pass. These are resolution/source checks, not a Rust compile or PDF run.

After the owner's push, run or reuse matching-source results for:

```sh
python3 -m unittest scripts.test_release_candidate.DependencyPinTests -v
cargo check --locked --tests -j1
cargo test --locked --lib -j1 svg_vector_pdf::tests -- --test-threads=1
```

Require nonzero counts, including the multipage/order/page-size and invalid
later-page output-preservation regressions. Return the tested SHA, lock hash,
toolchain/platform and raw commands/results. No local Rust execution, native,
package or release pass is claimed. Do not repeat the completed 7fa6b15a,
reader/BLAST or barcode campaigns; the new exon-locus checklist remains separate.

## Glen Routine Scan At 7fa6b15a, Received 2026-11-10

Separate routine scan, not the earlier completed focused campaign below.
Forwarded range: `a549626664fe88aeaa1dd5f634f51cc09d4b876d` through
`7fa6b15adc28d0adfe3e1edce0190a2621b54705`, two commits, clean detached
worktree; scan boundary recorded as `2026-11-10T04:40:35Z`.
Diff/whitespace and scoped formatting pass; renderer **6/6** passes. No
concrete regression or source patch is recommended. The broader focused Rust
attempt stops during compilation with **No space left on device**: incomplete
execution/environment limitation, not a passing test gate or a proven code bug.
No toolchain/profile/raw-stream/archive receipt accompanies this message.

The Python return is worded as "8 passed, 2 skipped". Source at the tested SHA
has eight test methods total: six fixture/source checks and two optional actual
CLI cases. Without its raw unittest summary, retain the reported wording as
unreconciled rather than claiming ten tests or eight successes plus two skips.
Both actual-process cases are skipped in this routine; none is newly accepted
by it. A retained summary is sufficient for count reconciliation; do not rerun
the old campaign merely to repair wording.

The return reports removal of its detached worktree and run-owned Cargo cache
(4519 files / 1.7 GiB). Its first renderer attempt had resolved to the shared
`/home/clawbio/GENtle/target`; that cache was correctly left untouched. Before
the next authorized build, verify the effective run-owned target and available
space on the actual target, Cargo-home and temporary-filesystem volumes. Abort
on a shared-target mismatch or insufficient space; do not clean shared caches
or an active build to make room. No change to Cargo configuration is requested.

This incomplete routine does **not** revoke or replace the separately attributed
[completed 7fa6b15a focused return](#glen-oriented-pileup-return-2026-11-09),
which reports 47/47 root/VM/MCP and 8/8 actual rehearsal with zero skips.
Do not sum the two runs, relabel skipped cases as passes, or infer new-source
acceptance. The exon-locus extension/Claude repairs are committed locally at
`244ad13952ac2772177028a6d53eb899932427ef`; this older scan does not contain or
verify them. After the owner's push, follow only the active repair scope below.
Native, Windows, release and scientific acceptance remain separately bounded.

## Exon-Locus Review Repairs, 2026-11-10

Authoring scope, with the later scoped macOS return recorded above. After the
owner's next push, verify only its remaining or affected boundaries. This supersedes the unexecuted
2026-11-09 exon-locus checklist below, not Glen's already completed transcript-
bound `7fa6b15a`, reader/BLAST or barcode campaigns. The owner supplied Claude's
read-only source review; Codex traced its structures/counterexamples and authored
the bounded repairs. It is **not** a reproduced Rust failure or passing Rust/VM/
native result. No local Claude invocation, Rust build/check/test, push or dispatch.
The unreleased locus v1 schemas and RULE_VERSION remain unchanged.

Local checks on Darwin arm64 / Python 3.15.0, `main@7fa6b15a` plus the working
tree (not a committed or runtime-certified candidate):

```text
python3 -m unittest discover -s tests -p 'test_*rna_read_end*.py' -v
Ran 54 tests; OK (53 passed, 1 actual-GENtle-CLI skip)
python3 -m unittest scripts.test_tutorial_checkouts -q
Ran 87 tests; OK
rustfmt --edition 2024 --check src/observed_transcript_structures/read_ends/locus.rs src/rna_read_ends/locus.rs src/rna_read_ends/locus_tests.rs
exit 0 (includes child tests/presentation; parsing/formatting only)
check_documentation_status.py: OK
manual_docs.py: OK
git diff --check: exit 0
```

Session-close: **4 OK, 2 warnings, 0 failures**: intentional dirty worktree and
manual plan-fidelity reminder. Scope was manually reconciled to repairs 1-9;
no counterexample remained untraced. The two withdrawn claims were not changed.
Work/capacity remains unmeasured locally; no Rust, actual VM/CLI, native GUI or
Windows result can be inferred from these Python/source checks. The explicit
actual-CLI skip remains open, not a pass. Cargo/dependencies/lockfile unchanged.

Repair files: `read_ends/locus.rs`, `locus/presentation.rs`, existing locus tests
and new `locus/review_tests.rs`; `rna_read_ends/locus.rs` and `locus_tests.rs`;
synthetic preparer/`tests/test_rna_read_end_locus.py`, discovery wording and the
linked decision/measurement/plan/interface/provenance docs. Earlier uncommitted
GUI/renderer/adapter changes remain part of the overall integration, not a newly
executed acceptance result. All unrelated `paper`, data, `output/`, untracked
prompts and recovery stashes remain untouched.

At authoring, thirteen new core test functions and three new service test
functions were **unrun**, alongside strengthened existing controls. The later
reported headless filter/overlay and renderer results are recorded above; the
following list still identifies the remaining/affected execution paths. Exact
main filters below select them with the earlier relevant locus tests. Check
test counts are nonzero; do not treat compile failure, a zero-match or a skipped
actual process/VM case as acceptance. Use one isolated checkout and explicit
run-owned target, not a shared target/cache. Keep the actual tested source and
binary hashes, platform/toolchain/features/profile and raw commands/results.

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
cargo test --locked --lib -j1 read_end_locus_ -- --test-threads=1 --nocapture
cargo test --locked -p gentle-render read_end_ -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 read_end_js_runtime -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 read_end_lua_runtime -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 read_end_mcp_discovery -- --test-threads=1
GENTLE_READ_END_BIN_DIR=/absolute/source-matched/bin python3 -m unittest tests.test_rna_read_end_locus -v
```

Specific targets under the locus filter: off-window softclip/internal-run plus
JSON-only service export; one/three-base overhang and intronic containment;
mirrored identical chains with/without a third model; indel-split full versus
partial junction chains; unavailable causes; query-disjoint collinear/chimeric
segments; positive novel compatibility, disagreeing hypotheses and CDS metrics;
run-variant unions, named budgets, two-sample/track filters, raw-only membership,
declared-absent annotation/CDS versus wrong hash, missing originals and exact
cancellation errors. Existing junction-backed support/RC controls must stay green.
Do not act on the withdrawn overlapping-query chimera or unknown-axis claims.

Retain actual `Report.work`/`Page.work`. The scale target
`read_end_locus_scale_30_models_12_exons_20_deletion_split_reads_records_work`
prints its work under `--nocapture`; it has 30 models, 12 exons, 10 minus-strand
models, 20 reads/~40 deletion-split blocks, no originals and a 100M-work ceiling.
Record this number and work for any independently authorized real locus run;
do not start a new private-data campaign merely to measure it. Work is shared
across samples, original-base scanning still charges per base, S2 still admits
at most 4096 units per analysis. No cap/charge/sweep tuning or capacity pass is
claimed by source checks. The later reported scale measurement above does not
establish real capacity; OD-L4 therefore remains pending.

The actual CLI test now executes both LF and CRLF input graphs, explicit model
projections, all five formats, deterministic output/no-clobber and a narrow-
window JSON-only export with unchanged off-axis run data. Fast Python fixture
generation/hash checks are **not** GENtle behavior/line-ending acceptance.
Native GUI (track position, scrolling, off-axis labels, counters, stale/cancel
receipts), native Windows paths/handles and real-data/scientific acceptance stay
separate. Use matching normal CI for Windows, without a new unauthorized dispatch.
Do not regenerate parity/catalogs speculatively: no new command/registry row;
retain ordinary freshness gates on the integrated source.

OD-L1..OD-L4 remain owner-pending beyond the selected defaults in the
[repair ledger](rna_read_end_integration_plan.md#exon-locus-review-repair-ledger-2026-11-10).
Return per-item pass/fail/not_run/blocked, source/binary/lock/input hashes,
commands/counts/skips, receipts/artifact paths and residual limits. No package,
release, provider/model call, publication or biological approval is requested.

## Exon-Locus And Scored Alternatives, 2026-11-09

Historical source-authoring checklist, superseded by the 2026-11-10 repair scope.

Active **new-source** scope after the owner's push; not an interim freeze or a
request to repeat completed 7fa6b15a/reader/BLAST/barcode campaigns. The canonical
prompt is bound below to `d72152e0...d79df6d`. Implementation is source-authored,
not Rust/VM/native accepted. Claude's subsequent review is
[prepared](exon_locus_claude_review_prompt.md), not performed.

Local authoring checks: normal Python read-end discovery **53 passed / 1 actual
CLI skip**, older rehearsal **6 passed / 2 actual CLI skips**, and tutorial/manual
LF/CRLF checkout suite **87/87**. These scopes are disjoint; reruns are not added
to the totals. The first package-style core/boundary invocation lacked the tests
import path; rerunning with the supported discovery/PYTHONPATH succeeded.
Documentation/manual freshness, scoped rustfmt parsing and whitespace pass.
Session-close: **4 OK, 2 warnings, 0 failures** (intentional dirty worktree and
manual plan-fidelity review). Unrelated `paper`, private data, `output/`, prompt
files and stashes are untouched. Dependencies/lockfile are unchanged.

Changed boundaries: pure locus records and scored hypotheses, shared stable
read keys/exon map, nested `RnaReadEnds::Locus` loader/actions, existing renderer
annotation tracks/sublanes, thin native controls and strengthened actual
JS/Lua/MCP fixtures. The gene-filtered v1/v2 path is preserved, not repurposed.
No new command path/registry row; do not hand-edit canonical projections.

Use one clean isolated checkout, explicit run-owned target and locked single-job
builds; record the actual tested integrated SHA, toolchain/features/platform.
Do not inherit/clean the shared target. Suggested focused commands (external,
not run by this authoring thread):

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
cargo test --locked --lib -j1 read_end_locus -- --test-threads=1
cargo test --locked -p gentle-render read_end_ -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 read_end_js_runtime -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 read_end_lua_runtime -- --test-threads=1
cargo test --locked --lib -j1 read_end_mcp_discovery -- --test-threads=1
GENTLE_READ_END_BIN_DIR=/absolute/source-matched/bin python3 -m unittest tests.test_rna_read_end_locus -v
```

Require nonzero counts, the actual CLI case (no skip-as-pass), the new GUI helper
test, actual VM execution and matching normal CI. Reuse existing registry
freshness evidence unless affected; normal CI still retains its mandatory gate.
The shared exporter/layout changed: retain PNG/PDF dimensions/pixel limits and
closed staged/no-clobber behavior, not only Rust source assertions.

Verify prompt tests 11-21: shared exon stays uninformative on both BAM directions;
sense/antisense-specific junctions outrank endpoints/A-T; RC/mirror invariance;
identical models remain ordered alternatives; unknown motifs are not fabricated
as canonical; supplementary flags alone are not chimeras. Model switch keeps
genomic blocks, read IDs and declared/locus counts unchanged. JSON/TSV export all
candidate/path scores, evidence/conflicts, rule version and annotation bindings.
Inspect missing annotation without losing genomic/original-base evidence; changed
binding recomputes scores, present wrong hash refuses; no stale-score reuse.
Check duplicate-source refusal, hypothesis/candidate/run/geometry budgets,
request-bound cancellation and input protection. Partial reference overlap stays
distinct from the full-coverage boolean; original A/T letters never rename.

Native GUI when available: genome and explicit transcript selection, both-strand
tracks, unique reads with hypothesis sublanes, off-window links/clips, reasons,
class/paging counters, long labels/scrolling, stale-view binding and cancellation
with terminal receipt. Native GUI is not inferred from headless helper tests.
Windows paths/handles/racing output and normal platform CI remain external;
Glen is Linux, so use the next owner-triggered matching Windows CI rather than
requiring a nonexistent Windows device or dispatching extra jobs.

Return per-item pass/fail/not_run/blocked, exact source/binary/lock/input hashes,
commands/counts/skips, profile/features/platform, artifacts and unresolved limits.
No private-data run, live model/provider call, release/package/dispatch or
scientific acceptance is requested. Candidate selection/publication is owner-gated.

Local repair-loop tooling now follows the
[recorded nextest workflow](testing.md#local-failure-collection-and-recorded-reruns).
The macOS development host has nextest 0.9.150 and recording enabled; this does
not install it on Glen's host or establish a Rust pass. Earlier plain-Cargo
results cannot be resumed as recordings. Keep failed-subset reruns diagnostic,
finish with a complete intended scope plus doctests, and retain normal CI gates.

## Read-End Integration Plan, 2026-11-08

**Current status:** Glen reports the [scoped oriented-pileup pass at 7fa6b15a](#glen-oriented-pileup-return-2026-11-09).
Do not repeat its completed Rust/VM/MCP/CLI campaign without affected source changes.
Native GUI and the still-running same-source CI remain separate.
The [pixel-guard Linux return](#glen-pixel-guard-linux-return-and-release-inventory-2026-11-09)
closes its earlier focused retry at `01d42b71`, not this later renderer/service.
Release-candidate selection/execution still requires the owner.
The earlier compile-failure/retry instructions below are historical; the
generator and listed downstream scopes have since executed at `7c5214ec`.

The [reconciled plan](rna_read_end_integration_plan.md) supersedes only the earlier
five-step organization. OD2 now permits only the separately authorized explicit
original-base input; other OD defaults remain. Do not execute
future-step commands merely because a plan authorizes development: each step's
implementation receipt/filter and runtime status will be recorded separately.
The prior core at `1df2db13` is still Rust-unverified, not accepted by the earlier
containment pass. No repeated BLAST/reader or private-data campaign is requested.

## Glen Oriented-Pileup Return, 2026-11-09

Glen's forwarded return binds exactly
`7fa6b15adc28d0adfe3e1edce0190a2621b54705` in an isolated checkout with an
owned target. Both locked checks pass (default tests and headless library),
renderer `read_end_` **6/6**, root library plus actual JS/Lua VMs and MCP
`read_end_` **47/47**, and the source-matched CLI build pass. Python rehearsal
is **8/8 with zero skips**, including both actual CLI process cases. The VM/MCP
checks are part of the 47, not additional test totals. No product defect reported.

The return confirms both transcript/read orientations, mRNA-relative groups and
Alternatives, original 12A/18T and query coordinates, literal 60/10 policy,
softclips and partial/complete alignment, unchanged denominators, all five
formats, repeat exports, complete labels, actual-layout pixel limits,
no-clobber and absence of a project file. Keep this scoped to the tested
transcript-bound contract, not physical tails or scientific interpretation.

Reported CLI SHA-256:
`1f6416b56f2f9ac631608f9cb3d6d571207de7bec56a51ef9890a4652b460258`.
The CLI self-reports the tested revision. Reported unchanged lockfile SHA-256:
`718c89649f21f2c8c437365fd33b256a40dce41eae7097ef1b4e7bc9a54e1314`;
Codex independently checked that hash against the local lockfile. This message
does not supply an audit archive/raw logs, toolchain/profile receipt or artifact
location; preserve those when supplied, without inventing an independent
manifest/binary verification or reopening completed tests merely for bookkeeping.

Codex read the existing same-source [CI run 42648565446](https://github.com/smoe/gentle_rs/actions/runs/42648565446)
after the forwarded status: Windows, headless and release safety now pass;
macOS remains in progress, full Linux and optional GUI/scroll/audit jobs are
skipped. This is a snapshot, not a final overall-CI verdict, native GUI or
Windows real-tool/package acceptance. No dispatch/retry was performed.

Remaining current-contract checks are native GUI (`not_run`) and completion of
the matching platform CI; release packages/container/readiness/owner approval
remain separately gated. Do not repeat completed reader/BLAST, containment,
private-barcode, parity or the above scoped execution campaign.

Glen explicitly identifies a **separate proposed extension**: exon-locus mode,
shared genomic track and scored gene/antisense `Assignments`. It is not in
`7fa6b15a`, and tests 11-21 are not satisfied by this pass. The expanded prompt
has now been [received via HTTPS](https://gentle.functional.domains/openclaw/canvas/codex-exon-locus-prompt-20261109-d72152e0/report.md).
Codex fetched it through the host's normal HTTPS route and verified its exact
SHA-256 against the owner's value and the `report.md` entry in the
[manifest](https://gentle.functional.domains/openclaw/canvas/codex-exon-locus-prompt-20261109-d72152e0/bundle-manifest.json):
`d72152e0197392e2f820312e5a338e860e56cc7a5737f94ef057fa610d79df6d`.
Only that Markdown entry was checked, not every bundle asset. No allowlist
bypass or credentials were used during retrieval. Implementation was subsequently
authorized and source-authored in the active scope above; no local Rust ran.
The [minimal draft plan](rna_read_end_integration_plan.md#exon-locus-extension-draft-2026-11-09)
separates this work from Step 7. TP73/TP73-AS may be examples, never hard-wired
gene choices. Do not relabel the verified transcript-bound feature as
incomplete merely because this later product scope was proposed.

## Oriented Original-Read Pileup, 2026-11-09

Historical implementation/verification checklist, fulfilled by Glen's scoped
`7fa6b15a` return above; native and platform boundaries remain explicit.

Verify the owner's next pushed integrated source containing this extension;
do not freeze development or substitute historical `01d42b71` acceptance.
No new operation/discovery rows exist, so do not hand-edit/regenerate the parity
matrix speculatively or reopen completed registry campaigns.

Changed boundaries: new pure `read_ends::oriented` module; optional authenticated
`Sample.original_bases` and `Build.terminal_runs`; v2 pileup/conservative v1
reading; shared renderer, thin GUI, strengthened JS/Lua/MCP fixtures; synthetic
rehearsal stages/tests and the linked request/manual/decision/plan documentation.
Original strings are explicitly supplied in sequencing 5'-to-3' orientation,
not reconstructed from SAM or implicitly joined by name. Binding checks do not
independently prove physical sequencing/trimming provenance. No actual molecule
identity, complete transcript, tail length, polyadenylation site or causality.
Local checks: nine targeted Python fixture/source tests passed, two real-CLI
cases explicitly skipped without a source-matched binary; documentation/manual
checks, scoped rustfmt parsing and whitespace passed. Session-close: four OK,
two expected warnings (intentional dirty worktree and manual plan-fidelity
review), zero failures. `paper`, private/untracked data/output and prompt files
remain untouched. Implementation is committed on main at the owner's request;
no Rust compilation/execution, push or model call is claimed.

Use one isolated checkout and explicit run-owned target, one Cargo job, locked
dependencies; retain the actually tested SHA/toolchain/features/platform:

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
cargo test --locked -p gentle-render -j1 read_end_ -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 read_end_ -- --test-threads=1
GENTLE_READ_END_BIN_DIR=/absolute/source-matched/bin python3 -m unittest discover -s tests -p test_read_end_rehearsal.py -v
```

Build/provide that source-matched CLI explicitly on the authorized external
runner. Both real-process cases must execute with nonzero counts, not skip.
The new oriented case covers plus/minus models, both read orientations, original
12A/18T versus independently reverse-complemented synthetic SAM, softclips,
JSON/TSV/SVG/PNG/PDF, repeat-export determinism, no-clobber and no project file.
The legacy process case remains a directly affected schema-compatibility control.
JS/Lua tests execute real VMs; MCP checks its additive confirmation boundary.

Require the 18 new Rust regressions plus affected existing controls: discovery
names the bound original input and defaults without bypassing typed inspection; terminal
60/10 defaults and internal/nonterminal runs; overlapping end windows; missing
versus bound-zero composition; original query coordinates unchanged by strand;
partial/complete alignment and hard-clipped original-axis refusal; distinct
sequence/mapping alternatives; stable unit keys and nonadditive A/T counters;
unchanged coverage/filter denominators; legacy hashes; LF/CRLF input hashes;
marker/input collisions, tampering, no-clobber, dense wrapping, both shared
lane/full-page heights, raster-budget refusal and actual converter dimensions.

Native review only with an available GUI: biological 5'-to-3' on both genomic
strands, same above/opposite below the divider, conflicts in Alternatives,
magenta A/turquoise T and hollow unaligned symbols, full labels, scrolling,
responsive cancellation and unchanged project/approval behavior. If unavailable,
report `not_run`; headless SVG is not native GUI acceptance. Normal Windows CI
must cover the affected input paths/closed handles/no-clobber; do not ask the
Linux-only auditor to invent Windows evidence or a GUI environment.

Return per-check pass/fail/not_run, exact source/binary/input hashes, commands,
counts/skips, features/profile/platform and retained logs/artifact locations.
No completed BLAST/reader, containment, private-barcode, biological, package or
release campaigns are requested. No local Rust, real CLI, VM or native acceptance
is claimed by Codex; Python fixture/source checks and rustfmt are not substitutes.

## Glen Pixel-Guard Linux Return And Release Inventory, 2026-11-09

Tested source: `01d42b71987fdc909d94e7900066f273d49c485d`; final fetch identical.
Glen used Debian Linux x86_64, Rust/Cargo 1.97.1, default dev/test features/profile,
an isolated checkout/target, one job and incremental disabled. Lockfile SHA-256:
`718c89649f21f2c8c437365fd33b256a40dce41eae7097ef1b4e7bc9a54e1314`.
Archive SHA-256: `74b6340d655364ffe8788055c580ae56d5288da25e7c12f7e492b2efa41eb043`;
all 17 manifest entries independently verify. The retained handoff is an exact
prefix, not a full copy; the retained release ledger matches the tested tree.

Raw logs confirm renderer **4/4**, budget **2/2**, actual small PNG/PDF conversion
**1/1**, and `cargo check --locked --tests -j1`, all exit 0 with no errors/skips.
The seven tests cover exact/over-limit/overflow refusal, underestimated wrapped
header/lane cases, shared dimension/SVG agreement and actual PNG dimensions.
The passing conversion test creates a PDF from the same SVG but does not retain
its temporary PDF for a fresh visual review. The archive's PDF reading
reference is byte-identical to the retained `1c4a4a5a` PDF, whose SHA-256 is
`3e61f5c153de70eeb7e3d0187a4f88a5790bc34bdf1cc07b6bdc13916ccf562d`.
Keep that carry-forward reference separately attributed, not a new candidate
artifact, native GUI, package or scientific verdict. No prior campaign repeated.

Codex also read the existing same-source [CI run 42614359386](https://github.com/smoe/gentle_rs/actions/runs/42614359386):
completed **success**, Windows/full suite, full Linux/full suite, headless,
release safety and aggregate selected-platform gate all pass; macOS and optional
GUI/scroll/audit jobs skipped. Windows [job 129150421407](https://github.com/smoe/gentle_rs/actions/runs/42614359386/job/129150421407)
passes its unfiltered workspace suite, covering the always-on raw-traversal
fixture. This is CI evidence at `01d42b71`, not a native GUI or Windows real-tool
acceptance. Do not ask the Linux auditor to repeat already covered Windows CI.

The point-release inventory agrees with the [.13 ledger](release_notes/release_notes_v0.1.0-internal.13.md#internal13-acceptance-ledger):
owner-selected final candidate; same-source Windows/macOS/full Linux CI; three
native packages/receipts/extracted smokes; headless container/digest/receipt and
no-network CLI/MCP smokes; remaining version/tutorial/parity/SCAD/native rack
gates; Glen readiness and explicit owner tag/publication approval. The report's
blanket missing Windows/full-Linux statement is narrowed by the existing CI
above; macOS is still not covered there. Keep all release rows Pending until
the candidate and its matching receipts are accepted under the existing owner
procedure. No candidate is selected or frozen by this evidence-recording change.

Biological UTR/Poly(A) verdicts, a raw-signal tail caller and improvements outside
the authorized .13 scope remain deferred/non-blocking, not silently accepted.
No new Rust/CLI campaign, dispatch, package/container build, model call, push,
tag, publication or scientific approval was performed or authorized by this
return. Continue independent development unless the owner starts release prep.

## Read-End Actual-Layout Pixel Guard, 2026-11-09

Historical authoring/retry checklist, completed by the `01d42b71` return above.

Glen's small focused bundle verifies `a549626664fe88aeaa1dd5f634f51cc09d4b876d`
on Debian Linux x86_64, Rust/Cargo 1.97.1 and Python 3.14.10, isolated target,
one job, incremental off. Archive SHA-256:
`7b19fb8fc622ad1ed3efc95cff64ab7ddf8a121e538f1febe171b84b43c19dd4`;
all seven manifest entries independently verify. Its raw logs pass renderer
3/3, dependency-policy 5/5 and full offline release-policy 54/54 (the five are
included). The lockfile remains `718c89649f21f2c8c437365fd33b256a40dce41eae7097ef1b4e7bc9a54e1314`.
This closes the historical Linux renderer-parser retry, not native Windows.
The separate scan's GitHub API query was incomplete (401); no CI pass is inferred.

The scan found a real resource-guard drift: export still estimated 54 pixels
per row after wrapping introduced variable header/lane heights. Ordinary
export passes at `1c4a4a5a` do not cover this oversized-page boundary.

The repair gives `gentle_render::read_end_pileup::dimensions` and SVG rendering
one shared checked layout, including wrapped headers, title/legend and all row
heights. Row/total height arithmetic refuses overflow. PNG/PDF requests now
check actual width times height with checked `u64` multiplication against the
unchanged 32,000,000-pixel ceiling before SVG serialization or conversion.
The whole mixed-format request refuses before publication. Text/vector-only
JSON/TSV/SVG requests keep their existing output-byte caps, not the pixel cap.
No biological values, report schemas, selection, fonts, paging or layout
constants change for ordinary pages; no dependency or lockfile changes.

New authored controls: one renderer overflow test; dimension/SVG equality added
to the existing 900/1400/8000-px case; two service regressions for exact-boundary/
overflow and long headers, wrapped lanes, or both. Each oversized fixture would
pass the old estimate but must refuse PNG and PDF without invoking converters;
SVG-only remains allowed. The existing small actual-converter service test now
compares the renderer's dimensions with both SVG and PNG metadata.
These Rust changes are **not compiled/run by Codex**.
Local authoring checks pass five Python source/rehearsal tests, Rust formatting/
parsing and documentation checks. Independent arithmetic confirms the three
counterexamples: old estimates 11,160,000 / 14,976,000 / 24,552,000 pixels versus
actual layout 56,784,000 / 46,000,000 / 101,152,000. This is not Rust execution.

After the next owner push, record the actual integrated SHA and a run-owned
target. Run only the affected scope, with nonzero counts:

```sh
cargo test --locked -p gentle-render -j1 read_end -- --test-threads=1
cargo test --locked --lib -j1 read_end_service_raster_budget_ -- --test-threads=1
cargo test --locked --lib -j1 read_end_service_png_pdf_derive_from_the_same_svg_with_explicit_dimensions -- --test-threads=1
cargo check --locked --tests -j1
```

Expected authored counts are four renderer tests, two budget tests and one
small PNG/PDF conversion test. Do not deliberately rasterize an over-limit
fixture; the refusal itself is the regression. Retain actual dimensions,
results/skips and logs. Reuse the completed bio/parity/47-test/VM/MCP/glossary/
supervisor/BAM/CRAM/CLI and private-data campaign evidence; do not repeat it.
Normal CI, native GUI/scrolling and native-Windows raw-path acceptance remain
separate. No scientific, release or publication approval follows.

## Compatible Bio Update Policy, 2026-11-08

Owner-authorized exception to the historical eight-pin repair: root
`workspace.dependencies.bio = "4.2.2"` permits deliberate compatible 4.x updates.
The seven other flagged pins, workspace inheritance and all other requirements
are unchanged. No Cargo update/resolution or local Rust execution occurred.
`Cargo.lock` remains byte-identical with SHA-256
`718c89649f21f2c8c437365fd33b256a40dce41eae7097ef1b4e7bc9a54e1314`.

The offline guard now verifies bio's stable bare/caret floor against a compatible
locked release and retains the other exact-pin checks. Regression controls cover
future patch/minor versions, old/next-major/prerelease/absent locked versions,
wildcards, exact pins and unbounded ranges, including LF/CRLF manifest/lock text.
Run the five `DependencyPinTests` and normal release-policy CI. There is no new
dependency version to reaccept and no reason to repeat completed RNA/BAM/CRAM/
CLI campaigns for this manifest-only policy change. Normal `--locked` CI still
must pass; a future deliberate lockfile update requires its own affected tests.
See [the decision](decisions.md#compatible-bio-updates). Keep the separately open
renderer/native/Windows boundaries separate; do not waive the existing guards.
Local verification: all five dependency-policy tests and the full 54-test
offline release-policy suite pass. Structured TOML/raw-byte comparisons confirm
only bio's requirement changed and the lockfile stayed identical. These checks
are not Cargo resolution, compiler or platform acceptance.

## Windows Specificity-Bundle Raw-Path Fixture, 2026-11-08

The existing Windows full-suite CI pass at `01d42b71` now covers this always-on
fixture; see the [current receipt](#glen-pixel-guard-linux-return-and-release-inventory-2026-11-09).
The instructions below describe the historical retry, not another Linux task.

The Windows full suite at `1c4a4a5a` failed only
`transcript_assay_panel_multi_destination_rollback_leaves_existing_evidence_alone`
in that invocation (4956 passed, 1 failed, 17 ignored). The final traversal
assertion built its input via `canonicalize().join("child").join("..").join("escape")`.
Windows verbatim-path joining can remove `..` before validation. The test now
constructs the native-separator raw string without normalization, asserts the
component remains present, and requires the exact unsafe-path error, absent
new output and unchanged foreign evidence. Production validation is unchanged;
do not remove or platform-skip the refusal assertion.

After the next owner push, run the focused test on native Windows and retain
the actual SHA/toolchain/count and raw refusal evidence, then normal CI:

```sh
cargo test --locked --lib -j1 transcript_assay_panel_multi_destination_rollback_leaves_existing_evidence_alone -- --test-threads=1
```

A Unix pass alone cannot certify the Windows fixture. No local Rust execution
by Codex. Keep this separate from the renderer-only retry and the offline bio
guard; no repeated RNA/BAM/CRAM or biological campaign is needed.

## Glen Linux Return At 1c4a4a5a And Renderer-Test Follow-Up, 2026-11-08

Historical retry: the three renderer tests now pass in Glen's `a5496266` bundle;
use the active pixel-guard scope above instead of repeating this checklist.

Tested source: `1c4a4a5aa4ac9c1515f683bb0c6fd0d91409b06d`; final fetch unchanged.
Debian Linux x86_64, Rust/Cargo 1.97.1, Python 3.14.10, locked default features
unless noted, isolated target, one job, incremental disabled. Archive SHA-256:
`660783235a9475e06b1657d1805878ebd8621577e96db0f1c1dedcd7053b5c9a`.
All 81 manifest entries independently verify. CLI SHA-256:
`b58ec855930185b53ccb55af0839f5013771c8b46308e7f8320acfae070116d1`;
lockfile SHA-256 remains `718c89649f21f2c8c437365fd33b256a40dce41eae7097ef1b4e7bc9a54e1314`.

**Passed at that source:** default test-target/headless-lib compilation;
`indexed_locus::tests` 9/9 with real `/usr/bin/samtools` 1.24, no skips;
separate synthetic CRAM/FAI/fallback 1/1 and actual BAM-to-CRAM 1/1. The exact
reruns overlap the nine-test group, not additional independent coverage.
Non-UTF-8 diagnostic preservation and strict malformed-SAM refusal pass.
The retained raw version stdout confirms byte `0xab` at offset 338 in compiler
flags (`-ffile-prefix-map=\xabBUILDPATH\xbb=.`), not alignment bytes. Raw stdout
SHA-256: `276ca13ee035861d0514bd32afb0b916c8fafd051e2c0fbf9168e32b521cd20c`.

The actual source-matched CLI/export process check passes 1/1. Fresh 1400x384
SVG has 16 text nodes within bounds; Glen visually checked its PNG and rendered
one-page PDF for complete, wrapped, non-overlapping text. JSON/TSV and the three
visual exports are retained; no-clobber refuses and no project state is created.
This is exported-page acceptance on Linux, not native GUI/scrolling or Windows.

**Still red:** renderer filter 2/3 at the tested SHA. The new long-label
regression panics in its test helper: searching `y="` finds the suffix of
`font-family="monospace"` and parses `monospace` as a float. The existing local
test-only repair requires a leading attribute separator; explicit `data-x`/
`data-y` collision controls are now added. Original layout, bounds, Unicode,
complete-text and determinism assertions remain intact. No production renderer,
selection, schema or export behavior changes in this follow-up.

After the next owner push, record the actual repaired SHA and use a run-owned
target. Repeat **only** the renderer filter (three tests, no skips):

```sh
cargo test --locked -p gentle-render -j1 read_end -- --test-threads=1
```

No new CLI binary or BAM/CRAM campaign is needed for this test-only repair.
Retain the prior macOS result with its uncommitted parser/fixture modifications
as attributed evidence, not a Linux pass for the next source. Native GUI lane
heights/scrolling/cancellation and Windows remain `not_run`; normal CI and
scientific acceptance remain separate. Do not repeat completed parity, 47-test,
VM/MCP, glossary, supervisor, Reader/BLAST, containment or barcode campaigns.
Codex performs no local Rust build/test, push or dispatch for this return.
An independent Python reproduction confirms the old parser failure and checks
the separator-based lookup against XML-parsed attributes/bounds on the 16
retained text nodes plus both collision controls. This is not Rust execution.

## Glen Linux Read-End Return And Bounded Repairs, 2026-11-08

Historical `7c5214ec` return; its retry commands are superseded by the active
`1c4a4a5a` follow-up above.

Glen tested `7c5214ec4aa0a28564edf4b7234a1f8386393f90` on Debian Linux x86_64,
Rust/Cargo 1.97.1, Python 3.14.10, locked default features unless noted, one job,
incremental off and an isolated target. The archive SHA-256 is
`560fbf7961d8a9519703bfd329236e1ef23ddb532220ac5e900956c4dde98f37`;
all 100 outer-manifest entries and the 38-entry tutorial submanifest verify.
Source-matched CLI SHA-256:
`696c40660f9dfca7aa06f12e341b16aacc108943f0f06c1760dc52ef345addd7`.
Lockfile SHA-256: `718c89649f21f2c8c437365fd33b256a40dce41eae7097ef1b4e7bc9a54e1314`.

Retain these scoped passes, not a new all-platform or scientific verdict:
canonical generator; registry 29/29 and freshness 1/1 **after generation**;
default test-target/headless compilation; read-end 47/47; renderer 2/2;
supervisor 4/4; async cancellation 2/2; glossary 16/16; script-interface actual
JS/Lua/MCP 29/29; tutorial 07.01 actual-process 6/6. The retained five exports
and no-clobber refusal created no project state. Filters overlap; do not sum.
The already-present local matrix diff is byte-identical to Glen's canonical
patch (only the two read-end rows and adapter counts).

Selection is **not accepted**: the default group is 6 passed/1 failed (including
the unconfigured real-tool early return); the explicit samtools 1.24 test is
0/1. The synthetic CRAM FASTA is 2016 bytes, exceeding its 1024-byte fingerprint
cap. The existing owner/Claude test edit raises that cap, not a production limit.
The native stub's existing test-only 30-second timeout is retained separately
from the production request's explicit timeout.

Source localization corrects the report's binary-input interpretation: in the
BAM branch, the exact `indexed_locus_invalid utf-8 ...` diagnostic can only come
from decoding `samtools --version` stdout, before either `view` command. SAM
decoder errors have the distinct `observed_sam_text_utf8` code. The archive does
not retain raw version stdout, so its invalid byte's upstream origin is unknown.
The new repair preserves such diagnostic bytes losslessly as explicitly labelled
ASCII escapes; valid UTF-8 stays unchanged. SAM and CRAM-header decoding remain
strict. Two new selection regressions cover version bytes, unchanged tool-name
requirements and refusal of malformed SAM through the native stub.

The retained PNG also demonstrates overflowing labels. Shared SVG text now
wraps at fixed, explicit glyph advances; row/header heights grow without
discarding text. Native lanes cache those renderer heights once per loaded
page and remain virtualized, with the same genomic frame. JSON/TSV, scientific
values and denominators are unchanged. Existing raster/output caps remain;
request a smaller explicit page if its taller output exceeds those caps.
One new renderer regression covers 900/1400/8000-px bounds, Unicode, complete
text and lane/full/model layout; the actual-CLI companion now checks SVG bounds.
These new Rust/native changes are authored, **not compiled/run by Codex**.
Local repair checks: five Python source/rehearsal tests and one scoped LF/CRLF
matrix checkout test pass; the actual-CLI companion is skipped without a new
source-matched binary. Documentation/manual/maturity and Rust parsing checks
pass, not a Rust/runtime verdict.

After the next owner push, record the actual integrated SHA and an explicit
run-owned target. Run only the affected scope (require nonzero test counts):

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
cargo test --locked --lib -j1 indexed_locus::tests -- --test-threads=1
GENTLE_TEST_SAMTOOLS=/absolute/configured/samtools cargo test --locked --lib -j1 indexed_locus_optional_real_samtools -- --test-threads=1 --nocapture
cargo test --locked -p gentle-render -j1 read_end -- --test-threads=1
GENTLE_READ_END_BIN_DIR=/absolute/same-source/bin python3 -m unittest tests.test_read_end_rehearsal.ReadEndRehearsalProcess -v
```

Build/provide a same-source CLI externally for that last check; the missing
binary skip is not acceptance. The selection group now authors nine tests on
Unix/eight on Windows, one still opt-in; the renderer filter authors three.
Retain BAM and CRAM results
separately, version-output bytes/hash and actual tool version; CRAM must reach
its reference/FAI/fallback assertions. Review fresh SVG/PNG/PDF at 1400 px and
native lane heights/readability/cancellation where a display is available.
Windows path/handle and full normal CI remain separate. Do not repeat the
completed 47-test/VM/MCP/glossary/supervisor/biological campaigns merely because
this return records them. No private data, live model, package or publication.

### Local macOS Check At 1c4a4a5a, 2026-11-08

Claude ran `1c4a4a5a` on macOS arm64 (rustc 1.101.0-nightly, `--locked`, the
same lockfile hash as above, incremental off). `cargo check --locked --tests`
and the headless `gentle_cli`/`gentle_mcp` check pass. Two failures were test
defects, fixed in uncommitted test-only edits:

- `read_end_long_labels_wrap_without_loss_or_crossing_page_or_lane_bounds`
  panicked with `ParseFloatError`: its attribute helper searched `y="`, which
  first matches the tail of `font-family="monospace"`. The pattern now requires
  a leading space. The layout assertions themselves are unchanged and pass.
- `supervised_process_caps_refuse_or_truncate_without_unbounded_capture` failed
  once and passed on three reruns. The changed `stub.rs` is a new binary whose
  first launch can exceed the fixture's 5 s default; `limits()` now uses 30 s.
  The timeout test keeps its explicit 50 ms limit.

After the fixes: renderer `read_end` 3/3, `supervised_process::tests` 4/4,
`indexed_locus::tests` 9/9 (the real-samtools test returned early without
`GENTLE_TEST_SAMTOOLS`), `read_end` 47/47, script-interfaces `read_end_` 29/29,
registry 29/29, freshness 1/1, and the actual-CLI rehearsal 6/6 with the new SVG
bounds assertion against a same-source `target/debug/gentle_cli`. Real BAM/CRAM,
fresh PNG/PDF review and native lane heights remain for Glen/native review.

## Read-End Compile-Gate Return And Retry, 2026-11-08

Historical return and retry, superseded by the successful generator run above.

Glen tested `1d5eed33bcd58460f93dd45f90d85588e9fb5859` on Debian Linux x86_64
with Rust/Cargo 1.97.1, Python 3.14.10, one Cargo job, incremental disabled and
an isolated target. The mandatory generator exited 101 **during compilation**:
nine diagnostics identify three causes, not nine independent defects. The
empty retained matrix diff does not establish freshness. All downstream
parity/read-end/service/VM/MCP/tutorial/native gates remain **not_run**.

The supplied archive's SHA-256 is
`c8e86e19fe373f647c98797ca11cc18181794d7ed29507bf4bb7cd4c0663e676`;
all 11 manifest entries were independently verified, and its `REPORT.md` matches
the separately supplied report byte-for-byte. The tested lockfile SHA-256 is
`7ad5b2dc9485675a75d272fa945c75342dceb60ef582f2885015f3e14df41e76`.

The owner's newer commit `6adee466` ("Claude fixing build failures") already
adds the GUI's `quote_shell_arg` import, declares coverage's accumulator as
`Option<Vec<Interval>>`, and clones the publication identity while cancellation
still borrows it in Build/Export. It also removes the reported unused import.
These source corrections have been inspected against the raw diagnostics;
**no successful compilation or acceptance of the repaired source is recorded**.
Do not substitute the historical failed SHA or infer a pass from these edits.

After the owner pushes the repaired integrated source, record its actual SHA
and use a fresh run-owned target. Retry this command **with default features**
first, so the GUI import is compiled too:

```sh
cargo run --locked --quiet --bin gentle_examples_docs -- parity-matrix-generate
```

If it fails, retain the primary compiler diagnostics and stop again. If it
succeeds, return only the actual canonical matrix diff, then run the pending
parity/freshness and new read-end scope below, including test-target compilation.
Keep Windows, native GUI and real configured samtools separately reported.
Do not repeat completed Python/reader/BLAST/containment or real-barcode analyses.

### Local macOS Execution, 2026-11-08

Claude ran the repaired source at `7c5214ec` on macOS arm64 (Darwin 27.2,
rustc 1.101.0-nightly, `--locked`, Cargo.lock SHA-256
`718c89649f21f2c8c437365fd33b256a40dce41eae7097ef1b4e7bc9a54e1314`, bio 4.2.2)
with two uncommitted test-only edits to
`src/observed_transcript_structures/indexed_locus_tests.rs`:

- fixture `timeout_ms` 2000 -> 30000: the first launch of each freshly compiled
  stub exceeded 2 s (`samtools version exceeded the 2 second timeout`), which
  also masked the expected `tool_refused`/`byte limit` messages;
- FASTA/FAI fingerprint cap 1024 -> 1 MiB: the 2016-byte synthetic FASTA
  deterministically failed `indexed_locus_hash_byte_budget` on every platform.

Results (exit 0 unless noted): `cargo check --locked --tests`; headless
`--no-default-features` lib and `gentle_cli`/`gentle_mcp` (26 pre-existing
dead-code warnings); `-p gentle-render --tests`; `--features script-interfaces
--lib --tests`. Tests: `read_end` 47/47, `supervised_process::tests` 4/4,
`indexed_locus::tests` 7/7 twice (the real-samtools test returned early without
`GENTLE_TEST_SAMTOOLS`), `async_blast_probe_cancellation` 2/2, `glossary_` 16/16,
`gentle-render read_end` 2/2, script-interfaces `read_end_` 29/29,
`capability_registry_parity` 29/29 and `parity_matrix_freshness` 1/1 after the
canonical `parity-matrix-generate` (two new rows: `rna-reads read-ends`,
`RnaReadEnds`). The actual CLI rehearsal
`GENTLE_READ_END_BIN_DIR=target/debug python3 -m unittest discover -s tests -p
'test_read_end_rehearsal.py'` passed 6/6 without a skip. The first parity run
failed only because the disk was full (linker); it passed after freeing space.

Still open: Linux and Windows CI at the pushed SHA, real configured samtools,
native GUI responsiveness/cancellation, Windows handles and full normal CI.

## Native Read-End Viewer And Migrated 07.01, 2026-11-08

Historical authoring checklist; use the active bounded-repair scope above rather
than repeating the now-completed Linux runtime campaigns.

Steps 1-6 are source-implemented, not runtime/native accepted. Preserve OD1-OD6.
No local Rust build/check/test, live model, private data or screenshot capture
was run. The GUI is a thin reviewed-request/virtualized shared-renderer tab;
the technical file/config bindings remain explicit in the JSON draft. It starts
without a selected transcript, sample or report. Channel toggles are in the
hashed view, never a combined score. Both denominators remain outside scrolling.
Local checks at this authoring boundary: 51 existing RNA Python tests, five new
rehearsal/source tests and 87 checkout tests pass (143 total); the one actual-CLI
test is explicitly skipped without a supplied binary. Manual/status/maturity
checks and scoped Rust parsing/formatting pass. These are not Rust execution.

The existing docs binary reports `2d2b8b42`; its generator sources are unchanged
against this integration. It generated/checked the current 70-entry catalog and
checked the unchanged 30-chapter manifest. **It cannot generate the newer embedded
registry/parity projection.** Run the canonical source-matched generation in the
Step 5 section first and return that minimal matrix diff. This pending gate is
intentional and must not be represented as CI-green or repaired by hand.

Then run only the new read-end scope, with an explicit temporary target:

```sh
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
cargo test --locked --lib -j1 read_end -- --test-threads=1
cargo test --locked -p gentle-render -j1 read_end -- --test-threads=1
cargo test --locked --lib --features script-interfaces -j1 read_end_ -- --test-threads=1
GENTLE_READ_END_BIN_DIR=/absolute/same-source/bin python3 -m unittest tests.test_read_end_rehearsal -v
```

The active core filter now includes the new cooperative hook regression (21,
not the historical 20-test count below). New actual VM/MCP filters are
`read_end_js_runtime_`, `read_end_lua_runtime_`, `read_end_mcp_discovery_`;
new GUI helpers use `read_end_gui_`. Require actual execution/nonzero counts,
not compile-only evidence or the one missing-binary process skip. Check the
unchanged generic MCP confirmation gate even for inspect, request-aware shell
approval, file-only managed work, exact page equality, no project/file mutation,
and explicit cancellation versus terminal result. No live model is required.

Replay 07.01 with synthetic inputs first. Native acceptance must separately
observe explicit selections, stacked hypotheses, both unchanged denominators,
all separate channels, paging/large-lane responsiveness and usable cancellation.
Verify report/view cache identity, export after filtering, exported report not
replacing the source pileup, save/reopen and native Windows closed handles/
spaced/apostrophe/Unicode paths. Renderers/converters are not preempted inside
their bounded calls. Keep actual configured-samtools controls, full normal CI,
packages, private-data accuracy and scientific acceptance separately attributed.
Do not repeat completed reader/BLAST/containment or real 15/16 campaigns.

## Shared Read-End Services And Export, 2026-11-08

Step 5 implements the shared operation/parser, file authentication, no-clobber
exports, indexed-selection S2.3 seam and managed artifact-only publication.
Six service and two managed-policy tests are authored, Rust-unrun locally.

First generate **only** the canonical projection on the integrated source:

```sh
scripts/regenerate_parity_matrix.sh
python3 scripts/generate_capability_maturity.py --check
cargo test --locked --lib -j1 rna_read_ends::tests -- --test-threads=1
cargo test --locked --lib -j1 read_end_policy_tests -- --test-threads=1
cargo test --locked --lib -j1 glossary_ -- --test-threads=1
cargo test --locked --test capability_registry_parity --test parity_matrix_freshness -j1
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
GENTLE_READ_END_BIN_DIR=/absolute/same-source/bin python3 -m unittest tests.test_read_end_rehearsal -v
```

Return the canonical generated matrix diff for integration; no local Rust
generator could run under the owner restriction, and no projection was invented
by hand. Maturity sidecars were generated via their Python generator. All Rust,
VM, MCP runtime and native behavior remain unverified at this source. Require
nonzero counts and actual CLI execution, not the intentional missing-binary skip.
The rehearsal helper creates only synthetic files and command drafts, no tool
launches. Check model/annotation raw hashes, saved-only inspection, explicit
sample/transcript selection, unchanged state/history, approval before I/O,
request-bound markers and late-completion retention, exports/view hashes and
Windows close/rename/no-clobber behavior. PNG/PDF must derive from the shared SVG;
test dimensions/derivation, not cross-platform raster equality. Keep unsupported
composition/tail evidence unavailable and incomplete region search explicit.
Exercise generic JS/Lua/Python and MCP routes on same-source binaries; do not
repeat prior unrelated campaigns or private-data/real-barcode analysis.

## Bounded Configured Alignment Selection, 2026-11-08

Step 4 implements the shared direct-child supervisor and owned private SAM/
manifest selection seam. No tool, Rust build or private input was run locally.
Authoring includes four supervisor tests, six portable selection tests (one
optional real-tool test) and one Unix symlink control:

```sh
cargo test --locked --lib -j1 supervised_process::tests -- --test-threads=1
cargo test --locked --lib -j1 indexed_locus::tests -- --test-threads=1
cargo test --locked --lib -j1 async_blast_probe_cancellation -- --test-threads=1
```

Use an explicit run-owned target. The stub is a synthetic native executable,
not binary BAM/CRAM evidence. On a compatible host set `GENTLE_TEST_SAMTOOLS` to
an absolute configured executable and rerun the optional test; it constructs
real synthetic BAM/CRAM/FAI/indexes locally. Without that variable its early
return is **not real-tool acceptance**. Verify declared versus verified input
hashes, incomplete unit manifests, cancellation/caps/timeout and closed handles,
missing/wrong/stale index refusal, CRAM reference/hash/FAI requirements and the
header-UR/environment fallback controls. Return actual counts/tool versions and
Windows handle/path results separately. Reuse existing narrow DuckDB process
tests if available; do not repeat completed biological or full BLAST campaigns.
The supervisor owns/reaps the direct configured child, not hostile descendants.
Index pairing is explicitly declared; a regional pass cannot prove every index
offset or verify an unrequested whole BAM/CRAM digest.

## Read-End Saved Pileup, 2026-11-08

Step 3 adds the pure engine report/view and headless SVG, not yet a file/GUI
route. Six root and one renderer tests are authored, Rust-unrun locally:

```sh
cargo test --locked --lib -j1 read_ends::pileup::tests -- --test-threads=1
cargo test --locked -p gentle-render -j1 read_end_pileup -- --test-threads=1
```

Check stacked hypotheses within one unit, both unchanged denominators, explicit
filter/paging budgets, stable ties, annotation-checked selected geometry and
distinct N/deletion/segment connectors. Inspect both transcript directions,
off-axis clip descriptors and unavailable composition, retained candidate/row
provenance and deterministic SVG/TSV. No source SAM/BAM is reread. PNG/PDF,
authenticated services and native acceptance remain later gates. No existing
data campaign or prior unrelated acceptance needs repeating.

## Read-End Coverage And Separate Evidence, 2026-11-08

Step 2 adds pure saved-projection coverage/denominators and separate composition/
tail-import records. It is not yet a file importer or command route. There are
five authored coverage tests and three evidence tests, all Rust-unrun locally:

```sh
cargo test --locked --lib -j1 read_ends::coverage::tests -- --test-threads=1
cargo test --locked --lib -j1 read_ends::evidence::tests -- --test-threads=1
```

Check gene-assigned versus selected-transcript denominators, all four exclusion
buckets, unique/shared/hypothesis-dependent membership, exact-only quantiles,
the two censoring fractions and conditional min/max coverage without duplicate
hypothesis counting. Empty denominators remain unavailable fractions, not zero
measurements. Refuse repeated sample/gene analyses and bad source bindings.

The A/T schema intentionally reports `no_bound_basecalls` under OD2, retaining
explicit orientation/trimming and required parameters; it never scans or claims
a tail. Tail imports require the exact configured caller/version/model/kit/
basecaller/parameter hash and manifest declaration, per-unit source and saved
QNAME witnesses. Unknown configurations and failed QC produce no accepted tail
value. Invalid unit/read/source/parameter bindings refuse. Pure functions do
not authenticate files or POD5/FAST5; that remains the later service boundary.
No UTR status changes, conversions between callers, BAM tag reads or private
analysis. Preserve OD1-OD6 and prior receipts; return actual tested source and
nonzero counts, not a claim that Python checks executed these Rust functions.

## Read-End Core Refinements, 2026-11-08

This is the active step-1 filter after the [integration-plan reconciliation](rna_read_end_integration_plan.md).
The original `1df2db13` receipt below is historical. The core now requires the
parsed normalized annotation; the caller must authenticate its exact file hash
against the saved binding. Confirm geometry/gene mismatch refusal, uncollapsed
hypothesis/candidate records, both gene-ID fields and additive clip descriptors.
There are 20 authored Rust tests, including the original 31-case geometry corpus
and six clip cases; no Rust compilation/execution occurred locally.

```sh
cargo test --locked --lib -j1 read_ends::tests -- --test-threads=1
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
python3 -m unittest discover -s tests -p 'test_*rna_read_end*.py' -v
```

Require 20 nonzero Rust tests; compare the unmodified collapsed literals as well
as new metadata, including minus/reverse clip sides and missing-CDS gene identity.
The original normalized-annotation hashing seam was inspected in both S2.3
branches before implementation. No completed BLAST/reader/private analysis or
containment campaign needs repetition. Return exact source/toolchain/features,
counts/skips/failures and logs; no new native or biological acceptance follows.

## Shared UTR Core Step 1, 2026-11-08

Historical authoring receipt at `1df2db13`: the active refined API/checklist is
above. Its 17-test count and authoring baseline below do not describe later code.

Verify the integrated revision containing this first Rust core after the owner's
normal push; do not create an interim freeze/bundle or repeat completed data
campaigns. Authoring baseline is `1b73f9dc` on main. The scoped change adds only
`observed_transcript_structures::read_ends`, its registration, synthetic tests
and documentation. No Operation, CLI/MCP/GUI route, Cargo dependency/lock change,
coverage, BAM/CRAM reader, terminal A/T scanner or physical-tail caller is added.
Registry/parity projections are intentionally unchanged, not hand-regenerated.

The pure core binds strict V2 models to existing admitted SAM assessments and
reuses the shared block-coordinate transform. All retained compatible candidates
and hypotheses must agree for each end. Missing CDS models stay visible; lower
bounds, unknowns and true zeros have distinct typed representations. Stop
intervals are explicit and raw CDS remains unchanged. V1 is refused by this new
Rust boundary, not silently interpreted as V2; the standalone Python legacy
reader and all its outputs remain unchanged. Source-file authentication remains
the responsibility of the later file/service boundary.

There are **17 authored Rust tests, not compiled or run locally**. One consumes
all 31 literal shared corpus cases; others cover strict shapes/bindings,
duplicate/cross-unit/nonoptimal/axis guards, exact/one-under work limits,
large coordinates, candidate/path ordering and actual existing SAM-kernel
assessments. The corpus's Python half checks both LF and CRLF against the
unchanged projector. Python/source checks do not establish Rust equivalence.

Local authoring checks at `1b73f9dc` plus this scoped diff: 50 Python read-end
tests (including the 31-case corpus in both encodings), 19 documentation/checkout
tests, documentation status and scoped rustfmt/whitespace checks pass on macOS
arm64 / Python 3.15.0. No Rust build, check or test was run. Existing owner
changes to `Cargo.toml`, `paper`, private data and untracked files are excluded.

Using an explicit run-owned target, run or reuse only matching normal-CI builds:

```sh
cargo test --locked --lib -j1 read_ends::tests -- --test-threads=1
cargo check --locked --tests -j1
cargo check --locked --no-default-features --lib -j1
python3 -m unittest discover -s tests -p 'test_*rna_read_end*.py' -v
```

Require nonzero Rust selection (currently 17 tests) and inspect the 31-case corpus
assertion, not just successful compilation. Confirm inputs remain unchanged,
missing candidates cannot disappear, both strands/split stops agree, and
observed zero/censored-zero bounds retain their distinct meanings.
Return exact tested SHA, platform/toolchain/features/profile, command counts,
failures/skips and retained logs. Preserve Glen's focused Linux containment pass
at `1cf7e120` below; it does not accept this new core and need not be repeated as
a standalone campaign. No BLAST/reader, private 15/16, GUI, package or biological
rerun is requested. Current-stage semantic mismatches must be reported rather
than fixed by relaxing the shared literals or existing Python contracts.

## Observed-Structure Containment Conflict, 2026-11-07

### Glen Linux Return, Received 2026-11-08

The owner supplied `REPORT (5).md` and
`gentle_verification_1cf7e120_20261107.tar.gz` for exact source
`1cf7e120352e46c2761ae77694892d87a48a4260`. Codex checked all 35 manifest-listed
files; `SHA256SUMS` itself is the only unlisted file. The separate report is
byte-identical to the archived report. Archive SHA-256:
`85f7d2155d3c8e991506ddde6d529c21843fe6fcabd830c29cde7dbe56d37a99`;
manifest SHA-256:
`0939d9290c6dce12beac0c36b7f429733f300ea82f0c79fef6d8448376817bad`.
The recorded lockfile hash and changed-Rust-path fingerprint below match Git
at that source. The test-binary hash is recorded in the evidence, not locally
recomputed from a supplied binary.

Raw commands, stdout and exit receipts establish **15 passed, zero failed or
skipped** on Debian forky/sid x86_64, Rust/Cargo 1.97.1, default dev/test
profile, one job and an explicit run-owned target: eight
`observed_structures_containment_`, one `sam_containment_`, four frozen
`observed_structures_s2a_`, one existing partial-chain/primer-boundary control,
and one legacy primer-payload byte/work preservation test. All five commands
exit zero. The clean detached checkout and unchanged final remote tip are
recorded before/after the run. The existing unused `lints.cargo` manifest-key
warning is not a containment failure.

This closes the focused Linux execution request below; do not queue those same
checks again as never run. The original authoring record is retained below as
history. Normal CI and any affected later-source tests remain separate, as do
Windows, native GUI, full-suite and scientific acceptance. No new BAM/CRAM,
UTR-projection or Poly(A) capability is established. No local Rust execution,
push, dispatch or private-data replay was performed to admit this receipt.

### Original Authoring And Verification Request

Add these authored regressions to the next already-planned consolidated Glen/CI
run after the owner's normal integration/push. No separate campaign, dispatch,
real-data replay or local Rust run is requested. This source is the uncommitted
follow-up over exact HEAD `6c1edcf3e4bf0e765f755e5b0cf5b64866925b25`; fetched
`origin/main` was `2d2b8b42e3610ef1874f63057b24b92f520a4e04`.
The supplied read-only review's core/comparison/DAG shapes at `e4bc7386`
were unchanged at that baseline. Do not attribute new execution to the bare
baseline or inherit historical S0/S1/S2 results as a pass of this fix.

The private shared core now rejects containment when a supported genomic block
in any same-record shorter-chain observation overlaps a missing longer-chain
skip. Earlier `ContainedPartial` labels in these conflicting cases were wrong;
they become `Alternatives`, without changing possible/definite chains, physical
chain IDs or Phase 2 behavior/bytes. Touching half-open ends and unsupported
`D` gaps remain nonconflicting. The record-local witness index also checks
nonoptimal observations of an already-possible chain without promoting them
into optimal hypotheses or counts. New witness/skip/block comparisons use the
existing total-work budget; observed-structure work counters can increase when
these checks run. Do not regenerate frozen preservation expectations blindly.

Nine inline synthetic tests are **authored but not run**, including compilation:

- `observed_structures_containment_` selects eight new cases: splice/retention
  conflict, right-flank zero/one-base overlap, reverse-query/minus-model genomic
  coordinates, left-flank zero/one-base overlap, record isolation/every witness,
  nonoptimal same-chain witness, deletion-gap control and exact/one-under budget.
- `sam_containment_` selects the new two-row, one-declared-unit test on both
  query orientations, keeping both equal-weight hypotheses possible/neither definite.
- Preserve the unchanged
  `observed_structures_partial_chain_and_primer_bounded_ends_do_not_become_full_rna`
  regression, the `observed_structures_s2a_` frozen bytes/IDs/work/order guards,
  and `primer_payload_s0_legacy_serialized_rows_and_work_are_byte_identical`.
  Static inspection found those frozen cases have one family or different-contig
  alternatives, so the new same-contig witness checks do not alter their work;
  this is reasoning, not an executed verdict.

Use the normal consolidated build/check scope and add these filters to its
existing run. If its current executable already includes this final source,
reuse it with the repository's configured test stack; do not use an older binary.
Record actual commit, platform, toolchain, scope/results and primary failures.
Windows/Linux/native and scientific acceptance remain unrun here. The local
checks for this change are formatting, documentation and whitespace only.

Changed-Rust-path fingerprint: `ab4beb68d50d28b412cc4fa048f87a2d7c9d2f76186d62e8b94d787059a3e5a3`.
Recreate by SHA-256 of the sorted UTF-8 concatenation
`<path> NUL <lowercase raw-file SHA-256> LF` for exactly
`src/observed_transcript_structures/{core,tests,sam_tests}.rs`.
This binds the authored source, not executed behavior or a final integration SHA.
Selected-gene relation clarification and the corrected future render/lane seam
are documentation-only; no public schema or relative-to field is added.

## Reported 15/16 Panel Analysis, 2026-11-07

The owner relayed Glen's already completed regular-cDNA comparison at
`7f667f4d1474fa6c4400e3cd99fb59e82240c4cf`, including `cf4f59d3` and a newly
prepared Ensembl-116 v2 model from the frozen GTF. Do not ask Glen to repeat
the initial 15/16 analysis. This was a prespecified 80-gene panel audit through
a run adapter using GENtle's unchanged UTR geometry on validated panel BAMs,
not an integrated full-data BAM/CRAM-to-structure-to-UTR execution.

These are owner-forwarded reported results only. No adapter source, commands,
input/model hashes, raw tables or result manifest accompanied this message;
Codex has not independently verified the biological data or adapter equivalence.
Retain those existing small provenance records and per-gene gate tables when
available, without requesting new large transfers or relabelling old results
as produced by the later ranking repair.

| Reported cohort | Observed 3-prime | Right-censored | Missing | Total |
|---|---:|---:|---:|---:|
| SK-MEL-28 / barcode15 | 946 | 29568 | 84381 | 114895 |
| WM793 / barcode16 | 2871 | 66704 | 181235 | 250810 |

Glen reports no gene passing the prespecified gate of at least ten observed
values per cell line and at most 25% censoring, without post-hoc relaxation.
Observed counts are TP73 6/0, MDM2 10/5 and MDM4 0/4 (15/16); E2F1, TP53,
TP63 and PATZ1 were reported zero. Observed is a metric status, not a complete
molecule/UTR or polyadenylation-site confirmation. A no-eligible-gene outcome
is not a product failure or permission to loosen the thresholds.

The reported 25.735% / 26.595% censored shares use all three statuses as the
denominator. The unchanged ranking gate uses **censored / (observed + censored)
per gene and sample**, excluding missing values. For the aggregate counts above
that second fraction is approximately 96.900% / 95.874%; these global fractions
cannot themselves prove every per-gene gate failed. Preserve the actual
per-gene gate table rather than substituting library-wide percentages.

Reported median whole-read lengths are 1549 / 1110 bp for the analyzed 15/16
cohorts. Keep the panel/read-selection scope explicit; do not extrapolate to
all barcode reads without the input-universe receipt or interpret a shift as
gene-specific UTR shortening. Glen's separate terminal A/T signal was weak and
threshold-dependent, not calibrated poly(A) presence or tail length. Causality,
biological replication and the later 14/13 and 18/17 contrasts remain separate.

Remaining bounded action: after the owner's push, run only the five ranking
tests below at a source containing `2d2b8b42`, then read-only validate the
retained detailed tables with the shared `load()` reader and preserve their
hashes. The repair changes malformed-input admission, not valid-data scoring.
If the original tables pass, do not reread BAM/FASTQ, rebuild models, reproject
or recreate the analysis merely for a new SHA. Keep the original analysis
attributed to `7f667f4d` plus a separate validation receipt. On refusal, report
the offending row/status and stop for a scoped correction decision; no silent
dropping, evidence edits or overwriting historical outputs. A scalable complete
intake and calibrated poly(A) filter are separate development proposals,
not work authorized or acceptance achieved by this report.

## UTR Ranking Row Validation, 2026-11-07

Glen's owner-forwarded scan at exact
`7f667f4d1474fa6c4400e3cd99fb59e82240c4cf` reports 45/45 read-end tests and
1/1 Rust `parity_matrix_freshness` passing. Preserve that source-bound Linux
return; it is not a Rust run of this subsequent Python-only repair. No new
logs/archive were supplied with this scan, so those execution results remain
attributed to Glen's report rather than independently revalidated here.

The scan found two malformed detailed-row cases accepted by the common
`rank_rna_read_end_genes.load()` reader: `missing/17` and a blank
`right_censored` bound. Codex reproduced both. The repaired loader validates
status/value before gene filtering or counts: observed/censored rows need a
nonnegative integer, missing rows an empty cell. Both ranking routes use that
one reader. True observed zero and censored zero remain valid; lower bounds
are not median observations. No threshold, contrast, nomination, model schema,
four-column comparator or generated-fixture output is changed.

Focused regressions exercise LF/CRLF, spaced/apostrophe/Unicode paths,
unavailable/multiple gene IDs, both public ranking writers, refusal before
output creation and unchanged input bytes. The ordinary Python CI discovery
already includes these tests. Local repair checks on macOS arm64 / Python
3.15.0 pass 48 read-end tests and 20 documentation/checkout tests (68 total,
no failures/skips), including unchanged synthetic generation under LF/CRLF.
Manual/status and whitespace checks pass; the roadmap stays at 400 lines.
These are worktree checks based on `7f667f4d`, not a new Rust or Windows pass.
Glen need not repeat completed reader, BLAST,
45-test or parity campaigns. For an independent focused check after the fix
is committed/published, only this small target is needed:

```sh
python3 -m unittest discover -s tests -p test_rank_rna_read_end_genes.py -v
```

Normal CI can run the whole read-end scope on Linux/Windows. The 15/16 panel
analysis is already reported completed above: validate its retained tables,
do not request a fresh full analysis. Do not silently discard invalid rows or
edit original evidence to make a run pass. This repair adds no private-data or
biological approval; existing explicitly authorized analysis stays separate.
No local Rust, push or dispatch is performed here.

## Slice B Panel Multi-Reference Specificity

Glen's initial 2026-11-06 return confirmed the projection defect at
`27a8de33d1b349f6c0dfd9eef75e799e066469f6`: 14 targeted Python checks passed,
but the locked headless `parity_matrix_freshness` run failed, 0 passed / 1
failed, because `docs/gui_cli_mcp_parity.md` was stale. Preserve that failure
as the pre-repair result, not a defect in the later UTR work.

The owner supplied Glen's subsequent `GENtle_27a8de33_parity_handoff.zip`,
23,333 bytes, SHA-256
`667043a8aaba7be3f5edbc8138ed8d75632b82840a25e7102a82c9690671cbd2`.
Codex verified exact membership and all 13 manifest hashes; the manifest hash is
`25a9516c6104b5b847d1de51ec4d534ad44a87dbd555b8f66443d36f0cf03774`.
The report and patch match the earlier compact verification ZIP byte-for-byte.
The source/tree/lock match local Git records. Protocol crate, generator, lock
and embedded glossary/override/collection-policy inputs have no diff from
`27a8de33` through pre-import `34660ecf`, so no new interim source freeze is
needed for this generated artifact.

Import only Glen's original generator patch (17 insertions / 5 deletions),
SHA-256 `23ef83413e1973bd4474d5563ebfbf7e380bfa71d5b0f998d0da79275b1328b3`.
The resulting tracked matrix is byte-identical to his supplied 299,119-byte
file, SHA-256
`03b3bd99e1342409fd6a4a0130f20355d25bd96fcf1b88de45053fdb1e728c55`.
No rows/counts were hand-edited. Generator/check logs and exit-zero receipts
are supplied. The retained Linux test stream records 29/29 capability parity
and 1/1 freshness after generation at `27a8de33` plus this Markdown-only diff.
This closes the missing-artifact request and imports the canonical repair;
it is not a local Rust rerun or acceptance of the later integrated source.

Glen's report additionally records locked check, 10 panel/script-interface,
10 multi-reference, 1 panel-history and 13 matrix passes; protocol was 119/1
before generation and 120/0 afterward. It reports 55 examples and 17 Python
maturity checks passing. Those broader execution results are attributed to
his report; the compact ZIP supplies raw logs for the generator/check and
29+1 parity scope, not the full broader campaign. No additional evidence or
repeat Linux campaign is requested. Windows/full current-source acceptance
waits for normal CI after the owner's push; no local build or dispatch occurs.

Local import checks: 37 Python metadata/documentation/checkout tests passed,
including exact LF/CRLF matrix preservation and its missing-policy control.
Manual/status and whitespace checks passed; the roadmap remains 400 lines.
The matrix still matches Glen's exact file hash. These checks do not execute
the Rust freshness test or replace the normal integrated-source CI gate.

Source work is based on `0d9909b5`, following the bounded single-pair Linux
return below. Use the final Slice B SHA/tree/lock, not the older reader binary.
Only synthetic panel evaluation is authorized: 24 assays, 8 references, 384
commands, 128 MiB cumulative TSV and 25,000 rows. These are input limits, not
measured runtime or RSS limits. No repeat ACTB/TP53 campaign or new BLAST run.

The Linux scopes below are now reported completed by Glen at that source.
Keep these reproduction commands for affected future source/normal CI, not a
request to repeat his completed checks. Retain platform/toolchain/features,
nonzero counts, skips and primary failures when execution is actually needed:

```sh
cargo check -q --locked --tests -j1
cargo test --locked --lib --features script-interfaces transcript_assay_panel_multi_ -- --test-threads=1
cargo test --locked --lib engine::primer_specificity_multi::tests:: -- --test-threads=1
cargo test --locked --lib transcript_assay_panel_specificity_ -- --test-threads=1
cargo test --locked --lib specificity_matrix_ -- --test-threads=1
cargo test --locked -p gentle-protocol
```

Exercise scripting cases with the JS/Lua features enabled; report zero matches
or platform skips separately. The mock BLAST end-to-end and historical panel
finalization fixtures are Unix-only; Windows must exercise pure admission,
manifest, byte/row, quoted-path and rollback cases as well as the outstanding
reader checks. Partial/cancelled is not an assay pass. Do not backfill sources
or attach to the panel to make a test pass.

Glen ran the canonical generator and the returned Markdown is now integrated.
Normal CI should check the **integrated** source without regenerating first:
`cargo run --locked --quiet --bin gentle_examples_docs -- parity-matrix-check`,
`cargo test --locked --test capability_registry_parity --test parity_matrix_freshness`
and the separate examples check `cargo run --locked --bin gentle_examples_docs -- --check`.
Do not use the old local `f3d8f86f` docs binary to certify the new registry.
Local source/format/metadata checks do not establish a Rust or scientific pass;
the imported matrix and historical Linux verdict do not claim a green new CI.

Earlier authoring checks on macOS arm64 (historical): 14 focused offline primer/tutorial/reader
checks passed. The wider 96-check selection had 89 passes, five skips and two
failures: the then-pending parity projection and an unrelated 10.04
RNA source-text assertion about command admission, reproduced on untouched
`0d9909b5`. Neither guard was weakened. The capability-maturity projections
were refreshed via the Python generator, not the Rust parity generator.
The shared pair reader still validates the database before interpretation and
again before retaining the report; its source guard follows the extracted call
path. The new import commits through the store serializer without the ordinary
writer's panel-source backfill. The five feature/adapter modules pass scoped
rustfmt checks; comparing 17 larger touched files against the base finds the
same 58 pre-existing formatting hunks and no new ones. Manual catalog,
documentation status, maturity freshness and whitespace checks pass. The parity
failure is superseded by the canonical return above; the separate 10.04 source
assertion is not automatically closed by this documentation-only import.

## Linux Reader And MCP Return At e4bc7386, 2026-11-05

Glen reports a scoped Linux pass at exact
`e4bc73868ee691aa87bf7968f431e6fb71c5d599`, tree
`02ae61eb49f55c8aaeee12b641205898bed47740`, Rust/Cargo 1.97.1,
ordinary dev/default features, one job, `CARGO_INCREMENTAL=0`. Default/test
and headless locked checks, CLI build, 9 reader + 2 equivalence + 10
multi-reference + 5 trace + 13 matrix tests passed without failures/skips.
Full capability parity 28/28 and freshness 1/1 passed. Separate framed MCP
discovery recorded 49 tools / 99 generic `op` paths with all four repaired
RT commands. The queued Linux MCP/parity request is therefore closed at this
source, not queued for another unchanged run.

The owner supplied [Glen's original report](audits/glen_followups_20261008/reader_verify_e4bc7386/REPORT.md):
6,219 bytes, SHA-256
`c8d83113058f4adfe31097db62bd361d9f2aee0d8bf57018afac353f6354f676`.
It is retained byte-identically with scoped LF/CRLF checkout protection.
Codex verified those report bytes and local Git source/tree/lock identities;
the sealed archive hash, binary identity, raw-run results and input/reference
checks remain attributed to Glen. The reported archive SHA-256 is
`7f23ce57c9938e2756d75a1feb36f5ced5c9163cf1fda0e3a4c47b8bf7810f15`;
that bundle and its raw members were not supplied locally or independently
rehash-verified. See [provenance and limits](audits/glen_followups_20261008/reader_verify_e4bc7386/README.md).

All four retained imports completed with exit 0 and reopened states: ACTB
untraced/traced 41.39 / 35.08 s, TP53 58.26 / 42.55 s. Both traces contain
52 begins and 52 completions, no unfinished/failed/truncated spans; the formerly
blocking forward-hit phases complete. Glen reports 213/213 unchanged input
checks and byte-identical reference receipts. Trace-mode output/state/reopen
comparisons match after excluding only `generated_at_unix_ms` and
`updated_at_unix_ms` and canonical JSON sorting. This is same-source tracing
equivalence, not real old/new report equivalence: `dfca3261` produced no report.

Do not infer an exact speedup ratio from censored old timeouts, faster tracing
from non-counterbalanced runs, or general throughput. TP53 untraced has only
1.74 s nominal deadline headroom; peak RSS is about 1.07-1.08 GiB, above the
earlier incomplete runs. Record the resource tradeoff without attributing all
of it to one allocation or opening a new optimization campaign automatically.
Exit-zero import completion is not a specificity verdict or scientific approval.

Remaining automated reader acceptance waits for normal Windows GitHub CI:
its full workspace suite includes the nine reader and two equivalence tests.
Glen's Linux-only host needs no Windows device or duplicate run. Retain the
actual Windows test outcomes for LF/CRLF, paths, timestamp restoration,
missing/corrupt inputs and handle closure; do not infer them from Linux.
No repeat Linux replay, new BLAST/reference preparation, package/GUI/live-model
work, push or dispatch is requested. Keep the published UTR stack separate:
this report neither reviews nor merges it, and no reader code was changed.

## Standalone UTR Boundary Repair, 2026-11-05

This repair is based on Glen / smoe-bot's three-commit UTR stack ending at
`3ff288027ba3bede9a3aeaadb951459f4d418832`. Preserve that authorship. Codex's
review reproduced four counterexamples; the owner-forwarded Claude review
independently confirmed the censoring, candidate-loss and anchor defects and
the stop-convention mismatch. Claude's reported 29 original passes and nine
probes are attributed evidence, not tests rerun by Codex on the repaired source.

Only the standalone Python boundary is changed: censored metrics are blank
for the unchanged comparator while detailed bounds remain numeric; full-path
candidate intersections precede model availability; connected coding-side
anchors preserve true zero; model/report v2 carries original CDS and explicit
stop evidence. Separate/included/overlapping/split stops normalize in transcript
coordinates. Partial/absent or legacy v1 stops give unknown 3′ reach, not `+3`.
See [the decision](decisions.md#standalone-transcript-read-end-projection) and
[model/recovery instructions](rna_read_end_metrics.md).

Authoring verification on macOS arm64 / Python 3.15.0: 45 focused read-end
tests, 86 documentation/checkout tests and 51 release-candidate helper tests
passed, zero failures or skips (182 total). This includes 16 new focused
regressions and the real Python CLI chain. Manual/status and whitespace checks
pass; the roadmap remains 400 lines. These are repaired worktree results based
on the named Glen tip, not an integrated-main or native Windows verdict.

Repair commit `cf4f59d3` is now owner-authorized for integration into main on
2026-11-06 together with Glen's original commits. Keep the branch authoring
results above distinct from the integrated-source checks recorded below.

Merge-result working-tree checks on 2026-11-06, macOS arm64 / Python 3.15.0:
45 read-end tests, 87 documentation/checkout tests and 51 release-helper tests
passed (183 total, no failures/skips), with manual/status and whitespace checks.
Both documentation conflicts retain main's newer panel/reader records and the
UTR request; no Rust source or dependency changed. Canonical Rust tutorial/
parity execution remains unrun locally; do not substitute the Python checkout
or synthetic-figure regeneration for those gates. Native Windows waits for
normal CI, and historical scientific/release limits remain unchanged.

Glen/CI: on the integrated source, run only the affected Python scope
and actual synthetic subprocess chain; do not repeat reader, BLAST, comparator
Windows-at-108cc5e5, gel or private six-barcode analysis campaigns:

```sh
python3 -m unittest discover -s tests -p 'test_*rna_read_end*.py' -v
python3 -m unittest scripts.test_tutorial_checkouts.TutorialCheckoutTests.test_rna_read_end_tutorial_hashes_and_generation_survive_both_checkout_modes -v
```

The boundary suite includes the real preparer/projector/comparator subprocesses,
LF/CRLF raw hashes, input-byte preservation and spaced/apostrophe/Unicode paths;
the checkout gate regenerates the synthetic teaching files after real policy
checkouts and has a missing-policy control. Normal Windows CI now runs this
same focused Python discovery; wait for that normal run rather than asking
the Linux-only Glen host to supply native Windows evidence or dispatching anew.
Require nonzero counts and report failures/skips with platform/Python/source.

Native Windows, large-reference preparation memory/performance, real frozen
annotation joins, mapping adequacy and biological interpretation remain unrun.
Do not reuse old exact 3′ values from v1 models as accepted v2 evidence. No Rust
changes/builds, GUI, live-model/provider calls, private-data analysis or release
acceptance follow. The owner authorized this merge, not a push or CI dispatch.
Existing closed Linux reader/MCP checks remain closed.

## Interpretation-Local Reference Reader Repair, 2026-11-05

Owner-authorized implementation follows the `dfca3261` diagnostic below.
`src/genomes/prepared_sequence_reader.rs` reuses only successfully parsed
FASTA indexes within one single-reference interpretation. Both primer
directions and all subject-window fallbacks share it. Per-read catalog
resolution, manifest/file checks and fresh FASTA reads are preserved through
the existing region reader; the public uncached and exact-only APIs remain.
No failed-load/sequence-window/global cache or annotation/alias reuse is added.

Raw index hashes bracket initial parsing; per-read metadata changes trigger
content rechecks. All reused indexes and bound manifests are rehashed before
returning the child report, including same-size/restored-mtime changes. Observed
drift is sticky and returns `primer_specificity_reference_changed_during_interpretation`
before a stale child can be persisted. Outer multi-reference admission and
post-interpretation validation are unchanged. These checks do not promise an
atomic filesystem snapshot, independently hash the whole FASTA on each hit or
authorize alternate references. Initial read failures keep their legacy
warnings/fallbacks; complete BLAST strings do not bypass the terminal fetch.

New synthetic Rust tests: nine `prepared_sequence_reader_` cases for one-parse
reuse, separate interpretations/resources, LF/CRLF and quoted/Unicode paths,
aliases/compatible-resolution refusal, ordinary missing/corrupt inputs, final
manifest/index drift and Windows handle closure; two `primer_specificity_local_reader_`
cases compare uncached/cached hit fields, raw/aggregated warnings, products and
summaries across both roles, full/partial/gapped/reversed HSPs and all three
alignment modes. A test-only uncached reader drives the existing path as the
comparison control, with literal alias/strand/terminal and no-truncation checks.
All fixtures are inline synthetic data; no external tools or private inputs.

The Linux commands below are now completed at `e4bc7386` as recorded above.
Retain them for reproduction or affected future source, not an unchanged repeat.
For unrun Windows checks use a source-matched isolated target and nonzero counts;
stop on primary failures. Do not repeat unrelated completed campaigns:

```sh
cargo check -q --locked --tests -j1
cargo check -q --locked --no-default-features --lib --bin gentle_cli -j1
cargo test --locked --lib -j1 prepared_sequence_reader_ -- --test-threads=1
cargo test --locked --lib -j1 primer_specificity_local_reader_ -- --test-threads=1
cargo test --locked --lib -j1 engine::primer_specificity_multi::tests:: -- --test-threads=1
cargo test --locked --lib -j1 primer_import_trace -- --test-threads=1
cargo test --locked --lib -j1 specificity_matrix_ -- --test-threads=1
```

The queued Linux MCP/parity check also passed at this source. No interim source
freeze or new bundle request follows.
Native Windows must execute LF/CRLF/path, timestamp restoration, rename/delete
and missing/corrupt cases; Unix passes alone are not Windows acceptance.

Historical executed Linux procedure, not a new replay request: after focused
checks passed, Glen repeated the authorized ACTB/TP53 retained-output
imports with the existing procedure: same host/dev/default profile and reference
identities, fresh absent state per attempt, one untraced and one traced run per
case, 60 s deadline/5 s grace, four launches/300 s total. Do not rerun BLAST,
reprepare references, relax hit/alignment/completeness policy or edit evidence.
Compare phases, wall/user/system/RSS and outcomes against `dfca3261`; timeout
is still incomplete, not a pass. If imports complete, retain/reopen reports and
compare scientific content/bindings between tracing modes, excluding only named
volatile fields. Completion or a faster phase is not a specificity/biology pass.
Return exact source/tree/lock/binary hashes, toolchain/platform/profile, commands,
counts/skips, failures and small manifest-bound evidence. Keep historical raw
trees unchanged. At authoring time Rust/native/performance acceptance was pending;
the Linux return above supersedes that status only within its named scope.
No local Codex Rust, push or dispatch occurred. No Claude review was invoked.

Authoring checks: 41 offline wiring/metadata/documentation checks passed, zero
skips, including source-negative controls and LF/CRLF policy checks. New Rust
files pass rustfmt; both changed parent modules parse under rustfmt without
rewriting unrelated formatting. Manual catalog/status and whitespace checks
pass; roadmap remains 400 lines. This is not Rust compilation or execution.

## Post-BLAST Phase-Trace Return, 2026-11-05

Glen completed the diagnostic at `dfca32611ad2b217382186038a4c3879d907fdde`,
tree `b5428910481eeaf23fe23ec0e04d885997fe62ac`, Linux x86_64 / Rust-Cargo
1.97.1, ordinary dev/default-feature CLI, one job, `CARGO_INCREMENTAL=0`.
The diff from `108cc5e5` contains only four documentation files. Retained
streams confirm 5 trace, 10 multi-reference and 13 matrix tests passed, plus
locked Cargo check. These results do not accept the later MCP-path repair.

The owner's 17,854-byte archive matches SHA-256
`978f85f6a18c19cc81475bfc9b68cc0aac882aea8af9ce8e09fc3569796926cc`.
Codex independently checked all 82 manifest entries and exact regular-file
coverage (83 files including the manifest). Manifest SHA-256:
`f87e592dd2cbb4ea798c453509a03cdc3730f42750b4d15a199c0e1577998d59`.
The [published report](https://github.com/smoe-bot/gentle_rs/blob/05055a36fb95528f9931fd1556498995c13005e0/docs/audits/glen_followups_20261008/post_blast_phase_trace_dfca3261/REPORT.md)
is byte-identical to the archive report. Historical evidence is unchanged.

All four ACTB/TP53 baseline/traced attempts ended at the prescribed deadline,
exit 124, with no stdout, saved state or surviving run-owned process. Each
trace contains 50 valid records, 27 begins and 23 completions, no truncation.
Both end with unfinished `forward_hits`: begin at 22.602560 s / 22.694039 s
from trace origin. Four earlier annotation decodes per run each process
465,769 records in 3.565-4.145 s; aliases take 5.212 / 5.260 s. Do not sum
nested spans or claim CPU stacks, completed hits, speedup or specificity.
The 37.397440 / 37.305961 s residuals in the supplied analyzer subtract trace
elapsed time from the external CLI deadline. Because those clocks start at
different boundaries, treat them as approximate deadline-budget residuals,
not independently established strict phase-duration lower bounds.

All 213 input-check lines and both reference-inspection receipts are identical
before/after. This verifies the retained receipts, not an independent rehash of
the absent large inputs or CLI binary. The diagnostic request is fulfilled;
real-import performance and full scientific acceptance remain open. Do not
repeat the unchanged four-run diagnostic or rerun BLAST/reference preparation.

**Proposal at the diagnostic return, now implemented separately above:**
`primer_specificity_hit_from_blast` unconditionally requests reference-backed
3-prime mismatch evidence before using complete aligned BLAST strings. The
prepared-window reader resolves/validates the resource and reloads the entire
FASTA index on each successful read. This is a source-level amplification
candidate inside `forward_hits`, not measured attribution to that subcall.
The proposal was an interpretation-local prepared reference/index reader, preserving
existing resolution, alias, coordinate/strand, warning/fallback and pre/post
reference checks. No global cache, hit truncation, mismatch-policy change or
skipped validation. Keep annotation-decode/alias reuse a separate later change.

Before a new real replay, require synthetic old/new result equivalence and
deterministic index-load-count checks, including partial/gapped/negative-strand
HSPs, missing/corrupt/stale references and full-alignment fallback. Glen/CI owns
Rust execution and the later bounded same-profile ACTB/TP53 comparison; no new
run, local Rust build, push or dispatch was performed for this review. The
plan has not been Claude-reviewed; a read-only second opinion is optional.
The documentation-only recording passes 19 documentation checks, manual
catalog/status and whitespace checks; roadmap remains 400 lines.

## RT Capture MCP Command-Path Repair, 2026-11-04

The reported `capability_registry_parity` failures are discovery drift, not a
missing engine operation or duplicate registry key. The registry already marks
`rna-reads primer-junctions` as MCP-prominent from `mcp-op`, but the generic
`op` tool omitted it from `tools/list.commandPaths`. The same omission affected
`primers assess-rt-candidates`, `primers compare-rt-candidates` and
`primers redesign-rt`. All four paths now appear on the existing tool; no
request, biological algorithm, permission, operation route or registry
classification changes. The parity/maturity projections remain unchanged.

The new offline source guard failed before the repair with exactly those four
paths and passes after it. A new Rust integration regression inspects the
actual tools/list projection, requires every glossary `mcp-op` path on the
single generic `op` tool and preserves its mutation/approval classification.
The two original failing parity assertions are unchanged.

The affected Linux discovery/freshness gate, including the new regression,
passed at `e4bc7386` in the return above. Keep these reproduction commands,
not a request to repeat completed comparator, gel, RT biology or Linux checks:

```sh
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo check -q --locked --tests -j1
```

Record actual tested source/toolchain/platform and nonzero counts. Local Rust
execution remains prohibited here; source-level passes are not a runtime pass.
Native, scientific and wider release acceptance remain separate.

Authoring checks on the working tree based at `dfca3261`: 36 offline metadata/
documentation checks passed, zero skips, including LF/CRLF policy checks.
Manual catalog/status, scoped regression formatting and whitespace checks
passed; roadmap remains 400 lines. The Rust integration test remains unrun.

## Read-End Metric Comparator Repair, 2026-11-04

The standalone Python utility is adapted from Glen / `smoe-bot`
`a1cf936d0cbccbbc6f1e2cbed8474006d740bc25`. It is not an engine operation,
`gentle_cli` command or MCP tool. Do not request a live-agent route that does
not exist, replay private FASTQ or claim transcript-aware UTR extraction.

Three review repairs: preserve CSV/TSV newline and exact-row semantics; separate
observed zero from blank censored/missing data (optional explicit
`--positive-metric` policy, no name-based inference); and reserve a fresh output
directory with closed/synced staging and atomic no-clobber publication, JSON
last. Numeric range validation prevents lossy integer inputs. Existing output
directories/symlinks and missing parents fail before output writes. Failed runs
can retain completed artifacts without a report; no automatic deletion/retry
or all-or-nothing/crash/hostile-parent-swap guarantee is claimed.

The offline CI step runs the synthetic process and publication tests:

```sh
python3 -m unittest discover -s tests -p test_compare_rna_read_end_metrics.py -v
```

Local verification on macOS arm64 / Python 3.15.0: all 16 focused tests and
86 offline release/metadata/documentation checks passed, zero skips. Rechecked
with the comparator working tree on main at `88cf72b6`, after resolving the
stash conflicts without reverting the equipment-pilot implementation. Manual
catalog/status and whitespace checks also passed. These counts do not certify
Rust, native Windows, upstream UTR extraction or biological interpretation.

Glen's focused native Windows check is closed at exact tested source
`108cc5e51d02b31f7af90d941821a2733ad50059`: Windows x64 / NTFS,
`Microsoft Windows NT 10.0.26100.0`, Python 3.15.2. The command above ran with
`python` rather than `python3`: 16 passed, zero failed/skipped, exit 0,
3.272 seconds. Both symlink cases executed. Spaced/apostrophe/Unicode paths,
LF/CRLF hashes, closed staging handles and rename/delete, existing/racing
destinations and fail-closed unsupported hard links are covered.

The successful [Windows run](https://github.com/smoe-bot/gentle_rs/actions/runs/41917134734)
uses wrapper commit `068ec42cdad5d6b9b65136bf58bea72e053c62d1`; its workflow
explicitly checks out and asserts the tested upstream SHA in `candidate/`.
The [retained receipts](https://github.com/smoe-bot/gentle_rs/tree/e22ca14b84dbef616678e8985781ff87ada8c97e/windows-evidence/108cc5e5)
confirm that identity, environment, command and all 16 test outcomes. These
were independently inspected here; no Rust or comparator rerun was performed.

Byte-accounting caveat: the retained Git logs use LF, while `SHA256.txt`
records their original CRLF Windows bytes. Both listed hashes match after
explicit LF-to-CRLF reconstruction, not the Git blobs as stored. Treat the Git
copies as declared line-ending derivatives, not byte-identical raw artifacts;
the original Actions artifact was not downloaded or independently revalidated.

Do not repeat completed TSS, CEL, primer/RT or sequence-design campaigns for
this standalone utility, or repeat this Windows check without an affected
change. This closes only the focused comparator gate: GUI, Rust, packages,
private-data accuracy and biological acceptance are not implied. Frozen
transcript models, CDS anchors, strand, clipping and cohort selection still
belong to the unimplemented upstream extraction; the input-table hash cannot
authenticate absent upstream provenance.

## RT Capture Integration Follow-Up, 2026-11-04

Integrate `fee72185` and `be6de439` from `gentle_rs_3` without inheriting their
worktree-local Rust passes as merged-source acceptance. The assessment,
comparison, junction-inspection and constrained-redesign contracts and earlier
verification limits remain in the [implementation ledger](rt_capture_primer_redesign_plan.md).
Main's microarray raw evidence and protein-length provenance repair are retained.

Merge-specific adaptation: all four new glossary JSON routes must obtain their
parser-only request from the corresponding discovery descriptor. Their strict
typed parsers reject the generic `{}` fallback. A focused Rust test exercises
the real smoke helper, parser and operation deserializer for each route, keeps
read-only admission and retains empty-object refusal. The existing full
glossary assertions are not weakened; the new offline guard is source evidence.

The four maturity projections were regenerated canonically in Python from the
merged glossary. Parity totals were reconciled from 965 distinct source/name
rows; row identities and ordering agree with merged declaration inputs. The
canonical Rust parity generator was not run. Glen/CI should execute affected
ledger scopes plus these integration gates with nonzero test counts:

```sh
cargo check -q --locked --tests
cargo test --locked --lib glossary_ -- --test-threads=1
cargo test --locked --test capability_registry_parity --test parity_matrix_freshness
cargo run --locked --bin gentle_examples_docs -- parity-matrix-check
```

Keep native/process, whole-reference performance and scientific acceptance
separate; no private data replay is authorized by this merge. At that merge
boundary the Q-Q utility repair had not started; its subsequent standalone
repair is recorded separately above. No local Rust, push or dispatch was
performed here.

Merge authoring checks: 37 offline metadata/documentation and targeted LF/CRLF
checks passed with zero skips. Manual/status, scoped regression formatting and
whitespace checks passed; roadmap remains 400 lines. These are not Rust/runtime
acceptance. Session-close reports four OK, two expected worktree/scope reminders
and no failures; unrelated `paper`, private data, prompts and `output/` are
excluded from the merge.

## Microarray Native Evidence Port, 2026-11-04

Glen / `smoe-bot` published the seven-file documentation commit
`37bc0ceab467b09050fe4e998477af8c460b9347` on
`glen/genome-transcript-protein-gui-20261103`. Port only that commit's
microarray changes, not the parent native overview or CD44/TSS tranche.
Both raw PNGs and `evidence.json` remain byte-identical to Glen's commit.
Retain main's reading card and protein-length provenance fix. Resolve the
backlog against current main, not by importing the parent's completed flags.

The receipt binds GUI source `503cb391` and CLI report source `3839426f`
separately. Both PNG hashes/dimensions and all four committed input hashes
match. The screenshots were visually inspected here; contrast and small-label
readability remain weak. They establish retained Linux/Xvfb observations at
their actual producer sources, not current-source, Windows/macOS or scientific
acceptance. Project/receipt hashes without retained local files are not claimed
independently verified. No genuine CEL processing was repeated.

Keep the high-contrast teaching-image task open. Do not edit/re-encode the raw
captures or rewrite their manifest. Any later crop/annotation needs a separate
file and derivation record; opacity/readability may need a new authorized
capture rather than a crop. Existing CD44/TSS source-review findings remain
separate from this microarray-only import.

Offline regressions bind canonical source/catalog graphics, both producer
identities, raw image/input hashes and the selected 121-bp beta-row coordinates.
Real LF/CRLF Git checkouts retain the pixels and text bytes; negative controls
remove the scoped LF rules and global fallback to expose conversion. No Rust
code, generated chapter, capture, provider call, push or CI dispatch is added.
Normal CI still owns the canonical Rust catalog/manifest freshness checks;
archive verification does not replace those checks or new native acceptance.

Authoring checks: 121 offline tests selected, 118 passed and three real-binary
cases skipped; the visual-only gate passed with 70 tutorial entries and 92
declared visuals. Manual/status and whitespace checks passed; roadmap remains
400 lines. Session-close: four OK, two intentional worktree/scope reminders,
zero failures. Unrelated `paper`, private data, prompts and `output/` remain
untouched.

## Tutorial Review Follow-Up, 2026-11-04

The owner's forwarded overnight report has SHA-256
`2efc2501a1b2143f849a187bf0de3c25f5d95eeecdcb0d17af0a81d32e584cea`.
It distinguishes upstream `069bd37b` from Glen's local `1e75d926`. A bounded
fetch at the time returned `89eec1b9`, whose foundation was already integrated;
the newer source/captures were then unavailable. The publication and bounded
microarray port above supersede that availability limitation. They do not fix
or independently verify the parent overview/CD44/TSS tranche.

The existing main table has the same expected-to-derived length fallback.
Repair only its presentation: derived length requires `comparison.derived`;
external expectations and unattributed legacy expectations remain separate.
No translation, coordinate, comparison-status, engine/schema or SVG algorithm
changes are made. Four synthetic regressions cover external-only records,
conflicting expected lengths, 742/139 label controls, zero versus unavailable,
missing source/length and semantic table painting without model changes.
These length controls are not new CD44 derivation evidence.

Tutorial 06.08 now distinguishes the native tables from the exported SVG and
teaches unavailable-versus-expected data, reverse-strand reading and separate
project/file persistence. The 07.04 reading card keeps authentic coordinates
apart from synthetic effects; 08.02's backbone link is checkout-relative.
All are hand-written canonical guides, not generated chapters. That earlier
table/prose repair kept catalog, manifest, historical review dates, raw images
and evidence manifests unchanged. The [image backlog](tutorial/coverage_backlog.md#overnight-review-follow-up-2026-11-04)
records the remaining overview/color mapping and readable-image work.

Run only the changed scopes before ordinary CI; require four regression matches:

```sh
cargo test --locked --lib isoform_protein_length_ -- --test-threads=1
cargo check -q --locked --tests
```

Native follow-up: inspect an external-only row and both genuinely derived
rows, horizontal scrolling and warnings. Reopen the same project and verify
that the display did not change its evidence. Keep Windows/native acceptance
separate from the semantic paint test. Do not repeat completed CD44 source
fetches, TSS DuckDB or array/CEL campaigns merely for these wording changes.
The existing offline source checks reproduce the old link/guide/fallback gaps;
no local Rust build, real-process replay, new capture, push or dispatch ran.

Authoring verification: 119 offline tutorial/checkout/screenshot-policy tests
were selected, 116 passed and three binary-dependent cases skipped. LF/CRLF,
manual/status, scoped rustfmt and whitespace checks passed. Session-close:
four OK, two expected worktree/scope reminders and no failures. Rust and native
acceptance remain pending; unrelated `paper`, private data and untracked files
are untouched.

## Primer Preflight Introspection Repair, 2026-11-04

The owner's local library run reported 4,838 passed, two failed and 16 ignored
(exact source/toolchain not supplied). Both failures name only the new MCP
`primers_preflight` row as `registry_only`; the shell-spelled command already
had its descriptor. Add the MCP descriptor using the existing shared read-only
inspection helper, its exact optional MCP argument names and matching empty
project-fact preconditions. Keep both whole-registry tests unchanged. This
does not probe Primer3 during discovery, change approvals or assert that a
backend is installed. The separate JS `u128` transport issue remains open.

Run only the affected Rust scopes before the ordinary suite, with nonzero counts:

```sh
cargo test --locked --lib execute_introspect_primer_preflight_mcp_contract_preserves_read_only_semantics -- --test-threads=1
cargo test --locked --lib execute_introspect_primers_preflight_readiness_is_not_tool_reachability -- --test-threads=1
cargo test --locked --lib execute_introspect_capabilities_projects_full_registry_with_fact_annotations -- --test-threads=1
cargo test --locked --lib execute_introspect_all_capability_section_has_full_fact_annotations -- --test-threads=1
```

The new regression compares both introspection routes with the shell contract,
checks the MCP registry origin/argument names and unchanged engine state. The
offline source guard fails before adding the descriptor; it is not Rust/runtime
acceptance. No completed scientific, Lua/MCP execution or publication campaign
is reopened by this metadata repair. No local Rust build, push or dispatch.
Authoring verification: 46 offline metadata/GUI/documentation checks passed,
including LF/CRLF documentation checkouts; manual/status and whitespace checks
passed. These counts do not establish the pending Rust regression results.

The attempted report transfer has no verified readable source path yet. A
bounded read-only SSH search through the owner's `gentle` alias received
`Permission denied` for `/home/clawbio/.openclaw/workspace/artifacts`; no sudo,
access-policy bypass or large download was attempted. Ask Glen for the small
report plus its manifest at an explicitly transfer-readable path/URL, rather
than treating his host-local working path as an established transfer path.

## Glen's Scoped Return At 89df315a, 2026-11-03

Owner-forwarded report for exact source
`89df315a80e0d3ae02e5685efd7a94e3895fbaa7`: dependency-pin guards 2/2,
panel-copy GUI regressions 3/3, annotation/SAM publication 40/40, Lua bindings
2/2, MCP 4/4, shell 3/3 and registry/matrix parity 2/2 passed. Default,
no-default and script-interface locked checks passed; the lockfile remained
byte-identical. These are Glen-reported results, not local reruns or independent
raw-log verification. Retain their scope; do not reissue the completed campaign
merely because older sections request it. Native GUI/Windows behavior is not
established by the counts alone. CI was still in progress when Glen reported;
no terminal full-platform verdict is inferred here.

The independent release-safety failure is the offline egui family parser, which
did not accept `=0.37.0`. It now accepts an optional leading `=` and preserves
all family comparisons. Fourteen compatibility/pin checks pass locally,
including a fail-before/pass-after parser regression. Do not loosen the pins.

Both JS binding tests failed at `load_project` with `u128 is not supported`,
before reaching their handoff assertions. Source inspection confirms that the
current `serde_v8` serializer/deserializer do not implement the 128-bit integer
methods used by populated typed projects. This is not a provider, permission,
compiler or Lua failure. The JS production transport has not been repaired in
this guard-only change; the original failing tests and populated fixture remain.

The next JS slice must define a lossless adapter-only transport repair before
editing it: preserve native byte buffers, engine JSON/schema/hash contracts and
read-only file/state boundaries; never round timestamps through a JS Number,
cast away high bits or remove lineage to make the fixture pass. Reuse the two
failing `js_protein_expression_handoff_` VM tests and add populated-project
load/inspect/save controls, buffer checks and integer-boundary/refusal cases.
Keep the options contract's existing unsafe-timestamp refusals. A generic codec
or dependency patch is not silently authorized by fixing the Python guard.
Offer an optional read-only Claude critique before finalizing that larger plan.

After the owner's push, rerun the repaired offline compatibility gate and inspect
normal CI. Repeat only the affected JS/transport scopes after a real JS repair;
retain the already passing Lua/MCP/shell and core scientific evidence separately.

## Flagged Dependency Pins, 2026-11-03

Historical exact-pin repair. For bio only, the 2026-11-08
[compatible-update policy](#compatible-bio-update-policy-2026-11-08) supersedes
this section's exact requirement; other flagged pins remain unchanged.

The owner authorized a separate manifest-only follow-up for Glen's eight
flagged constraints. Pin the versions already present in the lock:
bio `=4.0.1`; egui, eframe and egui_extras `=0.37.0`; mlua `=0.12.2`;
rustyline `=19.0.0`; lopdf `=0.46.0`; scraper `=0.29.0`. No feature, optional
flag, default-feature setting, workspace inheritance or other requirement is
changed. This does not pin the entire transitive graph or other wildcard
dependencies outside the finding.

`Cargo.lock` remains byte-identical, with SHA-256
`7ad5b2dc9485675a75d272fa945c75342dceb60ef582f2885015f3e14df41e76` before and
after the change. The two offline `DependencyPinTests` in
`scripts/test_release_candidate.py` check exact locked requirements and reject
wildcard/caret/unlocked controls; the actual-manifest guard failed for all eight
old constraints under both LF and CRLF before repair. CI already executes this
Python module. Semantic TOML checks are not Cargo resolution or Rust acceptance.

Local authoring verification: `python3 -m unittest scripts.test_release_candidate`
passes all 51 tests with no skips. A structured before/after TOML comparison
confirms only the eight version requirements changed; raw lockfile comparison
and the SHA-256 above confirm unchanged bytes. Documentation/maturity and
whitespace checks pass. Session-close reports four OK, two reminders (intentional
worktree changes and manual scope fidelity), zero failures.

Glen/CI should run, without changing the lockfile:

```sh
python3 -m unittest scripts.test_release_candidate.DependencyPinTests
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --tests
cargo check -q --locked --features script-interfaces --tests
git diff --exit-code -- Cargo.lock
```

No local Rust build/check/test, Cargo resolution/update, network operation,
push or CI dispatch was performed for this repair. Preserve unrelated local
files and the paused goal; no compiler-floor or scientific verdict changes.

## Native Panel Specificity Review, 2026-11-03

PCR Designer's **Review specificity** reads only the currently displayed saved
panel. Compare it to the unchanged **Copy specificity matrix command** on the
same project. The GUI has no separate matrix model, search or readiness update.
The worker uses the read-only snapshot operation; a busy frame defers delivery.
The GUI host disconnects a pending reply when the specialist closes even if
the parent DNA view remains alive.

Authoring began at `c6adcb83`, one workflow commit beyond the supplied
`e7778466`. Local main then advanced to
`89df315a80e0d3ae02e5685efd7a94e3895fbaa7` (dependency pins). Its status
entries are retained. Offline checks on macOS arm64 at that base plus the
viewer patch: eight specificity tutorial/source tests pass, the manual catalog
check passes, and the three requested GUI files pass scoped rustfmt. The
parent DNA check retains two pre-existing child-module hunks; the extra PCR
dialog close hook introduces none of the six existing `app.rs` formatting
hunks. Whitespace is clean; the roadmap stays at 400 lines. No Rust build or
test ran. After the owner pushes, record the actual SHA, platform, toolchain
and test counts; **zero matches are a failure**:

```sh
cargo check -q --locked
cargo test --locked --lib transcript_assay_specificity_review_ -- --test-threads=1
cargo test --locked --lib specificity_matrix_ -- --test-threads=1
```

The first filter has ten authored cases covering the renderer, both paging
axes, historical rows/citations, click admission, exact ID and state, engine
errors, contention, stale replies, post-acceptance edits, reload and close.
The second includes the existing engine cases and the three copy-action tests
from `da97696d`. Neither authoring slice ran those tests; Glen subsequently
reported their 3/3 pass at `89df315a` as retained above. That historical result
does not accept the new viewer; do not reissue it as an unfinished old campaign.

Record native GUI acceptance separately:

1. Open a saved panel with genomic, transcriptome and superseded receipts.
   Click **Review specificity**; check responsiveness and compare both grids,
   receipt details, separate citations and nonclaims with the copied command's
   JSON. Check that an unverified recorded status is not a verified current verdict.
2. While loading, separately switch panels, redesign under the same ID, undo
   and reopen another project. No stale result should appear. Change project
   state after acceptance and check the request to review again.
3. Page a large panel to its last assay and receipt in each dimension. Close
   PCR Designer mid-load, leave the DNA window alive, then reopen the specialist.
4. Check native Windows separately. A macOS run is not Windows acceptance.
   Use synthetic/disposable saved evidence: no BLAST, provider contact or private data.

Known limit: a panel already replaced externally **before** the click can
retain the same ID as the displayed report. The matrix header carries the
reader's current design digest; the GUI does not independently hash the cached
panel or claim to detect that prior replacement. No scientific verdict follows.

## Panel Specificity Review Shortcut, 2026-11-03

The transcript-panel result now offers `Copy specificity matrix command`, bound
to the displayed report ID through the shared shell's quoting helper. This is
clipboard-only: no searches, project reads/writes, reference selection or
readiness attachment. Missing evidence remains unassessed. The existing engine
matrix/CLI/MCP contracts and conditional multi-reference fan-out remain unchanged.
Authoring checks: 26 offline tutorial/source/documentation checks passed,
including LF/CRLF documentation checkouts; scoped rustfmt, manual/status and
whitespace checks passed. The three new Rust tests have not been executed.

After the owner's push, run the three new tests with nonzero counts:

```sh
cargo test --locked --lib specificity_matrix_copy_ -- --test-threads=1
```

Then natively open a saved transcript-panel result, copy/paste into `Shell...`
and run the printed command. Confirm the ID identifies the displayed report,
copying alone does not execute, and no-receipt results remain unassessed. Include
quoted IDs and native Windows clipboard/paste separately. Do not repeat completed
Primer3 or whole-reference campaigns. Offline/source checks are not Rust/native
acceptance; no local Rust build, push or dispatch is performed for this slice.

## Publication Fault-Test Port, 2026-11-03

The owner requested the two pending test fixes from
`codex/glen-six-gate-recheck-20261102` in worktree `1ee2`. Port only the changes
in `annotation_prepare/tests.rs` and `observed_sam/tests.rs`; the source
worktree's uncommitted drafts remain untouched. The already integrated
`3e1b7bda` glossary repair is not reapplied. Preserve the branch's useful
[historical audit evidence](glen_audit_followups_20261008.md#preserved-six-gate-evidence-2026-11-03)
without reviving its cleared goal or changing `.12` publication status.

The injected artifact-directory replacement may be refused on Windows while
identity handles remain open. Accept only Windows `PermissionDenied` at that
injection, then require the synthetic fault, no moved directory and complete
owned-output cleanup. Successful replacement retains the original orphan and
preserved-file assertions. Two additional tests inject an unowned file without
renaming, then require that file alone to survive cleanup, with no completion
report and the correct orphan diagnostic. Production publication and validation
rules are unchanged; no blanket platform skip or weakened fixture hashes.

Local validation passed on macOS arm64, Rust/Cargo `1.100.0-beta.7`, default
features plus `script-interfaces`, against base
`fa5de5530d8e34774d53d3aa89b112a8b78fec5f` plus this two-test-file port. The raw
`git diff --` for those two files has SHA-256
`9b121cf111a2707cc68c33620f1c87c73da8419567eefc2e1324a58d01e658cd`;
the unchanged `Cargo.lock` SHA-256 is
`7ad5b2dc9485675a75d272fa945c75342dceb60ef582f2885015f3e14df41e76`.
Formatting-normalized comparison matches both original drafts exactly.

Use `CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0`
for the test commands below; the check used the first two settings. All Cargo
commands used `--locked --offline --features script-interfaces` and `-j 2`:

| Scoped command / check | Result |
| --- | --- |
| `cargo test --lib annotation_prepare::tests:: -- --test-threads=1` | 28 passed, zero failures/ignored, including six adapter tests |
| `cargo test --lib observed_sam::tests:: -- --test-threads=1` | 12 passed, zero failures/ignored |
| Initial `cargo test --lib publication_preserves -- --test-threads=1` | 3 passed; overlaps the module runs, not three additional port regressions |
| `cargo check -q --tests` | Exit 0 |
| `scripts.test_frozen_024ed6_audit` and the 024ed6c0, 6d8b4db7 and non-publishing container cases in `scripts.test_tutorial_checkouts` | 8 passed, zero skips, including LF/CRLF and missing-policy controls |
| Scoped rustfmt, `scripts/manual_docs.py`, `git diff --check` | Passed |

The test build retains the existing macOS linker warning about `__eh_frame`
exceeding the compact-unwind table limit; compilation and execution succeed.
No unfiltered workspace/full-library suite, native Windows/Linux, package or
real-provider test was run for this test/documentation-only port.

Run the affected annotation-preparation and observed-SAM test modules before
ordinary CI. Native Windows must exercise both replacement/refusal and unowned
file preservation on the committed source; macOS cannot establish Windows
file-handle behavior. Full-suite, package and provider/scientific acceptance
remain separate. Uploads and CI dispatch remain with the owner.

## Protein Expression Named Scripting Bindings, 2026-11-03

Authored from `fa5de553`, preserving its preflight readers: named JS/Lua
`protein_expression_handoff(project, options)` returns the existing v2 report.
Shared inline admission also serves MCP; its saved-state/error boundary stays
intact. No new operation, provider action, write, agent-plan permission or
project codec is introduced. Lua nested typed decode preserves reviewed empty
lists/false. JS options, not the project, travel as JSON; unsafe timestamps are
refused. Per-interface glossary signatures and tutorial 09.03 describe these
first-class bindings, not a previously impossible execution route.

Rust checks/tests and JS/Lua VM execution have **not** run here. The owner's
exception covers only canonical generator builds/runs; their successful output
is not runtime or native acceptance. Current CI scripting jobs compile optional
features and smoke-test binaries but do not execute these new tests. Run each
adapter filter explicitly and require at least two tests, not a zero-test pass:

```sh
cargo test --locked --features js-interface --lib js_interface::tests::js_protein_expression_handoff_ -- --test-threads=1
cargo test --locked --features lua-interface --lib lua_interface::tests::lua_protein_expression_handoff_ -- --test-threads=1
cargo check -q --locked --features script-interfaces --tests
cargo test --locked --no-default-features --lib mcp_server::protein_expression:: -- --test-threads=1
cargo test --locked --no-default-features --lib engine_shell::protein_expression_handoff:: -- --test-threads=1
cargo test --locked --test capability_registry_parity protein_expression_handoff_named_scripting_bindings_preserve_canonical_usage
cargo test --locked --test parity_matrix_freshness
```

Use the same effective profile/provider catalog for the in-process oracle and
adapters. Compare every report field except `generated_at_unix_ms`; projected
provider status depends on that catalog, not a universal GeneArt guarantee.
The existing MCP positive control assumes the bundled catalog; run it with an
isolated configuration. Check saved and in-memory project preservation, saved
versus explicit-empty objectives, JSON/text, missing product, partial/withheld
requirements, unknown/null/path options and timestamp refusals. Read the named
signature back from capabilities and call the actual registered VM function.

Rehearse the one-line Lua/JS examples in tutorial 09.03 using source-built
REPLs; no `--script` flag exists. Retain manual/hybrid status and the four
historical macOS GUI captures at `9f01098b`; they prove no new scripting path.
Native/platform/live-agent/provider/biological acceptance remains separate.
No push or CI dispatch is authorized by this handoff.

Authoring verification at `fa5de553` plus this uncommitted diff on macOS:
127 offline tests passed, with seven expected binary-dependent skips, across
`scripts.test_specialist_tutorials`, `scripts.test_release_workflow_tutorials`,
`scripts.test_capability_maturity` and `scripts.test_tutorial_checkouts`.
The canonical generator-only build/run below succeeded; its parity and
tutorial-catalog check modes, maturity/manual/documentation checks and
whitespace checks pass. This build has the existing headless dead-code and
macOS linker warnings; it does not compile or execute optional scripting tests.
Whole-file JS/Lua rustfmt still reports pre-existing SAM-block formatting;
new hunks are formatted and those unrelated blocks were preserved.
Session-close reports four OK, two reminders (intentional dirty/generated
files and manual plan fidelity), zero failures. The old goal and unrelated
`outputs/` remain untouched; manifest/lock, agent prompt and screenshots are
unchanged.

```sh
cargo run --locked --no-default-features -j1 --bin gentle_examples_docs -- parity-matrix-generate
target/debug/gentle_examples_docs tutorial-catalog-generate
python3 scripts/generate_capability_maturity.py --write
python3 scripts/manual_docs.py --write
```

Changed code is in `src/engine_shell/protein_expression_handoff.rs`, its module
declaration, `src/{js_interface,lua_interface}.rs`,
`src/mcp_server/protein_expression.rs`, `src/test_support{.rs,/protein_expression.rs}`
and `tests/capability_registry_parity.rs`. Teaching checks are in
`scripts/test_specialist_tutorials.py`. Documentation covers README,
`docs/cli/{javascript,lua,planning}.md`, protocol/agent docs, changelog, roadmap,
this handoff and the status-prefixed original review prompt. Tutorial 09.03's
Markdown/source notes and generated catalog changed, along with the glossary,
parity matrix and glossary hashes in the four existing maturity catalogs.

## egui 0.37 And Ported Glossary Repair, 2026-11-02

The owner's updated root manifest/lock select egui/eframe 0.37, while the
vendored Markdown adapters still selected 0.36. The attached compiler output
also names removed screenshot events/UserData and private scroll offsets.
Align the two patched manifests/lock graph; use screenshot callbacks and public
scroll APIs. Preserve the owner's other dependency choices, capture approval,
exact viewport identity, one-shot completion, preview-before-transport and the
disabled CLI capture bridge. No actual screenshot, model call or local Rust
execution occurred; eleven offline API/manifest guards and source-matched Rust
regressions are authoring evidence only. The warning follow-up renames 186
`Id::new` calls to their `Id::unique` alias and four UI `id()` calls to
`scope_id()`, including the patched Markdown sources. A mechanical diff check
preserves every original argument/salt; no transient/positional ID or blanket
warning suppression is introduced.

The owner's test build after `1e46525c` exposed seven further compiler errors,
not just warnings: four wheel events lacked `source`, and three test Areas used
removed string-to-Id conversions. Use `MouseWheelSource::Unknown` for synthetic
events and explicit named `Id::unique` for those Areas. The two feature-gated
promoter events and GUI benchmark have the same wheel repair. Two source guards
failed before these edits and pass afterward; that is not a Rust compile verdict.
The later warning-only repair also removes the now-unused stored capture
transport viewport, retaining request identity, logical target/stale checks and
the actual approval-bound routing. A new `agent_screenshot_` regression checks
both native and embedded Help-window targets without completing any capture.
The owner's later full-library run reports 4816 passed, 1 failed and 16 ignored
(source/build identity not supplied). The sole failure is this new routing
fixture: the native child was listed in `RawInput` but not registered as live,
so egui removed it and its queued screenshot command at `end_pass`. The native
case now uses `show_viewport_deferred` before approval and asserts its retained
class/callback as well as the original one-command/target checks. Production
routing and consent are unchanged. The owner's subsequent compile rejected
`assert_eq!` on `ViewportClass` because the enum does not implement `Debug`;
the same class assertion now uses `matches!`, with an offline regression guard.
First rerun only this filter before another
broad suite; the supplied full-suite result remains failed, not accepted:

```sh
cargo test --locked --lib agent_screenshot_targets_native_or_embedded_viewport -- --test-threads=1
```

An earlier owner compile exposed three unqualified `ViewportId` references in
these regressions; they now use the test module's existing `egui::ViewportId`
path, with an offline guard. The native/embedded and one-shot assertions remain.
Do not infer warning-free compilation from source inspection: Rust/native
verification of this update remains pending.

Codex's warning-follow-up authoring checks: 91 offline passes, one explicit
binary-dependent skip, plus Rust syntax parsing (not type checking), manual/
documentation checks and clean whitespace. Packaging-policy tests use synthetic
or mocked executables; there is no new Rust/package/native verdict.

The owner's later local full-library run reports 4814 passed, 1 failed and
16 ignored (source/build identity not supplied). Its single failure reads a
pending edge-scroll request under a key derived from `ScrollAreaOutput.id`,
although egui 0.37 pre-hashes that scroll area's salt. The repair shares the
actual existing UI-scoped pending key between production and the test frame;
it does not remove the real queued-request or wheel/drag refusal assertions.
A new idle-frame control verifies application/consumption of a pending offset.
First rerun the two focused cases below, then the existing `tutorial_menu_`
filter before another broad suite. No new Rust pass is claimed here:

```sh
cargo test --locked --lib tutorial_menu_pending_edge_scroll_ -- --test-threads=1
```

Glen's `7f73a42a` glossary failure is fixed by the port of `3e1b7bda`.
`src/engine_shell/tests.rs` matches the already tested blob
`f9830fcf57aa578efd2ded8bc6104323b9c241c1`: both platforms pass glossary 12/12
at `43670b25` in [run 41375688160](https://github.com/smoe/gentle_rs/actions/runs/41375688160).
That run is terminal/failure overall: macOS passes (root 4818/0/17), Windows
fails two separate directory-replacement tests (root 4754/2/17); headless Linux
and safety pass, full Linux/optional audits skip. Do not repeat the old glossary
diagnosis or present it as a pass of the newer egui dependency graph.

After the owner commits/pushes these edits, use an isolated target/cache and
the actual published source/lock. Run only the changed scopes before ordinary CI:

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_mcp
cargo check -q --locked --all-features --tests
cargo check -q --locked --bin gentle_egui_window_repro --bin gentle_dna_viewer_repro
cargo check -q --locked -p gentle-benchmarks --bench gui_operations
cargo test --locked --lib engine_shell::tests::glossary_ -- --test-threads=1
cargo test --locked --lib egui_capture_ -- --test-threads=1
cargo test --locked --lib egui_compat::tests:: -- --test-threads=1
cargo test --locked --lib agent_screenshot -- --test-threads=1
cargo test --locked --lib tutorial_menu_ -- --test-threads=1
cargo test --locked --lib help_docs -- --test-threads=1
cargo test --locked --lib pcr_designer_specialist_polls_and_applies_completed_worker -- --test-threads=1
cargo test --locked --lib splicing_locus_interactive_inspector_renders_score_sequence_and_saved_region -- --test-threads=1
cargo test --locked --lib rna_read_support_and_preview_tables_reserve_practical_minimum_heights -- --test-threads=1
cargo test --locked --lib --features gui-test-support promoter_pair_button_is_reachable_without_horizontal_scrolling -- --test-threads=1
cargo test --locked --lib --features gui-test-support promoter_action_rows_are_readable_without_horizontal_scrolling -- --test-threads=1
python3 -m unittest scripts.test_egui_compatibility -v
```

Require nonzero counts. In native checks verify the selected window supplies
the image, decline/unsupported/stale requests do not capture, duplicate
completion produces no second attachment, and preview requires a later explicit
send. Include the new native-versus-embedded transport regression; preserve
window focus/drag state, Markdown scroll keys and the manual help-request ID.
Check last-row tutorial access, reversed edge hovering, wheel/drag priority,
Markdown links/images and lineage/rack scroll positions. Native captures require
the existing user consent; no live provider or headless capture is requested.
Report source/toolchain/platform/features, binary hashes, results and skips;
do not rerun completed gel/sequence-design campaigns or rewrite historical evidence.
The six-follow-up goal was cleared; this repair does not resume it or authorize
new pushes/dispatches. The separate publication-test drafts are now ported under
the [later owner request](#publication-fault-test-port-2026-11-03), not as part
of the egui update; their original worktree remains untouched.

## Discovery And Decision Status Follow-Up, 2026-11-02

Glen's [retained report](GENtle_discovery_verification_7f73a42a.md) confirms the
scoped discovery/decision repair at `7f73a42ab97d0aa3c11c37182da9e7b7c4b60c6a`:
canonicalization 4/4, capability parity 23/23, matrix freshness 1/1,
`parity-matrix-check` OK and metadata/status 3/3. The 952 capability rows and
164 gap keys are unique; 72/72 decisions are classified (70 active, 2 deferred).
Platform: Debian forky x86_64, Rust/Cargo 1.97.1, locked dev/test, headless root
gates. Preserve this positive verdict at that source/lock; it does not accept the
later egui update. Do not repeat this completed campaign unless its contracts
change. Native GUI/Windows Slice 6 checks remain separate and open.

The broader `glossary_` gate failed separately (12 passed, 1 failed; isolated
rerun 0/1) because the SAM inspector fixture supplied an untyped binding.
Glen's later scan locates its introduction in `afdbf8e4`; parent-byte equality
alone was too narrow to establish that provenance. The supplied report is
retained unchanged; this addendum records his clarification. The fixture repair
and its newer-source verification boundary are described above. CI at `7f73a42a`
remains failed, not retroactively green.

Supplied report SHA-256:
`8bf11a37d5ce7fe0b14d33dd6c5e3d57b254e9e2341d7b4d7e038f76bda7c307`.
Its trailing blank line is intentionally retained with the original bytes;
the corresponding whitespace-check warning is not a source-code defect.

The four redundant/conflicting introspection entries are removed. Readiness
retains an optional capability and `--arg`, `--seq-id`, `--readiness`,
`--evidence`, `--ui-host`. Registry identity is `(source, name/path)`, with
unioned routing metadata and primary catalog prose/order. Different JS/Lua
signatures survive in labelled usage and `x-gentle-usage-by-interface`;
same-interface usage/operation conflicts fail deterministically. No parser,
execution route or approval policy is changed.

Codex migrated the committed matrix by removing byte-identical duplicate rows
and deriving totals from unique rows; Glen subsequently verified its canonical
Rust generator and freshness at the reported SHA.
There are 164 distinct gap keys: MCP 42, JS 61, Lua 61. Lower counts do not
close adapter gaps. All 72 decisions now have explicit status: 70 active and
two deferred, with no new scientific approval or minimum Rust version.
Four maturity catalogs change only their glossary-source hash; evidence and
unassessed verdicts are unchanged.

Codex authoring checks: 37 targeted Python tests pass with zero skips, including
maturity/discovery metadata, all decision statuses, LF/CRLF matrix checkout,
CLI glossary coverage and existing parity-count/override checks. Maturity,
manual and documentation freshness, scoped rustfmt and whitespace pass;
`Cargo.lock` is unchanged. Session-close: four OK, two expected scope/worktree
warnings, no failures. These checks do not execute the registry or generator.
Unrelated `paper`, private data, `output/` and untracked prompts are untouched.

Glen executed the following commands with an isolated, run-owned target. These
are retained reproduction commands, not an outstanding rerun request:

```sh
cargo test --locked --offline -p gentle-protocol --lib glossary_canonicalization_tests -- --test-threads=1
cargo test --locked --offline --no-default-features -j1 --test capability_registry_parity --test parity_matrix_freshness -- --test-threads=1
cargo test --locked --offline --no-default-features -j1 --lib glossary_ -- --test-threads=1
cargo run --locked --offline --no-default-features --bin gentle_examples_docs -- parity-matrix-check
```

The report retains the actual source/tree, lockfile and binary hashes, commands,
counts and limits. Optional readiness, per-interface JS/Lua usages and
deterministic same-interface conflict refusals executed; neither remaining
adapter gaps nor native/platform/release acceptance is thereby closed.
Normal CI remains the platform gate. Codex did not rerun Glen's Rust checks,
push or dispatch; this evidence update does not resume the cleared goal.

## Headless Parity Test Gate, 2026-11-01

Codex executed the owner's narrow local-build exception after rebasing onto
`d9b1015cf27c7dfc2c8b73a675c11c232ded3bf4`. This repair's only Rust change is
`tests/capability_registry_parity.rs`, Git blob
`8ad3efa35252e068e4dbd509201fd9fc73bb78de` (SHA-256
`48dee34e494ee3848946549511b24b2617e3e781de4983ee1b17464ec01b7cbf`).
Platform: macOS 27.2, arm64; Rust `1.100.0-beta.7 (fd8e211ae 2026-10-17)`,
Cargo `1.100.0-beta.7 (72b77b8d9 2026-10-17)`, ordinary dev/test profile.

After these checks, rebase onto `3002836d18e3c65568b164d4c185b85998f9662b`.
The test blob is unchanged; that newer main's SAM implementation and Slice 6
harness changes were not rerun here. The results below stay bound to the
identified `d9b1015c` source plus the test patch, not the post-rebase integration.

The GUI-only import and two GUI-oracle assertion groups are now gated on
`desktop-gui`. Assertions/messages are preserved; production modules, capability
registry, GUI oracle, CI, fixtures and scientific behavior are unchanged.
Changed files: the test above, `CHANGELOG.md`, `roadmap.md`,
`operation_lifecycle_slice5_plan.md` and this handoff. The roadmap parks the
missing CI headless-test-target gate without changing CI/release-policy checks.

| Command | Result | Passed / Failed / Ignored |
| --- | --- | --- |
| `cargo test --locked -j1 --test capability_registry_parity` (before edit) | pass on retry | 21 / 0 / 0 |
| `cargo check -q --locked --no-default-features --tests` | pass | not a test run |
| `cargo check -q --locked --tests` | pass | not a test run |
| `cargo test --locked --no-default-features -j1 --test capability_registry_parity` | pass | 21 / 0 / 0 |
| `cargo test --locked -j1 --test capability_registry_parity` (after edit) | pass | 23 / 0 / 0 |
| `rustfmt --edition 2024 --check tests/capability_registry_parity.rs` | pass | not a test run |
| `git diff --exit-code -- Cargo.lock` | pass; unchanged | not a test run |

All test runs have zero filtered/measured tests. The initial baseline attempt
executed no tests: `apple-metal`'s Swift build was denied host-cache writes by
the sandbox. The same command passed with the required cache access. Existing
headless dead-code and macOS `__eh_frame` linker warnings remain unmodified.
Temporary local logs are under `/tmp/gentle-headless-parity-20261101/`; they
are not committed release/acceptance artifacts. Lockfile SHA-256 remains
`a393a1f3803b374bbaf5ef4f3e23813e974b1676f61a669dacb24063c059fb6a`.

Before/after lists: these same 21 names pass in the baseline and in both
post-edit configurations. Sorted name/outcome comparisons are identical after
excluding the two additional desktop-only tests below. All failure sets are
empty, so no pre-existing failure was repaired or hidden.

```text
arrays_parser_subcommands_are_glossary_paths_or_explicitly_allowlisted
cli_capabilities_reads_protocol_registry
collection_lift_policies_attach_only_to_existing_capabilities
engine_capabilities_project_protocol_registry
engine_operations_are_reachable_on_shell_route_adapters
gap_rows_are_listed_in_open_gap_section
gene_set_resolve_gui_affordance_is_prominent
glossary_commands_have_matching_registry_entries
logical_gene_set_is_not_silently_treated_as_a_physical_pool
mcp_tools_list_has_no_unregistered_prominent_tools
not_applicable_cells_have_curated_override_entries
parity_overrides_do_not_demote_prominent_or_shell_routes
physical_container_pool_actions_are_declared_as_context_agnostic_combines
prominent_projections_match_user_facing_listings
registry_for_adapter_includes_every_reachable_projection
registry_surfacing_covers_all_parity_adapters
sequence_design_registry_family_is_external_without_reclassifying_controls
serial_arrangement_gel_actions_are_declared_as_context_agnostic_arranges
shell_alias_registry_mutation_tags_match_contract
shell_alias_registry_resolves_to_known_capabilities_and_parity_policy
transcript_assay_panel_reachability_is_locked_across_adapters
```

The two additional tests both pass with desktop features and are not compiled
headlessly:

- `gene_set_resolve_gui_affordance_has_menu_palette_entry`: the original
  gene-set test's final menu/palette assertion, unchanged.
- `gui_prominent_projections_match_menu_palette_oracle`: the original prominent
  projection test's forward GUI membership assertion and unchanged reverse loop.

Keep the existing lifecycle `--no-default-features --tests` command lists.
This result supersedes their unrun compile status only at the identified source;
it does not execute their other Rust/VM/process suites or close native checks.
Windows, Linux and any CI run remain unverified here. Reuse these exact-source
results, verify changed-source/platform scopes normally, and handle the CI gap
separately. Unrelated `outputs/` and primary-checkout work are preserved; the
previous goal stays paused. No push, dispatch, model call or native capture.

## J1 Review Follow-Up, 2026-11-01

Verify the integrated revisions for Claude's A1--A4 follow-up:
`8f0c1938` (off-gene interval/shared-end semantics), `363e23ac` (overlapping
unassessed causes with unknown legacy fields) and `84b1845c` (JS/Lua VM tests).
The [plan](observed_junction_consequences_plan.md#j1-review-follow-up-2026-11-01)
and [contract](observed_splice_junction_assessment.md) explain the exact scope.
The B-i aligner policy is documentation-only: preserve `b6863907` admission and
all Phase 2/S1 formats/semantics. J2/J3 and real biological-sample recurrence
selection remain gated. The untracked Claude prompt is not modified/included.

Reuse your completed [b6863907 return](#glen-return-at-b6863907-2026-10-30)
for unchanged S1/primer/glossary/lifecycle scopes. It already verifies the
rebased completeness path, contrary to the older prompt's first-runtime claim.
New J1 semantics/DTOs and VM tests need their own nonzero counts; retain the
[reported Slice 6 pass at d9b1015c](#claude-reported-lifecycle-slice-6-results-at-d9b1015c-2026-11-01)
separately rather than reissuing its completed checks. No new Codex Rust/VM execution or
whole-suite acceptance is claimed from source checks.

```sh
cargo test --locked --offline --no-default-features -j1 --lib junction_assessment -- --test-threads=1
cargo test --locked --offline --no-default-features --features script-interfaces -j1 --lib junction_assessment_js_vm -- --test-threads=1
cargo test --locked --offline --no-default-features --features script-interfaces -j1 --lib junction_assessment_lua_vm -- --test-threads=1
cargo test --locked --offline -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo check -q --locked --offline --no-default-features --lib --bin gentle_cli --bin gentle_mcp
cargo check -q --locked --offline --tests
```

Retain same-contig/cross-contig off-gene and boundary-touch refusals, near hits
without gene events, exact shared donor/acceptor versus two-end near relations,
unoriented start/end labels and complete/limited/incomplete synthetic donors.
The three cause counts overlap; their union, not sum, is `unassessed_records`.
Missing fields in old v1 evidence stay unknown. Real aligner sources retain
limitations/possible-only support; definite zero is expected, not absent RNA.

Both VM tests must actually execute assess plus bound artifact-verifying
inspect and preserve exact project bytes. A feature baseline failure is a
primary failure, not permission to skip silently or repair unrelated code.
Run/reuse source/documentation, LF/CRLF, glossary and canonical parity freshness
checks without regenerating unchanged metadata. Do not repeat completed gel,
historical sequence-design or other unaffected audits.

Return tested SHA/tree, lockfile and binary hashes, platform/toolchain/features,
commands, counts, failures/skips and retained evidence via Steffen. Full-suite,
actual CLI/MCP, native Windows/Linux, real CD274/private sequencing and scientific
acceptance stay separate. No new live-model/provider call is required or
authorized; no push, dispatch, publication or goal resumption is requested.

Authoring checks: 50 source/documentation passes, one explicit real-CLI skip;
documentation/manual/maturity freshness, scoped rustfmt and whitespace pass.
The lockfile hash remains the `a393a1f3...` value in your b6863907 return below.
Session-close: 4 OK, 2 expected working-tree/plan-review warnings, 0 failures;
roadmap stays 400 lines. None of this executes the new Rust or VM regressions.

## Glen Return At b6863907, 2026-10-30

Owner-supplied `GENtle_verification_b6863907.md` (SHA-256
`5b601b74b7f747604657962c4784d57ddd6ec47d7fa3eca3bbf1628af4fd4da9`)
reports exact source `b68639078614fd9ada856e95073ed77598e06647`, Debian
forky/Linux x86_64, Rust/Cargo 1.97.1 and lockfile SHA-256
`a393a1f3803b374bbaf5ef4f3e23813e974b1676f61a669dacb24063c059fb6a`.
The supplied handoff matched that source byte-for-byte. These are Glen's
reported runs, not new Codex Rust execution or acceptance of later source.

Retain without repeating unchanged checks: J1 `junction_assessment` 27/27,
`observed_transcript_structures` 70/70, `primer_payload` 103/103,
`candidate_sidecar_is_missing` 2/2, `lifecycle_slice3_` 11/11, glossary 12/12,
legacy cryptic-protein projection 1/1 and CLI `candidate_store_unresolved_` 1/1.
Filters overlap; do not sum them. Source checks passed 8 observed-transcript,
12 maturity and 18 documentation/LF/CRLF cases, plus the reported freshness and
whitespace checks.

**Primary failure:** healthy candidate-store Save As stages its sidecar under a
missing parent before the project writer creates that parent. The recovery test
fails alone. The other three aggregate failures are mutex-poison cascades and
each passes alone. Aggregates: `candidate_store_` 8 pass/4 fail; JS and Lua
`candidate_store_unresolved_` each 7 pass/4 fail, not successful VM gates.

[Same-source CI run 40939624811](https://github.com/smoe/gentle_rs/actions/runs/40939624811)
reportedly passed Release-Safety and Linux/headless, but failed Linux/Windows
full suites; macOS/native jobs were skipped. Raw CI logs were unavailable:
Linux focused reproduction confirms the defect, but Windows failure attribution
remains unconfirmed. No release/full-suite acceptance follows.

Current repair moves existing parent creation before sidecar staging, still
after the unresolved-store guard. New synthetic regressions cover populated
and empty healthy saves, missing nested parents, LF/CRLF source bytes,
spaces/apostrophes/Unicode, unchanged source/state, no staging leftovers and a
non-directory parent. Existing degraded-store no-mkdir oracles remain unchanged.
At this authoring checkpoint the repair was not Rust-verified; the later
[d9b1015c report](#claude-reported-lifecycle-slice-6-results-at-d9b1015c-2026-11-01)
supersedes that status only for its specified scopes.
Local repair validation: 19 source/teaching checks and 20 documentation/checkout
checks passed; 15 binary-dependent tests explicitly skipped. Scoped rustfmt,
documentation/manual freshness and whitespace passed. Session-close: 4 OK,
2 expected working-tree/plan-review warnings, 0 failures; roadmap stays 400 lines.
No local Rust build/test/check was run.
Finalized on main on 2026-11-01 after repeating these same non-Rust checks.
The later report supplies scoped Rust/VM/process results; native and
newer-integrated-source verdicts remain separate.

J1 is implemented with scoped synthetic verification, not protein consequence
prediction. J2/J3 remain proposals: no current altered-ORF, transmembrane-anchor
or soluble-form assessment of CD274. Glen's proposed future evidence gates are
review input only, not authorization to implement them or select biological
thresholds. No real/private RNA analysis or scientific acceptance was performed.

## Lifecycle Slice 6, 2026-10-29

This reconciles the overlapping foundation committed on main as `1fdeb1d9`.
That commit separately recorded 37 Python/source passes and two new process
skips; retain those historical results at their source, not as execution of the
revised compatibility/tests below. Its old strict/count checklist is superseded.

Verify the [Claude-reconciled degraded-sidecar slice](operation_lifecycle_slice6_plan.md)
on the integrated source; do not repeat unrelated RNA or historical Slice 3/5
audits. The owner approved unknown schemas staying loadable but unavailable for
candidate use/save, legacy warning-only states failing closed, and no general
pre-mutation refusal. Existing strict hydration errors and count defaults stay
unchanged. Shared access/save guards, unavailable facts and agent notes, visible
GUI failures and one-shot persistence/non-rollback feedback are authored.

The authoring-only statements below are superseded within the
[d9b1015c reported return](#claude-reported-lifecycle-slice-6-results-at-d9b1015c-2026-11-01).
The commands are retained as the acceptance baseline, not a request to repeat
completed unchanged checks. Verify affected newer-source/platform/native gaps
with their own nonzero counts and source/binary/toolchain/features/platform
receipts; do not reissue the RNA/primer/lifecycle audit campaign:

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_mcp
cargo test --locked --no-default-features -j1 --lib candidate_store_ -- --test-threads=1
cargo test --locked -j1 --lib lifecycle_slice6_gui_save -- --test-threads=1
cargo test --locked --no-default-features --features js-interface -j1 --lib candidate_store_ -- --test-threads=1
cargo test --locked --no-default-features --features lua-interface -j1 --lib candidate_store_ -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
GENTLE_LIFECYCLE_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_operation_lifecycle_slice3 -v
```

Run the formerly failing recovery test alone first if an aggregate fails;
separate primary failures from poisoned-lock cascades. The `candidate_store_`
filter includes both new `candidate_store_save_as_` tests and the existing
unresolved-store controls; VM filters overlap, not a sum of independent totals.
Two GUI helper tests cover early refusal and an ordinary I/O failure; one CLI
and one MCP test exercise structured feedback. VM tests use public loaders on
warning-only files and try same-path save plus Save As. Three new process tests
cover CLI/export/import/shared shell, confirmed MCP, Python and survival of an
explicitly written report after a later refused autosave. No supplied binary
means those process checks skip, not pass; never auto-build from Python.

Verify exact retained reference/diagnostic, untouched derived folder and intact
sibling files, Save As with a foreign decoy and absent parents, legacy warning
and unknown-schema load compatibility, candidate unknown-versus-empty, history/
snapshot/serde status, restore-and-reopen recovery and intentional deletion as
the preservation oracle's negative control. Do not require newly rejected count
defaults or claim parseable truncation is detected. Preserve healthy payloads.

Native GUI/Windows remain separate: show save errors, keep dirty/path/pending
actions unchanged, inspect unaffected DNA and preserve correct path/handle
behavior. Retain the backslash-path refusal case without a separator migration.
Keep new-source full-suite, package, scientific and performance verdicts
separate. Check scoped formatting against the baseline; do not mass-format
unrelated files. No push, dispatch, model call, capture or goal resumption here.

Historical authoring validation before Glen's return: 19 source/teaching passes,
15 explicit optional binary skips, two real LF/CRLF checkout passes, documentation/manual freshness
and whitespace OK. Formatting is scoped to changed Rust lines; no build/type-
check or current runtime acceptance follows. Session-close: 4 OK, 2 expected
working-tree/plan-review warnings, 0 failures; unrelated `outputs/` is untouched.

## Claude-Reported Lifecycle Slice 6 Results At d9b1015c, 2026-11-01

Owner-forwarded feedback and its accompanying uncommitted note report the
[Slice 6 commands](#lifecycle-slice-6-2026-10-29) at
`d9b1015cf27c7dfc2c8b73a675c11c232ded3bf4`, which includes the `82fb250c`
Save As parent fix. Platform: macOS 27.2/ARM64, rustc 1.101.0-nightly
(41089d07a 2026-10-27). Cargo calls used `--locked`; test invocations used
`--test-threads=1`, not `-j1`. `Cargo.lock` was unchanged. Counts/toolchain are
attributed to Claude's supplied note, not independently rerun by Codex; no raw
logs or binary hashes were supplied with this return.

- Both `cargo check` commands: clean.
- No-default `candidate_store_`: 14 passed. GUI `lifecycle_slice6_gui_save`: 2
  passed. JS and Lua `candidate_store_`: 15 passed each; these overlap and
  include both VM tests.
- Parity: 21 + 1 passed. Extra controls: `candidate_sidecar_is_missing` 2
  passed, `--bin gentle_cli candidate_store_unresolved_` 1 passed.
- Process suite against same-source `target/debug` binaries, first run: 22 of 24
  tests passed. The degraded-store test failed in all six LF/CRLF × fault
  subtests, and the legacy-warning test errored. Both were harness defects.
- Corrected companion at the same binary SHA: **24/24 tests passed, zero
  failures**, including all six LF/CRLF × fault cases and legacy recovery.
  This is the combined source/process companion, not 24 independent process
  cases; do not add its counts to overlapping filters.

The harness now expects the specific `GentleCliError` and checks `.stderr`,
not `str(exception)`, which omits CLI stderr. Recovery Save As explicitly opts
out of JSON parsing for its successful plain-text message; all other successful
calls keep strict JSON parsing and exit-status assertions. No production code
changed. The final directory/file-byte preservation oracles, previously not
reached, now reportedly run and pass. Both script corrections were uncommitted
when Claude supplied this result; execution binds d9b1015c binaries plus that
harness overlay, not a claim of byte-identical committed test scripts.

Two subsequent offline mock regressions guard these harness contracts without
spawning binaries. They are not part of Claude's historical 24-test count.
Codex's local authoring check ran the companion with `GENTLE_LIFECYCLE_BIN_DIR`
unset plus `DocumentationStatusTests`: 33 passed, 11 explicit process skips.
Documentation/manual freshness and whitespace passed. Session-close: 4 OK,
2 expected working-tree/plan-review warnings, 0 failures; roadmap stays 400 lines.

Native GUI, native Windows, the full suite and packaging remain unrun. Do not
promote these counts to newer-main/SAM or J1 VM acceptance. Glen's remaining
native checks still require visible save refusal, unchanged dirty/path/pending
actions and Windows path/handle behavior. No renewed same-source build campaign
is requested solely to repeat these passes.

## Claude-Reported Lifecycle Rust Results At 284e5b70, 2026-10-29

Owner-forwarded older review, received 2026-10-29: Slice 5 is implemented in
`b0521123` with additional tests in `98f69995`, not a pending implementation
plan. Claude reports this command at
`284e5b7049016e42ea19dd4bd50d7e93516d7ca7`, macOS/default features:

```sh
cargo test --locked --lib -- lifecycle_slice5 lifecycle_slice3
```

**Reported: 28 passed, zero failed**: 17 Slice 5 tests, including three GUI
consumer unit tests, and 11 Slice 3 controls; `Cargo.lock` unchanged. Attribute
this result to Claude's supplied report. No raw log, exact toolchain or test
binary identity was supplied, and Codex did not rerun Rust. Counts are separate
from the 25-check binary rehearsal below; do not sum or promote them to native,
headless, whole-suite or newer-main acceptance.

Codex's current source inspection confirms the reported safeguards: project
snapshots for Preview, engine-read/job-lock publication order, owner/structural/
cancel/request checks, GUI result-owner gating and retained freshness stamps
for explicit sharing. The reviewed Preview branch, GUI consumer and test
bodies are unchanged since the reported source; shared files gained unrelated
junction support. That is source continuity, not execution of the newer source.
Claude's C1--C7 findings are closed in that scoped shipped-code review.

Formatting: Claude reports scoped rustfmt passing and whole-repository
`cargo fmt --check` failing with about 78 unrelated diffs at that source. Keep
the whole-repository check and its failure visible, distinguish baseline drift
from changed-file failures, and avoid unrelated formatting churn. Do not reuse
the historical diff count as a current measurement. Cargo accepts one filter
before `--`; the two filters above are passed after it to the test binary.

Remaining work uses the [existing Slice 5 commands](#lifecycle-slice-5-2026-10-28),
reusing unaffected source-matched receipts where available: headless test-target
check, uncovered sequence-design/lifecycle/agent suites, JS/Lua VM and parity
checks, native Windows and affected-source CI. Do not repeat the reported
28-test run solely because older authoring notes say unrun.

Native long-job concurrency/cancellation is still inconclusive. The rehearsal
below observed no History row, historical sharing, project replacement and
Apply/undo/redo, but did not complete the long-job/cancel or full structural-
versus-display native checks. Do not classify all native work as unrun or all
as accepted. Keep the brief fail-closed "Project is busy" display as a cosmetic
observation, not authority to weaken lock/owner checks or expand this repair.
The previous goal remains paused; no new build, model call, capture or push.

## Help And Observed-Transcript Smoke Fixtures, 2026-10-28

The owner reports 4706 library passes, three failures and 16 ignored tests;
the exact tested SHA/platform were not supplied. Inspection of integrated
`284e5b70` confirms two fixture defects, not authority to relax production
schemas or change a navigation destination assertion.

Help's pointer helper used the center of an inline galley's full bounds,
including the leading space occupied by earlier links. It now chooses a
nonempty row's rectangle without leading space. The existing embedded-book
click must still land at `docs/gui/documentation-automation-status.md`, not
`what-to-trust-today.md`. A new synthetic hyperlink-response test proves the
old center misses the link while the selected glyph-row point is inside it.
Browser, screenshot and unchanged-project assertions remain intact.

The four annotation/structure prepare/cancel smoke routes now use their existing
discovery request templates instead of generic `{}`. Inspection `--binding`
gets a typed, report-path-bound synthetic record. New controls require explicit
write confirmation, unchanged quoted POSIX/Windows path spelling, metadata-only
inspection defaults and rejection of empty/unlabelled/unknown-field requests
and malformed bindings. Discovery placeholders/zero budgets remain nonexecutable;
the parser, runtime validation and public contract were not broadened.

Local verification: 25 Python/source tests passed, including LF/CRLF manual
coverage and observed-transcript contract guards. Documentation status, manual
inventory freshness and whitespace checks passed. These do not execute the new
Rust regressions. Session-close reports four OK, two expected worktree/manual-
review warnings and no failures; unrelated dirty/untracked files are preserved.

On the next normally integrated source, run only these affected groups, or reuse
same-source CI with nonzero counts:

```sh
cargo test --locked --lib -j1 help_docs -- --test-threads=1
cargo test --locked --lib -j1 glossary_ -- --test-threads=1
```

They include the three original failures and the new geometry/template/binding
regressions. Do not repeat completed unrelated RNA/gel/lifecycle checks. Rust
and native execution of this repair remain unrun locally; return full tested
SHA, toolchain/platform, counts, failures/skips and retained evidence through
Steffen. No new audit freeze, push, package/model/provider call or goal resumption.

## Synthetic Lifecycle Rehearsal At 284e5b70, 2026-10-28

Owner authorizes this local rehearsal while Glen is unavailable, then explicitly
approves an isolated GUI profile and transient, unretained screen captures.
Source `284e5b7049016e42ea19dd4bd50d7e93516d7ca7`; all three existing main-checkout
binaries report that revision and build `1793215065252`. Hashes remain identical
before/after. This is binary self-report plus content identity, not a clean-build,
feature/profile/toolchain or package attestation. No Rust build or model call.

```sh
GENTLE_LIFECYCLE_BIN_DIR=/Users/u005069/GitHub/gentle_rs/target/debug \
GENTLE_TUTORIAL_BIN_DIR=/Users/u005069/GitHub/gentle_rs/target/debug \
python3 -m unittest scripts.test_operation_lifecycle_slice3 scripts.test_sequence_design_tutorial -v
```

Executed on macOS 27.0.1/ARM64 with disposable HOME/USERPROFILE/XDG paths:
**25 passed, zero skips** (13 source/teaching, 12 actual-process checks).
Includes both strategies, 18 exhaustive cases each, refusal/tampering/receipts,
project and sidecar byte preservation, LF/CRLF/quoted paths and macOS aliases.
The companions delete their per-test temporary state; do not describe this log
as an archive of every subprocess payload. The later alias teaching regression
is a separate source check, not a twenty-sixth original process-run result.

Scoped native observations using the same GUI bytes in a temporary app wrapper:
no-file preview completes without a History row or undo/redo; an unrelated edit
marks the retained output historical; explicit draft sharing includes that
warning without sending; project replacement clears the old result. Fresh
file preview yields exact `ATGTTTAAGTAA`/one base-9 edit/all seven GC windows;
reviewed Apply creates one sequence/container, same-session undo removes both,
and redo restores the 12-bp native viewer. No screenshots are retained, so these
are recorded observations, not a raw-image/semantic-snapshot acceptance bundle.

**Native cancellation/concurrent long-job responsiveness is inconclusive.**
An additional bounded synthetic workflow probe hit native-automation timeouts;
no active admission/cancel/terminal receipt was visibly established. A stack
sample found rendering rather than a captured design worker, and a later
SIGUSR1 dump had no active frames. Neither proves an engine defect or successful
cancel. Stop only the rehearsal process; do not waive or claim this acceptance.

Evidence outside Git:
`/private/tmp/gentle-lifecycle-rehearsal-284e5b70.R2QDjA/` (`REPORT.md`, runner
log, version outputs, source/binary bindings and synthetic input/report/sample).
`SHA256SUMS` SHA-256:
`c723f0e43a881c2e096de526ecc94c65621c1ec4ac1121d1b47ffe3d17506988`.
GUI binary hash `a96165d907b008731a3442ebff72249b00488f45f9638bbdfbc58519f9b1a71a`;
CLI/MCP hashes are in that report. Content hashes do not authenticate provenance.
Keep native cancellation/concurrency, Windows/Linux/package, full-suite, VM and
live-model/scientific/performance verdicts separate. The owner's successful local
CI report is not independently reconstructed here. Older authoring-time unrun
lists below are superseded only for the explicitly executed checks above.
Tutorial 06.07 now distinguishes GUI Shell history commands from Agent Assistant
`/history`, `/undo`, `/redo`; no runtime/approval/scientific contract changes.
No push, goal resumption, unrelated-output edits or new release gate.

## Lifecycle Compile And VM Oracles, 2026-10-28

Narrow repair of the owner-supplied Claude review at historical `550e6e0a`.
The missing Slice 4 test-module path was already fixed in `a52c7b74`.
`496830db` fixes Help navigation visibility and the one `FullOutput` command
access; other `PlatformOutput` accesses stay unchanged. The companion test
repair decodes Lua results as `ProjectState` before JSON comparison and feeds
JavaScript through the public `load_project` binding, not project JSON literals.
JS equality explicitly inspects typed-array and ArrayBuffer bytes. Production
JS/Lua adapter source is byte-identical before/after this test-only repair.

Run the following seven commands on the normally integrated full SHA, recording
toolchain/platform/features, counts and primary failures. A zero-match filter
is a failure. These commands have **not** been run in this Codex session; the
quoted review's request to build did not override the standing no-local-Rust
instruction. Older source-bound passes are not results for these repaired tests.

```sh
cargo check --locked -q
cargo check --locked -q --no-default-features --lib
cargo test --locked --lib lifecycle_slice -- --test-threads=1
cargo test --locked --lib --features js-interface lifecycle_slice -- --test-threads=1
cargo test --locked --lib --features lua-interface lifecycle_slice -- --test-threads=1
cargo test --locked --lib help_docs -- --test-threads=1
cargo fmt --check
```

Require all four repaired VM tests to execute. Their non-vacuous controls are:

| Test | Preservation oracle and counterexample |
| --- | --- |
| `lifecycle_slice3_js_runtime_preserves_state_and_legacy_controls` | Loader-backed input and independent baseline; covered operations/workflows preserve state, while mutation, empty and mixed workflows do not. Refusal keeps the input. A/C byte-view and buffer counterexamples exercise byte sensitivity. |
| `lifecycle_slice4_js_runtime_state_only_writer_is_fresh_only` | Exact output DNA and unchanged state; empty workflow must reconcile and differ under the same oracle. Existing destination must refuse without changing input or retained report bytes. |
| `lifecycle_slice3_lua_runtime_preserves_state_and_legacy_controls` | Decode each operation/workflow result and refused input as `ProjectState`; empty and mixed mutation workflows must still differ from the same typed baseline. |
| `lifecycle_slice4_lua_runtime_state_only_writer_is_fresh_only` | Typed state equality; empty workflow must differ under the same oracle. Existing destination must refuse with unchanged input and report bytes. |

These controls are authored, not an executed red/green demonstration. Glen/CI
must establish their runtime results; no production adapter conversion was
added to make them pass. JS project JSON-literal ingress remains an explicitly
deferred compatibility question; this repair uses the existing public loader.
Separate any pre-existing repository-wide formatting drift from this delta.
Python actual-binary tests, the full Rust suite, native GUI and native Windows
remain unverified here. Reuse unaffected existing receipts; do not launch a
parallel audit campaign. No push/CI dispatch, goal resumption or private data.

Local non-Rust checks: 13 lifecycle/tutorial source tests pass with 12 actual-
process skips; 19 documentation/real LF/CRLF-checkout checks pass. Manual
projection freshness, documentation status, scoped JS/Lua rustfmt and whitespace
pass. These checks do not execute a scripting VM or establish compilation.

## Lifecycle Slice 5, 2026-10-28

Include the owner's approved [Slice 5](operation_lifecycle_slice5_plan.md) in
the normal integrated-source lifecycle return, not a parallel audit campaign.
Only managed standalone no-path Plan changes: it retains the project snapshot
and execution capacity, but publishes a source-stamped command result with zero
live journal commits/History rows. Raw adapters, with-path Plan, S1, workflows,
opaque wrappers and arbitrary trusted work retain their contracts. No new
payload-resolution layer, scientific schema, journal merge or queue.

Changed runtime: `src/command_execution.rs`, its existing policy helper,
`src/app.rs` and `src/app/routine_and_agent_assistant_ui.rs`.
Inner-agent guidance in `src/agent_bridge.rs` retains full existing instructions
and adds only the snapshot/History/historical-evidence distinction.
New synthetic service tests: `src/command_execution/preview_tests.rs`; existing policy,
lifecycle and GUI test modules cover request binding, noncoverage and delayed
sharing. Test constructors in `src/app/tests.rs` supply no preview binding for
unrelated outputs. Tutorial 06.07, protocol, canonical GUI chapters, decisions/
architecture and lifecycle docs teach the bounded command audit versus History.

Local Python/source checks pass 25 focused cases with 12 process skips (37
total) and all 83 LF/CRLF checkout cases. Manual/maturity freshness, documentation
status and whitespace pass. Scoped rustfmt passes except the preserved existing
`app.rs` Help formatting groups: comparison to HEAD confirms no new format drift.
Whole-file `cargo fmt --check` may report those three baseline groups; keep that
separate from Slice 5 failures. No Rust/runtime/native acceptance is implied.
Session-close has 4 OK, 2 expected working-tree/manual-review warnings and no
failures. Unrelated `outputs/`, the paused goal and existing stashes are untouched.

The authoring-time checks above did not execute Rust/native/process/VM tests.
The [later Claude-reported return](#claude-reported-lifecycle-rust-results-at-284e5b70-2026-10-29)
and [process/native rehearsal](#synthetic-lifecycle-rehearsal-at-284e5b70-2026-10-28)
supersede that unrun status only for their named source/subsets. Use nonzero
counts for remaining checks at the integrated full SHA; do not repeat every
older evidence bundle:

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --tests
cargo test --locked -j1 --lib command_execution:: -- --test-threads=1
cargo test --locked -j1 --lib routine_and_agent_assistant_ui::command_tests -- --test-threads=1
cargo test --locked -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked -j1 --lib sequence_design -- --test-threads=1
cargo test --locked -j1 --lib operation_lifecycle -- --test-threads=1
cargo test --locked -j1 --lib agent_execution::tests::sequence_design -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo fmt --check
```

Exercise both service submit forms, named/raw inline/loaded previews and worker
`@file` input; hold real mutations and simultaneous previews with channels.
Require exact direct reports/approval hashes and unchanged live journal/state/
revisions/undo/redo. Source/replacement changes before publication reject; display/
metadata edits and Save As do not. Check cancellation before snapshot, during
search, after computed work but before publication, and after publication.
Exhausted/infeasible scientific reports are not feasible designs even if the
service completed. Oversize refusal has an infrastructure test, not a claimed
reachable full-solver case. Preserve the unchanged S1 and workflow stale guards.

GUI/native: consume after edit/replacement, then edit/replace again before
display/sharing. Historical warnings must reach the explicit inner-agent draft;
wrong-owner/busy-project sharing must fail without changing it. No auto-send,
model call or Apply approval follows. Check a long job remains usable beside a
preview and that Operation History gains no managed no-path row; with-path
reports/Apply still work. Keep existing raw CLI/MCP/JS/Lua/Python companions and
scientific parity passing. Return exact source/binary/platform/features/counts,
failures and unrun checks. No push/dispatch or previous-goal resumption here.

Additional source-pass tests authored after `b0521123` (unrun at authoring,
included in Claude's reported Slice 5 return at `284e5b70`):
`lifecycle_slice5_preview_tolerates_a_read_only_journal_append_after_snapshot`
must preserve the intervening S1 audit and the live engine/history, and
`lifecycle_slice5_published_preview_apply_revalidates_only_its_actual_constraints`
must accept an unchanged source after an unrelated edit but refuse changed
source DNA, changed annotations and an occupied output ID. They use the same
command-execution filter above; no new audit campaign or production behavior.
Tutorial 06.07 distinguishes historical context from actual Apply invalidity.

## Integrated Lifecycle And Help Follow-Up, 2026-10-28

After the owner's next push, use the normally integrated source containing
`1adad235` (Lifecycle Slice 4) and `7329f00f` (manual navigation/contextual Help).
Record the full tested SHA once; these commit IDs identify changes, not an
interim freeze or a request for another bundle. The Slice 4 and manual sections
below define this delta. Older sections retain their own source-bound results and outstanding
limits; they are not a blanket request to repeat completed checks.

Already completed, not new Glen work: at `7329f00f`, 164 Python/source/checkout/
synthetic-package checks passed, with six optional real-CLI skips (170 total).
Both manual migration checks retained all 7,171 GUI and 8,869 CLI nonblank lines
from `1adad235`, including indentation/multiplicity after documented link moves.
Manual projection freshness, documentation status and whitespace passed.
Slice 4's earlier checks are recorded below; none of these results executes
the pending Rust, real CLI/MCP/VM or native acceptance paths. Reuse same-source
CI receipts and run each missing command/feature group only once; rerun an
affected completed check only when its source or inputs have changed.

Return one concise report through Steffen to this main-maintaining chat:
per-item pass/fail/not-run/blocked, full source SHA, platform/toolchain/features/
profile, binary and input hashes, exact commands/counts/skips, primary failures
versus cascades, and retained artifact locations/checksums. No push, dispatch,
package build, tag, publication, live-model call or scientific approval is
requested. Keep the previous goal paused and unrelated outputs/stashes intact.

## Glen Scoped Scan At b8e0118f, Report Received 2026-10-28

Owner-relayed auto-report, not a new local run: Glen inspected
`ce4d33b5b990701f2ce51c85f6ce2821df9dc8ce..b8e0118feb44726392702398be6b6d52218eec66`
(six commits, 42 files) in a fresh worktree. The earlier source remains an
ancestor. He reports no concrete regression and recommends no patch.

Reported passes at exactly `b8e0118feb44726392702398be6b6d52218eec66`:

- Observed-transcript structure tests: 23/23; this predates S1 consumption.
- Lifecycle Slice 3: 11/11.
- Dotplot focus and resize: 2/2; not native window acceptance.
- Python lifecycle source tests: three passed; five binary tests explicitly skipped.
- Diff/whitespace check passed. Rust used one Cargo job and an explicitly
  isolated target. An initial inherited-shared-target attempt was stopped;
  all intended Rust groups were then rerun completely in that isolated target.

Glen reports removing the run-owned worktree and 7.5 GiB Cargo cache, leaving
the shared cache untouched. The supplied summary does not include raw streams,
artifact hashes, exact commands, platform, feature set or toolchain version;
retain those limits rather than reconstructing them. There is no full-suite,
cross-platform, real-CLI/VM, native GUI, package or scientific verdict here.

This SHA precedes `929c54e4` S1, `1adad235` Slice 4 and `7329f00f` manual Help.
Keep their missing checks open. Do not repeat these completed scoped groups
merely because an older section still lists them; use normal integrated CI for
affected-source regressions, keeping the new result separate from this receipt.

## Lifecycle Slice 4, 2026-10-27

Owner authorized the revised Plan-only fresh-file writer/host binding at the
`b8e0118f` baseline; implementation `1adad235` is now integrated beneath
`7329f00f`. Test the normally integrated source, including the existing
`a2637079` Slice 3 test-API repair. The [contract](operation_lifecycle_slice4_plan.md) records
support limits; tutorial 06.07 now requires fresh files in existing folders.
The owner-requested mid-implementation rebase advanced this branch to
`1804605f`; all Slice 4 working files restored cleanly, preserving RNA S1.
The scoped pre-rebase recovery stash `617a89ee` remains retained; unrelated
outputs and older stashes were not applied. No new verification follows from
the rebase or main's separately recorded RNA checks.
Post-rebase freshness checks required the four generated maturity projections
to adopt main's current glossary-source hash only; no capability definition,
maturity verdict or evidence changed.

Local non-Rust checks after the rebase: the focused lifecycle/tutorial/maturity
Python suite passes 22 checks with 12 explicit actual-process skips (34 total);
all 75 checkout tests pass, including LF/CRLF coverage. Maturity freshness,
documentation status, scoped rustfmt and whitespace checks pass. These execute
neither GENtle nor its Rust/VM/native tests; the checkout runner's synthetic
binary fixtures are not real-binary acceptance.
Session-close reports 4 OK, 2 warnings, 0 failures: intentional uncommitted work/
generated hashes and the manual plan-fidelity reminder. The only extra generated
refresh is reproduced by `python3 scripts/generate_capability_maturity.py --write`;
it repairs the post-rebase source-hash mismatch documented above.

Changed runtime: `src/engine.rs`, new
`src/engine/io/sequence_design_preview.rs`, the Plan dispatch arm,
`src/{engine_shell,agent_bridge,app,mcp_server}.rs` and `src/bin/gentle_cli.rs`.
The new direct `same-file` dependency reuses locked 1.0.6, changing only the root
dependency edge. Authored tests: new writer `sequence_design_preview/tests.rs`,
service/GUI/MCP/JS/Lua cases, the existing sequence-design tutorial regression
and `scripts/test_operation_lifecycle_slice3.py`. Documents include tutorial/
input provenance, protocol/CLI/GUI/architecture/decisions, lifecycle plans/status,
roadmap/changelog and this handoff. Private data and unrelated outputs are untouched.

Rust, actual-process/VM and native execution was not run locally. Reuse same-source
normal CI for groups already executed; run missing groups with nonzero counts:

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_mcp
cargo test --locked -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked -j1 --lib sequence_design -- --test-threads=1
cargo test --locked -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked -j1 --lib command_execution -- --test-threads=1
cargo test --locked --no-default-features --features js-interface -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked --no-default-features --features lua-interface -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo fmt --check
```

Build source-matched CLI/MCP once if needed, then run both actual-process
companions, including Python delegation/default state paths. Missing binaries
or prerequisites are unrun/blocked evidence, never parity:

```sh
cargo build --locked --no-default-features -j1 --bin gentle_cli --bin gentle_mcp
GENTLE_LIFECYCLE_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_operation_lifecycle_slice3 -v
GENTLE_TUTORIAL_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_sequence_design_tutorial -v
```

Require exact project/input/sidecar bytes, report equality, confirmation before
I/O, missing/default project refusal, derived-tree exclusion, existing/link/
special-file refusal, final-name race and injected staging errors. Keep both
strategies and distinct report paths; no expected scientific hash changes.
Verify host binding survives load/pool reconstruction, checkpoint import/fork,
rollback and newer live Save As commit, without serialization or dirty revisions.
Check native GUI load/Save As/new project/reset, undo/redo and responsiveness; Windows
handles/device/ADS/verbatim/UNC/reparse cases and macOS aliases/case/Unicode
remain actual native acceptance, not source or LF/CRLF proof. Retain inner-agent
ask/nested-execution guards and capability discovery; no model call is needed.

No overwrite mode, implicit mkdir, loader-fact context or payload resolver was
added. Same-parent non-ASCII protected-name comparisons intentionally over-refuse:
use a different existing directory. No hostile parent-swap, universal atomicity,
crash durability or survival against later explicit project save claim. Record
tested full SHA, feature/platform/toolchain, counts, failures and unrun checks;
keep the previous goal paused and coordination/push with the main thread.

## Manual Navigation And Contextual Help, 2026-10-28

`7329f00f` adds 73 GUI and 69 CLI authoritative topic files, a shared catalog,
legacy anchors, reciprocal links, complete-book search/copy, cross-book history
and 40 typed window Help destinations. Ordinary Help must remain separate from
Agent Help, screenshot capture, model calls and project mutation. See the
[implementation/acceptance ledger](manual_navigation_plan.md#external-acceptance-handoff).
The static checks recorded above are complete; actual Rust/GUI behavior remains
unverified. No engine/scientific contract change is claimed by the manual move.

Owner-reported compilation failed with E0624 at the parent toolbar's
`navigate_manual_history` calls and E0609 in the no-browser test. The scoped
repair gives the helper `pub(super)` visibility and reads
`FullOutput.platform_output.commands`, as defined by locked egui 0.36.2.
No navigation logic or assertion was removed; viewport screenshot and unchanged
project checks remain. Run the already-listed default-feature test-target check
first, then the missing Help regressions below. This repair is not a compilation
pass; no local Rust execution occurred.

Run these missing default-feature groups, or cite their same-source CI results.
Require nonzero selected counts. Shared Cargo checks/formatting are already in
the Slice 4 list; do not run them twice. Its default-feature `sequence_design`
filter also covers `sequence_design_window_gc_documented_walkthrough_is_replayable`,
so do not rerun that walkthrough separately. `help_docs` includes the real
topic-button dispatch, active-window registry, cross-book history and
package-layout cases:

```sh
cargo test --locked -j1 --lib help_manual -- --test-threads=1
cargo test --locked -j1 --lib help_docs -- --test-threads=1
cargo test --locked -j1 --lib cli_docs_mention_all_cli_glossary_paths -- --test-threads=1
cargo test --locked -j1 --lib genome_edit_intent_documented_operation_is_the_tested_synthetic_request -- --test-threads=1
```

Native checks still needed: narrow/wide Help windows, Contents/tree and
breadcrumbs, both directions of GUI/CLI links, Back/Forward, old bookmarks,
UTF-8/duplicate headings, whole-book Find/Copy, Reload, images and wheel/trackpad
scrolling. Check ordinary Help in sequence, specialist and detached windows,
including the active-window native menu route; no capture/model request or
project/dirty/undo change may result. Headless clicks are not native acceptance.

Keep installed-package checks on the normal authorized release path, not a new
package campaign. On macOS/Linux/Windows, launch from an unrelated working
directory with complete resources and with a missing chapter/incompatible
catalog; require one complete embedded-book fallback and a visible notice,
never an index-only or mixed-source book. Test on a disposable package copy,
not installed originals. Text embedding does not supply images/tutorial assets;
their installed-resource requirement and other guides' section-scroll limits
remain explicit. Record unavailable platforms/packages as pending, not passed.

## DNA Map And Dotplot Separation, 2026-10-27

The owner's screenshot shows the compact `Dotplot map` launcher in the main
DNA pane after a reasoning inspection opened a separate dotplot. Source review
found an explicit primary-mode assignment, not evidence of cross-window paint
corruption. Auxiliary analysis/feature/reasoning routes now preserve the primary
view; the explicit toolbar mode remains. The launcher fills its available
height and scrolls within narrow/short panes. Existing saved mode choices are
not migrated; select `Standard map` once if the old launcher is still selected.

Glen's [scoped scan](#glen-scoped-scan-at-b8e0118f-report-received-2026-10-28)
reports the focus/resize pair passing at `b8e0118f`. Do not call those two
checks unrun or launch another old-SHA audit; only affected-source regression
through normal integrated CI remains. The two contextual routes below have no
individually itemized result in that return; run them only if missing from
same-source CI. Local Rust execution remains prohibited:

```sh
cargo test --locked --lib -j1 open_dotplot_for_feature_opens_transcript_dotplot_on_explicit_request -- --test-threads=1
cargo test --locked --lib -j1 open_construct_reasoning_dotplot_action_focuses_repeat_region -- --test-threads=1
```

The resize regression drives the actual launcher for two frames at each of
1000x720, 360x180 and 780x940 logical body sizes; it must fill height without
horizontal body growth. Native macOS and hosted-window follow-up: begin in
Standard/Splicing, open a repeat/feature dotplot and focus stored analyses,
verify separate content and unchanged primary mode, resize both windows,
close/reopen Dotplot and toggle sequence text. Explicit `Dotplot map` must still
work, scroll and use the full remaining pane; restore Standard without closing
the auxiliary window. Retain existing repeat-region/provenance checks and
project DNA. Return source/binary/platform/profile and results with the normal
handoff; no private input, live-model run or new audit campaign is needed.

## Observed Transcript Structures S1, 2026-10-27

Include this in the GitHub-pushing thread's **single consolidated request**,
not a parallel per-branch acceptance run. S1 was committed as `a363f1bd` and
rebased onto main `b8e0118f` as `929c54e4`; the [post-rebase receipt](observed_transcript_structures_plan.md#s1-commit-and-rebase)
records repeated checks and separate source/binary fingerprints. For later
integration, record that exact SHA/binary before missing checks. No push/dispatch or new freeze is requested
by this note. [Contract and commands](observed_transcript_structures.md) and
[source ledger](observed_transcript_structures_plan.md#s1-source-implementation-ledger-2026-10-27)
separate this work from historical S0 acceptance. The owner subsequently
authorized local verification: [the executed receipt](observed_transcript_structures_plan.md#s1-executed-verification-2026-10-27)
records macOS S1 14/14, S0 24/24, primer-payload 46/46, lifecycle 11/11,
registry introspection 1/1 and canonical parity 1/1, plus headless CLI/MCP/library
and default-feature test-target Cargo check. Filters overlap; these passes
belong to the fingerprinted working diff, not later main or all platforms.

Reuse unchanged same-source evidence; rerun affected filters/checks after
integration or source changes, retaining the ten old Phase 2 byte/work cases.
The local tests cover source-bound
intake without original FASTQ/PAF/FASTA, configurable gene context, source/model
separation, whole-view unavailability, exact accounting, LF/CRLF/raw hashes,
cancel/fault/missing counts, no-overwrite output and saved-only inspection.
Raw apply retains its ordinary audit; standalone managed publication must keep
the live project/journal/redo unchanged and retain a published result after a
late host cancel. Generic MCP confirmation remains required.

The canonical producer/document equality and full-registry introspection now
pass; the standalone documentation-helper CLI wrapper was not built/run.
Exercise actual CLI/MCP processes with quoted/spaced/apostrophe paths and both
absent and existing project files; compare project/sidecar bytes. Recheck generic
JS/Lua state-out and all-covered versus mixed workflow behavior without treating
a managed workflow as standalone artifact-only publication. These actual-process
checks remain missing, not replaced by local in-process MCP/shell or the
separate lifecycle slice. Small managed responses retain typed metadata and
the exact binding in the job receipt; only oversized responses are compacted.
Return exact source/platform/toolchain/profile, command/count/failure records and
artifact hashes. No new S1 GUI, BAM/SAM, private-data replay, screenshots,
protein claim or scientific acceptance is requested. Native/Windows outcomes
outside an executed path remain unrun, even if the headless synthetic run passes.

## Lifecycle Slice 3, 2026-10-27

Glen reports 11/11 focused passes at `b8e0118f`; see the
[source-bound receipt](#glen-scoped-scan-at-b8e0118f-report-received-2026-10-28).
The remaining real-process, optional scripting-runtime and native checks are
not supplied by that result. The command list below retains the feature scopes
for missing/affected checks, not a request to repeat the completed old-SHA run.

Owner-authorized implementation starts from `ce4d33b5`; integrate normally and
record the exact tested full SHA. No new payload resolver or parser change.
Verify Plan/S1 request-aware construction/apply, raw no-autosave, generic
JS/Lua state-out and existing journal-publishing managed jobs. The
[reconciled plan](operation_lifecycle_slice3_plan.md) owns scope and limitations.
Keep the earlier goal paused; push/coordination remains with the main thread.

Changed runtime files: `src/engine.rs`, `src/engine/operation_lifecycle.rs`,
`src/engine/ops/operation_handlers.rs`, `src/bin/gentle_cli.rs`,
`src/mcp_server.rs`, `src/js_interface.rs`, `src/lua_interface.rs`.
Tests: new `src/engine/lifecycle_slice3_tests.rs`,
`src/command_execution/lifecycle_tests.rs` and
`scripts/test_operation_lifecycle_slice3.py`, plus the module registration in
`src/command_execution.rs` and adapter-local regressions. Documentation:
`docs/{architecture,decisions,protocol,cli,gui,roadmap,CHANGELOG}.md`, this
handoff, the lifecycle proposal/plan/historical review prompt, tutorial 06.07
and its input README. Shared shell/parser, Python client and Cargo files are
unchanged. Unrelated `outputs/` remains untouched.

Use same-source CI where already available; run only missing groups. Rust,
actual binaries/VMs and native acceptance were not executed by Codex:

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_mcp
cargo test --locked --no-default-features -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib sequence_design -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib genome_edit_intent -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib command_execution -- --test-threads=1
cargo test --locked --no-default-features --features js-interface -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --no-default-features --features lua-interface -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo fmt --check
```

Require nonzero selected counts; blocked optional runtimes are not parity. Build
CLI/MCP once at that source if needed, then supply them explicitly (never let the
Python client fall back to Cargo):

```sh
cargo build --locked --no-default-features -j1 --bin gentle_cli --bin gentle_mcp
GENTLE_LIFECYCLE_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_operation_lifecycle_slice3 -v
GENTLE_TUTORIAL_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_sequence_design_tutorial -v
```

Check absent/empty/noncanonical LF/CRLF projects and candidate sidecar bytes,
MCP denial before I/O, preview/denied Apply absence and approved exact DNA,
unreconciled MCP summary counts, loaded/inline targets, failure/cancel prefixes,
real JS/Lua/Python, redo and managed owner/history/structural guards. Existing
managed project-replacement tests remain essential controls. A concurrent
journal append is still allowed to make another job stale. No live model call
is needed: retain existing inner-agent approval/ask and nested-agent tests.
Glen owns native GUI dirty/redo behavior and Windows handle/path checks.

Do not relabel old approval passes as current-source evidence. Opaque wrappers
still repair on construction; permissive sidecar warnings can still be saved
by unmigrated commands. Report-path alias protection was deferred at Slice 3;
it is now authored in Slice 4 above, with acceptance still pending. Report these
boundaries, not a universal read-only guarantee. Return
changed-source SHA, feature/platform/toolchain, exact selected counts, failures
and unrun checks through the main-maintaining thread; no new audit campaign.

## Tutorial Agent-Parity Approval Count, 2026-10-26

The owner reports 4549 library passes, one failure and 16 ignored tests; tested
source/platform were not supplied. The remaining assertion expected 15 parser
review gates, before `98c9331a` made tutorial 06.07's `sequence-design plan
@request.json --path preview.json` require report-write approval. The count is
now 16; five tutorials, 27 cases and eight author-declared project mutations
remain unchanged. The production gate and tutorial's `execution: ask` stay intact.

On the next normally integrated source, run only this missing regression (or
reuse a same-source normal-CI result); do not repeat completed audit groups:

```sh
cargo test --locked --lib -j1 tss_tutorial_agent_drafts_parse_and_keep_mutations_reviewed -- --test-threads=1
```

It now also changes the real 06.07 preview to `auto` in memory and requires
exactly its `auto_blocked_by_runtime` finding. No command executes or fixture
changes. The Python companion pins the existing reviewed/non-project-mutating
declaration, not the Rust parser's runtime result. Rust verification is pending;
include it in the normal source-bound return below, not a new freeze or bundle.

## Primer Payload Phase 2, 2026-10-26

**Forward to Glen after Steffen's normal push.** Verify the integrated source
containing query preparation `72885fbf`, interpretation `07389910`, mask repair
`e622f12c`, source/exclusion binding repair `f1a60ae4` and CLI documentation
repair `320d5803`. Record the fetched full SHA once; no interim freeze, recurring
bundle or synchronization campaign. These are generic, input-driven contracts,
not a special case for the current panel. Existing open gates below stay separate.

**Do not repeat completed work.** Retain historical gel/labware/six-gate and
earlier tutorial acceptance at their actual source. The exclusion continuation
already passed 130 Python/source checks with one actual-CLI skip; the CLI-doc
repair passed 16 with one actual-CLI skip. Those are not Rust/runtime passes.
Use existing same-source CI logs for any command below that already passed with
the required feature/platform scope; run only the missing groups, not overlapping
narrow filters as well as the full `primer_payload_` filter.

Read [the contract/claim matrix](primer_payload_evidence.md),
[implementation boundary](primer_payload_partner_phase2_plan_20261026.md#implementation-boundary-2026-10-26)
and [tutorial 10.04's synthetic Phase 2 recipe](tutorial/10-04_targeted_nanopore_fastq_panel.md#phase-2-primer-free-payload-evidence-synthetic-only).
Use a clean isolated checkout and explicitly run-owned `CARGO_TARGET_DIR`, never
an inherited shared target. Use invented fixtures only. No private FASTQ,
downloads, aligner execution, provider/model calls, packages or publication.

### Rust And Freshness: Run Missing Groups Once

```sh
cargo check --locked --all-features --tests -j1
cargo check --locked --no-default-features --lib --bin gentle_cli --bin gentle_examples_docs -j1
cargo test --locked --lib -j1 primer_payload_ -- --test-threads=1
cargo test --locked --lib -j1 primer_trace_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 glossary_ -- --test-threads=1
cargo test --locked --lib -j1 cli_docs_mention_all_cli_glossary_paths -- --test-threads=1
cargo test --locked --test capability_registry_parity --test parity_matrix_freshness -j1
```

The payload filter includes the initial 20 authored regressions plus five mask
and nine exclusion/budget regressions, including MCP confirmation refusal,
named lifecycle and agent-outcome tests. Require nonzero selected counts and
retain test names/results. Rust compilation/execution remains pending for this
payload slice; the owner-reported 4535/1/16 run has no supplied source/platform.

Build source-matched binaries once if not already available from that source:

```sh
cargo build --locked --no-default-features -j1 --bin gentle_cli --bin gentle_examples_docs --bin gentle_mcp
```

Set `GENTLE_DOCS` to that run's absolute `gentle_examples_docs` binary path
(`.exe` on Windows). Check the canonical generators without accepting drift:

```sh
"$GENTLE_DOCS" parity-matrix-check
"$GENTLE_DOCS" tutorial-catalog-check
"$GENTLE_DOCS" tutorial-manifest-check
"$GENTLE_DOCS" tutorial-check
```

Eight payload parity rows/counts were initially source-projected, not verified
by the Rust generator. Do not overwrite projections just to pass: retain any
generator diff as a primary failure. Reuse same-source LF/CRLF checkout receipts;
if missing, run the existing `scripts.test_tutorial_checkouts` gate once.

### Real CLI And Approved Adapter Replay

Set `GENTLE_PAYLOAD_CLI` to the absolute source-matched `gentle_cli` binary path.
Execute only the previously skipped process-test class, not completed Python
source checks:

```sh
GENTLE_PAYLOAD_CLI=/absolute/path/to/same-source/gentle_cli python3 -m unittest -v scripts.test_primer_payload_contract.PrimerPayloadCliAcceptance
```

Require its one actual-process test to execute both LF and CRLF cases, with no
skip. It seals a synthetic panel, traces, prepares, summarizes and verifies
outputs; it tests spaced/apostrophe paths, rehashed false-exclusion refusal,
source-byte preservation, absent state creation and report rename/reopen.
Retain a separate tutorial 10.04 Phase 2 replay's reports/sidecars and command
streams; helper PAF is invented teaching input, not an aligner result.

Also exercise the same synthetic requests through shared shell and approved MCP
`op`: prepare, summarize, inspect and request-bound cancel. Compare scientific
payloads after accounting for run-specific paths/bindings. The existing MCP
test proves refusal without confirmation, **not** a successful approved replay.
Require `confirm:true` even for generic MCP inspection; confirm missing approval
causes no files, project/history/state bytes stay unchanged, and inspection does
not implicitly reopen original FASTQ. Cancel uses the exact submitted request,
not a reconstructed or edited one. Use fresh destinations for every replay.

### Contract And Failure Boundaries

Verify expected versus foreign payload, differing primer lengths/orientations,
repeated read names and competing anchors without gene/barcode/name shortcuts.
Check actual Phase 1 bindings and full query accounting, explicit unmapped versus
unassessed rows, hashes/annotation refusals, MAPQ 0/255, cg/cs indels/skips/
contradictions and coarse-mode unavailable claims. Require masked all-anchor
opt-in, distinct read versus anchor counts, union rather than overlapping-HSP
sums, equivalent partners and non-transitive edge clusters. Primer plus one
foreign payload is not cross-locus; two qualified ordered payload loci may be
candidates, never a confirmed mechanism. Keep original source bytes unchanged.

The mask and exclusion details below are acceptance criteria for the same
payload filter/replay, not additional repeat runs. Cover the plan's remaining
matrix where not already exercised: every cap, rehashed structural tampering,
ID renaming/permutation, partial coverage, marker faults/cancellation at streaming
and batch boundaries, verified-unit retention, fresh-output races and Unix
nonregular/symlink refusals. Native Windows must exercise LF/CRLF, spaced/
apostrophe paths and closed handles before rename/publication/cleanup; a Unix
run or headless check is not Windows acceptance. Mark missing cases not_run,
not inferred pass from a narrower test.

### CLI Documentation Coverage, 2026-10-26

The reported 4535-pass run's sole failure was missing full glossary paths in
`docs/cli.md`. `320d5803` replaces the abbreviation with all four complete
usages; the Rust coverage guard is unchanged. Its focused command above must
pass. The LF/CRLF-aware Python equivalent already passed and need not be rerun
as a substitute for Rust. No parser/confirmation/freshness guard was relaxed.

### Primer Payload Mask Geometry, 2026-10-26

`e622f12c` shares clipped mask geometry between preparation and verification,
avoiding eager out-of-range subtraction and empty-denominator division. Require
12,168 abstract interval/orientation cases, empty/full masks, coordinate extremes
and structured refusals rather than panic/wrap. Four all-anchor query copies in
opposite read orientations must have agreeing normalized masks/sequence hashes;
bound import must work after original FASTQ removal without altering Phase 1.
Run with ordinary overflow checks, not a profile workaround. The original
all-anchor regression remains; the independent Python oracle is not Rust proof.

### Primer Payload Exclusion Bindings, 2026-10-26

`f1a60ae4` derives each selected anchor's disposition from the bound source row
and captured policy, not just an anchor-ID set. First accept genuine eligible,
masked, short, wholly unassessed and unattributed cases, including valid list
permutations. Then rehash only deliberate copies of the query manifest/FASTA:
hiding eligible anchors, changing policy, swapping exclusions between repeated
read records, wrong primer groups/reasons and missing/duplicate/foreign
observations must fail. Out-of-range anchors return errors, not panic. Two
one-record units must exceed an import record cap of one. Preserve the anchor
identity encoding. Public summary import must accept the genuine fixture before
refusing its tampered copy, with no final output or changed source bytes; the
actual-CLI companion includes this refusal.

### Separate Gates And Return

Separately authorize any live inner-agent exercise: discover prepare/summarize/
show/cancel, draft reviewed requests, require write approval, inspect/cancel
while execution slots are busy, preserve completed/cancelled/control-error
outcomes and refuse recursive agent execution. Deterministic tests are not a
native or live-model responsiveness verdict. Phase 2 adds no GUI; older native
RNA/tutorial checks remain open independently, not a reason to repeat completed
CLI checks. No private-data, throughput or scientific acceptance follows.

Return one compact report with per-group pass/fail/not_run/blocked, full tested
SHA/tree, lockfile/binary/input hashes, platform/toolchain/features/profile,
commands/exit codes/test counts/skips, primary failures versus cascades and
retained artifact paths plus a SHA-256 manifest. Identify reused CI receipts and
their exact source/scope. Return through Steffen to this main-maintaining chat;
do not push, dispatch CI or modify sealed historical evidence. No release,
molecule/isoform assignment or scientific approval is implied.

## Sequence-Design Approval Facet, 2026-10-26

Fold this into the pushing thread's single integrated request, not a separate
branch freeze or audit campaign. Claude reviewed the old proposal at `adacb2bc`;
Codex rechecked the permission/construction claims at `b544a21b` and the owner
authorized only slices 1--2. The owner subsequently requested local parity
execution. Keep the [source/feature-bound results](operation_lifecycle_policy_proposal.md#executed-local-verification-2026-10-26)
separate: `ab80743e` default-feature registry 21/21, matrix freshness 1/1,
canonical matrix/catalog/manifest OK and configured tutorial CLI companion 9/9;
post-rebase `98c9331a` default Cargo check/formatting plus headless
sequence-design 32/32 and S1 14/14. No ignored tests in these Rust filters.
The earlier default library-test build was interrupted to rebase, not failed.
No model, private-data or native GUI run was performed.

At the ordinary integrated-source verification:

```sh
cargo check -q --locked
cargo test --locked --lib sequence_design_approval_policy_is_request_aware_and_legacy_conservative
cargo test --locked --lib sequence_design_ordinary_planner_requires_report_write_confirmation
cargo test --locked --lib sequence_design_forged_flags_cannot_authorize_wrapped_report_writes
cargo test --locked --lib sequence_design_approved_write_attempt_can_publish_a_cancelled_report
cargo test --locked --lib sequence_design_approval_distinguishes_optional_paths_and_static_family_metadata
cargo test --locked --test capability_registry_parity sequence_design_registry_family_is_external_without_reclassifying_controls
cargo test --locked --lib sequence_design -- --test-threads=1
cargo test --locked --lib genome_edit_intent -- --test-threads=1
cargo run --locked --bin gentle_examples_docs -- parity-matrix-check
cargo run --locked --bin gentle_examples_docs -- tutorial-catalog-check
cargo run --locked --bin gentle_examples_docs -- tutorial-manifest-check
python3 -m unittest scripts.test_sequence_design_tutorial scripts.test_tutorial_checkouts -q
```

The five shared-engine approval regressions pass in the headless 32-test filter;
the sixth registry-family regression passes in the earlier desktop-enabled run.
Preserve the existing scientific/report/apply equality assertions: only the with-path
approval assertion in `assert_sequence_design_shell_parity` flips. Use the
existing synthetic 06.07 fixture with its provenance, not private DNA. Verify
ordinary builtin-echo planner candidates and falsely read-only Shell/Op/Workflow,
raw-op and macro wrappers cannot write without confirmation. Confirm no-path
requests do not remove explicit candidate/tool gates, empty/malformed supplied
paths are still writes, and approved cancellation retains a written report's
`cancelled` status. File publication does not establish successful search.

Static Plan operation/command metadata is External with confirmation, while
request policy is engine-owned and typed. SaveFile/S1/Apply and unmigrated
operations retain their classifications and runtime fallbacks. Both existing
raw sequence-design IDs are now registered, with two shell-only matrix
rows/counts and no new operation. Matrix
columns encode surfacing, not mutation; canonical byte-freshness now passes at
the recorded source without regeneration. Existing LF/CRLF checkout coverage
is retained and no new byte-exact freshness gate or fixture is introduced.
Maturity attestations stay unchanged.

Constructor/reconciliation/history/autosave behavior, project/report aliases
and writer overwrite/atomicity are deliberately unchanged. S1's decision now
states its apply-time and CLI/MCP disk scope rather than promising unchanged
JS/Lua state-out. The `ced03a22` in-process/MCP control ran in the headless S1
filter; native Windows, GUI/default-feature library/full-suite and scripting
construction coverage remain separate. Carry this evidence into ordinary
integrated verification, not another branch campaign. No S2 authorization follows.

## Five RNA Fixture Repairs, 2026-10-26

The owner-reported full library run had 4504 passes, five failures and 16
ignored tests; its tested SHA/toolchain was not supplied. The five affected
fixtures are repaired without changing production matching, source admission,
confirmation, marker handling or publication guards. Verify on the next ordinary
integrated source, not a new interim freeze:

```sh
cargo test --locked --lib -j1 rna_panel_opposite_record_orientations_cannot_hide_same_strand_duplicate_exons -- --test-threads=1
cargo test --locked --lib -j1 primer_trace_other_template_probe_is_not_contiguous_core_evidence -- --test-threads=1
cargo test --locked --lib -j1 primer_trace_marker_fault_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_twenty_loci_two_hundred_shared_isoforms_have_no_sixty_four_posting_bias -- --test-threads=1
cargo test --locked --lib -j1 glossary_ -- --test-threads=1
```

The overlap case uses existing invalid FASTQ bytes, not an absent path: preflight
must pass, then overlap refusal must precede parsing/progress/publication. The
200-template panel must recompute sharing diagnostics before its content digest;
strict validation and the 200-member/order-independence assertions remain.

The original GAT probe cannot establish boundary divergence because an allowed
terminal core deletion exposes GCC (two initial mismatches). Keep this negative
control and the AAA positive control, which disagrees at all three initial
bases at every admitted boundary. An independent Python edit-distance oracle
checks these literal controls, not Rust-Bio execution or assay performance.

Marker faults on the second unit must retain only the first EOF/hash-verified
unit and report a streaming control error for file_0001; faults before the first
verified EOF must retain none. Preserve the occupied-output sentinel/refusal.
Glossary rows use a numeric value, retain read-only inspection, and still reject
literal range text, counts over 1000 and positive rows without a file selector.

Run the existing RNA/primer-trace filters, all-features/headless checks and normal
freshness/CI as already consolidated below; do not rerun unrelated completed gel
or labware audits. All Rust/native checks remain unrun by Codex. Return actual
source/toolchain/platform, command counts, failures/skips and retained evidence;
no private-data run or scientific acceptance is implied.

Local checks: 166 Python/source checks pass, five optional binary checks skip;
LF/CRLF and the independent boundary oracle pass. Scoped rustfmt, documentation
status, maturity freshness and whitespace pass. None executes the Rust repairs.

## RNA Genomic Exon And Strand Admission, 2026-10-26

Reviewed against `11a6cd9b` (README-only successor to `4acee734`). The main
checkout's untracked `rna_panel_overlap_strand_codex_prompt_20261026.md` was
read, not copied or edited. No fetch, push or CI dispatch is performed here;
the existing local `f6777b6c` repair still requires integrated Rust/CI evidence.

At the next ordinary verification, run:

```sh
cargo check --locked --all-features --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_examples_docs
cargo test --locked --lib -j1 rna_panel -- --test-threads=1
cargo test --locked --lib -j1 rna_legacy -- --test-threads=1
cargo test --locked --lib -j1 agent_rna_ -- --test-threads=1
```

Eleven new `overlap_tests` cases are authored, unrun. Check that reflected
records cannot hide same-genomic-strand exonic duplication, while opposite-
strand exonic overlaps are admitted independently of extraction orientation.
Keep exact assembly/chromosome comparison, adjacent one-base boundaries,
disjoint exons despite shared flanks, and unknown/unmappable conservative
fallbacks. Same-record/DNA-hash/declared-gene safeguards remain active even
for antisense records. Present assessments must be recomputed and reject
honestly rehashed tampering; absent historical fields retain old wire bytes,
hashes and admission, not an implicit upgrade. New builds/extensions get new
panel identities; do not pool their reports with historical-panel results.

Tutorial 10.04 now explains `source_overlap_assessment.rule=genomic_exons_v1`,
`antisense_overlapping_loci` and `conservative_source_pairs`. Inspect new
panels and summary JSON through GUI, shared shell/MCP and the live inner agent:
discover `source_overlap_contract`, preserve approved writes, distinguish
annotation strand from read/library orientation, and never collapse or reassign
ambiguous counts. Summaries retain source geometry even with zero verified
files; their inspection is structural, not fresh source authentication/mapping.
TSV counts and scientific nomination policy are unchanged.

Run canonical tutorial/parity checks and the existing no-shell/native Windows
path/handle acceptance at the actual integrated source. No live-model/native
or Rust pass is claimed here. Optional random-background characterization was
not added: this repair changes admission only, not sensitivity or budgets.

Offline verification: `python3 -m unittest scripts.test_rna_panel_tutorial
scripts.test_tutorial_checkouts scripts.test_capability_maturity` passes all
115 tests. Scoped rustfmt, documentation status, maturity freshness and
whitespace checks pass. Session-close reports four OK, two expected
worktree/manual-plan warnings and no failures; roadmap remains 400 lines.
These checks do not execute the new Rust regressions or establish native parity.

## RNA Work Profiles, Error Scripts And Summary Evidence, 2026-10-26

Use normal integrated source, not a new candidate freeze. The residual prompt
was read from the main checkout's untracked
`docs/rna_panel_residual_gaps_codex_prompt_20261025.md`; that draft is not copied
or changed here. Baseline local main was `adacb2bc`. The local `f6777b6c`
owned-ID repair is present; read-only search finds no explicit `GuiTestId`
callers outside its wrapper depending on `Copy` or a static `as_str` lifetime.
This is not compile/CI evidence, and no push/dispatch is performed here.

At the next ordinary verification:

```sh
cargo check --locked --all-features --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_examples_docs
cargo test --locked --lib -j1 rna_panel -- --test-threads=1
cargo test --locked --lib -j1 rna_legacy -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_eight_percent_acgt_edits -- --nocapture --test-threads=1
```

New `index_work_profile` measures maximum postings/templates/loci per word and
fixed fan-out counts. Each completed unit's `work_profile` retains total/peak
anchors, fixed histogram, peak record/ID, peak template-orientation chain calls
and file-local key bytes debited to the run budget. Long peak IDs are bounded,
explicitly truncated and full-ID hashed. No wall time is digest-bound, and
unfinished units publish no work claims. The no-profiling private test path
must produce byte-identical nomination TSV/FASTA and class/member counts.
Legacy omitted profiles retain old semantic hashes; panel digests are unchanged.

Summary JSON/TSV now retains exact sequence-sharing group records/citations
and a common-sequence flag, not a causal explanation or source-identity verdict.
Unknown historical evidence stays unknown, counts stay unchanged and no-source-
file-reverification claims remain. Check genuine differing-sequence cross-gene
sets, exact common sequences, unknown old panels, zero verified files and
re-digested tampering. Verify canonical generator/parity, inner-agent/MCP and
native projections in tutorial 10.04, including Windows cleanup/handles.

The deterministic substitution/short-indel scripts use actual A/C/G/T, not
skipped N windows. The 2%/8% software characterizations do not calibrate real
Nanopore sensitivity. Rust execution/native acceptance is not performed here.
On the separately authorized first complete chunk, retain profiles and external
time/memory, selected counts, panel/policy, exact SHA/toolchain and raw budget
failures before considering any index rewrite. Keep per-template equivalence,
hard budget failures and zero-length-record refusal unchanged.

Owner question, not an implementation change: does the facility/basecaller
actually emit zero-length FASTQ records? If yes, decide separately whether
they remain fatal or deserve an explicit digest-bound intake class. No such
class, scheduler, private-data run or paused-goal resumption is introduced.

Local verification: 114 tests pass with `python3 -m unittest
scripts.test_rna_panel_tutorial scripts.test_tutorial_checkouts
scripts.test_capability_maturity`; scoped rustfmt, documentation status,
maturity freshness and whitespace pass. Session-close: four OK, two expected
working-tree/manual-plan warnings, no failures; roadmap remains 400 lines.
Nineteen new Rust tests in `profile_tests`, `error_tests` and `sharing_tests`
are authored but unrun. Preserve their synthetic provenance and retain the
8% characterization output without converting it into a sensitivity verdict.

## Primer Trace Phase 0/1, 2026-10-26

Integrate this into the pushing thread's **single** candidate-bound request;
do not start a separate per-branch acceptance campaign. Implementation started
at `adacb2bc`; use the final integrated source/binary/lock/toolchain/profile in
your return. The original untracked primary-checkout Claude prompt is not
acceptance. Its [revisions](primer_trace_phase01_plan_20261026.md) and
[contract](protocol.md#primer-anchored-read-trace-phase-01) govern this slice.

No local Rust execution or private-read analysis was performed by Codex.
Execute, recording exact filters, discovered counts, failures and skips:

```sh
cargo check --locked --all-features --tests
cargo test --locked --lib -j1 primer_trace -- --test-threads=1
cargo test --locked --lib -j1 rna_panel -- --test-threads=1
cargo test --locked --lib -j1 rna_legacy_ -- --test-threads=1
cargo test --locked --lib -j1 glossary_cli_usage -- --test-threads=1
cargo run --locked --bin gentle_examples_docs -- parity-matrix-check
cargo run --locked --bin gentle_examples_docs -- tutorial-catalog-check
cargo run --locked --bin gentle_examples_docs -- tutorial-manifest-check
python3 -m unittest scripts.test_primer_trace_contract scripts.test_rna_panel_tutorial scripts.test_tutorial_checkouts -q
```

The primer-trace implementation is rebased onto local main `11a6cd9b`, including
`4acee734`'s separate RNA work-profile and summary-evidence changes. Verify the
integrated source, not the pre-rebase `adacb2bc` baseline. The `primer_trace`
Rust filter now includes typed-discovery, managed-worker and agent-feedback
regressions (25 new tests total, authored but not locally executed).
Full and compact results must preserve
cancelled/control-error outcomes in agent feedback; inspecting an old cancelled
report must still be a successfully completed inspection.

The typed-discovery regression remains unrun: all seven operation/shell
descriptor examples from the actual shared-shell introspection result must
deserialize through the canonical types/parser, preserve approval classification,
and leave sequence/digest/checksums/cancel path and completion unset. Inline
examples must be refused even by explicit sealing. Check the actual inner-agent capability
projection without treating example orientation/defaults as reviewed inputs;
generic MCP confirmation remains required. No model call, private-read use or
new GUI is necessary for the deterministic discovery regression. The cancellation
example must refuse marker publication, not merely fail after writing a dummy path.

Check the actual shell/MCP/managed-host boundaries, not just the matching
helper: writes and `--cancel` require approval, forged plan metadata cannot
remove that gate, validation/show leave state/journal alone, and the MCP
project file is byte-unchanged. Managed standalone trace keeps a typed
`primer_trace_result` without merging a live project journal. Test project
replacement, close, source-file changes, late host Cancel, explicit marker
Cancel, publication failure and reopening the bound report. Preserve completed,
cancelled and control-error outcomes; late host cancellation is not a new
scientific verdict. No GUI acceptance is claimed or needed for this slice.
The new managed tests use channels to hold the interval after publication, force
compact output, replace the live project and issue late Cancel; the published
completed binding must survive without any project commit. A separate test holds
the project write lock throughout partial-cancellation/control-error execution
and reopens the retained verified unit through its exact binding. These tests
are synthetic software checks, not native, platform or private-data acceptance.

Synthetic regressions cover four exact orientation geometries, substitutions/
indels, core limits, quality/nulls, shared handles, multiple anchors, template
ambiguity, EOF digests, artifacts, inspection binding and LF/CRLF intake.
Keep all retained alternatives and distinguish a payload-only template match
from contiguous core/probe evidence. The new TSV has no read sequence.
Confirm existing RNA-screen/report bytes and legacy nomination regressions
remain unchanged; the intake's added quality field is not serialized there.
Canonical matrix rows/counts were reconciled from registry rules without
running the Rust generator; its exact-byte check remains required. Generated
maturity hashes carry no new maturity attestations.

Only after separate authorization, select exact retained completed inputs
and caller-reviewed **complete ordered oligos**, not reconstructed September
core proposals. Confirm which oligos were present at RT/PCR and whether the
kit supplied a shared-handle oligo(dT) primer. Kit identity was not verified.
The `handle_followed_by_poly_t` flag remains a configurable sequence observation.

For authorized barcode13/14 analysis retain:

- exact full-hit counts separately by primer and read strand, with both read
  and hit denominators; do not equate the two or interpret them as enrichment;
- per-position start bins around 81..83 and remaining-tail bins around 25..29,
  without calibrating thresholds just to reproduce those observations;
- every expected/other/unresolved relation, all template alternatives, core
  terminal mismatch counts, and the earlier short off-locus nominations;
- cancellation/checksum/resource failures separately from observed zero hits;
- barcode15/16 controls when actually available, never invented local controls.

Report whether defaults are useful, not merely whether execution completes.
Return exact source/binary SHA-256, effective request and panel digests,
commands, platform/profile, typed outcome, denominator interpretation and
separate technical/scientific verdicts. Share private read-level evidence only
through the approved channel; no automatic publication is authorized.

Owner decisions remain: actual ordered oligos and RT/PCR/kit roles; acceptance
of uncalibrated defaults/axis model. FASTA stays off unless explicitly requested.
Phases 2--5 (segmentation, breakpoints, wider partners and later interpretation)
remain pending separate contracts; no precise unseen phase specification is
inferred from the external discussion draft.

## Shared-Sequence Ambiguity And Legacy Output Regressions, 2026-10-25

At ordinary integrated source, run:

```sh
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 rna_legacy_ -- --test-threads=1
cargo test --locked --lib -j1 multi_gene_sparse -- --test-threads=1
cargo test --locked --lib -j1 rna_sparse_target_gene_presence_ -- --test-threads=1
cargo check -q --locked --tests
cargo check -q --locked --all-features --tests
```

New panels bind optional `shared_transcript_sequences` over selected templates.
Exact mature-sequence equality across independent loci is a diagnostic, not an
overlap verdict. Review unknown source relations rather than merging/rejecting
genes by matching RNA. The synthetic tests cover differing source extents with
repeated declared gene identity (refused), distinct gene IDs with a reused label
and identical mature RNA/different exon boundaries (retained as two candidates,
one cross-gene set/read), deterministic order, honest re-digestion tampering and
legacy omission/nomination equivalence. Current-source `panel-show` must expose
the same groups as the native preview/Inspect and inner-agent/MCP shared route.
The two `rna_legacy_` tests compare entire scientific reports with only creation
time and documented mode notices normalized; candidate scores, supports,
histograms and other warnings are retained. These are regression comparisons,
not independently executed historical byte baselines or biological acceptance.

Follow tutorial 10.04's source-review note in the native/agent rehearsal.
New explicit builds change panel identities; absent legacy diagnostics remain
unrecorded without rewriting old panels. Existing budget failures, conservative
whole-source overlap guards and marker/cancellation outcomes are unchanged.
Rust/native execution and real-library calibration remain Glen/CI work; the
previous unrelated goal stays paused and no private inputs or push are involved.
Local checks: `python3 -m unittest scripts.test_rna_panel_tutorial
scripts.test_tutorial_checkouts scripts.test_capability_maturity` passes 113
tests; scoped rustfmt, maturity freshness and whitespace also pass. These are
source/checkout checks, not execution of the new Rust comparisons or GUI.

## Focused macOS Labware Evidence, 2026-10-25

The owner forwarded Claude's completed review/test report for exactly
`39ac4e4425d02d5374181a1c2a5e0fb9ae35d884`, including disclosure cleanup
`ab692041` and derived storage box `c217d0b0`. Reported results on macOS:

| Scope | Reported result |
| --- | --- |
| Protocol crate labware tests | 14 passed |
| Render crate labware tests | 4 passed |
| Root labware/rack/agent-parity/mutation/glossary filters | 289 passed |
| Python/OpenSCAD checks | 11 passed, none skipped |

Claude explicitly reports actual OpenSCAD rendering, including storage-box
mesh bounds and cut strips. These checks render repository templates with
literal parameters, not files exported by GENtle. They do not establish
physical fit, cap/lid clearance, material/temperature suitability or native GUI
acceptance. The earlier Codex-host OpenSCAD/Qt processor failure remains a
separate environment observation, not a failure of Claude's reported run.

The supplied summary does not include raw logs, exact command transcripts,
toolchain/profile or OpenSCAD version. Preserve those separately if available;
this is an attributed report, not an independently rerun or manifest-verified
receipt. It does not certify later source or close release/scientific gates.
The strip-test repair is no longer awaiting its first reported focused rerun.

Remaining Glen/CI work: normal CI and full `--lib`, including Windows; the nine
provenance-bound legacy carrier baseline files; OpenSCAD renders of fresh
source-matched GENtle exports for carriers/strips/boxes; and native GUI review.
Do not repeat the reported focused campaign solely because older sections below
still request it. Recheck affected boundaries when their source changes. Square
SVG cells and measured commercial-box profiles remain optional separate work.

## Historical RNA Review Follow-Up, 2026-10-25

Run focused `rna_panel_`, `rna_sparse_target_gene_presence_`, `rna_read` and
`multi_gene_sparse` library tests with `--locked -j1 -- --test-threads=1`, then
normal CI and `cargo check -q --locked --tests`. New tests cover exact source
DNA duplicates despite distinct annotations/no anchors, optional-field legacy
bytes, restricted competitors, default anchor exhaustion from 64 shared isoforms,
checked cumulative key-byte accounting, complete gzip-member loss and same-locus
opposite orientations. The key-budget unit test uses synthetic debits, not a
thousand-file end-to-end run. Existing noisy-read and publication/control tests
remain. Test aliases/current-source generator/parity and Windows handles too.

The attached review at `b1bae75e` is historical: tutorial restriction/sharding
guidance and several intake/error tests already existed. This follow-up does
not replace the per-template index, relax caps, add a resource-limited read class
or modify the legacy winner selector. Optional source-DNA hashes change newly
built panel identities, not old omitted-field digests. Rebuild deliberately;
missing historical duplicate coverage is unknown, not proof of uniqueness.

Before changing indexing/budgets, use the owner's separately authorized first
completed chunk with all annotated isoforms. Retain indexed postings, selected/
annotated counts, actual peak posting work (new report profiles above need
execution/measurement), key bytes, class counts, toolchain/SHA and memory/throughput evidence.
The authored default-cap regression exercises a legitimate synthetic shared
read exceeding the budget; it is unrun and establishes no real-library risk
frequency or performance. Require
per-template candidate equivalence before adopting a single-chain locus index.
No private data or Rust/native command was run locally for this follow-up.
Local source/checkout verification: 111 Python checks pass, as do scoped
rustfmt, maturity freshness and whitespace. This is not execution of the new
Rust tests. The earlier unrelated goal remains paused; `outputs/` is untouched.

## Derived Square-Cell Storage Box, 2026-10-25

The [reported focused macOS pass](#focused-macos-labware-evidence-2026-10-25)
supersedes the initial protocol/render/root-filter execution request at
`39ac4e44`; full `--lib`, normal CI/Windows and source-matched export/native
checks remain open. Retain `cargo check -q --locked --tests` in normal CI.
Authored storage-box tests include:
`labware_storage_box_geometry_has_square_cells_variable_grids_and_separate_lid`,
`labware_storage_box_refuses_invalid_grids_unknown_caps_and_colliding_collars`,
`labware_storage_box_exports_variable_grids_square_cells_and_all_disclosures`,
`labware_storage_box_refuses_compact_pitch_before_any_export_write`.
Extended carrier shell/generic-op, GUI label and legacy metadata guards remain.

Require 9x9 dimensions 104.95 x 104.95 x 47 mm and A1/I9 centres, plus 4x12
139 x 48.2 x 47 mm, complete physical cells despite blocked positions, the
separate lid and clear override refusals. The source-derived Python mesh checks
are not Rust exporter execution. A long synthetic cone propagates the existing
`collar_not_above_top` derivation refusal, not an invented cone-only error.

Generate source-matched 9x9 and 4x12 exports; render OpenSCAD box/insert/lid and
inspect all four SVG disclosures. Verify the native dropdown, round-marker
warning, variable grids and compact-override error. Commit the nine historical
carrier baselines from `1c5be317` with provenance/LF-CRLF coverage and compare
older spec/SCAD/simulation bytes. No existing output is intentionally changed;
byte identity, physical fit, cap/lid clearance, material/freezer/ice suitability,
Windows/native interaction and release acceptance are not established here.
Local verification: 19 focused source/maturity checks pass and existing SCAD
source-template projections match the pre-box source for collar/floor/PCR96.
On the earlier Codex host, `python3 -m unittest scripts.test_labware_models -v`
was attempted: installed OpenSCAD/Qt aborted before parsing with
`Incompatible processor ... neon`, also
for unchanged tube/carrier/strip models. That attempt did not accept any mesh
method. Claude subsequently reports all 11 checks passing on a compatible
macOS environment; GENtle-exported-file rendering remains separate. The four
maturity JSON projections changed only because the glossary now lists the
additive template choice.

## Compact Carrier Disclosure Construction, 2026-10-25

The optional cleanup builds disclosure vectors directly, preserving their
strings/order rather than rewriting `nonclaims[2]` or `assumptions[4]`. Run
`cargo test --locked --lib -j1 labware_compact_pitch_exports_geometry_persists_and_clears_without_drift -- --test-threads=1`
and the normal `labware_`/legacy-export checks. The extended assertion compares
the two common nonclaims against the conservative spec; existing clear/reset
SCAD equality remains. Rust/export-byte acceptance is unrun locally. No template
numbers, mutation policy or release scope changed.

## S1 SNV State-Preservation Follow-Up, 2026-10-25

At the next ordinary integrated verification, include the named
`genome_edit_intent_preserves_unreconciled_state_on_success_and_refusal` regression
and existing CLI/MCP boundary tests:

```sh
cargo test --locked --lib -j1 genome_edit_intent -- --test-threads=1
cargo test --locked --bin gentle_cli -j1 genome_edit_intent -- --test-threads=1
cargo check -q --locked --tests
```

The new fixture
demonstrates that each reconcile routine independently changes state, then
requires unchanged state/revisions/history on successful and refused SNV
inspection for both orientations and inline/loaded input. The case-only test
also distinguishes uppercase reported alleles from retained inline spelling.
These Rust changes have not been executed locally. Return the actual source,
platform/toolchain, executed counts and primary failures through Steffen; do not
reinterpret earlier S1 results as this patch's acceptance. No new capture,
release gate or standalone audit campaign. The [broader lifecycle-policy review](base_prime_editing_preparation.md#landed-s1-review-reconciliation-2026-10-25)
must precede S2 and must not grant file/confirmation exemptions from a broad
`ReadOnly` registry label. No editor feasibility or laboratory work is added.

## RNA Glossary And Bound-Rehearsal Fixtures, 2026-10-25

The owner's four library failures expose stale fixtures after result binding:
the glossary supplied untyped `demo` for `--binding` and omitted draft mode for
`--save`; MCP used a macOS temporary-path alias instead of the receipt path;
two summaries compared hashes across changed output/order request identities.
Tests now retain strict parser/location checks, compute independent expected
request/content hashes and compare every result field rather than masking hashes.
MCP/direct table verification uses the same bound summary; forged content must
fail on content binding. Production, historical digests and approvals are unchanged.

At the next ordinary integrated run, include `glossary_`, `rna_panel_summary_`
and `rna_panel_mcp_rehearsal_matches_direct_scientific_payload`, with locked
single-job library tests and `--test-threads=1`. Retain the macOS path-bound case,
LF/CRLF reopening, no-input-reread, tampering and incompatible-policy assertions.
Do not repeat completed campaigns or freeze development; return tested source,
platform/toolchain and primary failures through Steffen. Source-only checks do
not establish Rust/runtime acceptance. Continue the native rehearsal below only
with a verified source-matched binary.

## RNA Draft Admission And Cancellation Boundaries, 2026-10-25

At ordinary integrated main, verify deferred metadata-only `panel-screen
--request @FILE` loading. Native **Load request...** now uses parsed independent
observation admission, without execution capacity or a live-project baseline.
Read the file once on the worker (128 MB source cap; interactive result budget
still 16 MiB). Keep draft saving approval-bound. No private inputs, local Rust,
native/model execution or push was performed for this follow-up.

Run these authored regressions, then the normal RNA/command-service/agent
approval/parity suites and `cargo check -q --locked --tests`:

```sh
cargo test --locked --lib -j1 rna_screen_draft_load_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_screen_draft_ -- --test-threads=1
cargo test --locked --lib -j1 agent_rna_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_service_late_cancel_ -- --test-threads=1
cargo test --locked --lib -j1 replacing_gui_project_cancels_ -- --test-threads=1
```

The draft test fills all four execution slots and holds the live write lock
while invoking the real native Load handler; require loaded metadata, zero live
journal changes and no opening of referenced inputs. Test lazy missing/invalid/
oversized files, LF/CRLF/shebang JSON and quoted Unicode paths on Windows too.
The retry regression must preserve the current draft/review after missing or
invalid input, then clear the old error after a valid load without changing the
source file, live state/history or outputs. Loaded pins still require fresh review.
The agent tests hold a synthetic screen after its first verified file. Explicit
pending-command Cancel must retain exactly that unit with a bound cancelled
report, both with and without a marker path. Model Stop and project replacement
must leave both marker-backed and markerless artifact-only work running, while
the existing project-bound-work test still cancels its retired owner. Preserve
the held publication-to-retention late-cancel test: an already published
outcome cannot be relabelled.

Follow tutorial 10.04's **Which Stop Or Cancel Do I Use?** in native and inner-
agent acceptance. Runs marker-control failure stays visible with no implicit
service-cancel fallback; engine pre-publication cancellation remains valid.
No numeric cross-process job ID, worker restart recovery, queueing, automatic
result attachment or biological verdict is implied. Report actual SHA/platform,
commands and evidence via the main-maintaining thread, not a new source freeze.

## Native RNA Compile Repairs, 2026-10-25

The owner's build exposed two primary integration defects: missing explicit
paths for the native processing/evidence modules, and a numeric journal count
shadowing the optional command-admission baseline. The missing UI-method errors
were cascades of the module failure. Both are repaired without disabling UI
handlers or weakening project-commit/stale/cancel safeguards. Two new LF/CRLF
source regressions fail before repair and pass afterward; they are not Rust
verification. Include this repair at the next ordinary integrated source.

Before the native rehearsal, run the default `cargo check -q --locked --tests`
and headless check below, then the existing `command_execution`, `rna_screen_`
and `rna_native_` filters. Retain replaced-project refusal, retained-operation
counts, artifact-only no-commit behavior and the existing native handlers.
Report actual source/toolchain/platform, primary failures and cascades through
Steffen. Do not repeat completed campaigns or impose a new source freeze.

## Native RNA Processing And Tutorial 10.04, 2026-10-24

Integrate `4603d029` (execution/bindings), `a87fb9c6` (setup/monitor),
`6d7d6ef3` (inspector/export) and the subsequent rehearsal/parity commit at
ordinary current main; do not freeze another candidate. The earlier unrelated
goal remains paused. No local Codex Rust/native run, private sequencing access,
capture, provider call or push is claimed.

Run focused checks, then normal full CI on that exact integrated SHA:

```sh
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 rna_screen_ -- --test-threads=1
cargo test --locked --lib -j1 rna_native_ -- --test-threads=1
cargo test --locked --lib -j1 command_execution -- --test-threads=1
cargo test --locked --lib -j1 rna_artifact_writes -- --test-threads=1
cargo test --locked --test capability_registry_parity --test parity_matrix_freshness
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_examples_docs
cargo test --locked --lib --features gui-test-support -j1 rna_screen_native_semantic_targets -- --test-threads=1
cargo fmt --all -- --check
python3 -m unittest scripts.test_rna_panel_tutorial scripts.test_tutorial_checkouts scripts.test_capability_maturity
```

The `rna_screen_` filter includes draft-save agent confirmation. Preserve
existing RNA MCP, GUI intent/palette, approval and command-service contracts.
The authored app test holds the live project write lock across native screen,
monitor, bound reopening and summary export; require two synthetic samples,
four reads each, retained class/member sets and no live journal append. This
is not native GUI evidence. Keep the channel-controlled late-cancel/concurrent-
job tests below and ordinary mutator stale rejection.

With a verified source-matched docs binary, check catalog/manifest/parity
generation and tutorial replay using the repo's normal commands. The five
new native command cells and GUI totals were projected from the changed
glossary; exact Rust generator freshness still needs execution. Maturity
projections were refreshed by their Python generator. The synthetic demo now
contains stable toy `gene_id` qualifiers and a stored resolved gene set so
native preparation uses the real binding path. New demo/panel/report hashes
legitimately change; legacy omitted-field hashes must not.

Follow [the GUI-only rehearsal](tutorial/10-04_targeted_nanopore_fastq_panel.md#gui-only-rehearsal)
and its [typed checkpoints](tutorial/10-04_targeted_nanopore_fastq_panel.md#typed-checkpoints-and-capture-targets).
Retain permitted synthetic source/binary/snapshot-bound screenshots of setup,
progress/outcome, nomination prefix, summary ambiguity and table verification.
Test native and hosted narrow-window scrolling. Missing dialog/control automation
coverage is `harness_gap`, not guessed clicks or a biology verdict. This chapter
remains manual/hybrid, not a new generic acceptance-runner contract.

Use an observable source-bound synthetic job for scheduling checks: other windows
and agent conversation stay usable; close/reopen does not cancel; four occupied
execution slots do not block marker control/plain inspection; table verification
uses bounded execution capacity. Verify cancel before/during/after publication,
marker-write failure and post-start control_error retention. Reopen saved bindings
after restart without claiming worker recovery. On Windows test handles, quoted/
Unicode paths, raw traversal, fresh atomic draft save/no overwrite and LF/CRLF
request/cancel identity. Do not use private downloads to slow down the rehearsal.

Inner/outer-agent rehearsal covers `ui open rna-read-panel`, discovery/build,
metadata-only `panel-screen --request ... [--save ...]`, approved screening,
`introspect runtime`, request-bound cancel, bound `panel-show`, summary export
and table verification. Writes stay `ask`; draft saving via MCP uses
`agent_execute_plan`, scientific operations remain `op`/workflow/script accessible.
Compare canonical scientific payloads, not host-local job IDs or timestamps.
No model call was made here. Return SHA/toolchain/platform, exact commands,
ledger/hash locations and classified failures via the main-maintaining thread.
Queueing and short result-graph attachment stay separate.

Exercise existing size budgets explicitly: interactive result metadata is capped
at 16 MiB. An oversized artifact run retains its bound compact receipt/full
disk report; an oversized inspection can refuse rather than silently truncate.
Measure actual metadata size for the intended panel/file counts before claiming
large-library GUI coverage. Paging beyond that budget is not implemented here.

## File-Only RNA Execution And Result Bindings, 2026-10-24

At the ordinary integrated source, run the focused `rna_panel_` lib filter,
existing `command_execution::tests`, `agent_feedback`/approval/MCP parity tests,
glossary/capability/parity freshness and `cargo check -q --locked --tests`.
No new freeze, local Codex Rust build, private FASTQ run or native acceptance is
claimed. Regenerate/replay the synthetic tutorial with a source-matched binary;
new submitted-request fields legitimately change new report hashes, while
legacy omitted-field canonical bytes and cancel identities must remain stable.

New engine tests cover exact frozen legacy request bytes (LF/CRLF), reviewed
panel replacement before streaming, submitted/effective manifest normalization,
changed/deleted source files after publication, binding tampering, bound shared-
shell inspection and summary/legacy metadata. New service tests use channels to
hold the publication-to-retention window, request late cancellation and edit
live structure; require retained Completed evidence, zero committed operations
and no overwritten project edits. A separate held mutator must still commit
after artifact retention, proving there was no intervening live journal append.

Exercise callback cancellation and foreign-marker control errors after the
first verified file. Require Cancelled/Failed lifecycle states respectively,
retained bound partial reports and unchanged scientific outcomes. Fill all four
execution slots and hold the live project lock: marker control and metadata
inspection must still finish in their independent lanes; explicit expensive
verification must report execution-capacity refusal. Registered panel builds
remain normal project commits with correct mutation feedback. Oversized file-
only output must retain a small bound envelope and the full published report.

Preserve existing mutator stale/cancel rejection and bounded history tests;
do not infer file-only policy from `state_changed=false` or caller JSON. Review
`submit_parsed` freezing versus raw-text worker parsing, agent project-change
handling and compact-result feedback. Cross-process cancellation remains the
request-bound marker; these process-local job IDs are not portable handles.
Native setup/monitor/inspector/export and the rehearsal are implemented;
integrated runtime/native acceptance is described above and in the
[processing plan](rna_panel_gui_processing_plan.md).

## Gene-Set RNA Panel Preparation And Graph, 2026-10-24

At the ordinary integrated source, without another freeze, run:

```sh
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_mcp -- --test-threads=1
cargo test --locked -p gentle-shell ui_intent
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo check -q --locked --tests
```

Run the normal tutorial/catalog/maturity generators and freshness checks; the
source/metadata Python tests do not replace Rust execution. New authored tests
cover loaded-group discovery, immutable 10-to-20-locus extension, source-hash
drift, gene-set binding, registration/save/reload/undo/redo, GUI/SVG graph edges,
UI-intent lifecycle, stale workers and MCP read-only discovery. Require legacy
omitted-field digests and the synthetic screen/summary outputs to stay unchanged.

Follow tutorial 10.04 in a disposable native and hosted project: resolve an
explicit set in Gene Set Inspector; import its annotated loci; refresh, bind exact
stable IDs, review provenance/transcripts and build. Verify the set/source/panel
graph, reopen by right-click and double-click, save/restart, then extend with an
expanded set. Undo must remove only registration, not the external file; old
panels/screens must remain unchanged. A moved panel file must not prevent
snapshot inspection, but must prevent extension until explicitly resolved.
Missing/ambiguous IDs, unresolved members and stale source hashes are errors,
not zero transcripts. Other windows must remain usable during preparation.

Compare CLI `panel-build --list-sources [GENE_SET_ID]`, the copied build JSON,
MCP confirmed `op`, and inner-agent approved builds using the same scientific
payload. Discovery is read-only; file/project writes stay ask/confirmation-bound.
Live-model invocation and native semantics remain unverified locally. This is
panel preparation, not a new native FASTQ-screening form, isoform-discrimination
acceptance or a private sequencing run. Return exact source/platform/toolchain.

## PCR-Strip Test-Contract Repairs, 2026-10-23

Glen's two isolated failures from `fc75114a` are corrected in tests only:
16 successful exports append 16 external-operation receipts without changing
project state; invalid count types are refused during operation execution, not
generic shell parsing. Refusals must preserve journal length and existing model
bytes. Valid exports, quoted paths, legacy serialization and capability checks
remain covered. At the next ordinary integrated source, rerun only:

```sh
cargo test --locked --lib -j1 labware_pcr_strip_ -- --test-threads=1
```

Return actual source/toolchain/platform and results through Steffen. No repeated
mesh render, completed acceptance campaign or new source freeze is requested.
Codex ran no Rust verification; production journaling/parser behavior is unchanged.

## Mid-Run FASTQ Control-Error Retention, 2026-10-23

At the ordinary integrated SHA, run the focused `rna_panel_` and
`rna_panel_mcp` lib filters, agent-feedback/agent-execution and shared-shell
approval tests, existing glossary/parity/catalog freshness, and
`cargo check -q --locked --tests`. No new freeze or private chunk is requested.

New `rna_panel_mid_run_marker_errors_preserve_verified_units_with_explicit_outcomes`
injects faults through progress callbacks without sleeps. Require preserved
first-unit bytes/hashes and `status=control_error` at file boundary, during
intake and at finalization; malformed/nonregular markers follow the same policy.
Check zero verified units remain unassessed, matching cancellation remains
cancelled, and a later callback cancel cannot overwrite the control error.
The hard checksum/resource-error regression must still publish no result, even
after an earlier unit verified. Preflight mismatch remains a hard refusal.
The binding regression rejects missing diagnostics, mismatched paths and
relabelled success/cancel/legacy outcomes. MCP rehearsal now inspects the
retained control-error report and exports a partial summary with the same cause.

Replay tutorial 10.04 in a disposable project through CLI, MCP and the live
inner agent. Require a visible warning, failed screening feedback and retained
named missing inputs, while named inspection can complete. A summary keeps
`source_runs[].control_error` and remains partial even if finalization failed
after every file verified. The agent should explain fresh control/output paths
and explicit remaining-file retry, not silently restart or edit the marker.
Native Windows must verify closed handles before partial-unit cleanup and
publication; retain Unix nonregular/symlink coverage. Optional diagnostics are
omitted when absent so historical v1/v2 semantic digests remain unchanged.
No Rust/native checks were run locally; report exact source/toolchain/platform
and executed results, not a new biological or throughput verdict. The added
selection/cumulative-budget tutorial cautions do not implement the separate
index rewrite, resource-limited continuation or legacy-path follow-ups.

## Stored-Plan And RNA MCP Test Expectations, 2026-10-23

The owner's two latest library failures were obsolete approval expectations,
not evidence that either production gate should be relaxed. The shell test now
requires denial for a parameter write (including falsely read-only metadata),
unchanged denied state/history and one confirmed operation. A separate read-only
shell plan still executes without confirmation. The RNA MCP rehearsal first
denies raw `op` verification without `confirm=true`, checks untouched project/
report bytes, then explicitly approves it and retains scientific-payload parity
and unchanged semantic state. Shared-shell read-only exemptions are unchanged.

At the next ordinary integrated run, include these focused checks; do not
restart completed acceptance campaigns or impose an interim source freeze:

```sh
cargo test --locked --lib -j1 engine_shell::tests::execute_agents_execute_plan_ -- --test-threads=1
cargo test --locked --lib -j1 mcp_server::tests::rna_panel_mcp_rehearsal_matches_direct_scientific_payload -- --exact --test-threads=1
```

Use the normal run-owned target and return actual source/toolchain/platform,
failures/skips and retained logs through Steffen. Codex ran no Rust verification.

## RNA Table-Verifier Import Repair, 2026-10-23

The owner's compiler found E0432 in `rna_panel/verification.rs`: `sha2` was
not declared by the root crate. The verifier now uses the same already-declared
`ring::digest::Context` / SHA-256 as FASTQ intake, keeping bounded streaming and
all exact hash/size/path checks. No dependency or lockfile update is needed.
The previous repair's two unused `Engine` imports are removed. Include the
existing `rna_panel_summary_verification_is_explicit_bounded_and_hash_bound`
Rust regression in the ordinary `rna_panel_` filter below; source/format checks
are not compiler acceptance. No new campaign or repeated historical checks.

## RNA Confirmation, Starter And Tutorial Freshness Repairs, 2026-10-23

The reported five failures have three causes: the stored-plan executor trusted
shell candidate flags instead of shared command mutation policy; the demo saved
unreconciled lineage/containers; and 10.04 lacked a review entry while its catalog
projection had noncanonical indentation. The repair checks parsed commands for
every candidate kind, saves an engine-reconciled synthetic starter and registers
10.04 explicitly as unreviewed. It preserves V2 co-nomination, cancellation and
table-verification contracts from `874ff657`.

Do not relax confirmation, unchanged-state assertions or freshness checks.
New Rust tests exercise mislabeled shell/raw-op writes denied before artifacts,
successful explicit approval, unchanged blocked state/history and read-only
inspection without approval. The MCP rehearsal now also checks that loading
the saved starter is idempotent before its unchanged-state comparison.

At the next ordinary owner-pushed source, run the existing RNA/MCP checks below
plus `cargo test --locked --lib -j1 agent_execution::tests:: -- --test-threads=1`.
Re-run the three originally failing `workflow_examples::tests` cases:
`tutorial_source_units_generate_committed_catalog`,
`tutorial_catalog_check_passes_on_committed_tree` and
`tutorial_check_passes_on_committed_tree`. Keep canonical generator checks;
do not repeat completed gel or earlier glossary/protocol acceptance campaigns.
Report failures/skips and actual source/toolchain/platform through Steffen.

Codex used the existing `f37f7039` documentation binary for catalog regeneration
and freshness diagnostics; its generator source is unchanged at the integrated
revision. This is not new-source Rust compilation or acceptance of the revised
agent/demo code. Both existing-binary catalog and full tutorial freshness checks
pass after the repair; the two historical stale-review warnings remain. All 109
Python/source/LF-CRLF checkout checks pass without skips, alongside formatting
and documentation/maturity checks. No local Rust build/test harness, push,
dispatch or private-read run is performed. Existing scientific artifacts and
tutorial report stay intact.

## Narrowed FASTQ Repairs, 2026-10-23

The follow-up keeps the existing seed/chaining algorithm and fail-closed
budget policy. V2 records retain class plus members; the renamed counts are
partial co-nomination, never whole-read compatibility. The demo_1 regression
requires A's junction `[40]` and B's `[]` while both remain nominated. V1
screen aggregation uses `legacy_class_unresolved`; v1 summary JSON must keep
its original semantic digest. Named sample coverage distinguishes incomplete
and wholly unassessed libraries. Progress counts processed reads and compressed
read-ahead bytes, not verified read fractions. Matching pre-start cancellation
works; stale/mismatched markers fail, and polls are throttled. Explicit table
verification reads three allowlisted tables only, capped at 256 MB; ordinary
display and bounded nomination prefixes never silently rehash large sidecars.

At the integrated SHA, using the usual run-owned target, Glen/CI should run:

```sh
cargo test --locked --lib -j1 rna_panel_ -- --test-threads=1
cargo test --locked --lib -j1 rna_panel_mcp -- --test-threads=1
cargo test --locked --test glossary_cli_usage
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo check -q --locked --tests
```

Replay all of 10.04 synthetically, through CLI and MCP/shared-agent commands,
including `panel-show summary/report.json --verify-artifacts`. Native Windows
must cover quoted paths, closed handles, TSV LF/CRLF tampering and publication;
Unix adds symlink refusal. Check named missing sample B after cancellation
following sample A, and wholly unassessed samples with blank counts. Keep the
original `artifact_content_verified=false` for bounded nomination prefixes.
The live inner agent must discover the new flag and explain co-nomination,
unknown coverage and pending integrity separately; no OS helper is required.
Retain exact SHA/toolchain/platform/test counts and distinguish failures from
cascades. No local Rust or native execution is claimed here; no private reads,
performance/biological verdict, new freeze, Glen contact, push or dispatch.

## FASTQ Supporting-Read Summary Follow-Up, 2026-10-23

The new shared `ExportRnaReadPanelSummary` / `rna-reads panel-export-summary`
consumes only stored screen reports and retained panels; no original FASTQ or
large nomination TSV is reopened. It exports sample/transcript support and
exact ambiguous sets, preserves cancelled/zero-verified-file uncertainty,
refuses duplicate input hashes and incompatible panel/policy/sample bindings,
and publishes a content-bound JSON marker after closed, hashed LF TSVs.
`panel-show` reopens summaries without journaling. Tutorial 10.04 now supplies
an engine-created summary request and teaches non-additive ambiguous support.
The demo's two library headers are deliberately distinct; identical byte
inputs would otherwise correctly be refused as duplicate evidence.

Run the existing `rna_panel_` filter, glossary, capability/registry parity,
canonical catalog/manifest freshness and `cargo check -q --locked --tests`
listed below at the ordinary integrated source. Seven new `rna_panel_summary_`
tests cover actual aggregation, source deletion/no reread, LF/CRLF report JSON,
stable ordering/table hashes, chunk/sample pooling, duplicates, cancellation,
tampering/count overflow, traversal/no-overwrite and shared-shell inspection.
The existing MCP/agent-approval rehearsal now covers export and reopening.
These Rust tests are authored, not run locally; neither source checks nor
rustfmt establish their compilation or execution.

Native Windows must exercise spaced/quoted paths and handle closure before
publication. In the live inner agent, discover the exact typed request,
require approval for export, reopen it via `panel-show` and explain that the
numbers are supporting reads, not molecule copies or confirmed isoforms.
Do not silently redistribute or add ambiguous per-member support. Keep prior
native/background/real-library acceptance open; do not process private chunks
or contact Glen from this worktree without explicit authorization.

## Provider Preview Fixture Follow-Up, 2026-10-22

After integrating the owner's Windows LF commit `1cfea9e3`, three failures in
the supplied macOS log are already addressed in source: the numeric glossary
fixture (`254091c7`) and optional-source discovery in two provider tests
(`bfe685ed`). Keep their focused regressions; no new parser broadening or
missing-explicit-source tolerance is intended.

The remaining preview-reopen test mixed its synthetic preview/apply Context
with host-wide discovery in the writer. It now injects the same fixture source
list into the shared writer implementation. The public entry point still
discovers and checks every configured source; raw traversal, symlink, source-
destination and no-overwrite refusals are unchanged, not bypassed for macOS.

On the next ordinary verification, use a clean checkout and run-owned target:

```sh
cargo test --locked --lib -j1 external_service_providers::updates::tests:: -- --test-threads=1
```

Include macOS for the original reopen failure and the synthetic symlink-source
ancestor guard; retain native Windows coverage for the existing drive/verbatim
tests. Report executed test counts, failures/skips, actual source/lockfile/
toolchain/platform and logs through Steffen. The completed glossary/all-features
return below is historical evidence, not a second full campaign. No interim
freeze, bundle, push, dispatch, native capture, provider call or package build
is requested here. Local source checks do not certify the Rust paths.

## Targeted FASTQ Panel, 2026-10-22

Main integrated `9a2594fc` on `2ea67b1b`. Its recursive formatting check exposed
unresolved bare child-module paths under the path-loaded `rna_panel.rs`;
explicit `rna_panel/{input,tools,tests}.rs` paths now resolve. A new LF/CRLF
source guard fails before that repair and passes afterward. Recursive rustfmt
is not compiler/runtime acceptance; use the integrated source for the checks
below, retaining the original feature handoff's limits.

The owner's subsequent local compiler exposed E0425/E0433 for `Value` in
`tools.rs` and `tests.rs`. Repair `f37f7039` imports `serde_json::Value` in both
modules. Six source tests and recursive formatting pass, not a Rust check;
include that repair in the next ordinary RNA verification.

Changed areas: `src/engine/analysis/rna_panel.rs` and its `input.rs`, `tools.rs`,
`tests.rs`; `src/engine/protocol/rna_panel.rs`; thin enum/result/parser/MCP/
feedback/introspection wiring (including required `OpResult` literal fields),
legacy FASTA format refusal, tutorial 10.04/source/catalog/index, glossary/parity
and derived maturity source hashes. The only dependency change makes the
already-locked `md-5` version directly available. No external binary added.

Authored, not Rust-verified: `rna_panel_` covers independent source/member
identities, 20 loci / 200 ambiguous members without a 64-posting cap, both
orientations, junction flanks, strict FASTQ/gzip including truncation, quality
`>`, changed files, same-pass MD5/SHA-256, overlapping sources, cancellation,
file-boundary continuation, no overwrite, report binding and shared-shell/MCP
scientific payloads. Run normal checks at the integrated source, not a new
intermediate acceptance freeze:

```bash
cargo test --locked --lib rna_panel_
cargo test --locked --lib glossary_
cargo test --locked --lib multi_gene_sparse
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo run --locked --bin gentle_examples_docs -- tutorial-catalog-check
cargo run --locked --bin gentle_examples_docs -- tutorial-manifest-check
cargo check -q --locked --tests
```

Confirm the catalog/parity projections with the canonical Rust generators;
their Python/source checks are not generator execution. Follow [10.04](tutorial/10-04_targeted_nanopore_fastq_panel.md)
on synthetic data before any explicitly authorized real chunk. Native Windows
must exercise spaced paths, gzip failure and handle closure before publication.
For live inner-agent acceptance: each tutorial stage must use the shared route,
approve external writes/project loading, preserve compatibility sets and
support runtime inspection/cancellation without inventing completion. Native
responsiveness is unestablished; no new persistent job scheduler was added.
Any throughput/background run must identify exact source/policy/input hashes;
twenty-target capacity is a synthetic functional fixture, not a measured claim.
Do not process ongoing downloads, POD5 or private libraries without approval.

## Glossary And All-Features Follow-Up, 2026-10-22

The two repairs are independently verified at
`5ab0b02aef2c05edc15e8434be2cf1429ebf1ce5`, tree
`ac970706ac5ee4804b5022050f33479f15cb4871`; do not repeat the completed campaign
or transfer this scoped pass to later source. The previous failures were the
glossary fixture and E0521 at `70969635` / CI `39739876557`.

Retained owner-supplied archive:
`/tmp/gentle-two-repair-5ab0b02a-20261022.tar.gz` (128,815 bytes), SHA-256
`ed724dbe04adfe480872e2e072ba76a5f792f98fdc1979bdbd3ef7d7093dbed0`.
Codex validated all 26 manifested payloads in place; `SHA256SUMS` SHA-256 is
`5e50695983ef47ac6d1db6635306421fa35c7693b2f820951439d77a2be72ead`.
`REPORT.md`, `COMMANDS.md`, four command streams and CI JSON receipts preserve
the scope; binaries are identified by hash but not included.

- Linux: glossary 8/8; protocol UI 2/2, including seven widgets with one static
  step ID, distinct expected pseudonymous scopes and unchanged engine state.
- Linux: `cargo check --tests --features gui-test-support` and `--all-features`
  both returned 0. Checks compile targets; they do not execute tests.
- Toolchain/profile: Rust/Cargo 1.97.1, Debian forky x86_64, locked dev/test,
  one build job/test thread and repository-configured 16 MiB `RUST_MIN_STACK`.
  This is not default-stack evidence or native GUI acceptance.
- Same-source [CI 39774820033](https://github.com/smoe/gentle_rs/actions/runs/39774820033)
  passed both repaired gates on Windows/macOS. Windows completed all later
  gates; Linux headless and release safety passed. Full Linux was skipped.
- Overall CI remains **FAIL**: macOS later failed `Full test suite` (exit 101).
  Retained public logs returned HTTP 403 and disclose no failing test. Do not
  attribute it to either repair or claim the later provider-fixture fix
  resolves it without the named failure/results.

Next action is normal current-source CI/RNA verification plus the separate
provider-fixture follow-up above, not another glossary/protocol survey. The
RNA glossary filter covers newly added routes. Obtain the named macOS failure
from authorized logs if it remains reproducible. Disk-full warnings and a
failed final fetch are retained environment limits; confirm adequate free space
before Glen's next run, without cleaning a shared cache. A write-free remote
comparison established that the tested source was an ancestor of `1cfea9e3`,
not evidence for that descendant. No local Rust, push, dispatch or rerun was
performed when reconciling this return.

Preserve Glen's reported `70969635` passes separately: root test-target check,
protocol labware 11/11, renderer labware 2/2, root labware 16/16, rack mutation
2/2 and planning 20/20. Do not repeat those surveys solely to relabel them.
Legacy raw baselines, source-matched meshes and native carrier controls remain
open; missing OpenSCAD skips are not acceptance. The proposed commercial
product/offer/protocol slice and further bench models are separate work, not
part of these CI repairs. Historical procedures below remain source-scoped.

## Catalog-Parity Repair And Gel Return, 2026-10-20

Preserve Glen's [scoped gel pass and partial native return at bb1686bd](gel_workbench.md#integrated-linux-return-at-bb1686bd-2026-10-20).
The integration and normal [CI 39346038531](https://github.com/smoe/gentle_rs/actions/runs/39346038531)
remain failed. Both Linux/Windows primer-matrix logs show the same invalid
`catalog doctor` override assertion; do not invent a second primer defect.

After an ordinary owner-approved push, use a clean checkout and run-owned target
to verify the narrow canonical-reachability repair. Keep default GUI features
for the root test check, record the actual SHA/toolchain/platform/profile/stack
and retain complete unfiltered results:

```sh
cargo check -q --locked --tests
cargo test --locked -p gentle-protocol -j1 -- --test-threads=1
cargo test --locked --lib -j1 specificity_matrix_ -- --test-threads=1
cargo test --locked --lib -j1 introspection_distinguishes_drafts_from_validated_reports -- --test-threads=1
cargo test --locked --lib -j1 catalog_mcp_read_only_parity_and_import_confirmation -- --test-threads=1
cargo test --locked -j1 --test parity_matrix_freshness
```

Require normal Linux/Windows CI beyond the formerly failing matrix step;
later skips and a suite with the six poisoned-registry tests filtered are not
acceptance. No new interim source freeze or workflow dispatch is requested.
The repaired registry keeps all four catalog commands shell-reachable on
MCP/JS/Lua without fake engine operations or n/a exceptions; provider update
confirmation remains unchanged. Do not re-pin matrix counts merely to pass.

Review reconciliation, 2026-10-24: the repair is already in `91070ae4`, not a
new gel defect to implement. The [density rebase return](gel_workbench.md#comparison-density-verification-2026-10-20)
records all six protocol `gel_workbench` tests passing on `09639034` plus the
density slice, including the canonical parity projection. The earlier failure
above remains evidence for its older source, not today's gate result. That
filtered pass does not certify a current unfiltered protocol/CI run; keep the
commands above in normal consolidated verification, without repeating a full
campaign merely to relabel the SHA.

For gel native follow-up, actual project replacement, explicit selection
switch and separate Cancel-button behavior remain from Glen's `bb1686bd` return.
At the next consolidated native run, also check the later density UI explicitly:

- Select a saved comparison: summary and QQ must wait for **Verify and display
  saved comparison**, rather than appearing immediately. Successful replay
  reveals the shared profile panel, summary and QQ together. A deliberately
  mismatching copy must show the failure reason without a trusted curve,
  summary or QQ; raw CLI `inspect_comparison` remains diagnostic, not validation.
- Select **Signal per size width**, verify, then explicitly select **Legacy bin
  integrals** and verify again. Mode/selection/revision changes must remove the
  previous view; exercise project replacement, Cancel and close during work so
  no late result returns. Confirm the selector also drives PDF export. Use
  explicit density in the CLI for GUI/PDF comparison; omission remains legacy.
- Preserve saved integrals/statistics, zero versus gaps and the common-axis
  acquisition caveat. The separate unscaled custom-only reference is not an
  algebra coefficient or normalization. The [owner permits cross-ladder residuals](decisions.md#cross-ladder-residual-policy)
  as of 2026-10-25; do not ask again whether different ladder IDs alone forbid
  them. Retain both calibrations. Coefficients and background/marker masks
  remain undecided; no new real-data subtraction or rerun is requested here.

Retain honest platform/X11/Xvfb versus physical-host scope and original captures;
green automated guards do not close these paths. Do not repeat completed gel
arithmetic, six-command 10.03, JSON/PDF or five-page legacy comparisons just
to relabel them; recheck affected paths if their source changes. Any native
capture requires the existing owner-approved scope, not private desktop data.
Return pass/fail/not-run, command logs, source/binary/input/artifact hashes and
locations through Steffen. No live model, provider contact, package or scientific
approval is requested. The carrier checks below remain separately outstanding.

## Carrier Follow-Up, 2026-10-20

After the next ordinary owner push, verify the new carrier changes without
repeating completed sequence-design/expression/gel campaigns. Record the source
SHA, toolchain/profile, lock/binary/input hashes, commands, pass/fail/not-run
and artifact locations. No new interim freeze, publication or scientific approval.
Use a clean checkout and explicit run-owned `CARGO_TARGET_DIR`.

```sh
cargo check -q --locked --tests
cargo test --locked -p gentle-protocol labware
cargo test --locked -p gentle-render labware
cargo test --locked --lib labware_ -- --test-threads=1
cargo test --locked --lib planning_wetlab_ -- --test-threads=1
cargo test --locked --lib rack_openscad_all_templates_preserve_bottom_floor_and_open_top -- --test-threads=1
cargo check -q --locked
python3 -m unittest scripts.test_labware_models -v
```

Keep default GUI features enabled for the test-target check: library-only
`cargo check` does not compile `src/app/tests.rs`. The missing carrier-template
import there is repaired, but the new source still needs external Rust checks.
Claude's reported 20/20 ranking run used a temporary import fix; it is not a
pass for the unpatched `28d964e6` or an acceptance claim for this branch.

Check both Safe-Lock support modes, exact 4x9 dimensions/A1/D9, refusal reasons,
generic-op/shared-shell/GUI labels, unverified disclosures and unchanged project
state. Verify the canonical parity projection: only template arguments change,
no operation or command rows are added. Retain native GUI tube/carrier export
control acceptance separately; authored GUI tests are not live acceptance.

The `labware_` filter also covers these new, unrun export regressions:

Claude reports three root `labware_` failures at `09639034`: the two tests
below and `labware_carriers_preserve_legacy_json_and_export_both_support_modes`
used no placements, but carrier labels require placement scope. Each fixture
now has one unoccupied A1 guard placement. Claude's temporary macOS probe passed
12/12 after this change and was reverted; this committed repair has not been
run here. Keep the engine's empty-scope refusal intact and rerun the root
`labware_` filter at the integrated SHA; no new acceptance is claimed.

- `labware_pcr96_carrier_exports_preserve_floor_and_disclose_nonclaims`: exactly
  96 floor-preserving/open-top cuts, fixed height/opening/floor and SCAD
  assumptions/nonclaims; each of the four SVG kinds has exactly one description
  and one visible unverified disclosure. Project state remains unchanged.
- `labware_legacy_carrier_exports_omit_drawing_disclosures`: all three legacy
  templates retain the original SCAD first line, none of the new SCAD comments
  and neither disclosure in any of the four SVG kinds. This guards disclosure
  leakage, not full byte identity; the historical baselines remain required.

The earlier `rack_openscad_all_templates_preserve_bottom_floor_and_open_top`
already resizes its rack to 8x12 and counts 96 PCR96 floor-preserving cuts.
Retain that coverage; the new test adds total-cut, header and SVG-count guards.
All new Rust checks remain Glen/CI work, not a locally executed pass.

Capture and commit the nine original carrier baselines from `1c5be317` source-
matched binaries: SCAD, simulation JSON and the embedded physical spec for each
legacy template on one fixed synthetic rack. Document exact recreation and raw
hashes; add scoped LF rules, LF/CRLF coverage in `scripts.test_tutorial_checkouts`
and Rust byte-comparison tests. Compare against new-source raw exports, without
normalizing or regenerating away differences. Glen's earlier 9/9 comparison is
historical evidence, not a committed fixture or this source's pass.

Render both new Safe-Lock carriers in OpenSCAD from real source-matched GENtle
exports. Require finite nonempty meshes, expected bounds, 2-mm floors and open
bores. Retain original SCAD/spec/mesh hashes and logs; optional tool skips are
not acceptance. Physical fit, cap clearance and thermal/material suitability
remain unverified even with successful mesh renders.

Include the conditional PCR slice, whose source gate is now met: Eppendorf
`0012 595 670-00/2026-05`, page 4, raw SHA-256
`67b84dea22bcb2cbea26907b8f1116c506210f1917744feb64b2fafb6745b40d`.
The [source record](../data/labware/README.md#pcr-02-ml-drawing-and-96-position-carrier)
explains maximum OD 6.0 versus main-cylinder 5.8 mm, unknown seating collar,
6.6-mm opening, 13-mm floor bore and checked nominal standard coordinates.
The same filters cover new footprint/A1/H12, 8x11/12x8/16x24 refusals,
no-overwrite, shell/generic-op and GUI labels. Run the renderer's complete
tests too (`cargo test --locked -p gentle-render`): its older source-binding
test name does not include `labware`. Verify clean `InvalidInput` when exporting
the PCR tube visualization (unknown collar), not a crash or invented dimension.
Render a third, PCR96 carrier mesh using an exact 8x12 grid. Keep all legacy
and Safe-Lock checks; inspect new-template SVG disclosures without inferring
full ANSI/SLAS compliance. No dense spacing, boxes, 0.1 ml profile or thermal
contact was implemented in the original carrier tranche. The offline mesh
harness uses literal parameters and
actual SCAD format strings, not a freshly built Rust exporter.

Owner requirements clarified 2026-10-20: offer **adjustable** compact spacing
below 20 mm in a separate extension without changing existing templates.
Storage prioritizes dense transport/freezer organization rather than treatment
access. Commercial paper/cardboard boxes should be modeled as possible storage
destinations, not presumed fabrication targets. Box dimensions/dividers/lids,
closed-cap clearance, fit and material/temperature qualification still need
their own evidence/contract. Box formats must be variable, not fixed to 9x9.
Compact spacing is now authored separately; no physical box model or new release
gate follows.

### Compact Carrier Follow-Up

Run the existing protocol/root `labware` filters and `cargo check --tests` on the
integrated commit, then canonical parity freshness and glossary/interface checks.
New authored tests (not run here):

- `labware_compact_pitch_bounds_preserve_conservative_edges`
- `labware_compact_pitch_decimal_and_legacy_rack_bytes_are_exact`
- `labware_compact_pitch_exports_geometry_persists_and_clears_without_drift`
- `labware_compact_pitch_refuses_invalid_modes_and_other_templates_without_writes`
- `labware_compact_pitch_shell_and_generic_operation_share_mutation_contract`
- `labware_compact_pitch_keeps_variable_grid_formats_independent_of_spacing`

Existing GUI-label coverage also checks pitch-control visibility. Review the
numeric field, exact mm/micrometre conversion, Clear, seating-mode refusal and
unsupported-template message natively. Source-matched OpenSCAD renders should
include one compact layout. Retain conservative/legacy byte baselines; none was
regenerated here. The parity Markdown projection was updated without running
its Rust generator: Glen/CI must confirm byte-exact canonical freshness.
The new operation is shared by GUI, CLI, generic MCP/scripting/workflow routes;
no agent auto-execution/approval exemption is added.

### Rack-Mutation Confirmation Follow-Up

Claude reports protocol 11/11 and render 2/2, but root `labware_` 15/16 at
`17f846eb`: `labware_compact_pitch_shell_and_generic_operation_share_mutation_contract`
correctly rejects the missing mutation classification. All eleven rack mutators,
not just compact pitch, were absent from `ShellCommand::is_state_mutating`, so
`--allow-auto-exec` could bypass the assistant's mutation-confirmation gate.
The repair adds exactly those eleven variants; no geometry, export strings,
non-rack classification or MCP tool permission changes. The shared assistant,
GUI and planner admission now use the corrected flag. Claude's reverted
267-test probe is not execution of this committed source or current acceptance.

New Rust tests, authored but not run here:

- `rack_mutators_are_state_mutating_and_block_agent_auto_exec`: parse all eleven
  commands; retain read-only `racks show`; block an actual auto-intent pitch
  suggestion with `AGENT_MUTATION_CONFIRMATION_REQUIRED` and unchanged rack,
  state, journal and revision; explicit index selection can still execute it.
- `rack_mutation_admission_cases_stay_in_sync_with_executor_dispatch`: bind the
  per-variant cases and classifier to the actual rack mutation handler. This
  source-set guard plus adjacent comments catches a future dispatch-list drift;
  it does not substitute for the execution regression above.

After ordinary integration, use the same clean checkout/run-owned target rules
above, default GUI features and the exact tested SHA. In addition to the
protocol/render and carrier checks, run the **full root library suite**, not
only the `labware_` filter:

```sh
cargo check -q --locked --tests
cargo test --locked --lib rack_mutat -- --test-threads=1
cargo test --locked --lib labware_ -- --test-threads=1
cargo test --locked --lib -- --test-threads=1
cargo check -q --locked
```

Keep agent-parity/mutation/GUI results and primary failures separate from any
mutex-poisoning cascade. The source-only Python LF/CRLF membership/negative
control checks are not Rust, live-agent or native acceptance. OpenSCAD renders,
nine legacy baseline files and native GUI review remain outstanding. Include
11.75..12.049 mm: setter admission uses the collar limit; floor-seated export
must clearly refuse it. Optional positional disclosure cleanup was deferred.

Export classification inventory (unchanged): all seven named rack SVG/SCAD/JSON
exports return `is_state_mutating=false`, as do named `SaveProject`, `ExportPool`
and `LaddersExport`. Named guide oligo/protocol exports return `true`, as does
the conservative generic `op`/`workflow` route. Thus rack exports match some
file-writing routes, not a uniform export policy. `false` means no declared
project mutation here, not absence of filesystem side effects; no broader
export-policy repair is part of this commit.

The clarified future box requirement is floor-supported 1.5/2.0 ml tubes with
snug divider contact in openings that are square in top view. Cells are
separated by wall thickness; literal tube-to-tube contact is not required.
The [Semadeni example](../data/labware/README.md#different-storage-box-formats)
supplies a source lead, not interior/divider/lid dimensions or Safe-Lock fit.
Variable formats and all physical qualification remain separate future work.

## Historical Covering Request, 2026-10-19

Updated 2026-10-19: after his push, Steffen should forward the
**[concise parity/expression follow-up](glen_parity_expression_followup_20261019.md)**.
Only that message is active for this request. The execution checklists and
earlier capture requests below are retained background, not renewed authority
or a request to repeat completed checks. No release freeze is selected.

The original covering request was prepared 2026-10-18. Its source-bound return
and historical procedures follow; continue development independently.

## Returned Results And Remaining Request

Glen executed the pre-inspector checks below at
`394dd5d0523ea91c19a8eb94d5ed8ad7bd5299b2`; see the
[verified evidence and six-item status](glen_audit_followups_20261008.md#focused-verification-at-394dd5d-2026-10-18).
README 4/4, sequence-design Rust 67/67 plus actual CLI suite 22/22, starter
suite 6/6, workflows 95/95 and rack/labware 127/127 pass. The nine historical
carrier outputs are byteidentical. Keep those checks closed at that source;
do not reissue this entire original checklist as unfinished work.

One reported defect is repaired locally: the parity matrix's two missing
tube-operation rows/counts. Codex applies Glen's canonical generated delta and
preserves newer inspector surfacing without running local Rust. On the next
ordinary owner-pushed revision, run `parity-matrix-check` with its source-matched
docs binary, `cargo test --locked -j1 --test parity_matrix_freshness`, and inspect
normal CI. Run the seven inspector regressions and `gentle-shell` UI-intent
checks in original item 3: `5893ebb0` postdates the tested source. Its updated
09.03 tutorial/catalog projection also needs source-matched freshness checking.
No new bundle or intermediate freeze is requested; final-candidate affected
checks remain on the normal release path. OpenSCAD mesh execution remains
`not_run` on Glen's host, not a failed Rust exporter. Native, agent, scroll and
package boundaries below receive no new acceptance from this return.

## Source And Isolation

Fetch published `smoe/gentle_rs main` once and record the resolved full SHA and
tree. Test that source in a clean detached checkout with an explicitly set,
run-owned absolute `CARGO_TARGET_DIR` outside the checkout. Record toolchain,
platform, features/profile and raw lockfile SHA-256. Use only binaries built
from this checkout; preserve the primary checkout, shared caches and historical
evidence. A later main commit does not invalidate this scoped report or require
another synchronization request, but it is not automatically accepted by it.

No push, workflow dispatch, tag, release publication, provider contact, live model
or package/container build is requested. The owner's follow-up additionally
requests only the four synthetic inspector tutorial screenshots described below;
no unrelated GUI/desktop/private-data capture is authorized. Existing normal CI
can be inspected without restarting it. Keep `package-opt1-fallback` separate;
the selected .13 packaging recipe remains opt2. Codex runs no local Rust.

## Original Focused Checks And Inspector Follow-Up

1. **README metadata repair:** verify the `57a2403a` repair is present. Run the
   release-version regression below and inspect the new normal push CI when
   available. Historical run `39066178748` at `1c5be317` failed the README
   label guard on Linux and Windows; it did not test the repair. Both root
   libraries and all five stripped-child families passed at that older SHA.
   Preserve those bounded results rather than calling the failed run green.

   ```sh
   cargo test --locked -j1 --test release_version_consistency
   ```

2. **Sequence-design follow-up:** verify `8a0e2670` and `958b2238` are included.
   The first retains a freshly validated best after edited-candidate match-cap
   overflow, marking search incomplete without optimality/infeasibility proofs
   instead of returning Invalid and discarding it. The second adds independent
   180 single-symbol and 10,800 paired-symbol IUPAC cases without changing the
   production matcher. Eighteen offline oracle/harness passes are not Rust
   acceptance; execution subsequently passed at `394dd5d` as recorded above.
   Execute the core/protocol tests that the failed older workspace run did not
   reach, plus the shared engine and CLI boundaries:

   ```sh
   cargo test --locked -p gentle-sequence-design -j1
   cargo test --locked -p gentle-protocol -j1 sequence_design
   cargo test --locked --lib --no-default-features -j1 sequence_design -- --test-threads=1
   cargo test --locked --bin gentle_cli --no-default-features -j1 sequence_design -- --test-threads=1
   cargo check -q --locked -j1
   cargo build --locked --no-default-features -j1 --bin gentle_cli --bin gentle_examples_docs
   ```

   Replay tutorial 06.07 and the four actual CLI acceptance tests using the
   [functional procedure](sequence_design_acceptance.md#exact-candidate-functional-checks).
   For these functional checks the dev binaries just built are sufficient;
   record that profile, without claiming package or performance acceptance.
   Set `GENTLE_TUTORIAL_BIN_DIR` to the absolute run-owned target's `debug`
   directory, `GENTLE_SEQUENCE_DESIGN_CANDIDATE` to the recorded full HEAD, and
   `GENTLE_SEQUENCE_DESIGN_EVIDENCE_DIR` to a fresh absolute evidence directory:

   ```sh
   python3 -m unittest scripts.test_sequence_design_audit_oracle scripts.test_sequence_design_tutorial -v
   ```

   All four actual CLI tests must execute; skips and the four Rust CLI unit
   tests are not substitutes. Check retained-best behavior after match-cap
   overflow, incomplete status without false proof flags, cancellation clearing
   output, lexical ties and independent IUPAC coordinates/strands/palindromes.
   Preserve both strategies' joint motif/global/window-GC, tampering,
   materialization, save/reopen and same-session undo/redo checks. Historical
   `de9df8af` prototype acceptance stays separately attributed.

3. **Integrated .13 expression/storage/tutorial paths:** follow the focused
   commands, generator freshness/parity checks and two real-CLI starter replays
   in [the preparation handoff](internal_13_preparation_handoff.md#external-completion-before-candidate-selection).
   Its owner dispatch/publication section is not part of this request. Use the
   same freshly built binaries with `scripts.test_release_workflow_tutorials`;
   no optional-binary skip is acceptance. Retain LF/CRLF and quoted-path results,
   deterministic container order, conservative expression refusals and unchanged
   review state/files. Report drift; do not regenerate historical evidence.

   The merged `gentle_rs_3` inspector at `5893ebb03b4915ba1af2f77f1d7681904266b6f9`
   supersedes the older preparation note deferring all inspector work. Include:

   ```sh
   cargo test --locked --lib -j1 protein_expression_inspector_ -- --test-threads=1
   cargo test --locked -p gentle-shell -j1 ui_intent
   cargo test --locked --lib -j1 mcp_server::protein_expression::tests:: -- --test-threads=1
   cargo test --locked --lib -j1 protein_expression_inner_agent_suggestion_ -- --test-threads=1
   cargo test --locked --lib -j1 protein_expression_agent_card_ -- --test-threads=1
   cargo test --locked --lib --features gui-test-support -j1 protein_expression_capture_targets_ -- --test-threads=1
   ```

   Require all seven inspector regressions to execute. Verify shared-report
   parity, explicit/absent inputs, typed withholding, copy-only suggestions,
   stale/closed-worker rejection and unchanged project/history/files. Check
   canonical tutorial/catalog/parity freshness at this source; static metadata
   edits are not generator acceptance. Native window review stays at release
   preparation. The later owner request additionally asks for the four
   [09.03 screenshot checkpoints](tutorial/09-03_protein_expression_handoff_cli.md#required-tutorial-screenshots)
   from a verified binary and clean synthetic profile. Four manual macOS PNGs
   at `9f01098b` are now retained with binary/input/image hashes and limits in
   the linked archive. They do not contain the uncompiled parity follow-up.
   Review them; retain typed GUI reports and semantic snapshots at the integrated
   candidate separately, without treating these images as automated acceptance.
   No provider/model contact, unrelated capture or automatic capture is allowed.
   Also execute the named MCP `protein_expression_handoff` discovery/report/
   refusal regressions, the deterministic inner-agent command comparison and
   the `gui-test-support` status/locator test. The new control card is guidance,
   not evidence of a model's behavior; compare actual reports on saved state.

4. **Merged rack safety and tube reference models:** verify the integrated
   `gentle_rs_2_main` commit `5b1da7d3981e7a2a4280cadc4b86b9345d82896c`.
   Run its protocol/renderer, headless engine/shell and default GUI regressions,
   then existing rack export and capability/parity checks:

   ```sh
   cargo test --locked -p gentle-protocol -j1 labware::tests::
   cargo test --locked -p gentle-render -j1 labware::tests::
   cargo test --locked --lib --no-default-features -j1 labware_ -- --test-threads=1
   cargo test --locked --lib -j1 labware_ -- --test-threads=1
   cargo test --locked --lib -j1 rack_ -- --test-threads=1
   python3 -m unittest scripts.test_labware_models -v
   ```

   Require positive, overflow-checked admission through 1,536 slots, refusal
   above the limit, unchanged unsupported saved racks on reopen, and refusal
   before overwriting any existing export. Exercise generic `op` inspection
   and all three tube exports with quoted paths: deterministic SCAD, legacy
   result decoding, explicit artifact-write confirmation and unchanged project
   state/history. Generator/parity checks must include the new operation IDs.
   OpenSCAD mesh checks alone do not execute the Rust exporter; missing
   OpenSCAD is a recorded skip, not mesh acceptance. Keep source URL/digest/
   revision/page and exact/approximate/unknown dimensions intact.

   Follow [the provenance and legacy-baseline boundary](../data/labware/README.md#verification-boundary).
   In separate isolated checkouts, compare raw SCAD, simulation JSON and spec
   JSON for the same fixed rack and all three existing carrier templates,
   before this feature at `1c5be3175323eec41359b980b5a41e22c7edaec1` and at
   the tested integrated source. Retain input and before/after hashes outside
   Git; no source-only or template-only result proves byte identity. Native
   Rack controls remain a release-stage check, not new capture authority here.
   Tube models are body/collar visualizations with cap/hinge/cavity omitted,
   not printable vessels or fit evidence. New carriers, support/pitch,
   96-position interpretation and storage boxes remain deferred owner choices.

Use `--locked` throughout. If dependencies are already cached, `--offline` may
also be used; otherwise record retrieval failures as environmental. Do not
change the lockfile or widen resource limits silently to obtain a pass. Missing
tools, zero selected tests, timeouts and unrun platforms remain explicit.

## Remaining Six-Gate Boundaries

Use the [six-item ledger](glen_audit_followups_20261008.md) as the status owner,
and the [release-stage request](glen_final_candidate_execution_handoff.md) only
when release preparation is authorized. Audit preservation, 08.04 contracts,
live-agent acceptance, honest visual evidence, scroll benchmarks and package
acceptance retain their individual evidence boundaries. These focused checks
do not close all six. Existing tutorial/link repairs and historical audits
must not be reissued as unfinished work. Native/physical-host, model, scroll,
package and scientific verdicts stay separate; this request authorizes none
of their additional execution.

## How To Return Results

Send **one concise response to Steffen for forwarding into the main Codex
chat**. Do not rely on direct Glen-to-Codex messaging or assume a path on one
host exists on the other. Include:

- Tested full SHA/tree, platform, Rust/Cargo versions, features/profile,
  lockfile and binary hashes, and any same-source CI URL/conclusion.
- Per-check `pass`, `fail`, `not_run` or `blocked`, actual selected test counts,
  skips, exact commands/exit codes and primary failure versus cascades.
- Remaining six-gate boundaries, distinguishing new results from historical
  passes and software execution from native/scientific/release acceptance.
- Evidence location and SHA-256 of its manifest, plus any concrete defect and
  minimal proposed repair. Do not modify or push source as part of this check.

Supply a small folder/archive containing `REPORT.md`, command/exit-code logs,
raw stdout/stderr, CI/CLI receipts and a per-file SHA-256 manifest. Steffen can
copy it to `/tmp/glen-verification-<short-sha>-<date>/` on the Codex host and
relay the path and manifest hash here. Keep bulky binaries, raw captures and
reference inputs on Glen's host with hashes/locations; transfer them only on
request. Never include credentials or private projects. Record content hashes
as identity checks, not authenticated provenance or human scientific approval.

## Preserved Prompt Records: Not Active Instructions

The following received files are committed verbatim for history:

- [Release-link repair prompt](release_handoff_link_repair_prompt.md): repair
  already closed; its preparation-time failure and old runner/version status
  are historical, not a request to restore the deleted .12 planning ledger.
- [Sequence-design closure prompt](sequence_design_closure_codex_prompt_20261014.md):
  closure was reconciled in the plan/acceptance records. The current
  [compiler decision](decisions.md#workspace-compiler-compatibility) supersedes
  its broader proposed workspace minimum/publication gate; select no MSRV.
  The [issued handoff's SHA correction](glen_sequence_design_acceptance_handoff_20261014.md)
  remains intact.
- [Original Claude VKORC1 prompt](vkorc1_acceptance_codex_prompt_20261014.original-claude.md):
  retain the original separately from the [reconciled task](vkorc1_acceptance_codex_prompt_20261014.md).
  Its in-progress CI, undecided .13 aim and other snapshot claims are not
  current status. `--offline` is not a substitute for `--locked`, nor proof
  that a lockfile was rewritten. No old embedded dispatch is authorized here.
