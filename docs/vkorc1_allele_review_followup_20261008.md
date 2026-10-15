# VKORC1 Allele Review Follow-Up

## Scope And Baseline

Reconcile Claude's read-only review of `91038fe0` against this checkout at
`d8a2a1c1` on `gentle_rs_2_main`. At the initial source reconciliation,
local `main` was `e04f6044` and already contained `f14255cb`.
No remote fetch, private-data run or release acceptance
is implied. The review's M1/M2/L1 baseline omissions have landed in `f14255cb`;
this follow-up closes a remaining M2 mixed-allele case and strengthens GUI M3
execution atomicity.

## Findings And Implementation

| Item | Current implementation | Regression coverage |
| --- | --- | --- |
| M1: reference and assembly evidence | `DbsnpResolvedPlacement::reference_check` compares the complete resolvable SPDI reference span with the extracted forward DNA. Missing, cropped or ambiguous reference bases are `unavailable`. `FetchDbSnpRegion` records `dbsnp_reference_check`, a separate `dbsnp_assembly_check` and fallback provenance; mismatch/unavailable and every fallback produce warnings. A matching base never authorizes an unverified assembly. | `dbsnp_fetch_checks_reference_and_assembly_and_retains_non_snv_spdi`; `dbsnp_reference_check_validates_the_full_reference_and_keeps_unknown_explicit` |
| M2: lossless non-SNV evidence | `build_dbsnp_variant_marker_feature` emits `vcf_ref`/`vcf_alt` only when every retained SPDI allele has canonical single-base deleted/inserted strings. Mixed SNP/indel placements must not become SNPs merely because an indel has a different reference span and is excluded from normalized alternates. All selected-placement raw SPDI entries remain under `dbsnp_spdi`, including empty insertion/deletion strings. Non-SNV focal markers are not full variant intervals or materializable SNVs. | Fetch regression covers insertion, deletion, multi-base replacement, ambiguous bases and mixed SNP/deletion or SNP/insertion cases, checking raw SPDI equality, absence of VCF allele qualifiers and refusal to materialize. |
| L1: explicit handoff alternate | `reporter_construct_handoff_commands` uses the shared materialization validator to supply one reported, reference-distinct A/C/G/T alternate from the loaded source. Missing, multiallelic or invalid evidence leaves the field unset and explicitly requires review; execution still revalidates. | `reporter_construct_handoff_binds_unique_loaded_alternates_or_requires_review`, including parser/typed-operation round-trip and unchanged state. |
| M3: all-or-nothing GUI pair | Preserve shared preflight for both choices. Apply the two canonical operations under one write lock using `with_rollback_on_error`; do not publish success feedback or new output names until both succeed. The rollback restores project, journal and undo/redo history, while revisions and consumed operation identities remain invalidated/non-reused. Successful pairs still have two undo entries. | `variant_followup_allele_pair_preflights_before_creation_and_retries_without_orphans` additionally checks the active view; `variant_followup_allele_pair_rolls_back_late_refusal_and_preserves_history` exercises a successful reference followed by a refused alternate, prior undo/redo preservation and an unsuffixed corrected retry. |

The GUI transaction is thin composition, not another allele-selection algorithm
or a new operation/protocol. The late-refusal test calls the same private
execution seam as the GUI directly to reach the rollback boundary, rather
than treating a preflight-only refusal as transaction coverage.

Local macOS checks: formatting and whitespace pass; the Python tutorial
checkout/walkthrough suite runs 29 cases, with 26 passes and three explicit
real-binary replay skips. These checks do not execute the new Rust regressions.
Session-close has no failures; its warnings concern the intentional edits and
manual scope confirmation. No Rust build/test or remote operation was run.

## Request To Glen / CI

Rust execution remains deferred under the owner's instruction. Use the final
committed candidate, record its full SHA, dirty status, lockfile/toolchain,
platform and commands, and run:

```sh
cargo test --locked --lib dbsnp_ -- --test-threads=1
cargo test --locked --lib materialize_variant_allele -- --test-threads=1
cargo test --locked --lib reporter_construct_handoff_ -- --test-threads=1
cargo test --locked --lib variant_followup_allele_pair_ -- --test-threads=1
cargo check -q --locked
```

Report individual pass/fail/skip counts and retain failure output. At exact
`df5c8244787d8c920ee171ee55167de55e8a7294`, Windows job `115778176141`
executes the file-backed reference/assembly fetch, annotated extraction, SPDI
GenBank round-trip and late pair-rollback regressions; their named `ok` lines
are in the [retained CI scope](audits/glen_followups_20261008/candidate_df5c8244/README.md).
Those dbSNP tests are not Unix-only. Keep this historical unit-test evidence
separate from newer-source/default-stack, real-data and native GUI acceptance;
Unix-only mock-BLAST tests elsewhere retain their own platform limits.

The named Windows passes are
`engine::tests::dbsnp_fetch_checks_reference_and_assembly_and_retains_non_snv_spdi`,
`engine::tests::test_fetch_dbsnp_region_operation_extracts_annotated_slice_and_provenance`,
`dna_sequence::tests::dbsnp_spdi_genbank_round_trip_preserves_wrapped_json` and
`main_area_dna::tests::variant_followup_allele_pair_rolls_back_late_refusal_and_preserves_history`.
The mock payload uses `file://`, not Unix-only mock-resource machinery.

## Still Open

### 2026-10-14 Default-Stack Follow-Up

Glen's retained `024ed6c0520d38c4221517219d44c792084de4e5` Linux/Rust 1.96.1
logs show the four original filters selecting 11/2/7/2 tests, then aborting with
test-thread stack overflow (exit 101). Explicit `RUST_MIN_STACK=16777216`
reruns pass all four families; these are labelled diagnostic, not replacement
passes. The `dbsnp_spdi_` filter passes 2/2 without that override.

The local repair moves only these heavyweight fixture bodies onto an explicit
16-MiB test worker, preserves their names/assertions and propagates panics.
Helper regressions exercise a small caller stack and a worker failure. Native
CI initially used a shell `unset`, which does not prevent Cargo from restoring
`RUST_MIN_STACK` through `.cargo/config.toml`. That step did not establish
default-stack acceptance. No parser, allele-selection or production stack policy
changes are implied by the test-worker repair.

The corrected CI uses a Cargo runner that removes the override at the final
executable-spawn boundary, preserving Cargo's library paths and normal full-suite
settings. In a clean checkout of the final committed candidate, run:

```sh
runner_config="$(python3 scripts/default_stack_test_runner.py --cargo-config)"
for filter in dbsnp_ materialize_variant_allele reporter_construct_handoff_ \
  variant_followup_allele_pair_ test_support::heavy_stack_tests; do
  cargo test --config "$runner_config" --locked --lib "$filter" -- --test-threads=1 || exit
done
cargo test --locked --lib tss_profile_export::tests::svg_dimensions_ -- --test-threads=1
cargo check -q --locked
```

Use `python` on native Windows. Retain every runner JSON receipt and confirm
`child_rust_min_stack:null`; an absent receipt is not a default-stack pass
(a target-specific Cargo runner may take precedence). Keep the original Linux
11/2/7/2 selections and two helper tests visible, report platform-specific skips,
and retain exact SHA, lockfile/binary hashes, toolchain, commands and exit codes.
The two TSS parser tests cover the owner-observed quick-xml 0.41-to-0.42 boundary.
Rust verification remains external and pending, not inherited from `024ed6c0`
or from the offline Python runner tests.

### Nested dbSNP Worker Repair (2026-10-14)

At exact `de9df8af0a56923df01a7492f61005b1fe079d8c`, Linux job
`115975878354` in run `38599805415` confirms `child_rust_min_stack:null`,
selects 11 `dbsnp_` tests and aborts with stack overflow/exit 134 in
`dbsnp_dialog_fetch_extracts_region_and_opens_window`. The remaining four
filters do not execute after that first failure. This is a failed acceptance
run, not a missing runner or a passed diagnostic replacement.

The same [terminal CI record](audits/glen_followups_20261008/candidate_de9df8af/README.md)
shows Windows failing at that same first fixture with the stripped-child receipt,
child exit `0xffffffff` and step exit 127. Both later four-filter loops and
full suites remain unexecuted; selected-platform failure is a cascade. macOS
is skipped, not a default-stack pass. Neither platform accepts the repair.

The test's fixture thread is explicitly named and sized at 16 MiB; the dbSNP
dialog still spawned an unnamed default-sized production worker around detached
engine execution. The local repair gives that worker the adjacent GenBank
ingress policy: explicit 16 MiB, named `gentle-dbsnp-fetch`, with visible creation
failure and no pending task left behind. Engine, allele/provenance rules,
fixture assertions, runner and global stack settings are unchanged.

The new `dbsnp_worker_stack_is_explicit_and_named` uses a hand-crafted 4-MiB
frame through the production spawning seam. Expect 12 Linux `dbsnp_` tests at
the repaired source, including the original dialog boundary; keep per-platform
counts and all five stripped-child receipts. Rust execution of the repair is
still pending externally. The failed SHA/logs remain immutable, and no skip,
larger ambient stack or ordinary CI summary substitutes for the rerun.

### Focused Linux Repair Verification (2026-10-15)

Glen verifies exact `6b7109dbb7fcc419866e1537dc9aa4cfa1ab70be`, tree
`8d6e0ee93579c55873f3bb5cf53a1b6e6673018c`, in a clean detached checkout.
On Debian x86_64/Rust-Cargo 1.96.1, the default-feature test profile passes all
12 `dbsnp_` tests with an actual runner receipt: inherited stack `16777216`,
child stack `null`. Named passes include the formerly overflowing GUI fetch
and `dbsnp_worker_stack_is_explicit_and_named`. Exit is zero; tests take 2.47 s
after a 6:06 cold build. Runner Python self-tests pass 6/6; formatting and
whitespace checks pass. This closes the reproduced Linux failure at that source.

Producer evidence is retained at
`/home/clawbio/.openclaw/workspace/artifacts/gentle-6b7109db-verification-20261015/`;
the owner-supplied copy is `/private/tmp/gentle-6b7109db-verification-20261015/`.
Codex verifies all 34 inventory entries, safe exact membership (36 files including
the inventory/companion), candidate tree and runner/test/exit records. Inventory
SHA-256: `93ba9a64fa20ea318aeef3dd98a744eb8c86c92f529a93128f75980e0ba58c78`.
The included bundle matches the original transfer. Test-binary SHA-256
`1716c9c53cfe0f6591ccd121cff75eae61c0888d16b89441ad75d986ac76ca58`
is a producer record; the binary was removed with run-owned scratch, not supplied
for independent rehashing. The runner's executable path confirms the isolated
target, despite the earlier environment snapshot showing an ambient shared target.

No Windows, other four Rust filters, broad suite or GUI replay ran in this
focused verification. The six Python runner checks are not those Rust filters.
Keep the failed `de9df8af` CI unchanged; no newer-source, cross-platform or
release acceptance is inferred. Remaining consolidated gates follow the
[release-preparation cadence](decisions.md#development-and-release-verification-cadence),
without another interim freeze, bundle or dispatch. Codex ran no local Rust.

- Tutorial 08.04 live native GUI: confirm reference/alternate refusal, retry,
  unchanged source and honest genomic-forward C/T evidence. Reverse-strand
  gene nomenclature must not silently complement the selected allele.
- Raw refSNP archive/checksum and selected assembly/placement verification
  against the actual prepared reference. Synthetic mocks and a retained raw
  response alone are not this scientific acceptance.
- Package/release acceptance and real reporter-study interpretation remain
  independent owner/Glen decisions; no network retrieval or WIP artifact
  regeneration was performed for this follow-up.
- L2/L3 script-template, field-label and path ergonomics, and the successful
  pair's two-Undo behavior remain separate cosmetic follow-ups, not hidden
  additions to this minimal robustness change.

Return the exact tested candidate and retained evidence paths/hashes with a
typed verdict distinguishing product defects, missing/incompatible data,
environmental failures and pending human scientific review. Do not combine
passes from the old `91038fe0` review or nearby revisions into a new acceptance.

## Public Audit And Remaining Acceptance Preparation (2026-10-15)

The Python public audit now requires the saved rs9923231 marker's reference
and assembly checks to match, assembly fallback to be false, SNV class,
genomic-forward C with a reported T, GRCh38 and NC_000016.* at one-based
31096368. It retains the actual complete alternate list and raw SPDI entries
before validating them; an invalid marker stops the run with unverified evidence
retained. No A,G,T set is assumed. A multiallelic bare choice must refuse with
unchanged persisted state; a single alternate is explicitly not applicable,
not a fabricated refusal. Thirty helper tests pass using hand-crafted data
and mocked ingress, without executing GENtle or contacting a provider.

Tasks 1/2 of the [bounded acceptance request](vkorc1_acceptance_codex_prompt_20261014.md)
were already represented in the retained Windows correction and terminal
`de9df8af` records. Rehashing the original Linux/Windows logs locally confirms
their existing digests; no historical audit file is changed or fetched again.
Do not replace an aborted filter's absent result summary with invented
pass/fail/ignored totals:

| Stripped-child Filter At `de9df8af` | Linux | Windows |
| --- | --- | --- |
| `dbsnp_` | 11 selected; stack abort/134; no completed P/F/I summary; null child-stack receipt present | 11 selected; child 0xffffffff/step 127; no completed P/F/I summary; null child-stack receipt present |
| `materialize_variant_allele` | Not run; no receipt | Not run; no receipt |
| `reporter_construct_handoff_` | Not run; no receipt | Not run; no receipt |
| `variant_followup_allele_pair_` | Not run; no receipt | Not run; no receipt |
| `test_support::heavy_stack_tests` | Not run; no receipt | Not run; no receipt |

macOS was skipped. The primary failure is a product stack overflow, not an
environmental isolation failure; selected-platform failure is a cascade.
Glen's older Linux selections were 11/2/7/2 plus two helper tests. The repair
adds one dbSNP probe, explaining the [later 12/12 focused Linux pass](#focused-linux-repair-verification-2026-10-15)
at `6b7109db`. That closes only the reproduced Linux defect there, not Windows,
the four unrun filters, newer source, or the quick-xml `svg_dimensions_` gate.

The [ledger's prepared owner request](glen_audit_followups_20261008.md#prepared-vkorc1-next-actions-2026-10-15)
and [physical-host checklist](glen_final_candidate_execution_handoff.md#vkorc1-0804-physical-host-checklist)
carry the residual work. No dispatch or physical run has occurred here.

| Item | Typed Status / Scope |
| --- | --- |
| Windows synthetic-coverage claim | **closed** at `df5c8244`; not physical Windows acceptance |
| Public audit assertions | **closed** as a Python implementation/offline regression repair; not live-data acceptance |
| Remaining Windows/default-stack and public replay | **pending owner dispatch**, with successful same-source normal CI required first |
| Physical 08.04 GUI and native Windows | **pending Glen**, with approved isolation and capture consent |
| .13 aim and release acceptance | **owner decision**; package/extracted-installer acceptance stays separate |
