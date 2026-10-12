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

Report individual pass/fail/skip counts and retain failure output. The synthetic
fetch tests use Unix-only mock-resource machinery; a green native-Windows job
does not prove those cases. Keep source inspection, focused tests, full-suite
and native GUI verdicts separate. No focused Rust results for this diff are
claimed here.

## Still Open

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
