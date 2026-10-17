# Codex Prompt: VKORC1 / Tutorial 08.04 Remaining Acceptance (2026-10-14)

Prepared by Claude from a read-only reconciliation of local `main` at
`de9df8af` (merge of `gentle_rs_3`). Re-verify HEAD before you start; Codex
worktrees merge into `main` mid-session.

## Goal

Close or precisely re-scope the four VKORC1 items still open after the allele
review: live dbSNP/GRCh38 replay, native 08.04 GUI, Windows, and release
acceptance. Most of the remaining work needs an owner dispatch or a human on a
physical host. Your job is to:

- record what CI already proves;
- fix one stale claim;
- close one real gap in the public-audit script;
- prepare, but **not execute**, the single dispatch and the Glen checklist that
  remain.

## Hard Boundaries

- No push, tag, publish, release, model call or network retrieval.
- No `workflow_dispatch`. Prepare the exact inputs for the owner instead.
- No local Rust build or test unless the owner explicitly allows it. If allowed,
  always pass `--locked`. A bare `--offline` rewrites `Cargo.lock` and produces
  fake compile errors.
- Python unit tests and the maintenance scripts are allowed offline.
- Do not change engine allele-selection, dbSNP resolution, GUI transaction or
  `reporter_ops` code. Their review findings (M1-M3, L1, SPDI GenBank
  round-trip) are closed at `60fad40d` and `dfadf2a5`. Do not re-plan them.
- Preserve every historical receipt and audit byte. Add new evidence subsets;
  never rewrite old ones.
- Do not choose the `.13` release aim (see Task 6).

## Verified Facts (do not re-derive; re-check only if HEAD moved)

| Fact | Evidence |
| --- | --- |
| Normal CI at `df5c8244` passed | Run `38541824374`: macOS and Windows success; full Linux skipped. |
| Windows already executes the synthetic dbSNP fetch tests | Windows job `115778176141` at `df5c8244` logs `ok` for `engine::tests::dbsnp_fetch_checks_reference_and_assembly_and_retains_non_snv_spdi`, `engine::tests::test_fetch_dbsnp_region_operation_extracts_annotated_slice_and_provenance`, `dna_sequence::tests::dbsnp_spdi_genbank_round_trip_preserves_wrapped_json` and `main_area_dna::tests::variant_followup_allele_pair_rolls_back_late_refusal_and_preserves_history`. These tests have no `cfg(unix)`; the mock is a `file://` URL. |
| The corrected default-stack runner first reaches CI at `de9df8af` | `fcb3a867` adds `scripts/default_stack_test_runner.py` and the CI step "Allele and reporter fixtures without ambient stack override". It is not in `df5c8244`. Run `38599805415` at `de9df8af` was **in progress** at preparation time: Linux and Windows running; macOS, the 08.04 synthetic GUI job, scroll and audit-CLI skipped. |
| The viewport repair is unverified | `5909a5dd` ("wrap promoter actions within visible viewport") is not in `df5c8244`. Its `promoter_action_rows_` tests run only in the opt-in `audit_followups` job (`.github/workflows/ci.yml` around lines 131-141), not in normal CI. |
| Last synthetic native 08.04 replay predates the viewport repair | Audit run `38476197008` at `024ed6c0` passed with explicit T, eight steps, three-`Seq` oracle parity and network isolation. Glen's physical-host GUI replay is policy-blocked (interactive `sudo unshare`). |
| The public refSNP pilot predates the dbSNP checks | Pilot `a1fb305a` (run `37778513095`) is not a descendant of `f14255cb` (dbSNP evidence validation). Its markers therefore cannot carry `dbsnp_reference_check` or `dbsnp_assembly_check`. |
| The public audit script does not assert the new checks | `scripts/prepare_public_vkorc1_audit.py` and `scripts/ci_vkorc1_gui_audit.py` never mention `dbsnp_reference_check`, `dbsnp_assembly_check`, `dbsnp_assembly_fallback`, `vcf_alt` or the placement. |
| The release aim is undecided | `docs/roadmap.md:24` reads "Active aim: undecided for `.13`... Owner selection is required". |

## Tasks

### 1. Correct the stale Windows claim

[`docs/vkorc1_allele_review_followup_20261008.md`](vkorc1_allele_review_followup_20261008.md)
says "The synthetic fetch tests use Unix-only mock-resource machinery; a green
native-Windows job does not prove those cases." This is wrong at `df5c8244`.

- Replace it with the job/test evidence above, bound to that SHA.
- Keep the rest of the paragraph's caution: unit tests are not native GUI or
  real-data acceptance.
- Apply the same correction in the
  [audit ledger](glen_audit_followups_20261008.md) wherever it is repeated.

### 2. Record the default-stack result from run `38599805415`

Wait for the run to finish on its own. Do not poll in a tight loop, re-run or
dispatch. Then, for Linux and Windows:

- Confirm the step "Allele and reporter fixtures without ambient stack
  override" executed every filter: `dbsnp_`, `materialize_variant_allele`,
  `reporter_construct_handoff_`, `variant_followup_allele_pair_` and
  `test_support::heavy_stack_tests`.
- Record pass/fail/ignored counts per filter.
- Confirm each runner JSON receipt shows `child_rust_min_stack: null`. A
  missing receipt is not a default-stack pass.
- Compare the counts with Glen's original Linux selections (11/2/7/2 plus two
  helper tests). Explain any difference.
- Retain original-log SHA-256s in a new evidence subset under
  `docs/audits/glen_followups_20261008/candidate_de9df8af/`, following the
  existing `candidate_df5c8244` pattern.
- macOS was skipped. Say so; do not infer a macOS default-stack pass.
- If anything failed, separate the primary failure from cascades (for example
  mutex poisoning). Report it as a typed product defect or environmental
  failure, and stop there. No speculative fixes.

If the steps pass, mark the default-stack item closed for Linux and Windows at
`de9df8af` in the follow-up doc, the ledger and `docs/roadmap.md:24`. Leave the
quick-xml `svg_dimensions_` item as is unless the same run shows it executed.

### 3. Make the public audit assert the dbSNP checks (Python only)

Extend `scripts/prepare_public_vkorc1_audit.py` so that a future public run
fails closed unless the fetched marker in the saved project state shows all of:

- `dbsnp_reference_check=match`;
- `dbsnp_assembly_check=match`;
- `dbsnp_assembly_fallback=false`;
- `dbsnp_variant_class=snv`;
- `vcf_ref=C`;
- `vcf_alt` containing `T`;
- a GRCh38-family assembly name;
- chromosome 16 (`NC_000016.*`) at 1-based position `31096368`.

Also:

- **Record, don't assume, the alternates.** Store the full fetched `vcf_alt`
  list and the raw `dbsnp_spdi` entries in the receipt. Do not hard-code
  `A,G,T`: dbSNP may change, and the receipt must show what was served.
- **Bare `alternate` refusal.** Keep the existing refusal step. It must refuse
  when more than one alternate is present. If the record ever has exactly one
  alternate, record that explicitly instead of failing on a missing refusal.
- **Unit tests.** Add tests to `scripts/test_prepare_public_vkorc1_audit.py`
  with synthetic project JSON for each failing qualifier (mismatch, fallback,
  non-SNV, wrong position) plus one passing case. Mark the fixtures as
  hand-crafted.
- **Run offline:**
  `python3 -m unittest scripts.test_prepare_public_vkorc1_audit scripts.test_ci_vkorc1_gui_audit -v`

### 4. Prepare (do not execute) the one owner dispatch

Write the exact request for the owner to approve, in the ledger's next-actions
section:

- workflow: `.github/workflows/ci.yml`, `workflow_dispatch`;
- inputs: `platform=linux`, `audit_followups=true`, `audit_public_vkorc1=true`;
- target: the final committed tip after Tasks 1-3, with its full SHA;
- precondition: normal CI at that SHA completed successfully.

This single run would cover:

- the live public dbSNP/GRCh38 replay with the new assertions;
- the synthetic native 08.04 replay on source that includes `5909a5dd`;
- the `promoter_action_rows_` and `promoter_pair_` viewport tests.

List the acceptance criteria per job and the artifacts to hash. State that it
needs network access to NCBI and Ensembl and that it makes no model call.

### 5. Draft the Glen physical-host checklist

Add a short section to
[`glen_final_candidate_execution_handoff.md`](glen_final_candidate_execution_handoff.md)
under "Remaining native and human checks". Glen runs it at the final tip, with
the packaged or dev build named explicitly, on physical macOS or Linux and on
native Windows. For each platform, record:

1. A multiallelic marker with an empty alternate field. Refused; project state,
   history and active view unchanged; no `_reference` output.
2. Retry with `T`. Exactly two outputs with unsuffixed IDs, identical length,
   one C→T difference at the marker, and an unchanged source.
3. The pair action and the analysis/export rows are reachable at default and
   minimum window size without horizontal scrolling.
4. Alleles are labelled genomic-forward. Nothing complements T into
   gene-strand or clinical nomenclature (−1639G>A is not entered or derived).
5. Platform, build profile, binary SHA-256, screenshots (with their hashes),
   and an explicit pass/fail per step.

Physical-host isolation must not be bypassed. If the policy block persists,
record it as environmental, not as a product failure.

### 6. Release acceptance: record inputs only

Do not pick or propose a `.13` aim. Add one line to the roadmap's open gates
listing the VKORC1 inputs available to the owner:

- default-stack result (Task 2);
- public replay (Task 4, pending dispatch);
- native checklist (Task 5, pending Glen).

Keep package/extracted-installer acceptance separate, as today.

## Deliverables And Handoff

- Edits to:
  - `docs/vkorc1_allele_review_followup_20261008.md`;
  - `docs/glen_audit_followups_20261008.md`;
  - `docs/glen_final_candidate_execution_handoff.md`;
  - `docs/roadmap.md`;
  - `docs/CHANGELOG.md` (completed items go here, never as roadmap "Done").
- `scripts/prepare_public_vkorc1_audit.py` and its tests.
- A new `candidate_de9df8af` evidence subset, if Task 2 has a completed run.
- Before handoff, run
  `scripts/maintenance_chore.py session-close --plan docs/vkorc1_acceptance_codex_prompt_20261014.md`
  and report its warnings.
- The handoff note lists the files changed, tests run (with counts), unrun
  checks and remaining risks.
- End with a typed verdict per open item: closed, pending owner dispatch,
  pending Glen, or owner decision.
- Commit locally on a branch if the owner wants that; do not push.
