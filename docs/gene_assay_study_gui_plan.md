# Gene Assay Study Workspace

Status: implementation slices for internal .11, 2026-09-18. Not Claude-reviewed;
an optional read-only consultation was offered. Work continues on
`gentle_rs_2_main`, fast-forwarded to local `main` at `8597a6ee`. No .10 release
action, tag or publication is part of this work.

## Existing Contracts

At the audit baseline, Splicing Expert already opened PCR Designer's transcript-panel mode. That mode
uses `DesignTranscriptAssayPanel`, renders the saved product matrix, and builds
the shared experimental handoff. It was not a gene-study workspace: the
`PlanGeneIsoformAssayStudy` request/plan and canonical publication request
required shell continuation. The panel candidate table's historical
`order_ready_primers` field does not establish procurement readiness.

Use the existing PCR Designer rather than another designer. The added study
mode must inspect typed source JSON, never reverse-engineer publication PDFs.
The transcript-capture pool is a different programme, not this primer-pair aim.

## Bounded Slices

1. Add a Gene assay study entry in Splicing Expert and PCR Designer. Inspect
   a saved plan's question, transcript scope, annotation identity, declared
   evidence, missing assessments, recommendation, override and exact operations.
   Link matching persisted panels to a selected-pair inspector showing stored
   binding footprints, per-transcript products and selection evidence.
2. Normalize a supplied typed request through the engine; explicitly review it
   before planning. Allow a short-product ceiling edit, retain full JSON for
   expert review, and invalidate review on edits. Execute only a separately
   reviewed exact plan/workflow through managed background commands, with
   cancellation, stale-result protection and fresh output paths. Saved-plan
   inspection by itself grants no execution approval.
3. Expose the existing canonical dossier export from an explicitly supplied
   publication request. Retain its content/digest and order-readiness checks;
   do not generate an independent GUI narrative or submit orders.
4. Add one catalogued PATZ1 tutorial and a forwardable screenshot/acceptance
   request for Glen. Reuse the existing synthetic minus-strand fixture, internal
   backend, and declared evidence. Separate automated checks from live GUI and
   human scientific acceptance.

## Acceptance Boundaries

Test request invalidation, foreign-locus refusal, exact command execution,
pair-to-transcript presentation, stale/cancelled results and canonical export.
Run focused tests and Cargo check before handoff. Retain existing digest checks;
missing evidence is not a zero or pass. Do not change scientific fixtures to
force a successful design.

No original paper-facing gene-assay source bundle was found among the repository's
public examples/fixtures. Glen must supply an original typed publication request
and its digest-bound study plans, evidence, panels and handoffs for that separate
acceptance; published images cannot substitute for the missing inputs.

## Verification Handoff

Follow-up, 2026-09-19: the previously interrupted ten study tests and eleven
`command_execution::tests` all passed locally at `8284c4ca` using
`cargo test --locked --offline -j 1 --lib FILTER`. This closes the old
unverified receipt/stale-result fix, not Glen's live or real-data acceptance.
The earlier generated tutorial mismatch was also closed by the recorded TSS
rebase verification in the changelog; it is not a waived failure.

The current user-requested continuation is a presentation-only completion:
replace matrix/band truncation with independent, bounded pages; page the
selected-pair outcomes too; expose the engine handoff's existing coverage
accounting in both panel and study modes. Do not alter scientific selection,
the coverage universe, report schemas or exported rows. The new headless
navigation regression uses explicitly synthetic display records (95 transcripts,
25 assays, 81 band rows), not evidence that those are feasible assays. This
scope was offered for optional Claude consultation but was not Claude-reviewed.

Verification of this working-tree continuation: 60 `transcript_assay_` tests
(including four new coverage/navigation tests), ten `gene_assay_study_ui` tests
and eleven `command_execution::tests` passed. The large-panel test clicks the
live headless egui page controls and verifies the last transcript, assay and
band row, missing-cell wording and unchanged report JSON. Python tutorial
discovery ran 88 tests with 12 declared optional skips. Default offline locked
Cargo check, the `gui-test-support --lib` check, formatting and whitespace
checks also passed. Session-close hygiene reported only the intentional dirty
tree and manual plan-fidelity reminder. These results do not replace native
GUI or real-reference acceptance at a frozen candidate SHA. The full workspace
and generated tutorial replay were not rerun for this presentation-only slice.

Implementation files: `src/main_area_dna/transcript_assay_report_ui.rs`, its
module declaration in `src/main_area_dna.rs`, and the existing
`gene_assay_study_ui.rs` / `primer_design_ui.rs` views. Documentation changes
are limited to this plan, the GUI guide, tutorial 04.08, Glen's acceptance
request and the roadmap/changelog status.

### Original Handoff Record

The implementation was prepared on the baseline above and is now committed at
the user's request. Acceptance must name the exact commit containing this work,
not the baseline SHA alone.
Final Rust verification was delegated to Glen at the user's explicit request
after the focused rebuild spent almost 30 minutes competing with other builds.
Only this task's rebuild and queued semantic-feature check were terminated.

- The initial seven Rust tests passed four cases and exposed three failures:
  planning advanced the project revision and its own result was rejected as
  stale. Receipts now capture the committed structural revision, and result
  attachment reads it against the project under one read lock. The final ten
  tests, including post-commit edits, reused IDs and canonical export, required
  the later completed run recorded above.
- An intermediate default `cargo check -q --locked --offline -j1` passed before
  the final receipt/cache changes. Final default and semantic-feature checks
  remain required. Formatting and whitespace checks passed.
- Python tutorial discovery ran 86 tests successfully with 12 optional skips.
  Both helper tests also passed with `GENTLE_TUTORIAL_BIN_DIR=target/debug`,
  including the real offline preparation/export replay using the pre-existing
  CLI. That replay is not proof of the new GUI or a newly built binary.
- Catalog and manifest checks passed (57 entries, 28 generated chapters).
  The full generated tutorial check reported `report.json` drift after the new
  tutorial was added without a review-manifest entry. An explicitly unreviewed
  entry was added and the catalog regenerated; the later TSS rebase verification
  ran the full check successfully without waiving the earlier failure.
- No live GUI screenshot, inner-agent run, real-data study, order approval or
  publication was performed. Pre-commit session-close hygiene reported only
  the intentional dirty/generated-catalog and manual plan-fidelity warnings.
