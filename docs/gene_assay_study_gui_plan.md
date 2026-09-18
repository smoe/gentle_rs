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
  tests, including post-commit edits, reused IDs and canonical export, still
  require a completed rebuild/run. No passing final Rust result is claimed.
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
  entry was added and the catalog regenerated; the costly full check needs a
  rerun, not a presumed pass or unrelated artifact refresh.
- No live GUI screenshot, inner-agent run, real-data study, order approval or
  publication was performed. Pre-commit session-close hygiene reported only
  the intentional dirty/generated-catalog and manual plan-fidelity warnings.
