# Request to Glen: Gene Assay Study GUI and Tutorial Screenshots

Please test the [04.08 tutorial](tutorial/04-08_gene_assay_study_gui.md) at the
exact merged candidate SHA, using matching GUI and CLI binaries. This is .11
primer-PAIR study work, not transcript capture and not an additional .10 gate.
Do not change design constraints or fixtures to obtain a green result.

The implementation was committed at the user's request on `gentle_rs_2_main`,
based on `8597a6ee35cd484e0452dde358742347a4b0cb0b`; that base SHA alone does
**not** contain this GUI. Use the exact implementation commit, or a later
explicitly selected merged candidate, for acceptance. A local commit is not
evidence that it has been pushed or merged. The user delegated final Rust verification
to you after a memory-constrained local rebuild. See the
[verification handoff](gene_assay_study_gui_plan.md#verification-handoff),
including the initially failing stale-plan tests and the not-yet-rerun fix.

Before live capture, run the focused implementation checks from that checkout:

```bash
cargo test --locked --lib gene_assay_study_ui
cargo test --locked --lib command_execution::tests
cargo check -q --locked
cargo check --locked --features gui-test-support --lib
cargo fmt --check
git diff --check
python3 -m unittest discover -s scripts -p 'test_*tutorial*.py'
```

Use the tutorial helper with your exact built CLI for the real replay. Also run
the usual tutorial catalog/manifest and generated-document checks; a manual
guide's catalogue entry alone is not execution or screenshot acceptance.
In particular rerun `gentle_examples_docs tutorial-check` after the review-entry
correction; do not accept the earlier `report.json` mismatch as a waived gate.

Please annotate the tutorial with screenshots at these checkpoints:

| ID | Capture | Scientific/UI check |
| --- | --- | --- |
| G1 | Splicing Expert launch and saved study header | Gene, source feature, transcript scope, recommendation/override, missing evidence; no implied execution approval |
| G2 | One selected pair and its transcript outcomes | Both primers 5' to 3', mature-cDNA coordinates on minus locus, shared/no-product/missing distinctions |
| G3 | Experimental readiness card | Specificity absent means blocked; candidate oligos are not an approved order |
| G4 | Effective request and new plan | Changed maximum, new plan ID/output directory, first review; old files preserved |
| G5 | Second review and execution/cancellation | Exact workflow digest; no design before approval; tamper/stale/cancelled outcomes cannot become success |
| G6 | Canonical pending dossier | Same scientific content as CLI; no comparison panel falsely attributed to the plan |

Retain untouched raw captures separately from cropped/annotated versions. Record
source SHA, binary hashes, OS/profile, input/preparation receipt hashes, action,
effective parameters, expected and observed result, and image hashes. Public
screenshots must use the synthetic fixture only. Capture only passed checkpoints
for tutorial illustrations; failed checkpoints need separate diagnostic evidence,
not edited-to-look-passing illustrations. Add concise captions for wet-lab readers,
and keep previous image revisions recoverable.

Semantic controls are `splicing.gene_assay_study.open`,
`splicing.transcript_panel.open`, and PCR Designer's `assay_study.open_plan`,
`assay_study.open_request`, `assay_study.normalize`, `assay_study.plan`,
`assay_study.execute`, `assay_study.publish`. Use a `gui-test-support` build for
those controls. File dialogs, review checkboxes, pair picker and remaining
controls are manual/hybrid checkpoints, not claims of complete scripted coverage.

For paper-facing acceptance, please supply one original gene-assay source bundle:
the typed publication request, exact study plan and request, evidence ledger,
input hashes/reference/annotation identities, persisted panels, handoffs and any
reviewed order forms. No such original bundle was found in the public fixtures.
Do not reconstruct it from a PDF. Keep private real-data receipts external and
checksum-bound; distinguish product defects, missing/incompatible input,
environment problems and needs-human-review verdicts.

Please return the exact tested revision and G1-G6 results before publication.
Tutorial captures are evidence of GUI behaviour, not approval to order or proof
of biological/experimental suitability. This file is a forwardable request;
writing it does not mean Glen has already been contacted or accepted the work.
