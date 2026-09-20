# Request to Glen: Gene Assay Study GUI and Tutorial Screenshots

**Result:** the six public-reference checkpoints were completed manually at the
exact `.11` code candidate `3c1c32bc`; see the
[2026-09-20 acceptance report](gene_assay_study_glen_acceptance_20260920.md)
and its [hash-bound screenshot evidence](screenshots/gene_assay_study_gui/README.md).
This records GUI presentation, not biological validation or order approval.

Please test the [04.08 tutorial](tutorial/04-08_gene_assay_study_gui.md) at the
exact merged candidate SHA, using matching GUI and CLI binaries. This is .11
primer-PAIR study work, not transcript capture and not an additional .10 gate.
Do not change design constraints or fixtures to obtain a green result.

The implementation was committed at the user's request on `gentle_rs_2_main`,
based on `8597a6ee35cd484e0452dde358742347a4b0cb0b`; that base SHA alone does
**not** contain this GUI. Use the exact implementation commit, or a later
explicitly selected merged candidate, for acceptance. A local commit is not
evidence that it has been pushed or merged. The previously interrupted ten
study tests and eleven managed-command tests now pass locally at `8284c4ca`.
This does not replace your exact-candidate rerun, live capture or real-data
acceptance. See the [verification handoff](gene_assay_study_gui_plan.md#verification-handoff).

Before live capture, run the focused implementation checks from that checkout:

```bash
cargo test --locked --lib gene_assay_study_ui
cargo test --locked --lib transcript_assay_report_ui
cargo test --locked --lib command_execution::tests
cargo check -q --locked
cargo check --locked --features gui-test-support --lib
cargo fmt --check
git diff --check
python3 -m unittest discover -s scripts -p 'test_*tutorial*.py'
```

Use `scripts/prepare_real_patz1_tutorial.py` with your exact built CLI for the
authentic-reference replay. The pinned public Ensembl 116/RefSeq files replace
the synthetic example for this acceptance. Also run
the usual tutorial catalog/manifest and generated-document checks; a manual
guide's catalogue entry alone is not execution or screenshot acceptance.
In particular rerun `gentle_examples_docs tutorial-check` after the review-entry
correction; do not accept the earlier `report.json` mismatch as a waived gate.

Please annotate the tutorial with screenshots at these checkpoints:

| ID | Capture | Scientific/UI check |
| --- | --- | --- |
| G1 | Locus figure, DNA-map and Structure source-comparison panels | 13 Ensembl and four RefSeq records; exact shared full chain ENST00000266269.10 / NM_014323.3; versions, hashes and distinct CDS rows retained; minus-strand coordinates not reversed twice; shared filter/highlight, DNA viewport synchronization, stale-binding rejection; comparison navigation leaves 13 actionable design targets and does not assign PCR selection |
| G2 | Transcript-panel request and scope | All 13 loaded Ensembl transcripts, minimal-discrimination objective, explicit best-effort policy; RefSeq comparison is not silently added to the design universe |
| G3 | Both primers, product matrix and unresolved distinctions | Both primers 5' to 3', mature-cDNA coordinates, coverage versus discrimination, shared/no-product/missing distinctions; a partial panel is not labelled a complete solution |
| G4 | Separate effective study request and new plan | No synthetic expression evidence; new plan ID/output directory, first and second review; exploratory panel not attributed to an unexecuted study |
| G5 | Experimental readiness card | Missing authentic genomic and whole-transcriptome specificity blocks readiness; no design before approval; tamper/stale/cancelled outcomes cannot become success |
| G6 | Canonical pending dossier | Same scientific content as CLI; no comparison panel falsely attributed to the plan |

The development replay (Primer3 2.6.1) yielded seven pairs covering 13 exact-cDNA
classes but nine unresolved class-pair distinctions, correctly `partial`. This
is an observation, not a prescribed pass count or order approval. Please rerun
at the exact candidate and retain actual output, including any unresolved cases.

At G2/G3 also review the coverage-scope explanation, distinguishing transcript
records from exact-cDNA groups, uncovered from unassessed, and coverage from
isoform discrimination. For a real panel exceeding 40 transcripts, 24 assays
or 80 band rows, verify the independent page controls reach the final entries
without changing the report or export. This 13-transcript authentic fixture
does not establish performance at those larger sizes.

Retain untouched raw captures separately from cropped/annotated versions. Record
source SHA, binary hashes, OS/profile, input/preparation receipt hashes, action,
effective parameters, expected and observed result, and image hashes. Public
screenshots must use the pinned public reference fixture, not private samples.
Capture only passed checkpoints
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
