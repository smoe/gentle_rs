# Glen: Scrutinize the Microarray and Measured-Gel Tutorials

This is an independent review request, **not acceptance**. The owner asked for
an immediately mergeable first tutorial slice. Preparation baseline is local
`150a339074e7a2d282567d791a8f50834b70ff9c`; audit the exact commit containing
this note, record its full SHA, and report any later drift. No direct message
has been sent to Glen; Steffen can forward this document.

The slice contains navigation hygiene plus manual guides 07.04 and 10.03.
It changes no engine/GUI biology. Local Python source/fixture/checkout checks
are not an executed GENtle workflow or screenshot acceptance. New repeat and
imported-DuckDB execution tests, and FLNA public inputs, belong to later slices;
they are deliberately not prerequisites for merging these two manual guides.

Local preparation on macOS: 79 Python tests pass; catalog/manifest freshness
checks pass using the existing `150a3390` generator binary (SHA-256
`906131ff6153c1cb76b2bbc6b5f56819a410b9caaa538e63776a81362a18b5e2`).
Catalog: 64 entries; executable manifest: unchanged 30 chapters; GUI contracts:
unchanged six. No new Rust build/check or learner-command execution is claimed.

## Common Review Contract

Please use a clean isolated worktree and one exact-source CLI/GUI pair. Keep
your private study untouched. Return source SHA, binary hashes/self-identities,
Cargo.lock hash, toolchain, effective profile/features, OS/display environment,
exact commands and per-check outcomes. No push, tag, publication or replacement
of old tutorial evidence is requested. A build/test pass is not scientific
approval or a numerical performance verdict.

First read the learner prose as an experienced wet-lab biologist. Check every
command against the shared parser, then actually execute it. Report product
defect, missing/incompatible data, environment problem and human-scientific-
review-needed separately. Missing evidence is not zero; known warnings remain
visible. Keep original raw screenshots, semantic snapshots and hashes separate
from annotated teaching copies. Do not claim macOS/Windows from a Linux run.

## 07.04: Clariom D on TP73

Read `docs/tutorial/07-04_microarray_evidence_gui_cli.md` and the existing
`test_files/fixtures/microarray_tracks/README.md`. Check real vendor IDs and
coordinates versus **synthetic** contrast statistics. Do not use private CEL
files or describe this as E-MTAB-14704 expression acceptance.

```sh
cargo test --locked --lib project_microarray_track_uses_vendor_subset_on_tp73_genbank_anchor
```

Execute the lesson's inspect/open/project commands in one fresh state file.
Expect five imported contrast rows, not five independent probesets; verify
alpha PSR0100145779.hg.1 spans genomic 3652520..3652564 and local `[4,49)`.
Review both contrast values and the meaning of the cutoffs.

In the native GUI, capture Configuration -> Microarrays with file-role guidance,
the anchored TP73 Array track with readable hover/summary, and nearest-feature
navigation from a region with no visible feature. Verify that nearest means
nearest **projected** feature, not a whole-genome vendor lookup. Try the existing
hg19-projected synthetic manifest only in a disposable copy: distinguish
unsupported mapping from a supported target with zero TP73 overlap. Preserve
assembly warnings. Save/reopen and verify features/contrasts remain intact.

## 10.03: Measured Image Calibration

Read `docs/tutorial/10-03_measured_gel_sizing_gui_cli.md` and
`docs/tutorial/inputs/gel_sizing/README.md`. The PNG is drawn, not experimental.

```sh
python3 docs/tutorial/inputs/make_gel_sizing_fixture.py docs/tutorial/inputs/gel_sizing --check
cargo test --locked --lib gel_image
```

Execute import, both analyze requests, inspect and SVG/TSV exports in one fresh
state. Compare the import digest with the PNG's actual SHA-256 and requests.
DNA inside band: approximately 707.106781 bp, not 750. SDS appendix: approximately
35.355339 apparent kDa; invented prestained ladder/gel-system labels are not a
commercial calibration. Both outside bands must remain null with
`outside_calibrated_range`, never extrapolated or treated as zero. Localization
bounds are not statistical confidence intervals or amount measurements.

Native screenshots: original image with marked lanes/ladder; result table with
inside and outside outcomes; optional apparent-kDa appendix. Check zoom/scroll
does not change original zero-based marks, no-clobber export, and project
save/reopen retains original image/digest, assignments and reports. Select real
GUI-assigned IDs instead of forcing CLI fixture IDs into an unrelated image.
Automatic band detection is out of scope.

## Catalog and Return

```sh
python3 -m unittest scripts.test_tutorial_gap_lessons scripts.test_tutorial_checkouts scripts.test_tutorial_acceptance
gentle_examples_docs tutorial-catalog-check
gentle_examples_docs tutorial-manifest-check
gentle_examples_docs tutorial-check
cargo check -q --locked
cargo fmt --check
git diff --check
```

Confirm the new guides are discoverable manual/hybrid entries, not fabricated
typed GUI contracts. Existing executable chapters, old receipts and historical
human review dates must remain unchanged. The README no longer duplicates
review badges; generated catalog/Help still obtain reviews from the manifest.

Reply through Steffen with a concise finding list, exact SHA and artifact paths,
per-check pass/fail/not-run/blocked, screenshots/hash manifest and any proposed
wording corrections. Keep native platforms and human scientific review pending
unless actually performed; do not convert this review into release approval.
