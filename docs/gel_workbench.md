# Measured-Gel Workbench

Development implementation, not native/scientific acceptance. This extends the
[manual sizing editor](gel_image_analysis.md), not the simulated gel renderer.
One engine operation (`GelWorkbench`) and shared command (`gel-image workbench`)
serve GUI, CLI, MCP `op` and scripting. Keep one saved project to compare any
explicitly selected lanes across its imported gels; no adjacent-well or cell-line
pairing is guessed. Save Project retains sources, calibrations, profiles and
immutable comparisons. Current limits: 64 profiles, 64 comparisons, 20,000 bins
per profile/common grid and eight analyses/eight comparisons per report export.

## Import Sources Without Inventing Raw Measurements

The normal image importer still accepts PNG/JPEG/TIFF. The workbench additionally
extracts an original embedded raster from a supported PDF, retaining both the
original document and extracted image with separate SHA-256 values. It does not
render a screenshot, execute PDF actions, digitize the electropherogram or read
the Peak/Region Tables as a distribution. Multiple supported rasters require
an explicit `pdf_object` number rather than guessing which one is a gel.

PDF admission is deliberately bounded: single-revision classic xref tables,
unencrypted sources, supported 8-bit grayscale/RGB raster/filter combinations,
no unsupported masks/transforms, at most 16 MiB, 32 pages and 4,096 objects.
Unsupported PDFs fail explicitly instead of silently changing their image.
The existing public TapeStation PDFs exercise embedded image object 5.

Native `.D5000` **numerical decoding is not implemented**. The supplied archive's
measurement/model XML members are encrypted. Explicit `allow_native_thumbnail`
imports only its available `Thumbnail.png`, with a persistent display-only
warning, not the numerical measurement or operational logs. Without that opt-in
the import refuses. Do not infer assay calibration, samples or fluorescence
units from this preview. Request exported numerical traces from the laboratory;
Agilent documents CSV/XML export in its
[software description](https://www.agilent.com/en/p/electrophoresis/software-laptops/tapestation).

```sh
gentle_cli --state gels.json gel-image workbench '{"action":"import_source","image_id":"before","path":"before.pdf","pdf_object":5}'
gentle_cli --state gels.json gel-image workbench '{"action":"import_source","image_id":"after","path":"after.pdf","pdf_object":5}'
gentle_cli --state gels.json gel-image workbench '{"action":"import_source","image_id":"native-preview","path":"run.D5000","allow_native_thumbnail":true}'
```

These are local, mutating imports. They do not upload a source or establish
permission to publish its embedded identifiers or sample data.

## Propose Bands, Then Review

Draw the lanes on original pixels. Set polarity (dark/light bands), smoothing,
minimum prominence and minimum spacing. `detect_bands` averages across each
specified lane and finds local peaks. Its returned marks are **unconfirmed
proposals**. Accepting proposals into the editable sample assignments is a
separate user action. Ladder identities/sizes are never auto-assigned.
Display extraction normalizes native 8/16-bit channels and uses the fixed
stored-channel weights `(2126 R + 7152 G + 722 B) / 10000`
(`stored_rgb_rec709_v1`), without implicit color-space/transfer conversion.
Native alpha is composited over white (`alpha_white_v1`). Both policies are
recorded in provenance; original bytes are not changed. Non-finite/out-of-range
floating-point channels are rejected, not clipped into display measurements.
Transparent previews use the same white background, rather than displaying
hidden RGB values after discarding the PDF/raster alpha channel.
Smears, labels, saturation and shoulders need review; a detected image maximum
is not a specific molecular product or concentration measurement.

The GUI uses the same `GelDetectRequest` as CLI:

```json
{"action":"detect_bands","image_id":"before","migration":"down",
 "dark_bands":true,"minimum_prominence":0.08,"minimum_distance_px":5,
 "smoothing_radius_px":1,
 "lanes":[{"id":"B1","label":"Cell 1 custom","min":{"x":100,"y":40},"max":{"x":125,"y":650}}]}
```

The geometry above is only a request-shape example, not a confirmed assignment
for the public TapeStation gel. Use its retained lane/calibration geometry or
review the actual image. Existing `gel-image analyze` validates accepted marks
and confirmed ladder calibration normally; detection alone creates no sizing
report and does not modify the project.

## Retain Profiles

`extract_profile` requires a **saved valid calibration**, a chosen sample lane,
explicit background and polarity. It integrates original image display intensity
over each calibrated migration-pixel interval. The calibrated size interval
retains that weight; it is not falsely integrated again as density per bp.
No whole-lane extrapolation, background fit or molecular-weight conversion is
invented. The resulting `image_display` profile is a contrast-dependent shape
diagnostic, not fluorescence, DNA mass, molarity, molecular counts or yield.
Restrict alignment markers explicitly; optional `masked_ranges` lists reviewed
size intervals and makes each intersecting pixel bin unavailable. No automatic
marker recognition is claimed. A retained source digest distinguishes the same
lane name in two runs.

```sh
gentle_cli --state gels.json gel-image workbench '{"action":"extract_profile","profile_id":"before-custom-1","analysis_id":"before_A1","lane_id":"B1","run_id":"2026-10-06","dark_bands":true,"background":0}'
```

Numeric imports are separately `vendor_trace` evidence. `import_profile` accepts
a complete explicit `GelProfile` in request JSON; source identity and calibration
there are caller-declared, not independently verified detector evidence.
Invented teaching inputs instead use explicit `synthetic` evidence; they cannot
be pooled with vendor/image profiles by merely matching display units.
`import_profile_file` additionally retains exact original CSV bytes/digest.
Its strict, documented **integrated-bin** CSV has these three columns:

```csv
left_edge,right_edge,weight
100,200,2.5
200,300,
400,500,1.2
```

This is synthetic. The blank bin and the absent 300..400 interval are missing,
not zero. Raw point-sampled FU-vs-time/FU-vs-bp exports cannot be relabelled as
this format: their integration measure/calibration must first be established.
No vendor CSV layout or Jacobian is guessed. `size_kind`, `signal_units` and
`measure` are explicit, reviewed caller inputs:

```sh
gentle_cli --state gels.json gel-image workbench '{"action":"import_profile_file","path":"trace.csv","profile_id":"trace1","label":"Calibrated bin areas","run_id":"run1","lane_id":"B1","size_kind":"linear_dna_bp","signal_units":"vendor_integrated_area","measure":"reviewed_vendor_area_bins"}'
```

## Compare Any Compatible Pair

Choose two stored profiles, an explicit size range and optionally a short-size
threshold. Units, size kind, integration measure and evidence kind must match;
image brightness is never pooled with vendor/calibrated signal. GENtle
area-conservingly rebins onto the union of their edges, assuming uniform weight
**within** each recorded size bin. Intervals missing in either operand are
excluded from both; there is no interpolation across missing bins or markers.
Each distribution is independently unit-area normalized on that common valid
support. Exclusions, original profiles and unnormalized totals remain in the
immutable report. Interpret a fraction/quantile as conditional on that support.

```sh
gentle_cli --state gels.json gel-image workbench '{"action":"compare","comparison_id":"cleanup-cell1","left_profile":"before-custom-1","right_profile":"after-custom-1","range":[100,5000],"short_threshold":500}'
gentle_cli --state gels.json gel-image workbench '{"action":"inspect"}'
gentle_cli --state gels.json gel-image workbench '{"action":"inspect_comparison","comparison_id":"cleanup-cell1"}'
```

The result includes q10/median/q90, the fraction below the chosen threshold,
99 weighted QQ pairs, the maximum cumulative-distribution gap (`ks_distance`)
and Wasserstein-1 distance in bp or kDa. No positive common weight means no
distribution, not a row of zeros. Signed algebra results cannot be probability
distributions.

### What A QQ Plot Says

The x-axis is the left profile's fragment-size quantile and y the right one's
quantile at the **same cumulative weight**, not sorted brightness pixels or
paired peak indices. Identical shapes lie on y=x. With before on x and after
on y, a curve above the diagonal indicates larger after quantiles. Deviation
only near the lower tail is compatible with short-fragment removal; a broadly
shifted diagonal may also reflect calibration. Shape agreement does not show
equal concentrations or recovery. Mass-weighted, molar-weighted and
display-intensity-weighted quantiles answer different questions.

### KS Distance Is Not An Automatic KS Significance Test

The classical two-sample KS test compares distributions of **independent
observations**; see the
[SciPy contract](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html).
Smoothed neighboring detector bins and gel-image pixels are not independent
fragment observations. Their number depends on scanning/export resolution,
not the number of molecules or biological replicates. Therefore the report
retains `p_value: null`. A larger image or denser export must not manufacture
stronger biological significance. `ks_distance` here is a descriptive maximum
gap between weighted cumulative curves, with no rejection threshold.

The complementary
[Wasserstein distance](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wasserstein_distance.html)
describes how far signal weight must move on the size axis; unlike the KS gap
it carries the size unit. Quantile shifts and below-threshold fractions are
often the most direct laboratory summaries. Later inferential work should use
independent library/biological replicates and a declared paired permutation or
replicate-bootstrap design. Tail-sensitive Anderson-Darling/Cramer-von Mises
comparisons also need that sampling model; they do not fix pseudoreplication.
Multiple pairwise comparisons additionally need a prespecified multiplicity
policy. These inferential tests are not implemented in this slice.

## Signed Lane Algebra

Explicit linear combinations use the same valid shared size support and
declared units/measure, **without automatic unit-area normalization or fitting**.
Missing active operands stay missing; negative residuals are retained.

```sh
gentle_cli --state gels.json gel-image workbench '{"action":"algebra","output_id":"cell1-residual","range":[100,5000],"terms":[{"profile_id":"mixed-cell1","coefficient":1},{"profile_id":"manufacturer-cell1","coefficient":-1}]}'
```

Choose coefficients from reviewed loading/normalization information. A residual
near a custom-only band is exploratory support, not a gene assignment or proof
that mixed reverse transcription equals a physical sum of separate reactions.
Ratios, fitted mixture models and replicate inference remain separate follow-ups.

### Verify And Display A Saved Residual

`inspect_algebra` is read-only. It strictly parses the recorded v1 request and
replays it against saved operands, including nested residuals, before exposing a
plot series. It does not fetch data, redetect bands, re-extract lanes, fit scales
or change project history. The statuses are `verified`, `mismatch` and
`unverifiable`; verification means **internal consistency**, not authenticity or
biological acceptance. Original/imported non-residual profiles still carry their
declared evidence limitations.

```sh
gentle_cli --state gels.json gel-image workbench '{"action":"inspect_algebra","profile_id":"cell1-residual"}'
```

The GUI's **Verify and display saved algebra** uses the same operation in a
background worker. That worker also rasterizes the shared SVG with GENtle's
configured font path; the GUI only installs/caches the prepared texture. Idle
frames do not decode SVG, replay or hash operands. Project revision/selection changes, cancellation,
closing and project replacement discard stale inspections. Reopen and inspect
again after an edit.

The standalone graph shows **saved bin integral divided by bin width**, labelled
with the original signal units **per bp** (or **per kDa**). The area under a step
therefore retains the signed integral even on unequal-width bins. This display
transform does not change stored weights or imply a calibrated density of DNA
mass/molecules. The zero line is explicit, negative values remain negative,
and hatched intervals remain unavailable, not zero. All-zero and all-missing
results have different messages. Neither defines a new inferential test.

The derived result is identified by its profile ID and **Algebra input digest**,
not by legacy run/lane fields inherited from its first operand. Its presentation
retains every dependency's identity, evidence kind, units, source/input digest,
processing/calibration and warnings. Zero-coefficient operands bind provenance
and grid edges, but contribute neither signal nor missing-value exclusions.
Different ladder calibrations need alignment review; subtraction does not fix
registration uncertainty.

Replay is bounded to depth 16, 64 distinct profiles and 1,280,000 bin visits,
with memoization of shared dependencies and cycle detection. Exceeding a budget
is unavailable verification, not a successful zero result. The original algebra
v1 calculation, digest field order and receipt bytes are frozen. Its 2,048-byte
receipt limit still applies; long requests fail explicitly. Old saved residuals
remain readable and can verify without migration. New `import_profile` requests
cannot claim `signed_residual`; use the algebra operation. Legacy/unverifiable
records remain inspectable, but **report export refuses them**.

Reference overlays, fitted/manual display factors and changes to the existing
comparison-page bin-total plot are deferred. The latter is labelled NOT a density;
its bin-width-sensitive appearance is a separate display follow-up, not changed
scientific weights or a silent revision of existing report pages.

## Export The Saved Before/After Report

Select saved analyses in explicit before/after order, plus saved comparison IDs.
The GUI PDF export and CLI use the same operation. The PDF has a paired overview,
the saved image/calibration details, and profile/CDF/QQ panels with units,
source identities, limitations and descriptive metrics. It does not redetect
bands, re-extract profiles, refetch sources or pick new comparison settings.
Integrity validation checks saved results against their frozen inputs.
Source validation also replays the recorded PDF/thumbnail extraction and
integrated-bin CSV parsing to verify source-to-result binding; this does not
redetect bands or extract a new quantitative lane profile. Original document
and numeric-source summaries retain full printable hashes without raw payloads.
SVG hover remains available in existing per-analysis SVG exports; PDF is static.

```sh
gentle_cli --state gels.json gel-image workbench '{"action":"export","path":"gel-comparison.pdf","title":"Before and after cleanup","analysis_ids":["before_A1","after_L1"],"comparison_ids":["cleanup-cell1"]}'
gentle_cli --state gels.json gel-image workbench '{"action":"export","path":"gel-algebra.pdf","title":"Exploratory signed residual","analysis_ids":[],"algebra_ids":["cell1-residual"]}'
```

Use `.json` instead for complete portable machine-readable results. Exports
never overwrite existing files and are external writes, not undoable project
edits. The new report does not append the manufacturer's complete PDFs or copy
its Region Tables automatically; originals remain retained in the project.
An unavailable encrypted raw trace stays unavailable in every presentation.
`algebra_ids` is optional and ordered (at most eight, like each other record
kind). Selected results must verify again; no stale GUI cache authorizes export.
Residual pages and complete wrapped operand/warning appendices follow the
existing pages. Algebra-only exports need no empty before/after overview.
JSON adds `algebra` only when results are selected, with the prepared series and
all saved dependency snapshots. Requests omitting `algebra_ids` retain their
old serialization and existing report pages. There are no residual CDF/QQ/KS,
Wasserstein, quantile or significance summaries.

Algebra export is bounded to 1,280,000 retained input/result bins in total,
64 provenance pages per result and 128 total PDF pages. Oversized metadata is
refused rather than silently clipped; use fewer results or JSON. The new pages
do not assign public cell-line subtraction coefficients: `1/-1` in the
[teaching example](tutorial/10-03_measured_gel_sizing_gui_cli.md#optional-synthetic-signed-algebra)
is illustrative only.

## Verification Boundary

Deterministic Rust regressions cover synthetic band proposals, conservation of
integrated weights, explicit gaps, identity/shifted QQ, analytically known
KS/Wasserstein distances, signed subtraction, CSV LF/CRLF identity, source
import refusal/retention and shared shell/storage. Renderer tests cover plot
semantics and bounded escaped output. The public PDF regression imports both
original documents, uses the existing reviewed assignments, extracts before B1
and after C1 and reopens the serialized project before exporting five pages.
It verifies invoked PDF text/font mappings, unchanged state and no-clobber.

Owner-authorized local verification on 2026-10-19: base
`03e1b80041163ea9044cfc263a59307517bd948e` **plus the uncommitted workbench changes**,
Darwin 27 arm64, Rust `1.100.0-beta.5 (8ea38001b 2026-10-13)` and Cargo
`1.100.0-beta.5 (7552a885e 2026-09-29)`. Default `desktop-gui` and
`screenshot-capture` features were compiled; no screenshot capture was invoked.

| Check | Result |
| --- | --- |
| Entire `gentle-protocol` suite | 90 passed, including four new workbench contracts and canonical parity freshness |
| `gentle-render` / `gel_workbench` | 15 passed |
| Root library / `gel_workbench` | 48 passed; one opt-in visual writer ignored |
| Root library / `gel_image` | 41 passed, including the unchanged real sizing, GUI-adapter and MCP-consent paths |
| Explicit opt-in PDF writer | One passed; five-page report rendered with Poppler and visually inspected |
| Glossary parser and CLI-documentation checks | One passed in each filter |
| `cargo check -q --locked --offline` | Passed with the default dev profile |
| Python fixture and LF/CRLF checkout suite | 79 tests, one explicit existing-CLI-binary skip |

Filters overlap; their counts must not be added as independent coverage. The
check uses the normal dev profile. Tests only disable debug symbols via
`CARGO_PROFILE_TEST_DEBUG=0`; optimization and assertions are unchanged. Builds
use one Cargo job. No production profile, dependency version or fixture hash
was altered by the verification repairs.
Formatting, documentation status and whitespace checks pass. Session-close
reports four OK, two expected warnings (intentional dirty worktree and manual
plan-fidelity reminder), and no failures; the deferred subset is explicit here.

Reproduce the targeted commands from the repository root after fetching the
locked dependencies:

```bash
export CARGO_BUILD_JOBS=1
export CARGO_PROFILE_TEST_DEBUG=0
# macOS Swift bridge cache; not needed on Linux/Windows.
export CLANG_MODULE_CACHE_PATH="$PWD/target/local-module-cache"
export SWIFTPM_MODULECACHE_OVERRIDE="$CLANG_MODULE_CACHE_PATH"
cargo test --locked --offline -p gentle-protocol -- --test-threads=1
cargo test --locked --offline -p gentle-render gel_workbench -- --test-threads=1
cargo test --locked --offline --lib gel_workbench -- --test-threads=1
cargo test --locked --offline --lib gel_image -- --test-threads=1
cargo test --locked --offline --lib glossary_cli_usage_smoke_commands_parse -- --test-threads=1
cargo test --locked --offline --lib cli_docs_mention_all_cli_glossary_paths -- --test-threads=1
cargo check -q --locked --offline
python3 -m unittest scripts.test_gel_tapestation scripts.test_tutorial_checkouts -q
```

The first macOS build attempts were blocked by the outer sandbox's cache and
nested `sandbox-exec` restrictions in the existing `apple-metal` Swift bridge.
Cache redirection plus the authorized unsandboxed build resolved that host issue;
no dependency/profile workaround was committed. The large test executable also
emits the linker's existing 16-MiB unwind-table-size warning; linking succeeds.

The deterministic runs found and repaired implicit color-conversion overflow,
discarded preview alpha, permissive unknown fields on the `inspect` action,
headless texture-delta disposal and incorrectly constructed synthetic xref
fixtures. Opaque previews retain their old bytes/warnings. Pixel conversion is
now explicit and versioned, with invalid float channels rejected rather than
clipped. PDF text tests follow Form XObjects instead of assuming lopdf's basic
text extractor traverses them.

The opt-in writer is documented in the
[fixture README](../test_files/fixtures/gel_images/tapestation_26r266/README.md#workbench-pdf-regression).
The inspected local report has SHA-256
`bb3faad1a27b44501c061af8821df7766767b9b34cada34b5feaa8888a3679ce`;
this is a receipt, not a cross-platform golden. It illustrates image-display
shapes only, not measured fluorescence or cleanup yield.

Native GUI input/responsiveness, Windows/Linux execution, real CLI process
replays, the full root/workspace suite and human scientific acceptance remain
open. No speed, release, raw-D5000-decoding or significance verdict is claimed.
Do not substitute the older size-only CLI evidence for this new path. Next
native review should inspect both cell lines independently, explicit object 5,
marker masks, source hashes and actual saved PDF/QQ output.

### Signed-Algebra Verification (2026-10-20)

This separate record uses local main
`5a8bd8cec472f68bc73b88448ed135b1131eb792` plus the signed-algebra working patch
on `gentle_rs_2_main`, Darwin 27 arm64, Rust/Cargo `1.100.0-beta.5` as above.
Default desktop/screenshot features are compiled, not exercised as capture.
Tests use one Cargo job, `CARGO_PROFILE_TEST_DEBUG=0` and the same local Swift
module-cache redirects shown above. No dependency or Cargo profile changed.
The desktop build also requires one missing `use eframe::egui` import in the
inherited planning-protocol UI; this prerequisite changes no planning behavior.

- `gentle-render` / `gel_workbench`: 20 passed, including signed density,
  gaps/zero, complete hashes, escaped provenance and rendering budgets.
- The normal protocol gate exposes an inherited canonical-parity failure:
  `catalog doctor` has a `ShellPassthrough` MCP override where the canonical
  contract requires `NotApplicable`. Its producer, glossary and parity inputs
  are unchanged from the rebase base. This is not repaired or counted as green.
  The remaining full protocol suite passes 90 tests with exactly
  `shared_workbench_route_keeps_canonical_parity_projection_current` excluded.
- Root library / `gel_workbench`: 61 passed, two opt-in visual writers ignored.
  This final run includes worker-side rasterization, stale/canceled/closed/project
  guards, work/dependency budgets, the narrow Shell inspection exemption and all
  six tutorial commands. Existing before/after report regressions also pass.
- Root library / `gel_image`: 42 passed, one failed. The introspection test
  `introspection_distinguishes_drafts_from_validated_reports` calls the same
  canonical parity builder and fails on `catalog doctor`. It is a cascade of
  the primary protocol failure above, not a second independent gel defect.
- Catalog and manifest freshness: one test passed in each filter, executed from
  the freshly built root test executable. The optional algebra PDF writer also
  passed explicitly, outside the default ignored-test set.
- `cargo check -q --locked --offline` and `cargo fmt --check`: passed on the
  final sources. Python fixture and LF/CRLF checkout suite: 80 tests, one documented
  existing-CLI-binary skip; the final focused fixture/tutorial rerun passes eight
  tests with the same skip. Counts overlap and must not be added.
- Documentation status, whitespace and roadmap's 400-line bound pass.
  Pre-commit session-close reports four OK, two warnings (intentional working
  changes/catalog notes and manual plan fidelity), no failures. The only
  out-of-gel source prerequisite is the named planning UI import above.

The exported synthetic report has two vector-PDF pages, both rendered with
Poppler and visually inspected: signed density, a printed zero label, unavailable
hatching, complete hashes and word-wrapped processing/warnings. Its SHA-256 is
`9fd0ba6c05423465c9031bfb1f6839e9a6eb0e14550a5a8b87eff4bfc2e30926`.
This is a local receipt, not a cross-platform golden or real-lane subtraction.
The report uses the explicit invented recipe from the optional 10.03 exercise.

Reproduce the checks with the common environment above. The known failing gates
are deliberately visible, not weakened:

```bash
cargo test --locked --offline -p gentle-protocol -- --test-threads=1
# Remaining protocol tests only; not a green canonical parity verdict:
cargo test --locked --offline -p gentle-protocol -- --test-threads=1 --skip shared_workbench_route_keeps_canonical_parity_projection_current
cargo test --locked --offline -p gentle-render gel_workbench -- --test-threads=1
cargo test --locked --offline --lib gel_workbench -- --test-threads=1
cargo test --locked --offline --lib gel_image -- --test-threads=1
cargo test --locked --offline --lib tutorial_catalog_check_passes_on_committed_tree -- --test-threads=1
cargo test --locked --offline --lib tutorial_manifest_check_passes_on_committed_tree -- --test-threads=1
GENTLE_GEL_ALGEBRA_TEST_REPORT=/absolute/new-synthetic-report.pdf cargo test --locked --offline --lib write_verified_algebra_pdf_for_visual_review -- --ignored --test-threads=1
cargo check -q --locked --offline
cargo fmt --check
python3 -m unittest scripts.test_gel_tapestation scripts.test_tutorial_checkouts -q
git diff --check
```

The freshness tests and opt-in writer were invoked directly from Cargo's built
test executable; their equivalent Cargo commands are shown above. Direct test
invocation must preserve `.cargo/config.toml`'s `RUST_MIN_STACK=16777216`.
The first standalone writer invocation omitted it and overflowed the default
test-thread stack before creating a report; rerunning with the repository's
existing setting passed. No stack configuration, dependency or production
profile was changed. Linking retains the existing 16-MiB unwind-table warning.

No full workspace, native GUI input/performance, Windows/Linux, live-agent,
release or scientific acceptance is implied by these scoped local checks.
