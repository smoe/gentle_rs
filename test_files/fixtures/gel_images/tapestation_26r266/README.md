# Public 26R266 TapeStation gel-sizing regression

## Origin and permission

Real experimental reports supplied by GENtle's owner, who explicitly approved
public use as regular GENtle tests on 2026-10-18. These are not synthetic gels.
The original experimental identifiers are retained; cell lines remain the two
anonymous sample identifiers in the reports. No sequence-level gene enrichment,
primer identity or full-length cDNA claim can be inferred from their sizes.
Agilent/TapeStation/ScreenTape product names and report artwork identify the
source, not GENtle manufacture or endorsement.

| Retained original | Original supplied filename | Pages |
| --- | --- | --- |
| `before.pdf` | `LibCheck_26R266 - 2026-10-06 - 12.46.57.pdf` | 7 |
| `after.pdf` | `26R266 ONT Libs nach extra CU - 2026-10-13.pdf` | 9 |

`SHA256SUMS` binds original PDFs, decoded PNGs and the raw coordinate/expectation
manifest. Never normalize those bytes or regenerate to conceal a mismatch.
Scoped Git attributes plus the LF/CRLF checkout regression protect these inputs.
`region_summaries.json` separately transcribes all ten vendor Region Table rows,
including page references, original measurement fields, source PDF digests and
explicit cleanup/priming pairings. It is not a full distribution or a GENtle
measurement. Its checksum is retained alongside the original fixture hashes.
Recreate it by reading the Region Tables on the listed pages of these PDFs;
`pdftotext -layout before.pdf -` and `pdftotext -layout after.pdf -` provide a
text cross-check. Preserve field labels and units and review against the rendered
pages; do not substitute the separate Sample Table or Peak Table concentrations.

## Reproduction

The gel on page 1 is a 1375 x 731 RGB raster in PDF image object `5 0`, not a
screenshot. Extracted with Poppler `pdfimages` 26.06.0 without resizing, cropping,
contrast adjustment, recompression to JPEG or redrawing:

```sh
pdfimages -f 1 -l 1 -png before.pdf before-extracted
pdfimages -f 1 -l 1 -png after.pdf after-extracted
```

The `*-000.png` outputs are `before.png`/`after.png`; later extracted objects
are the separate soft mask and warning icons, not additional gel lanes. The
retained RGB rasters were visually checked against the rendered PDF pages on
white backgrounds. Both PDFs and PNGs are retained; these are vendor-processed
gel displays, not original detector/acquisition data. Another extractor version
may encode PNG bytes differently: review any difference, do not update hashes
automatically.

`manifest.json` records all input dimensions, page/table references, source
warnings, contrast, sample groups and coordinates. Pixel positions were read
from TapeStation's small peak-indicator glyphs next to the lanes, then reviewed
against the gel and electropherograms. Thus this is **annotation-assisted manual
sizing**, not an independent peak detector or independent assay validation.
Half-pixel centers describe the raster glyph extent; they do not establish
half-pixel experimental accuracy. The 1-pixel localization span in the request
is not a confidence interval or the total measurement uncertainty.

GENtle's existing `piecewise_log10_size_v1` engine alone calculates sizes.
Vendor sample sizes live separately from the engine request and never select
coordinates, calibration coefficients or expected return values. The 5% relative
comparison tolerance was declared before the first GENtle run; it is a software
consistency criterion, not a D5000 accuracy specification.

## Cases and limitations

- `before_A1`: 13 vendor sample-peak values from pages 3/4/6/7, using the A1
  ladder values from page 2. Upper/lower marker positions are explicit negative
  controls: they return null outside the 100-5000 bp calibration, not sample sizes.
- `after_L1`: C1/D1/E1; `after_L2`: F1/G1/H1. These are explicit test groupings
  by display block, not reconstructed instrument-internal calibration settings.
  The manufacturer-primer cell-line pair E1/F1 straddles the second ladder;
  adjacent displayed lanes must not be blindly grouped into pairs.
- Each after-case checks seven interior bands of the *other* ladder against
  their nominal sizes, for 14 cross-ladder registration checks. Those are
  controls, not biological samples. Do not count fitted calibration points as
  independent agreement observations.
- The after-PDF contains no sample Peak Tables. Its 15 in-range selected peak
  positions receive GENtle estimates, but **no vendor numerical peak comparison**
  is claimed. The additional D1 mark above the 5000 bp reference returns null.
  Region means are not substituted for absent peak sizes.
- ScreenTape expiry warnings apply to every before-lane and after L1/C1/D1/E1.
  They remain in the original PDF and manifest, and in the selected ladder's
  source note. Absence of a warning for L2/F1/G1/H1 is not proof of calibration.
- Display contrast is 130% before and 119% after. Region windows change from
  400-5000 to 500-7000 bp. The 7000 bp region boundary exceeds the D5000 nominal
  100-5000 bp sizing range ([manufacturer guide](https://www.agilent.com/cs/library/usermanuals/public/D5000_QuickGuide.pdf)).
  A 10000 bp alignment marker does not extend that range. No fluorescence,
  concentration, region mean, cleanup yield or sequencing read length is
  calculated or validated here.

## Execution and use in GENtle

The optional real-data companion in
[tutorial 10.03](../../../../docs/tutorial/10-03_measured_gel_sizing_gui_cli.md#real-tapestation-companion)
teaches replay, GUI report inspection, sample pairing and the limits of the
available region summaries. Targeted priming and physical length enrichment
are separate dimensions. The
[profile-comparison plan](../../../../docs/gel_lane_profile_comparison_plan.md)
keeps full distributions and lane algebra pending numeric trace admission;
neither is implemented by this fixture or its Python sizing replay.

- `src/gel_image/tests.rs::tapestation_real_images_preserve_sources_and_match_vendor_peak_sizes`
  exercises real PNG import, engine analysis, persisted-project reopen, shared
  Shell inspection and SVG export. It is an ordinary root-library test.
- `scripts.test_gel_tapestation.TapeStationFixtureTests` checks source identity,
  input/result separation, missing-data semantics and comparison failures without
  invoking a binary. `scripts.test_tutorial_checkouts` additionally exercises real
  LF/CRLF checkouts, with a missing-JSON-policy negative control.
- `scripts.test_gel_tapestation.TapeStationBinaryTests` runs the real CLI only
  when `GENTLE_GEL_TEST_CLI` is explicitly supplied. Native CI supplies the CLI
  already built by its ordinary smoke; no extra Cargo build or network is added.
  Without a supplied binary that optional local test explicitly skips, not passes.

To retain a fresh replay, use an existing binary and a **new** output directory:

```sh
python3 scripts/gel_tapestation_acceptance.py \
  --gentle-cli /absolute/path/to/gentle_cli \
  --output-dir /absolute/path/to/new-audit-directory
```

The runner isolates HOME/config/cache/temp, performs real direct-CLI operations,
then reloads each report through the shared Shell in another process. It exports
JSON, TSV and SVG, and retains the project, requests, comparisons and a hashed
receipt. Existing output directories are refused. It checks the binary digest
before/after, requires its full source identity and verifies that the four gel
implementation files match that source. It never compiles, runs an inner agent,
downloads a resource, submits data or modifies an existing project.

## Observed first replay, 2026-10-18

macOS arm64; supplied CLI source `1c5be3175323eec41359b980b5a41e22c7edaec1`,
binary SHA-256 `f38780893ca3ad2087ff4ee15224d4e70b17217ec9dd3cfec53ae94f89826561`.
The gel implementation files were unchanged at the working source `78e43576`.
The binary does not embed effective profile settings; no speed/profile or new
whole-candidate acceptance claim is made.

All 13 vendor sample-peak comparisons passed: maximum relative difference
0.890%. All 14 cross-ladder checks passed: maximum 3.671%. All three deliberately
out-of-calibration marks remained null. The other 15 after-estimates correctly
had no vendor peak value to compare. No engine sizing defect was demonstrated;
the missing real-image, grouped-ladder and truthful-comparison regression
coverage is added without changing calibration mathematics.

| Before lane | TapeStation bp | GENtle bp, rounded for this summary |
| --- | --- | --- |
| B1 | 215, 280, 683 | 215, 280, 680 |
| C1 | 203, 278, 475, 673, 1013, 1219, 2005, 2435 | 204, 279, 473, 673, 1010, 1218, 2007, 2434 |
| D1 | 2047 | 2053 |
| E1 | 1149 | 1139 |

These observations are not committed runtime goldens and are not expected output
inputs. The tests calculate from the retained geometry with the current engine;
the runner retains unrounded results and actual binary identity for each run.
