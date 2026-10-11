# Estimate a PCR Band's Size from a Gel Image

> Which size is consistent with this band's migration between two confirmed
> ladder bands, and where must I refuse an estimate?

This 15-minute exercise uses a drawn image, not an experiment. It teaches the
manual calibration path shared by the GUI, CLI and agents. Automatic band
detection, identity, concentration and experimental validation are not included.

## Before You Start

Run from the repository root. Keep a fresh writable output directory and one
disposable project, for example `/tmp/gentle-gel-lesson/state.json` (create the
parent directory first). Use the same `--state` on every terminal invocation.
A separate CLI does not inherit the GUI's unsaved project.

Inputs and recreation are documented in the
[fixture README](inputs/gel_sizing/README.md). No download or external image
library is needed:

- [Drawn image](inputs/gel_sizing/gel_sizing_synthetic.png).
- [DNA assignments](inputs/gel_sizing/dna-analysis.json).
- [Exact checksums](inputs/gel_sizing/SHA256SUMS).

![Synthetic ladder and two sample bands](inputs/gel_sizing/gel_sizing_synthetic.png)

The picture is deliberately small. Zoom for marking; zoom must not change its
original zero-based pixel coordinates. The left lane is the ladder and the
right lane contains a bracketed band and an unbracketed band.

## 1. Import the Original Image

**GUI:** choose **Patterns -> Gel Image Analysis...**, then **Import image...**.
Open the PNG above. DNA selection is not required. Do not check TIFF page 0 for
a PNG. The shared `ui open gel-image-editor` opens the same independent editor
when applied by the GUI host; a headless receipt is not an opened window.

**CLI / shared engine:**

```sh
gentle_cli --state /tmp/gentle-gel-lesson/state.json gel-image import '{"image_id":"gel_teaching","path":"docs/tutorial/inputs/gel_sizing/gel_sizing_synthetic.png"}'
```

Inspect `result.gel_image.sha256` and compare it with the image digest in
`dna-analysis.json` and `SHA256SUMS`. An image mismatch means stop, not edit the
expected digest. The project retains original bytes; a thumbnail is display only.

**Ask the inner agent:** "Open the measured-gel editor. Import the retained
synthetic PNG only after showing the exact command. Keep its original pixel
coordinates and verify its digest; do not infer bands or send the image away."

## 2. Confirm Lanes and Ladder Assignments

**GUI:** draw the ladder rectangle first, then the sample rectangle. Select bp
and downward migration. Use numeric correction in **Assignments** to match:

| Item | Original pixels | Assignment |
| --- | --- | --- |
| Ladder lane | x=0..40, y=0..219 | synthetic marker |
| Sample lane | x=60..130, y=0..219 | synthetic sample |
| Ladder centers | (20,20), (20,80), (20,140), (20,200) | 2000, 1000, 500, 250 bp |
| First sample | (90,110) | unknown; localization half-width 2 px |
| Second sample | (90,210) | unknown; outside reference range |

Enter the invented ladder's source explicitly. In a real gel, use the supplier,
lot and confirmed visible reference bands; do not shift missing bands down a
predefined list. Marks are user assignments, not automatically detected peaks.

**CLI:** the committed request supplies exactly these assignments:

```sh
gentle_cli --state /tmp/gentle-gel-lesson/state.json gel-image analyze @docs/tutorial/inputs/gel_sizing/dna-analysis.json
gentle_cli --state /tmp/gentle-gel-lesson/state.json gel-image inspect synthetic_dna_1
```

**GUI:** choose **Analyze and save report** after confirming all assignments.
The GUI may assign its own image/report IDs; do not paste the fixture's IDs over
an unrelated imported image. Ask an agent to draft the request using the actual
imported ID/digest and confirmed marks, then explicitly approve the analysis.

## 3. Explain the Answer

The first band is halfway in migration between 1000 and 500 bp:

```text
log10(size) = (log10(1000) + log10(500)) / 2
size = sqrt(1000 * 500) = 707.106781... bp
```

It is **not 750 bp**. GENtle interpolates log size between adjacent references,
not arithmetic size and not a universal global gel curve. The second band lies
beyond y=200: expect `outside_calibrated_range` and a null estimate. Do not
extrapolate or treat null as zero bp.

The optional +/-2 px bounds reflect mark localization only. They exclude ladder
uncertainty, gel smiling, sample-dependent migration and photography distortion;
they are not a statistical confidence interval. Straight comparable lanes and
correct calibration remain operator responsibilities.

## 4. Export and Reopen Without Replacing Evidence

```sh
gentle_cli --state /tmp/gentle-gel-lesson/state.json gel-image export '{"report_id":"synthetic_dna_1","path":"/tmp/gentle-gel-lesson/sizing.tsv","format":"tsv"}'
gentle_cli --state /tmp/gentle-gel-lesson/state.json gel-image export '{"report_id":"synthetic_dna_1","path":"/tmp/gentle-gel-lesson/sizing.svg","format":"svg"}'
```

**GUI:** select SVG/TSV/JSON and **Export...**, then **File -> Save Project**.
Reopen the project and the editor: original image and assignments should remain
available even if the external source is moved. Editing a mark makes the old
report outdated; analyze under a new report ID. Existing exports are refused,
not overwritten. Use fresh output names for the second run.

Ask an agent to inspect/export the saved report, not to claim a screenshot or
an experimental result from a successful headless command. A defensible summary
is "approximately 707 bp under these confirmed synthetic ladder assignments;
the lower band is outside calibration", not "the PCR product is identified".

## Appendix: Apparent SDS Protein Mass

Use the separate [protein request](inputs/gel_sizing/protein-analysis.json) on
the same imported `gel_teaching` image, or confirm its assignments in the GUI.
It sets `sds_protein_kda`, synthetic 100/50/25/12.5 kDa markers and an explicitly
synthetic prestained gel-system label. The midpoint is sqrt(50*25), about
35.355339 **apparent kDa**. The outside band is still unestimated.

```sh
gentle_cli --state /tmp/gentle-gel-lesson/state.json gel-image analyze @docs/tutorial/inputs/gel_sizing/protein-analysis.json
gentle_cli --state /tmp/gentle-gel-lesson/state.json gel-image inspect synthetic_sds_1
```

Do not use these invented masses for a real Western blot. Prestained marker
calibrations depend on gel/buffer system; record the applicable vendor values.
Switching units does not convert existing ladder assignments. Apparent mass is
not protein identity or quantity, and separate marker exposures are not silently
registered. See the [full method and limitations](../gel_image_analysis.md).

## Verification and Acceptance

Fixture recreation, source checks and independent arithmetic are covered by
`scripts.test_tutorial_gap_lessons`; binary/LF preservation by
`scripts.test_tutorial_checkouts`. Existing `gel_image` Rust tests exercise the
shared algorithm and editor. Run those and the exact commands above at one
source SHA before calling this an executed learner pass. Native screenshots,
Windows/macOS GUI and real-image scientific acceptance remain pending with
Glen/the owner; this manual guide is not an automated GUI contract.
