# Measured-Gel Profiles And Lane Comparisons

Status: initial implementation added and locally Rust/PDF-verified on 2026-10-19;
native and scientific acceptance remain open. The [workbench contract](gel_workbench.md) records the
delivered subset and remaining native numerical-decoding/vendor-integration
limits. The numbered sequence below retains the scientific design rationale. The
owner requested a cleanup-distribution comparison and lane algebra on
2026-10-18. An optional read-only Claude review was offered; this plan has not
been Claude-reviewed; the owner authorized proceeding without another review
on 2026-10-19. Supplied Claude follow-up reviews were reconciled against
`bb1686bd` on 2026-10-20 as described below; that is separate from the original
unreviewed implementation plan. Do not make this an additional `.13` release blocker.

## Questions And Available Evidence

Two different questions must remain separate:

1. Does physical cleanup change the fragment-length distribution of the same
   library, including its short-fragment component?
2. Does combined custom/manufacturer priming show structure compatible with
   the custom-only preparation above the manufacturer-only smear?

The [public 26R266 fixture](../test_files/fixtures/gel_images/tapestation_26r266/README.md)
retains original PDFs, extracted gel rasters, annotation-assisted band positions
and separately transcribed vendor region summaries. These support size checks
and descriptive tables, not full numerical distributions. Before/after region
windows differ (400-5000 versus 500-7000 bp), contrast differs, and expiry warnings
are present. The 7000 bp bound is beyond the nominal 100-5000 bp D5000 range.
Never label a region mean as a peak, reconstruct a smear from a sparse peak
table, or divide unmatched-window percentages to report cleanup enrichment.

The owner also supplied the original before-cleanup `.D5000` file. On
2026-10-18 its SHA-256 was
`fd2150e15daf0fb32ad11028bc8ad4355ecaf4b64237ef980aaf8a1de66ac51d`
(1,367,616 bytes). ZIP-directory inspection found encrypted `ModelData.xml` and
`RunData.xml`; only `Thumbnail.png` and `Logs.zip` are unencrypted. No numerical
curves were decoded. The new native-source importer can retain an explicitly
opted-in thumbnail, but cannot decode those numerical members. Do not brute
force, treat the thumbnail as raw detector data, or publish embedded operational
logs as a scientific fixture. The second file/export is pending; development
need not wait for it.

Ask for vendor-exported numerical electropherograms of both measurements,
including the actual x-axis and signal units, calibration/normalization,
marker handling, sample associations, dilution and loading metadata. CSV/XML
export is supported by the [manufacturer's software](https://www.agilent.com/en/p/electrophoresis/software-laptops/tapestation).
Verify the supplied export actually contains curve samples, not only tables.
Do not prescribe a menu or assume a trace schema before inspecting the export.

## Pairing Contract

Source/run plus lane is the identity; lane names alone repeat across runs.
Use the fixture's `region_summaries.json` as the single recorded pairing source.

| Question | Cell 1 | Cell 2 |
| --- | --- | --- |
| Cleanup, custom only | before B1 vs after C1 | before C1 vs after D1 |
| Cleanup, manufacturer only | before D1 vs after E1 | before E1 vs after F1 |
| Combined vs components, after cleanup | custom C1, manufacturer E1, combined G1 | custom D1, manufacturer F1, combined H1 |

L1/L2 are ladders, not samples. The manufacturer pair E1/F1 straddles L2; neither
adjacency nor a shared well name establishes a pair. No combined-before sample
exists in the supplied report. Pairing does not prove identical preparation,
loading or recovery; preserve unknowns until the owner supplies those facts.

## Minimal Implementation Sequence

### P1: Retain And Overlay Numeric Profiles

After inspecting a real export, introduce one shared engine-owned portable
profile/report contract under the existing measured-gel boundary, not a GUI
calculator or a permanent TapeStation-only Python analysis path. Keep band
sizing and dense trace evidence separate; do not change `GelImageAnalysisReport`
v1 into an abundance report. Expose new operations through the existing typed
`op` and shared Shell routes before adding GUI affordances or agent guidance.

Bind original bytes/hash, run/sample/lane, assay, extraction/import version,
coordinates, axis units, signal kind/units, calibration, marker masks,
normalization and warnings. Distinguish raw detector traces, calibrated vendor
traces, digitized electropherograms and image-derived display intensity. Missing
values and unsupported calibration remain typed unavailable, never zero.
Keep native archives local unless their contents are reviewed for publication.

First deliverable: per-cell overlays with explicit common x-axis and independent
source provenance. A picture-derived brightness curve may be useful later as a
clearly labelled display diagnostic, but must never be silently interchangeable
with calibrated fluorescence or concentration. Native `.D5000` import and PDF
curve digitization are not prerequisites or authorized format guesses.

### P2: Distribution Summaries And Cleanup

Use the explicit intersection of valid calibrated coverage, with no
extrapolation, interpolation across missing intervals or markers. Do not crop
away short fragments merely because the existing reports' region tables begin
at 400/500 bp: their loss is the very question. For these D5000 data the candidate
comparison envelope is 100-5000 bp, subject to actual exported coverage.

Offer unnormalized traces and an explicit area-normalized shape comparison;
disclose the normalization interval and factors. Only report cumulative
fractions, quantiles and short/long-region areas when the integration measure
is justified. Detector intensity sampled against migration time does not become
a mass density per bp just by relabelling x: retain the calibration/Jacobian or
vendor area representation and use area-conserving binning. An arbitrary FU-vs-bp
integral must not be presented as ng, molecule count or mass fraction.

Keep fluorescence-, mass- and molar-weighted distributions distinct. A 500 bp
analysis threshold may be chosen explicitly; it is not inferred to be the
laboratory's cleanup cutoff. Total recovery additionally needs dilution and
input/output volumes and appropriate quantitative calibration. Unknown values
must prevent a yield claim, not prevent an honest shape-only overlay.

### P3: Explicit Lane Algebra

Only after axis, coverage, signal kind and normalization admission succeeds,
allow a declared linear combination over compatible profiles. Start with
explicit coefficients, provenance and signed output; never clip negative
residuals to zero. Sum/difference require compatible units; division additionally
requires a stated denominator floor, with zero/missing denominators unavailable.
Do not silently average runs or choose scale factors to maximize an apparent
custom contribution.

Useful exploratory expressions on the after-cleanup profiles are:

```text
R1 = G1 - alpha1 * E1     # cell 1: combined minus scaled manufacturer
R2 = H1 - alpha2 * F1     # cell 2: combined minus scaled manufacturer
G1 approximately alpha1 * E1 + beta1 * C1
H1 approximately alpha2 * F1 + beta2 * D1
```

These are mathematical shorthand, not GENtle command syntax or a validated
additive chemistry model. Choose coefficients by an explicit reviewed loading
or normalization policy. A later fitted decomposition must disclose its fit
region, constraints, uncertainty and lack of identifiability, and must not be
the default first implementation. Reverse transcription/amplification with two
primer sets need not equal a physical mixture of the two separate preparations.

Compare residual locations with custom-only structure, including registration
uncertainty; a residual is not gene-specific enrichment. The annotation-assisted
first replay has custom short marks near 697/700 bp and combined marks near
644/635 bp, so those must not automatically be declared identical bands.
Assignment to PATZ1, E2F1 or the TP53 family needs independent sequence evidence.
Biological replicates and controls are needed before inferential claims.

## Tests And Acceptance

### Signed-Algebra Follow-Up (2026-10-20)

Claude's supplied read-only review of `3a4e458b` was checked against unchanged
gel sources on main. The revised bounded slice adds replay-verified standalone
residual plots and report selection, not another statistical test. Freeze v1
replay and expose full operand provenance before plotting; divide integrals by
bin width for explicitly labelled display density. Cache the shared result in
the GUI worker path. Fail closed on mismatching/unverifiable report export while
keeping legacy inspection. Leave existing comparison pages and reference
overlays unchanged; use `1/-1` only in a synthetic lesson, not as an approved
public biological subtraction. Scope/verification are recorded in
[the workbench contract](gel_workbench.md#verify-and-display-a-saved-residual).

### Comparison Display Follow-Up (2026-10-20)

The older algebra review is already implemented; do not reopen frozen algebra
v1 or add deferred reference/fitted overlays. The new ordinary comparison view
uses the replay-verified `GelComparisonReport`, not a parallel engine-prepared
result. Density is a pure renderer projection alongside CDF geometry. Extract
the existing full-report verification for export and a GUI worker; raw plain
inspection remains available. This existing equality includes warning prose:
algorithm/warning changes need a future compatibility policy, not a hidden fix.

Add optional PDF `comparison_display`, omitted/explicit legacy integrals or
explicit signal per size width. Keep old request/report schemas and omitted-mode
SVG bytes unchanged. Reject a mode for JSON or no selected comparison. Use one
common numeric density axis with acquisition/calibration disclosures; different
source hashes do not prove different exposure settings. Preserve zeros/gaps and
refuse unsafe numeric geometry. No normalization/fitting or new statistics.

The GUI defaults to density per session, verifies/rasterizes only the shared
profile panel in the worker, and retains its native summary/QQ. Cache by revision,
comparison and mode, under the existing project-instance binding. Discard stale,
canceled, closed and replaced-project results. Mode changes invalidate immediately;
prepare explicitly when the worker is free, rather than dropping a mode request.
Use the same session selector for PDF export and no persisted display field.

Extend the **hand-written** 10.03 Markdown with its own density markers, leaving
the six algebra commands alone. Source metadata declares `hand_written_markdown`;
it is not a generated chapter. Correct the supplied example: A/C's union is
`[100,125,200,400]`, weights `[25,75,200]`, not four bins. Independent tests check
flat density, preserved area, legacy bytes, strict requests, failure/no-write,
read-only state, stale GUI guards and vector PDF output. Native/platform/real-data
acceptance remains open; no new Glen request, release gate or performance claim.

### Existing Acceptance Boundary

- Preserve the existing size-only real-gel regression without changing its
  interpretation or passing expected sizes into engine requests.
- For numeric import: bounded/malformed data, unknown units, unsorted/duplicate
  axes, mismatched runs, hashes, markers, missing intervals and real LF/CRLF
  exports; native encrypted input must fail as unsupported without guessed data.
- Independently calculable synthetic profiles: identical subtraction is zero
  on shared valid support, missing remains missing, signed residuals survive,
  zero denominators stay unavailable, area-conserving rebinning preserves its
  declared measure, and scaling/normalization are recorded rather than hidden.
- Test engine persistence and GUI/CLI/MCP/agent reachability over the same
  request/report. UI selection, zoom and export never recompute the analysis.
- Extend tutorial 10.03 only as each capability exists. Retain original vendor
  data and independent expected values; compare both cell lines separately.
  Glen's native checks and the owner's scientific interpretation remain distinct
  from import/math tests, and no result here grants release acceptance.

The next implementation decision is the smallest importer/profile contract
supported by the actual numerical export. Current work may finish the tutorial,
pairing and evidence boundaries without inventing an unavailable distribution.
