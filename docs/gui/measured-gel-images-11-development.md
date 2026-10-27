### Measured Gel Images (`.11` Development)

<!-- manual-counterparts:start -->
Related CLI reference: [Measured Gel Images](../cli/gel-images.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Open **Patterns -> Gel Image Analysis...**, the command palette, or
`ui open gel-image-editor`. Import a PNG/JPEG or explicitly selected TIFF page 0,
draw lane rectangles, enter known ladder sizes and click band centers. Select/move
marks or correct their coordinates in the assignment editor. Zoom and scroll
only change the view. Drafts are undoable project data and survive Save Project.
**Analyze and save report** validates the assignments through the same engine as
the CLI; changed assignments mark the old result as outdated and disable export
until reanalysis. Results include a calibration plot and JSON/TSV/SVG exports.
See [the walkthrough and limitations](../gel_image_analysis.md#graphical-walkthrough).
Image loading and analysis run in background workers; canceled or stale mutation
results do not commit. Existing virtual gel previews remain simulations.
The [development workbench](../gel_workbench.md) additionally exposes PDF/native
preview import, reviewable band proposals, stored profiles from selected sample
lanes, cross-gel selectors, QQ/shape summaries, explicit signed combinations
and saved before/after PDF reports. These use the shared engine in background
workers; [test scope and evidence](../gel_workbench.md#verification-boundary) stay
separate from pending native responsiveness and real-gel acceptance. Numerical
`.D5000` decoding is unsupported; a native thumbnail is
display-only. This development is not an added `.13` release requirement.

For saved signed results, **Verify and display saved algebra** prepares a
shared SVG and its raster preview in the background and caches it by revision/result ID; the GUI only installs the texture. The graph
retains negative values, a zero line and missing intervals with units per
bp/kDa; complete operand processing/calibration and warnings remain accessible.
Thin styled curves show the recorded coefficient times each operand's rebinned
signal per size unit, on the same axis as the thick teal residual. The numbered
legend follows request order, including duplicate IDs; zero terms have no curve.
Each operand breaks at its own gaps, and a visible operand never fills hatched
residual coverage. No coefficients are chosen by the display.
**Optional custom-only reference** adds a caller-selected, unscaled panel on
the same size range with its own labelled Y axis, units, native bins and gaps.
No role/normalization or amplitude comparability is inferred. Changing/clearing
this selection invalidates cached and late worker pictures. Export with this
reference requires exactly the inspected algebra ID; clear it for unpaired or
multiple-result export. Full reference provenance remains separate from operands.
Ordered algebra export appends verified result/provenance pages. Legacy or
mismatching results stay diagnosable but cannot be report-exported. Existing
comparison pages retain their legacy bin-total display when export mode is omitted.
**Verify and display saved comparison** prepares a verified profile-only panel
in the worker. Selecting a saved comparison alone shows neither its summary
nor QQ: both appear with the panel only after successful verification. A failed
check shows the reason, not trusted numbers or a previous plot. For diagnosis,
CLI `inspect_comparison` still exposes the raw, unverified saved report; it does
not authorize display or export. The GUI's session-only selector defaults to
**Signal per size width** and also controls
PDF comparison export; **Legacy bin integrals** is available explicitly. Mode,
selection or revision changes invalidate the picture; canceled/closed/project
changes reject late results. If busy, wait and prepare the selected mode again.
Equal units/common axes do not prove acquisition comparability. No automatic
scale, inferred reference pairing, statistical significance or enrichment is inferred.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: microRNA Target Scan](microrna-target-scan.md) | [Next: Pool Distribution (Engine Ops)](pool-distribution-engine-ops.md)
<!-- manual-page-nav:end -->
