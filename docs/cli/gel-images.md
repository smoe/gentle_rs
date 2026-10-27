# Measured Gel Images

<!-- manual-counterparts:start -->
Related GUI reference: [Measured Gel Images (`.11` Development)](../gui/measured-gel-images-11-development.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## Measured Gel Images (`.11` Development)

The following routes share typed engine operations with GUI Shell and MCP `op`:

- `gel-image import REQUEST_JSON_OR_@FILE`
- `gel-image save-draft REQUEST_JSON_OR_@FILE`
- `gel-image analyze REQUEST_JSON_OR_@FILE`
- `gel-image inspect REPORT_ID`
- `gel-image export REQUEST_JSON_OR_@FILE`
- `gel-image workbench REQUEST_JSON_OR_@FILE`

The original five routes size manually confirmed bands in bp or apparent kDa, not band
intensity or concentration. See [measured-gel sizing](../gel_image_analysis.md)
for requests, original-image retention, calibration limits and exports.
`save-draft` retains incomplete `GelImageAnalysisRequest` assignments without
creating a measurement. The GUI editor uses the same operation. Its shared
window intent is `ui open|focus|close gel-image-editor`; headless hosts return
the intent rather than opening a native window.

Introspection describes these routes and their confirmation requirements.
Only validated analyses contribute `report.exists` (`gel_image_analysis`);
imported images and saved drafts are not sizing reports. Report readiness is
advisory: execution still validates the original image hash and calibration.

The [development workbench](../gel_workbench.md) adds typed `action` requests for
PDF extraction/explicit native thumbnails, read-only band proposals, persisted
image/numeric profiles, QQ/KS-distance/Wasserstein comparisons, signed linear
algebra and saved before/after PDF/JSON export. It does not decode encrypted
`.D5000` detector data, infer concentration from images or supply pixel-based
significance p-values. Every command also works in GUI Shell and MCP `op`.

`inspect_algebra` returns read-only replay verification and the shared signed
density presentation; `export` accepts optional ordered `algebra_ids`. Mismatching
or unverifiable residuals cannot be report-exported. See
[the contract and limits](../gel_workbench.md#verify-and-display-a-saved-residual).

For an unscaled, separately labelled custom-only reference panel, pass optional
`reference_profile_id` to `inspect_algebra`. PDF/JSON export uses
`algebra_reference_profile_id` with exactly one selected `algebra_ids` entry.
The caller selects the role; no priming identity, subtraction coefficient or
amplitude match is inferred. Native reference bins retain their own units,
calibration and gaps; unsafe or stale sources refuse. See the
[reference contract](../gel_workbench.md#unscaled-custom-only-reference-panel).

PDF `export` accepts optional `comparison_display`: `signal_density_per_size`
plots saved common-bin weight / size width; `bin_integral` or omission retains
the byte-identical legacy comparison SVG. This field requires at least one
`comparison_ids` entry and is refused for JSON. It changes no stored integrals,
statistics or confirmation requirements. See [density display](../gel_workbench.md#comparison-density-display).


TATA-box evidence is available through
`promoters tata-screen REQUEST_JSON_OR_@FILE [--output REPORT.json]` and the
separate, digest-gated `promoters tata-materialize REQUEST_JSON_OR_@FILE`.
See [TATA-box evidence](../tata_box_evidence.md) for the three evidence sources,
EPD input files and exact request examples. These routes also work through
the shared Shell, MCP `op`, and normal operation/workflow adapters.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: MicroRNA Target Scans](microrna.md) | [Next: Racks, Pools and Ladders](racks-pools.md)
<!-- manual-page-nav:end -->
