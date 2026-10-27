# Multi-Gene Publication Bundles

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

## Multi-gene publication bundle

`gentle_publication_report REQUEST.json OUTPUT_DIR [--pdf]` resolves one
portable gene-set manifest into a responsive HTML report, complete primer
tables, copied assets, and an optional printable companion. HTML and print
derive from the same `resolved-report.json`; see
[`gene_set_publication_reports.md`](../gene_set_publication_reports.md).
`bundle-manifest.json` records the request digest and the digest/media type of
every finalized bundle artifact. `generation-report.json` records that
manifest's digest while retaining the legacy `copied_files[]` field.

For content-bound isoform-assay dossiers, the same binary accepts
`gentle.gene_isoform_assay_publication_request.v1` plus optional
`--profile ID`, `--blocks ID,ID`, and `--pdf`. The shared-shell equivalent is:

```text
primers publish-gene-isoform-study REQUEST.json OUTPUT_DIR [--profile ID] [--blocks ID,ID] [--pdf]
```

This form writes one HTML page per gene and uses `print.html` as the sole PDF
source. GENtle searches `PATH` and conventional platform locations for
Chromium, Google Chrome, or Microsoft Edge on Linux, macOS, and Windows before
falling back to a `chromium` command. Set `GENTLE_BROWSER_BIN` to select another
Chromium-compatible executable; an explicit override always takes precedence.

This page documents command-line entry points for GENtle.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Tracks and Arrays](tracks-arrays.md) | [Next: RNA Read Panel](rna-read-panel.md)
<!-- manual-page-nav:end -->
