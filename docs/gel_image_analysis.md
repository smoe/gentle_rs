# Measured Gel Image Sizing

This is `.11` development on `codex/gel-image-dev`, not a `.10` release gate.
The first slice is manual calibration through the shared engine and Shell.
There is no click-based image editor or automatic band detector yet.

## What Is Measured

For linear double-stranded PCR/digest DNA, the result is fragment length in bp.
For denaturing SDS-PAGE or a corresponding Western blot, it is apparent
molecular mass in kDa. Neither establishes molecular identity. Intensity,
abundance, and concentration are not calculated by this workflow.

The method interpolates log10(size) between adjacent, confirmed ladder bands.
For example, halfway in migration between 1000 and 500 bp gives approximately
707 bp, not 750 bp. It does not assume that the whole gel has one globally
linear calibration. A reference pair must bracket an estimated position;
outside positions are retained with `outside_calibrated_range` and null size.

Select the correct ladder and confirm which visible band has which size. At
least two distinct reference bands are needed; three or more allow the operator
to review consistency. Missing, unresolved or merged ladder bands must not be
assigned by automatically shifting a list of expected sizes. Piecewise
interpolation passes through the reference points, so its zero residual at those
points is not independent evidence of accuracy.

SDS-PAGE size estimates conventionally use log molecular weight versus migration
or relative mobility ([Bio-Rad method](https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6210.pdf)).
Prestained protein ladders can have different apparent masses in different
gel/buffer systems ([manufacturer guidance](https://www.thermofisher.com/order/catalog/product/LC5625/faqs)).
For `prestained=true`, GENtle requires `sds_protein_kda` and an explicit
`gel_system`; use the corresponding vendor calibration values, not a generic
virtual-gel preset. `ladder.source` records that document or an explicit custom
or synthetic origin. The first slice accepts explicit values; catalog picking
will be added with the graphical editor.

## Import And Coordinates

PNG, JPEG and TIFF are accepted, up to 16 MiB compressed and 32 megapixels per
image. The self-contained project currently allows 128 MiB of original images.
Original bytes preserve 16-bit data and acquisition metadata. A PNG thumbnail
is generated once for display, never for coordinate measurement. JPEG produces
a lossy-input warning. Image decoding is bounded; damaged images fail cleanly.

TIFF requires `tiff_page: 0` explicitly. Later pages/exposures are not currently
supported or silently combined. Omit `tiff_page` for PNG and JPEG.

Coordinates always refer to zero-based pixel centers in the original decoded
image: x increases right, y down. Valid centers are `0..width-1` and
`0..height-1`, inclusive. No EXIF rotation, crop, mirror, lens correction or
blot registration is inferred. Migration can be `down`, `up`, `right`, or
`left`. Smaller ladder sizes must occur later along that direction.

This initial calibration assumes straight, locally comparable lanes. A distant
ladder, smiling gel, tilted photograph or separately cropped Western marker
exposure may violate that assumption. Do not treat the output as reliable until
the geometry has been reviewed. Separate marker exposures cannot currently be
registered through this workflow; use one correctly registered image.

## Manual Workflow

All commands below also work in GENtle's GUI Shell. With the CLI, use the same
`--state PROJECT.json` for each call so images and reports stay together.
`@FILE` may replace the inline request JSON.

```sh
gentle_cli --state gel-project.json gel-image import '{"image_id":"gel1","path":"gel.png"}'
```

Read `result.gel_image.sha256` from the response. Assign bands using original
image coordinates, not coordinates from a scaled screenshot. The following
request is a **synthetic coordinate example for a 160x220 image**, not a ladder
assignment for your photograph. Replace the digest, geometry and ladder values
with your own confirmed data and put the request in `analysis.json`:

```json
{
  "report_id": "measurement1",
  "image_id": "gel1",
  "image_sha256": "REPLACE_WITH_IMPORTED_SHA256",
  "size_kind": "linear_dna_bp",
  "migration": "down",
  "lanes": [
    {"id":"marker","label":"Ladder","min":{"x":0,"y":0},"max":{"x":40,"y":219}},
    {"id":"sample","label":"PCR","min":{"x":60,"y":0},"max":{"x":130,"y":219}}
  ],
  "ladder": {
    "lane_id":"marker",
    "label":"Synthetic reference bands",
    "source":"synthetic: known band sizes halving every 60 pixels",
    "bands":[
      {"id":"r1","center":{"x":20,"y":20},"size":2000},
      {"id":"r2","center":{"x":20,"y":80},"size":1000},
      {"id":"r3","center":{"x":20,"y":140},"size":500},
      {"id":"r4","center":{"x":20,"y":200},"size":250}
    ]
  },
  "sample_bands":[
    {"id":"b1","lane_id":"sample","label":"Unknown band","center":{"x":90,"y":110},"position_half_width_px":2}
  ]
}
```

```sh
gentle_cli --state gel-project.json gel-image analyze @analysis.json
gentle_cli --state gel-project.json gel-image inspect measurement1
gentle_cli --state gel-project.json gel-image export '{"report_id":"measurement1","path":"measurement1.json","format":"json"}'
gentle_cli --state gel-project.json gel-image export '{"report_id":"measurement1","path":"measurement1.tsv","format":"tsv"}'
gentle_cli --state gel-project.json gel-image export '{"report_id":"measurement1","path":"measurement1.svg","format":"svg"}'
```

The optional localization half-width gives min/max sizes by moving the marked
center along migration within that span. It is not a confidence interval and
does not include ladder uncertainty, gel distortion, transfer effects, or
protein-specific migration. Bounds crossing the calibration edge remain absent
with a warning. Unknown quality remains for human review; no automatic saturation
or smear detection is claimed.

Image/report ids cannot be overwritten. To revise an assignment, use a new
report id; engine undo/redo restores previous project state. Stored images use
shared immutable allocations so undo does not copy their pixel payloads.
Saving/reopening a project preserves the original image even if its input file
has moved. Portable reports include its basename and digest, not its absolute
workstation path. Import itself does not transmit data to an external agent.

Exports fail rather than overwrite existing files, and validate stored source
bytes and derived measurements first. SVG embeds the preview, lane boundaries,
ladder/sample legend, measured ladder scale, estimated-size table and warnings.
Dense ladder ruler labels are thinned to remain readable; JSON contains all
confirmed ladder assignments. JSON/TSV retain unrounded numerical estimates.

MCP uses the normal typed `op` tool with `ImportGelImage`, `AnalyzeGelImage`,
`InspectGelImageAnalysis` or `ExportGelImageAnalysis`; existing confirmation
guardrails still apply. JS/Lua/Python use their existing `op`/Shell wrappers,
not a second calibration implementation.

## Verification

`cargo test --lib gel_image` exercises synthetic, in-memory fixtures defined in
`src/gel_image/tests.rs`. They are explicitly drawn rectangles at independently
chosen coordinates, not outputs of GENtle's virtual migration algorithm.
Tests cover geometric interpolation, all four migration directions, protein
units, two-reference warnings, missing references, out-of-range rows, invalid
assignments, native 16-bit PNG/TIFF retention, source/report tampering,
project save/reopen, undo/redo, no-clobber exports, SVG semantics and shared
Shell execution. These are not real-gel accuracy or GUI responsiveness verdicts.

## Remaining Work

1. Dedicated image editor with lane/marker selection, zoom, calibration curve,
   ladder catalogs and clear original-to-display coordinate transforms.
2. Background decode/detection with progress/cancellation, cached previews,
   editable peak suggestions, quality checks and semantic GUI test identifiers.
3. Reviewed rotation/crop and two-sided ladder corrections; explicitly registered
   Western marker exposures; optional restricted-range log-linear fitting.
4. Reference-image validation with documented acquisition/provenance, sample
   association and expected-versus-observed comparison without identity claims.
5. Glen's image-size/memory/repaint audits and human wet-lab acceptance, plus
   richer image/report inventory and readiness discovery for the inner agent.

Densitometry, concentration estimates and automatic ladder identification remain
separate later capabilities, not hidden interpretations of this size report.
