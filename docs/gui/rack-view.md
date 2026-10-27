## Rack View

<!-- manual-counterparts:start -->
Related CLI reference: [Racks, Pools and Ladders](../cli/racks-pools.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Rack view is the first physical-placement layer on top of arrangements.

- `Open Rack` opens one saved rack/plate draft linked to one or more
  arrangements.
- The rack header shows:
  - rack name / rack id
  - current profile
  - arrangement chips for each contiguous arrangement block on the rack
- The rack window also includes a compact in-place help strip so the sample
  move, block move, preview, and keyboard-cancel modes stay visible while you
  work.
  - the strip now also shows tiny arrangement-color chips and can be collapsed
    with `Hide help` / `Show help`
  - shortcut hints are rendered as small keycap-style labels for `Cmd`,
    `Ctrl`, and `Esc`
  - the collapsed/expanded state is remembered in project metadata, so a saved
    project reopens with the same rack-help visibility you last used
  - when the strip is not pinned open, GENtle auto-minimizes it after a few
    successful rack moves; `Pin open` keeps it visible
- At compact window sizes, rack/profile/template controls wrap instead of
  disappearing beyond the right edge. The four help cards become one concise
  interaction summary, and the rack grid receives its own bounded scroll
  viewport so the label and physical-carrier export controls remain visible
  below it; scrolling a large plate no longer hides the export block.
- Built-in profiles:
  - `Small tube rack (4 x 6)`
  - `Plate 96`
  - `Plate 384`
- `Custom` profile authoring is now available directly in the rack window:
  - pick `Custom`
  - edit `rows` / `columns`
  - click `Apply Custom`
  - occupied positions are reflowed under the same order-preserving rules
- Rack snapshots now also carry:
  - `Template`
    - `Bench rows`
    - `Plate columns`
    - `Plate edge avoidance`
    - `Apply Template`
  - `Fill`
    - `Row-major`
    - `Column-major`
  - `Blocked`
    - comma- or whitespace-separated reserved coordinates
    - `Apply Blocked` / `Clear Blocked`
- Blocked slots render directly in the rack grid and are excluded from
  placement/drop targets.
- Coordinates are A1-style everywhere, and row labels continue beyond `Z` as
  `AA`, `AB`, ...

Interaction model:

- Drag one occupied sample position onto another slot to move that sample
  within its arrangement block.
- `Command`/`Ctrl`-click occupied sample positions to build a same-arrangement
  multi-selection, then drag one selected sample or click a target slot to move
  that selected sample group together.
- Drag one arrangement chip onto a target slot to move the whole arrangement
  block.
- `Command`/`Ctrl`-click arrangement chips to build a multi-selection, then
  drag one selected chip to move those selected arrangement blocks together.
- Sample multi-selection and arrangement-chip multi-selection are intentionally
  separate modes: clicking one mode clears the other so rack intent stays easy
  to read.
- Click-select plus click-target still exists as a fallback, but drag/drop is
  now the primary interaction.
- While dragging, the rack grid now paints a ghost preview of the reordered
  sample/block positions before you release the mouse.
- Ghost previews now pulse subtly while dragging, and successful drops briefly
  fade-highlight the changed coordinates so the final reorder is easier to
  spot.
- While dragging near the rack viewport edge, the rack view now autoscrolls so
  larger racks remain editable without dropping and restarting the move.
- All rack moves use shift-neighbor semantics:
  - later occupied positions move to keep the rack block contiguous
  - free-floating holes are not the primary editing model in this baseline
- Changing `Fill` or `Blocked` reflows occupied positions through the same
  shared order-preserving engine path instead of editing coordinates ad hoc.
- `Template` is a shared shortcut layered on that same engine path:
  - `Bench rows` restores row-major fill with no blocked slots
  - `Plate columns` restores column-major fill with no blocked slots
  - `Plate edge avoidance` uses column-major fill and blocks the outer ring

Exports:

- The compact footer separates label preview from physical-carrier export
  settings and keeps both groups visible at the documented minimum Rack-window
  size.
- `Preview Labels...` in the rack window opens one in-app preview of the
  deterministic rack label sheet using the selected label-sheet preset.
  - `Export SVG...` from that preview writes the exact same label sheet once
    you are happy with it
- `Preview Labels` from the arrangement table opens the same label style but
  filtered to the selected arrangement only.
- `Physical template` in the rack window selects the first printable carrier
  family layered on top of the same rack snapshot:
  - `Storage PCR tube rack`
  - `Pipetting PCR tube rack`
  - `Cell-culture plate` (use an explicit plate profile for the intended well format)
  - `Safe-Lock 1.5/2.0 ml - collar seated (unverified)`
  - `Safe-Lock 1.5/2.0 ml - floor seated (unverified)`
  - `PCR 0.2 ml - 8x12 carrier (unverified)`
- `Fabrication SVG...` exports a top-view planning/fabrication sketch for the
  currently shown rack using the selected physical template.
- `Isometric SVG...` exports a presentation-grade pseudo-3D rack figure for
  the currently shown rack using that same physical template. This is the
  intended README/documentation-facing hero export for the current rack layer.
- `OpenSCAD...` exports one parameterized 3D model draft for the currently
  shown rack using that same physical template. Wells open above the declared
  bottom floor. Rendering is not fit, manufacturing or cell-culture validation;
  corner-radius metadata does not currently round the cube, and occupancy text
  belongs to the separate SVG/label outputs, not embossed model labels.
- `Carrier labels SVG...` exports one carrier-matched front-strip plus
  module-label SVG sheet for the currently shown rack using that same physical
  template and the currently selected carrier-label preset.
- `Carrier preset` in the arrangement table and rack window selects:
  - `Front strip + cards`
  - `Front strip only`
  - `Module cards only`
- `Simulation JSON...` exports one machine-readable rack geometry/placement
  JSON for downstream simulation adapters using that same physical template.
- Built-in label-sheet presets:
  - `Compact cards`
  - `Print A4`
  - `Wide cards`
- The `Label preset` selector above `Arrangements` and in the rack window is
  shared, so arrangement-scoped and rack-wide exports use the same engine-owned
  preset choice.
- Physical-rack exports intentionally remain rack-wide in this first baseline;
  arrangement-scoped occupancy labels/front strips still flow through the
  shared label-sheet export path instead of being baked into the 3D geometry.
- `Tube reference models (visualization only)` is a separate collapsed section
  in the Rack workspace. It offers source-backed Safe-Lock 0.5/1.5/2.0 ml
  outer-body/collar OpenSCAD exports through the shared engine operation, without
  requiring a rack or assigning a tube type to any sample. Source URLs are in
  the hover; complete dimensions/provenance accompany the exported file.
  Cap, hinge and cavity are omitted and assumptions are labelled. These are
  not printable vessels or proof of compatibility with an existing template.
- The same section includes **BRAND PCR strip 0.2 ml - 1-8 tubes (781320)**.
  **Retained tubes** defaults to eight and admits one through eight contiguous
  tubes, representing a cut strip. Export goes through `ExportTubeOpenScad`
  with `tube_count`; the selected count is included in the suggested filename.
  Capacity is 0.2 ml per tube. Only internal bridges are retained; separate caps,
  cavities and cut stubs/end tab are omitted. Bridge and taper details remain
  illustrative, not physical-fit evidence. This creates no sample inventory.
- Grid processing is bounded to 1,536 positions with overflow-checked dimensions.
  Unsupported saved racks still open unchanged; a visible explanation replaces
  the grid and processing refuses, rather than expanding or truncating it.
  The physical-template combo additionally offers Safe-Lock 1.5/2.0 ml collar-
  and floor-seated carriers, both explicitly unverified at 20-mm pitch.
  These two selections expose **Carrier pitch (mm)** and **Clear**, backed by
  the shared `SetRackCarrierPitch` operation. The optional saved override uses
  exact whole micrometres; compact exports recheck seating-specific limits and
  show **Compact spacing: cap clearance not checked**. Other physical templates
  refuse an override rather than ignore it; the field is hidden for them.
  New exports retain support mode, assumptions and nonclaims. Selecting a
  template neither assigns tube identity to samples nor qualifies fabrication,
  cap clearance, ice/freezer use or commercial-rack compatibility.
  **Safe-Lock 1.5/2.0 ml - derived storage box (unverified)** uses the chosen
  row/column grid, square 10.15-mm cells, 1.2-mm dividers and a 2-mm floor.
  OpenSCAD exports the box, insert and separate lid; SVG positions remain round
  and visibly say **Square cells drawn as round markers**. Clear any compact
  pitch override first. This is a derived design target, not measured fit or
  suitability; caps/hinges may collide. Commercial-box profiles still require
  independent dimensions/dividers/lid evidence.
  Grid formats are variable, not restricted to 9x9; changing rows/columns keeps
  the selected pitch and does not qualify a commercial-box fit.
  `PCR 0.2 ml - 8x12 carrier (unverified)` requires exactly eight rows and twelve
  columns. Nominal plate-format coordinates imply neither thermocycler use nor
  holder fit. The separate PCR reference profile retains an unknown seating
  collar; its tube-body export currently reports a structured refusal, not a
  fabricated collar. Safe-Lock tube-body exports remain available.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Serial Gel (Ladder View)](serial-gel-ladder-view.md) | [Next: Catalog Inspection And Provider Imports](catalog-inspection-and-provider-imports.md)
<!-- manual-page-nav:end -->
