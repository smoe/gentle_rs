## Map interactions

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

The DNA map supports mouse interactions:

- Hover: highlights/inspects map elements
- Click: selects a feature
- Double-click: creates a sequence selection from the clicked feature or restriction site

### Unified scroll/zoom policy

Across GUI panes, GENtle uses one interaction policy:

- Default wheel/trackpad scroll pans or scrolls content.
- `Shift + wheel/trackpad` zooms on zoomable canvases.
- `Option` (`Alt`) + drag enables hand-pan mode (`Grab`/`Grabbing` cursor).
- Arrow keys pan the active hovered region.
- `Shift + ArrowUp/ArrowRight` zoom in, `Shift + ArrowDown/ArrowLeft` zoom out
  on zoomable canvases.
- Non-zoom scroll panes (help/table/list areas) support keyboard scrolling with
  `Arrow`, `PageUp/PageDown`, and `Home/End` when hovered and no text field has
  keyboard focus.

### Linear DNA map: zoom and pan (mouse/touchpad)

When the sequence is in linear mode and your pointer is over the map:

- wheel/trackpad scroll pans the map:
  - horizontal motion pans the bp viewport left/right
  - vertical motion pans the rendered map lanes up/down
- `Shift + wheel` zooms around the current cursor position
- `Option + drag` pans the map by hand

Toolbar alternatives (linear mode):

- `-`: zoom out
- `+`: zoom in
- `Fit Seq`: reset to full sequence span and recenter vertically
- `Fit Features`: keep the current subsequence span and recenter the feature lanes vertically
- `Zoom to selection`: fit the complete current selection, including one entered
  by formula outside the visible range
- `Go to selection`: centre the selection at the current zoom level; if the
  selection is wider than the view, use `Zoom to selection` to see all of it
- `Pan` slider: move the current viewport left/right
- Right-side vertical slider (`V pan`) in the map panel: move feature lanes up/down
  directly; `0` applies vertical feature fit for the current subsequence

For a promoter and nearby first-exon region, select the intended coordinates
with `Selection formula`, press `Apply Sel`, then `Zoom to selection`. Numeric
selection formulas use 0-based, end-exclusive boundaries, whereas the `Go`
viewport fields use 1-based inclusive coordinates. For feature-relative
selection, use the intended annotated transcript/TSS rather than treating
`CDS.start` as a transcription start. A negative resolved coordinate means the
requested interval extends outside the loaded sequence: retrieve/extend the
sequence or explicitly choose a smaller interval. GENtle does not silently
truncate that biological request.

### Lineage graph: zoom and pan (mouse/touchpad)

In `Main window -> Graph` view:

- Zoom:
  - `Shift + wheel` over the graph (primary)
  - `Cmd/Ctrl + wheel` also zooms (legacy alias kept for compatibility)
  - `-`, `+`, `Reset`, `Fit`, and `Zoom` slider remain available
  - `Reset Layout` restores default node placement after manual moves
- Pan:
  - scrollbars, wheel, or touchpad scrolling
  - `Option + drag` pans by hand
  - `Space + drag` on empty background remains available as a legacy alias
- Node layout:
  - drag a node with the mouse to reposition it
  - double-click a node to open it (`pool` nodes open pool view)
  - single-click selects a node for visual focus
  - hover shows node details at pointer; pool nodes include pool range and
    auto-selected ladder hints
  - moved positions are persisted in project metadata and restored when the
    project is reopened
  - automatic graph layout uses DAG layering (parents left, children right)
    with crossing-minimizing ordering between layers
- Dense graph readability:
  - `Compact labels` reduces label density when the graph is crowded
  - operation labels are shortened and node detail labels are reduced
  - operation transitions are rendered as square glyph nodes between tube nodes
    (Petri-net style), with compact operation symbols inside the square
  - operation squares now use a consistent symbol/color set by operation family
    (for example digest, ligation, PCR, extract, filter)
  - when collapsed groups project multiple hidden operations onto one visible
    transition, the square shows a count symbol (for example `2`) and grouped
    edge text summarizes the merged operation set
  - collapsed group representative nodes display a compact hidden-internals
    badge: total hidden operation count plus top hidden operation-family chips
    (same symbol/color set as transition squares)
  - edge text labels remain secondary annotations for detailed operation names
- Workspace persistence:
  - graph node positions, zoom, compact-label toggle, graph scroll/pan offset,
    and preferred lineage/container split are stored in project metadata
  - the lineage page itself is vertically scrollable; oversized graph areas no
    longer hide the containers/arrangements section
  - in both `Table` and `Graph` view, drag the lower edge of the lineage pane
    to resize it against the containers pane

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Toolbar buttons](toolbar-buttons.md) | [Next: Linear map conventions](linear-map-conventions.md)
<!-- manual-page-nav:end -->
