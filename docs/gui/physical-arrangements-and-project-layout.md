### Physical Arrangements And Project Layout

<!-- manual-counterparts:start -->
Related CLI reference: [Racks, Pools and Ladders](../cli/racks-pools.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Built-in rack/plate profiles include small tube racks, six-well
  cell-culture plates, 96-well plates, and 384-well plates. Physical-template
  export includes PCR tube rack templates plus a cell-culture 6-well plate
  template for SVG/PDF documentation assets.
- Every new serial arrangement now gets one default draft rack automatically.
- `Open Lanes` now opens a chooser menu rather than immediately opening every lane window:
  - choose one lane container, or `All`
  - lane windows opened from arrangements start in a compact map-first layout with the sequence text panel hidden by default
- The lineage (`Table`/`Graph`) area is split from `Containers` with a
  draggable horizontal divider in the main window.
- In `Table` view, the lineage grid supports both horizontal and vertical
  scrolling; `Node`/`Op` cells use compact IDs with full values on hover.
- The Lineage graph drawing surface uses a dedicated science-canvas frame with
  mode-aware light/dark colors, separate from surrounding controls.
- In DNA sequence windows, the same canvas treatment applies only to the
  graphical map canvas; the sequence text panel and toolbar/control rows remain
  on their normal control surfaces.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Project Overview And Collections](project-overview-and-collections.md) | [Next: Project And Sequence Window Examples](project-and-sequence-window-examples.md)
<!-- manual-page-nav:end -->
