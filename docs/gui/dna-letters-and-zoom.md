### DNA Letters And Zoom

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

Linear map zoom detail:

- In linear mode, GENtle uses adaptive letter routing based on viewport width
  and estimated glyph width:
  - low density: standard one-row base letters
  - medium density: helical-compressed letters (when compressed mode is enabled)
  - high density: condensed 10-row letters (when compressed mode is enabled)
  - over-capacity density: letters are hidden (`OFF`) to preserve readability
- In `HELICAL` mode, base letters now use a true oscillating phase path
  (sin/cos projection) with endpoint pinning, so the first/last visible bases
  remain anchored to viewport edges while interior overlap is distributed.
- Helical X routing uses sub-column interpolation to keep adjacent base
  positions evenly distributed and avoid local crowding plateaus that can
  distort symmetry near curve extrema.
- Helical mode also applies cycle-level arc-length compensation so consecutive
  base positions stay visually more uniform along the curve.
- When `Hide DNA backbone line when letters are shown` is enabled and bases are
  visible, the backbone line and bp tick marks are both suppressed.
- `Reverse-strand letter opacity` controls reverse-strand emphasis in both the
  linear map and the sequence-panel reverse-complement row.
- Normal DNA-map labels are independent of the Cargo build profile. Internal
  `dev` packages do not add routing-tier debug text to the map or SVG exports;
  explicit tracing/profiling remains opt-in.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Map Modes And Dotplot Workspace](map-modes-and-dotplot-workspace.md) | [Next: Feature Tree And Selection](feature-tree-and-selection.md)
<!-- manual-page-nav:end -->
