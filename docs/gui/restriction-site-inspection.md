### Restriction-Site Inspection

<!-- manual-counterparts:start -->
Related CLI reference: [Restriction Cloning](../cli/restriction-cloning.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Selecting a restriction-site marker keeps the inline restriction expert view
  and now shows the active enzyme's raw catalog metadata
  (`recognition_iupac`, `enzyme_cut_offset_0based`, `overlap_bp`, optional
  `note`) plus a direct REBASE enzyme URL
  (`https://rebase.neb.com/rebase/enz/<Enzyme>.html`).
- The pinned restriction-site expert view includes `Copy Summary` and
  `Copy Detail JSON` actions so the SnapGene-style tooltip summary or full
  shared expert record can be pasted into protocols, bug reports, or external
  agent/MCP handoffs.
- Right-clicking a hovered or selected restriction-site label/tick on the map
  exposes the same actions directly: `Pin Restriction Site Details`,
  `Copy Restriction Summary`, and `Copy Detail JSON`. This only changes the
  interaction menu; restriction-site rendering itself stays unchanged.
- Restriction-site markers now preserve end geometry in the UI:
  - blunt cuts use aligned markers in a neutral slate color
  - `5'` overhang cuts use blue staggered markers
  - `3'` overhang cuts use amber staggered markers
  - when the linear DNA view is zoomed far enough to show sequence letters, the
    top and bottom strand cuts are drawn at their actual offset positions

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Feature Tree And Selection](feature-tree-and-selection.md) | [Next: Circular Map Labels](circular-map-labels.md)
<!-- manual-page-nav:end -->
