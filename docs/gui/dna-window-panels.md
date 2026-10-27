### DNA Window Panels

<!-- manual-counterparts:start -->
Related CLI reference: [Window and Display Intents](../cli/display.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

A DNA window is split into two visual areas:

- DNA map panel: graphical map (circular or linear)
- Sequence panel: text-oriented sequence rows

Both panels can be shown/hidden from the toolbar.
- On macOS-hosted windows, the sequence-text toggle now keeps a stable
  bottom-panel layout slot even while hidden/auto-hidden, reducing repaint
  churn and transient relayout artifacts when switching the text panel on/off.
- Hovering a visible restriction-site label or tick shows a compact shared
  restriction-site tooltip with enzyme/site count, recognition span, cut
  geometry, caret-marked top/bottom cut positions, and REBASE/note metadata
  when available; clicking the same site still opens the richer description
  panel expert view and the hover reminds users that click pins the details.
  The same presentation-ready lines are serialized as
  `tooltip_lines[]` by `inspect-feature-expert ... restriction ...` for
  scripts and external tools.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Loading and creating sequences](loading-and-creating-sequences.md) | [Next: Strand Placement](strand-placement.md)
<!-- manual-page-nav:end -->
