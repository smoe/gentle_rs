## RNA Structure (ssRNA)

<!-- manual-counterparts:start -->
Related CLI reference: [RNA Secondary Structure](../cli/rna-structure.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

When the active sequence is single-stranded RNA (`molecule_type` `RNA`/`ssRNA`),
the DNA window shows an `RNA Structure (rnapkin)` panel in the top area.

Features:

- `Refresh RNA Structure`
  - Runs `RNAfold --noPS` through shared engine APIs to fetch the dot-bracket/MFE text output and refreshes the SVG preview through `rnapkin`.
- `Export RNA SVG`
  - Saves a chosen SVG path through engine operation `RenderRnaStructureSvg`.
- Textual report
  - Shows command metadata, dot-bracket structure, optional MFE, and `stdout`/`stderr` from `RNAfold`.
- SVG preview
  - Displays rendered structure image in-panel.

Runtime dependency:

- `RNAfold` and `rnapkin` must be installed and reachable in `PATH`, or
- set `GENTLE_RNAFOLD_BIN` / `GENTLE_RNAPKIN_BIN` to explicit executable paths.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: RNA-Read Mapping Workspace](rna-read-mapping-workspace.md) | [Next: microRNA Target Scan](microrna-target-scan.md)
<!-- manual-page-nav:end -->
