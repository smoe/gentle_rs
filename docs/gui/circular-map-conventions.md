## Circular map conventions

<!-- manual-counterparts:start -->
Related CLI reference: [Figure and Layout Exports](../cli/rendering.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Current circular map conventions include:

- Features arranged around the circular backbone
- Restriction enzyme labels around the perimeter
- Optional overlays for GC content, ORFs, and methylation sites
- Feature geometry is internally 0-based and end-exclusive. The last base is
  retained when a feature ends at the sequence boundary; hover coordinates
  are shown as 1-based inclusive positions (`1-1 bp` denotes the first base).
- An origin-crossing ORF is drawn across the origin on either strand, not
  around the opposite arc. This is an ORF prediction, not evidence of expression.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Linear map conventions](linear-map-conventions.md) | [Next: Anchored Region Extraction (Engine Ops)](anchored-region-extraction-engine-ops.md)
<!-- manual-page-nav:end -->
