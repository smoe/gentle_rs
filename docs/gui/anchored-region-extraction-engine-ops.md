## Anchored Region Extraction (Engine Ops)

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

The Engine Ops panel includes an `Extract Anchored` form for promoter-like
design constraints:

- fixed anchor:
  - absolute position (0-based), or
  - feature boundary (`kind`, optional `label`, `Start/End/Middle`, `occurrence`)
- direction:
  - `Upstream` or `Downstream`
- flexible target length:
  - `target len` with `tolerance`
- hard constraints:
  - required restriction enzyme sites
  - required TF motifs (ID/name or IUPAC)
  - optional forward/reverse primer constraints
- candidate controls:
  - `Unique`
  - `max candidates`
  - output prefix

Execution calls engine operation `ExtractAnchoredRegion` and creates one or
more candidate sequences as operation results.

Note:

- This anchored extraction uses local in-sequence anchor resolution. It is
  separate from genome-provenance anchoring used by genome extract/extend
  workflows.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Circular map conventions](circular-map-conventions.md) | [Next: Gene Set Inspector And Collection Operations](gene-set-inspector-and-collection-operations.md)
<!-- manual-page-nav:end -->
