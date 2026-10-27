## Pool Distribution (Engine Ops)

<!-- manual-counterparts:start -->
Related CLI reference: [Racks, Pools and Ladders](../cli/racks-pools.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

When a strict engine operation creates multiple products (for example digest or
ligation), the Engine Ops panel shows a molecular-weight proxy distribution
based on sequence length in base pairs (`bp`).

Bin breaks are computed as follows:

- Compute `span = max_bp - min_bp` across created products.
- Choose `bucket_size` adaptively:
  - `50 bp` if `span <= 500`
  - `100 bp` if `span <= 2000`
  - `250 bp` otherwise
- Anchor buckets at `min_bp` (not around a center value).
- For each product length `bp`, compute:
  - `idx = (bp - min_bp) / bucket_size`
- Display bucket range:
  - `lo = min_bp + idx * bucket_size`
- `hi = lo + bucket_size - 1`

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Measured Gel Images (`.11` Development)](measured-gel-images-11-development.md) | [Next: Serial Gel (Ladder View)](serial-gel-ladder-view.md)
<!-- manual-page-nav:end -->
