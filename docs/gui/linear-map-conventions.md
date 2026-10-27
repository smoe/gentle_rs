## Linear map conventions

<!-- manual-counterparts:start -->
Related CLI reference: [Figure and Layout Exports](../cli/rendering.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Current linear map conventions are:

- Forward-strand features are shown above the DNA backbone
- Reverse-strand features are shown below the DNA backbone
- Regulatory features (for example TFBS / `regulatory_region`) are grouped into
  dedicated upper lanes above forward-strand coding features to keep dense loci
  readable
- `REG@TOP` / `REG@DNA` toggle in the map toolbar switches regulatory-feature
  placement between dedicated top lanes and a single near-baseline/GC-strip
  lane
- lane packing uses tighter spacing and stronger non-overlap padding in dense
  views
- strict `REG@DNA` mode keeps regulatory tracks in one near-baseline lane with
  enforced clearance from coding/gene lanes (no regulatory intrusion into gene
  lanes)
- Directional features use arrow-shaped ends
- Feature labels are lane-packed to reduce overlap
- Coordinate fallback labels are suppressed for unlabeled regulatory features
  to avoid clutter in dense tracks
- Restriction enzyme labels are lane-packed to reduce overlap
- Adaptive DNA-letter mode routing (linear map):
  - default mode is `Auto adaptive`
  - `Show DNA letters on linear map` is the master visibility switch for the
    base-letter layer; disabling it hides DNA letters regardless of the adaptive
    or forced routing mode while keeping feature tracks and the backbone visible
  - explicit override modes are available:
    - `Force standard`
    - `Force helical`
    - `Force condensed 10-row`
  - in `Auto adaptive`, density tiers are deterministic:
    - `density <= 1.5`: standard letters
    - `density <= 2.0`: helical letters (if compressed letters enabled)
    - `density <= 10.0`: condensed 10-row letters (if compressed letters enabled)
    - `density > 10.0`: letters hidden (`OFF`)
  - `Enable compressed DNA letters` applies only to `Auto adaptive` mode:
    disabling it forces `Auto adaptive` to use standard-only behavior.
  - row mapping uses `row = (bp + offset_bp) mod 10` where `offset_bp` is
    `Helical phase offset` (`0..9`)
  - increasing phase offset shifts the top-to-bottom seam position without
    changing base order
  - in `Condensed 10-row` mode, DNA letters replace/suppress the black
    backbone line when visible
  - condensed rows stack outward from the baseline with a deterministic minimum
    row step, so the full 10-row band remains readable in dense views
  - feature lanes are pushed outward in condensed mode to preserve readability
    around the DNA text band
  - reverse-strand letter visibility is configurable (`Show reverse-strand DNA letters`)
  - helical strand geometry is configurable:
    - `parallel` (same slant direction for forward/reverse)
    - `mirrored` (cross-over slant)
  - optional 180° reverse-letter rotation is available independently
  - the top-right map diagnostics show active mode, route policy, density,
    estimated columns fit, glyph width, and a deterministic reason string
  - Sequence-window edits are runtime-local; restart persistence still follows
    `Configuration -> Graphics -> Apply`

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Map interactions](map-interactions.md) | [Next: Circular map conventions](circular-map-conventions.md)
<!-- manual-page-nav:end -->
