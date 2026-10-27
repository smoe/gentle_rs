### Cryptic Splicing Screen

<!-- manual-counterparts:start -->
Related CLI reference: [Cryptic Splicing](../cli/splicing.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- `Cryptic screen` is a thin background adapter over
  `gentle.cryptic_splicing_screen.v1`:
  - choose a 1-based inclusive source span and orientation, or copy the current
    splicing-group span; optional insert boundaries classify vector/insert
    context and an optional CDS feature id enables coding-consequence projection
  - the scan runs on an immutable engine snapshot, so the DNA viewer remains
    responsive and project state is not mutated
  - candidate rows distinguish structural GT-AG configuration,
    branchpoint/polypyrimidine heuristics, optional MaxEnt donor/acceptor scores,
    RNA-evidence status, projected-protein-feature status, and optional CDS
    consequences. `model absent` is displayed explicitly and candidates are
    never labelled as observed junctions merely because a motif/model score
    exists
  - `use active MaxEnt snapshot` is opt-in. The panel shows the active model
    fingerprint and accepts an optional expected-fingerprint lock. If the model
    is replaced after a run, the cached report remains copyable as the exact old
    JSON, but recomputing exports/evidence joins is disabled until the screen is
    rerun
  - the collapsed evidence section accepts saved RNA-read report ids and one
    saved UniProt projection id. Both joins run as portable read-only operations;
    legacy or digest/coordinate-mismatched evidence is shown as
    `not assessable`, not silently joined
  - pair accounting shows admissible, evaluated, reported, and ranking-complete
    state. Structural caps retain a positional prefix; scored caps retain the
    globally ranked top set
  - the same compact controls are available under Sequence Tools; both surfaces
    share one result cache keyed by the effective request, source, and optional
    model content identity
  - JSON copy/save and SVG rendering use the shared portable operations and
    preserve operation/run ids in the typed report

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Conservation Workspace](conservation-workspace.md) | [Next: Isoform Architecture Panels (Engine Ops)](isoform-architecture-panels-engine-ops.md)
<!-- manual-page-nav:end -->
