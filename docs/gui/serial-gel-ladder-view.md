## Serial Gel (Ladder View)

<!-- manual-counterparts:start -->
Related CLI reference: [Racks, Pools and Ladders](../cli/racks-pools.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

When pool members are available, Engine Ops shows a virtual agarose-gel preview:

- one or two DNA ladders are auto-selected to span the pool bp range
- ladders are drawn in lanes next to a pool lane (bands from pool members)
- pool-band labels show bp, with multiplicity when multiple members have the
  same length

Serial gel export is available in two places:

- Engine Ops (`Export Pool Gel SVG`) for current sequence-id inputs.
- Main lineage page:
  - `Arrangements` table: `Preview Gel` opens an in-app serial-gel preview with
    live left/right ladder selectors, a rendered gel image, `Save to
    Arrangement`, and `Export SVG...`.
    Preview and gel-export controls are enabled only for `Serial`
    arrangements; `Plate` arrangements remain available for rack/plate work
    and show the engine-owned serial-only readiness explanation on the disabled
    gel controls.
  - `Containers` table: `Gel SVG` exports one lane for that container;
    `Export Pool` calls `collections run export-pool` and writes one
    `gentle.pool.v1` artifact. The action is disabled for empty or
    non-exclusive containers and for in-silico `Selection` containers, matching
    the shared engine preflight. After the native save dialog closes,
    validation, serialization, and writing run in a
    detached engine snapshot; the GUI reports success only after decoding and
    validating the typed export report against the requested container and
    destination path.
  - `Arrangements` table: `Export Gel` exports all lanes defined in that
    arrangement, using the arrangement's saved ladder pair when present or auto
    ladder selection otherwise. The canonical collection policy classifies
    this as context-agnostic `arrange`, preserving the declared lane order.
  - `Arrangements` table: `Open Lanes` lets you open one chosen lane or all lanes as compact DNA windows instead of always opening every lane at once.
  - `Arrangements` table: `Open Rack` opens the linked physical rack draft for
    that arrangement, creating the default rack on demand for older/legacy
    arrangements that do not yet have one.
  - `Arrangements` table: `Preview Labels` opens one in-app preview of the
    deterministic label sheet scoped to that arrangement on its linked rack
    draft.
    - the same preview window then exports the SVG once the cards look ready
      to print
    - arrangement-scoped label cards now prioritize:
      - rack position
      - role
      - sequence or ladder identity
      - bp length/topology when sequence-backed
    - long arrangement/op provenance no longer overflows each card; the scope
      is shown in the sheet header instead
  - `Arrangements` table: `Carrier SVG` exports one carrier-matched front-strip
    plus module-label SVG for that arrangement from its linked rack draft,
    using the currently selected physical template and carrier-label preset.
  - `Arrangements` table: `Place on Existing Rack...` appends the arrangement
    as one contiguous block onto another saved rack.

Tutorial companion:

- [`docs/tutorial/03-07_gibson_arrangements_gui.md`](../tutorial/03-07_gibson_arrangements_gui.md)
  walks through how a Gibson apply step creates singleton output containers,
  records one reusable serial arrangement, and then exports the arrangement as
  one ladder-flanked gel.

Ladder source:

- built-in ladder catalog: `assets/dna_ladders.json` (derived from historical
  GENtle ladder data; upstream legacy files used "marker" naming)
- built-in RNA ladder catalog: `assets/rna_ladders.json`
- historical references:
  - `https://github.com/GENtle-persons/gentle-m/blob/main/src/marker.txt`
  - `http://en.wikibooks.org/wiki/GENtle/DNA_markers`

Controls:

- `ladder preset`: quick selection for common ladder pairs (or Auto)
- `gel ladders`: optional comma-separated ladder names
  - if empty, ladder selection is automatic
- `Export Pool Gel SVG`: writes the current ladder + pool band view to SVG via
  shared engine operation `RenderPoolGelSvg`
- `Arrangement Preview Gel`:
  - left/right ladder selectors can be changed independently
  - changing one side while the other stays `Auto` coerces to the same named
    ladder on both sides, so single-ladder selection remains simple
  - `Auto` on both sides previews/export uses shared engine auto ladder
    selection
  - arrangement gel preview now also exposes one shared gel-condition bundle:
    - agarose %
    - buffer (`TAE` / `TBE`)
    - topology-aware circular migration on/off
  - topology-aware mode now distinguishes explicit circular-form hints like
    `supercoiled`, `relaxed circular`, and `nicked/open circular` when those
    terms appear in sequence names/definitions/comments; otherwise circular DNA
    stays on the generic circular heuristic
  - nearby fragments that co-migrate now collapse into one merged-band
    annotation instead of drawing several indistinguishable bars
  - when the lane roles look like `vector` / `insert` / `product`, the
    right-hand detail panel adds `Comparison hints` for insert sizing and
    vector-vs-product interpretation
  - matching role badges are drawn under the lane labels, so even generic lane
    names like `Tube A/B/C` still read as vector/insert/product immediately
  - the right-hand detail panel adds `Merged-band notes` so grouped bands are
    explained as observed-vs-actual lane readouts
  - `Save to Arrangement` persists the currently selected ladder pair for later
    `Export Gel` reuse
  - condition changes refresh the preview immediately and are reused by
    `Export SVG...` from that preview window

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Pool Distribution (Engine Ops)](pool-distribution-engine-ops.md) | [Next: Rack View](rack-view.md)
<!-- manual-page-nav:end -->
