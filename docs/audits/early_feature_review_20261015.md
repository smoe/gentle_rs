# Early Feature Review For Internal.13

Reviewed 2026-10-15 from local `22b06a9b` plus the scoped rack-floor correction
`bb18cb44`. This is a source/documentation and OpenSCAD review, **not** executed
Rust, full native acceptance, an exhaustive bug audit or physical validation.
The inventory was cross-read against the release notes for
[.2](../release_notes/release_notes_v0.1.0-internal.2.md),
[.3](../release_notes/release_notes_v0.1.0-internal.3.md),
[.4](../release_notes/release_notes_v0.1.0-internal.4.md),
[.5](../release_notes/release_notes_v0.1.0-internal.5.md),
[.6](../release_notes/release_notes_v0.1.0-internal.6.md),
[.7](../release_notes/release_notes_v0.1.0-internal.7.md),
[.8](../release_notes/release_notes_v0.1.0-internal.8.md) and
[.9](../release_notes/release_notes_v0.1.0-internal.9.md), plus the older cloning
chapters 02/03. Those notes describe historical scope, not current acceptance.

## Entry Points, Coverage And Verdicts

All rows retain CLI/GUI-Shell/inner-agent reachability through the shared parser;
MCP `op`/`workflow` and deterministic `agent_execute_plan` reach the corresponding
engine/shell contracts. This does **not** claim one named MCP tool or a prominent
GUI button per operation. Mutations and artifact writes still require review.
Source `src/engine/tests.rs` is abbreviated **E**, `src/engine_shell/tests.rs`
**S**, `src/app/tests.rs` **A**; line references are review-time locations.
"Current" means routes/tests/docs remain connected by source inspection only.
Every Rust test and native interaction listed here remains unrun in this review.

| Feature | GUI; shared shell / operation | Deterministic coverage | Tutorial / docs | Recent substantive evidence; verdict |
| --- | --- | --- | --- | --- |
| Lineage containers | Main container table, DNA operations; `state-summary`, `DigestContainer`, `SelectCandidate` | E:37281 container creation; E:37347 selection; A:3887 split bounds | 02.01, 03.01, 03.07; `gui.md` Containers | `13377b12` single-apply ownership; **current**, native review open |
| Pools and derived aliquots | Main pool merge/select/export; `MergeContainers`, `ExportPool`, `gene-sets create-pool` | E:8065 merge copies; E:64191 source/derived lineage; S:21908 explicit preview/apply | 03.02, 03.07; `cli.md` pool contracts | collection pooling/export added after the original rack baseline; **current**, not physical stock evidence |
| Serial arrangements | Main Arrangements; `arrange-serial`, `arrange-set-ladders` | E:33758 serial; E:33816 ladder update; A:9161 order/pool description | generated 03.06, manual 03.07/03.08; `protocol.md` | serial contract retained through `ac614625`; **current**, automatic draft rack is intentional |
| Arrangement gels | Arrangement Preview Gel; `render-pool-gel-svg --arrangement` | E:35366 arrangement inputs; E:35556 comparison hints; A:9287 preview | 03.06/03.07; `cli.md` label/topology controls | `3246e774` redundant isoform identity; **drifted from early plain bands**, current docs explain markers |
| Racks and plates | Open Rack window; `racks create-from-arrangement`, `place-arrangement`, `move*`, `set-profile` | E:33864 default rack; E:33927/34024 plate selection; E:34758 appended blocks; A:9897 drag preview | 03.08 manual/hybrid; `gui.md` Rack Workspace | `4f3eb9c9` six/12/24/48-well profiles, `de7114c5` compact controls; **current**, native drag/drop untested |
| Templates and blocked slots | Rack profile/template controls; `racks apply-template`, `set-blocked`, `set-fill-direction` | E:34584 edge avoidance; E:34678 reservations; S:5993 template | 03.08; CLI/protocol templates | plate-edge policy and explicit family selection retained; **drifted docs** omitted cell-culture choice, corrected |
| Labels and carrier labels | Preview Labels / exports; `racks labels-svg`, `carrier-labels-svg` | E:34758 scoped labels; E:35190 carrier/simulation; E:35351 presets; A:9358 preview | 03.08; `cli.md` label presets | shared scope/preset path retained through `ac614625`; **current**, printer-scale check untested |
| Fabrication SVG | Rack export; `racks fabrication-svg` | E:34889 marker/clearance regression; S:6374 physical exports | 03.08; `figures/README.md` | template dimensions share the rack spec; **current**, no fit/manufacturing qualification |
| Isometric SVG | Rack export; `racks isometric-svg` | E:34889 tubes; E:35055 culture geometry; Python `test_early_feature_review` | 03.08 and retained figures | `01225ccc` culture-well representation; **current**, Rust re-export freshness untested |
| Hero SVG | Shared `racks hero-svg`; Python presentation transform | E:35131 generic slots; E:35055 plate labels; Python surface-coordinate preservation | `figures/README.md`, README figures | `01225ccc` native top-down hero supersedes using the Python transform as sole exporter; **current**, styles are intentionally different |
| OpenSCAD | Rack export; `racks openscad` | E:34973 all-template floor regression; Python semantic floor checks; six actual renders | 03.08; `figures/README.md` | **defect found and corrected** in `bb18cb44`: inverted well-floor cut. Rust regression still unrun |
| Simulation JSON | Rack export; `racks simulation-json` | E:35190 dimensions/arrangement links; S:6374 | 03.08; CLI/protocol | consumes the same physical spec, not dynamics; **current**, no external simulation acceptance |
| Digest and ligation | DNA restriction tools, routine setup; `digest`, `DigestContainer`, ligation workflow | E:8243 deterministic chain; E:8154 sticky ends; E:37447 reject unknown enzyme lists | executable 03.01/03.02; CLI | `618b15cf` reverse/IUPAC/circular repair supersedes palindromic-only behavior; **current**, rerun regressions in CI |
| Gibson | Patterns/Gibson Specialist; `gibson preview/apply` | S:6829 preview, S:6983 apply; A:5064 ordered inserts | 03.03–03.08; `gibson_planning` docs | later multi-insert/UI/lineage refinements exist; **current**, old release limitations are historical, not universal current failures |
| Protocol cartoons | Shared exports/previews; `protocol-cartoon list/render-svg/render-with-bindings` | `crates/gentle-render/src/protocol_cartoon.rs`:4992/5507 | 03.05/03.09, CLI templates | `448aac3b` explicit PCR cartoon and shared renderer extraction; **current**, visual/native review untested |
| Dotplots | DNA Dotplot map; `ComputeDotplot`, overlay render routes | E:32975 SVG; E:44360 stored view; E:44853 exon anchors | 02.03 manual, 02.04 online; CLI | .4/.5 multi/exon-anchor, .7 context-rail extensions retained; **current**, online resources optional |
| Sequencing confirmation | dedicated dialog; `sequencing-confirmation`/`ConfirmConstructReads` | E:46115 junctions; E:46908 trace; A:8001 dialog | 10.01/10.02 manual/hybrid | `063c7463` primer guidance, .4/.5 trace changes; **current**, not raw-read assembly or biological approval |
| Guide RNAs | Patterns gRNA routines; `guides filter/generate-oligos/export-csv` | E:44250 design/filter/export; S:18880 end-to-end | executable 04.04/04.05; CLI | filtering/export remains distinct from genome specificity; **current**, no efficacy/order claim |
| Construct reasoning | analysis node/window; `construct-reasoning build/show-graph/run-inspection-action` | E:58413 spans; S:4021 handoff; A:8238 opener | 03.09 manual/hybrid; CLI/protocol | .5 helper/host/adapter interpretation grew beyond early cartoons; **current**, planned constructs are not validated assays |
| Early feature experts | feature double-click/context menu; `inspect-feature-expert`, `render-feature-expert-svg`, `ui open splicing-expert` | E:28341 TFBS, E:28470 splicing, E:32422 SVG | 06.02/06.03, 07.03; GUI/protocol | `44208314`/`a02ef777` readable exact-bound UniProt presentation; **current**, not an inferred preferred transcript |

No inspected family is orphaned. The standalone set-to-arrangement link was
absent at the baseline and is the bounded Phase-6 follow-up, not a hidden pool
or GUI loop. Absence of a new native result is **untested**, never a green pass.
Last-change evidence above names substantive changes identified in the reviewed
path history/release notes, not a guarantee that every leaf last changed there.

## Actual OpenSCAD Review

Installed `/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD` exists (universal
arm64/x86_64). No installation was attempted. For each of the three committed
SCAD files, run `-o TEMP/NAME.png --imgsize=1000,750 --viewall --autocenter
--render=true INPUT.scad`. Original and corrected versions each exit **0**,
with `Status: NoError`, manifold geometry and **1000 x 750** PNGs. There are no
OpenSCAD geometry warnings. Two non-product invocation failures are retained:
sandbox Qt CPU detection aborts with required `neon` (134); the first unrestricted
command consumes the filename as `--render`'s value (usage/1). Explicit
`--render=true` outside that restriction produces the six successful renders.

Temporary images: `/private/tmp/gentle13-early-review.SmFlVn/`, suffixes
`.before.png` and `.after.png`; no images or native screenshots are committed.
Visual inspection confirms the old plate was capped and the correction opens
the wells above a solid bottom. All 54 cutter centres/dimensions remain the
same. Storage/pipetting use 24 openings, plate uses six. The examples match the
current physical-spec values (10/14-mm PCR pitch, 39.1-mm six-well pitch,
9.2/34.8-mm openings, 1.2/2.0/1.0-mm floors).

The isometric exporter and Python transform retain matching template markers
and three surface polygons; the two offline semantic tests pass. This verifies
the Python transform against retained SVG inputs, not regeneration by an unrun
Rust binary. Existing `.scad` corner-radius metadata is not used to round the
cube; label-strip depth is descriptive, while the current top recess uses
front-top clearance. These are disclosed draft-model limitations, not a new
fillet/label geometry implementation. Occupancy labels remain SVG/label-sheet
outputs, not embossed text. Rendering success does not establish printability,
tube fit, cell-culture suitability, material safety or physical inventory.

## Required External Checks

Glen/CI: run rack/container/arrangement, cloning, gel, protocol-cartoon,
dotplot, guide and expert focused tests plus full release CI, then produce fresh
exports. Native checks are listed in the [release handoff](../glen_final_candidate_execution_handoff.md#internal13-early-features).
Keep drag/drop, printer scaling, normal-size 3D/label readability and human
visual approval separate from source tests and OpenSCAD image generation.
