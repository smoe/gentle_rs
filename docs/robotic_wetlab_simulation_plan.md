# Robotic Wet-Lab Simulation And Equipment Digital Twins

Status: future post-0.1 roadmap track. This document records the intended
ownership and evidence boundaries; it does not authorize physical robot
execution or add a release gate to the current internal series.

Pilot scope, 2026-11-04: the owner authorized the narrowed test-only emission
regression and retained the existing row-A / label-strip orientation. One new
rack test and the unchanged reference-tube tests have scoped source-emission
verification, recorded [below](#test-only-emitter-pilot-2026-11-04). This adds no
exporter, simulator, CAD language or runtime report and establishes no native
CAD, physical or scientific acceptance.

## Goal

Make selected GENtle protocols executable as provenance-bound simulations in a
robotics environment such as Gazebo. The first benchmark should be an
automated CUT&RUN workflow in a synthetic reference laboratory. That is the
first end-to-end protocol benchmark, not the first CAD-conformance gate. The
simulation must expose material flow, equipment use, timing, robot reach and collision
constraints, checkpoints and recoverable failures without claiming that a
virtual run validates biological performance.

The longer-term aim is to test protocol and laboratory-automation designs in a
digital world before committing scarce samples, equipment time or physical
robot motion.

## Ownership Boundaries

- **GENtle** owns versioned protocol semantics, material and sample identity,
  dependencies, quantities and timing claims, approval state, expected
  checkpoints, simulation requests and returned receipts.
- **A separate shareable equipment-asset project** should own reusable models
  of labware and instruments. GENtle references exact catalog records and
  content digests instead of embedding a second geometry library.
- **Gazebo or another simulator** owns world loading, physics, collision,
  sensors, robot/controller execution and simulator-native traces.
- **The real laboratory world** remains private local infrastructure. A public
  synthetic reference world and redistributable equipment fixtures provide CI
  and tutorial coverage without exposing room geometry, labels, screens,
  access paths or other sensitive details.

Conversation and planning agents may propose a simulation, but they do not
gain authority to move a physical robot, mutate GENtle state or publish a
private world. Existing preview, confirmation, freshness and rollback rules
remain in force.

## Equipment Asset Contract

Parametric CAD tools such as SolidCraft may author editable source models, but
the shared contract must be tool-neutral. Each released equipment or labware
asset should bind at least:

- stable asset and revision identifiers, licence, author/source provenance and
  content hashes;
- source CAD plus named exchange artifacts where appropriate (for example STEP,
  glTF/mesh and Gazebo SDF or URDF);
- physical dimensions, scale and coordinate frames; collision/inertial
  geometry, joints and safety envelopes only when separately supported by
  reviewed evidence, otherwise explicitly unavailable;
- semantic interaction points such as deck slots, tube or plate positions,
  lids, doors, grippers, pipette tips, liquid ports and observable states;
- supported actions, preconditions, resource requirements and known fidelity
  limits, with no inferred action capability when these are unknown;
- calibration or local override records kept separately from the public base
  asset.

Photogrammetry or descriptive room models are useful visual references, but
must not be used directly as authoritative physics geometry. Derive simplified,
reviewed semantic and collision objects for simulation.

### Canonical Abstract Equipment Representation

Use a versioned, human-readable intermediate representation as the canonical
source for equipment families instead of making any CAD application's native
file authoritative. Deterministic generators may rewrite that representation
into SolidCraft- and FreeCAD-targeted source or projects and then produce
exchange and simulator artifacts. Generated files bind back to the source
description, generator version and parameters by digest.

Parametric primitives, reusable subassemblies, constraints, materials and
interaction points are longer-term requirements, not a mandate to build a
general CAD language in the next slice. Reuse existing source-backed profiles
and template derivations first; add a new abstraction only for a demonstrated
consumer. Richer families may eventually describe:

- tubes describe body, base, rim, opening, cap or lid, nominal capacity,
  material, grasp region, liquid frame and observable open/closed state;
- racks describe outer envelope, rows, columns, pitch, offsets, slot geometry,
  orientation and identifiers, nesting rules, grasp regions and the transform
  from every slot to its contained tube;
- variants change parameters or compose referenced parts rather than copying
  opaque meshes;
- exporters must preserve and test dimensions, scale, handedness, frame names,
  slot count, pitch and semantic interaction points across supported targets.

Round-trip editing in every CAD tool is not a prerequisite. When a generated
model is manually refined, the change must either be expressed back in the
abstract source or recorded as a derived, tool-specific override with explicit
provenance; silent divergence is not acceptable.

The shipped `gentle.rack_equipment_description.v1` remains a closed geometry
description. Preserve its schema/projection constants, existing bytes, unknown
mass/inertia and empty permitted actions. Future physics, interaction frames
or capabilities require a separately reviewed new schema/version or companion
record; do not populate v1 placeholders or reinterpret slot origins as poses.

## Protocol-To-Simulation Boundary

Do not translate free-form protocol prose directly into robot commands. Add a
versioned, simulator-neutral execution plan that binds:

1. the reviewed GENtle protocol and exact inputs;
2. typed material transfers, incubations, washes, measurements and human
   checkpoints;
3. required equipment capabilities and resolved asset revisions;
4. spatial frames and consumable placements supplied by the selected world;
5. permitted robot actions, timing assumptions and explicit unresolved steps;
6. a preview digest requiring approval before simulator execution.

The simulator adapter returns a receipt bound to that plan and world. It should
record simulator/version, robot/controller and asset hashes, initial state,
executed and refused actions, simulated time, collisions or limit violations,
material-lineage transitions, checkpoint outcomes, injected faults and output
artifact hashes. Simulator-reported facts remain observations; GENtle derives a
separate assessment with explicit reasons.

Before any receipt importer is implemented, its contract must also define:

- stable plan-step and execution-attempt identifiers, exact input/plan/world/
  asset bindings, and explicit terminal outcomes for every required step;
- execution status separately from evidence completeness and GENtle's
  `pass`, `fail` or `incomplete` assessment; a completed simulator process is
  not automatically a successful protocol;
- refusal or non-passing findings for missing/conflicting evidence, duplicate
  step identifiers within an attempt, mismatched bindings, invalid dependency
  order and failed processes. A valid empty result must be distinguishable
  from absent execution evidence;
- idempotent imports and explicit retry/checkpoint rules: identical repeated
  evidence cannot transfer material twice, and conflicting retries cannot
  silently replace accepted evidence;
- explicit units and material balances, including declared losses/waste.
  Nominal tube capacity is not an observed liquid volume;
- bounded validation and a separate assessment/preview before any confirmed,
  freshness-checked, atomic engine update. Retain partial-run observations;
  do not imply that physical or simulated transfers were undone just because
  receipt acceptance failed. Missing evidence cannot mean zero collisions or
  completed lineage transitions.

These are requirements for the later protocol contract, not new public enums,
receipt schemas or lineage mutations authorized by the CAD-source pilot.

### Determinism And Replay Levels

- GENtle-generated CAD source and canonical input bindings must reproduce
  byte-for-byte for pinned inputs, generator and parameters. Declare any
  volatile diagnostic fields separately; do not erase scientific input fields
  or change existing byte contracts to make a comparison pass.
- Geometry conformance uses named, predeclared tolerances and explicit units,
  accounting for current decimal formatting and f32 widening. Numerical
  tolerance is not manufacturing clearance or a physical-fit claim.
- Later simulation receipts must bind random seeds, timestep, physics backend,
  plugins, controller/configuration digests and relevant runtime/platform
  versions. State which logical outcomes must agree and which numerical
  measurements use tolerances. Do not promise bitwise-identical trajectories
  across unspecified platforms or runtime versions.
- Exact plan/binding checks remain deterministic even when numerical replay
  is unavailable. Source generation, native CAD acceptance, simulator replay
  and physical/scientific acceptance are separate verdicts.

## Staged Development

1. **One-family CAD-source conformance.** Execute only the bounded proposal
   below: one existing rack family and independent tube-emitter checks, one
   OpenSCAD-source backend and ordinary test diagnostics. No new engine record,
   schema, CAD representation, backend, catalog, execution plan or world is
   required.
2. **Shareable catalog and justified extensions.** Resolve repository,
   licensing, review and local-override policy before new asset publication.
   Add another family or CAD backend only for a named consumer and its own
   conformance gate. SolidCraft and FreeCAD remain candidates, not simultaneous
   prerequisites. Reconcile any abstract representation with existing models.
3. **Minimal synthetic reference cell and protocol contract.** After a reviewed,
   quantity/sample-bound synthetic protocol input exists, define the neutral
   execution/receipt contract above and test one bounded operation before a
   whole protocol. Bind one world/robot and the minimum equipment needed;
   loading, placement and collision checks use the declared replay levels.
4. **CUT&RUN benchmark.** Map a reviewed CUT&RUN protocol to typed operations,
   including consumable and sample lineage, timing and explicit manual or
   unsupported steps. Exercise success, missing-equipment, unreachable-pose,
   collision, depleted-consumable and interrupted-run cases.
5. **Private real-lab world.** Bind the locally reconstructed room and equipment
   inventory without moving its geometry or sensitive metadata into public
   fixtures. Compare planned paths and resource schedules with the synthetic
   benchmark.
6. **Physical rehearsal later.** Only after separate robotics safety review,
   calibrated hardware identities, emergency-stop and human-presence policy,
   execute bounded dry motions or hardware-in-the-loop tests. Physical wet-lab
   execution remains a separately approved capability.

## Acceptance Boundaries

The first CAD-source gate is the one-family acceptance checklist below. It
does not require a public catalog, a second CAD tool or a simulator.

The later end-to-end simulation benchmark separately requires:

- one public equipment catalog and synthetic Gazebo world reproducible from
  exact hashes;
- one reviewed tube-and-rack family with exporter conformance against the
  same bound source; a second backend requires separate justification;
- one CUT&RUN simulation plan that can be previewed, approved, replayed and
  refused under the declared exact-binding and numerical-replay rules;
- preserved sample/material lineage and protocol checkpoints across the run;
- negative tests for stale assets/worlds, scale or frame mismatch, missing
  capabilities, collision and interrupted execution;
- no private laboratory geometry or identifying metadata in public artifacts;
- a receipt that clearly distinguishes simulated completion from physical
  execution and biological validation.

Native GUI visualization is optional. A conversational or MCP caller may drive
the workflow, while Gazebo and ordinary 3D viewers provide the spatial
inspection surface.

## Open Maintainer Decisions

- repository, licence and review policy for the shared equipment-asset catalog;
- schema and implementation home for the canonical abstract equipment
  representation and its SolidCraft/FreeCAD generators;
- canonical exchange formats and the minimum fidelity required for collision,
  manipulation and visual models;
- protocol-operation granularity and the boundary between GENtle, a workflow
  scheduler and robot-specific controllers;
- first robot, pipetting and instrument interfaces represented by the public
  reference cell;
- retention, access and publication policy for private world and run artifacts.
- **D1 - Pilot orientation resolved, 2026-11-04:** the owner confirms retaining
  the existing row-A / label-strip edge as GENtle's "front" for the test-only
  pilot. Future ANSI/SLAS deck mapping or legacy-orientation migration remains
  separately reviewed; no frames or labels change. See the
  [durable decision](decisions.md#equipment-descriptions-do-not-grant-execution-authority).
- **D2 - Published file mode:** tempfile-based no-clobber exporters publish 0600
  on Unix, unlike ordinary `fs::write` exports. Set repository-wide publication
  policy separately; do not change one exporter in isolation.

## Next Slice: One Family, One Exporter

Status: implemented as one test-only emitter regression, with the scoped
[verification receipt](#test-only-emitter-pilot-2026-11-04). This is not
geometric validation or completion of a simulation milestone. The original
[read-only Claude review prompt](robotic_wetlab_plan_review_prompt_20261104.md)
covers the source-review follow-up and remaining decision gates.
The separately authorized [native inspection](#native-openscad-inspection-2026-11-07)
records successful rendering and two unresolved parameter/geometry discrepancies.

### Fixed Working Set

- Family: the existing `safe_lock_tube_rack_floor_seated` Safe-Lock carrier,
  with `RackProfileSnapshot::from_kind(RackProfileKind::SmallTube4x6)`, fixed
  `rack_id="conformance_safe_lock_floor_4x6"`, no placements, no compact-pitch
  override and `blocked_coordinates=[]`. There is no sample data.
- Reference tube: `eppendorf_safe_lock_1_5ml`. The carrier is already derived
  jointly from the 1.5 and 2.0 ml profiles; retain both as bound derivation
  inputs rather than silently changing its geometry. The 2.0 ml profile is not
  a second rendered model or an expansion to another carrier family.
- Backend: OpenSCAD source only, reusing `ExportTubeOpenScad` /
  `gentle_render::labware::render_tube_openscad` and
  `racks openscad` / `export_rack_openscad`. A tube file and rack file are two
  artifacts from one backend, not two new exporters or command vocabularies.
- Result: test assertions and ordinary test-runner diagnostics only. Defer a
  portable conformance record, runtime reader, report ledger and schema until
  a real runtime consumer needs them; do not overload equipment v1.

The tube and rack checks are independent. The rack source does not reference
or place the exported tube; selecting the 1.5 ml reference does not change the
carrier. No cross-artifact pairing, assembly, fit or insertion invariant is
being tested. The family relationship is shared source provenance, not a
demonstrated physical tube/rack relationship.

The zero-blocked baseline must emit 24 hole cuts. Preserve the existing
distinction outside that baseline: rack OpenSCAD omits holes at blocked
coordinates, while equipment v1 retains those physical slots with
`blocked=true`. Do not test or describe those outputs as equivalent for a
nonempty blocked set, or change either behavior to force agreement.

The reference tube remains the existing outer-body/collar visualization.
Caps, hinges, cavities, liquid/grasp frames and a calibrated slot-to-tube
insertion transform remain unavailable. Do not make a printable vessel, infer
fit from combined renderings or add an assembly pose to fill these gaps.
PCR96, the undocumented individual PCR seating collar, strips, collar-seated
carriers, storage boxes and other tube variants are outside this pilot.

### Admission And Compatibility Gates

- D1 is owner-confirmed on 2026-11-04: retain and pin the existing row-A /
  label-strip-edge convention for this pilot. This is not an ANSI/SLAS deck
  mapping or approval to migrate existing frames. A different orientation
  requires a separate decision and implementation scope.
- Generate test artifacts only in owned temporary directories; commit no new
  `.scad`, drawing, mesh or JSON snapshot. This avoids a new fixture-publication
  gate, not the existing source-use/provenance obligations. Document the inline
  synthetic setup and existing drawing sources. Any later committed artifact
  needs separate rights review; drawing digests are not a licence grant.
  No catalog publication, private lab geometry or sample data.
- Reuse the existing source/profile/template ownership and existing command
  routes; do not introduce a generic constraint solver or CAD geometry DSL.
- D2 remains a separate repository-wide decision. Do not change legacy
  overwrite/file-mode behavior or introduce a new publication path in this
  pilot; use fresh owned temporary destinations for conformance execution.
- Preserve equipment v1, its projection and pinned hashes, legacy output bytes,
  project content/history and existing adapter approval rules. Any needed
  geometry, writer-policy or shipped-schema change is a separate reviewed slice.

### First-Gate Acceptance And Verification

1. Keep the fixed inputs above and use the existing tube/rack operation paths;
   artifact writes keep their normal confirmation requirements. Reuse the
   existing tube dimension/provenance/repeatability tests, not a duplicate suite.
2. **Emission check, not independent derivation verification:** parse the
   controlled generated rack source and compare its scalar declarations and
   all 24 row-major hole translations with the emitted values expected from
   the current spec. Expected spec and exporter share derivation helpers;
   agreement alone proves only emission/formatting consistency. Do not build a
   general OpenSCAD parser or accept arbitrary user CAD input.
3. **Independent literal oracle:** additionally pin the constants and arithmetic
   below in test code without calling the production derivation for expected
   values. This catches shared-helper drift against the declared source/rules,
   not physical correctness. Keep existing derivation tests rather than adding
   a second derivation implementation.
4. Assert byte-identical repeated rack source with the fixed `rack_id`, profile
   and generator inputs, using fresh temporary destinations. Reuse the existing
   tube repeatability/profile-JSON checks. No committed artifact freshness check
   is proposed; if introduced later, scoped LF policy, LF/CRLF coverage and
   fixture provenance are mandatory. The full embedded tube `profile_json`
   remains part of any future artifact-publication review.
5. Retain existing missing-dimension, invalid-geometry and export-refusal tests;
   do not add a runtime validator just to create new refusal/staleness cases.
   Missing physical data remains unavailable. Equipment-v1 unit conversion,
   hash/validation and blocked-slot behavior keep their existing regressions.
6. Probe installed CAD availability/version only when execution is separately
   authorized. A controlled native parser/load check is a separate verdict;
   absent tooling is unverified, not a pass or permission to install/download.
   On macOS check the direct OpenSCAD app binary before suggesting installation.
7. Preserve legacy bytes, equipment v1 and existing project/history/approval
   behavior through existing regression tests. No GUI, backend vocabulary,
   command route, engine report or screenshot campaign is required.

**Literal oracle and fixed numerical rules.** All values below are millimetres.
The existing profile data gives body diameter 9.85, lengths 38.9/40.0 and
cylindrical lengths 20.0/32.3. Under the existing engineering assumptions:

- opening: `9.85 + 2 * 0.30 = 10.45`;
- bore: `floor(max(38.9 - 20.0, 40.0 - 32.3) + 3.0) + 1 = 22.0`;
- rack height: `22.0 + 2.0 = 24.0`; pitch: `20.0`;
- edge: `(20.0 - 10.45) / 2 = 4.775`;
- width: `2 * 4.775 + 10.45 + 5 * 20.0 = 120.0`;
- depth: `10.0 + 4.775 + 10.45 + 3 * 20.0 = 85.225`;
- hole coordinates: row-major Cartesian product of
  `x=[10, 30, 50, 70, 90, 110]` and `y=[15.225, 35.225, 55.225, 75.225]`;
  iterate rows (`y`) outside columns (`x`). A1 is `(10, 15.225)` and D6 is
  `(110, 75.225)`.

Pin the canonical three-decimal declarations and hole translations, not just
their presence. In particular:

```text
outer_width = 120.000;
outer_depth = 85.225;
rack_height = 24.000;
opening_diameter = 10.450;
floor_thickness = 2.000;
front_top_clearance = 10.000;
front_label_strip_depth = 6.000;
front_label_strip_recess = 0.800;
corner_radius = 2.000;
```

`SCAD_EMISSION_ABS_TOL_MM = 0.000500001` covers three-decimal rounding plus a
`1e-9` mm floating-point comparison guard when comparing parsed emission to
the already widened f32 spec. Do not add the spec's f32 error a second time.
`SPEC_LITERAL_ABS_TOL_MM = 0.0001` separately compares spec fields to the
hand-calculated literals, matching the existing carrier-test allowance for
f32 arithmetic. Exact string pins remain exact; neither tolerance weakens
those checks or represents fabrication clearance. No adaptive tolerance or
new unit conversion is needed; changing either constant requires review.

**Increment beyond existing coverage.** In `src/engine/tests.rs`,
`labware_carriers_preserve_legacy_json_and_export_both_support_modes` already
pins 4x9 JSON derivation values and SCAD disclosures;
`labware_legacy_carrier_exports_omit_drawing_disclosures` protects legacy
headers; `export_rack_fabrication_isometric_svg_and_openscad_include_template_markers` checks basic
SCAD/module/floor output. `rack_openscad_all_templates_preserve_bottom_floor_and_open_top`
already covers the floor/open-top structure. The render crate's
`tube_models_bind_their_source_dimensions_and_disclose_omissions` already pins
tube literals, repeated bytes and embedded profile JSON. Reuse those unchanged.
The new value is the fixed 4x6 rack's parsed numerical declarations, exact
hole count/order/coordinates and repeatability; do not relabel existing tests
as newly implemented coverage. Native CAD and full geometric validation remain
unverified by these source assertions.

Passing this gate establishes only the specified emitter regression checks.
Native CAD acceptance, mesh geometry, collision, inertia, actions, protocol
execution, Gazebo/worlds, physical fit and biological performance are separate
or deferred. A reviewed CUT&RUN input is needed for the later protocol stage,
not for this sample-free geometry pilot.

### Test-Only Emitter Pilot, 2026-11-04

Owner D1 decision: retain the existing row-A / label-strip-edge convention
for this pilot; no legacy coordinates, labels or production geometry change.
The durable decision is in `docs/decisions.md`.

`src/engine/state/lineage_containers/equipment/tests.rs` adds only
`rack_equipment_safe_lock_4x6_scad_emission_matches_spec_and_literals`.
The existing equipment test module can call the actual private geometry
helpers for the emission/spec comparison without exposing another API; the
separate literal assertions remain independent of those helpers. The test
uses the normal `ExportRackOpenScad` operation, parses nine declarations and
all 24 hole translations, pins the agreed constants/tolerances and checks
repeat bytes plus unchanged project/undo/redo state. Outputs are temporary;
no CAD, drawing, mesh or JSON fixture is committed.

Source: `c03cdca8` plus this working-tree pilot on macOS arm64,
`rustc 1.100.0-beta.7 (fd8e211ae 2026-10-17)`. Verification used:

```sh
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
cargo test -q --locked --offline -p gentle-render --lib tube_models_bind_their_source_dimensions_and_disclose_omissions -- --test-threads=1
cargo test -q --locked --offline --no-default-features --lib rack_equipment_safe_lock_4x6_scad_emission_matches_spec_and_literals -- --test-threads=1
cargo test -q --locked --offline --no-default-features --lib -- rack_equipment_ labware_ export_rack_ rack_openscad_ --test-threads=1
cargo test -q --locked --offline -p gentle-render --lib labware::tests -- --test-threads=1
cargo check -q --locked --offline --no-default-features
rustfmt --edition 2024 --config skip_children=true --check src/engine/state/lineage_containers/equipment/tests.rs
python3 scripts/manual_docs.py
python3 scripts/check_documentation_status.py
git diff --check
```

The new test passes alone; the 43 selected engine regressions and four existing
render-labware tests all pass without skips. The separately run tube and new
rack tests overlap those groups, not additional unique coverage. Headless
Cargo check, scoped formatting, documentation and whitespace checks pass.
Existing headless dead-code and linker compact-unwind warnings remain.
The implementation was delivered in `88cf72b6`; the source description above
records its pre-commit execution, not a later committed/merged-source rerun.

The tube tests, production helpers/exporters, equipment v1, schema/projection
constants, dependency manifests and lockfile are unchanged. No native CAD
application or simulator was launched, and no physical fit, pairing or
manufacturing claim follows. Native Windows, default-feature/GUI and full
workspace execution were not run. D2 and all later simulation/physical gates
remain separate. The original implementation ledgers below are historical;
their results are not relabelled as this pilot's verification.

### Native OpenSCAD Inspection, 2026-11-07

The owner authorized a fast-forward to local main and independent OpenSCAD
inspection, without installing tooling or changing the exporters. The branch
advanced to `6c1edcf3`, then to `1cf7e120` when main advanced during inspection;
the latter changes unrelated observed-structure code, not these labware paths.

Execution: macOS 27.2 arm64, installed
`/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD`, version `2026.04.04`,
explicit `Manifold` backend and `--hardwarnings`. Sandboxed Qt probes abort
with a NEON error while `sysctl` is denied; the same binary works outside the
sandbox. This is not evidence of an absent installation or unsupported CPU.

No Rust build was run. The existing `target/debug/gentle_cli` embeds source
`89df315a80e0d3ae02e5685efd7a94e3895fbaa7` and has SHA-256
`f2eb91da15d57acdcea70e396f836786b3c9d77ecaeb1b5edf6bef8f400a0332`.
Byte comparison confirms that `src/engine/state/lineage_containers.rs`,
`crates/gentle-render/src/labware.rs` and
`crates/gentle-protocol/src/labware.rs` are unchanged from that source; the
`data/labware` sources are also unchanged. This is an existing-binary CAD
inspection, not a compilation or test receipt for current main.

The real CLI `ExportRackOpenScad` and `ExportTubeOpenScad` operations generated
the fixed unblocked, unoccupied 4x6 floor-seated rack above and the independent
1.5 ml reference tube. The synthetic state starts from CLI `export-state`,
with only that rack added; no experimental project is loaded. Repeated exports
are byte-identical and leave the parsed synthetic project unchanged. SCAD hashes:

- Rack: `cb4d463c9b9b90eecabfabd855641d07a691fd163c2776a22355c0c96765d81c`.
- Tube: `3733d578bd4eb438e5a33ab51bafa3b641cd42e3f1c66dd37500f47359e5bc60`.

Native mesh generation, three rendered PNG views and seven mesh-derived SVG
sections complete without OpenSCAD warnings/errors. Representative commands,
from the owned inspection directory, are:

```sh
SCAD=/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD
"$SCAD" --hardwarnings --backend Manifold -o rack.off rack.scad
"$SCAD" --hardwarnings --backend Manifold -o tube.off tube.scad
"$SCAD" --hardwarnings --backend Manifold --render --projection=o --imgsize=1400,1000 --camera=0,0,0,55,0,35,300 --autocenter --viewall -o rack.png rack.scad
"$SCAD" --hardwarnings --backend Manifold -o rack-section-bores.svg rack-section-bores.scad
```

Each section wrapper uses `projection(cut=true) translate([0,0,-Z])
import("rack.off");`: the actual exported mesh, not a duplicate rack model.
Independent inspection finds:

- Rack bounds `120 x 85.225 x 24` mm; tube bounds `10.8 x 10.8` mm with
  top at `38.9` mm and faceted tip at `0.000722726` mm. Each mesh is finite,
  one connected surface, with every edge incident to two faces.
- All 24 bores have diameter `10.45` mm and the literal Cartesian-product
  centers above. A section at `z=1.999` has only the outer contour; at
  `z=2.001` it has that contour plus 24 holes, bracketing the 2 mm floor.
- The label pocket begins between `z=23.199` and `23.201`, consistent with
  its 0.8 mm recess. Its bounds are `x=2.388..117.613`, `y=0.8..9.2`.
  Row A is adjacent to that low-Y strip, preserving D1. The isometric/top
  rack PNGs and separate tube PNG were visually inspected; SCAD emits no
  embossed row/column labels and no tube/rack assembly.

**Two existing discrepancies remain open.** `corner_radius=2.000` is declared
but unused: mesh sections have square outer corners. Likewise,
`front_label_strip_depth=6.000` is declared but unused: the actual pocket is
`8.4` mm deep in Y (`0.8..9.2`), derived from `front_top_clearance - 1.6`.
Thus native parsing/rendering and the listed geometry checks pass, but full
agreement between the declared scalars and the solid is not established.
Do not change legacy geometry/bytes or the pinned test to conceal this;
resolve the intended semantics in a separate owner-reviewed geometry change.

Local evidence is temporary, not a committed fixture or a runtime report:
`/private/tmp/gentle-openscad-20261107.XZiKSW/` contains SCAD/OFF/PNG/SVG,
`inspect_models.py`, per-command arguments/stdout/stderr and `inspection.json`
with bindings/hashes. The first local mesh reader incorrectly assumed OFF
counts occupied a separate line; it was corrected for valid same-line OFF
headers before the successful checks. No GENtle code was changed.
No new fixture publication, installation, screenshot capture, GUI interaction,
native Windows/Linux, physics, manufacturing, physical-fit or scientific
acceptance is claimed; D2 remains open. Existing Rust receipts stay separate.

## First Implementation Slice

2026-11-03: the three forward-ported Glen commits (`8d91d696`, `5ae0b514`,
`4a633418`) were imported intact from his checksum-verified incremental bundle
based on `888e49ce`. The first implementation deliberately reuses existing
GENtle tube/rack drawings and derivations rather than choosing a new catalog
repository/licence, universal CAD language or simulator/controller prematurely.

`racks equipment-json RACK_ID OUTPUT.json [--template TEMPLATE]` exports a
separate, typed `gentle.rack_equipment_description.v1`. It projects the exact
existing carrier geometry into explicit right-handed local metre frames,
carrier-top slot origins and an unverified body envelope. Every physical
slot remains present, including blocked slots. Embedded source geometry,
drawing provenance/unknowns, template revision, source and report hashes make
the result independently inspectable. Sample identity, occupancy and private
world geometry are not exported; rack ID/blocking are explicit instance data.
The same description returns through the operation result. Shared-shell,
generic MCP and agent approval use existing engine paths. Existing simulation
JSON and CAD/SVG exports remain unchanged; publication refuses overwrites.

This is a **geometry-description foundation only**, not completion of staged
milestone 1 or its acceptance criteria. It exports no collision/inertial model,
grasp/liquid/insertion pose, simulator world, material transfer, robot action,
protocol execution plan, reviewed CUT&RUN instructions or physical approval.
No synthetic recipe is promoted into a reviewed laboratory protocol. The
existing synthetic cloning dependency diagram cannot supply those missing
material/protocol semantics. Mass/inertia remain unknown and actions empty.

Here "front" means the row-A / label-strip edge, not the side facing an
operator; it does not define a calibrated ANSI/SLAS deck pose. D1 is an owner
decision, not a request to rotate or relabel existing exports in this slice.

The [one-family emitter pilot](#test-only-emitter-pilot-2026-11-04) now provides
only its bounded source-level regressions. Further CAD or simulator work needs
its own proposal. Reviewed, quantity/sample-bound protocol input and catalog
ownership/licensing remain prerequisites for the later protocol/world stage,
not invented inputs to the CAD-source check. A shareable CAD representation
must reconcile with existing source-backed labware rather than duplicate it.
Physical execution requires its own later approval.

See [CLI usage](cli/rendering.md), [wire semantics](protocol.md)
and the [review prompt](robotic_wetlab_equipment_review_prompt_20261103.md).

## Verification Ledger

Authoring is on macOS arm64, based on the preserved Glen tip `4a633418`.
All Cargo commands use the unchanged owner lockfile with `--locked --offline`,
one build job and `CARGO_INCREMENTAL=0`. Transient development/test debug level
is zero; repository profiles, dependencies and the lockfile are not edited.

The first check exposed an incorrect extracted-module path, corrected before
runtime testing. The first eight-test run passed six and found two issues:
the PCR96 fixture must respect the existing 8x12 requirement, and direct f32
text serialization must match the JSON-value projection. The fixture and new
description serialization were corrected; old carrier admission/geometry and
legacy serialization were not changed. The initial discovery script also
correctly flagged the not-yet-regenerated parity matrix; it is not a runtime
failure or an accepted stale projection.

Final executed checks, with Rust `1.100.0-beta.7 (fd8e211ae 2026-10-17)`:

```sh
# All Cargo invocations also use these transient environment settings:
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
cargo test -q --locked --offline --lib -- rack_equipment_ export_rack_ parse_racks_ --test-threads=1
cargo test -q --locked --offline -p gentle-protocol --lib
cargo test -q --locked --offline --test capability_registry_parity --test mcp_capability_surface --test parity_matrix_freshness
cargo run -q --locked --offline --bin gentle_examples_docs -- parity-matrix-generate
cargo run -q --locked --offline --bin gentle_examples_docs -- parity-matrix-check
cargo check -q --locked --offline
python3 scripts/generate_capability_maturity.py --write
python3 scripts/generate_capability_maturity.py --check
python3 -m unittest scripts.test_capability_maturity scripts.test_tutorial_checkouts -q
python3 scripts/manual_docs.py
python3 scripts/check_documentation_status.py
git diff --check
```

Results: 29 focused rack/equipment/parser tests (all eight new tests plus 21
existing regressions), 120 protocol tests, 26 capability-registry tests,
six MCP discovery/error-contract tests, one parity-freshness test and 98
Python maturity/checkout tests pass, with no ignored/skipped tests in those
runs. Canonical matrix generate/check, the four required glossary-hash-only
maturity updates, default-feature Cargo check, manual/status and whitespace
checks pass. The LF/CRLF coverage is synthetic checkout/byte-preservation
coverage, not a native Windows execution receipt. The linker reports its
large existing `__eh_frame` compact-unwind warning; no performance claim follows.
The unchanged lockfile SHA-256 is
`7ad5b2dc9485675a75d272fa945c75342dceb60ef582f2885015f3e14df41e76`.

A deliberately broader 30-test run (the 29 above plus
`execute_introspect_all_capability_section_has_full_fact_annotations`) passed
29 and failed that global guard only: `primers_preflight` is `registry_only`.
At baseline `4a633418`, `src/mcp_server.rs` already contains that MCP tool,
while `src/engine_shell.rs` annotates only `primers preflight` (space form).
Those preflight records are untouched here; this is a separate inherited
annotation gap, not a reason to weaken the guard or expand this robotics diff.
No separate baseline rebuild or full-suite pass is claimed. Both new equipment
descriptors are explicitly tested as `fact_annotated` and confirmation-required.

None of this work establishes native Windows/Linux publication semantics, physical fit, simulator
or CAD-tool conformance, private-world acceptance, scientific success, a live
agent/GUI verdict, standalone CLI/MCP export replay, optional JS/Lua VM execution,
or full-workspace-suite acceptance. The unrelated untracked
`output/` directory is untouched. No push, CI dispatch, install, large download,
world export or hardware action is performed.

## Verification Follow-Up, 2026-11-03

The owner supplied Claude's concrete review of `602d7fe6` and authorized the
same locked/offline, single-job, zero-debug Rust verification recipe above.
The expected baseline was `c6adcb83`; actual local main had advanced once to
`89df315a` (dependency requirements pinned to already locked versions). This
worktree was fast-forwarded to that main before editing. Since `602d7fe6`,
the cited shell file also gained an unrelated protein-expression module include;
the equipment implementation/geometry remained unchanged. Both
`git branch -r --contains 602d7fe6` and `git tag --contains 602d7fe6` were empty,
permitting the optional, absent-by-default v1 pitch-input field. No stale
REBASE_HEAD was treated as an active rebase.

Publication changed during implementation: a repeated remote-ref check showed
`origin/main` now containing `602d7fe6`. Work stopped before any commit or main
update, and Steffen explicitly approved: "Keep v1 and explicitly approve this
additive optional field." This is a recorded compatibility exception: old v1
readers reject a populated new field, while absent-override records keep the
same bytes. Both schema/projection constants remain v1 by that owner decision.

The implemented follow-up is byte/value-level validation plus separate local
catalog comparison, not new geometry. Unknown embedded keys/lossy f32 values
are retained long enough to reject, exact input bytes receive their own hash,
and diagnostics never include parser text. Fully rehashed edited sources may
be consistent but cannot claim a catalog match. Old override records without
the recorded override fail the catalog comparison safely. Validation is shared
by shell/agent/API and generic confirmed MCP `op`, with no project-file save.

Equipment export now uses the generalized sequence-design JSON destination
and writer; the sequence-design wrapper keeps its original error namespace.
Existing staging-fault tests remain unchanged. A closed staging file is
published without clobber; bound project/sidecar, raw-path and native checks
are shared rather than copied. D1 and D2 above remain owner decisions. V1
non-null mass/inertia remain invalid; no component order or frame is invented.
The typed compatibility verifier also rejects non-null fields before JSON
serialization, so `Some(NaN)` cannot silently become an accepted null mass.
The existing tamper test now covers this typed-API serialization edge case.

### Executed Verification

Source: `89df315a` plus this scoped change on macOS arm64 (`Darwin arm64`),
Rust `1.100.0-beta.7 (fd8e211ae 2026-10-17)`. Cargo.lock remains unchanged at
`7ad5b2dc9485675a75d272fa945c75342dceb60ef582f2885015f3e14df41e76`.
No dependency, repository profile or feature policy was changed.

```sh
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
cargo test -q --locked --offline --lib -- rack_equipment_ export_rack_ parse_racks_ lifecycle_slice4_ --test-threads=1
cargo run -q --locked --offline --bin gentle_examples_docs -- parity-matrix-generate
cargo run -q --locked --offline --bin gentle_examples_docs -- parity-matrix-check
cargo test -q --locked --offline -p gentle-protocol --lib
cargo test -q --locked --offline --test capability_registry_parity --test mcp_capability_surface --test parity_matrix_freshness
cargo test -q --locked --offline --lib -- rack_equipment_mcp_export_preserves_state_and_matches_saved_json rack_equipment_mcp_validate_preserves_existing_and_missing_state_files --test-threads=1
cargo check -q --locked --offline
python3 scripts/generate_capability_maturity.py --write
python3 scripts/generate_capability_maturity.py --check
python3 -m unittest scripts.test_capability_maturity scripts.test_tutorial_checkouts -q
python3 scripts/manual_docs.py
python3 scripts/check_documentation_status.py
rustfmt --edition 2024 --config skip_children=true --check crates/gentle-protocol/src/rack_equipment.rs src/engine/state/lineage_containers/equipment.rs src/engine/state/lineage_containers/equipment/tests.rs src/engine/io/sequence_design_preview.rs
git diff --check
scripts/maintenance_chore.py session-close --plan docs/robotic_wetlab_simulation_plan.md
```

Final results: 58 focused tests pass, including both MCP tests, the unchanged
applicable `lifecycle_slice4_` tests, seven-template hash pins and orientation
characterization. The two named MCP tests also pass separately (overlapping
coverage, not two additional unique tests). Protocol tests pass 120/120;
capability/MCP/freshness suites pass 27/27, 6/6 and 1/1; Python maturity and
LF/CRLF checkout tests pass 98/98. None are ignored or skipped in those runs.
Canonical matrix generation/check, glossary-hash-only maturity projections,
Cargo check, the four-file scoped rustfmt check, manuals/status and whitespace
checks pass. No maturity verdict was promoted. Tutorial inspection found only
the unchanged 06.07 sequence-design wrapper path, with no equipment tutorial
claims to revise. All seven hash triples are established as code literals
using `producer_version="pinned"`, not Python-generated expectations.

Authoring diagnostics are retained separately: the first 58-test run passed
56, exposed the deliberately unpopulated hash-pin constants, and was denied
permission for the existing Unix-socket fixture under the sandbox. After
recording the emitted hash constants, the exact focused command was rerun
with approved unsandboxed permissions and passed all 58. No sequence-design
test body or assertion was changed. An intermediate 33-test Python subset
correctly rejected the not-yet-regenerated parity matrix; the final full
98-test run above passes after canonical regeneration. Linker compact-unwind
warnings remain non-fatal and are not performance evidence.

The separate command below was also executed and **failed**, unchanged, with
only `("primers_preflight", "registry_only")`; zero tests passed in that run.
It is not included in any pass count and was not weakened or fixed here:

```sh
cargo test -q --locked --offline --lib execute_introspect_all_capability_section_has_full_fact_annotations -- --test-threads=1
```

### Changed Files

- Contract/discovery: `crates/gentle-protocol/src/rack_equipment.rs`,
  `crates/gentle-protocol/src/lib.rs`, `docs/glossary.json`,
  `docs/gui_cli_mcp_parity.md`, `docs/agent_interface.md`.
- Shared implementation and tests:
  `src/engine/state/lineage_containers/equipment.rs`,
  `src/engine/state/lineage_containers/equipment/tests.rs`,
  `src/engine/io/sequence_design_preview.rs`, `src/engine/operation_lifecycle.rs`,
  `src/engine/ops/operation_handlers.rs`, `src/engine.rs`,
  `src/engine/protocol.rs`, `src/engine_shell.rs`, `src/mcp_server.rs`.
- Result-literal/description plumbing only: `src/app.rs`, `src/app/tests.rs`,
  `src/main_area_dna/tests.rs`, `src/engine/analysis/rna_reads.rs`.
- Documentation: `docs/protocol.md`, `docs/cli/rendering.md`,
  `docs/robotic_wetlab_simulation_plan.md`, `docs/CHANGELOG.md`, `docs/roadmap.md`.
- Required glossary-hash-only rollups: `assets/capability_maturity.json`,
  `integrations/clawbio/skills/gentle-cloning/CAPABILITY_MATURITY.json`,
  `integrations/clawbio/skills/gentle-pcr-primer-design/CAPABILITY_MATURITY.json`,
  `integrations/clawbio/skills/gentle-tss-collection/CAPABILITY_MATURITY.json`.

### Remaining Boundaries

Native Windows
publication, JS/Lua VM execution (not run by CI), full-suite/native GUI and
physical/scientific acceptance remain unverified. The inherited
`primers_preflight` global fact-annotation failure is neither fixed nor counted
as a pass. No push, tag, CI dispatch, installation or hardware action is allowed.
The paper submodule, untracked prompt notes, data and output are excluded from
staging/committing.

Delivery decision: after overlapping uncommitted work appeared on local main,
Steffen requested committing only this follow-up on `gentle_rs_2_main`, leaving
main integration pending. No main update, merge, rebase or push is performed.
Session-close reports four OK, two warnings and zero failures: the expected
working-tree/generated-artifact warning and the manual plan-fidelity reminder.
The unrelated untracked `output/` directory remains untouched.

### Status After Integration, 2026-11-04

The delivery note above is historical. `10c3545b` has since been integrated
into `main` and reached `origin/main`. The later `069bd37b` annotates the MCP
`primers_preflight` descriptor, so the inherited global fact-annotation failure
recorded above is no longer open.

Claude's review of `602d7fe6` is closed: `10c3545b` addresses every finding.
At `069bd37b` on macOS arm64 (`rustc 1.101.0-nightly (3e92e1a2d 2026-11-01)`),
with the unchanged lockfile and the transient settings above, this command
passed 59/59 tests (the 58 focused tests plus the global guard):

```sh
cargo test -q --locked --offline --lib -- rack_equipment_ export_rack_ parse_racks_ lifecycle_slice4_ execute_introspect_all_capability_section_has_full_fact_annotations --test-threads=1
```

The protocol, parity and Python suites were not rerun for this note. D1, D2,
native Windows publication, JS/Lua VM execution and physical/scientific
acceptance remain open.
