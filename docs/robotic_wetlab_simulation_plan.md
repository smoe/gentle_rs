# Robotic Wet-Lab Simulation And Equipment Digital Twins

Status: future post-0.1 roadmap track. This document records the intended
ownership and evidence boundaries; it does not authorize physical robot
execution or add a release gate to the current internal series.

## Goal

Make selected GENtle protocols executable as provenance-bound simulations in a
robotics environment such as Gazebo. The first benchmark should be an
automated CUT&RUN workflow in a synthetic reference laboratory. The simulation
must expose material flow, equipment use, timing, robot reach and collision
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
- physical dimensions, scale, coordinate frames, simplified collision and
  inertial geometry, movable joints and safety envelopes;
- semantic interaction points such as deck slots, tube or plate positions,
  lids, doors, grippers, pipette tips, liquid ports and observable states;
- supported actions, preconditions, resource requirements and known fidelity
  limits;
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

The representation must support parametric primitives, reusable subassemblies,
constraints, materials and named coordinate frames as well as the semantic
interaction points above. The first conformance family is tubes and racks:

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

## Staged Development

1. **Catalog, abstract representation and generators.** Define the
   equipment/labware manifest, canonical parametric representation and the
   simulator-neutral protocol execution plan. Demonstrate deterministic
   SolidCraft and FreeCAD generation with a small tube-and-rack family. Resolve
   ownership, licensing and private-local overrides before adding
   simulator-specific operations.
2. **Synthetic reference cell.** Build a small redistributable lab world with
   one robot, bench, racks, plates/tubes, pipetting and representative
   instrument interfaces. Prove deterministic loading, placement and collision
   checks.
3. **CUT&RUN benchmark.** Map a reviewed CUT&RUN protocol to typed operations,
   including consumable and sample lineage, timing and explicit manual or
   unsupported steps. Exercise success, missing-equipment, unreachable-pose,
   collision, depleted-consumable and interrupted-run cases.
4. **Private real-lab world.** Bind the locally reconstructed room and equipment
   inventory without moving its geometry or sensitive metadata into public
   fixtures. Compare planned paths and resource schedules with the synthetic
   benchmark.
5. **Physical rehearsal later.** Only after separate robotics safety review,
   calibrated hardware identities, emergency-stop and human-presence policy,
   execute bounded dry motions or hardware-in-the-loop tests. Physical wet-lab
   execution remains a separately approved capability.

## Acceptance Boundaries

The first useful milestone requires:

- one public equipment catalog and synthetic Gazebo world reproducible from
  exact hashes;
- one canonical tube-and-rack family reproducibly generated for SolidCraft and
  FreeCAD, with exported geometry and semantic invariants checked against the
  same abstract source;
- one CUT&RUN simulation plan that can be previewed, approved, replayed and
  refused deterministically;
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

The next useful slice needs maintainer decisions on catalog ownership/licensing
and an explicitly reviewed, quantity/sample-bound protocol input. Only then
should a neutral plan and synthetic-world adapter be designed. A shareable CAD
representation must reconcile with the existing source-backed labware model,
not duplicate it. Physical execution requires its own later approval.

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
