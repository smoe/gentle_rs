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
- **D1 - Orientation:** retain and explicitly name the current row-A / label-strip
  edge as GENtle's "front" (not the operator-facing side), or adopt ANSI/SLAS
  plan-view orientation in a separate legacy-export migration. See the
  [current-output characterization](protocol.md); frames and labels are unchanged here.
- **D2 - Published file mode:** tempfile-based no-clobber exporters publish 0600
  on Unix, unlike ordinary `fs::write` exports. Set repository-wide publication
  policy separately; do not change one exporter in isolation.

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
