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
