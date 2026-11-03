# Claude Review: Offline Rack Equipment Foundation

Please review the GENtle changes in
`/Users/u005069/.codex/worktrees/47dd/gentle_rs` read-only. Do not implement,
commit, fetch, launch simulators/robots, install tools, or contact another agent.
Locate the implementation commit with subject
`feat(labware): export source-bound rack equipment descriptions`, whose parent
is `4a633418`. Review that exact commit against its parent, not an arbitrary
later HEAD or unrelated changes already integrated from main.

## Context And Original Plan

Steffen requested merging Glen's forward-ported robotics/wet-lab roadmap and
implementing the reasonably bounded next slice. Glen's bundle SHA-256 is
`5083c5f1131228236468fc1f35e46202adba5741486e80495dd31f2e41caeaa8`;
its exact base is `888e49ce` and its three preserved commits are `8d91d696`,
`5ae0b514`, `4a633418`. This plan is post-0.1 and grants no physical authority.

Codex's implementation plan: reuse the existing source-backed rack/template
geometry, leave historical simulation JSON/CAD/SVG unchanged, and add one
engine-owned typed equipment-description projection. Bind canonical source
and output content, make local frames/units explicit, preserve unknowns, use
existing shared-shell/operation routes and approval, publish without overwrites,
and test geometry/parity/refusal deterministically. Do not invent an equipment
catalog, universal CAD language, reviewed protocol, Gazebo adapter, physics or
robot commands before ownership/licensing and reviewed material inputs exist.

## Review Targets

- `docs/robotic_wetlab_simulation_plan.md`, especially First Implementation Slice.
- `crates/gentle-protocol/src/rack_equipment.rs` and the registry additions.
- `src/engine/state/lineage_containers/equipment.rs` and `equipment/tests.rs`.
- `src/engine.rs`, `src/engine/operation_lifecycle.rs`, operation handler and
  optional `OpResult.rack_equipment` wiring (other result literals only add None).
- `src/engine_shell.rs`, `docs/glossary.json`, introspection descriptors,
  `docs/agent_interface.md` and the generated parity/maturity projections.
- `docs/protocol.md`, `docs/cli/rendering.md`, changelog and verification ledger
  in the plan. Check exact executed-source/test limitations there; source review
  does not turn unrun native/physical/full-suite checks into acceptance.

## Questions

Known separate check failure: the expanded global introspection guard finds
`primers_preflight` as `registry_only`. Baseline `4a633418` already advertises
that MCP name but annotates only the space-form `primers preflight`. The 29
equipment/rack/parser regressions, 120 protocol tests, 26 + 6 + 1 parity/MCP
tests, 98 Python tests and Cargo check pass. See the ledger for exact scope;
do not count the global guard as passed or fold an unrelated preflight repair
into this equipment review. Confirm the source-based baseline attribution if
needed; a separate baseline Rust build was not performed.

1. Do source, asset and content hashes bind the documented inputs independently
   of JSON map ordering? Can modified derived geometry/actions pass verification
   merely by recomputing the top-level hash? Distinguish internal consistency
   from authentic or scientifically/physically validated inputs.
2. Are local right-handed metre frames, embedded millimetres, slot order,
   blocking, body-envelope limits and quaternion order unambiguous? Verify
   existing geometry helpers are reused, not duplicated. Is an absent physics
   model clearly absent rather than a guessed default?
3. Are legacy formats and project content/undo preserved? Do actual shell/API
   results contain the same record as the file, including embedded f32 values
   through both typed JSON text and JSON-value serialization? Are source drawings/unknowns
   retained without leaking occupancy, sequences or private room geometry?
4. Are no-clobber publication/cleanup and file-handle lifetimes sound on Windows
   as well as Unix? Do not treat macOS execution as native Windows proof.
5. Are CLI, generic MCP, agent and optional scripting paths genuinely shared,
   with request-derived approval that forged planner flags cannot bypass?
   Is discovery/introspection/parity fresh without promoting maturity?
6. Is this a useful minimal foundation, rather than prematurely claiming the
   full milestone? Identify concrete blockers before any simulator/world slice.

Return concise findings by severity with file/line evidence, missing-test
boundaries and a revised minimal next plan. No speculative robotics redesign.
