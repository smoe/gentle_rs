# Claude Review: Equipment Emitter Pilot Follow-Up

Review the revised plan read-only in
`/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
The verified starting HEAD is `60bae970`; review the documentation changes
against that source, whether committed or in the working tree, and report any
drift before relying on it. If your read-only tools cannot inspect Git state,
say so and review the named files without widening permissions. Do not
implement, commit, fetch, build/test Rust, install/download tools, run CAD or
simulators, take screenshots, contact another agent or launch nested Codex.

## Task And Original Codex Review

Steffen asked to revise Claude's robotics roadmap, narrow the first slice to
one tube/rack family and one exporter, and prepare this review prompt. The
previous attachment matched the tracked plan exactly. This is not a request
to reopen the completed rack-equipment implementation/review or add a release
gate to the current internal series.

The original broad plan required a canonical parametric representation,
SolidCraft and FreeCAD generators and a neutral execution plan in Stage 1;
its first useful milestone additionally required a public world and CUT&RUN
simulation. Codex's read-only review raised three P2 concerns:

1. The first gate combined several independently useful stages.
2. Determinism did not distinguish exact source/binding reproduction from
   numerical replay with declared tolerances and runtime configuration.
3. Interrupted-run receipts needed stable step/attempt IDs, evidence
   completeness, idempotent import/retry and explicit partial-state acceptance.

Compatibility guardrail: shipped equipment v1 is geometry-only. The validator
rejects non-null mass/inertia; richer physics/interaction semantics need a
separately reviewed schema/version or companion record, not filled placeholders.

## Supplied Claude Review And Codex Response

The owner relayed your source-only review at `60bae970`; no runtime evidence
was claimed. This revision incorporates its concrete findings:

- P2-1: name the check as emission consistency, not independent derivation
  verification. Add hand-calculated literal oracles alongside parsed-spec
  comparisons, without using the production helper to compute expectations.
- P2-2: say tube and rack checks are independent, with no pairing, placement
  or fit invariant. The selected reference tube does not alter the carrier.
- P2-3: make the pilot unblocked. Document that SCAD omits blocked holes while
  equipment v1 retains blocked slots; do not force equivalent behavior.
- P3: fix the rack identity, numeric tolerances and incremental test scope.
  Existing tube tests already pin dimensions, bytes and profile JSON. Existing
  carrier tests pin 4x9 JSON geometry but not all 4x6 SCAD numerical emissions.
- Defer the proposed engine-owned conformance record. With no runtime consumer,
  the useful first slice is test-only, with temporary generated artifacts and
  no new schema, report writer, ledger or parser exposed to users.

One numerical clarification: parsed source versus an already widened f32 spec
needs the three-decimal rounding allowance, not a second f32 allowance. The
separate spec-versus-hand-literal comparison retains the existing `1e-4` mm
carrier-test tolerance. Please check these as distinct assertions.

## Revised Codex Proposal

- First gate: test-only emission regression, not geometry-engine redesign,
  independent physical derivation proof or simulator execution. One existing
  `safe_lock_tube_rack_floor_seated` family, a synthetic sample-free 4x6 rack,
  fixed `rack_id="conformance_safe_lock_floor_4x6"`, no compact pitch, zero
  blocked coordinates, no placements and the independent
  `eppendorf_safe_lock_1_5ml` reference visualization.
- Keep both 1.5 and 2.0 ml profile inputs to the existing carrier derivation.
  Selecting one rendered reference tube must not silently turn that shared
  carrier into a different geometry family.
- One backend: existing OpenSCAD source through `ExportTubeOpenScad` and
  `racks openscad`. Two model files remain one backend. Preserve the existing
  outer-body/collar limitations; no cap/cavity/grasp/liquid/insertion pose.
- Parse the fixed generated rack declarations and 24 hole cuts; compare to
  the spec and to the plan's separately hand-calculated literal oracle.
  `SCAD_EMISSION_ABS_TOL_MM=0.000500001` and
  `SPEC_LITERAL_ABS_TOL_MM=0.0001` serve separate comparisons. Canonical string
  pins remain exact. Neither tolerance represents manufacturing clearance.
- Reuse existing tube/provenance/repeatability and carrier/legacy tests; add
  only the missing 4x6 rack numerical/coordinate/repeatability assertions.
  Use test diagnostics, no new runtime conformance record or command route.
- Keep D1's limited legacy-orientation choice explicit. Generate only temporary
  outputs; no new committed CAD fixture or publication gate is required for
  those outputs, without waiving existing source-use/provenance obligations.
  Any later committed artifact needs rights/LF/CRLF review. D2 remains
  repository-wide; no writer or legacy mode/overwrite changes.
- Separate source checks from native CAD acceptance, then catalog/backend
  expansion, a minimal synthetic operation/world, CUT&RUN, private worlds and
  separately approved physical rehearsal. No reviewed CUT&RUN input is needed
  merely to compare sample-free CAD source.
- State exact/numerical replay levels and later fail-closed receipt assessment,
  completeness, retries and partial observations explicitly. These do not
  authorize implementing an execution/receipt schema in the pilot.

## Grounding And Files To Read

- `AGENTS.md`, `docs/architecture.md` (rack/profile/export ownership),
  `docs/decisions.md` (Equipment Descriptions Do Not Grant Execution Authority).
- `docs/robotic_wetlab_simulation_plan.md`: revised staged development,
  determinism/replay, receipt requirements and Next Slice. Historical
  implementation/test ledgers must remain historical and unchanged.
- `docs/roadmap.md`: only the robotics parking-lot/later-release pointer changes;
  `docs/CHANGELOG.md`: documentation-only planning entry.
- `crates/gentle-protocol/src/labware.rs`: drawing-backed profile fields,
  unknowns, `derive_tube_carrier_geometry` and profile IDs.
- `src/engine/state/lineage_containers.rs`: `safe_lock_carrier_spec` consumes
  both Safe-Lock profiles; `rack_physical_template_spec` and
  `export_rack_openscad` already own carrier derivation and source export.
  SCAD skips blocked holes; the header includes `rack_id`; scalars and hole
  translations use three decimal places.
- `crates/gentle-render/src/labware.rs`: `render_tube_openscad` is an existing
  limited visualization, not a printable vessel or fit proof.
  Reuse `tube_models_bind_their_source_dimensions_and_disclose_omissions`.
- `src/engine/tests.rs`: compare the planned increment with
  `labware_carriers_preserve_legacy_json_and_export_both_support_modes`,
  `labware_legacy_carrier_exports_omit_drawing_disclosures`,
  `export_rack_fabrication_isometric_svg_and_openscad_include_template_markers`
  and `rack_openscad_all_templates_preserve_bottom_floor_and_open_top`.
- `crates/gentle-protocol/src/rack_equipment.rs`,
  `src/engine/state/lineage_containers/equipment.rs` and `equipment/tests.rs`:
  current v1, validation, pinned hashes and orientation characterization.
- `docs/cli/rendering.md`: existing tube/rack routes and the individual PCR
  profile's undocumented-collar refusal. Do not substitute it just to make the
  pilot appear single-profile.
- `docs/testing.md`: byte-exact source fixtures need scoped LF/CRLF coverage;
  fixtures need origin/recreation/usage provenance. No runtime result is
  established by this planning review.

## Questions

1. Does this really stay within one family and one CAD backend? Is OpenSCAD
   the smallest credible existing backend, with SolidCraft/FreeCAD later,
   rather than another mandatory implementation campaign?
2. Are the shared two-profile derivation and independent tube/rack checks now
   accurately described, including SCAD's blocked-hole behavior and the fixed
   unblocked case? Flag any implied pairing or missing-dimension inference.
3. Are the hand-computed oracle, exact declarations and both numerical allowances
   correct for the current emitters? Does this avoid circular verification and
   double-counting f32 error without claiming native CAD or geometric validity?
4. Do D1/D2, fixture rights, schema boundaries and writer compatibility remain
   explicit decisions rather than silently approved defaults? Is anything
   incorrectly made a prerequisite for this sample-free pilot?
5. Do later receipt rules prevent duplicate lineage application, missing-evidence
   success and false rollback claims, without designing that whole subsystem now?
6. Is the stated test increment genuinely absent from existing coverage, and
   is deferring a runtime result record the right minimal scope? If the checks
   already exist, identify them rather than requesting duplicate tests.

Return concise findings by severity with repository file/line evidence,
remaining owner decisions and a revised minimal plan only if needed. Distinguish
planning suggestions from verified code behavior. Do not treat source review,
old receipts or generated hashes as new-source/native/scientific acceptance.
