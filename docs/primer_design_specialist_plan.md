# Primer Design Specialist Window Plan

Last updated: 2026-03-03

## Current Gap

Engine and shared-shell support for primer-pair design is already available via
`DesignPrimerPairs` plus report inspection/export (`primers design`,
`primers list-reports`, `primers show-report`, `primers export-report`).

The remaining gap is GUI presentation: there is no dedicated, report-centric
primer design workflow window for designing pairs, inspecting alternatives, and
filtering/ranking outcomes before PCR handoff.

## Scope (This Plan)

- Add a dedicated primer specialist GUI window.
- Support primer-pair design execution and report browsing in that window.
- Support alternative-result filtering with persisted filtered views.
- Add a UI-derived single-primer inspection pane (derived from pair reports).
- Add explicit handoff actions into PCR, PCR Advanced, and PCR Mutagenesis
  workflows.
- Preserve UI-intent and shared-shell parity expectations for agent/MCP/CLI
  pathways.

## Interfaces

### Engine

- Add persisted report-view model for filtered primer alternatives
  (metadata-backed, deterministic replay of view criteria).
- Keep `DesignPrimerPairs` as the canonical pair-design operation.
- Keep ranking/scoring logic engine-owned (no GUI-local biology forks).

### Shell

- Extend primer shell contract with filtered-view commands:
  - `primers view-create`
  - `primers view-list`
  - `primers view-show`
  - `primers view-export`
  - `primers view-delete`
- Keep existing report commands unchanged for backward compatibility.

### UI intent

- Add `primer-design` target to `ui open|focus` intent routing.
- Allow optional context prefill (`seq_id`, ROI, report/view id) while
  preserving deterministic argument validation.

## Phases

### Phase 1: specialist window + pair reports

- Add dedicated primer specialist window entry points (menu + command palette).
- Implement pair-design run form and report list/detail rendering.
- Expose backend provenance in GUI result detail.

### Phase 2: persisted filtered views

- Add filter controls for pair alternatives.
- Add save/load/export/delete flows for filtered views.
- Persist views in project metadata with deterministic IDs and payload schema.

### Phase 3: UI-derived single-primer pane

- Derive single-primer candidates from selected pair report/view.
- Add ranking/filter inspection table for forward/reverse primer records.
- Keep this pane explicitly UI-derived (engine-native single-primer op deferred).

### Phase 4: PCR handoff polish + parity tests

- Add explicit actions to push selected pair into PCR/PCR Advanced/PCR
  Mutagenesis controls.
- Add parity tests for shell/ui-intent behavior and deterministic filtered-view
  payload semantics.

## Acceptance Criteria

- Users can run pair design in a dedicated GUI specialist window without using
  strict engine operation forms.
- Users can inspect alternative pair results and persist filtered views.
- Persisted views survive project save/load and can be exported deterministically.
- Selected primer pairs can be forwarded to PCR flows with explicit user action.
- Shared-shell and UI-intent surfaces expose equivalent capabilities for agent
  routes.

## Test Matrix

- Engine unit tests:
  - filtered-view creation/validation/persistence/export/delete
  - deterministic ordering and tie-break behavior under filter/sort constraints
- Shell tests:
  - parse/execute primer `view-*` commands
  - stable error contracts for invalid view payloads/options
- UI-intent tests:
  - `ui open primer-design` and context prefill semantics
  - target-specific option validation failures are deterministic
- GUI behavior tests:
  - report loading, filtering, view persistence lifecycle, PCR handoff actions

## Risks / Follow-ups

- Scope risk: mixing pair-design, single-primer UX, and specificity BLAST fanout
  in one implementation pass can increase regressions; phase boundaries are used
  to keep rollout deterministic.
- Data-model risk: filtered-view persistence schema must remain stable across
  adapter surfaces; schema/version tests are required.
- Follow-up: introduce an engine-native single-primer design contract after this
  UI track stabilizes.
