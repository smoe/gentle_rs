# AGENTS.md

Repository-local instructions for coding agents working in this project.

## Scope

- This file applies to the whole repository rooted at this directory.
- Follow these rules in addition to platform/tool-level safety policies.

## Read First

Before non-trivial edits, read:

1. [`docs/architecture.md`](docs/architecture.md)
2. [`docs/roadmap.md`](docs/roadmap.md)
3. Relevant interface docs for your change:
   - GUI: [`docs/gui.md`](docs/gui.md)
   - CLI: [`docs/cli.md`](docs/cli.md)
   - Protocol: [`docs/protocol.md`](docs/protocol.md)

## Core Principles

- Keep one deterministic engine contract across GUI/CLI/JS/Lua.
- Do not duplicate business/biology logic in frontend adapters.
- Prefer extending shared parser/executor paths over adding one-off command paths.
- Keep behavior deterministic and machine-readable (structured errors/results).

Reference: [`docs/architecture.md`](docs/architecture.md)

## Planning and Status Discipline

- Treat [`docs/roadmap.md`](docs/roadmap.md) as the source of truth for:
  - known gaps,
  - priority tracks,
  - execution order.
- If you materially change project status or close/open a gap, update
  `docs/roadmap.md` in the same change.

## Testing Expectations

- Run targeted tests for changed code paths first.
- Run `cargo check -q` before handoff.
- If tests are missing for the changed behavior, add them or document the gap in
  `docs/roadmap.md`.
- Do not claim behavior is complete without at least one deterministic test path
  (unit or integration) that exercises it.

## Test Data Provenance (Mandatory)

- Every committed fixture under `test_files/` or `tests/` must document:
  - origin (exact source or explicit synthetic/hand-crafted status),
  - deterministic recreation/retrieval steps,
  - where it is used in GENtle (runtime/parser/tests).
- Add or update a nearby fixture README/manifest as part of the same change.

Reference: the test-data provenance rule in
[`docs/architecture.md`](docs/architecture.md).

## GUI and Agent Safety

- Keep screenshot capture paths disabled unless explicitly approved by project
  policy.
- Preserve agent guardrails:
  - no recursive/nested `agents ask` execution from agent-suggested commands,
  - avoid creating command routes that can silently bypass that rule.
- For UI-intent features, use shared `ui ...` shell command contracts and keep
  GUI handlers thin.

References:

- [`docs/architecture.md`](docs/architecture.md)
- [`docs/roadmap.md`](docs/roadmap.md)

## Documentation and Handoff

- Update docs when user-visible behavior changes.
- In handoff notes, list:
  - files changed,
  - tests run,
  - any remaining risks or missing tests.

## Source Documentation Plan

- Goal: make `cargo doc` output useful at module/API level for every file under
  `src/`.
- Baseline requirement:
  - each module/file should have a top-level `//!` summary,
  - each `pub mod` in `src/lib.rs` should carry a short `///` description,
  - key public structs/enums/functions should explain purpose + invariants.
- Pilot module for style and scope:
  - `src/genomes.rs` should document:
    - catalog/source resolution model,
    - preparation/indexing lifecycle and manifest side effects,
    - annotation parsing behavior (formats + malformed-line policy),
    - extraction/BLAST contracts and failure modes.
- Rollout sequence:
  1. Complete module-level docs in infrastructure-heavy modules
     (`genomes`, `engine`, `engine_shell`).
  2. Cover rendering and UI orchestration modules.
  3. Cover adapter modules (`cli/js/lua`) and remaining utilities.
  4. Keep docs updated in the same change whenever public behavior changes.
