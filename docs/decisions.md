# GENtle Decisions

This file records durable implementation constraints and architectural
decisions. It is intentionally short: if an entry becomes a backlog item, move
that work to [`roadmap.md`](roadmap.md); if it describes completed history,
move the outcome to [`CHANGELOG.md`](CHANGELOG.md).

Status values:

- `active`: currently governs implementation.
- `superseded`: kept for history; no longer governs new work.
- `deferred`: accepted direction, but not currently active work.

## DEC-001: Single Shared Engine

Status: active

GENtle uses one deterministic engine contract across GUI, CLI, shared shell,
JS, Lua, Python wrapper, MCP, and agent-facing routes. Adapters must not fork
biology or business logic.

## DEC-002: Protocol-First Machine Contracts

Status: active

Machine-facing routes use structured, versioned request/result records. Durable
schemas and wire contracts live in [`protocol.md`](protocol.md), not in the
roadmap.

## DEC-003: Thin Adapter Boundary

Status: active

Frontend and integration adapters should parse, route, present, and collect
artifacts. Core decisions, reports, findings, and explanation records should be
engine-owned and portable whenever practical.

## DEC-004: Roadmap/Changelog/Decision Split

Status: active

`roadmap.md` is for next work only. Completed work goes to
[`CHANGELOG.md`](CHANGELOG.md). Durable constraints go to this file. A session
should be able to start by reading only `roadmap.md` in under two minutes.

## DEC-005: No Roadmap `Done` Entries

Status: active

Do not add a `Done` entry to `roadmap.md`. Move completed work directly to
`CHANGELOG.md` in the same session. This prevents the roadmap from growing back
into a mixed history/planning document.

## DEC-006: Workspace Extraction Order

Status: active

The intended crate extraction order remains:
`gentle-protocol -> gentle-engine -> gentle-render -> gentle-shell -> gentle-gui`.
First-wave extraction avoids per-feature micro-crates; related analysis logic
should remain together until the engine boundary is stable.

## DEC-007: Root-Crate Compatibility During Extraction

Status: active

While production code still partly lives in the root crate, extracted crates
must expose stable contracts without forcing all callers through a one-shot
import rewrite. Root modules may re-export extracted contracts during the
transition.

## DEC-008: Shell Dispatch Split Pattern

Status: active

Shared shell parsing and execution should remain one behavior path for GUI
Shell, `gentle_cli shell`, MCP, JS/Lua wrapper helpers, and agent execution.
Large command families may split into helper dispatch functions, but those
helpers must preserve shared parser/executor semantics.

## DEC-009: Stack-Hardening Policy

Status: active

Stack-overflow fixes should reduce dispatcher frame depth without forking
behavior. Prefer narrow helper dispatch and expanded-stack workers for
confirmed small-stack failures in shell/engine routes.

## DEC-010: `#[inline(never)]` Helper Dispatch

Status: active

Until the monolithic root dispatchers are fully decomposed, stack-sensitive
command or operation families may dispatch through dedicated `#[inline(never)]`
helpers before entering broader match frames. Inner match branches should
delegate back to the same helper so there remains one implementation path.

## DEC-011: UI-Intent Shared Catalog

Status: active

UI-intent targets and discoverability metadata are shared catalog data.
Menus, command palette, shell/agent routes, MCP discovery, and ClawBio/OpenClaw
handoff should consume the same target catalog instead of hard-coding parallel
target lists.

## DEC-012: GUI Intent Handlers Stay Thin

Status: active

GUI `ui open|focus|close|selection ...` handlers may open/focus/close host
windows, dialogs, and active viewer selections, but must route target
resolution through shared `ui ...` shell command contracts where possible.

## DEC-013: Mutating Agent/MCP Safety

Status: active

Mutating `op` and `workflow` execution through agent/MCP paths requires
explicit confirmation. Agent-suggested commands must not recursively invoke
`agents ask`, `agents plan`, or `agents execute-plan`, including through macro
expansion paths.

## DEC-014: Screenshot Capture Policy

Status: active

Screenshot capture paths stay disabled unless project policy explicitly
approves them. Documentation should prefer deterministic SVG/export routes
while screenshot execution is policy-disabled.

## DEC-015: ClawBio/OpenClaw Boundary

Status: active

The ClawBio/OpenClaw skill scaffold wraps deterministic `gentle_cli` routes and
writes reproducibility bundles. It should surface shared capabilities,
suggested actions, and UI-intent handoffs without creating ClawBio-only biology
logic.

## DEC-016: External AI/Automation Deployment

Status: active

For tool-driven external AI deployment, prefer the published GHCR image in
headless `mcp` mode with explicit project/state mounts and stdio communication
instead of depending on browser GUI containers.

## DEC-017: Fixture Provenance

Status: active

Every committed fixture under `test_files/` or `tests/` must document origin,
deterministic recreation/retrieval steps, and where GENtle uses it.

## DEC-018: Source Documentation Policy

Status: active

Public modules and key public records should carry concise rustdoc explaining
purpose, invariants, and non-obvious behavior. Prefer why/invariant/edge-case
documentation over line-by-line narration.

## DEC-019: Primer/Oligo Material Identity

Status: active

A planned primer/probe or in-silico design artifact is not physical stock.
GUI, CLI, JS, Lua, Python, MCP, and agent routes must preserve the distinction
between design, reviewed order, received material, local availability, and
consumption.

## DEC-020: Inline/Stateless Sequence Inspection

Status: active

Read-only operations that need only sequence letters, optional topology, and an
optional span should be state-optional. Promoting such results into project
state remains an explicit second step.

## DEC-021: Helper-Construct Terminology Migration

Status: deferred

The legacy "helper genome" vocabulary should eventually move toward
"helper construct" semantics. The migration should be atomic across contracts,
docs, GUI, and examples to avoid mixed terminology.

## DEC-022: Primer3 Backend Parity

Status: deferred

Primer3 is the intended external backend parity layer for primer design, but
deeper constraint mapping, fixture coverage, and backend equivalence checks
should land through shared engine/report contracts rather than GUI-only code.

## DEC-023: Presentation i18n Does Not Localise Machine Contracts

Status: active

GUI labels and dialog text may be translated at runtime, but shared shell
commands, protocol schema fields, saved project records, adapter payloads,
agent response schemas, and scientific identifiers remain deterministic
English. Localised text is presentation only.
