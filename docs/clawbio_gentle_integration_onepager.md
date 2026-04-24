# ClawBio × GENtle Integration Sketch (Technical One-Pager Attachment)

## Objective

Define a deterministic integration boundary where ClawBio orchestrates strategy and GENtle executes biology operations through stable contracts.

## Why This Is Different From Standalone Cloning Tools

Desktop cloning tools (including prior GENtle-style usage) already provide local sequence editing/design.
The ClawBio x GENtle integration adds a separate orchestration and policy layer:

- policy-aware intent handling before execution
- explicit non-mutating-first planning loop
- deterministic run gating with role-based control
- machine-readable provenance bundles for each orchestrated run

## Role Split

- ClawBio: intent interpretation, strategy selection, constraint policy, fallback policy.
- GENtle: deterministic execution (`shell`, `op`, `workflow`), schema-tagged outputs, reproducibility artifacts.
- Shared middle ground: GENtle may compile prose into typed plan candidates,
  but ClawBio still owns whether those candidates should be executed.

## Governance and Control Model

Control is intentionally split into three planes:

- Intent plane (ClawBio):
  - Proposes strategies and candidate command/workflow payloads.
  - Does not get implicit execute authority.
- Execution plane (GENtle):
  - Executes only explicit, schema-valid commands.
  - Keeps deterministic outputs and run artifacts.
- Policy plane (project/lab policy owner):
  - Enforces allow/deny scope and approval requirements before mutating runs.

Minimum control baseline for integration:

- explicit non-use scope for human/animal/plant genome-editing workflows
- default advisory/read-only orchestration mode
- explicit approval gate before mutating execution
- local-first data boundary and explicit export only
- immutable run-level audit artifacts and checksums

## Runtime Boundary

- Headless runtime: `gentle_cli` (GUI-independent).
- ClawBio wrapper scaffold: `integrations/clawbio/skills/gentle-cloning/`.
- Wrapper request schema: `gentle.clawbio_skill_request.v1`.
- Wrapper result schema: `gentle.clawbio_skill_result.v1`.

## Request/Result Contract

Request (`gentle.clawbio_skill_request.v1`):
- `mode`: `capabilities | state-summary | shell | op | workflow | agent-plan | agent-execute-plan | raw`
- Optional controls: `state_path`, `timeout_secs`
- Mode payload:
  - `shell`: `shell_line`
  - `op`: `operation`
  - `workflow`: `workflow` or `workflow_path`
  - `agent-plan`: `system_id`, `prompt`, optional planner/runtime overrides
  - `agent-execute-plan`: `plan` or `plan_path`, `candidate_id`, optional `confirm`
  - `raw`: `raw_args[]`

Result (`gentle.clawbio_skill_result.v1`):
- `status`: `ok | command_failed | timeout | failed | degraded_demo`
- Includes executed command, exit code, stdout/stderr, resolver metadata, and artifact paths.
- `stdout_json` is populated when wrapped `gentle_cli` stdout parses as JSON.
- `chat_summary_lines[]` is populated for `gentle.sequence_context_view.v1`
  results so chat layers can answer with the compact DNA-window summary first.

## Deterministic Execution Pattern

1. Inspect (non-mutating): `capabilities`, `state-summary`, `features query`, `primers preflight`.
2. Either plan in ClawBio directly or ask GENtle to compile prose with `agent-plan`.
3. Execute one stored candidate explicitly with `agent-execute-plan`, or run direct `op` / `workflow`.
4. Collect artifacts for audit/replay.

## Artifacts and Provenance

Expected reproducibility bundle per run:
- `report.md`
- `result.json`
- `reproducibility/commands.sh`
- `reproducibility/environment.yml`
- `reproducibility/checksums.sha256`

These artifacts are the handoff unit for ranking, review, and replay.

## Safety and Guardrails

- No implicit mutating execution: ClawBio must choose mutating modes explicitly.
- Keep data local by default; export is explicit.
- Maintain recursion guardrails for agent-invoked commands.
- ClawBio should not treat `agents ask` as its primary integration boundary;
  the preferred shared AI-facing boundary is the typed planner plus the
  read-only transport/preflight metadata.
- Prefer schema-tagged JSON outputs over prose when ClawBio consumes results.
- Treat policy denials as first-class deterministic outcomes (not hidden failures).

## Minimal Wire Example

```json
{
  "schema": "gentle.clawbio_skill_request.v1",
  "mode": "shell",
  "shell_line": "features query tp73 --kind CDS --within --range 61784..63000 --label TP73 --limit 200"
}
```

```json
{
  "schema": "gentle.clawbio_skill_result.v1",
  "status": "ok",
  "chat_summary_lines": [
    "1 CDS feature in view"
  ],
  "stdout_json": {
    "schema": "gentle.sequence_context_view.v1"
  }
}
```

## MVP Integration Scope

- Orchestrated feature/region inspection
- Primer-design planning and execution via deterministic workflows
- Structured report export for downstream strategy ranking in ClawBio
