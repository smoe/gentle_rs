# Maintenance Chore Plan

This plan defines a small set of recurring repository chores that keep GENtle
efficient for humans, coding agents, and external inspectors. These chores do
not replace feature work. They reduce drift, keep evidence discoverable, and
make the project easy to audit without creating a second process treadmill.

Grounding rule: every finding must cite concrete repository evidence such as a
commit SHA, diff, file path, test result, CI signal, generated artifact, or
documented command. If the evidence is weak, record no finding.

## Output Contract

Every chore should produce the same compact handoff shape:

- Scope: files, commits, time window, or release gate inspected.
- Evidence: exact paths, commands, SHAs, tests, or artifacts.
- Findings: only concrete problems, ordered by release or maintenance risk.
- Minimal fix: the smallest safe patch or documentation move.
- No-finding note: what was checked when no issue is found.

Completed chore outcomes go to [`CHANGELOG.md`](CHANGELOG.md) only when they
change project state, documentation, policy, or release status. The roadmap
should receive only changed next-work priorities or one-line parking-lot items.

## Chore 1: Session-Close Hygiene

Cadence: at the end of any non-trivial coding or documentation session.

Purpose: prevent planning entropy from creeping back into the roadmap, confirm
that the implemented work matches the plan that authorized it, and keep the
working tree understandable for the next session.

Inputs:

- [`roadmap.md`](roadmap.md)
- [`CHANGELOG.md`](CHANGELOG.md)
- [`decisions.md`](decisions.md)
- [`../AGENTS.md`](../AGENTS.md)
- Any plan, issue, PR, or user request that authorized the session work.
- `git status --short`

Checks:

```bash
scripts/maintenance_chore.py session-close --plan PATH_OR_LABEL
wc -l docs/roadmap.md
rg -n "^\s*[-*]\s*(Done|shipped|baseline is now available|already implemented)" docs/roadmap.md
rg -n "CHANGELOG.md|decisions.md|Done entry|completed work directly" docs/roadmap.md docs/architecture.md AGENTS.md
rg -n "^(<<<<<<< .+|=======|>>>>>>> .+)$" docs AGENTS.md
git status --short
```

The script is the preferred first pass because it reports the same checks in
the shared output-contract shape. The raw commands remain documented so a human
or external inspector can verify the script's behavior without trusting it.
Its plan-fidelity item is deliberately a reminder, not an automatic pass/fail
check, because only the session owner can confirm whether the final diff stayed
inside the authorized plan.

Pass criteria:

- `docs/roadmap.md` stays between 200 and 400 lines unless explicitly accepted.
- No completed-history bullets remain in `docs/roadmap.md`.
- Completed work from the session is in `docs/CHANGELOG.md`.
- New durable rules are in `docs/decisions.md` or `docs/architecture.md`.
- Implemented changes match the plan, issue, PR, or user request that authorized
  them; any divergence is recorded in `docs/CHANGELOG.md` or handoff notes.
- Handoff notes distinguish intentional edits from unrelated dirty files,
  untracked generated artifacts, and dirty submodules.
- Generated artifacts referenced by docs are reproducible from documented
  commands or explicitly called out as release attachments/local evidence.

## Chore 2: Daily Bug Scan

Cadence: daily, or after a merge-heavy development day.

Purpose: catch likely regressions introduced by recent commits without
inventing issues.

Inputs:

- Commits since the last scan or last 24 hours.
- Diffs, touched tests, CI signals, and relevant docs.
- Optional Codex automation memory at
  `$CODEX_HOME/automations/daily-bug-scan/memory.md`, which records the last
  automation run boundary for the coding agent and is not required for human
  review.

Pass criteria:

- Findings cite commit SHAs, file paths, diffs, failing tests, or CI output.
- Proposed fixes are minimal and avoid unrelated cleanup.
- Weak suspicions are skipped rather than inflated into bugs.

## Chore 3: Weekly Drift, Test, And Provenance Scan

Cadence: weekly, after large feature merges, and whenever new fixtures or
parser inputs are added.

Purpose: identify maintainability drift and test/provenance gaps before they
become release blockers. This is the companion to the daily bug scan: it asks
what is quietly becoming harder to maintain, not only what broke.

Preferred first pass:

```bash
scripts/maintenance_chore.py drift-scan --base-ref REF
```

The first implemented scan is intentionally narrow: it warns when shared
engine/protocol/shell-contract code changed in `REF...HEAD` without obvious
separate test-file evidence in the same diff. Broader adapter-parity,
fixture-provenance, generated-artifact, pure protocol-doc, and inline Rust-test
heuristics remain manual until they can be made low-noise.

Evidence to inspect:

- Adapter files in GUI, CLI, shared shell, MCP, JS, Lua, and Python wrapper paths.
- Changes to `src/engine*`, `src/*shell*`, `crates/*`, and protocol records.
- Changed files under `tests/`, `test_files/`, and nearby fixture manifests.
- Diffs touching `docs/protocol.md`, `docs/gui.md`, `docs/cli.md`, and
  `docs/architecture.md`.
- Tutorial source files, generated chapters, tutorial catalog/manifest metadata,
  and workflow examples that define user-facing tutorial behavior.
- New command routes, operation records, report records, workflow examples,
  generated artifacts, or runbooks.

Findings to report:

- Adapter-owned business or biology logic that should be engine-owned.
- Duplicate GUI/CLI/MCP/shell paths for the same behavior.
- Protocol or schema changes without matching docs and deterministic tests.
- New engine behavior exposed in one adapter but invisible in other intended
  machine-facing surfaces.
- Behavior changes without unit or integration coverage.
- Fixtures without origin, recreation/retrieval steps, or usage notes.
- Tests depending on live network resources without an offline or skipped
  deterministic path.
- Generated artifacts committed without documented source commands.
- Tutorials whose prose drifted from the implementation, whose workflow steps
  omit important biological interpretation, or whose wording is unnecessarily
  hard for an experienced wet-lab biologist to follow.
- Tutorial review/sign-off metadata that is stale relative to the tutorial
  source-version date, or missing when a tutorial was materially rewritten.
- Stack-hardening workarounds that should move from active to superseded in
  [`decisions.md`](decisions.md).
- Deferred probe-region structure drift: once PM-probe import/projection/
  interpretation/backend-execution behavior stops moving, split
  `src/engine/io/probe_regions.rs` into import, inspect/render, projection,
  interpretation, and backend-execution modules without changing behavior.

Pass criteria:

- Each finding names the exact duplicated path, missing contract, uncovered
  behavior, fixture gap, or drifted adapter.
- Suggested fixes preserve the shared engine contract and avoid refactors unless
  the drift blocks release confidence.
- Every committed fixture has provenance.
- Every completed behavior claim has a deterministic test path or an explicit
  documented test gap.
- Tutorials touched by the scan are either confirmed current/readable or have a
  minimal follow-up/fix recorded with evidence. Codex readability reviews may be
  recorded as `codex` sign-offs, but human scientific approval remains distinct.

## Chore 4: Release-Gate Readiness Scan

Cadence: before tags, before release-candidate handoff, and weekly while the
release gate is active.

Purpose: keep the roadmap release gate aligned with runnable evidence.

Preferred first pass:

```bash
scripts/maintenance_chore.py release-gate
```

The first implemented scan verifies compact roadmap release-gate sections,
local Release Gate links, core release-facing docs, and `cargo check` alignment
between the roadmap and `docs/release.md`. It does not validate generated TP73
proof artifacts yet; those remain manual unless a future artifact-manifest
option is added.

Evidence to inspect:

- [`roadmap.md`](roadmap.md) release gate and next-session priorities.
- TP73 pancreatic benchmark runbooks and generated proof artifacts.
- CUT&RUN release-smoke document and any referenced reports.
- Release notes, install/run docs, and container/headless deployment docs.
- Targeted tests, `cargo check -q`, GUI smoke notes, and known skips.

Findings to report:

- Roadmap claims that are not backed by artifact paths or commands.
- Release notes describing features that are not demonstrated by current proof
  artifacts.
- Missing or stale install/run instructions for the demonstrated path.
- GUI smoke gaps that affect demonstration confidence.

Pass criteria:

- A new session can name the release blocker or continue condition after reading
  only the roadmap.
- Release claims cite existing files, commands, or generated artifacts.

## Chore 5: Monthly Decisions And Parity Review

Cadence: monthly, after major refactors, and after adding an operation or
command family.

Purpose: keep durable decisions current and adapter parity explicit instead of
folkloric.

Evidence to inspect:

- `Status: active | superseded | deferred` entries in [`decisions.md`](decisions.md).
- `docs/glossary.json`, `docs/protocol.md`, `docs/cli.md`, and `docs/gui.md`.
- MCP tool names, shared shell command families, and JS/Lua/Python wrapper
  surfaces where applicable.
- Current code paths for stack-hardening, shell dispatch, adapter boundaries,
  screenshot policy, ClawBio/OpenClaw boundaries, and fixture provenance.
- Changelog entries that may have made old constraints obsolete.

Findings to report:

- Active decisions that are no longer true.
- Deferred decisions that have become current roadmap priorities.
- Workarounds that should become `superseded`.
- Important capability exposed in one surface with no intended support or
  explicit non-support note in the others.
- Divergent command names, result fields, or confirmation wording.
- UI-intent targets missing from menu, command palette, shell, agent, MCP, or
  ClawBio/OpenClaw discovery.

Pass criteria:

- Active decisions are few, current, and implementation-guiding.
- Obsolete constraints remain available for history but no longer govern new work.
- Missing parity is either fixed, documented as intentionally unsupported, or
  added to the roadmap as one-line future work.

## Implementation Order

1. Keep the Daily Bug Scan as the defect-oriented daily chore.
2. Add Session-Close Hygiene as the default end-of-session checklist.
3. Add Weekly Drift, Test, And Provenance Scan as the first weekly maintenance
   chore.
4. Add Release-Gate Readiness Scan for pre-tag periods.
5. Add Monthly Decisions And Parity Review after the weekly scan is stable.

Inspector-friendly success condition: an external reviewer can trace current
priorities, completed work, durable rules, tests, artifacts, and known gaps
without reading a monolithic planning document or trusting undocumented claims.
