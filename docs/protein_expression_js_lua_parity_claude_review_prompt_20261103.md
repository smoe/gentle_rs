# Claude Review Prompt: Protein Expression Handoff JS/Lua Parity

## Status, 2026-11-03

The owner forwarded this proposal to Claude, supplied its read-only review and
authorized the revised implementation on `fa5de553`. First-class JS/Lua
bindings, shared inline admission, discovery and tutorial 09.03 changes are
now authored; Rust/VM/native verification remains with Glen/CI. The owner
permits generation-only builds, not local Rust checks/tests or VM execution.
Do not push or resume unrelated goals. See the
[execution checklist](glen_next_verification_handoff.md#protein-expression-named-scripting-bindings-2026-11-03).

Correction to the historical proposal below: the existing agent-plan executor
already provided indirect reachability. This slice adds named bindings and
accurate discovery, not a new permission or a general shell API. Also,
JSON-stringifying JS options cannot restore Number precision already lost;
unsafe reviewed timestamps are refused instead of claiming arbitrary u128
transport. Lua uses typed nested decoding for reviewed empty lists and false.

## Original Review Prompt

Prepared 2026-11-03 for Steffen to forward manually, as a read-only proposal
review rather than implementation authorization. No local Claude invocation
was performed by Codex. The original base and diagnosis below are historical.

## Task And Repository

Please critique the minimal Codex plan below before implementation. Verify the
source rather than assuming that a parity-matrix gap proves runtime absence.
Separate routing, discovery, authored tests and executed acceptance.

Repository: `/Users/u005069/.codex/worktrees/faae/gentle_rs`.
Branch: `gentle_rs_3`. Inspected HEAD:
`9c07188d95bbb742502db5af2ac5fe90594b7253`.
Recheck drift using read-only source inspection; this SHA is evidence, not an
audit freeze or permission to reset/rebase any checkout. At preparation the
worktree contains unrelated untracked `outputs/`; preserve it.

Use Read/Grep/Glob only. No shell, builds, tests, VM/native execution, network,
provider/model calls, screenshots, file writes, commit, rebase, push or CI
dispatch. Do not inspect credentials, private projects/sequencing inputs,
`outputs/`, `output/`, `data/` or `paper/`. Read tracked source and the explicitly
named public synthetic examples only. Return your review in your answer, not
a new file. Rust/VM verification belongs to Glen/CI.

Roadmap context: `.13` priority B retains the provider-neutral expression draft
and tutorial 09.03, with explicit provider selection and no submission/contact.
This proposed slice closes only the expression handoff's JS/Lua reachability
and discovery gaps. It adds no release gate or provider capability.

Read `AGENTS.md`, relevant parity/optional-adapter sections of
`docs/architecture.md`, DEC-001/002/003/016 in `docs/decisions.md`,
`docs/roadmap.md`, `docs/testing.md`, and relevant expression sections of
`docs/protocol.md`, `docs/agent_interface.md` and `docs/CHANGELOG.md`.

## Evidence Already Gathered

- `src/engine_shell.rs::execute_planning_command` handles
  `ShellCommand::PlanningProteinExpressionHandoff`, delegates to
  `execute_planning_protein_expression_handoff`, and returns a `ShellRunResult`
  with `state_changed=false` and the shared report. This is a shared-shell
  report, **not an existing workflow/engine Operation**.
- `src/mcp_server/protein_expression.rs` already provides a named read-only
  MCP tool. Its `Arguments` are strict; objective/requirements are inline JSON
  objects, not paths. Explicit null fields are rejected. `command` prepares
  tokens for the shared parser; `call` uses the non-mutating shell-tool route.
  The saved `state_path` is an MCP boundary, not a product-selection mechanism.
- That MCP module has deterministic shared-report/inner-agent equality and
  no-write tests, plus invalid-input tests. The report comparison removes
  only `generated_at_unix_ms`. These are existing test definitions, not a new
  claim that they ran on today's source.
- `src/js_interface.rs` and `src/lua_interface.rs` already expose named
  read-only helpers such as `show_construct_reasoning_graph` through
  `execute_shell_command`, returning `run.output`. No expression/planning
  wrapper or registration was found in either file. Confirm no alternate
  public route already closes this exact gap.
- `docs/glossary.json`'s `planning protein-expression-handoff` row currently
  has `engine_operations: []` and interfaces `cli-shell`, `gui-shell`,
  `gui-menu`, `mcp`. The generated matrix records JS and Lua as `gap`.
  `crates/gentle-protocol/src/lib.rs::glossary_command_surfacing` recognizes
  explicit `js`/`lua` interfaces as prominent; generic operation reachability
  is a separate mechanism. Discovery metadata must describe real exports.
- The canonical glossary now merges duplicate identities and preserves
  per-interface usage variants. Do not create conflicting duplicate rows or
  overwrite the CLI syntax with a scripting signature to remove these gaps.
- `docs/tutorial/09-03_protein_expression_handoff_cli.md` already teaches
  missing product, neutral draft, withheld provider and explicit local provider
  projection through CLI/GUI/MCP/inner-agent paths. It has no JS/Lua examples.
  Its `Interface Parity` table currently has four rows.
- `src/bin/gentle_js.rs` and `src/bin/gentle_lua.rs` are interactive REPLs.
  Their `--project PATH` option loads global `project` via `load_project`.
  They do not advertise a script-file execution option. Scripting is optional
  source-build functionality, not part of the current native/container inventory.

## Original Codex Plan

### A. Add Two Thin Read-Only Bindings

Provisional public shape, for critique rather than a claimed existing API:

```text
JavaScript: protein_expression_handoff(project, options = {}) -> report
Lua:        protein_expression_handoff(project, options)      -> report
```

Use the same public name and report semantics in both languages. Prefer the
existing report-only getter convention over inventing a new response envelope.
Accept an explicit in-memory project; do not silently load/save a state path,
choose the active sequence, or replace the script's `project` variable.

Options should expose only the current read-only operands: optional `seq_id`,
inline `objective`, inline `requirements`, explicit `provider`,
`profile_scope` (currently `effective` only), and `format` (`json` or `text`).
Omitted values retain shared defaults and unanswered questions. The prototype
does not authorize new fields or biological defaults.

Route through the existing shared command/parser/executor and report builder.
Use one narrow transport-neutral admission/token builder if otherwise the two
adapters would duplicate validation; do not depend on a GUI or MCP module from
the engine. Decide whether existing helpers suffice before adding another one.
Do not require a new Operation/schema merely to satisfy the registry, expose
a general arbitrary-shell executor, or refactor the whole planning subsystem.

Match the existing inline-object boundary where applicable: reject unknown
option fields, wrong types, `@FILE`/path strings for objective/requirements,
unsupported format/scope and accidental approval/output arguments. Review the
exact omission/null/empty-object policy before implementation; preserve the
distinction between absent requirements and explicitly reviewed empty topics.
Do not change existing CLI `@FILE` support or MCP validation as a side effect.

Register both the native bindings and their public VM/REPL functions, following
the existing feature gates, help and error conventions. Registering a helper
only in Rust is not scripting reachability. No new dependencies/default
features, package binaries or provider access should be needed.

### B. Preserve The Scientific And Safety Contract

- Return the same `gentle.protein_expression_handoff.v2` report as the shared
  command, including request/product provenance, readiness, missing questions,
  CDS/tag assessment, neutral draft, provider binding and withholding reasons.
  Preserve legacy readability; do not invent a third report family.
- With no product or an unresolved product ID, preserve the shared outcome;
  do not invent a new refusal or select a loaded product automatically.
- Provider omission stays neutral. Naming a provider does not supply missing
  product/requirements, outsourcing permission or catalog capability. The
  existing eligibility rules remain authoritative, including withheld cases.
- Return text-bearing reports only under the same format policy. Suggested
  actions remain data: never execute preflight, quote, optimization, synthesis,
  cloning or ordering. Do not contact providers or models or write artifacts.
- Preserve the supplied project, sequence/annotation content, metadata,
  saved-project bytes and any existing approval/history semantics. Do not add
  journal entries or autosave for this read-only binding.
- Leave the inner agent's existing shared command and instructions intact.
  This is not a reason to shorten its skills or broaden auto-execution.

### C. Close Only Proven Discovery Gaps

Once public exports exist, bind the scripting signatures to the existing
canonical capability using the established interface/usage metadata pattern.
Keep command identity, CLI aliases, GUI/MCP behavior and mutation policy stable.
Update generated projections through their canonical generators, not manual
Markdown relabelling. Derive counts from the typed registry; do not hard-code
today's gap totals or lower maturity/acceptance thresholds.

Check glossary/help/introspection, JS/Lua runtime help, parity freshness and any
necessary glossary-hash projections. Regenerate only outputs required by the
changed source; do not refresh unrelated catalogs, evidence or screenshots.
`planning consult cloning` and other adapter gaps remain outside this slice.

### D. Extend Tutorial 09.03, Not Its Biological Scope

Use the existing synthetic starter and reviewed-requirements examples:

- `docs/examples/workflows/protein_expression_review_starter.json`
- `docs/examples/plans/protein_expression_requirements_in_house.json`
- `docs/examples/plans/protein_expression_requirements_review_only.json`
- `docs/examples/plans/README.md` for their provenance.

Teach why the additional bindings matter: reproduce a reviewed handoff from a
script without creating a second expression planner or granting action consent.
Add JS/Lua rows to the parity table and runnable, clearly language-labelled
REPL examples for missing/neutral product, in-house refusal and permitted local
provider projection. State explicit expected fields and how to inspect/report
them; no implicit provider or example requirements may be loaded.

Show optional-feature/source-build prerequisites without asking users to build
inside the read-only review. Use existing `--project`/`load_project` entrypoints;
do not invent `--script`, script-file execution or a plain-JSON JS project codec.
Keep multiline examples usable with the actual line-oriented REPL. Explain
that the input is the script's loaded snapshot, not unsaved GUI edits. Comparing
all fields requires the same effective project/profile/catalog, operands and
format; exclude only the known generation timestamp.

Preserve the chapter's `manual/hybrid` status. The existing starter workflow
is executable; the shared-shell handoff is not thereby a workflow operation.
Update its source/catalog metadata and teaching checks only as needed. Keep the
four existing screenshot receipts at their actual source; they do not prove
new scripting execution, new-source GUI behavior or biological acceptance.

### E. Author Deterministic Tests; Leave Execution To Glen/CI

Use existing explicit synthetic examples. Test both callable Rust adapters and
the actual registered JS/Lua functions, not registration strings alone:

1. Missing product, usable product with unanswered requirements, neutral draft,
   in-house/selected-provider withholding, and explicitly permitted projection.
2. JSON default versus text-bearing report, saved/effective objective/profile,
   and exact provider/product identity. Compare complete structured reports to
   the shared-shell oracle, excluding only `generated_at_unix_ms`.
3. Unknown/wrong-typed options, null/omission/empty-object boundaries, path-like
   operands, unsupported format/scope and unknown nested requirements. Preserve
   shared scientific outcomes rather than replacing them with adapter guesses.
4. Exact input-state and saved-byte preservation on success and refusal; no new
   files or action execution. Use independent state-preservation oracles.
5. Public export discovery and canonical capability/usage consistency. Keep
   existing GUI/CLI/MCP/inner-agent tests as controls; no live-model test.
6. Tutorial/source/catalog agreement and correct feature prerequisites. Reuse
   existing teaching tests and checkout policy; any new raw-byte freshness
   contract needs its own scoped LF/CRLF coverage.

For JS, load synthetic projects through the public loader/native buffers;
plain JSON literals are not equivalent ingress. For Lua, decode returned
values to the typed Rust representation for equality where appropriate.
Explicitly review empty Lua table versus JSON array/object and false-versus-nil
behavior: `outsourcing.allowed=false` must never become an omitted answer.
Do not alter production codecs to simplify a fixture. Isolate host/catalog
inputs using existing test patterns, not ambient credentials or global-env races.

Implementation bookkeeping: reconcile the existing roadmap line only if status
changes; put completed outcomes in `docs/CHANGELOG.md`, never roadmap Done
entries. Record names/features of authored-but-unrun tests. Rust/VM/native
execution and current-source acceptance stay with Glen/CI; no push/dispatch.

## Questions To Resolve In Your Review

1. Is the missing-binding diagnosis correct? Is there an existing public path
   we missed, and what minimal change makes both execution and discovery true?
2. Is `protein_expression_handoff(project, options)` with report-only output the
   best local convention? Specify defaults, null/nil handling and typed errors.
3. What is the smallest shared admission seam? Can it reuse the current parser
   without duplicating MCP rules or changing its existing behavior? Justify any
   requested new Operation instead of treating it as a parity-table requirement.
4. Which source metadata must change so the two cells can honestly stop being
   gaps? Which generators/tests prove it without falsely promoting maturity?
5. Are project/profile/catalog bindings, timestamp exclusion and state oracles
   strong enough? Identify adapter-specific serialization or error pitfalls.
6. Are the proposed tutorial examples actually runnable in today's REPLs? What
   minimal teaching changes cover input through output without bloating 09.03?

## Requested Reply

Return a concise verdict (approve/revise), severity-ordered findings grounded
in file/function references, then a revised minimal plan with expected changed
files. Distinguish essential fixes from optional follow-ups. Include precise
Rust/VM commands/test filters for Glen/CI, with feature gates and nonzero-count
requirements, but do not run them. Leave packaging, native screenshots, live
providers/models and biological usefulness outside this review's verdict.
