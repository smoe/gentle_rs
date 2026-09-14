# Review Request: Contextual Menu And Action Availability

Historical prompt: the user supplied Claude's review on 2026-09-14. See the
[revised plan and response](contextual_action_availability_plan.md) for the
current proposal; the text below preserves the original review brief.

Please review this GENtle implementation proposal against the actual source,
then recommend the smallest sound abstraction and implementation sequence.
This is a read-only architectural review, not authorization to implement,
commit, rebase, run workflows, or change the release scope.

## Repository And Context

Repository: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Full proposal: `docs/contextual_action_availability_plan.md`.
Initial source inspection: `860d7ef6`; current local HEAD when preparing this
prompt: `98f74ee8`. Record the revision you actually inspect and recheck paths
rather than assuming that either revision is still current. The draft may be
uncommitted; read the working file. If using a different checkout, use its
equivalent files. If source is unavailable, explicitly distinguish a conceptual
review from verified repository findings.

The user wants menus, submenus and individual items enabled/disabled according
to current project state and selection, with helpful explanations. Their
example is a gRNA action offered when no DNA is selected. They also request
correct `gRNA` spelling instead of `Grna`. GUI/inner-agent parity and a responsive
GUI matter. We want to discuss the plan before implementation, without adding
another `.10` release requirement.

## Important Findings To Verify

- The pictured **Patterns > Crispr > Guides > Candidate Scans** leaves currently
  **import macro templates**, not execute guide design. See
  `src/app.rs::render_cloning_pattern_catalog_menu_entries` and its caller.
  Thus "disable everything without DNA" would incorrectly disable useful
  imports. Clarify labels and distinguish import, setup and execution.
- `src/app.rs::humanize_catalog_label` title-cases filenames, causing `Grna`.
  `assets/cloning_routines.json` already supplies correct gRNA titles. Prefer
  canonical display metadata and an acronym-aware fallback; retain stable IDs,
  filenames, command syntax and saved workflow bindings.
- The anchor-scan template needs `seq_id`, `anchor_a_pos`, `anchor_b_pos`, but
  its routine ports use `anchor_a` and `anchor_b`. Do not guess bindings.
  The oligo routine needs a guide set rather than arbitrary selected DNA.
- Existing shared machinery includes `FactExpression`, `FactTruth` and
  `FactEvaluationResult` in `src/engine/protocol.rs`, the evaluator in
  `src/engine/state/sequence_ops.rs`, and capability binding/readiness in
  `src/engine_shell.rs::introspection_readiness_for_capability`.
- `src/app/routine_and_agent_assistant_ui.rs` already gates agent suggestions
  with supplied fact expressions and has separate routine preflight/execution
  checks. `src/app/collection_operations_ui.rs` distinguishes collection policy,
  missing bindings and unavailable adapters. Reuse rather than replace these.
- `src/app.rs::active_dna_window_context` differs from selection in the main
  project view. Some launchers fall back to the first sequence; do not silently
  adopt that behavior. Palette Enter/click also needs guarding, not just menus.
- The pictured menus are egui; `src/about.rs` contains limited native macOS
  bridges. There is no reason to recreate all menus in native code.

## Original Codex Proposal

1. Fix labels and distinguish actions that import/discover, configure, launch
   with a subject, preview/execute, or inspect/export. Opening a useful input
   form need not require every parameter needed for execution.
2. Use a small shared `ActionDescriptor`, keyed by stable action identity,
   referencing existing capabilities/routines, explicit binding policy and
   launch prerequisites. Do not introduce another general rule language.
3. Produce one ephemeral `ActionAvailability`: `ready`, `blocked` or `unknown`,
   with typed reasons, bound subjects, evaluation scope, freshness token and
   optional recovery actions. Availability is not permission or scientific
   success. Use headless helpers; no engine dependency on egui.
4. Bind context explicitly: clicked object, initiating DNA viewer or selected
   project object(s). Preserve launch origin through popup focus changes;
   invalidate removed/replaced subjects. Never silently choose the first DNA.
5. Project that result consistently to menu leaves, palette, context actions
   and matching agent/UI-intent invocations. Recheck on activation; keep actual
   engine validation, resource checks and confirmation authoritative.
6. Aggregate submenus bottom-up: enabled if an actionable descendant is usable,
   otherwise dimmed. Keep reasons accessible even when native disabled menus
   cannot be hovered. Import/help/configuration must remain reachable as recovery.
7. Keep evaluation bounded and cached by project instance/relevant revision,
   binding/selection, catalog and resource/job generations. No disk walks,
   whole-project fact projection per menu item, hashing, external probes or
   scientific preflight in paint. Stale background results cannot enable actions.
8. Migrate a closed pilot set before broader adoption. Do not globally disable
   historical commands lacking metadata or manufacture readiness for them.
   Proposed first delivery: truthful labels/gRNA casing, then the evaluator on
   existing launchers, then explicit gRNA Routine Assistant bindings.

## Questions To Challenge

- Is a new descriptor justified, or can the existing UI-intent/capability and
  routine catalogs express this with fewer additions? Identify the minimal
  types and their ownership, not a large speculative action framework.
- Are import, subject-bound launch and execution separated at the right points?
  What is a useful enabled action in an empty project?
- How should ambiguous selections and focus changes be handled without sticky
  accidental targets or disabling a form that could collect its own input?
- How can canonical requirements coexist with optional agent expressions and
  incomplete legacy metadata without bypasses or widespread regressions?
- Are disabled-parent aggregation and reason/recovery accessibility sound on
  egui and the native surfaces actually present?
- Which existing caches/revisions/preflight contracts can be reused safely?
  Where would the proposal accidentally add menu latency or stale readiness?
- Which first slice gives visible improvement without a new gRNA algorithm or
  an application-wide refactor? What should explicitly remain deferred?

## Requested Response

Return a concise verdict, prioritized concrete findings with source references,
a revised minimal plan, and acceptance tests. In particular cover empty project,
DNA versus protein/RNA/guide set, missing anchors, multiple selections, popup
focus, project replacement/Undo, stale preflight, disabled keyboard activation,
agent parity and no expensive work in steady frames. Preserve uncertain facts
as uncertain. Distinguish repository evidence from your design recommendations.

No runtime changes have been implemented for this proposal. Do not describe
it as landed or imply that an earlier Claude review approved it.
