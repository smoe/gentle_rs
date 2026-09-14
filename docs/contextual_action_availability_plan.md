# Contextual Menu And Action Availability

Status: slices 0 and 1 are committed as `89050bac`. Slices 2 and 3 now have
the shared readiness/menu/palette implementation and explicit gRNA bindings;
verification results are recorded below. This is a closed pilot, not completed
application-wide availability.
Claude's user-forwarded review inspected `98f74ee8`, 2026-09-14. Codex checked
the feedback against that source and rechecked the relevant app paths after
HEAD advanced to `a5b7028157493d44c1812fc8e1f39e1ccbea6f6a`.
This does not expand the `.10` release gate. The original proposal and review
response are retained below; the revised plan has not had a second Claude
review. No further direct Claude connection is part of this task.

## Evidence And User Intent

The user wants menus, their nested submenus and individual actions to reflect
the current project/selection, explain unavailable actions, and spell gRNA
correctly. The screenshot is evidence of confusing presentation, not a request
to run a guide-design workflow. Line references below describe the reviewed
`98f74ee8` source; prefer the named functions if subsequent edits move them.

- `src/app.rs:15701` title-cases filenames: `grna` becomes `Grna`.
  `assets/cloning_routines.json` already has the correct gRNA titles.
- `src/app.rs:15782` / `16930` render the Patterns hierarchy as template-import
  buttons. Clicking a pictured candidate scan imports a macro; it does not
  execute it. Importing is useful even without DNA. Keep that distinction.
- `src/app.rs:15823` routine-catalog entries also import templates.
- The anchor-scan template requires `seq_id`, `anchor_a_pos`, `anchor_b_pos`;
  its routine ports use `seq_id`, `anchor_a`, `anchor_b`. Port-to-parameter
  binding is not safely inferable from equal names alone.
- `src/app.rs:8725` exposes active DNA-window context; it is not the same as
  selecting a sequence row in the main project view. Several launchers fall
  back to the first ID from `project_sequence_ids_for_blast`, which lists all
  molecule kinds. This can prefill the wrong subject for a DNA-only action.
  Audit each consumer: the same helper also supports non-DNA workflows.
- `src/app.rs:22763` paints command-palette rows and accepts Enter/click without
  a common availability result. Disabling menu buttons alone is insufficient.
- `src/engine/protocol.rs:2642` and `src/engine/state/sequence_ops.rs:1931`
  already provide `FactExpression`, `FactTruth` and deterministic evaluation.
- `src/engine_shell.rs:31878` binds capability prerequisites and returns
  `ready`, `blocked` or `unknown`; it currently projects state inside each
  capability evaluation. Do not call it repeatedly from paint.
- `src/app/routine_and_agent_assistant_ui.rs:1613` already blocks agent
  suggestions on unmet/unknown supplied expressions. Routine execution also
  has a separate preflight. Neither is a complete shared menu policy today.
- `src/app/collection_operations_ui.rs` already separates collection-policy
  rejection, missing bindings and unavailable GUI adapters. Preserve this.
- `src/about.rs` owns limited native macOS menu bridges; the pictured Patterns
  hierarchy is egui. Native bridges must not become another policy engine.
- `collect_cloning_pattern_catalog_entries` performs recursive directory I/O
  from the Patterns menu closure (`src/app.rs:16965`), with a relative default
  catalog path. This is existing work in paint, not merely a prospective risk.
- `capability_precondition_expr_value` (`src/engine_shell.rs:31620`) allocates
  and searches the static descriptor registry on each lookup.
  `introspection_project_graph` also loads the agent catalog and probes
  executable availability. These are not pure project-state lookups.

## Original Codex Plan

Historical proposal, superseded where the review response below says otherwise:

1. Clarify action purpose and fix display names first. Keep import/browse/setup
   accessible where they can genuinely work without project DNA.
2. Add a small shared action descriptor and availability result using existing
   capabilities, fact evaluation and routine preflight, not a second rules DSL.
3. Resolve the invocation's explicit subject/context once, evaluate cheap
   prerequisites, and project the same decision into menus, palette, context
   actions and agent/UI-intent dispatch.
4. Aggregate nested menus from their actionable descendants and explain why
   disabled choices cannot run. Recheck at invocation; engine validation and
   confirmation remain authoritative.
5. Pilot with gRNA candidate scans/oligos, template imports, and a small set of
   existing subject-bound launchers; add deterministic transitions, parity and
   no-work-in-paint tests before migrating additional families.

## Claude Feedback And Codex Response

Adopt Claude's main simplification: generalize the existing app-layer
`CollectionLauncherReadiness` presentation instead of introducing protocol
DTOs and a parallel action registry in the pilot. Existing capability facts,
collection policies and routine preflight remain authoritative. Also move
catalog I/O and incorrect subject fallback fixes ahead of the abstraction.

Adopt canonical routine titles keyed by `template_path`, separate project
facts from host-resource probes, and cover both command-palette render paths.
Do not repeat the agent's per-suggestion fact-graph construction per menu row.

Three qualifications matter:

- **Focus is not already solved.** Claude missed the helper call chain:
  `set_active_window_viewport` reports focus, and
  `note_viewport_focus_if_active` calls it (`src/app.rs:3309`). The separate
  palette calls that helper (`src/app.rs:22879`), and the root frame explicitly
  selects `ViewportId::ROOT` when focused (`src/app.rs:25636`). Therefore
  "only DNA reports focus" and "the main window never clears DNA context"
  are not supported. Keep explicit invocation-origin binding and test both
  accidental target loss and stale target retention, in both window modes.
  No live focus regression has yet been demonstrated in this review.
- **Filter per action, not globally.** The sequence-list helper is shared
  with BLAST, UniProt and input choosers. Do not turn it into a DNA-only list.
  For subject-bound actions, explicit selection wins; otherwise use a real
  chooser or explain the missing input. Kind filtering alone does not make
  an alphabetically first sequence the user's intended subject.
- **A dim row is not an invocation guard.** Put availability on collected
  palette entries for presentation, but also recheck at the shared dispatcher.
  That covers both Enter/click paths and state changes after collection.

Claude ran no build or tests. Source inspection supports these revisions,
not a claim of measured latency or a reproduced GUI focus failure.

## Revised Abstraction

### An Action Is More Than A Capability

Reuse existing capability/routine IDs, `CommandPaletteAction` and UI-intent
targets, not translated labels or menu paths. Add only the small app binding
adapter needed to associate a migrated action with its purpose, explicit
subject and existing prerequisites. Do not create another capability registry.
Keep execution mode/confirmation separate from availability. A capability ID
alone is insufficient to distinguish importing a template from executing it.

Distinguish purposes explicitly:

| Purpose | Example | Meaning of enabled |
| --- | --- | --- |
| Discover/import | Import gRNA Anchor Window Scan template | The template can be inspected/imported; no claim of runnable design |
| Configure | Open Configuration or choose a reference | User can supply missing input; no DNA required just to open setup |
| Launch with subject | Use routine with selected DNA | The explicitly selected subject is suitable for the parameter form |
| Preview/execute | Generate candidates | Required bindings and applicable engine preconditions/preflight are current |
| Inspect/export | Export an existing result | Required result exists, is usable for this action and has the required provenance |

Do not give every action every stage or require fully bound execution parameters
just to open a form that collects them. Conversely, an empty form that cannot
select an input is not a useful enabled action.

Generalize `CollectionLauncherReadiness` into a reusable app-layer presentation
type, retaining its typed collection-policy rejection reasons and existing
binding/adapter distinctions. Add missing-subject, wrong-kind, ambiguous-input
or checking/unknown cases only where the pilot needs them. Preserve the
existing collection tests; do not flatten domain reasons into generic strings.

The presentation carries the bound subject/purpose and snapshot token needed
for stale-state checks. Reasons have localized text and a recovery route where
available. Keep capability support, adapter support and permission separate.
Do not use translated labels as action IDs.
Availability is not authorization, an experimental-success prediction, or a
proof that external files remain unchanged.

### Shared Evaluation, Small Host Context

No new `gentle-protocol` DTO in this pilot. Keep biological rules, fact
evaluation and collection policies headless; the app layer resolves GUI
context and presents their results. Preserve the root-engine compatibility
boundary and never add an engine-to-egui dependency. Reuse/extract only the
needed capability argument binding from the shared shell; do not execute shell
commands or serialize/parse JSON merely to paint a menu. A portable availability
result can follow when a concrete non-GUI consumer needs that wire contract.

The host supplies a bounded context: explicit clicked row or launch-window
subject, selected range/feature IDs if needed, pending-load state, and known
host/resource/job status. Validate subject identity and molecule kind against
engine state. Project inventory is not an implicit selection.

- A context-menu action uses its clicked object. A DNA-viewer action uses that
  viewer. A main-window action uses its explicit selected project object(s).
- Opening a menu/palette must preserve the launch origin; taking focus must
  not accidentally remove the DNA context. Closing/removing/changing the
  originating object must invalidate that context.
- Opening from the main window must not borrow a previously focused DNA
  viewer. Capturing an origin for a palette invocation is not a persistent
  global "last DNA" fallback.
- Do not silently choose the first sequence or resolve multiple selections to
  one. Offer a chooser where supported; otherwise disable with a reason.
- Distinguish a selected DNA molecule from a highlighted DNA range, an anchor,
  protein, logical set, pool, report or guide set. Use existing collection
  policies for multi-subject operations rather than treating all as DNA.
- CLI callers provide explicit bindings; they never inherit a GUI selection.
  Missing GUI context is unknown/unavailable for a GUI-only intent, not a
  reason to disable ordinary explicit headless operations.

Capability prerequisites are authoritative where annotated. Required routine
ports establish missing bindings, not complete biological validity. For the
pilot, add explicit, validated bindings/kind constraints only where existing
metadata is insufficient; never infer them from filenames or descriptions.
Old/custom templates remain importable and inspectable; absent readiness
metadata means unknown, not ready. Do not execute macros to discover eligibility.
This is not a blanket new gate on all historical shell commands. Keep unmigrated
routes on their existing validation path and report registry coverage honestly;
enforce the new invocation policy only for explicitly migrated actions. Unknown
prerequisites of a migrated execution action block that action, not its separate
discovery/import/configuration route.

### Menus And Feedback

Keep unavailable actions visible, dimmed and named consistently. Explain the
reason on disabled hover where supported, in palette detail rows, and through
context help/accessibility descriptions. A disabled native submenu may not
deliver hover: retain a reachable action-help/browser surface for its reasons.

Compute submenu state bottom-up from actionable children, excluding separators
and captions. Enable a submenu if at least one child can perform its advertised
action now; disable an empty/all-blocked/all-unknown subtree. A browse/import
child can legitimately keep a category enabled even if all execution actions
are blocked. Parent summaries should describe this rather than assert that a
whole scientific family is unsupported.

Recovery actions must remain reachable outside a disabled parent. Do not make
"Configure reference" inaccessible because reference-dependent analysis is
disabled. Keep root File/Help/Configuration paths useful in an empty project.

Menu, palette Enter/click, keyboard shortcut and GUI-intent routes that reach
a migrated action call the same invocation guard. Re-evaluate against fresh
bound state; reject stale/unready actions with typed reasons. Availability on
`CommandPaletteEntry` controls display, not permission to skip dispatch checks.
The engine still validates its actual operation and resource/preflight inputs.
Agent-supplied expressions may add constraints, never override canonical
requirements; no UI readiness result grants auto-execution or permission.
Unifying the agent's expression evaluator and native-menu presentation is
deferred, not a claim of completed whole-application parity.

### GUI Cost And Freshness

No filesystem walks, template parsing, sequence hashing, full-project graph
projection, BLAST/Primer3 probes or preflight execution in paint. The existing
catalog walk/routine listing should move behind an explicit cached refresh for
the migrated hierarchy. Evaluate only requested actions and necessary subjects.

Separate three kinds of data rather than caching the existing introspection
function wholesale:

- Immutable capability descriptors: one initialized lookup of the existing
  registry, such as `OnceLock`, without changing its serialized contract or
  keeping a second independently maintained definition.
- Project-derived facts: a bounded snapshot keyed by relevant project state,
  shared across requested actions rather than rebuilt for each row.
- Host/catalog/resource checks: explicit refresh with its own generation and
  pending/error state. Resource checks cannot be cached indefinitely by an
  engine revision, and a static descriptor cache must not capture probe results.

Cache presentation using project instance identity + relevant mutation revision,
selection/origin generation, catalog generation and resource/job generation.
An engine revision alone is insufficient across project replacement or focus
changes. Do not use read-only execution counters as invalidation triggers.
If a needed lock/snapshot is unavailable, report checking/unknown rather than
blocking paint or treating missing data as ready. Background results must match
their input token before acceptance; no cached scientific preflight becomes an
evergreen Boolean. Busy state is action-specific, not a global GUI lock.

## Minimal Implementation Slices

### Slice 0: Existing Bugs, No New Framework

Move catalog discovery out of the paint closure. Resolve the default through
the existing resource-location mechanisms, independent of launch directory;
preserve explicit path overrides. Cache a snapshot with explicit refresh and
visible loading/error states, not a background refresh on every menu frame.

Inventory the `.first()` consumers and their actual molecule/input contracts.
Fix unsafe subject-bound launches with explicit selection, a useful chooser,
or a missing-input reason. Preserve legitimate multi-kind choosers. Add direct
regressions before touching the readiness abstraction.

### Slice 1: Truthful Labels

Prefer routine titles joined by normalized `template_path`; use an acronym
fallback for directory segments and uncatalogued imports (`gRNA`, `CRISPR`,
`PCR`, `DNA`, `RNA`, `TFBS`). Test all current catalog entries. Do not rename
commands, files, template IDs or workflow keys. Label leaves "Import ...
template" and the subtree as template import. Do not change import into run.

### Slice 2: One Readiness Presentation, Two Surfaces

Generalize the existing app readiness enum and reuse static descriptor lookup,
pure project facts and separately refreshed host facts for a closed pilot.
Start with the Patterns hierarchy and the corresponding palette actions plus
the audited subject-bound launchers from slice 0. Include explicit invocation
origin, collection-policy compatibility, submenu aggregation, disabled reasons
and a fresh dispatcher guard. Exercise both palette rendering paths.

Document migrated coverage. Where agent/UI-intent routes invoke these same
actions, use their guarded dispatcher; do not expand this slice into a rewrite
of the agent's expression evaluation or all shell capability prerequisites.

### Slice 3: gRNA Subject Binding

Reuse the Routine Assistant for a distinct "Use with selected DNA..." action
with validated bindings for both candidate scans. Resolve the anchor port-name
mismatch explicitly. Anchor coordinates are collected before Run, not required
just to open a useful form. The oligo route requires a guide set, not arbitrary
DNA. No new guide-design, PAM, specificity or off-target biology is in scope.

Defer protocol DTOs, native-menu presentation migration, agent-expression
unification and broader capability/toolbar migration. Preserve their existing
validation and permission rules. Update affected interface docs and changelog
as slices land, without expanding the release gate or calling a pilot complete
application-wide parity.

## Acceptance Cases

- Empty project: import/open/help/configuration work; "Use with selected DNA"
  is disabled with a recovery path. No phantom enabled execution action.
- DNA loaded but not selected: no implicit binding; after explicit selection,
  a subject-bound launcher enables. RNA/protein/pool/report are not DNA by
  filename or availability of an open window.
- Protein-only and mixed projects: a DNA-only launcher cannot silently bind
  the alphabetically first protein or choose an unrelated DNA. Multiple
  selected subjects require the advertised collection policy or a chooser.
- Anchor scan setup opens on selected DNA; Run stays unavailable until both
  anchors and other required parameters pass the existing bound validation.
- An existing guide set enables its oligo route even without a selected DNA
  viewer. Selecting unrelated DNA alone does not enable that route.
- Nested menus correctly aggregate mixed and empty descendants; recovery
  remains reachable. Disabled leaf/submenu reasons are inspectable.
- Deleting/replacing DNA, Undo/Redo, switching projects, changing selection,
  receiving stale background checks or changing resource status invalidates
  readiness. An open menu/palette cannot target an unintended window.
- Palette launched from DNA retains its explicit initiating subject; palette
  launched from the main window requires that window's selection. Closing the
  origin invalidates it. Test embedded and separate window paths, without
  assuming that focus is always sticky or always lost.
- Mouse, Enter and shortcuts agree on bound prerequisites for migrated actions,
  as do GUI intents that invoke their dispatcher. In both palette paths,
  attempted disabled/stale execution produces no mutation or job launch.
- Unknown legacy metadata is not mislabeled as a biological rejection or a
  verified pass. Missing permissions remain separate and cannot be bypassed.
- gRNA spelling comes from shared presentation metadata; old IDs round-trip.
- Counter-based tests prove no repeated catalog I/O, per-item whole-project
  projection or external preflight in steady menu frames. Glen measures actual
  responsiveness using the existing audit profile, not a new noisy CI threshold.
- Launch-directory changes do not break the packaged/default catalog; explicit
  invalid overrides yield an honest error. Refresh invalidates the snapshot.
- Descriptor lookup does not rebuild the registry per action. Resource refresh
  changes readiness even without a project edit; stale probe results are ignored.

Testing order: catalog/subject regressions, shared readiness/lookup tests,
host-dispatch guard tests, headless egui frames (including disabled Enter),
existing routine/agent
regressions, `cargo check`, formatting/whitespace, then Glen's live GUI check.
Use inline synthetic fixtures or documented fixture provenance. No real
guide design, private project upload or screenshot capture is needed for
these deterministic tests.

## Claude Consultation

An earlier direct read-only attempt failed authentication and produced no
review. The user then forwarded the
[original review prompt](contextual_action_availability_review_prompt.md)
and supplied Claude's source-based critique on 2026-09-14. The summary above
separates that feedback from Codex's corrections and revised scope.

The user approved starting implementation after this review. Slices 0 and 1
also include explicit palette-origin capture because removing inventory
fallbacks must not break subject-bound launches when the palette takes focus.
This is not the general readiness-presentation migration in slice 2. The
revised proposal remains available for user-forwarded Claude review if desired;
do not contact Claude directly.

## First-Slice Implementation Record

- Catalog snapshot/refresh and labels: `src/app/pattern_catalog_ui.rs` and
  `src/app.rs`, using shared validation in `src/engine_shell.rs`.
- Explicit subjects and palette origin: `src/app/subject_selection.rs`,
  `src/app.rs`, and the existing molecule-kind helper exposed through
  `src/engine/state/sequence_ops.rs`.
- Shared asset lookup: `src/runtime_assets.rs`, `src/lib.rs` and
  `src/bin/gentle.rs`; explicit import paths do not use default-asset fallbacks.
- Regression coverage: inline tests in the new modules and `src/app/tests.rs`.
  Existing PCR/confirmation dispatch tests now check rejection without selection
  and successful opening after explicit selection, rather than requiring the
  removed inventory fallback.
- Necessary feedback fix: `src/app/routine_and_agent_assistant_ui.rs` preserves
  subject-launcher messages alongside the UI-intent identifier. Otherwise the
  generic dispatch summary hid the new missing-subject explanation.
- User-facing/status documentation: this plan, `docs/gui.md`, `docs/CHANGELOG.md`
  and `docs/roadmap.md`. No dependency, fixture artifact, template ID, assay
  algorithm or release-gate change.

Verification on 2026-09-14: the initial focused set passed 21/21; the final
broader command below passed 595/595 app/runtime-asset tests, including both
palette-origin modes, stale origins, explicit subject transitions, catalog
worker/cache behavior and packaged-path lookup.

```sh
cargo test --lib --locked -j 1 -- app:: runtime_assets::tests --test-threads=2 --quiet
cargo check -q --locked -j 1
cargo fmt --all --check
git diff --check
```

All four commands passed. Session-close hygiene reported 4 OK, 2 warnings and
0 failures: intentional uncommitted changes and the manual plan-fidelity
reminder. The palette-origin and UI-intent feedback additions above are the
explicitly documented dependencies of removing unsafe launch fallbacks.

This is headless egui/dispatch coverage, not live macOS focus acceptance or a
GUI performance measurement. Full workspace/release gates were not rerun.
The beta toolchain also reported the existing nightly-only `lints.cargo`
manifest warning and a large debug-test unwind-table linker warning.

## Readiness And Binding Slices

- `src/app/action_readiness.rs` lifts the collection readiness presentation
  without changing its typed rejection reasons. One bounded context per surface
  supplies DNA/sequence/guide-set readiness; project and viewer locks are
  nonblocking. The pilot never calls whole-project introspection or host probes
  from paint. The broader agent/introspection host-probe path is not rewritten.
- Menus aggregate children bottom-up; imports remain separate. Both palette
  paths use the same disabled-row/Enter policy and a fresh dispatch guard.
  `docs/gui.md` enumerates migrated actions and the already-open-window and
  empty-project retrieval exceptions.
- `src/app/grna_routine_ui.rs` opens the existing Routine Assistant with the
  explicitly selected DNA or a guide-set chooser. Catalog port/template checks
  run in the background snapshot and again on explicit invocation. Missing
  anchors block preflight, not setup. A previous gRNA preflight cannot authorize
  execution after its project revision or bound parameters change.
- `src/engine_shell/routine_bindings.rs` explicitly maps anchor port names,
  validates through the engine's existing anchor geometry, and is reused by
  macro preflight. No new PAM/off-target or guide-design algorithm is implied.
  Static descriptor lookup now uses `OnceLock`; it does not cache mutable host
  availability. Default preflight catalog lookup also uses packaged asset
  resolution, required for the same workflow outside the checkout directory.
- Existing collection launchers, engine validation, confirmation and non-pilot
  actions retain their prior contracts. Native menus, broader capability
  migration and agent-expression unification remain deferred. Glen's live GUI
  and timing acceptance remains separate from headless regression tests.

Verification on 2026-09-14: 634 app/runtime/macro/guide-design tests passed,
followed by 79 introspection, guide-shell and anchor-operation regressions.
The two commands below cover 713 tests, with no failures or ignored tests.
Initial focused testing caught a cache semantic regression (duplicate descriptor
IDs must keep the former first-match behavior) and an incorrect test assumption:
rejected shell macro execution intentionally preserves a failed lineage receipt,
whereas a disabled GUI launch and `--validate-only` do not mutate the project.
The fixes preserve both existing contracts rather than suppressing failure audit.

```sh
cargo test --lib --locked -j 1 -- app:: runtime_assets::tests engine_shell::routine_bindings::tests engine_shell::tests::execute_macros engine_shell::tests::parse_macros engine_shell::tests::execute_introspection engine::tests::test_guide_design --test-threads=2 --quiet
cargo test --lib --locked -j 1 -- engine_shell::tests::execute_introspect engine::tests::test_generate_candidate_set_between_two_sequence_anchors engine_shell::tests::execute_guides --test-threads=2 --quiet
```

The first command's `execute_introspection` filter matches no current test names;
the second deliberately covers the actual `execute_introspect` names.
No full-workspace, release packaging or measured latency claim is made.

Final `cargo check -q --locked -j 1`, `cargo fmt --all --check` and
`git diff --check` passed. Session-close hygiene reported 4 OK, 2 warnings,
0 failures before commit: intentional edits and the manual plan-fidelity
reminder. The small subject-lock, preflight asset-resolution and stale-status
changes above are dependencies of the pilot, not unrelated refactoring.
