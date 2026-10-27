### Routine Assistant

<!-- manual-counterparts:start -->
Related CLI reference: [Macros and Routines](../cli/macros-routines.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- `Patterns -> Routine Assistant...`
  - opens a dedicated staged workflow window for routine application:
    1. goal + candidate search
    2. routine alternative comparison + disambiguation answers
    3. typed binding form from routine input ports
    4. shared-engine preflight preview (`--validate-only`)
    5. transactional run + process run-bundle export
  - compare stage includes explicit disambiguation question prompts from
    routine explain/compare payloads and records user-entered answers in
    routine decision traces.
  - when planning context is active, the goal stage now shows a compact
    helper-aware planning context strip plus candidate planning metrics
    (guardrail pass state, score/fit/time/cost, family-alignment bonus).
    Goal and comparison stages display the shared estimate's hands-on,
    unattended and elapsed ranges, with unknown breakdowns for unmodeled
    routines and an explicit legacy-scalar label. Stock-draw value,
    missing-material purchasing proxy, machine/base cost and known-unit total
    are separately labelled. Only the restriction/ligation routine has a
    synthetic step model; none of these figures is a validated protocol time
    or an orderable quote. The GUI derives no second estimate.
    A separate shared label shows the legacy ranking/time-guardrail proxy
    (family/base duration plus delays). Bench-model availability/duration does
    not change ranking or guardrails, and the proxy is not measured lab time.
    Missing historical policy metadata displays as unknown, not reconstructed.
  - comparison has `Explanation` / `Protocol Overview` tabs. The latter is a
    read-only dependency graph from the exact shared explanation, with zoom,
    fit, scroll/drag navigation, selectable steps, source/range details and
    `Copy SVG` / `Copy JSON` (full explanation). It shows the recorded typical
    critical path; columns are dependency depth, not time. Independent roots
    express prerequisites only, not resource-feasible scheduling. The current
    restriction example is synthetic; other routines/legacy records show
    unavailable, never fabricated bench or material flows. Stale/tampered
    diagrams fail rebinding. The content-keyed view cache never trusts pointer
    identity or a supplied digest alone. Semantic IDs `routine.protocol.open`,
    `.fit`, `.copy_svg`, `.copy_json` and static `.step` support later native
    acceptance; step identity is a separate pseudonymous subject scope, not a
    dynamically allocated widget ID (`.explanation` returns to the text tab).
    Both the graph and its detail panel scroll in a small window. No lineage
    entry, run, supplier call or project edit results.
  - when Routine Assistant is opened from a sequence window with construct
    reasoning, that same strip now keeps the active `construct reasoning`
    sequence id plus variant effect tags / suggested assay ids visible, and
    candidate ranking reuses the same sequence-aware routine-preference
    context instead of falling back to planning-objective context alone.
  - compare/parameter stages now surface engine-owned macro-template
    suggestions derived from the selected routine plus the synthesized
    `routine_preference_context`, so users see reusable next-step automation
    ideas without adapter-local heuristics.
  - routine explanation reload also carries the active construct-reasoning
    sequence id, so the compare stage can show the selected routine's own
    sequence-aware planning score / fit / family-alignment bonus instead of
    only the candidate-list snapshot.
  - explainability and comparison data are loaded via shared shell commands:
    - `routines explain ROUTINE_ID [--seq-id SEQ_ID]`
    - `routines compare ROUTINE_A ROUTINE_B`
  - execution is routed through shared macro paths:
    - `macros template-import PATH` (auto-import linked template)
    - `macros template-run TEMPLATE_NAME --bind ... --validate-only`
    - `macros template-run TEMPLATE_NAME --bind ... --transactional`
  - Gibson topology guard:
    - sequence bindings show live topology (`sequence | linear|circular | N bp`)
    - if a Gibson input is circular, preflight and execution are blocked
    - `Linearize Vector...` creates a branched copy, sets topology to linear,
      and auto-rebinds that routine input
  - export stage uses shared process artifact route:
    - `export-run-bundle OUTPUT.run_bundle.json`
  - persisted routine decision traces now retain the helper-aware planning
    context, sequence-aware variant routine preferences, candidate
    planning-score snapshot, and macro suggestions that were visible during
    the assistant flow, so exported run bundles preserve why a routine and
    follow-on automation were favored.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Patterns And gRNA Setup](patterns-and-grna-setup.md) | [Next: Gibson Workspace](gibson-workspace.md)
<!-- manual-page-nav:end -->
