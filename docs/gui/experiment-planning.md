### Experiment Planning

<!-- manual-counterparts:start -->
Related CLI reference: [Experiment Planning](../cli/planning.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- `Patterns -> Planning...`
  - opens a dedicated standalone `Planning` window for the planning meta-layer.
  - supports editing and applying:
    - `global` planning profile
    - `confirmed_agent_overlay` profile
    - `project_override` profile
    - planning objective
  - shows merged effective profile (`global -> confirmed_agent_overlay -> project_override`).
  - includes a searchable `Host Profile Browser` backed by the shared starter
    host-profile catalog, so propagation/expression host traits remain
    inspectable with genotype/phenotype/source-note details.
  - includes a searchable `Helper Construct Browser` backed by the shared
    helper catalog + normalized interpretation layer, so planning work can
    inspect helper semantics without reparsing free text.
  - the same browser surfaces the shared helper vector-card projection and
    deterministic catalog-doctor issues used by `helpers show-card` and
    `helpers doctor-catalog`, preserving GUI/CLI parity for metadata-only
    vector candidates.
  - the Agent Assistant and GUI Shell can ask the same read-only consultation
    route (`planning consult cloning`) for cloning strategy/vector advice backed
    by the effective profile, host/helper catalogs, and routine estimates;
    v1 keeps `seq_id` traceability-only and turns marker/promoter/MCS gaps into
    explicit follow-up questions instead of hidden heuristics.
  - high-yield protein-expression requests can use the separate read-only
    shell route (`planning protein-expression-handoff`) to surface product
    context, optional reviewed `gentle.protein_expression_requirements.v1`,
    chassis/route candidates, missing yield/folding/chassis questions,
    and a readiness-gated GeneArt preflight scaffold without creating or
    ordering constructs; provider-ready products use a project-bound draft,
    while undefined products and explicitly in-house-only requirements expose
    no provider action.
  - supports registering pending planning sync suggestions (`pull`/`push`) from
    JSON payload and explicit `Accept`/`Reject` resolution in-window.
  - exposes sync status (`pending count`, latest pull/push timestamps, last
    source/snapshot, last error).
- `Patterns -> Routine catalog`
  - routine discovery is grouped by `family` and `status`.
  - selecting a routine imports its linked template file when the routine
    manifest provides `template_path`.
  - Gibson baseline is available as:
    - `gibson.two_fragment_overlap_preview`
    - template: `gibson_two_fragment_overlap_preview`
  - current execution path for imported routines is the shared shell panel:
    - `macros template-run TEMPLATE_NAME --bind ... --validate-only`
    - then run again without `--validate-only` for a mutating execution.
  - manifest source: `assets/cloning_routines.json`
    (`gentle.cloning_routines.v1`).

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Sequencing Confirmation Workspace](sequencing-confirmation-workspace.md) | [Next: Candidate-Set Workflow (Engine Ops + GUI Shell)](candidate-set-workflow-engine-ops--gui-shell.md)
<!-- manual-page-nav:end -->
