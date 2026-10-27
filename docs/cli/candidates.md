# Candidate Sets

<!-- manual-counterparts:start -->
Related GUI reference: [Candidate-Set Workflow (Engine Ops + GUI Shell)](../gui/candidate-set-workflow-engine-ops--gui-shell.md) | [Design Constraints Filter (Engine Ops)](../gui/design-constraints-filter-engine-ops.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Candidate-set capability status:

- `gentle_cli`: supported as first-class `candidates ...` commands and shared-shell `candidates ...` commands, backed by shared engine operations (`GenerateCandidateSet`, `GenerateCandidateSetBetweenAnchors`, `DeleteCandidateSet`, `ScoreCandidateSetExpression`, `ScoreCandidateSetDistance`, `ScoreCandidateSetWeightedObjective`, `TopKCandidateSet`, `ParetoFrontierCandidateSet`, `FilterCandidateSet`, `CandidateSetOp`, `UpsertCandidateMacroTemplate`, `DeleteCandidateMacroTemplate`)
- `gentle_js`: supported via `apply_operation` with the same candidate-set operations
- `gentle_lua`: supported via `apply_operation` with the same candidate-set operations


- Shared shell:
  - Commands:
    - `candidates list`
    - `candidates delete SET_NAME`
    - `candidates generate SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite] [--limit N]`
    - `candidates generate-between-anchors SET_NAME SEQ_ID --length N (--anchor-a-pos N|--anchor-a-json JSON) (--anchor-b-pos N|--anchor-b-json JSON) [--step N] [--limit N]`
    - `candidates show SET_NAME [--limit N] [--offset N]`
    - `candidates metrics SET_NAME`
    - `candidates score SET_NAME METRIC_NAME EXPRESSION`
    - `candidates score-distance SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite]`
    - `candidates score-weighted SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]`
    - `candidates top-k INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
    - `candidates pareto INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
    - `candidates filter INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]`
    - `candidates set-op union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET`
    - `candidates macro [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
    - `candidates template-list`
    - `candidates template-show TEMPLATE_NAME`
    - `candidates template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]`
    - `candidates template-delete TEMPLATE_NAME`
    - `candidates template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]`

Candidate-set commands (`gentle_cli candidates ...` and `gentle_cli shell 'candidates ...'`):

- `candidates list`
  - Lists available candidate sets from project metadata.
- `candidates delete SET_NAME`
  - Removes one candidate set.
- `candidates generate SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite] [--limit N]`
  - Creates a candidate set from fixed-length windows on `SEQ_ID`.
  - `--feature-kind` can be repeated to constrain nearest-feature context.
  - `--feature-geometry`, `--feature-boundary`, and `--strand-relation`
    control how directed feature distance targets are selected.
- `candidates generate-between-anchors SET_NAME SEQ_ID --length N (--anchor-a-pos N|--anchor-a-json JSON) (--anchor-b-pos N|--anchor-b-json JSON) [--step N] [--limit N]`
  - Creates a candidate set from fixed-length windows constrained to the
    interval between two local sequence anchors on `SEQ_ID`.
  - Anchor JSON accepts the same local `SequenceAnchor` schema used by
    `ExtractAnchoredRegion` (`Position` or `FeatureBoundary` with
    `Start|End|Middle`).
- `candidates show SET_NAME [--limit N] [--offset N]`
  - Pages records (`sequence`, `coordinates`, `metrics`) from one set.
- `candidates metrics SET_NAME`
  - Lists available metric names in one set.
- `candidates score SET_NAME METRIC_NAME EXPRESSION`
  - Computes a derived metric expression for all records in a set.
- `candidates score-distance SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite]`
  - Computes nearest-feature distance metric with optional feature filters.
- `candidates score-weighted SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]`
  - Computes one weighted objective metric from existing metrics.
  - `--term` can be repeated; default direction is `max`.
  - `--normalize` is enabled by default (min-max scaling per term).
- `candidates top-k INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
  - Materializes explicit top-k selection for a metric.
  - deterministic tie-break policy avoids unstable ordering.
- `candidates pareto INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
  - Materializes a Pareto frontier for multi-objective optimization.
  - optionally truncates with deterministic tie-break (`--max-candidates`).
- `candidates filter INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]`
  - Creates `OUTPUT_SET` by value and/or quantile thresholds.
- `candidates set-op union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET`
  - Creates set algebra output from two sets.
- `candidates macro SCRIPT_OR_@FILE`
  - Runs multiple candidate statements in order (semicolon/newline separated).
  - Nested macro calls are rejected.
- `candidates template-list`
  - Lists persisted candidate macro templates.
- `candidates template-show TEMPLATE_NAME`
  - Shows one persisted template definition.
- `candidates template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]`
  - Creates/updates a named template in project metadata.
  - placeholders in script use `${param_name}` and must be declared via `--param`.
  - Optional `--details-url URL` records external protocol/reference details for
    catalog display.
- `candidates template-delete TEMPLATE_NAME`
  - Deletes one persisted template.
- `candidates template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]`
  - Expands a named template with provided bindings/defaults, then executes it as
    a candidate macro script.


Notes:

- In-memory candidate sets persist in
  `ProjectState.metadata["candidate_sets"]` (`gentle.candidate_sets.v1`).
- On save, candidate sets are externalized into a sidecar index + JSONL files;
  project metadata stores a reference schema (`gentle.candidate_sets.ref.v1`).
- On load, sidecar-backed candidate metadata is hydrated automatically.
- `list/show/metrics` are read-only commands.
- `macros run/template-put/template-delete/template-run` and
  `delete/generate/generate-between-anchors/score/score-distance/score-weighted/top-k/pareto/filter/set-op/macro/template-put/template-delete/template-run`
  mutate state and are available through CLI shell and GUI shell.

Recommended suffixes:

- pool artifacts: `*.pool.gentle.json`
- project/state artifacts: `*.project.gentle.json` (or existing `*.gentle.json`)

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Experiment Planning](planning.md) | [Next: Guide Design](guides.md)
<!-- manual-page-nav:end -->
