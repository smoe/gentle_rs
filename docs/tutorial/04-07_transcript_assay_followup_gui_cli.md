# Follow Up Transcript Assays: RT Primer Pools, Differential Junctions and Specificity

> Type: `GUI + shared-shell walkthrough`; status: `manual/hybrid`.
> Offline fixture/replay checks are not human scientific approval or GUI acceptance.
> Synthetic teaching data only: no orderable human PATZ1 primers, BLAST search,
> vendor submission, or experimentally validated result is produced here.

After [04.03, qPCR across exon junctions](04-03_qpcr_exon_junctions_gui.md)
and [04.06, PATZ1 transcript panels](generated/chapters/04-06_patz1_transcript_assay_panels_cli.md),
the next question is not simply "which primer scored highest?" It is "what
does this design establish, what remains unknown, and what must I approve next?"

This follow-up connects four decisions without conflating them:

| Decision | Engine artifact | What it does not establish |
| --- | --- | --- |
| Where to prime reverse transcription | `gentle.terminal_exon_rt_primer_pool.v1` | Successful RT or isoform-specific capture |
| Which junction to investigate | `gentle.primer_pair_summary.v2` inside a panel | Statistical significance or independent replicated observations |
| Whether selected pairs satisfy a reference-specific policy | Aggregate specificity handoff and acceptance | A pass before execution/finalization, or wet-lab specificity |
| How much follow-up to plan | `gentle.gene_isoform_assay_study_plan.v1` | Approval to execute, validated biology, or an order |

## 1. Start With Known Inputs

Run from the repository root with an already built `gentle_cli`, Python 3
and `jq`. No build, reference download, Primer3, BLAST, or external agent is
needed. Use a fresh temporary directory so old state cannot masquerade as a
new result:

```bash
GENTLE=/absolute/path/to/gentle_cli
OUT=$(mktemp -d "${TMPDIR:-/tmp}/gentle-assay-followup.XXXXXX")
STATE="$OUT/followup.gentle.json"
WF=docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json

"$GENTLE" --state "$STATE" op \
  '{"LoadFile":{"path":"test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_minus_strand.gb","as_id":"patz1_transcript_assay_demo"}}'
"$GENTLE" --state "$STATE" op \
  '{"SetParameter":{"name":"primer_design_backend","value":"internal"}}'
```

The [existing fixture provenance](../../test_files/fixtures/transcript_assay_panel/patz1/README.md)
describes an invented 240 bp minus-strand locus. Its engine gene feature index
is `1` (`n-2` in the GUI); the preceding source feature is index `0`.
PATZ1-201, PATZ1-202 and PATZ1-203 have mature cDNA lengths 120, 80 and
100 nt respectively. Those lengths are annotation-derived, not RNA measurements.
The GRCh38-like anchor is synthetic and must not be used as proof of human
reference identity.

## 2. Design an RT Pool, Not a PCR Panel

The [small RT request](inputs/transcript_assay_followup_rt_pool.json) names
PATZ1-201 first and PATZ1-202 second. This order is deliberate pool priority,
not measured abundance. It retains five candidates per target, searches the
first 40 bases of each transcript-oriented terminal exon, and uses an exact
22 nt sequence-specific segment plus a fixed 22 nt test adapter.

```bash
"$GENTLE" --state "$STATE" primers design-terminal-exon-rt-pool \
  @docs/tutorial/inputs/transcript_assay_followup_rt_pool.json
"$GENTLE" --state "$STATE" primers export-report \
  transcript_assay_followup_rt_pool "$OUT/rt_pool.json"
jq '{tm_policy, pool_selection_policy, targets:[.targets[] | {
  priority_1based, resolved_transcript_id, strand,
  terminal_exon_source_start_0based, terminal_exon_source_end_0based_exclusive,
  terminal_exon_transcript_start_0based, terminal_exon_transcript_end_0based_exclusive,
  evaluated_candidate_count, selected:[.candidates[] | select(.selected)]
}]}' "$OUT/rt_pool.json"
```

Check the geometry before reading the oligo scores. Both terminal exons are
source interval `[20,60)` in zero-based half-open coordinates, not the
rightmost exon on the genomic page. Their mature-transcript intervals are
`[80,120)` for PATZ1-201 and `[40,80)` for PATZ1-202. Each 40 nt canonical
window permits 19 overlapping 22-mers; five are retained and exactly one is
selected per target. The report has one pairwise selected-pool interaction.
Every full oligo is adapter followed by the reverse complement of its target
segment, written 5-prime to 3-prime. Tm is descriptive, excluded from ranking.

These two transcripts share this terminal exon. Naming two targets therefore
does **not** make their RT primers isoform-specific, and pool interaction
screening is not a genome search. `genomic_specificity` is absent in this
request/report. Do not infer a pass from its absence or silently deduplicate
the targets.

RT primers initiate cDNA synthesis upstream of their binding sites; this
design does not simulate RT yield, full-length synthesis, or inclusion of
sequence downstream of the site. In particular, designing this pool does not
change the `oligo_dt` synthesis model in the existing PATZ1 PCR operations
below. Comparing a gene-specific RT strategy requires a separately reviewed
assay request and appropriate controls, not relabelling this panel.

## 3. Inspect Junction Selection Without Inventing Differential Evidence

Reuse the exact three operation objects from 04.06, changing only output
paths. Their permissive synthetic primer constraints are not recommended
laboratory settings. The RT pool is a separate report in the same project.

```bash
jq --arg p "$OUT/endpoint.json" \
  '.workflow.ops[2] | .DesignTranscriptAssayPanel.path=$p' "$WF" > "$OUT/endpoint.operation.json"
jq --arg p "$OUT/sybr.json" \
  '.workflow.ops[3] | .DesignTranscriptAssayPanel.path=$p' "$WF" > "$OUT/sybr.operation.json"
jq --arg p "$OUT/common.json" \
  '.workflow.ops[4] | .DesignTranscriptAssayPanel.path=$p' "$WF" > "$OUT/common.operation.json"
"$GENTLE" --state "$STATE" primers design-transcript-assay-panel \
  "@$OUT/endpoint.operation.json" --backend internal
"$GENTLE" --state "$STATE" primers design-transcript-assay-panel \
  "@$OUT/sybr.operation.json" --backend internal
"$GENTLE" --state "$STATE" primers design-transcript-assay-panel \
  "@$OUT/common.operation.json" --backend internal
jq '{operation_sha256, junction_evaluations, selected:[.selected_assays[] |
  .primer_pair_summary | {selection_operation_sha256, selection_audit_status,
  selection_audit_method, selection_evidence_observation_count,
  selection_evidence_projection_count, selection_evidence_matched_target_count,
  retained_because_of_differential_junction_evidence,
  retained_despite_zero_marginal_discrimination, selection_evidence}],
  specificity_followups}' "$OUT/sybr.json"
```

Known observations from the committed 04.06 reports:

- The SYBR panel is complete, has three transcript classes and three selected
  primer-only assays. Its required JUC targets PATZ1-202 exon ordinals `1-2`
  at transcript position `40`, with `selected_spanning_assay` status.
- The synthetic JUC effect is `-1.2` for `synthetic_case_vs_control`.
  Nevertheless, its `differential_eligibility` is `not_assessed` and its
  disposition is `incomplete_missing_threshold`. Requiring the geometry did
  not supply a content-bound differential threshold.
- The JUC-bearing pair has one source observation and one transcript
  projection in this fixture. Several projections of one source observation
  would still not be independent measurements. Read the observation,
  projection and matched-target counts separately.
- `selection_audit_status` is `computed`; the method is
  `binary_detection_leave_one_out_v1`. Each pair's
  `selection_operation_sha256` matches the panel's `operation_sha256`.
  A legacy `not_computed_legacy` audit is unknown, not a measured zero.
- Every specificity follow-up still says `genomic_confirmation_status: not_run`.
  The common-region panel derives commonality from annotation; its PSR support
  cannot manufacture commonality. Endpoint band lengths are predictions, not
  gel measurements or quantitative abundance.

For a **new** differential-qualified design, first approve and run a new
`InterpretProbeRegionEvidence` request against the actual projected array
features with `min_abs_logfc`, `threshold_source` and `policy_sha256`.
The shared `arrays interpret-probe-region-evidence` route exposes those fields;
the engine records `interpretation_request_sha256`. Bind that new report to
the subsequent design. Do not edit a saved report to add missing provenance.
The source fixture here has no projected array feature track to regenerate
that positive branch, so this exercise deliberately stops at the missing
threshold result.

An eligible JUC can become a `required_validation_obligation` or
`preferred_differential`, according to its requested priority. A value below
the declared threshold is `ineligible_below_threshold`; missing contrast or
measurement has its own incomplete disposition. A required evidence obligation
can explain retention even when the pair adds zero exclusive binary transcript
distinctions. Inspect the explicit retention fields, not just a redundancy
count. An absolute-effect gate is not a significance test or isoform validation.

## 4. Normalize, Review, Then Plan a Study

The [study request](inputs/transcript_assay_followup_study.json) reuses the
unchanged digest-bound isoform ledger and JUC report. It explicitly chooses
`targeted_junction_validation` for this bounded planning exercise; GENtle
retains its automatic recommendation separately. The 30..120 bp policy is
chosen for this tiny fixture, not as a general qPCR recipe.

```bash
"$GENTLE" --state "$STATE" primers plan-gene-isoform-study \
  @docs/tutorial/inputs/transcript_assay_followup_study.json \
  --normalize-only --normalized-request "$OUT/study.normalized.json"
jq '{policy, policy_sha256, expected_isoform_evidence_sha256,
  junction_evidence, coverage_universe, fallback_submission, profile_override}' \
  "$OUT/study.normalized.json"
```

**Approval 1:** review the complete normalized request, including effective
defaults, hashes, coverage universe, missing evidence and override reason.
In Agent Assistant, approval is a human decision about this exact planning
basis, not permission for the agent to change it. Direct CLI commands do not
prove that a human reviewed anything. For this synthetic exercise, explicitly
choose to proceed with planning before running:

```bash
"$GENTLE" --state "$STATE" primers plan-gene-isoform-study \
  "@$OUT/study.normalized.json" --path "$OUT/study.plan.json" \
  --workflow "$OUT/study.workflow.json"
jq '{request_sha256, recommended_profile, selected_profile, profile_override,
  evidence_summary, decision_factors, planned_operations,
  operation_batch_sha256, approved_workflow_sha256}' "$OUT/study.plan.json"
```

Normalization/planning do not execute primer design or modify the project;
their requested files are outputs. Expect two ordered operations, a
`pan_transcript` common control and a `minimal_discrimination_panel` junction
panel, both strict `require_all`. The planner uses its own default primer
constraints and `preferred` junction priority: these are **not** the required
JUC operation or relaxed constraints replayed in section 3. Success of that
earlier design does not prove feasibility of the planned operations.

Expect `junction_interpretation_incomplete_provenance` among the decision
factors. Setting the study policy threshold to 0.5 does not repair old evidence.
The ledger supplies its own splicing target feature (index `2` here); retain
the emitted operation rather than replacing it with the gene index `1`.

**Approval 2:** separately review the exact ordered `planned_operations[]`,
`operation_batch_sha256` and exact `approved_workflow_sha256`. The bounded
offline exercise ends here. Only after deliberate execution approval is the
continuation:

```bash
"$GENTLE" --state "$STATE" primers execute-gene-isoform-study-workflow \
  "$OUT/study.plan.json" "$OUT/study.workflow.json"
```

Do not reformat, regenerate or alter approved workflow bytes. Changed inputs
or constraints require a fresh plan and approval, not silent best-effort
fallback. The companion test changes the workflow bytes and expects refusal
before any project mutation; it does not claim the untouched plan was executed.

For several genes, the existing `primers compose-gene-isoform-study-workflow-batch`
route accepts ordered `{plan_path, workflow_path}` entries and produces one
reviewable batch. After separate approval, use
`primers execute-gene-isoform-study-workflow-batch`, rather than serially
reusing stale single-state agent approvals. See the [CLI contract](../cli.md).
Do not add a second gene merely to complete this tutorial.

## 5. Optional External Specificity: Plan Is Not Pass

Do **not** submit the synthetic PATZ1 sequences to a real human reference and
interpret the result as PATZ1 validation. This continuation is for a separately
designed real panel with authentic assembly/transcript identity, reviewed
specificity policy, prepared catalogued `genomic_dna` or `transcriptome_cdna`
resources, and an external runner. It is not required by the offline exercise.

The following are shared-shell **templates**, not ready-to-run fixture commands:

```text
primers transcript-assay-specificity-plan PANEL_REPORT_ID --target-genome GENOME_ID --output-dir NEW_OUTPUT_DIR --catalog CATALOG_PATH --cache-dir CACHE_DIR
primers transcript-assay-specificity-finalize HANDOFF.json @EXECUTION_MANIFEST.json --path ACCEPTANCE.json
```

Planning writes an aggregate handoff and an execution-manifest template covering
every selected assay's two primer searches; it does not launch those searches.
Read the returned paths instead of guessing filenames. In a new output
directory, the external runner must execute every exact declared `program`
and `args[]`, then record every command's exit status and output byte length
and SHA-256, including failures. A nonempty TSV is not proof of success; an
empty TSV with successful execution can be a legitimate completed no-hit result.
Never mark the unexecuted manifest template complete.

Always finalize after the runner finishes, even when a command failed:

| Outcome | Meaning and action |
| --- | --- |
| `not_run` | No completed confirmation attached; keep it unknown. |
| `incomplete` | Missing/duplicate/failed/stale/mismatched execution evidence; no partial acceptance is attached. Repair execution/provenance, not biology by assumption. |
| `not_assessed` | Execution completed but intended-target geometry cannot be assessed; do not accept. |
| `specificity_fail` | Complete evidence failed the biological policy; review/redesign, not a scheduler retry labelled success. |
| `pass` | Complete all-assay acceptance for that declared target and policy only; not experimental validation. |

Finalization binds the current panel, assays/primer sequences, policy,
reference content, handoff and execution outputs. Genomic-DNA and whole-cDNA
assessments remain distinct. The optional RT-pool genomic screen is also a
different operation: it checks variable RT-primer segments and does not
replace transcript-panel finalization. No vendor submission follows automatically.

## 6. Continue in the GUI or Inner Agent

Open the saved project through `File -> Open Project...` in a separate GUI
session; do not concurrently overwrite its state from CLI. The current code
supports these controls, but this walkthrough is not a live GUI sign-off:

1. Open the PATZ1 sequence, select its gene feature and use `Open Splicing Window`.
2. In Splicing Expert, use `Design all-transcript panel` to open PCR Designer's
   `Transcript panels` mode. Review the source group and transcript matrix.
3. Switch to `RT primer pool`. With the Splicing Expert group open, select
   PATZ1-201 then PATZ1-202 and use `Add current transcript` once each in that
   order. Start with an empty target list; check the feature label is `n-2`.
4. Enter the RT request's adapter, variable length `22`, terminal-exon window
   `40`, retain-per-target `5`, and report ID. `Design RT primer pool` runs
   the shared operation; `Open saved report`, `Export report JSON...` and
   `Copy selected oligos TSV` inspect/export it, not order it.
5. Use the GUI Shell for exact JSON-based continuation. For example,
   `primers show-transcript-assay-panel patz1_sybr_juc_panel` inspects the
   persisted panel. The same `primers plan-gene-isoform-study` and specificity
   commands above are Shell routes; this guide does not claim dedicated study
   approval or aggregate-finalization buttons. Replace `$OUT` with the actual
   absolute path inside the GUI Shell, which is not Bash.

After catalog generation, Help or Agent Assistant can open this teaching text
with `ui open tutorial-guide transcript_assay_followup_gui_cli`. It opens a
guide, not a populated tutorial project. An interest prompt such as "Help me
review terminal-exon RT pools, differential junction selection and transcript
assay specificity before planning a gene isoform study" matches its title and
catalog notes through the existing interest-guided retrieval. Retrieval score
is not scientific confidence. Keep provider configuration, model suggestions,
execution confirmation and human scientific approval separate.

## Reproducibility and Review

```bash
python3 scripts/test_transcript_assay_followup_tutorial.py
GENTLE_TUTORIAL_BIN_DIR=/absolute/path/to/binary/directory \
  python3 scripts/test_transcript_assay_followup_tutorial.py
```

Without the environment variable, source/fixture checks run and real-binary
tests explicitly skip. With it, tests invoke only `gentle_cli`, use the
internal backend and temporary state, compare direct/shared-shell RT reports,
replay the three existing panels, inspect/export them, normalize/plan a study
and reject changed workflow bytes. They do not build GENtle, launch BLAST,
execute the approved study workflow, contact a model, or submit an order.

Retain the request, exported reports, exact design operations, normalized
study basis, plan/workflow bytes, and any later specificity manifest/acceptance
with your review notes. Record unresolved evidence and failed checks as such.
Human GUI acceptance, positive threshold-qualified replay, actual external
specificity and wet-lab controls remain separate reviews, not implied passes.
Use Help's `Copy Feedback Context` when reporting a tutorial mismatch.
