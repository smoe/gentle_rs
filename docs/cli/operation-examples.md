# Operation JSON Examples

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

### Example operations


## Load and Create Sequences

Load a file:

```json
{"LoadFile":{"path":"test_files/pGEX-3X.gb","as_id":"pgex"}}
```

Create a persistent project sequence from inline bases:

```json
{"CreateSequenceFromText":{"sequence_text":"ATGGAATTCGGGCCCTAA","output_id":"paper_candidate","name":"Candidate from publication","circular":false}}
```

Digest:

```json
{"Digest":{"input":"pgex","enzymes":["BamHI","EcoRI"],"output_prefix":"pgex_frag"}}
```

Ligation:

```json
{"Ligation":{"inputs":["pgex_frag_1","pgex_frag_2"],"circularize_if_possible":true,"output_id":"re_ligated","protocol":"Sticky","output_prefix":"lig","unique":true}}
```

Merge multiple containers/pools into one candidate pool namespace:

```json
{"MergeContainers":{"inputs":["tubeA_1","tubeB_3","tubeC_2"],"output_prefix":"merged_pool"}}
```

Note:

- Ligation is protocol-driven and requires explicit `protocol`
  (`"Sticky"` or `"Blunt"`).


## PCR and Primer Mutagenesis

PCR:

```json
{"Pcr":{"template":"pgex","forward_primer":"ATGGCT","reverse_primer":"CGTACC","output_id":"amplicon1","unique":true}}
```

Advanced PCR (tail insertion / mismatch-tolerant annealing):

```json
{"PcrAdvanced":{"template":"pgex","forward_primer":{"sequence":"GGATCCATGGCT","anneal_len":6,"max_mismatches":1,"require_3prime_exact_bases":4},"reverse_primer":{"sequence":"CGTACC","anneal_len":6,"max_mismatches":0,"require_3prime_exact_bases":4},"output_id":"amplicon_site_added","unique":true}}
```

Advanced PCR with degenerate/randomized primer library:

```json
{"PcrAdvanced":{"template":"pgex","forward_primer":{"sequence":"ATNAAA","anneal_len":6,"max_mismatches":1,"require_3prime_exact_bases":3,"library_mode":"Sample","max_variants":10,"sample_seed":42},"reverse_primer":{"sequence":"AAACCC","anneal_len":6,"max_mismatches":0,"require_3prime_exact_bases":4},"unique":false}}
```

Restriction-site cloning PCR handoff from a saved primer-pair report:

```json
{"PrepareRestrictionCloningPcrHandoff":{"template":"insert_tpl","primer_report_id":"insert_tpl_pairs_v1","pair_index":0,"destination_vector_seq_id":"pgl3_mcs","mode":"directed_pair","forward_enzyme":"EcoRI","reverse_enzyme":"HindIII","forward_leader_5prime":"GC","reverse_leader_5prime":"AT"}}
```

PCR mutagenesis (explicit SNP intent):

```json
{"PcrMutagenesis":{"template":"pgex","forward_primer":{"sequence":"ATCAAA","anneal_len":6,"max_mismatches":1,"require_3prime_exact_bases":3},"reverse_primer":{"sequence":"AAACCC","anneal_len":6,"max_mismatches":0,"require_3prime_exact_bases":4},"mutations":[{"zero_based_position":2,"reference":"G","alternate":"C"}],"output_id":"snp_product","unique":true,"require_all_mutations":true}}
```

Extract region (`from` inclusive, `to` exclusive):

```json
{"ExtractRegion":{"input":"pgex","from":100,"to":900,"output_id":"insert_candidate"}}
```

Derive splicing references from a sequence window (DNA window + mRNA isoforms +
exon-reference):

```json
{"DeriveSplicingReferences":{"seq_id":"tp73_window","span_start_0based":1200,"span_end_0based":2600,"seed_feature_id":null,"scope":"target_group_target_strand","output_prefix":"tp73_refs"}}
```


## Alignment and Region Extraction

Pairwise sequence alignment (global/local) with structured result payload:

```json
{"AlignSequences":{"query_seq_id":"tp73_refs_mrna_NM_001204186","target_seq_id":"tp73_refs_exon_reference","mode":"local","match_score":2,"mismatch_score":-3,"gap_open":-5,"gap_extend":-1}}
```

The same operation can run statelessly over inline ASCII operands:

```json
{"AlignSequences":{"query":{"kind":"inline_sequence","sequence_text":"TTTACGTAA","id_hint":"inline_query","span_start_0based":3,"span_end_0based_exclusive":7},"target":{"kind":"inline_sequence","sequence_text":"GGGACGTCCC","id_hint":"inline_target","span_start_0based":3,"span_end_0based_exclusive":7},"mode":"global"}}
```

Shared shell:

```bash
gentle_cli shell 'align compute --query-sequence-text TTTACGTAA --query-id-hint inline_query --query-range 3..7 --target-sequence-text GGGACGTCCC --target-id-hint inline_target --target-range 3..7 --mode global'
```

Extract anchored region with flexible 5' boundary and fixed anchor side:

```json
{"ExtractAnchoredRegion":{"input":"pgex","anchor":{"FeatureBoundary":{"feature_kind":"CDS","feature_label":null,"boundary":"Start","occurrence":0}},"direction":"Upstream","target_length_bp":500,"length_tolerance_bp":100,"required_re_sites":["EcoRI"],"required_tf_motifs":["TATAAA"],"forward_primer":"GAATTC","reverse_primer":"CGTACC","output_prefix":"promoter","unique":false,"max_candidates":20}}
```

Note:

- `ExtractAnchoredRegion.anchor` is a local sequence anchor (in-sequence
  coordinate resolver), not genome-assembly provenance anchoring.

`required_tf_motifs` accepts either IUPAC motif strings or motif IDs/names from
the local JASPAR snapshot.


## TFBS and Promoter Evidence

Annotate TFBS features using log-likelihood ratio thresholds:

```json
{"AnnotateTfbs":{"seq_id":"pgex","motifs":["MA0139.1","SP1","TATAAA"],"min_llr_bits":0.0,"min_llr_quantile":0.95,"per_tf_thresholds":[{"tf":"SP1","min_llr_bits":-1.0,"min_llr_quantile":0.80}],"clear_existing":true,"max_hits":500}}
```

`AnnotateTfbs` writes generated TFBS features with score qualifiers:

- `llr_bits`: absolute log-likelihood ratio score (base 2)
- `llr_quantile`: empirical score quantile in the scanned sequence region (both strands)
- `true_log_odds_bits`: smoothed true log-odds score (base 2)
- `true_log_odds_quantile`: empirical quantile of `true_log_odds_bits` in the scanned region
- compatibility aliases are also written: `log_odds_ratio_bits`, `log_odds_ratio_quantile`

Motif-selection shortcuts:

- `motifs: ["ALL"]` or `motifs: ["*"]` scans all motifs from the local JASPAR registry.
- `max_hits` controls safety capping of generated TFBS features:
  - omitted: default cap of `500` hits
  - `0`: unlimited (no cap)
  - `N > 0`: stop after `N` accepted hits

Summarize grouped TFBS occupancy in one focus window relative to a wider
context through the shared shell:

```bash
cargo run --quiet --bin gentle_cli -- --state STATE.json shell \
  'features tfbs-summary vkorc1_rs9923231_context --focus 2900..3100 --context 0..6001 --limit 25'
```

The JSON payload groups TFBS features by factor name (preferring
`bound_moiety`, then `standard_name`, then `tf_id`) and reports:

- `focus_occurrences`
- `context_occurrences`
- `outside_focus_occurrences`
- `focus_density_per_kb`
- `context_density_per_kb`
- `outside_focus_density_per_kb`
- `focus_vs_context_density_ratio`
- `focus_vs_outside_density_ratio`

Notes:

- `--context` is optional; when omitted, the full sequence is used as the wider
  comparison window.
- `--min-focus-count` defaults to `1`, so the command naturally lists TFs that
  are actually present in the requested focus window.

The same summary is also available as a first-class JSON operation:

```json
{"SummarizeTfbsRegion":{"seq_id":"vkorc1_rs9923231_context","focus_start_0based":2900,"focus_end_0based_exclusive":3100,"context_start_0based":0,"context_end_0based_exclusive":6001,"min_focus_occurrences":1,"min_context_occurrences":0,"limit":25}}
```


## Promoter-SNP Follow-Up

Promoter-SNP follow-up commands for the VKORC1 / `rs9923231` handoff:

```bash
cargo run --quiet --bin gentle_cli -- \
  variant annotate-promoters vkorc1_rs9923231_context \
  --gene-label VKORC1 \
  --upstream-bp 1000 \
  --downstream-bp 200

cargo run --quiet --bin gentle_cli -- \
  variant promoter-context vkorc1_rs9923231_context \
  --variant rs9923231 \
  --gene-label VKORC1 \
  --path docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/variant_promoter_context.json

cargo run --quiet --bin gentle_cli -- \
  variant reporter-fragments vkorc1_rs9923231_context \
  --variant rs9923231 \
  --gene-label VKORC1 \
  --retain-downstream-from-tss-bp 200 \
  --retain-upstream-beyond-variant-bp 500 \
  --path docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json

cargo run --quiet --bin gentle_cli -- \
  variant materialize-allele vkorc1_rs9923231_promoter_fragment \
  --variant rs9923231 \
  --allele alternate \
  --alternate-base T \
  --output-id vkorc1_rs9923231_promoter_alternate
```

These direct shell/CLI routes keep the promoter classification, fragment
suggestion, and allele materialization on the same shared engine contracts used
by operation JSON, workflows, and the updated tutorial path.

`variant annotate-promoters --collapse tss-cluster` uses a 50 bp tolerance;
append a tolerance such as `tss-cluster:25` to make it explicit. A TSS class
requires the same gene and strand plus overlapping transcript-oriented first
exons, and retains its complete transcript membership and grouping reason.

The dedicated `variant reporter-fragments` route remains the compatibility
surface for variant-anchored designs in this engine-first slice. Local
motif-anchored or explicit-interval planning is available through the same
`SuggestPromoterReporterFragments` operation contract, for example:

```json
{"SuggestPromoterReporterFragments":{"input":"tgfb1_promoter_context","gene_label":"TGFB1","retain_downstream_from_tss_bp":200,"retain_upstream_beyond_variant_bp":500,"max_candidates":5,"fragment_policy":{"anchor":{"kind":"motif_hit","label_or_id":"TP73","occurrence":1},"collapse_mode":{"tss_cluster":{"tolerance_bp":50}},"promoter_upstream_baseline_bp":500,"anchor_flank_bp":150,"max_fragment_length_bp":5000}}}
```

Motif selection uses only local TFBS annotations. Its report is explicitly
sequence-motif evidence, not evidence of occupancy or functional regulation;
overlength geometries remain visible as typed rejections.

Compare transcript-aware promoter/reporter geometries without creating a
construct through the shared `promoters` command family:

```bash
cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters compare-architectures \
  @docs/examples/assets/serpine1_promoter_reporter_architecture_request.json \
  --path analysis/serpine1_architecture_comparison.json \
  --svg-path analysis/serpine1_architecture_comparison.svg
```

The loaded state must contain the annotated locus named by `seq_id`. The
example request selects three pinned Ensembl-116 SERPINE1 transcript models
representing the main and rare-upstream TSS classes and compares
TSS-proximal, spliced-leader, and genomic-leader models. It deliberately
supplies neither initiation evidence nor CUT&RUN sources, so those layers stay
explicitly unevaluated rather than being inferred from annotation. Add
`endogenous_atg_retained_fusion` only together with a matching
`fusion_requests` row that states retained endogenous CDS length, fusion mode,
vector-ATG removal, and junction sequence.

`--path` and `--svg-path` are independent. Omitting either keeps that artifact
in memory while returning the same portable report. MCP callers use the same
`ComparePromoterReporterArchitectures` operation through the confirmed generic
`op` route; there is no MCP-only analysis or scoring model. The corresponding
request and report schemas are
`gentle.promoter_reporter_architecture_comparison_request.v1` and
`gentle.promoter_reporter_architecture_comparison.v1`; their field-level
contract is documented in [protocol.md](../protocol.md).

Inspect and then build a complete promoter-reporter panel proposal through the
distinct `promoters` command family (`reporters` remains the reporter-gene
catalog family):

```bash
cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters compose-study @study_request.json \
  --path study_report.json
```

For read-only, variable-length boundary proposals from a composed locus, use
`promoters fragment-candidates REQUEST_JSON_OR_@FILE --path candidates.json`.
It supports explained adjustments at either end, a soft compact-length target,
fixed-envelope signal comparisons, explicit provenance-bound called-peak BED
sources and optional validated-vector MCS checks. Coverage alone never seeds.
See [Evidence-Guided Reporter Boundaries](../reporter_fragment_selection.md).
This is upstream of exact-ROI planning; it does not save regions or approve constructs.
For fixed-window similarity geometry instead, use
`gentle_cli shell 'op @geometry_op.json'` with `ComputeTssWindowGeometry`.
The [operation example](../reporter_fragment_selection.md#fixed-window-similarity-geometry)
returns bounded, strand-aware windows/unions/intersections without a project;
the preparation script requires `--gentle /path/to/gentle_cli`, with no Python fallback.

`promoters compose-study` is the multi-gene front-half route. Its
`gentle.regulatory_reporter_study_request.v1` input accepts exactly one of a
gene-set `source` or an existing `gentle.gene_set_resolution.v1`, a prepared
genome/build, a TSS policy, an optional exact per-member evidence-kind gate,
independent evidence rows, exact vector identity, native-versus-engineered
mutation policy, and an output directory. The default TSS policy selects the
outermost annotated 5' transcript; `explicit_per_member` requires one transcript
id for every resolved member. The default fragment anchor is the selected TSS.

The operation extracts one annotated promoter locus per resolved gene into the
project and writes content-hashed `gentle.promoter_reporter_candidates.v1`
files plus one `gentle.promoter_reporter_panel_request.v1`. The returned
`gentle.regulatory_reporter_study.v1` embeds the gene-set/promoter cohort and a
fresh context-bound panel-readiness report. It does not create constructs or
infer occupancy, causality, promoter sufficiency, or reporter activity. A
successful composition can still report `blocked`; correct the named readiness
checks before proceeding. Supplying a reviewed co-regulated resolution adds
its retrieval provenance as `perturbation_response`, separate from motif or
occupancy evidence.

For manually assembled candidate files, or after reviewing the composer's
panel request, use the existing readiness and planning steps:

```bash
cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters panel-readiness request @panel_request.json \
  --path panel_readiness.json

cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters panel-plan @panel_request.json --path panel_proposal.json
```

When resuming from the saved proposal, ask GENtle to check its currency before
materialization. Omit `--approve` to stop at `ready_for_approval`:

```bash
DIGEST=$(jq -r .proposal_digest panel_proposal.json)
cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters panel-readiness proposal @panel_proposal.json \
  --approve "$DIGEST" --path panel_approval_readiness.json
```

The first form reports `ready_to_plan` only after the exact vector, selected
candidate/TSS geometry, supplied evidence coordinates, and complete detached
panel simulation pass. The proposal form additionally recomputes the canonical
proposal against current project and file contents; it reports
`ready_for_approval` without a digest and `ready_to_materialize` only when the
supplied digest matches exactly. Evidence may be `not_supplied` without making
the mechanics invalid, and evidence compatibility is not a claim of occupancy
or promoter function. `evidence_kinds_by_member_id` exposes the normalized
per-member evidence vocabulary without requiring callers to parse check prose.
The wrapper/report schemas are
`gentle.promoter_reporter_panel_readiness_request.v1` and
`gentle.promoter_reporter_panel_readiness.v1`.

The request names one loaded, exactly catalog-validated vector and one or more
saved `gentle.promoter_reporter_candidates.v1` members. Use
`mutation_policy: "native_only_v1"` for an unchanged regulatory-fragment panel.
Use `p53_family_core_disruption_v1` only for explicitly requested p53-family
motif controls; that policy retains the existing wild-type/mutant pair and
requires a motif-hit anchor. A member may override the request policy, allowing
selected controls without mutating every fragment. Omitted or `unspecified`
effective policies fail before sequence editing. Optional study-specific
response, occupancy, or other selection support belongs in each member's
`evidence[]` ledger with a stable id, kind, source, summary, provenance, and
interpretation tags. Study-wide interpretation limits belong in request-bound
`scientific_caveats`. Both are preserved in the proposal rather than supplied
by the generic operation.

Planning is read-only: GENtle simulates fragment extraction, optional
policy-selected controls, shared directional restriction selection or explicit
Gibson fallback, cloning primers, circular products, feature transfer, and
GenBank/SVG/manifest paths in a detached engine. Candidate anchor provenance and
overlapping evidence rows, including motif scores and track provenance, are
bound into the proposal without turning motif evidence into occupancy evidence.
Each product must preserve all catalog-required vector annotations outside the
replaced MCS and carries the assembly junction-validation basis plus a final scan
from the shared restriction-enzyme machinery. Any site count exceeding the
combined vector-plus-insert count is flagged for review. The proposal binds
those results and every source/state hash under `proposal_digest`.

Core members retain the ranked candidate's fragment geometry. An extended
member must explicitly bind one annotated transcript and the canonical CDS
start-codon boundary:

```json
{
  "candidate_set_path": "tgfb1_candidates.json",
  "candidate_id": "tgfb1_core",
  "fragment_role": "extended",
  "extended_boundary": {
    "kind": "canonical_cds_start_exclusive",
    "transcript_id": "ENST00000221930"
  },
  "label": "TGFB1 extended"
}
```

GENtle derives the strand-aware genomic endpoint from that transcript's CDS
annotation. `canonical_cds_start_exclusive` retains the complete annotated 5'
UTR but excludes every genomic base supplying the start codon; use
`canonical_cds_start_codon` only when the source ATG itself belongs in the
insert. The proposal's `extended_boundary_audit` exposes the exact codon source
ranges, spliced 5'-UTR exon ranges, and typed upstream-ATG, uORF, and
5'-UTR-intron warnings. Missing, duplicate, strand-mismatched, or CDS-free
annotations fail instead of falling back to an inferred genomic ATG.

After reviewing the proposal, materialize only that exact digest:

```bash
DIGEST=$(jq -r .proposal_digest panel_proposal.json)
cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters panel-materialize @panel_proposal.json --approve "$DIGEST"
```

Materialization recomputes the proposal before writing. A changed project,
candidate file, helper-vector catalog, vector sequence, proposal field, output
path, or approval digest fails closed; existing artifacts are not overwritten.
Successful runs commit the simulated sequences and write one GenBank and SVG
per circular product plus a construct/primer TSV manifest. Planning does not
order oligos or execute an experiment, and motif/PWM changes are not occupancy
or functional-ablation evidence. See the executable offline synthetic tutorial
`08-09_promoter_reporter_panel_planning_offline`.

Plan an exact multi-fragment regulatory comparison from persisted genomic
regions through the separate read-only contract:

```bash
cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters regulatory-panel-plan @regulatory_panel_request.json \
  --path regulatory_panel_plan.json

cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters regulatory-panel-render @regulatory_panel_plan.json \
  --path regulatory_panel_plan.svg
```

The request binds each selected region, its containing set digest, current
sequence projection, reference release, exact vector identity, declared order,
orientation, spacer DNA, and requested comparison. Candidate A and a minimal
promoter are required; partner B and a reference control are optional. Partner,
order, orientation, and spacing questions require an exact partner binding and
matching declared geometry. The planner evaluates five local sequence/context
lanes independently. An evidence binding can supply `report_path` to an exact
typed locus/reporter JSON, `report_sha256` over its bytes, `report_id` equal to
the nested locus panel ID, and an optional feature/track/lane `row_id`.
Current sequence, anchor, assembly and gene-annotation release are checked.
Opaque citations and unavailable content still cannot become an evaluated lane.

The SVG command accepts only an unchanged plan with a valid embedded digest.
Both commands are non-mutating. Planning labels describe testable comparisons rather than
regulatory sufficiency, activity, or causation. See the executable offline
walkthrough `08-11_regulatory_fragment_panel_planning_offline` for a bounded
A/B panel with explicit order, orientation, and spacing variants.

Prepare exact designed products, inspect the proposal, and only then approve
that exact object. These aliases return the shared operation-result envelope:

```bash
gentle_cli --state STATE.json promoters regulatory-products-plan \
  @regulatory_panel_plan.json --output-prefix design > product_result.json
jq '.result.regulatory_fragment_materialization_proposal' product_result.json > products.json
# Review all product sequences, instances, spacers and omitted annotations first.
DIGEST=$(jq -r .proposal_digest products.json)
gentle_cli --state STATE.json promoters regulatory-products-materialize \
  @products.json --approve "$DIGEST"
```

The second command is an explicit mutation. It creates all exact products or
none, never overwrites IDs, and records a receipt with final-product audit
`not_evaluated`. This is sequence design by vector-context replacement, not
restriction/Gibson simulation or biological validation. Source or annotation
changes require a fresh proposal and review. Use normal sequence export for
GenBank/FASTA files after inspecting the products.

The Promoter design GUI's `Regulatory-fragment panel > Exact design products`
controls consume these same two operations. Product approval is separate from
panel review, and the copied GUI proposal JSON can be inspected through these
shell routes without changing its digest.

Select one candidate in-silico (explicit provenance step):

```json
{"SelectCandidate":{"input":"pgex_frag_1","criterion":"band_size_range:450-550bp","output_id":"picked_band"}}
```

Filter candidates by molecular-weight proxy (bp length) with tolerance and uniqueness:

```json
{"FilterByMolecularWeight":{"inputs":["pgex_frag_1","pgex_frag_2","pgex_frag_3"],"min_bp":450,"max_bp":550,"error":0.10,"unique":true,"output_prefix":"mw_pick"}}
```

Apply practical design-constraint filters (GC bounds, homopolymer cap, U6 `TTTT`
avoidance, forbidden motifs):

```json
{"FilterByDesignConstraints":{"inputs":["g1","g2","g3"],"gc_min":0.30,"gc_max":0.70,"max_homopolymer_run":4,"reject_ambiguous_bases":true,"avoid_u6_terminator_tttt":true,"forbidden_motifs":["GAATTC"],"unique":false,"output_prefix":"design_pick"}}
```

Notes:

- `reject_ambiguous_bases` defaults to `true` if omitted.
- `avoid_u6_terminator_tttt` defaults to `true` if omitted.
- `forbidden_motifs` supports IUPAC and is checked on both strands.


## Candidate Sets

Generate a candidate set from windows on one sequence:

```json
{"GenerateCandidateSet":{"set_name":"sgrna_windows","seq_id":"grch38_tp53","length_bp":20,"step_bp":1,"feature_kinds":["gene"],"feature_label_regex":"^TP53$","max_distance_bp":5000,"limit":10000}}
```

Generate a candidate set only between two local anchors:

```json
{"GenerateCandidateSetBetweenAnchors":{"set_name":"sgrna_between","seq_id":"grch38_tp53","anchor_a":{"FeatureBoundary":{"feature_kind":"gene","feature_label":"TP53","boundary":"Start","occurrence":0}},"anchor_b":{"FeatureBoundary":{"feature_kind":"gene","feature_label":"TP53","boundary":"End","occurrence":0}},"length_bp":20,"step_bp":1,"limit":10000}}
```

Add a derived metric expression:

```json
{"ScoreCandidateSetExpression":{"set_name":"sgrna_windows","metric":"gc_balance","expression":"100 * (gc_fraction - at_fraction)"}}
```

Add distance-to-feature score:

```json
{"ScoreCandidateSetDistance":{"set_name":"sgrna_windows","metric":"distance_to_cds_bp","feature_kinds":["CDS"],"feature_label_regex":null}}
```

Add a weighted multi-objective score:

```json
{"ScoreCandidateSetWeightedObjective":{"set_name":"sgrna_windows","metric":"priority_score","objectives":[{"metric":"gc_fraction","weight":0.7,"direction":"maximize"},{"metric":"distance_to_cds_bp","weight":0.3,"direction":"minimize"}],"normalize_metrics":true}}
```

Select top-k by one metric:

```json
{"TopKCandidateSet":{"input_set":"sgrna_windows","output_set":"sgrna_top20","metric":"priority_score","k":20,"direction":"maximize","tie_break":"seq_start_end"}}
```

Compute Pareto frontier for multiple objectives:

```json
{"ParetoFrontierCandidateSet":{"input_set":"sgrna_windows","output_set":"sgrna_frontier","objectives":[{"metric":"gc_fraction","direction":"maximize"},{"metric":"distance_to_cds_bp","direction":"minimize"}],"max_candidates":100,"tie_break":"seq_start_end"}}
```

Filter by absolute value and quantiles into a new explicit set:

```json
{"FilterCandidateSet":{"input_set":"sgrna_windows","output_set":"sgrna_windows_top","metric":"gc_balance","min":-20.0,"max":20.0,"min_quantile":0.10,"max_quantile":0.90}}
```

Intersect/union/subtract explicit candidate sets:

```json
{"CandidateSetOp":{"op":"intersect","left_set":"sgrna_windows_top","right_set":"near_tss","output_set":"sgrna_final"}}
{"CandidateSetOp":{"op":"union","left_set":"set_a","right_set":"set_b","output_set":"set_union"}}
{"CandidateSetOp":{"op":"subtract","left_set":"set_a","right_set":"set_b","output_set":"set_a_minus_b"}}
```

Delete a candidate set:

```json
{"DeleteCandidateSet":{"set_name":"sgrna_windows"}}
```

Create transformed or branched candidates:

```json
{"Reverse":{"input":"pgex","output_id":"pgex_rev"}}
{"Complement":{"input":"pgex","output_id":"pgex_comp"}}
{"ReverseComplement":{"input":"pgex","output_id":"pgex_revcomp"}}
{"Branch":{"input":"pgex","output_id":"pgex_branch"}}
```


## Display and Engine Parameters

Set visibility of a GUI-equivalent display toggle (example: features):

```json
{"SetDisplayVisibility":{"target":"Features","visible":false}}
```

Set linear-map viewport (zoom/pan state):

```json
{"SetLinearViewport":{"start_bp":1000,"span_bp":5000}}
```

Set an in-silico engine parameter (example: cap fragment/product combinatorics):

```json
{"SetParameter":{"name":"max_fragments_per_container","value":80000}}
```

Enable strict anchor-verification policy for genome-anchor extension:

```json
{"SetParameter":{"name":"require_verified_genome_anchor_for_extension","value":true}}
```

Set prepared-genome fallback policy for anchor extension/verification:

```json
{"SetParameter":{"name":"genome_anchor_prepared_fallback_policy","value":"single_compatible"}}
{"SetParameter":{"name":"genome_anchor_prepared_fallback_policy","value":"always_explicit"}}
{"SetParameter":{"name":"genome_anchor_prepared_fallback_policy","value":"off"}}
```

Set feature-details font size used in the feature tree/details panel (valid range `8.0..24.0`):

```json
{"SetParameter":{"name":"feature_details_font_size","value":10.5}}
```

Set sequence text-panel maximum length (`200000` default, `0` means unlimited):

```json
{"SetParameter":{"name":"sequence_panel_max_text_length_bp","value":200000}}
```

Set regulatory-overlay max linear view span threshold (`50000` recommended for anchored genome maps):

```json
{"SetParameter":{"name":"regulatory_feature_max_view_span_bp","value":50000}}
```

Set primer-design backend controls (internal baseline + optional Primer3 backend):

```json
{"SetParameter":{"name":"primer_design_backend","value":"auto"}}
{"SetParameter":{"name":"primer3_executable","value":"primer3_core"}}
```

Set adaptive linear DNA-letter routing parameters (shared GUI/runtime semantics):

```json
{"SetParameter":{"name":"linear_sequence_letter_layout_mode","value":"auto"}}
{"SetParameter":{"name":"linear_show_sequence_bases","value":true}}
{"SetParameter":{"name":"linear_sequence_helical_letters_enabled","value":true}}
{"SetParameter":{"name":"linear_sequence_helical_phase_offset_bp","value":3}}
{"SetParameter":{"name":"linear_show_reverse_strand_bases","value":true}}
{"SetParameter":{"name":"linear_helical_parallel_strands","value":true}}
{"SetParameter":{"name":"reverse_strand_visual_opacity","value":0.55}}
```

Supported `linear_sequence_letter_layout_mode` values:

- `auto`, `adaptive`, `auto_adaptive`
- `standard`, `standard_linear`
- `helical`, `continuous_helical`
- `condensed_10_row`, `condensed-10-row`, `condensed`

Compatibility notes:

- `linear_show_sequence_bases` is the master visibility switch for DNA letters
  in the linear map. Setting it to `false` hides base letters even when an
  adaptive or forced letter mode would otherwise draw them; feature tracks and
  the DNA backbone remain available. Aliases accepted by `SetParameter` are
  `linear_show_sequence_letters`, `show_linear_sequence_bases`, and
  `show_linear_sequence_letters`.
- Legacy fixed-threshold parameters are still accepted for compatibility but are
  deprecated no-ops under adaptive routing:
  - `linear_sequence_base_text_max_view_span_bp`
  - `linear_sequence_helical_max_view_span_bp`
  - `linear_sequence_condensed_max_view_span_bp`
- `linear_sequence_helical_letters_enabled` applies to auto mode only
  (`linear_sequence_letter_layout_mode = auto*`).
- reverse-strand visibility aliases map to the same control:
  - `linear_show_double_strand_bases`
  - `linear_show_reverse_strand_bases`
- `linear_helical_parallel_strands=true` keeps forward/reverse helical slant
  in parallel; `false` uses mirrored slant.
- `reverse_strand_visual_opacity` (aliases:
  `linear_reverse_strand_visual_opacity`,
  `linear_reverse_strand_letter_opacity`) controls reverse-strand letter
  emphasis across linear map and sequence panel (`0.2..1.0`).


## TFBS, Restriction and VCF Display

Set TFBS display filtering parameters shared by GUI and SVG export:

```json
{"SetParameter":{"name":"show_tfbs","value":true}}
{"SetParameter":{"name":"tfbs_display_use_llr_bits","value":true}}
{"SetParameter":{"name":"tfbs_display_min_llr_bits","value":0.0}}
{"SetParameter":{"name":"tfbs_display_use_llr_quantile","value":true}}
{"SetParameter":{"name":"tfbs_display_min_llr_quantile","value":0.95}}
{"SetParameter":{"name":"tfbs_display_use_true_log_odds_bits","value":false}}
{"SetParameter":{"name":"tfbs_display_min_true_log_odds_bits","value":0.0}}
{"SetParameter":{"name":"tfbs_display_use_true_log_odds_quantile","value":false}}
{"SetParameter":{"name":"tfbs_display_min_true_log_odds_quantile","value":0.95}}
```

Supported TFBS display parameter names:

- `show_tfbs`
- `tfbs_display_use_llr_bits`
- `tfbs_display_min_llr_bits`
- `tfbs_display_use_llr_quantile`
- `tfbs_display_min_llr_quantile` (quantile in range `0.0..1.0`)
- `tfbs_display_use_true_log_odds_bits`
- `tfbs_display_min_true_log_odds_bits`
- `tfbs_display_use_true_log_odds_quantile`
- `tfbs_display_min_true_log_odds_quantile` (quantile in range `0.0..1.0`)

Set restriction-enzyme display parameters shared by GUI and SVG export:

```json
{"SetParameter":{"name":"show_restriction_enzymes","value":true}}
{"SetParameter":{"name":"restriction_enzyme_display_mode","value":"preferred_only"}}
{"SetParameter":{"name":"preferred_restriction_enzymes","value":["EcoRI","BamHI"]}}
```

Supported restriction display parameter names:

- `show_restriction_enzymes`
- `show_restriction_enzyme_sites` (bool alias)
- `restriction_enzyme_display_mode`
- `restriction_display_mode` (string alias)
- `preferred_restriction_enzymes` (CSV string or string array)
- `preferred_restriction_enzymes_csv` (CSV string alias)
- `restriction_preferred_enzymes` (CSV string or string array alias)

Supported `restriction_enzyme_display_mode` values:

- `preferred_only`
- `preferred_and_unique`
- `unique_only`
- `all_in_view`

Set VCF display filtering parameters shared by GUI and SVG export:

```json
{"SetParameter":{"name":"vcf_display_pass_only","value":true}}
{"SetParameter":{"name":"vcf_display_use_min_qual","value":true}}
{"SetParameter":{"name":"vcf_display_min_qual","value":30.0}}
{"SetParameter":{"name":"vcf_display_required_info_keys","value":["AF","DP"]}}
```

Supported VCF display parameter names:

- `vcf_display_show_snp`
- `vcf_display_show_ins`
- `vcf_display_show_del`
- `vcf_display_show_sv`
- `vcf_display_show_other`
- `vcf_display_pass_only`
- `vcf_display_use_min_qual`
- `vcf_display_min_qual`
- `vcf_display_use_max_qual`
- `vcf_display_max_qual`
- `vcf_display_required_info_keys` (string CSV or string array)


## Export and Rendering

Export selected pool members (engine operation):

```json
{"ExportPool":{"inputs":["frag_1","frag_2"],"path":"digest.pool.gentle.json","pool_id":"digest_1","human_id":"BamHI+EcoRI digest"}}
```

Export built-in DNA ladder catalog (engine operation):

```json
{"ExportDnaLadders":{"path":"dna_ladders.snapshot.json","name_filter":null}}
{"ExportDnaLadders":{"path":"dna_ladders.neb.json","name_filter":"NEB"}}
{"ExportRnaLadders":{"path":"rna_ladders.snapshot.json","name_filter":null}}
{"ExportRnaLadders":{"path":"rna_ladders.neb.json","name_filter":"NEB"}}
```

Render sequence SVG (engine operation):

```json
{"RenderSequenceSvg":{"seq_id":"pgex","mode":"Linear","path":"pgex.linear.svg"}}
{"RenderSequenceSvg":{"seq_id":"pgex","mode":"Circular","path":"pgex.circular.svg"}}
```

Example: linear DNA-window export with explicit restriction display state:

```json
{"SetParameter":{"name":"show_restriction_enzymes","value":true}}
{"SetParameter":{"name":"restriction_enzyme_display_mode","value":"preferred_only"}}
{"SetParameter":{"name":"preferred_restriction_enzymes","value":["EcoRI","BamHI"]}}
{"RenderSequenceSvg":{"seq_id":"pgex","mode":"Linear","path":"pgex.linear.restriction.svg"}}
```

Example: linear DNA-window export with explicit JASPAR/TFBS display state:

```json
{"AnnotateTfbs":{"seq_id":"pgex","motifs":["SP1","CTCF","TATAAA"],"min_llr_quantile":0.95,"clear_existing":true,"max_hits":200}}
{"SetParameter":{"name":"show_tfbs","value":true}}
{"SetParameter":{"name":"tfbs_display_use_llr_bits","value":false}}
{"SetParameter":{"name":"tfbs_display_min_llr_quantile","value":0.98}}
{"SetParameter":{"name":"tfbs_display_use_true_log_odds_quantile","value":true}}
{"SetParameter":{"name":"tfbs_display_min_true_log_odds_quantile","value":0.98}}
{"RenderSequenceSvg":{"seq_id":"pgex","mode":"Linear","path":"pgex.linear.tfbs.svg"}}
```

Render RNA secondary-structure SVG (engine operation):

```json
{"RenderRnaStructureSvg":{"seq_id":"rna_seq","path":"rna.secondary.svg"}}
```

Render lineage SVG (engine operation):

```json
{"RenderLineageSvg":{"path":"lineage.svg"}}
```

- Gibson apply operations export as the same dedicated `Gibson cloning` hub
  shown in the GUI, instead of a raw parent-to-each-output multi-edge graph.

Render pool gel SVG with ladder selection (engine operation):

```json
{"RenderPoolGelSvg":{"inputs":["frag_1","frag_2","frag_3"],"path":"digest.gel.svg","ladders":["NEB 100bp DNA Ladder","NEB 1kb DNA Ladder"]}}
```

Render pool gel SVG with automatic ladder selection:

```json
{"RenderPoolGelSvg":{"inputs":["frag_1","frag_2","frag_3"],"path":"digest.auto.gel.svg","ladders":null}}
```

The shared-shell route `render-pool-gel-svg` returns the SVG result plus
`gel_band_rows[]` and `gel_summary_lines[]`, giving CLI/Telegram/ClawBio callers
a text-first lane/band description even when they do not display the SVG.

Persist a ladder pair on an existing arrangement:

```json
{"SetArrangementLadders":{"arrangement_id":"arrangement-2","ladders":["Plasmid Factory 1kb DNA Ladder","GeneRuler 100bp DNA Ladder Plus"]}}
```

Clear an arrangement back to automatic ladder selection:

```json
{"SetArrangementLadders":{"arrangement_id":"arrangement-2","ladders":null}}
```

Prepare a whole reference genome once (download/copy sequence + annotation and
build local FASTA and BLAST indexes):

```json
{"PrepareGenome":{"genome_id":"Human GRCh38 Ensembl 116","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```

Extract a genomic interval from the prepared cache (1-based inclusive
coordinates):

```json
{"ExtractGenomeRegion":{"genome_id":"Human GRCh38 Ensembl 116","chromosome":"1","start_1based":1000000,"end_1based":1001500,"output_id":"grch38_chr1_1000000_1001500","annotation_scope":"core","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```


## Genome Extraction and Projection

Extract by gene query (name/id) from the prepared gene index:

```json
{"ExtractGenomeGene":{"genome_id":"Human GRCh38 Ensembl 116","gene_query":"TP53","occurrence":1,"output_id":"grch38_tp53","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```

Extract a promoter slice directly from transcript TSS geometry:

```json
{"ExtractGenomePromoterSlice":{"genome_id":"Human GRCh38 Ensembl 116","gene_query":"TERT","transcript_id":"ENST00000310581","output_id":"grch38_tert_promoter","upstream_bp":1000,"downstream_bp":200,"annotation_scope":"core","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```

Extend an anchored sequence on the contextual 5' side by 250 bp:

```json
{"ExtendGenomeAnchor":{"seq_id":"grch38_tp53","side":"five_prime","length_bp":250,"output_id":"grch38_tp53_ext5","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```

Import BED intervals (plain `.bed` or gzipped `.bed.gz`) onto a genome-anchored
sequence:

```json
{"ImportGenomeBedTrack":{"seq_id":"grch38_tp53","path":"data/chipseq/peaks.bed.gz","track_name":"H3K27ac","min_score":10.0,"max_score":null,"clear_existing":false}}
```

Import BigWig signal tracks (`.bw` / `.bigWig`) onto a genome-anchored
sequence:

```json
{"ImportGenomeBigWigTrack":{"seq_id":"grch38_tp53","path":"data/chipseq/signal.bw","track_name":"ATAC","min_score":0.2,"max_score":null,"clear_existing":false}}
```

Import VCF variants (`.vcf` / `.vcf.gz`) onto a genome-anchored sequence:

```json
{"ImportGenomeVcfTrack":{"seq_id":"grch38_tp53","path":"data/variants/sample.vcf.gz","track_name":"Variants","min_score":20.0,"max_score":null,"clear_existing":false}}
```

Project prepared microarray contrast tracks onto a genome-anchored sequence:

```json
{"ProjectMicroarrayTrack":{"seq_id":"grch38_tp73","manifest_path":"data/publication_resources/rostock_p73_clariomd_e_mtab_14704/analysis/clariomd_probe_level/clariomd_microarray_track_manifest.json","contrasts":["AdTAp73alpha-AdGFP","AdTAp73beta-AdGFP"],"level":"probeset","min_abs_logfc":0.5,"max_adj_p":0.05,"max_features":5000,"clear_existing":true}}
```

Project one interval through an explicit build-to-build interval map:

```json
{"ProjectGenomeInterval":{"source_genome_id":"hg19","target_genome_id":"hg38","projection_path":"data/projections/hg19-to-hg38.interval-map.tsv","chrom":"chr1","start_1based":1010,"end_1based":1020,"strand":"+"}}
```

Notes:

- `PrepareGenome` is intended as a one-time setup step per genome and cache
  location.
- During prepare, GENtle also attempts to create a BLAST nucleotide index
  (`makeblastdb`) under the genome install directory.
- HTTP-based source downloads are resumable (Range requests with retry/backoff),
  and completed installs persist legacy SHA-1 fields for sequence/annotation
  files only when a platform checksum tool is available.
- A catalog entry can either define explicit URLs (`sequence_remote` /
  `annotations_remote`) or define `ncbi_assembly_accession` +
  `ncbi_assembly_name` to derive direct NCBI GenBank/RefSeq FTP downloads.
- Helper/vector entries can also define `genbank_accession`; when explicit
  URLs are absent, GENtle derives direct NCBI EFetch sources
  (`rettype=fasta` + `rettype=gbwithparts`) for one-time prepare/index.
- Catalog entries may also include optional physical parameters:
  `nucleotide_length_bp` and `molecular_mass_da`.
- If `molecular_mass_da` is omitted but nucleotide length is available
  (directly or from prepared FASTA index), GENtle estimates dsDNA molecular
  mass and labels it as `estimated_from_nucleotide_length`.
- `ExtractGenomeRegion` expects the genome to have been prepared already.
- `ExtractGenomeGene` also expects prepared cache and gene index.
- `ExtractGenomePromoterSlice` also expects prepared cache and transcript/gene index.
- `genomes/helpers blast` expects prepared cache and a BLAST index.
  If index files are missing, GENtle tries to build them on demand.
- BLAST progress:
  - BLAST+ does not provide a native percent-progress CLI flag for `blastn`.
  - GENtle surfaces deterministic progress at orchestration level (query counts,
    running elapsed time, and invocation context in adapters that support live status).
  - Async BLAST jobs additionally expose queue/scheduler state:
    `queued|running|completed|failed|cancelled` and scheduler counters.
  - Scheduler concurrency defaults to available CPU parallelism and can be
    overridden with `GENTLE_BLAST_ASYNC_MAX_CONCURRENT` (`1..256`).
- BLAST options layering (engine contract):
  - built-in defaults (`task=blastn-short`, `max_hits=25`)
  - optional defaults file (`assets/blast_defaults.json` or path set by parameter)
  - optional project override metadata (`set-param blast_options_override ...`)
  - quick flags (`--max-hits`, `--task`)
  - per-command JSON override (`--options-json` / `--options-file`)
  - project defaults file path can be set via
    `set-param blast_options_defaults_path '"path/to/blast_defaults.json"'`
- BLAST executable overrides:
  - `GENTLE_MAKEBLASTDB_BIN` (default: `makeblastdb`)
  - `GENTLE_BLASTN_BIN` (default: `blastn`)
- Legacy SHA-1 verification override:
  - `GENTLE_SHA1_TOOL` (optional explicit `sha1sum`/`shasum`/`certutil` path)
  - `GENTLE_DISABLE_LEGACY_SHA1=1` to disable legacy SHA-1 verification
- When BLAST hits are imported as features (`ImportBlastHitsTrack` /
  shell `genomes|helpers blast-track`), operation payload now includes
  `blast_provenance` with invocation details (`blastn_executable`, `blast_db_prefix`,
  raw `command` args, and compact `command_line`) so sequence history can
  trace how hits were generated.
- `ImportGenomeBedTrack` expects `seq_id` to be a sequence created by
  `ExtractGenomeRegion`, `ExtractGenomeGene`, or `ExtendGenomeAnchor`
  (genome-anchored provenance).
- `ImportGenomeBigWigTrack` expects the same genome-anchored `seq_id`.
- `ImportGenomeVcfTrack` expects the same genome-anchored `seq_id`.
- `ProjectMicroarrayTrack` expects the same genome-anchored `seq_id` and a
  `gentle.microarray_track_manifest.v1` manifest. Projection is refused unless
  the manifest `coordinate_system` / `supported_genome_ids` match the sequence
  anchor `genome_id` or the manifest declares an explicit
  `coordinate_projections[]` map into that anchor build; GENtle does not guess
  array coordinate builds.
- `arrays probe-regions` is currently a plan/preflight command for arbitrary
  Affymetrix CEL inspection. It emits `gentle.probe_region_plan.v1` JSON with
  CEL file size/mtime cache keys, parsed TSV/CSV/SDRF-style metadata previews,
  default condition contrasts, explicit output/cache path status, annotation
  source checks, resource-registry platform hints such as
  `Clariom_D_Human -> pd.clariom.d.human`, backend-candidate readiness for
  `r_oligo`, Affymetrix Power Tools, and legacy 3' IVT/CDF `r_affy_cdf`
  candidates, an advisory
  `scripts/probe_regions_oligo.R` command for explicit RMA/CEL requests, an
  advisory `apt-probeset-summarize -a rma-sketch` command when user-supplied
  PGF/CLF and optional MPS files are detected, and local `Rscript` / APT /
  R-package dependency status; it does not run CEL summarization itself and
  never downloads or installs missing files/packages. Platform recognition now
  comes from `data/resources/affymetrix/platform_registry.json`
  (`gentle.affymetrix_platform_registry.v1`): verified entries may report
  concrete local support-file expectations, while provisional historical chips
  identify their family/backend/CDF needs without claiming local validation.
  With `--output DIR`, the preflight also writes `DIR/plan.json`, a pretty
  versioned copy of the same `gentle.probe_region_plan.v1` report returned on
  stdout.
  Repeat `--r-library-path PATH` to make an agent sandbox, user library, or
  system package tree explicit. The normalized request records those paths,
  `r_library_paths_checked` records R's complete effective search path, and
  every generated R helper command receives the same values. The preflight
  checks all direct helper dependencies in one bounded R process without
  attaching package namespaces and reports their exact versions. A timeout or
  probe failure is distinct from a missing or unchecked package; each
  diagnostic names the paths inspected and tells the user to verify
  `--r-library-path` if the result conflicts with their interactive R session.
  With `--dataset E-MTAB-14704`, the preflight resolves the publication
  resource's declared local CEL paths and any locally present SDRF metadata, then
  reports missing raw files as ordinary file-status errors.
  For registry entries such as `Clariom_D_Human`, the same preflight reports
  the expected local Thermo Fisher na36 hg38 support ZIPs under
  `annotation_source.vendor_support_files[]`; place them manually in
  `data/resources/affymetrix/clariom_d_human_na36_hg38/` when probe/probeset
  coordinate development needs vendor CSV annotations. The preflight accepts
  both the concise canonical ZIP names documented in that directory and the
  browser-preserved `TFS-Assets_LSG_Support-Files_...` download names.
  Legacy 3' IVT arrays such as HG-U133 / Mouse 430 / Rat 230 families are
  recognized as provisional CDF-backed platforms; their `r_affy_cdf` candidate
  renders an explicit `scripts/probe_regions_affy.R` command when CEL inputs
  are present, while readiness still depends on local R/`affy`, `limma`,
  `Biobase`, CDF, and annotation resources.
- `arrays run-probe-region-backend PLAN.json --allow-external-execution`
  (or `arrays run-probe-region-backend --plan PLAN.json --allow-external-execution`)
  reads a persisted `gentle.probe_region_plan.v1`, checks the recorded
  preflight/backend readiness, and only then runs the rendered R/oligo or APT
  command. The `--allow-external-execution` gate is required; without it GENtle
  refuses external R/APT execution. The command captures stdout/stderr/exit
  status, validates the resulting four-file helper-output contract, and writes
  hardened `gentle.probe_region_backend_provenance.v1` provenance. It never
  downloads or installs CEL files, vendor resources, R packages, or APT.
- `arrays inspect-probe-region-output OUTPUT_DIR` validates the explicit helper
  outputs after the user has run R themselves. It does not project features;
  it prepares a shared GUI/CLI-readable summary for projection triage.
- `arrays render-probe-region-output-svg OUTPUT_DIR OUTPUT.svg` draws the
  inspected helper output natively as a deterministic SVG plot, without running
  R/APT or changing project state.
- `arrays render-probe-region-evidence-svg REPORT.json OUTPUT.svg` draws an
  already-exported probe-region evidence interpretation report as a
  deterministic SVG. It visualizes transcript lanes, report-provided exon
  overlaps, report-provided junction spans, parent probeset spans, and PM probe
  intervals. It draws only geometry present in the report and remains a
  review-only constraint view, not an isoform-support call.
  The same report is what the GUI Splicing Expert uses to label
  junction-spanning array probes, so CLI-generated reports are the parity path
  for checking those rows outside the GUI.
- `arrays project-probe-region-output SEQ_ID OUTPUT_DIR` projects selected or
  all `log2FC_*` helper-output rows into genome-anchored array features when
  the helper output's declared `coordinate_system` or `genome_build` directly
  matches the sequence anchor genome id, or when a declared
  `coordinate_projections[]`/`projection_maps[]` interval map bridges the
  helper-output build to the anchor build.
- `arrays interpret-probe-region-evidence SEQ_ID` reads those projected
  features back from the sequence and summarizes transcript/exon geometry
  compatibility as a read-only report, including structured exon ordinals,
  exon ranges, junction spans, overlap base counts, conservative geometry
  scores, score-basis guardrails, and review-only transcript labels for each
  evidence-to-transcript mapping. It does not run R/APT and does not claim
  isoform support until probe/probeset-to-transcript ambiguity has been
  audited.
  When the effect threshold participates in a study decision, also supply
  `--threshold-source` and `--policy-sha256`. Legacy or ad-hoc reports without
  those bindings remain geometry evidence but cannot establish
  response-associated assay retention.
- The built-in genome catalog includes both `Human GRCh38 Ensembl 116` and
  `Human GRCh37 Ensembl 87` (`hg19`/`GRCh37.p13` aliases), so direct native
  extraction can use either build when the corresponding cache has been
  prepared.
- BED import accepts local `.bed` and `.bed.gz` files.
- Concatenated gzip members are accepted for `.bed.gz` track input.
- BigWig import accepts local `.bw` and `.bigWig` files and uses
  `bigWigToBedGraph` (override with `GENTLE_BIGWIG_TO_BEDGRAPH_BIN`).
- VCF import accepts local `.vcf` and `.vcf.gz` files.
- Concatenated gzip members are accepted for `.vcf.gz` track input.
- For VCF import, `min_score` / `max_score` filter on VCF `QUAL`.
- `ExtractGenomeRegion`, `ExtractGenomeGene`, and `ExtendGenomeAnchor` append extraction provenance
  records into `ProjectState.metadata["provenance"]["genome_extractions"]`
  (genome id, coordinates/query, source descriptors, and checksums when present).
- If `catalog_path` is omitted, engine default catalog is `assets/genomes.json`.
- Bundled `assets/genomes.json` currently includes Human GRCh38 (Ensembl 113 and 116),
  Mouse GRCm39 Ensembl 116, Rat GRCr8 Ensembl 115, Danio rerio GRCz11
  Ensembl 115, Pan troglodytes Pan_tro_3.0 Ensembl 115, Canis lupus familiaris
  ROS_Cfam_1.0 Ensembl 115, Drosophila melanogaster BDGP6.54 Ensembl Metazoa 62,
  Caenorhabditis elegans WBcel235 Ensembl 115, Saccharomyces cerevisiae S288c
  (Ensembl 113 and 115), and `LocalProject` (backed by
  `test_files/fixtures/genomes/AB011549.2.fa` +
  `test_files/fixtures/genomes/AB011549.2.gb`).
- `cache_dir` is optional. If omitted, catalog/default cache settings are used.
- `PrepareGenome` now validates that gene-bearing contigs parsed from the
  prepared annotation are present in the prepared FASTA index; truncated or
  mismatched installs fail during preparation instead of only later during
  extraction.
- Catalog rows with `ensembl_template` metadata can be refreshed on demand via
  `update-ensembl-specs`, which rewrites explicit pinned URLs without preparing
  any genomes.
- `chromosome` accepts exact contig names and also tolerates `chr` prefix
  differences (`1` vs `chr1`).
- Missing-contig extraction errors now also report the prepared `sequence.fa`
  and `sequence.fa.fai` paths for cache debugging.
- For `ExtractGenomeGene`, `occurrence` is 1-based among matching records.
- `ExtractGenomeGene.extract_mode=coding_with_promoter` resolves the CDS span
  from transcript `CDS` annotation and applies `promoter_upstream_bp`
  strand-aware on the gene's 5' side.
- For `ExtendGenomeAnchor`, `side` is contextual to anchor strand.
  On anchor strand `-`, `5'` increases physical genomic position.
- For anchor-extension reads, genome ids can resolve through assembly-family
  compatibility (for example `GRCh38.p14` using a prepared `Human GRCh38 ...`
  cache) when unique; ambiguous matches are rejected with options.

Available `target` values:

- `SequencePanel`
- `MapPanel`
- `Features`
- `CdsFeatures`
- `GeneFeatures`
- `MrnaFeatures`
- `Tfbs`
- `RestrictionEnzymes`
- `GcContents`
- `OpenReadingFrames`
- `MethylationSites`


## Save Sequence Files

Save as GenBank:

```json
{"SaveFile":{"seq_id":"pgex_frag_1","path":"frag1.gb","format":"GenBank"}}
```

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Independent Command Examples](command-examples.md) | [Next: Operation Limitations](limitations.md)
<!-- manual-page-nav:end -->
