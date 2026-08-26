# GENtle Engine Protocol (Draft v1)

## Gene-set publication report

`gentle.gene_set_publication_request.v1` is the portable multi-gene publication
input. It resolves to `gentle.gene_set_publication_report.v1`, which normalises
the complete primer table and copied asset paths. The generation receipt is
`gentle.gene_set_publication_generation.v1`. HTML and printable output are
projections of this record, not independent sources. Every generated legacy
bundle also carries `gentle.gene_set_publication_bundle_manifest.v1`, which
binds the exact request SHA-256 and every finalized user-facing artifact. The
manifest excludes itself and `generation-report.json` to avoid recursive
hashes; the generation receipt records the manifest path and SHA-256.

Gene-centred isoform-assay dossiers use
`gentle.gene_isoform_assay_publication_request.v1`. Every study plan,
experimental handoff, and oligo order form is supplied with an exact SHA-256;
GENtle embeds each source JSON once in
`gentle.gene_isoform_assay_publication.v1`. Its `content_blocks[]` contain only
JSON pointers and named deterministic projections, while `profiles[]` contain
ordered block ids. The renderer emits one meta-page, one page per gene, one
print HTML document, optional browser-derived PDF, order-sheet TSVs, and a
`gentle.gene_isoform_assay_publication_projection.v1` receipt. Selecting a
profile or explicit declared blocks changes presentation scope only. MCP
clients invoke the same generator through `gene_isoform_assay_publication`;
its required `confirm=true` authorizes artifact writes, not another scientific
decision. Experimental handoffs are admitted only when their panel report id,
sequence id, and source feature match an exact operation payload in the bound
study plan. Every standard profile includes a `quality_assurance` projection
that exposes each handoff's named readiness policy, source-panel hash,
per-assay required/optional gates, recorded status and GENtle summary,
evidence-report ids, blockers, and warnings. Annotated-transcript product
testing and whole-transcriptome/cDNA specificity remain separate gates; an
absent, incomplete, or `not_evaluated` gate is never rendered as a pass.
Every standard profile also includes an explanatory `provenance` projection.
It identifies the role, schema, report id, path, and SHA-256 of the study plan,
experimental handoff, and readiness-bound order form, plus their bound panel,
policy, operation-batch, and handoff identities where present. These rows
describe what GENtle processed and content-bind the exact source bytes; they do
not promote user-supplied/external claims into GENtle results or make approval
a biological validation.

A dossier may be published while part of a study is still running or has
failed, so each request gene carries a `status` of `resolved`, `pending`, or
`unresolved` plus an optional `status_reason`. An omitted status is derived
from the gene's contents — a gene with at least one experimental handoff is
`resolved` and a gene without one is `pending` — so a request written before
this field existed still reports honestly instead of presenting a gene without
results as finished. A declared status is checked rather than trusted:
`unresolved` requires a `status_reason` stating why the established automatism
could not address the gene, `resolved` without any handoff is rejected, and
`pending`/`unresolved` with handoffs present is rejected. The canonical report
adds `resolved_gene_count`, `pending_gene_count`, `unresolved_gene_count`, and
`complete`, and gains a warning naming the outstanding genes while `complete`
is false. The index lists every gene with its status and reason, and each
incomplete gene page and the print document open with a visible notice, so the
finished genes are usable immediately and the dossier is regenerated from the
same request once the outstanding results arrive. The natural producer of an
`unresolved` reason is a `gene_failures[]` entry from
`primers execute-gene-isoform-study-workflow-batch --on-gene-failure continue`.

This document defines the draft machine-facing protocol for operating GENtle
through a shared core engine.

Goal:

- GUI, CLI, JavaScript, Lua, and Python wrappers call the same core routines.
- AI tools can run deterministic cloning workflows with reproducible logs.

## Design principles

- Protocol-first: versioned JSON request/response shapes
- Capability negotiation: clients discover supported operations and formats
- Deterministic operation log: each operation emits a stable op id and result
- Structured errors: machine-parseable error code + message

## Capabilities

`gentle_cli capabilities` returns:

- `protocol_version`
- `supported_operations`
- `supported_export_formats`
- `deterministic_operation_log`
- `capability_registry`

`capability_registry` is the shared discovery surface projected by CLI, MCP,
JavaScript, and Lua adapters. Each row has a stable `name`, glossary-sourced
`description` where applicable, JSON-schema `input_schema`/`output_schema`,
`mutating` (`false`, `true`, or `external`), per-adapter surfacing
(`prominent`, `shell-only`, `gap`, or `n/a`), required
not-applicable justifications backed by
[`parity_matrix_overrides.json`](parity_matrix_overrides.json), and optional
`inline_operand_ok`. Glossary commands with declared engine operations project
as shell-pass-through on adapters that expose a generic operation/workflow
route.

Adapter error payloads that cross machine boundaries use the shared
`EngineError` shape (`code`, `message`, optional `cause_chain`) so adapters can
preserve lower-level string failures without changing their transport.

## Protocol-first workflow examples

Canonical, adapter-independent examples are defined in:

- `docs/examples/workflows/*.json`
- schema: `gentle.workflow_example.v1`

Each example includes:

- metadata (`id`, `title`, `summary`)
- test policy (`test_mode`: `always|online|skip`)
- required local files (`required_files`)
- canonical `workflow` payload

Adapter snippets (CLI/shared shell/JavaScript/Lua) are generated on demand from
those canonical files:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Validation only:

```bash
cargo run --bin gentle_examples_docs -- --check
```

Tutorial manifest + generated outputs:

- discovery catalog: `docs/tutorial/catalog.json`
- discovery schema: `gentle.tutorial_catalog.v2`
- shared tutorial source units:
  - `docs/tutorial/sources/catalog_meta.json`
  - `docs/tutorial/sources/*.json`
- source-unit schemas:
  - `gentle.tutorial_catalog_meta.v2`
  - `gentle.tutorial_source.v4`
  - `gentle.tutorial_source.v3`
  - `gentle.tutorial_source.v2`
- generated runtime manifest: `docs/tutorial/manifest.json`
- runtime manifest schema: `gentle.tutorial_manifest.v2`
- review freshness manifest: `docs/tutorial/review_manifest.json`
- review freshness schema: `gentle.tutorial_review_manifest.v1`
- committed generated outputs: `docs/tutorial/generated/`

Catalog/manifest split:

- `docs/tutorial/catalog.json` is the canonical discovery layer for all
  tutorials, including hand-written walkthroughs and agent/reference guides.
- Catalog/manifest v2 add optional group placement fields. `decimal_id` is
  nullable/absent so reference/navigation units such as the generated hub and
  tutorial landscape can stay catalogued without receiving a tutorial number.
- Catalog entries and generated chapter reports expose derived review routing
  fields: `review_stale_reason`, `review_issue_template`, and
  `review_issue_template_path`. These are derived from
  `docs/tutorial/review_manifest.json` and from dependency freshness checks;
  they are not authored in tutorial source JSON.
- Tutorial source units may declare `graphics[]` entries with `kind`
  (`generated` or `screenshot`), `path`, `caption`, `illustrates_step`, and
  either `regen_command` or `capture_date`. Generated chapters carry those
  declarations through the manifest/report, embed declared generated figures at
  the step they illustrate, and keep the trailing artifact list as an audit
  appendix. Declared graphics also participate in review-staleness checks.
- `docs/tutorial/sources/` is the authoring layer for both the discovery
  catalog and the executable tutorial runtime manifest.
- `docs/tutorial/manifest.json` is a generated runtime contract used for
  chapter output and tutorial runtime checks.
- `docs/tutorial/review_manifest.json` stores non-executable review metadata
  keyed by tutorial id. The validated tutorial source JSON remains the
  executable contract; review metadata tracks `tutorial_kind`, active or
  deprecated status, optional replacement, and optional `codex_reviewed_at` /
  `human_reviewed_at` dates.
- GUI help/tutorial discovery may consume the catalog directly for curated
  ordering and metadata, while executable tutorial project materialization still
  resolves through the manifest/workflow example path.

Review-manifest checks are warnings, not hard failures:

- missing entries for known tutorial ids warn
- entries for unknown tutorial ids warn
- `human_reviewed_at` older than `warn_after_months` warns
- if a human review date predates the tutorial source JSON, workflow JSON,
  catalog Markdown page, or declared graphics file, the tutorial is marked
  stale with a dependency-specific reason
- `deprecated` tutorials, or entries with `replaced_by`, do not escalate
  execution failures during tutorial checks

Generate/check tutorial outputs:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_examples_docs -- tutorial-catalog-generate
cargo run --bin gentle_examples_docs -- tutorial-catalog-check
cargo run --bin gentle_examples_docs -- tutorial-manifest-generate
cargo run --bin gentle_examples_docs -- tutorial-manifest-check
```

## TFBS region summary contract

GENtle now exposes a portable grouped TFBS summary contract for comparing one
focus window against a wider context window on the same sequence.

Current shared-shell route:

```bash
gentle_cli shell 'features tfbs-summary SEQ_ID --focus START..END [--context START..END] [--min-focus-count N] [--min-context-count N] [--limit N]'
```

First-class operation route:

```json
{"SummarizeTfbsRegion":{"seq_id":"SEQ_ID","focus_start_0based":2900,"focus_end_0based_exclusive":3100,"context_start_0based":0,"context_end_0based_exclusive":6001,"min_focus_occurrences":1,"min_context_occurrences":0,"limit":25}}
```

Portable schema:

- `gentle.tfbs_region_summary.v1`

Request fields:

- `seq_id`
- `focus_start_0based`
- `focus_end_0based_exclusive`
- optional `context_start_0based`
- optional `context_end_0based_exclusive`
- `min_focus_occurrences`
- `min_context_occurrences`
- optional `limit`

Result fields:

- sequence/focus/context bounds and widths
- total TFBS hit counts in the focus and context spans
- grouped rows keyed by TF name with:
  - `motif_ids`
  - `focus_occurrences`
  - `context_occurrences`
  - `outside_focus_occurrences`
  - focus/context/outside densities per kb
  - focus-vs-context and focus-vs-outside density ratios

Grouping policy:

- prefer `bound_moiety`
- otherwise `standard_name`
- otherwise `gene`
- otherwise `name`
- otherwise `tf_id`

## TFBS score-track contract

GENtle also exposes a portable continuous TFBS/PSSM score-track contract for one
selected sequence span.

Boundary note:

- This is GENtle's generic DNA motif-statistics surface today.
- It is the right contract for:
  - local JASPAR/TF motif scoring
  - score-track figures and JSON reports over one DNA span
  - deterministic motif benchmark/report flows such as
    `SummarizeJasparEntries`
- It is intentionally distinct from splice-aware ATtRACT/RBP evidence:
  - TFBS/PSSM score tracks describe score landscapes over one selected DNA span
  - ATtRACT evidence describes transcript-aware splice-region hits with exon /
    donor flank / acceptor flank / intron classification
- If future ATtRACT PWM scoring reuses PSSM-like math, that reuse should live
  in a lower shared helper layer while keeping the adapter-facing contracts
  distinct.

Current shared-shell route:

```bash
gentle_cli shell 'features tfbs-score-tracks-svg SEQ_ID OUTPUT.svg --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative]'
```

```bash
gentle_cli shell 'features tfbs-score-tracks-svg --sequence-text DNA --output OUTPUT.svg [--topology linear|circular] [--id-hint TEXT] --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative]'
gentle_cli shell 'features tfbs-score-track-correlation-svg SEQ_ID OUTPUT.svg --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--correlation-metric pearson|spearman] [--signal-source max_strands|forward_only|reverse_only] [--allow-negative]'
```

First-class operation routes:

```json
{"SummarizeTfbsScoreTracks":{"target":{"kind":"seq_id","seq_id":"SEQ_ID","span_start_0based":15564,"span_end_0based_exclusive":16764},"motifs":["TP53","TP63","TP73","PATZ1","SP1","BACH2","REST"],"score_kind":"llr_background_tail_log10","clip_negative":false}}
```

```json
{"RenderTfbsScoreTracksSvg":{"target":{"kind":"inline_sequence","sequence_text":"ACGTACGTACGT","topology":"linear","id_hint":"inline_promoter","span_start_0based":0,"span_end_0based_exclusive":12},"motifs":["TP53","TP63","TP73","PATZ1","SP1","BACH2","REST"],"score_kind":"llr_background_tail_log10","clip_negative":false,"path":"docs/figures/tp73_upstream_tfbs_score_tracks.svg"}}
```

```json
{"RenderTfbsScoreTrackCorrelationSvg":{"seq_id":"SEQ_ID","motifs":["TP53","TP63","TP73","PATZ1","SP1","BACH2","REST"],"start_0based":15564,"end_0based_exclusive":16764,"score_kind":"llr_background_tail_log10","correlation_metric":"spearman","correlation_signal_source":"forward_only","clip_negative":false,"path":"docs/figures/tp73_upstream_tfbs_score_tracks_correlation.svg"}}
```

Portable schema:

- `gentle.tfbs_score_tracks.v1`

Behavior notes:

- `SummarizeTfbsScoreTracks` returns the structured per-position forward/reverse
  score arrays over a `SequenceScanTarget`, so the same report path works for
  stored `seq_id` spans and inline ASCII DNA.
- motif tokens in `motifs[]` resolve through the same shared TF-query layer:
  - exact motif ids / TF names
  - aliases such as `OCT4`
  - catalog-backed functional groups such as `Yamanaka factors` / `stemness`
  - family-like queries such as `KLF family`
- compact motif entries without PFM rows are supplemented from the bundled
  full-PFM table when the match is unambiguous. Consensus-derived matrices are
  therefore a last-resort fallback, not the normal source for TFBS scoring or
  sequence-logo display.
- when the shell command omits an explicit range, the full selected
  `SequenceScanTarget` span is used.
- when the target is a stored `seq_id`, the same shared report can also carry
  `overlay_tracks[]` for imported genome-track evidence already materialized on
  that sequence (for example projected BED / CUT&RUN peak intervals):
  - `source_kind`
  - `track_name`
  - `display_label`
  - optional `source_file_name`
  - `interval_count`
  - optional `max_score`
  - `intervals[]` with clipped local `start_0based`, `end_0based_exclusive`,
    optional `label`, optional numeric `score`, and optional `strand`
  These overlays reuse the same deterministic report/renderer path instead of
  forcing GUI-only evidence lanes.
- the same shared report also carries correlation sidecars:
  - `correlation_summary` remains the default `max_strands` summary for
    compatibility
  - `correlation_summaries[]` can carry:
    - `max_strands`
    - `forward_only`
    - `reverse_only`
  Each summary records:
  - `signal_source`
  - `smoothing_method`
  - `smoothing_window_bp`
  - `pair_count`
  - `rows[]` with `left_tf_id`, optional `left_tf_name`, `right_tf_id`,
    optional `right_tf_name`, `overlap_window_count`, `raw_pearson`,
    `smoothed_pearson`, `raw_spearman`, `smoothed_spearman`, and optional
    `signed_primary_peak_offset_bp`
- SVG rendering includes a compact per-track sequence-logo preview derived from
  the resolved JASPAR matrix. Logo glyphs are scaled and clipped like a standard
  stacked letter logo while the columns remain readable; when a motif is too
  narrow for legible letters, the renderer switches to colored stacked bars
  instead of drawing misleading microscopic text.
- the same report now also carries one explicit
  `cross_strand_correlation_summary` for the correlation-export overview:
  - `pair_count`
  - `smoothing_method`
  - `smoothing_window_bp`
  - `rows[]` keyed by one TF pair
  - each row carries `cells[]` in deterministic `F-F`, `F-R`, `R-F`, `R-R`
    order with the same raw/smoothed Pearson/Spearman fields plus one optional
    strand-specific peak offset
  - the shared SVG export now uses that overview when
    `correlation_signal_source=max_strands`, but it renders the result as one
    strand-expanded all-vs-all matrix:
    - every motif appears twice on each axis, ordered `F` then `R`
    - each TF-pair therefore reads naturally as one adjacent 2x2 block
      (`F-F / F-R / R-F / R-R`)
    - no prior `max(forward, reverse)` aggregation is used for those block
      values
- each returned track now also carries one `directional_summary` describing how
  its own forward and reverse curves travel together over the selected span:
  - optional `forward_primary_peak_position_0based`
  - optional `reverse_primary_peak_position_0based`
  - optional `signed_primary_peak_offset_bp`
  - `raw_pearson`
  - `smoothed_pearson`
  - `raw_spearman`
  - `smoothed_spearman`
- `score_kind` selects which per-window value is exported:
  - `llr_bits`
  - `llr_quantile`
  - `llr_background_quantile`
  - `llr_background_tail_log10`
  - `true_log_odds_bits`
  - `true_log_odds_quantile`
  - `true_log_odds_background_quantile`
  - `true_log_odds_background_tail_log10`
- the background-normalized score kinds suppress everything below the `0.95`
  modeled random-background quantile, so the plot can focus on the unusual
  tail instead of the visually noisy body of the background
- each returned track now also carries a deterministic normalization reference
  computed from two shared sources:
  - a larger deterministic random sample (`100000 bp`) used for summary
    moments such as `mean_score`, `stddev_score`, and `positive_fraction`
  - a quantized IID random-DNA window model used for upper-tail calibration,
    theoretical bounds, and the non-saturating background percentile/tail
    views
  with the same underlying score family and the same clipping semantics as the
  exported track itself:
  - `chance_model`
  - `mean_score`
  - `stddev_score`
  - `p95_score`
  - `p99_score`
  - `positive_fraction`
  - `observed_peak_empirical_quantile`
  - `observed_peak_modeled_quantile`
  - `observed_peak_modeled_tail_probability`
  - `observed_peak_modeled_tail_log10`
  - `observed_peak_delta_from_p95`
  - `observed_peak_delta_from_p99`
  - `theoretical_min_score`
  - `theoretical_max_score`
- each returned track now also carries up to three highlighted motif windows
  (`top_peaks`) chosen from local maxima on the rendered strands:
  - windows at or above the deterministic random-background `p99` threshold are
    preferred
  - if none cross that threshold, the single best local peak is still kept as a
    fallback anchor
  - each peak records start/end, strand, score, empirical quantile, and
    `delta_from_p99`
- each returned track now also carries a compact `motif_logo_columns` payload
  with the same per-position base counts/fractions/information-content fields
  used by the JASPAR expert, so static renderers and specialist surfaces can
  reuse one shared motif/logo payload without doing adapter-only motif lookups
- `RenderTfbsScoreTracksSvg` reuses the same shared report payload and writes a
  deterministic stacked SVG figure suitable for GUI/CLI/agent/README parity.
- `RenderTfbsScoreTrackCorrelationSvg` reuses that exact same report and writes
  a second deterministic SVG view:
  - left panel: smoothed selected-metric heatmap
  - right panel: raw selected-metric heatmap
  - footer: top synchronized pairs with signed primary-peak offsets
  - `correlation_metric` accepts `pearson` or `spearman`
  - `signal_source` accepts:
    - `max_strands` for the current orientation-agnostic neighborhood view
    - `forward_only`
    - `reverse_only`
  - `spearman` is often the safer interpretation mode when promoter score
    tracks are clearly non-normal or contain many tied values from clipping
- the SVG figure now prints score-family-aware normalization labels per motif:
  - raw bit views keep the compact `p99 / Δp99 / bg+` summary
  - background-normalized views instead show `theory max / peak q / -log10 tail`
    so a visually busy tail plot still carries its calibration context
- The shared report labels the source span through `target_kind`,
  `target_label`, `source_sequence_length_bp`, and `scan_topology` so adapter
  surfaces do not need live GUI state to explain what was scored.
- The shared report also carries transcription-start markers for covered or
  directly adjacent starts, and the SVG renderer shows them as short hooked
  arrows so strand direction survives figure-oriented exports; the TSS label
  now also carries a small outline box so it survives README/downsampled
  figure use better.
- `SummarizeTfbsScoreTracks`, `RenderTfbsScoreTracksSvg`, and
  `RenderTfbsScoreTrackCorrelationSvg` now also emit `OperationProgress::Tfbs`
  events with `task_kind=score_tracks`, `stage_label`, and `stage_percent` so
  GUI/agent callers can show progress while deterministic background
  calibration and target scanning are still running.
  Those progress callbacks are now also honored cooperatively for
  `SummarizeTfbsScoreTracks`, so returning `false` can cancel a long-running
  score-track job during background calibration or target scanning.
- `clip_negative=true` is the default presentation mode for promoter-design
  inspection on the bit-based score kinds because it suppresses negative-only
  windows and leaves only positive support.
- `clip_negative=false` keeps the full continuous score landscape, which is
  useful for figure/report contexts where one wants to show that only some
  factors cross into positive support while the others remain below zero.

## Multi-gene promoter TFBS contract

GENtle now also exposes one promoter-aligned multi-gene TFBS comparison report
for cases where users want one shared upstream figure across several genes
instead of repeating a single-gene promoter workflow manually.

Current shared-shell routes:

```bash
gentle_cli shell 'genomes promoter-tfbs-summary GENOME_ID [--gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] ...|--gene-json JSON|--gene-set GROUP_OR_GO|--gene-set-resolution RESOLUTION.json] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--gene-group-catalog PATH] [--cache-dir PATH] [--path FILE.json]'
gentle_cli shell 'genomes promoter-tfbs-svg GENOME_ID [--gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] ...|--gene-json JSON|--gene-set GROUP_OR_GO|--gene-set-resolution RESOLUTION.json] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--gene-group-catalog PATH] [--cache-dir PATH] OUTPUT.svg'
gentle_cli shell 'genomes promoter-cohort-comparison GENOME_ID --cohort-label LABEL --cohort-kind manual|co_regulated|anti_co_regulated --gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] [--gene ...|--gene-json JSON] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--source-seq-id SEQ_ID ...] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--expression-json JSON] [--expression-source-label LABEL] [--cutrun-dataset ID ...] [--cutrun-read-report ID ...] [--catalog PATH] [--cache-dir PATH] [--path FILE.json]'
```

First-class operation routes:

```json
{"SummarizeMultiGenePromoterTfbs":{"genome_id":"Human GRCh38 Ensembl 116","genes":[{"gene_query":"TERT"},{"gene_query":"TP73"}],"motifs":["Yamanaka factors","SP1"],"upstream_bp":1000,"downstream_bp":200,"score_kind":"llr_background_tail_log10","clip_negative":true,"path":"artifacts/tert_tp73_promoter_tfbs.summary.json"}}
```

```json
{"SummarizePromoterCohortComparison":{"genome_id":"Human GRCh38 Ensembl 116","cohort_label":"TERT_TP73_manual","cohort_kind":"manual","genes":[{"gene_query":"TERT"},{"gene_query":"TP73"}],"motifs":["Yamanaka factors","SP1"],"upstream_bp":1000,"downstream_bp":200,"score_kind":"llr_background_tail_log10","clip_negative":true,"expression_source_label":"reviewed_table","expression_rows":[{"gene_label":"TERT","sample_id":"case_1","condition":"case","value":9.8,"unit":"TPM"}],"cutrun_dataset_ids":["toy_ctcf"],"path":"artifacts/tert_tp73_promoter_cohort.json"}}
```

```json
{"RenderMultiGenePromoterTfbsSvg":{"genome_id":"Human GRCh38 Ensembl 116","genes":[{"gene_query":"TERT"},{"gene_query":"TP73"}],"motifs":["Yamanaka factors","SP1"],"upstream_bp":1000,"downstream_bp":200,"score_kind":"llr_background_tail_log10","clip_negative":true,"path":"artifacts/tert_tp73_promoter_tfbs.svg"}}
```

Portable schema:

- `gentle.multi_gene_promoter_tfbs.v1`
- `gentle.promoter_cohort_comparison.v1`

Behavior notes:

- each `genes[]` entry uses the same shared TF-query and promoter-resolution
  machinery as the single-gene promoter route, including exact/fuzzy gene
  selection and optional transcript pinning
- shorthand shell gene tokens support:
  - `QUERY`
  - `QUERY::OCCURRENCE`
  - `QUERY@TRANSCRIPT_ID`
  - `QUERY#DISPLAY_LABEL`
  - combinations such as `QUERY::2@ENST...#short_label`
- every returned `genes[]` row carries:
  - resolved gene/transcript identifiers
  - promoter genomic interval
  - transcript TSS position
  - one nested `gentle.tfbs_score_tracks.v1` report
- minus-strand promoters are reverse-complemented before scoring, and
  `sequence_orientation="transcription_aligned"` makes that explicit so
  upstream/downstream comparisons stay meaningful across mixed-strand genes
- `summary_rows[]` flatten one per-gene/per-factor comparison table with:
  - `max_score`
  - `peak_position_0based`
  - `peak_position_promoter_relative_bp`
  - `peak_genomic_position_1based`
  - `positive_fraction`
- the SVG export renders one small-multiples promoter panel per gene with a
  shared promoter-relative axis convention and explicit TSS markers
- `SummarizePromoterCohortComparison` is a deterministic first-slice cohort
  report for `manual`, `co_regulated`, and `anti_co_regulated` cohorts only.
  It resolves each requested gene/transcript to a strand-aware promoter window,
  reuses the same multi-gene TFBS scoring and TFBS similarity ranking helpers,
  reports pairwise motif-track similarity, shared/common score peaks,
  cohort-specific or outlier peaks, optional expression-row association, and
  warnings for unresolved genes/transcripts. For declared `co_regulated` and
  `anti_co_regulated` cohorts, additive `relationship_flags[]` compare the
  available TFBS-track similarity evidence with that expectation and surface
  unexpected divergence or concordance as review cues, never as regulatory
  verdicts. `source_seq_ids`,
  CUT&RUN dataset ids, and saved read-report ids are retained for traceability
  in this slice.

## Ortholog promoter resources and comparisons

GENtle’s ortholog promoter v1 is offline-first and engine-owned. It does not
call Ensembl, GO, or any live orthology service during normal execution.

Portable schemas:

- `gentle.ortholog_resource.v1`
- `gentle.ortholog_promoter_cohort.v1`
- `gentle.ortholog_promoter_comparison.v1`

Behavior notes:

- `gentle.ortholog_resource.v1` is a local mapping table with
  `source_species`, source gene id/symbol, `target_species`, target gene
  id/symbol, `orthology_type`, `confidence`, `source`, `evidence[]`, and
  `species_aliases[]`. Additive context binding uses a flattened
  `contexts[]` registry plus optional `source_context_id` and
  `target_context_id` row references; old resources without those fields
  remain valid.
- `orthology_type` and `confidence` remain JSON strings for compatibility but
  deserialize into open typed vocabularies. Canonical one-to-one, one-to-many,
  many-to-one, and many-to-many spellings expose typed cardinality, and
  canonical high/medium/low confidence values expose typed levels. Unknown
  provider or fixture values round-trip exactly rather than being rejected or
  silently normalized.
- `ResolveOrthologPromoterCohort` resolves the anchor gene first, maps each
  target species through the local ortholog table, and derives promoter
  windows with the same prepared-genome promoter/TSS resolver used by
  `genomes extract-promoter`.
- Explicit context references must exist. Their organism and optional
  `genome_id` must agree with the row and request; conflicts fail before
  prepared-genome catalog access. When a legacy row has no context ids, the
  resolver creates deterministic report-local contexts from the requested
  species/genome pair. A preserved mapping with only a species declaration
  receives an organism-only context and does not fabricate a `genome_id`.
- Species aliases are normalized for matching. Ambiguous target mappings are
  unresolved by default; `ambiguity_policy=first` chooses the stable first
  candidate and records a warning. `ambiguity_policy=preserve` also leaves the
  target unresolved, but adds ordered `candidate_mappings[]` rows carrying
  each candidate's gene identity, genome/context references, orthology
  type/confidence, provider source, and evidence. Preserved candidates are not
  inserted into the resolved promoter `rows[]`, so downstream comparison does
  not silently treat alternatives as accepted cohort members. The shipped
  `reject` and `first` behavior remains unchanged. Candidate labels include
  provider source text when available, but labels and ranks are local to one
  report and must not be persisted as stable candidate identities. Ambiguity
  policy is a closed operation-control enum: unknown values are rejected
  rather than defaulted. The additive `preserve` value remains in v1 because
  it is opt-in and does not alter existing policy meanings; older v1 readers
  cannot read reports that explicitly select it.
- Resolved rows carry species, genome id, gene id/symbol, transcript id,
  strand, TSS, promoter span, transcription-aligned promoter sequence, and
  orthology evidence/provenance. They also refer to a context copied into the
  cohort report; target rows retain the mapping's oriented source/target
  context references. Unresolved rows make missing or ambiguous mappings
  explicit, and preserved ambiguity candidates refer only to contexts copied
  into the same portable report.
- Symbol-only mappings remain supported for legacy/local resources. A symbol
  is a lookup key, not evidence of orthology or functional equivalence; those
  claims remain bound to the row's declared type, source, and evidence.
- Ortholog promoter cohort and comparison reports may carry an additive
  `relationship` expectation (`manual`, `co_regulated`,
  `anti_co_regulated`, or `unspecified`) plus non-blocking
  `relationship_flags[]`. Co-regulated expectations flag unexpected TFBS or
  CUT&RUN divergence; anti-co-regulated expectations flag unexpected
  concordance. These rows are evidence triage only, not regulatory proof.
- `SummarizeOrthologPromoterComparison` keeps evidence channels separate:
  promoter-sequence identity, TFBS score-track similarity, motif peak
  presence, optional expression assignment, and CUT&RUN/occupancy status are
  reported in distinct arrays.
- CUT&RUN source ids can assign species/genome-matched occupancy support as
  `confirmed`, `nearby`, `motif_only`, `occupancy_only`, or `no_data`.
  Rows without matching provenance are `not_comparable`; raw peak intensity is
  never compared across species by default.
- `SummarizeOrthologPromoterComparison` accepts an optional additive
  `cutrun_normalization` record. It must name a normalization method, unit,
  shared comparison reference, provenance statement, and one unambiguous
  normalized value for every resolved promoter. Each value identifies its
  contributing selected CUT&RUN dataset/read-report ids and may add
  value-specific provenance.
- The optional `cutrun_quantitative_comparison` report section is
  `comparable` only when every normalized value matches one resolved promoter,
  every cited source was selected and evaluated for that species/genome, and
  all required normalization metadata is present. It then reports neutral
  pairwise values and `right - left` deltas without inferring activation or
  conservation.
- Selected CUT&RUN sources without that complete record produce
  `cutrun_quantitative_comparison.status=not_comparable` while retaining the
  qualitative support rows. GENtle does not derive cross-species numbers from
  BED scores, peak height, signal intensity, fragment counts, or read depth.

## TFBS score-track similarity contract

GENtle now also exposes one anchor-vs-candidate TFBS score-track similarity
report for the same selected DNA span.

First-class operation route:

```json
{"SummarizeTfbsTrackSimilarity":{"target":{"kind":"seq_id","seq_id":"SEQ_ID","span_start_0based":15564,"span_end_0based_exclusive":16764},"anchor_motif":"TP73","candidate_motifs":["ALL"],"ranking_metric":"smoothed_spearman","score_kind":"llr_background_tail_log10","clip_negative":true,"species_filters":["Homo sapiens"],"include_remote_metadata":true,"limit":25}}
```

Portable schema:

- `gentle.tfbs_track_similarity.v1`

Behavior notes:

- the report ranks candidate motifs against one selected anchor motif using the
  same displayed score-track vectors already used by the continuous TFBS score
  track route
- initial ranking metrics are:
  - `raw_pearson`
  - `smoothed_pearson`
  - `raw_spearman`
  - `smoothed_spearman`
- rows retain the full metric family plus signed primary-peak offset so
  adapters can later re-sort without recomputing biology
- `candidate_motifs=["ALL"]` or `["*"]` expands to the full local JASPAR
  registry
- otherwise `anchor_motif` and every `candidate_motifs[]` token resolve
  through the same shared TF-query layer used by score tracks and JASPAR
  summaries, so aliases/groups/family-like queries expand deterministically
- optional `species_filters[]` currently use the persisted local
  `jaspar.remote_metadata.json` snapshot rather than live network refresh:
  motifs without cached metadata are skipped when such a filter is active
- remote metadata remains additive:
  - the ranking can run without it
  - `include_remote_metadata=true` adds compact per-row summaries when cached
    rows exist locally
- unlike the pairwise synchrony heatmap, this route sorts by the selected
  metric itself rather than its absolute value because the intent is “most
  similarly scoring”, not “most strongly related in either direction”

## Promoter Evidence Matrix Contract

GENtle now exposes a conservative promoter-region evidence ledger for one
sequence/locus. This is deliberately not a final biological verdict; it records
promoter candidates and evidence items in a versioned structure that can later
carry orthology, co-regulation, anti-co-regulation, and CUT&RUN evidence.

Shared-shell route:

```bash
gentle_cli shell 'features promoter-evidence-matrix SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--no-feature-overlaps] [--path FILE.json]'
gentle_cli shell 'features promoter-isoform-comparison SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--no-feature-overlaps] [--path FILE.json]'
gentle_cli shell 'features promoter-expression-evidence SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--expression-json JSON] [--source-label LABEL] [--path FILE.json]'
gentle_cli shell 'features promoter-artifact-manifest SEQ_ID --artifact-json JSON --path FILE.json [--gene-label LABEL]'
```

First-class operation route:

```json
{"SummarizePromoterEvidenceMatrix":{"input":"SEQ_ID","gene_label":"TP73","promoter_upstream_bp":1000,"promoter_downstream_bp":200,"include_feature_overlaps":true,"path":"tp73_promoter_evidence.json"}}
```

```json
{"SummarizeIsoformPromoterComparison":{"input":"SEQ_ID","gene_label":"TP73","promoter_upstream_bp":1000,"promoter_downstream_bp":200,"include_feature_overlaps":true,"path":"tp73_isoform_promoter_comparison.json"}}
```

```json
{"SummarizePromoterExpressionEvidence":{"input":"SEQ_ID","gene_label":"TP73","promoter_upstream_bp":1000,"promoter_downstream_bp":200,"expression_source_label":"rna_seq_table_1","expression_rows":[{"transcript_id":"ENST...","sample_id":"case_1","condition":"case","value":12.4,"unit":"TPM","artifact_path":"rna_expression.tsv"}],"path":"tp73_promoter_expression.json"}}
```

```json
{"ExportPromoterArtifactManifest":{"input":"SEQ_ID","gene_label":"TP73","artifacts":[{"artifact_id":"evidence_matrix","artifact_kind":"promoter_evidence_matrix","path":"tp73_promoter_evidence.json","schema_hint":"gentle.promoter_evidence_matrix.v1","required":true}],"path":"tp73_promoter_artifacts.json"}}
```

Portable schema:

- `gentle.promoter_evidence_matrix.v1`
- `gentle.isoform_promoter_comparison.v1`
- `gentle.promoter_expression_evidence.v1`
- `gentle.promoter_artifact_manifest.v1`

Behavior notes:

- transcript-derived promoter windows reuse the same strand-aware TSS geometry
  and exact-span collapse as Alternative Promoter Comparison
- existing DNA-level `promoter` annotations are merged into the same candidate
  row when their span matches a transcript-derived promoter window
- evidence rows are string-kind based so unknown future evidence can be
  preserved by older adapters
- current evidence kinds include:
  - `promoter_geometry`
  - `transcript_support`
  - `promoter_annotation`
  - `tfbs_annotation`
  - `variant_overlap`
  - `repeat_context`
  - `enhancer_overlap`
  - `terminator_overlap`
  - `external_interval_overlap`
  - `cutrun_peak_overlap` when projected interval metadata identifies CUT&RUN
- the current `ranking_mode` is only deterministic display ordering:
  transcript support first, then evidence count, then coordinates and label
- richer scoring is expected to be additive and method/provenance-tagged rather
  than replacing the evidence ledger
- `SummarizeIsoformPromoterComparison` consumes the same evidence rows and
  reports common versus differential evidence signatures between promoter
  groups for isoforms of the same gene; geometry/support bookkeeping is kept
  out of the comparison signatures so differences focus on annotated/regulatory
  features
- `SummarizePromoterExpressionEvidence` accepts expression rows as external
  evidence with transcript/gene/promoter labels plus optional artifact paths.
  GENtle assigns them to promoter groups and reports unassigned rows, but it
  does not infer causality or perform normalization. The Promoter design GUI
  can paste or load the same row JSON, run this shared operation, cache the
  returned report, and export the same portable
  `gentle.promoter_expression_evidence.v1` payload; GUI consumers display
  promoter-associated evidence rather than a causal or validation verdict
- `ExportPromoterArtifactManifest` is a lightweight component index, not a
  narrative bundle. It records which JSON/SVG artifacts exist and whether
  required components are missing so downstream tools can compose their own
  presentation

## TF query resolution contract

GENtle also exposes one lightweight TF-query resolution report so CLI, GUI,
and ClawBio/OpenClaw users can audit what a natural-language transcription
factor query will map to before running TFBS/JASPAR analysis.

Current shared-shell route:

```bash
gentle_cli shell 'resources resolve-tf-query "Yamanaka factors" OCT4 "KLF family" --output tf_queries.json'
```

First-class operation route:

```json
{"ResolveTfQueries":{"queries":["Yamanaka factors","OCT4","KLF family"],"path":"tf_queries.json"}}
```

Portable schema:

- `gentle.tf_query_resolution.v1`

Behavior notes:

- each query row reports:
  - the original query string
  - whether it matched exactly, by alias, as a built-in group, or by a
    family-like local-registry expansion
  - the expanded motif ids/names GENtle will actually use downstream
- this is the intended bridge between user-facing names such as
  `Yamanaka factors`, `stemness`, or `OCT4` and the concrete local motif
  registry rows that power:
  - `SummarizeJasparEntries`
  - `SummarizeTfbsScoreTracks`
  - `SummarizeTfbsTrackSimilarity`
  - `ScanTfbsHits`

## Gene group catalog contract

Implemented baseline:

- schema: `gentle.gene_group_catalog.v1`
- purpose: deterministic, local, catalog-extensible gene groups that can map to
  public ontologies without being limited to them
- built-in catalog: `assets/gene_groups.json`
- default discovery: built-in `assets/gene_groups.json`, built-in
  `assets/gene_groups.d`, system/user/project `catalogs/gene_groups.json`, and
  system/user/project `catalogs/gene_groups.d`

Rationale:

- Gene Ontology and related public resources are valuable external mappings,
  but GENtle also needs lab-facing groups such as reprogramming-factor panels,
  clinical/research cohorts, family aliases, and project-local gene sets.
- The authoritative runtime input should be a reviewed catalog record, not a
  hidden agent prompt or one-off GUI list.
- AI-assisted generation is useful as a drafting path, but draft memberships
  must carry provenance and review state before downstream analyses treat them
  as trusted facts.
- Gene Ontology is represented as an external resource namespace (`GO`) in the
  catalog metadata, analogous to other GENtle external resources. A group may
  map to GO when useful, but the catalog can also hold lab terms without GO
  coverage.
- `resources status` also reports `gene_ontology` as a declared external
  database/mapping namespace. This is intentionally not a GO download/index
  implementation yet.

Expected record shape:

- `id`, `label`, `aliases[]`, `short_description`, `long_definition`
- `organism` / taxon and gene-symbol namespace
- optional reference-genome or transcript-source expectations
- `members[]` with canonical gene symbol/id, optional external identifiers,
  evidence notes, confidence/status, and provenance
- `external_mappings[]` for GO, Reactome, MSigDB, HPO, MeSH, UniProt keyword,
  and lab/project namespaces
- `curation_status` such as `draft`, `reviewed`, `curated`, or `deprecated`

Current shared-shell routes:

```bash
gentle_cli shell 'gene-groups list [--catalog PATH] [--filter TEXT] [--output OUTPUT.json]'
gentle_cli shell 'gene-groups show GROUP_ID [--catalog PATH] [--output OUTPUT.json]'
gentle_cli shell 'gene-groups resolve TOKEN [--catalog PATH] [--output OUTPUT.json]'
gentle_cli shell 'gene-groups doctor [--catalog PATH] [--output OUTPUT.json]'
gentle_cli shell 'gene-groups draft --description TEXT [--member SYMBOL] [--candidate SYMBOL=EVIDENCE] [--go GO:NNNNNNN] [--output GROUP.json]'
```

The initial built-in rows seed the previous hard-coded TF-query groups:

- `yamanaka_factors`
- `p53_family`
- `regulation_of_alternative_splicing`, a lab-facing group anchored to
  `GO:0000381` as an external reference while remaining reviewable local
  knowledge.

`resources resolve-tf-query` now consults catalog-backed groups marked with
`usages=["tf_query", ...]` before falling back to family-like motif-registry
matching. Existing built-in rows still report `resolution_kind=builtin_group`
for compatibility.

AI/user-assisted drafting is an explicit review-gated route, not a hidden
mutation of trusted catalogs:

```bash
gentle_cli shell 'gene-groups draft --description TEXT [--id ID] [--label LABEL] [--short-description TEXT] [--organism NAME] [--taxon-id N] [--namespace NAMESPACE] [--alias TEXT] [--tag TEXT] [--usage TEXT] [--member SYMBOL] [--members A,B,C] [--candidate SYMBOL=EVIDENCE] [--unresolved-candidate TEXT] [--go GO:NNNNNNN] [--agent-provider NAME] [--agent-model NAME] [--agent-generated-at UTC] [--provenance TEXT] [--output GROUP.json]'
```

The draft route writes a `gentle.gene_group_catalog.v1` catalog fragment plus a
`gentle.gene_group_draft.v1` report with `review_required=true`, the input
description hash, user/member counts, agent-suggested candidate members,
unresolved candidates, warnings, and the proposed local record. It does not
mutate trusted catalogs by default. A later import/promote step can require
explicit review before a drafted group becomes authoritative.

## Gene set analysis operands

Implemented baseline:

- schemas:
  - `gentle.gene_set_resolution.v1`
  - `gentle.gene_set_promoter_cohort.v1`
  - `gentle.gene_set_cutrun_regulatory_support.v1`
  - `gentle.gene_set_pool_creation.v1`
- engine operations:
  - `ResolveGeneSet`
  - `ProduceGeneSetDirectList`
  - `ProduceGeneSetOntologyAssignment`
  - `ProduceGeneSetCoRegulatedCohort`
  - `BuildGeneSetPromoterCohort`
  - `CreateGeneSetPool`
  - `InspectCutRunGeneSetRegulatorySupport`
- `gene-groups resolve` answers "which catalog entry matches this token";
  `gene-sets resolve` answers "which genes are in this analysis operand after
  expansion, gating, deduplication, and provenance recording".
- Choose `ExternalMapping` (including `--go`) for reviewed local catalog
  groups that cite a term; choose the ontology-assignment producer for genes
  assigned to that term by a specified offline provider/cache.
- V1 is offline-first: GO/external mappings resolve only from local
  `external_mappings`; no live GO download or Ensembl ortholog/paralog API call
  is made.
- Promoter cohorts may carry an additive declared relationship expectation
  (`unspecified`, `manual`, `co_regulated`, or `anti_co_regulated`).
  Downstream CUT&RUN support reports derive non-blocking occupancy
  relationship flags over evaluated members only; unevaluated members are not
  coerced into support or absence.
- Resolved gene-set reports, promoter cohorts, and gene-set CUT&RUN support
  reports are persisted as logical gene-set artifacts. Shared lineage SVG
  export renders them as `GeneSet`/analysis nodes linked from producer
  operation to gene set to downstream analysis, not as sequence or pool nodes.

Current shared-shell routes:

```bash
gentle_cli shell 'gene-sets resolve [GROUP_ID|--group GROUP_ID|--members A,B|--go GO:NNNNNNN|--neighbors GENE --flank-genes N|--random-size N --seed N] [--genome GENOME_ID] [--catalog PATH] [--genome-catalog PATH] [--allow-draft] [--allow-deprecated] [--output OUTPUT.json]'
gentle_cli shell 'gene-sets produce direct-list --cache CACHE.json_or_tsv [--query LIST_ID] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]'
gentle_cli shell 'gene-sets produce ontology-assignment --cache CACHE.json_or_tsv --term GO:NNNNNNN [--ontology-namespace GO] [--evidence-code CODE] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]'
gentle_cli shell 'gene-sets produce co-regulated --cache CACHE.json_or_tsv --dataset DATASET_ID --contrast LABEL --score METHOD --threshold RULE --direction both|positive|negative [--relationship co-regulated|anti-co-regulated|manual] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]'
gentle_cli shell 'gene-sets create-pool RESOLUTION_ID --member-container MEMBER_ID=CONTAINER_ID [--member-container MEMBER_ID=CONTAINER_ID ...] [--output-prefix PREFIX] [--container-name NAME] [--preview | --apply --expected-plan-fingerprint-sha256 SHA256] [--path REPORT.json]'
gentle_cli shell 'gene-sets promoter-cohort GENOME_ID [--resolution RESOLUTION.json|--group GROUP_ID|--members A,B|--go GO:NNNNNNN|--neighbors GENE --flank-genes N|--random-size N --seed N] [--relationship manual|co-regulated|anti-co-regulated] [--upstream-bp N] [--downstream-bp N] [--catalog PATH] [--genome-catalog PATH] [--output OUTPUT.json]'
```

Resolution notes:

- group-level `curation_status` and member-level `status` are independent
  gates; draft groups require `--allow-draft`, deprecated groups require
  `--allow-deprecated`, draft members are included with a warning, excluded
  members are skipped, and unknown member statuses are included with warnings.
- deterministic random sets use a stable prepared-genome universe and record
  genome id/build, gene-index source, seed, universe size, and foreground
  exclusion count.
- genomic neighbors are same-chromosome only and require a uniquely resolved
  anchor.
- duplicate collapse prefers resolved `gene_id`, then normalized symbol, and
  keeps merge provenance auditable.
- `--relationship` records a user-declared expectation (`manual`,
  `co_regulated`, or `anti_co_regulated`) on the cohort. It is not a computed
  verdict about regulation; downstream evidence reports can surface
  expectation-consistent or divergent CUT&RUN/TFBS signals as non-blocking
  flags.
- existing jq recipes that approximate membership with `.status // "included"`
  may silently drop or miss draft-member warnings; the engine now reports draft
  members explicitly in gene-set resolution warnings.

Biological-context binding:

- `gentle.gene_set_resolution.v1` owns a `contexts[]` registry plus an optional
  `default_context_id`. Each context can identify organism, taxon, GENtle
  genome-catalog entry, assembly accession, annotation source/release, and
  symbol namespace.
- A resolved member may select a registry row through `context_id`. Effective
  context resolution uses member context, then report default, then legacy
  report-level `genome_id`/organism/taxon/namespace fields. Lower-precedence
  fields may fill missing values but may not contradict an already selected
  context.
- Existing V1 reports remain readable. When they contain legacy report-level
  context fields, GENtle promotes those fields into one deterministic default
  registry row when the report is produced, stored, or read.
- Context-sensitive collection consumers must explicitly declare
  `context_requirement="homogeneous"` in the lift-policy registry. An absent
  declaration deserializes as `not_reviewed`, never as permission to treat
  mixed contexts as safe.
- The complete declarative vocabulary is `not_reviewed`, `context_agnostic`,
  `homogeneous`, `partitionable`, and `explicit_cross_context`. This slice
  implements the first two current behaviors (`context_agnostic` and
  `homogeneous`); partitioning and explicit cross-context comparison require
  their own consumers and are not inferred automatically.
- Promoter derivation and gene-set primer-specificity mapping require one
  resolvable homogeneous context whose `genome_id` equals the requested target
  genome. They reject before coordinate lookup or BLAST with
  `missing_biological_context` or `mixed_biological_context` otherwise.

Retrieval producer metadata:

- `gentle.gene_set_resolution.v1` carries additive, defaulted metadata for
  retrieval-backed producers that supply candidate members before resolution.
  Existing payloads without these fields remain valid and default to
  `review_status="unreviewed"` with no producer metadata.
- The resolver request remains the normalized analysis source. For example,
  `gene-sets produce direct-list` and `gene-sets produce ontology-assignment`
  resolve via `request.source_kind="explicit_members"` after cache lookup, while
  `producer`, `query_metadata`, and set-level
  `organism`/`taxon_id`/`symbol_namespace` fields record where those candidate
  members came from.
- Producer metadata is not a trust shortcut. Retrieved sets should carry
  provider/cache version, query/filter metadata, and review status before
  downstream analyses treat them as reviewed operands.
- `gene-sets produce direct-list` is the first offline retrieval producer. It
  reads a local JSON cache marked with
  `gentle.gene_set_direct_list_cache.v1`, or a TSV/CSV with `symbol` or
  `gene_id` rows and optional `list_id`, then resolves the selected members
  through the existing `explicit_members` resolver. The route requires
  provider provenance plus provider/cache version metadata and at least one
  interpretation context field (`organism`, `taxon_id`, or
  `symbol_namespace`) from the cache or CLI flags.
- `gene-sets produce ontology-assignment` is also offline-only. It reads local
  JSON caches marked with `gentle.gene_set_ontology_assignment_cache.v1` or
  TSV/CSV assignment rows, selects genes assigned to `--term` under optional
  evidence/filter constraints, and then resolves those candidate members
  through the existing `explicit_members` resolver. A term with no matching
  assignment rows returns a gene-set report with an unresolved term row and a
  warning, not a silent empty success.
- `gene-sets produce co-regulated` is the first evidence-derived retrieval
  producer. It reads local JSON caches marked with
  `gentle.gene_set_co_regulated_cache.v1` or TSV/CSV rows, filters rows by
  dataset/contrast/condition, applies an explicit numeric threshold rule
  (`score>=N`, `score>N`, `score<=N`, `score<N`, `abs>=N`, or `abs>N`) plus a
  sign-direction rule (`both`, `positive`, or `negative`), and resolves the
  selected members through `explicit_members`. The report carries
  `co_regulated_metadata` with dataset, normalization, scoring, threshold,
  direction, and declared relationship expectation. It is a retrieval result,
  not proof of regulation.
- Local cache import routes are explicit shell resource operations:
  `resources import-gene-list-cache`,
  `resources import-ontology-assignment-cache`, and
  `resources import-co-regulated-cache`. They normalize local TSV/CSV files
  into the corresponding offline cache schemas. They do not resolve members,
  mutate trusted gene-group catalogs, download live GO annotations, or infer
  regulation.
- Current producer families represented by the protocol metadata are
  `direct_gene_list`, `ontology_assignment`, and `co_regulated_cohort`. The
  ontology-assignment producer is distinct from existing local
  `external_mapping` resolution: `external_mapping` asks which local GENtle
  catalog groups cite a term, while ontology-assignment producer metadata
  records a provider/cache membership lookup for that term.
- Co-regulated producer metadata records dataset, contrast, normalization,
  scoring, threshold, sign-direction rule, and declared relationship
  expectation. It must not be read as proof of regulation; the default
  interpretation note says the evidence-derived cohort is a retrieval result
  and does not prove regulation.
- Gene-set lineage artifacts are graph-visible when reports are produced or
  promoter/CUT&RUN analyses are built. `render-lineage-svg` shows logical
  gene-set nodes with resolved/unresolved counts, producer kind, taxon, and
  namespace, then links those nodes to promoter-cohort and CUT&RUN support
  analysis nodes.

### Sequence collection subjects

Gene sets are one source of a broader protocol subject: an auditable collection
of sequence-like members. The collection subject is engine-owned so GUI,
CLI/shared-shell, MCP, JS, Lua, Python, and agent routes can expose the same
readiness, errors, results, and provenance instead of each adapter looping over
single-sequence behavior independently.

Detailed implementation plan:
[`gui_gene_set_collection_operations_plan.md`](gui_gene_set_collection_operations_plan.md).

Collection semantics stay explicit:

| Collection subject | Meaning | Typical producers | Typical projections |
|---|---|---|---|
| Logical set | Members plus provenance; no physical mixing implied | `GeneSetResolutionReport`, project selection, query result | per-member operations, multi-record export, promoter-window derivation |
| Physical pool/container | One sample/tube/lane may contain multiple molecules | `ExportPool`, operation-created product containers, imported pools | pool gel, storage, split/filter follow-ups |
| Arrangement | Ordered semantic layout across containers | `CreateArrangementSerial`, Gibson/bench handoff workflows | serial gel, rack placement, arrangement-scoped labels |
| Alignment | Ordered members plus aligned-column correspondence | future multiple-sequence alignment reports | alignment inspection, consensus/variant projection |
| Derived collection | New members derived from source members | promoter windows, restriction fragments, amplicons, neighboring loci | materialized sequences, downstream analysis sets |
| Storage collection | Physical placement of members or arrangements | rack/plate placement, future inventory/freezer views | labels, carrier templates, inventory reports |

Operation lifting rule:

- New or promoted operations that accept one sequence must state their
  collection behavior before gaining a prominent collection GUI affordance.
- Valid lifting modes are:
  - `map`: run the operation independently for each member and return
    per-member reports plus a collection summary.
  - `combine`: intentionally treat all members as one sample/pool.
  - `compare`: consume multiple members together, as in alignment or cohort
    comparison.
  - `arrange`: preserve or create explicit member order/placement.
  - `derive`: create one or more descendant members per source member while
    retaining source-member provenance.
  - `reject`: decline the collection operand with a typed reason when the
    operation is biologically or contractually single-sequence only.
- The selected lifting mode is part of the shared engine/shell contract, not a
  GUI implementation detail. GUI controls should collect the collection
  operand, show per-member and aggregate readiness/errors/results, and call the
  same operation path available to CLI/MCP/agent routes.

Initial high-value lifting expectations:

| Single-sequence operation family | Collection behavior |
|---|---|
| Sequence export/inspection | `map` per member or write one explicit multi-record artifact when requested |
| Restriction digest / PCR / primer design | `map` per member, with per-member products/reports and an aggregate warning summary |
| Pool/gel rendering | `combine` into one pool lane only when explicitly requested; otherwise `arrange` one lane per member/container |
| BLAST / feature scans / TFBS scans | `map` per member, preserving one result set per source member |
| Promoter or neighboring-sequence derivation | `derive` one or more windows per resolved member, preserving gene-set/source provenance |
| Multiple sequence alignment | `compare` members together and return an alignment report with member order and column correspondence |
| Rack/freezer/inventory placement | `arrange` into a storage projection from containers/arrangements, without changing logical set identity |

Implemented collection-lifting baseline:

- `gentle.collection_lift_policy_registry.v1` is loaded from
  `docs/collection_lift_policies.json`. Each curated row is keyed by capability
  source/name and collection subject kind. It declares either one supported
  lifting mode and result payload kind, or a typed rejection reason. The
  additive `context_requirement` field is `not_reviewed` by default and can
  require a homogeneous biological context for coordinate- or
  reference-sensitive operations.
- Capability descriptors project their applicable `collection_lift_policies`;
  adapters should consume that field rather than inventing local collection
  behavior.
- Dynamic readiness remains in the normal fact/precondition and operation
  error machinery. A temporary missing input is not encoded as a static
  collection-policy rejection.
- `RenderPoolGelSvg` makes the lane geometry explicit in that registry:
  physical containers are `combine`, serial arrangements are ordered
  `arrange`, and project selections plus logical gene sets are rejected with
  `requires_physical_pool`. A plate-mode arrangement is instead a dynamic
  serial-only readiness error. The existing shared shell operation is
  executable by arrangement id even though the generic collection launcher
  does not yet have a typed arrangement adapter or collection-report wrapper.
- `gentle.collection_operation.v1` records the selected subject, capability,
  lifting mode, policy, source-membership lock, per-member outcomes and typed
  errors, produced report ids, warnings, and provenance. `dry_run` and
  `applied` distinguish previews from committed operations. Because this report
  may outlive or be transported without its source report, it owns a portable
  copy of the source `contexts[]` registry and `default_context_id`.
- Domain reports that embed their source resolution by value, such as
  `GeneSetPromoterCohortReport`, inherit the registry through that embedded
  resolution and do not duplicate a second domain-level registry. Derivation
  preserves an existing source resolution `op_id` / `run_id`; only an
  anonymous inline resolution receives the deriving operation's identity when
  it is first persisted.
- Derived members are additional `per_member_status` rows. Their
  `parent_member_id` points to the source member, while
  `produced_report_ids` links both source and descendants to the domain report.
- `BuildGeneSetPromoterCohort` is the first proving consumer. Its normal
  promoter-cohort result embeds the generic collection report, and the same
  report is returned in `OpResult.collection_operation`.
- `AssessPrimerPairSpecificityCollection` is the first `map` consumer. It
  applies the existing `AssessPrimerPairSpecificity` interpretation separately
  to each resolved primer-design report, persists successful child reports,
  and returns their ids through the collection report. The aggregate wrapper
  does not reimplement hit, amplicon, completeness, or biological pass/fail
  logic. Its operation JSON is accepted through the normal CLI `op`/workflow,
  MCP `op`, JavaScript, and Lua paths; `collections run primer-specificity` is
  the shared-shell convenience form.
- `FindRestrictionSitesCollection` proves the same `map` contract for a second
  operation family. `collections run restriction-scan` calls the existing
  single-sequence scanner once per loaded project sequence or explicitly
  sequence-bound gene-set member, preserving sequence topology and effective
  enzyme disclosure. Its `gentle.collection_restriction_site_scan.v1` wrapper
  owns successful child reports because restriction scans have no persisted
  report store; generic `produced_report_ids` therefore remains empty.
- `InspectConstructReasoningCollection` maps fresh construct-reasoning
  inspection over project sequences or explicitly sequence-bound gene-set
  members. Its additive
  `gentle.collection_construct_reasoning_inspection.v1` wrapper owns each
  successful member's graph identity, input fingerprint, and recommended
  inspection actions. Graphs are assembled on a cloned engine and are not
  inserted into project metadata, so member rows declare
  `graph_persisted=false` and generic `produced_report_ids` remains empty.
  Project-sequence selections are context-agnostic; gene-set resolutions must
  have homogeneous biological context. Containers and arrangements are
  explicitly unsupported because physical occupancy and order do not alter
  independent per-sequence construct reasoning.
- Collection domain-result wrappers are additive wire contracts: A5
  restriction-scan, A6 TFBS-scan, and A7 digest reports default missing fields
  and ignore unknown future output fields. Their member-binding request
  records remain strict and reject unknown fields.
- Restriction scanning is biological-context agnostic because it has no
  genome-scoped parameter. Missing gene-set member bindings are typed dynamic
  member failures, while unknown/duplicate bindings and unknown enzymes reject
  the request. The membership fingerprint covers membership/order only, not
  mutable sequence content.
- `ScanTfbsHitsCollection` is the next `map` consumer. The shell convenience
  route is `collections run tfbs-scan`; both paths invoke `ScanTfbsHits` once
  per loaded project sequence or explicitly sequence-bound gene-set member.
  Motif tokens and all default/per-TF thresholds have the same meaning as the
  direct operation, and motifs are never silently defaulted.
- `gentle.collection_tfbs_hit_scan.v1` owns successful child
  `gentle.tfbs_hit_scan.v1` reports because those reports have no persisted
  store. It exposes requested motifs, canonical effective motif ids,
  thresholds, optional per-member cap, retained hit totals/counts by TF,
  incomplete member ids, truncated member count, and
  `aggregate_counts_complete`. Aggregate counts are complete only when every
  member succeeds without truncation and scans every effective motif; a child
  cap may stop motif iteration, so retained zero/count values must not imply
  complete absence when that flag is false.
- `DigestCollection` is the first mutating `map` lift. Project-sequence
  subjects resolve directly, while logical gene-set members require explicit
  `DigestCollectionMemberBinding` rows. The operation reuses the existing
  topology-aware `Digest` implementation, rejects unknown enzymes before
  member execution, and retains typed member failures for member-local digest
  errors.
- Digest convergence is guarded by deterministic fragment-count, round-count,
  repeated-state, and no-progress checks. No wall-clock deadline participates
  in direct or collection digest acceptance, so machine load cannot make the
  same member succeed directly and fail when mapped over a collection.
- Preview is the default for `collections run digest`. The portable
  `gentle.collection_digest.v1` report owns per-member source snapshot hashes,
  effective output prefixes, planned fragment ids, lengths, topology,
  fragment snapshot hashes, and materialization state. Aggregate fields expose
  planned/created counts, completeness, and incomplete member ids.
- `sha256_collection_digest_plan_v1` binds the collection subject and
  membership fingerprint, exact bound source snapshots, effective enzymes,
  normalized output prefix, deterministic sequential id reservation, planned
  fragment snapshots, and typed failed-member details. Apply requires that
  exact preview fingerprint and rejects drift before mutation. Applied
  fragments receive lineage only from their own source member; logical
  collection digestion does not create a physical container. Existing
  `DigestContainer` remains the explicit physical-container route, while
  generic `Digest` collection policies reject container and arrangement
  subjects.
- A logical gene-set member has no inherent primer-assay identity. Its
  `PrimerSpecificityCollectionMemberBinding` must therefore name the exact
  `stable_member_id` and persisted `primer_report_id`. Project-sequence
  collections may resolve that binding automatically only when exactly one
  primer-design report names the sequence as its template. Missing,
  ambiguous, mismatched, duplicate, and unknown bindings become typed member
  or request errors rather than symbol-based guesses.
- `CollectionMemberOutcome::Succeeded` records successful execution and
  persistence, not a biological specificity pass. A child report with
  `summary.status = fail|incomplete|not_assessed` remains a successful
  collection execution row and contributes an aggregate warning. This keeps
  infrastructure failures distinct from scientific conclusions.
- Logical gene sets are explicitly rejected for `ExportPool` and
  `RenderPoolGelSvg` with `requires_physical_pool`; resolving genes never
  silently asserts that their products occupy one tube or gel lane.
- `ExportPoolCollection` is the first physical-container `combine`. It accepts
  only `CollectionSubjectRef::Container`, requires context-agnostic policy and
  declared exhaustive contents, rejects in-silico `selection` containers,
  validates every member, then delegates one atomic write to the existing
  `ExportPool` implementation. Physical `singleton` and `pool` containers are
  eligible. Its additive
  `gentle.collection_pool_export.v1` report embeds the generic collection
  report, identifies the `gentle.pool.v1` artifact/path and source container,
  and links every successful member to the same wrapper report id. It does not
  mutate project state. Non-exclusive containers fail before writing because
  the current pool artifact has no field for incomplete physical membership.
- `CreateGeneSetPool` is the explicit logical-to-physical `combine`. It accepts
  only a complete `GeneSetResolution` and strict
  `GeneSetPoolMemberBinding` rows from every stable member id to one distinct,
  exclusive singleton source container. The context requirement is
  `context_agnostic`; the physical bindings, rather than gene labels, identify
  the material being pooled.
- Preview returns additive `gentle.gene_set_pool_creation.v1` with the generic
  collection report, source-resolution id, per-member source container and
  source sequence, source snapshot hash, reserved aliquot sequence id,
  reserved pool-container id, and `sha256_gene_set_pool_plan_v1` lock. Apply
  requires that exact lock, rejects all drift before mutation, derives one
  lineage-linked aliquot per source, and places only those aliquots in one new
  exclusive `ContainerKind::Pool`. Source containers and their latest mappings
  are retained. The report/journal owns gene-set provenance; the base
  `Container` wire shape is unchanged.

Collection membership fingerprints use
`sha256_canonical_collection_members_v1`. GENtle hashes the UTF-8 canonical
JSON projection of the source members:

- project sequence selections, containers, and resolved gene sets are set-like;
  stable member ids are sorted and de-duplicated
- arrangements are order-bearing; rows are sorted by numeric
  `ordering_index`, duplicates are retained, and each index travels with its
  stable member id
- the subject kind and order semantics are part of the hashed projection

This deliberately avoids lexicographic position ordering (`10` before `2`) and
ensures that swapping two arrangement lanes changes the membership lock.
The membership fingerprint intentionally covers membership and ordering only.
Two collections with the same stable member ids but different assemblies or
annotation releases therefore have the same membership fingerprint. Consumers
must inspect the biological-context registry; the membership digest must not be
used as a proxy for identical biological input.

## Stateless sequence-scan contract

Implemented additive contract:

- shared non-mutating sequence-inspection routes can accept either:
  - a stored project sequence, or
  - inline ASCII DNA text without first materializing a project-state record
- intended consumers:
  - GUI selection-first actions (`current selection`, `visible span`,
    `whole sequence`)
  - CLI/shared shell
  - JS/Lua/Python wrappers
  - MCP/agent wrappers including ClawBio/OpenClaw-style low-latency helpers

Implemented shared operand shape:

- `SequenceScanTarget`
- stored-sequence form:

```json
{
  "kind": "seq_id",
  "seq_id": "tp73_region",
  "span_start_0based": 700,
  "span_end_0based_exclusive": 1200
}
```

- inline-sequence form:

```json
{
  "kind": "inline_sequence",
  "sequence_text": "ACGTACGTACGT",
  "topology": "linear",
  "id_hint": "inline_tp73_window",
  "span_start_0based": 0,
  "span_end_0based_exclusive": 12
}
```

Semantics:

- inline `sequence_text` is normalized through the same parsing path used for
  in-memory DNA sequences
- `topology` matters for full-span scans where circular wrap-around changes hit
  detection
- span fields remain local to the supplied operand:
  - for `seq_id`, local sequence coordinates on the stored record
  - for inline text, coordinates on the submitted text after normalization
- inspection remains non-mutating unless the caller explicitly chooses a second
  promote/export/materialize step

Implemented first-class operation on top of that operand:

- `FindRestrictionSites { target, enzymes?, max_sites_per_enzyme?, include_cut_geometry?, path? }`
  - purpose:
    - non-mutating REBASE-backed scan for recognition sites and cleavage
      geometry directly on one operand
  - result schema:
    - `gentle.restriction_site_scan.v1`
  - report fields:
    - `target_kind`
    - `target_label`
    - `source_sequence_length_bp`
    - `scan_start_0based`
    - `scan_end_0based_exclusive`
    - `scan_length_bp`
    - `scan_topology`
    - `enzyme_filters`
    - `enzymes_scanned`
    - `matched_site_count`
    - `skipped_enzyme_names_due_to_max_sites`
    - `rows[]`
  - row fields:
    - `enzyme_name`
    - `recognition_sequence`
    - local scan coordinates:
      - `recognition_start_0based`
      - `recognition_end_0based_exclusive`
    - source-sequence coordinates:
      - `source_recognition_start_0based`
      - `source_recognition_end_0based_exclusive`
    - optional cut geometry:
      - `forward_cut_0based`
      - `reverse_cut_0based`
      - `opening_start_0based`
      - `opening_end_0based_exclusive`
      - source-space equivalents for the same cut/opening positions
    - `end_geometry`
  - `enzymes=[]` means:
    - use the current shared preferred restriction-enzyme list
    - if that list is empty, fall back to the default preferred enzyme set
- `ScanTfbsHits { target, motifs, min_llr_bits?, min_llr_quantile?, per_tf_thresholds?, max_hits?, path? }`
  - purpose:
    - non-mutating thresholded JASPAR/IUPAC hit scan directly on one operand
      without creating `TFBS` features by default
  - result schema:
    - `gentle.tfbs_hit_scan.v1`
  - report fields:
    - `target_kind`
    - `target_label`
    - `source_sequence_length_bp`
    - `scan_start_0based`
    - `scan_end_0based_exclusive`
    - `scan_length_bp`
    - `scan_topology`
    - `motifs_requested`
    - `motifs_scanned`
    - `default_min_llr_bits`
    - `default_min_llr_quantile`
    - `per_tf_thresholds`
    - `max_hits`
    - `truncated_at_max_hits`
    - `matched_hit_count`
    - `rows[]`
  - row fields:
    - `tf_id`
    - `tf_name`
    - `motif_consensus_iupac`
    - `motif_length_bp`
    - local scan coordinates:
      - `match_start_0based`
      - `match_end_0based_exclusive`
    - source-sequence coordinates:
      - `source_match_start_0based`
      - `source_match_end_0based_exclusive`
    - `wraps_origin`
    - `forward_strand`
    - `matched_sequence`
    - `llr_bits`
    - `llr_quantile`
    - `true_log_odds_bits`
    - `true_log_odds_quantile`
  - threshold semantics:
    - `min_llr_bits` and `min_llr_quantile` provide default gates across all
      motifs
    - `per_tf_thresholds[]` can override those defaults by token, resolved TF
      id, or resolved TF name
    - `motifs=["ALL"]` or `motifs=["*"]` expands to the full local motif
      registry before scanning

Implemented follow-up:

- `SummarizeTfbsScoreTracks` and `RenderTfbsScoreTracksSvg` now accept the same
  `SequenceScanTarget` operand used by stateless restriction-site and TFBS-hit
  scans

UX parity expectations:

- GUI:
  - DNA-window actions can now call the shared restriction-site scan for:
    - current selection
    - current visible span
    - whole active sequence
- CLI/shared shell:
  - one command family accepts either stored `SEQ_ID` or inline
    `--sequence-text`
- testing:
  - deterministic parity tests compare inline-target results against stored
    `seq_id` results for the same sequence/span

## JASPAR entry presentation contract

GENtle also exposes a portable JASPAR-entry presentation contract that answers
three practical questions for one local motif snapshot:

- what sequence should maximize this motif score,
- what sequence should minimize it,
- and what score range should one expect on one deterministic random DNA
  background.

Current shared-shell route:

```bash
gentle_cli shell 'resources summarize-jaspar --motif SP1 --motif REST --random-length 10000 --seed 123 --output jaspar.summary.json'
```

First-class operation route:

```json
{"SummarizeJasparEntries":{"motifs":["SP1","REST"],"random_sequence_length_bp":10000,"random_seed":123,"path":"jaspar.summary.json"}}
```

Portable schema:

- `gentle.jaspar_entry_presentation.v1`

Behavior notes:

- When `motifs` is omitted or empty, GENtle summarizes every entry currently
  visible through the active local JASPAR registry.
- `motifs[]` accepts the same query forms described above for
  `ResolveTfQueries`, so callers can use aliases, built-in groups, or
  family-like queries instead of only exact motif ids.
- `maximizing_sequence` and `minimizing_sequence` are derived column-wise from
  the same local matrix/consensus information GENtle already uses for TFBS
  scoring; they are not hand-curated examples.
- The background distribution uses one deterministic pseudorandom
  `uniform_acgt_lcg` DNA sequence so shell/CLI/agent runs are reproducible.
- Per-entry output includes:
  - local entry identity (`motif_id`, optional `motif_name`,
    `consensus_iupac`)
  - motif length
  - maximizing/minimizing sequences
  - their exact `llr_bits` and `true_log_odds_bits`
  - their empirical quantiles against the same random-background distribution
  - min/max/mean/stddev plus `p01/p05/p25/p50/p75/p95/p99` summaries for both
    bit-based score families

## JASPAR registry benchmark contract

GENtle also exposes one routine registry-wide JASPAR benchmark artifact for
drift/regression work.

Current shared-shell route:

```bash
gentle_cli shell 'resources benchmark-jaspar --random-length 10000 --seed 123 --output data/resources/jaspar.registry_benchmark.json'
```

First-class operation route:

```json
{"BenchmarkJasparRegistry":{"random_sequence_length_bp":10000,"random_seed":123,"path":"data/resources/jaspar.registry_benchmark.json"}}
```

Portable schema:

- `gentle.jaspar_registry_benchmark.v1`

Behavior notes:

- this route always benchmarks the full active local JASPAR registry
- it reuses the same deterministic internal math as
  `SummarizeJasparEntries`, rather than introducing a separate
  benchmark-only scoring path
- the report includes:
  - all per-entry presentation rows
  - aggregate score-family summaries for `llr_bits` and
    `true_log_odds_bits`
  - top motifs by maximum score
  - top motifs by positive random-background hit rate
- the intended artifact is a cached/exportable JSON snapshot suitable for
  future scoring-drift comparisons

## JASPAR catalog contract

GENtle also exposes a lighter-weight portable JASPAR catalog report for
registry browsing without having to compute full per-entry score
distributions:

- the active local JASPAR snapshot remains the source of truth for which motif
  entries exist
- each catalog row includes:
  - motif id
  - optional TF name
  - consensus IUPAC sequence
  - motif length
  - and optionally a compact remote metadata summary
- when remote enrichment is requested, GENtle adds compact collection/class/
  family/tax-group fields plus species-count previews for only the returned
  subset instead of requiring every consumer to hand-roll JASPAR REST lookups

Current shared-shell route:

```bash
gentle_cli shell 'resources list-jaspar --filter TP --limit 50 --output jaspar.catalog.json'
```

First-class operation route:

```json
{"ListJasparCatalog":{"filter":"TP","limit":50,"include_remote_metadata":false,"path":"jaspar.catalog.json"}}
```

Portable schema:

- `gentle.jaspar_catalog.v1`

Behavior notes:

- when remote metadata is requested for the catalog, GENtle first reuses any
  persisted JASPAR remote-metadata snapshot rows that already exist locally
- explicit refresh is a separate operation so GUI/CLI/agent consumers can
  choose between:
  - cache-first registry browsing
  - and deliberate network refresh/persistence for selected subsets

## JASPAR remote-metadata snapshot contract

GENtle also exposes one persisted reusable JASPAR remote-metadata snapshot so
species/collection/class/family assignments can survive across sessions instead
of living only in one GUI/app run:

- the local motif snapshot still governs which entries exist
- the remote snapshot is an additive enrichment keyed by motif id
- each persisted row stores:
  - motif id / optional TF name
  - consensus IUPAC sequence
  - motif length
  - full remote metadata payload
  - and one compact derived summary for catalog/list views
- the default runtime snapshot path is:
  - `data/resources/jaspar.remote_metadata.json`
- `resources status` surfaces this snapshot separately from the local motif
  registry status so operators can tell whether cached species metadata exists,
  is valid, and how many persisted rows it currently holds

Current shared-shell route:

```bash
gentle_cli shell 'resources sync-jaspar-remote-metadata --filter TP --limit 50 --output jaspar.remote_metadata.json'
```

First-class operation route:

```json
{"SyncJasparRemoteMetadata":{"filter":"TP","limit":50,"path":"jaspar.remote_metadata.json"}}
```

Portable schema:

- `gentle.jaspar_remote_metadata_snapshot.v1`

## UCSC rmsk resource contract

GENtle now treats the UCSC RepeatMasker `rmsk` table as a portable external
resource rather than a GUI-only track import convention.

Current shared-shell routes:

```bash
gentle_cli shell 'resources install-ucsc-rmsk --assembly hg38'
gentle_cli shell 'resources sync-ucsc-rmsk https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz data/resources/ucsc.rmsk.hg38.json --assembly hg38'
gentle_cli shell 'resources prepare-ucsc-rmsk-index data/resources/ucsc.rmsk.hg38.json data/resources/ucsc.rmsk.hg38.interval-index.json'
gentle_cli shell 'resources suggest-ucsc-rmsk-index --assembly hg38 --output rmsk.indexing.json'
```

Portable schemas:

- `gentle.ucsc_rmsk_resource.v1`
- `gentle.ucsc_rmsk_install_report.v1`
- `gentle.ucsc_rmsk_interval_index.v1`
- `gentle.ucsc_rmsk_descriptor.v1`

Behavior notes:

- the normalized snapshot preserves the UCSC 17-column `rmsk` table:
  `bin`, `swScore`, `milliDiv`, `milliDel`, `milliIns`, `genoName`,
  `genoStart`, `genoEnd`, `genoLeft`, `strand`, `repName`, `repClass`,
  `repFamily`, `repStart`, `repEnd`, `repLeft`, `id`
- genomic coordinates remain UCSC 0-based half-open in the resource; callers
  convert to GENtle feature coordinates only when materializing display
  features
- snapshots include source URLs, field specs, and index recommendations so
  future GUI/CLI/projectors can inspect the resource without hard-coding UCSC
  table details
- `resources status` reports the default hg38 runtime snapshot path
  `data/resources/ucsc.rmsk.hg38.json`, the matching interval-index path,
  validation state, row/index counts when available, and the same index
  recommendations
- `resources install-ucsc-rmsk` is the preferred runtime setup route: it
  downloads or reads the source table, writes the normalized snapshot, and builds
  the interval index in one reportable step
- `resources prepare-ucsc-rmsk-index` builds the prepared interval sidecar at
  `data/resources/ucsc.rmsk.<assembly>.interval-index.json` by default; it
  refuses truncated snapshots so fixtures cannot accidentally masquerade as a
  full assembly index
- prepared interval indexes carry chromosome aliases such as `1` <-> `chr1`, so
  genome-anchor projection does not depend on catalog/source naming style
- `--limit N` creates a deliberately truncated snapshot for fixtures/smoke
  checks and marks `truncated=true`

Indexing guidance:

- primary: build a per-assembly interval index keyed by
  `(genoName, bin, genoStart, genoEnd)` so genome-region/gene extraction can
  project only overlapping repeats; GENtle's prepared sidecar stores intervals
  grouped by chromosome and sorted by `(start,end,row_offset)` for deterministic
  overlap queries, plus an alias dictionary for common contig spelling variants
- secondary: keep class/family partitions keyed by normalized
  `(repClass, repFamily, repName)` for display filtering, feature-tree
  grouping, and legends
- optional expert lookup: add a `repName -> row offsets` dictionary for queries
  such as `AluY` or `L1PA2`
- optional display-quality sidecar: store `swScore`, `milliDiv`, `milliDel`,
  and `milliIns` columnar by row offset so divergence/quality shading can be
  added without changing the interval index

## Publication dataset resource contract

GENtle now has a small built-in catalog for publication-associated external
datasets whose raw files should not be committed to the repository but should
be easy to inspect and fetch reproducibly.

Current shared-shell routes:

```bash
gentle_cli shell 'resources list-publication-datasets --filter p73'
gentle_cli shell 'resources status-publication-dataset E-MTAB-14704'
gentle_cli shell 'resources prepare-publication-dataset E-MTAB-14704'
gentle_cli shell 'resources prepare-publication-dataset PXD058816 --download-files --max-files 1'
```

Portable schemas:

- `gentle.publication_dataset_list.v1`
- `gentle.publication_dataset_status.v1`
- `gentle.publication_dataset_prepare.v1`
- `gentle.publication_dataset_download_manifest.v1`

Behavior notes:

- default catalog: `assets/publication_resources.json`
- default prepared-manifest/cache root: `data/publication_resources`
- the built-in catalog ships with checked-in, no-byte-download prepared
  manifests/download scripts for the current Rostock p73 entries
- `resources status` includes a compact `publication_datasets` block that
  reports catalog validity, dataset count, declared file count, and how many
  datasets already have prepared manifests
  - `resources prepare-publication-dataset` is manifest-first:
    - always writes `manifest.json`, `download_manifest.tsv`, and `download.sh`
      with declared file sizes and optional archive checksums when the catalog
      records them
    - does not download raw files unless `--download-files` is explicit
  - accepts `--max-files N` for smoke checks and small staged fetches
  - accepts `--category NAME` or `--categories CSV` to restrict planned files
    before download; this is the preferred way to fetch Clariom D CEL files
    without touching multi-GB PRIDE or unresolved CUT&RUN payloads
- built-in Rostock p73 multimodal entries:
  - `E-MTAB-14704`: Affymetrix Clariom D transcriptome CEL files and MAGE-TAB
    metadata; concrete BioStudies FTP file URLs are declared
  - `E-MTAB-15709`: CUT&RUN sequencing accession; resolved through SRA/ENA as
    `PRJEB100610` / `ERP182066`, with concrete paired FASTQ URLs and byte counts
    declared while keeping full raw-read download explicit
  - `PXD058816`: PRIDE Co-IP/MS raw/search files; concrete PRIDE archive URLs
    are declared
- entries with no concrete file URLs still produce a metadata manifest and
  warning, so agents can report unresolved accessions cleanly and refresh the
  catalog later instead of failing with an opaque missing-path error
- `scripts/analyze_p73_clariomd_probe_level.R` is the first external-analysis
  bridge for `E-MTAB-14704`; it keeps R/Bioconductor use outside the core
  engine while producing deterministic TSV outputs suitable for later GENtle
  panels. It expects `oligo`, `limma`, and the Clariom D platform design
  package `pd.clariom.d.human`.
- `scripts/probe_regions_oligo.R` is the generic R/oligo backend helper for
  `arrays probe-regions`; it consumes explicit CEL paths, optional TSV/CSV/SDRF
  metadata, selectors, a platform design package, and repeatable
  `--r-library-path PATH` values, then writes
  chromosome-ordered intensity CSVs plus expression/feature TSVs, limma
  contrast TSVs, provenance JSON, and a normalized matrix manifest. It
  currently supports `--normalization rma`; the Rust preflight emits an advisory
  command only for compatible explicit CEL requests, and users still run it
  explicitly. The helper records `r_version`, exact `package_versions`,
  requested and effective R library paths, `analysis_method_version`, and input
  fingerprints in provenance.
  `arrays inspect-probe-region-output` exposes these fields through the existing
  `gentle.probe_region_output_inspection.v1` contract. GENtle reports missing
  R/Bioconductor dependencies but never installs them. The `r_oligo` contract
  checks the helper's direct packages `oligo`, `limma`, `Biobase`, `DBI`, and
  `RSQLite`, plus the selected platform package, in one bounded, non-loading R
  probe and records every detected package version. If an agent-local package
  installation and a system installation disagree, the diagnostic lists the
  exact R library paths checked and explicitly asks the user to verify
  `--r-library-path PATH`. A direct R invocation
  records input path, size, and modification time; execution through
  `arrays run-probe-region-backend` enriches the finalized provenance with
  Rust-computed SHA-256 input fingerprints. The selected pdInfo package is
  installed outside GENtle and normally bundles its own
  `extdata/<package-name>.sqlite` database. The helper resolves that
  package-specific database, accepts a differently named file only when it is
  the sole package-local SQLite candidate, and refuses probeset output when no
  unambiguous database is available. Transcript-cluster-only output may
  continue from `netaffxTranscript.rda`, but emits an explicit warning instead
  of silently implying that probeset coordinates were loaded.
- `scripts/probe_regions_affy.R` is the matching explicit R/`affy` helper for
  legacy 3' IVT/CDF arrays. It consumes explicit CEL paths, a local CDF package
  or CDF name, optional metadata, optional Bioconductor annotation package, and
  optional user-supplied probeset coordinate table, accepts the same repeatable
  R-library-path flag, and writes the same helper-output contract including
  package versions and effective library paths. A CDF supplies probe grouping for RMA; genome
  placement for region displays still requires local annotation coordinates.
- Affymetrix/Thermo Fisher platform knowledge for `arrays probe-regions` is
  now a resource specification:
  `data/resources/affymetrix/platform_registry.json`
  (`gentle.affymetrix_platform_registry.v1`). The registry records aliases,
  family (`clariom`, `gene_st`, `exon_st`, `legacy_3prime_ivt_cdf`), species,
  genome-build hints, backend candidates, Bioconductor/CDF package names, and
  manually staged support-file expectations. Verified entries, currently led by
  Clariom D Human na36 hg38, may report concrete vendor ZIP paths; provisional
  historical entries recognize the chip family and CDF/backend needs while
  explicitly requiring local verification before external execution.
- Thermo Fisher Clariom D Human na36 hg38 probeset/transcript support ZIPs are
  login-walled vendor files. GENtle never downloads them automatically; users
  place them manually under
  `data/resources/affymetrix/clariom_d_human_na36_hg38/` as documented in that
  directory's README. Both concise canonical ZIP names and browser-preserved
  `TFS-Assets_LSG_Support-Files_...` download names are accepted through the
  Affymetrix platform registry, and
  `arrays probe-regions --platform Clariom_D_Human` reports their local status
  under `annotation_source.vendor_support_files[]`.

## RNA structure executable resources

`resources status` and `services status` also expose the executable resources
used by RNA secondary-structure rendering so ClawBio and other automation
layers can distinguish missing tools from missing biological data snapshots.

Resource fields:

- `vienna_rna`
  - display name: `ViennaRNA RNAfold`
  - resolves `GENTLE_RNAFOLD_BIN`, then `RNAfold` on `PATH`
  - version probe: `RNAfold --version`
  - used to compute dot-bracket structure and MFE
- `rnapkin`
  - display name: `rnapkin RNA structure renderer`
  - resolves `GENTLE_RNAPKIN_BIN`, then `rnapkin` on `PATH`
  - version probe: `rnapkin --version`, falling back to `rnapkin -V`
  - used to render RNAfold structure input as SVG/PNG
- `legacy_sha1`
  - display name: `Legacy SHA-1 verification tool`
  - resolves `GENTLE_SHA1_TOOL` first, then platform defaults
    (`sha1sum`, `shasum -a 1`, and Windows `certutil`)
  - `GENTLE_DISABLE_LEGACY_SHA1=1` reports `support_status=disabled`
  - used only for legacy upstream download/manifest SHA-1 fields; if missing
    or disabled, SHA-1 verification is skipped and basic corruption, size,
    gzip, or parser checks are the remaining safeguards

Most executable resources report `support_status`, `available`,
`resolved_executable`, `version_command`, optional `version_output`, optional
`error`, and explanatory `notes[]`. In service handoff/readiness views these
appear as `external_tool:*` rows rather than as normalized data-snapshot rows.
The `legacy_sha1` row is a compatibility verifier rather than a versioned
program dependency; it reports `env_var`, `disable_env_var`,
`resolved_executable`, `tool_kind`, `probe_digest`, `available`, `error`, and
`notes[]`.

## External service contracts

GENtle exposes vendor/CRO integrations through provider-neutral records owned
by the engine. Provider behavior is loaded from overlay-discoverable config
catalogs so GUI/CLI/MCP and ClawBio callers depend on shared schemas rather
than GeneArt-, metabion-, or portal-specific fields.

- `gentle.external_service_provider_catalog.v1`
  - returned by `services providers list` and embedded in `services status`
    under `external_providers`
  - top-level fields: `schema`, `generated_at_unix_ms`, `providers[]`,
    `summary_lines[]`
  - each provider row includes `provider`, `display_name`, `support_status`,
    `website_url`, `dashboard_url`, optional `api_documentation_url`,
    `capabilities[]`, `account_enablement_notes[]`, and `warnings[]`
  - each capability includes `service_kind`, `track`, `display_name`,
    `quote_handoff_supported`, `direct_api_documented`,
    `direct_api_implemented`, `supported_submission_modes[]`,
    `status_tracking`, `artifact_kinds[]`, and `notes[]`
- `gentle.external_service_provider_config.v1`
  - loaded by the shared services layer to build provider catalog/preflight/
    quote behavior
  - discovery order follows the catalog overlay policy:
    `assets/external_service_providers.json`, `assets/external_service_providers.d/`,
    `/etc/gentle/catalogs/external_service_providers.json`, matching system
    `.d/`, `~/.config/gentle/catalogs/external_service_providers.json`,
    matching user `.d/`, and project-local
    `.gentle/catalogs/external_service_providers.json` plus `.d/`
  - later sources override earlier provider ids; the doctor route reports
    source provenance and schema/validation errors
  - provider rows include metadata/capabilities, channel definitions such as
    `wop`, `email_excel`, or `oci_contact`, product-template mappings,
    validation rules, default delivery/purification/QC hints, and required
    follow-up policies
- `gentle.external_service_provider_config_doctor.v1`
  - returned by `services providers doctor [--catalog PATH] [--output PATH]`
  - reports inspected sources, parsed source/provider counts, SHA1 checksums,
    warnings, and deterministic validation errors
- `gentle.external_service_delivery_route_request.v1`
  - accepted by `services delivery-route REQUEST_JSON_OR_@FILE`
  - also produced and immediately classified by
    `services route-project-source --kind sequence|oligo-form|primer-report-rows ...`
  - fields: `schema`, `source_target`, optional `optimization_target`,
    optional `vector_spec`, optional `delivery_options`, optional
    `commercial_context_ref`, `return_spec`, and optional `request_metadata`
  - this is the provider-neutral input for generic wording such as "deliver
    this sequence"; callers should fill it from the active sequence or selected
    line items before asking GENtle to choose a provider route
  - `route-project-source` v1 supports selected project sequence/span,
    persisted oligo order forms, and primer report pair ranks; primer rows are
    always persisted as an oligo order form before routing so duplicate-review
    provenance is not bypassed
- `gentle.external_service_delivery_route.v1`
  - returned by `services delivery-route`
  - fields include `status`, `molecule_type`, `sequence_kind`,
    `sequence_count`, optional sequence-length fields,
    `recommended_provider`, `recommended_service_kind`, `candidates[]`,
    `summary_lines[]`, `rationale[]`, `clarification_questions[]`, and
    `warnings[]`
  - `candidates[].request` is a concrete
    `gentle.external_service_request.v1` that can be passed to
    `services project-preflight` or `services project-quote` after human review
  - current routing treats short DNA oligos/primers/probes as Metabion
    `dna_oligo_single_tube`, medium DNA fragments as Metabion `dna_fragment`,
    cloned/vector-backed DNA as GeneArt `cloned_gene`, long synthetic DNA as
    GeneArt `dna_fragment`, and protein products as GeneArt
    `protein_expression`; RNA or unclear molecule context asks for
    clarification rather than guessing
- `gentle.external_service_request.v1`
  - accepted by `services project-preflight REQUEST_JSON_OR_@FILE` and
    `services project-quote REQUEST_JSON_OR_@FILE [--output-dir DIR]`
  - fields: `schema`, `provider`, `service_kind`, `source_target`,
    optional `optimization_target`, optional `vector_spec`, optional
    `delivery_options`, optional `commercial_context_ref`, `return_spec`,
    and optional `request_metadata`
  - `commercial_context_ref` is a caller/session reference only; GENtle does
    not persist PO, shipping, credential, or account-secret values in project
    state
  - `return_spec` contains `requested_payloads[]`, optional
    `inline_max_bytes`, `redact_commercial_fields`, and
    `prefer_artifact_bundle` so ClawBio/MCP can request only the payload it
    needs, such as an amino-acid sequence, adjusted GenBank, quote metadata,
    or a full handoff bundle
- `gentle.external_service_preflight.v1`
  - returned by `services project-preflight`
  - fields include provider/service capability status, `eligible`,
    `quote_handoff_available`, `direct_submission_available`,
    `supported_submission_modes[]`, `blocking_issues[]`, `warnings[]`,
    optional turnaround/cost hints, `required_followup[]`, `dashboard_links[]`,
    and `request_summary[]`
  - in the first GeneArt slice this is local and deterministic; it does not
    contact Thermo Fisher/GeneArt
  - in the metabion slice this is also local and deterministic; WOP/email/Excel
    routes are represented as handoff channels, not automated submission
  - requests carrying `source_target.duplicate_review.status =
    review_required` are preflight-blocked until the source oligo order form is
    reviewed with `primers oligo-order review-dedup FORM_ID`; delivery-route
    classification remains allowed before that review
- `gentle.external_service_quote.v1`
  - returned by `services project-quote`
  - fields include `quote_status`, `quote_mode`, the embedded preflight
    report, `dashboard_links[]`, `required_followup[]`,
    `service_ready_bundle`, `return_spec`, and `warnings[]`
  - `service_ready_bundle.schema` is
    `gentle.external_service_artifact_bundle.v1`; V1 uses inline payloads for
    human-reviewable handoff markdown and redacted request JSON, and configured
    providers may also add normalized line-item JSON/CSV, email-draft markdown,
    and guided portal/WOP checklists
  - no vendor order, cart submission, direct API upload, status polling, or QAD
    retrieval is performed until a later confirmed direct-integration phase
  - Metabion v1 supports `provider = metabion` for DNA single-tube oligos and
    m-block DNA fragments/libraries as quote/handoff artifacts only; template
    URLs are advisory, and missing local templates produce warnings rather than
    hard quote-handoff failures
  - when `--output-dir DIR` is supplied, inline handoff payloads are also
    materialized as deterministic local files under `DIR`; the returned quote
    report records those generated files in
    `service_ready_bundle.local_files[]` and writes a complete
    `quote_report.json` into the same bundle directory

## JASPAR expert contract

GENtle also exposes a portable single-entry JASPAR expert contract for the
GUI, shell/CLI, and agent surfaces when they need a richer inspection view than
the registry-wide summary:

- which local entries are available is still governed by the active local
  JASPAR snapshot
- one selected entry can be expanded into:
  - its count matrix
  - a simple sequence-logo payload derived from per-column base frequencies and
    information content
  - multiple score-family distributions (`llr_bits`, `true_log_odds_bits`)
  - optional remote JASPAR metadata such as assigned species, collection, tax
    group, class, and family

Current shared-shell route:

```bash
gentle_cli shell 'resources inspect-jaspar SP1 --random-length 10000 --seed 123 --fetch-remote --output jaspar.expert.json'
```

First-class operation route:

```json
{"InspectJasparEntry":{"motif":"SP1","random_sequence_length_bp":10000,"random_seed":123,"include_remote_metadata":true,"refresh_remote_metadata":true,"path":"jaspar.expert.json"}}
```

Portable schema:

- `gentle.jaspar_entry_expert.v1`

Behavior notes:

- remote metadata enrichment is optional and best-effort; the expert still
  works from the local registry when offline or when the REST request fails
- when remote metadata is requested, the expert first reuses any matching row
  already persisted in `jaspar.remote_metadata.json` before falling back to a
  live refresh
- score panels include compact histogram bins so GUI/agent consumers can show
  actual distribution shape instead of only percentiles
- count-matrix and logo payloads come from the same local matrix counts GENtle
  uses for motif scoring, rather than a separate presentation-only transform

## Draft design resources

### `gentle.gibson_assembly_plan.v1`

Purpose:

- describe one Gibson cloning project in a destination-first way,
- separate user-specified plan inputs from derived design consequences,
- provide one canonical JSON artifact that future routines, primer design, and
  protocol-cartoon rendering can all read from.

Status:

- destination-first single-insert and ordered multi-insert plans are now
  consumed by the shared Gibson preview/apply path (`gibson preview ...`,
  `gibson apply`, and the `Patterns -> Gibson...` specialist window)
- current limit:
  - multi-insert execution currently assumes a defined destination opening
  - `existing_termini` remains the single-fragment handoff path for now

Canonical examples:

- `docs/examples/plans/gibson_destination_first_single_insert.json`
- `docs/examples/plans/gibson_destination_first_multi_insert.json`

Top-level structure:

- `schema`, `id`, `title`, `summary`
- `destination`
  - destination molecule (`seq_id`, prior topology)
  - explicit opening definition (`mode`, label, resulting left/right ends)
- `product`
  - intended output topology and output-id hint
- `fragments[]`
  - participating inserts or non-destination fragments
  - orientation plus per-end adaptation strategy
  - optional source-coordinate hints for deterministic ordering
- `assembly_order[]`
  - explicit left-to-right order of destination ends and inserts
  - supports future multi-fragment Gibson plans without changing the model
- `junctions[]`
  - one record per adjacent join
  - required overlap length
  - explicit overlap partition across the left/right adjacent members
  - whether overlap is derived from destination context or user-specified
  - explicit `distinct_from` constraints for terminal junctions
- `validation_policy`
  - hard requirements:
    - unique/unambiguous destination opening
    - distinct terminal junctions
    - adjacency-consistent overlaps
  - advisory checks:
    - overlap-length design range
    - fragment-count-aware overlap targets
    - overlap Tm
    - destination/fragment/reference uniqueness heuristics
  - optional design request:
    - `desired_unique_restriction_site_enzyme_name`
      asks the shared preview/apply path to try introducing one new unique
      REBASE cutter site on one terminal overlap if the assembled product can
      still remain uniquely cut there
- `derived_design`
  - derived overlap sequences
  - primer design suggestions
  - advisory notes and validation outcomes

Current draft value vocabulary:

- `destination.topology_before_opening`
  - `linear`
  - `circular`
- `destination.opening.mode`
  - `existing_termini`
    - use the current termini of an already-linear destination sequence
  - `defined_site`
    - a user-selected opening site/window on an existing destination molecule
  - reserved future values:
    - `restriction_digest`
    - `pcr_linearization`
    - `inverse_pcr`
- `destination.opening.uniqueness_requirement`
  - `must_be_unambiguous`
    - opening ambiguity is a hard validation error
  - `advisory_only`
    - ambiguity is surfaced, but not automatically fatal
- `product.topology`
  - `linear`
  - `circular`
- `fragments[].role`
  - `insert`
  - reserved future values:
    - `backbone_fragment`
    - `bridge_fragment`
- `fragments[].orientation`
  - `forward`
  - `reverse`
- `fragments[].source_span_1based`
  - optional source-coordinate hint for plans that preserve source order
  - shape:
    - `source_seq_id`
    - `start`
    - `end`
- `fragments[].left_end_strategy.mode` / `right_end_strategy.mode`
  - `native_overlap`
    - fragment terminus is already expected to satisfy the required overlap
  - `primer_added_overlap`
    - overlap is expected to be introduced via a primer tail
    - this is the current draft bucket for overlap-extension / primer-stitching
      style adaptation
  - reserved future values:
    - `synthetic_terminal_sequence`
    - `library_defined_overlap`
- `assembly_order[].kind`
  - `destination_end`
  - `fragment`
- `junctions[].overlap_source`
  - `derive_from_destination_left_flank`
  - `derive_from_destination_right_flank`
  - `derive_from_adjacent_fragment_ends`
  - `designed_bridge_sequence`
    - internal junction overlap chosen as a synthetic bridge/adaptor sequence
  - reserved future value:
    - `user_specified_sequence`
- `junctions[].overlap_partition`
  - explicit contribution of the overlap region from the adjacent members
  - shape:
    - `left_member_bp`
    - `right_member_bp`
  - invariant:
    - `left_member_bp + right_member_bp == required_overlap_bp`
  - examples:
    - left-member only overlap: `30 + 0`
    - right-member only overlap: `0 + 30`
    - split overlap: `20 + 20`
- `validation_policy.adjacency_overlap_mismatch`
  - `error`
  - `warn`
- `validation_policy.uniqueness_checks.*`
  - `off`
  - `warn`
  - `error`
- `validation_policy.reference_contexts[].severity`
  - `warn`
  - `error`

Input-vs-derived boundary in the draft model:

- Intended user/planner inputs:
  - `destination`
  - `product`
  - `fragments`
  - `assembly_order`
  - `junctions[].required_overlap_bp`
  - `junctions[].overlap_partition`
  - `junctions[].distinct_from`
  - `validation_policy`
- Intended normalized/derived outputs:
  - `derived_design.junction_overlaps`
  - `derived_design.primer_design_suggestions`
  - `derived_design.notes`
- Transition fields that may begin as user hints and later become resolved
  values:
  - `junctions[].overlap_source`
  - fragment end strategies (`native_overlap` vs `primer_added_overlap`)

Interpretation:

- Gibson plans are modeled as explicit assembly junctions around an opened
  destination, not merely as an unordered bag of fragments.
- The destination opening defines two terminal junctions, and therefore two
  required overlap regions.
- Inserts may already satisfy those terminal overlaps or may require primer-tail
  adaptation.
- The overlap at one junction should be treated as a selection around the
  in-silico junction, not merely as one scalar length:
  - it may come entirely from the left member
  - entirely from the right member
  - or be split across both members
- For plans that preserve an existing source order, fragment ordering should
  follow ascending bp coordinates (low bp positions first) unless an explicit
  alternative order is requested.
- For multi-fragment plans, internal fragment-fragment junctions may be created
  through primer-added bridge overlaps rather than relying on pre-existing
  native overlap.
- Uniqueness is best treated in layers:
  - destination opening uniqueness: hard validation
  - left/right terminal overlap distinctness: hard validation
  - destination/fragment/genome uniqueness heuristics: advisory checks

Design intent:

- make the same JSON artifact useful for:
  - preflight Gibson validation
  - primer design derivation
  - workflow/macro instantiation
  - factual protocol-cartoon generation
  - reproducible AI-facing project context

Practical overlap heuristics (draft defaults):

- single-insert / two-fragment style assemblies often fit comfortably in the
  20-40 bp range
- multi-fragment assemblies should usually move toward longer overlaps
- a practical starting rule for the draft model is:
  - 1-2 assembled fragments: 20-40 bp overlaps
  - 3-5 assembled fragments: 40 bp overlaps
  - 6+ assembled fragments: 50-100 bp overlaps
- internal multi-fragment junctions introduced by primer-added bridge overlaps
  should follow the same fragment-count-aware guidance rather than being treated
  as exempt from overlap design heuristics

Primer design conventions (draft):

- Gibson primer suggestions should be modeled as two-part primers:
  - `overlap_5prime`
    - non-priming 5' overlap segment used for assembly of adjacent fragments
  - `priming_3prime`
    - gene-specific 3' priming segment used for PCR amplification from the
      source template
- Primer design should start from an in silico assembled product/junction view,
  then work backward to fragment-specific PCR primers.
- Overlap choice is best treated as Tm-aware rather than length-only:
  - simple PCR-fragment-into-vector assemblies may use shorter overlaps when
    overlap Tm is already adequate
  - more complex or multi-fragment assemblies often justify longer overlaps
- The overlap region may lie entirely within one adjacent member or be split
  across the two members around a junction.
- Two primers that implement the same overlap sequence can still belong to
  different PCR reactions, because each primes a different template fragment.

Assembly setup heuristics (draft advisory layer):

- linearized destination can be prepared by PCR amplification or by restriction
  digestion
- PCR cleanup is not always required, but carryover should stay modest
  relative to the final assembly reaction volume
- column purification is especially worth recommending for:
  - assemblies of three or more PCR fragments
  - assemblies involving fragments longer than ~5 kb
- direct vector + insert assemblies often benefit from insert concentration
  above vector concentration
- multi-fragment vector assemblies should generally move toward equimolar
  fragment usage
- some constructs may validate in silico but still perform poorly because of
  biological burden or instability in the propagation host
  (for example repeats or toxic products)

Planning implication:

- these factors should usually surface as `derived_design` advisories or future
  execution/setup guidance, not as hard failures in the core Gibson junction
  model

### `gentle.gibson_assembly_preview.v1`

Purpose:

- provide one deterministic, non-mutating preview response for the current
  Gibson specialist flow,
- keep GUI, shared shell, and direct CLI on the same overlap/primer/cartoon
  derivation path.

Current shared entry points:

- `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
- `gibson apply PLAN_JSON_OR_@FILE`
- GUI specialist window: `Patterns -> Gibson...`

Top-level structure:

- `schema`, `plan_id`, `title`, `summary`
- `can_execute`
- `destination`
  - resolved opening mode/span or cutpoint and actual topology
- `fragments[]`
  - resolved ordered insert rows (fragment id, template seq id, orientation,
    length)
- `insert`
  - compatibility mirror of the first insert row for older single-insert
    consumers
- `resolved_junctions[]`
  - overlap bp, left/right member contributions, overlap Tm, resolved overlap
    sequence, source note
- `primer_suggestions[]`
  - full primer sequence plus explicit `overlap_5prime` and
    `priming_3prime` segments
- `warnings[]`, `errors[]`, `notes[]`
  - includes the shared Tₘ-model note used by GUI/CLI so the assumptions stay
    visible to the user
  - notes also carry explicit design-review guidance that separates:
    - overlap-side success/failure
    - PCR 3' priming-side success/failure
    so adapters can explain when the current blocker is priming rather than
    Gibson overlap derivation
- `suggested_design_adjustments[]`
  - optional structured next-step relaxations when overlap derivation already
    succeeds and the remaining blocker is only the 3' priming window
  - current v1 targets:
    - increasing `priming_segment_max_length_bp`
    - lowering `priming_segment_tm_min_celsius`
  - intended for adapters to offer deterministic “apply and rerun preview”
    actions without parsing prose notes
- `unique_restriction_site`
  - optional structured outcome for a requested
    `validation_policy.desired_unique_restriction_site_enzyme_name`
  - reports whether the requested site was:
    - already unique in the assembled product
    - newly engineered on one terminal overlap
  - carries the enzyme name, terminal side/junction, engineered overlap
    sequence, motif offset, mutation count, and user-facing message so adapters
    do not have to infer this from notes/error prose
- `cartoon`
  - built-in protocol id plus template bindings for single-insert flows
  - multi-insert previews may instead carry one fully resolved
    `ProtocolCartoonSpec` directly
  - intended to stay mechanism-first:
    - show resolved fragment flow and achieved homology/overlap relationships
    - preserve strand-specific 5' chew-back / exposed-tail geometry rather
      than flattening the mechanism to duplex-only blocks
    - avoid drawing full primer objects or low-level PCR parameterization inside
      the cartoon itself
    - keep primer sequences, priming segments, Tm assumptions, and related PCR
      details in adjacent textual/review payloads instead
- `routine_handoff`
  - best-effort Routine Assistant handoff metadata for existing execution paths

Current v1 scope and limits:

- one or more insert fragments in an explicit ordered chain
- destination-first order:
  `destination_left -> insert_1 -> ... -> insert_n -> destination_right`
- the shared preview derives `n + 1` explicit Gibson junctions for `n` inserts
- terminal overlaps are derived from destination context; internal junctions
  are normalized from the adjacent fragment ends / partition rules
- user influence over PCR design stays high-level and Gibson-specific:
  overlap bp range, minimum overlap Tm, priming-segment Tm window, and
  priming-segment length window
- current execution limitation:
  - multi-insert apply currently requires `destination.opening.mode=defined_site`
  - `existing_termini` remains the single-fragment path used by the current
    Routine Assistant handoff
- current unique-site engineering limitation:
  - only the single-insert `defined_site` path is supported
  - only palindromic cutter recognition sequences are currently handled
  - overlap windows must be non-wrapping in the displayed destination sequence
- current Tₘ fields use the shared GENtle Thermo Fisher-style modified
  Allawi/SantaLucia nearest-neighbor estimate with fixed assumptions:
  - exact complement
  - Allawi/SantaLucia 1997 Watson-Crick nearest-neighbor table and terminal
    initiation terms
  - 215 mM effective monovalent salt term, matching the Thermo Fisher
    high-fidelity calculator path
  - 500 nM primer concentration
  - Thermo Fisher high-fidelity empirical adjustment
  - no mismatch/dangling-end/Mg correction
  - fallback to the simple 2/4 estimate for ambiguous or very short sequences
- generic PCR/qPCR request editing is intentionally out of scope for this
  specialist flow
- mutating execution now exists as engine operation
  `ApplyGibsonAssemblyPlan`:
  - consumes the same plan JSON
  - creates deterministic sequence outputs for:
    - left insert primer
    - right insert primer
    - assembled product
  - creates one shared serial arrangement for downstream gel review:
    - original destination vector
    - ordered insert lane(s)
    - assembled product
    - recommended DNA ladders carried with the arrangement for flanking export
  - transfers destination and insert features onto the assembled product
    deterministically through the shared engine path
  - destination features intersecting the consumed opening are now projected
    when a truthful rewrite is available:
    - one-sided overlaps are trimmed to the surviving product span
    - simple spanning features can survive as multipart remnants
    - MCS-like annotations are projected to the edited locus and revalidated
      against actual restriction-enzyme sites on the assembled product
  - the MCS cross-check is product-aware:
    - `mcs_expected_sites` is rewritten to the currently unique cutter set for
      that annotated region on the assembled product
    - `mcs_expected_sites_original`, `mcs_region_sites`,
      `mcs_nonunique_sites`, `mcs_gained_unique_sites`, and
      `mcs_lost_or_nonunique_sites` preserve the cross-check result
    - insert-derived sequence may introduce new sites, and those new sites are
      considered during the same validation pass
  - records one operation-log row so GUI lineage/CLI state replay can reopen
    the specialist from the saved plan without silently re-running it

Normalization/derivation phases:

1. Resolve the destination opening into explicit `dest_left` / `dest_right`
   terminal context.
   - for cutter-derived openings, the resolved coordinates represent the actual
     cleavage window between the recessed termini rather than the whole
     recognition span
   - equal start/end is therefore valid and means a blunt cutpoint
2. Normalize `assembly_order[]` into one adjacency chain.
   - when fragments carry compatible `source_span_1based` hints for one source
     context, default normalization should preserve ascending bp order
3. Materialize one `junction` per adjacent pair in that chain.
4. Derive required overlap sequences from destination flanks and/or adjacent
   fragment termini.
   - respect the junction-specific overlap partition when choosing the final
     overlap sequence around that adjacency
   - internal multi-fragment junctions may instead use designed bridge
     sequences introduced by primer-added overlaps
5. Detect whether each fragment end already satisfies its required overlap or
   requires adaptation (for example primer-added tails).
6. Run hard validation and advisory design checks.
7. Expose derived overlaps, primer design suggestions, and cartoon-ready event
   semantics through `derived_design`.
8. Attach reaction/setup advisories (cleanup, stoichiometry, host-risk notes)
   without conflating them with the hard overlap/junction logic.

Current invariants for the draft model:

- `assembly_order[]` defines the intended adjacency order explicitly.
- `junctions[]` should cover every adjacent pair in `assembly_order[]`.
- terminal junctions are the ones adjacent to the opened destination ends.
- `junctions[].overlap_partition.left_member_bp +
   junctions[].overlap_partition.right_member_bp` should equal
  `junctions[].required_overlap_bp`.
- terminal junction distinctness is a hard validation rule for opened
  destination-vector Gibson plans.
- when source-order hints are present and no contrary manual order is given,
  low bp positions should precede high bp positions in normalization.
- destination-opening uniqueness is a hard validation rule.
- broader destination/fragment/genome uniqueness checks are advisory unless a
  stricter policy is requested.
- `derived_design` may contain unresolved/null sequences at pure planning time;
  this allows the same schema to exist before sequence extraction or primer
  design has been run.
- `junctions[].distinct_from` is currently intended primarily for terminal
  destination-defined junctions, not as a requirement that every internal
  fragment-fragment junction be globally unique.
- `native_overlap` is an expectation about the fragment terminus; it still
  requires sequence confirmation in validation/derivation.
- `designed_bridge_sequence` should be treated as a designed internal overlap,
  suitable for primer-stitching style workflows; GENtle should validate it for
  distinctness/design heuristics rather than treating it as biologically
  privileged just because it was user-supplied.

### `gentle.prepared_cache_inspection.v1`

Purpose:

- provide one deterministic, non-mutating inspection payload for prepared
  reference/helper cache roots,
- let GUI, shared shell, and direct CLI report exactly what GENtle created
  locally before any deletion happens.

Current shared entry points:

- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
- GUI specialist window: `Genome -> Clear Caches...`
- `Prepared References... -> Clear Caches...`

Top-level structure:

- `schema`, `cache_roots[]`
- `entries[]`
- `entry_count`, `total_size_bytes`, `total_file_count`

Entry structure:

- `entry_id`
- `classification`
  - `prepared_install`
  - `orphaned_remnant`
- `cache_root`, `path`
- `artifact_stats[]`
  - `group`
    - `cached_sources`
    - `derived_indexes`
    - `blast_db`
  - `total_size_bytes`
  - `file_count`
- `total_size_bytes`, `file_count`

Inspection rules:

- inspection stays rooted in the selected cache roots only
- default roots are adapter-facing conventions:
  - references: `data/genomes`
  - helpers: `data/helper_genomes`
- orphaned remnants are inspectable even when they are not backed by a
  manifest
- catalog JSON, project state files, MCP/runtime files, backdrop/runtime
  caches, and developer build artifacts are out of scope

### `gentle.prepared_cache_cleanup.v1`

Purpose:

- provide one deterministic cleanup result payload for conservative prepared
  cache deletion workflows,
- keep partial rebuild/reindex cleanup and full prepared-install deletion on
  the same shared contract across GUI/CLI/shell.

Current shared entry points:

- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache ...`
- GUI specialist window: `Genome -> Clear Caches...`

Top-level structure:

- `schema`, `mode`, `cache_roots[]`
- `selected_prepared_ids[]`
- `selected_prepared_paths[]`
- `include_orphaned_remnants`
- `results[]`
- `entry_count`, `removed_item_count`, `removed_bytes`, `removed_file_count`

Per-entry result structure:

- `entry_id`
- `classification`
  - `prepared_install`
  - `orphaned_remnant`
- `cache_root`, `path`
- `removed`
- `removed_artifact_groups[]`
- `removed_bytes`, `removed_file_count`
- `skipped_reason?`

Cleanup modes:

- `blast_db_only`
  - remove only BLAST DB sidecars for selected manifest-backed installs
- `derived_indexes_only`
  - remove BLAST DB sidecars plus `sequence.fa.fai` and `genes.json`
  - cached sources and manifests remain so reindex-from-cached-files still
    works
- `selected_prepared_installs`
  - remove only explicitly selected prepared installs
  - optional `include_orphaned_remnants` also removes orphaned remnants under
    the same selected roots
- `all_prepared_in_cache`
  - remove all prepared installs under the selected roots
  - optional `include_orphaned_remnants` extends that deletion to orphaned
    remnants

Cleanup rules:

- `blast_db_only` and `derived_indexes_only` apply only to manifest-backed
  prepared installs
- selective cleanup modes accept either `selected_prepared_ids[]` or
  `selected_prepared_paths[]`
- `selected_prepared_paths[]` are the precise selector when duplicate prepared
  ids exist across multiple selected cache roots
- orphaned remnants can only be deleted through the full-delete modes
- cleanup never scans the whole workspace; it only touches the selected roots
- cleanup does not treat catalog JSON, `.gentle_state.json`, MCP/runtime files,
  backdrop/runtime caches, or `target/` as cache

## Core entities

### ProjectState

```json
{
  "sequences": {"seq_id": "DNAsequence object"},
  "metadata": {"any": "json"},
  "display": {"ui_visibility_tfbs_and_linear_viewport_state": "..."},
  "lineage": {"nodes": {}, "edges": []},
  "parameters": {"max_fragments_per_container": 80000},
  "container_state": {
    "containers": {},
    "arrangements": {},
    "racks": {},
    "seq_to_latest_container": {}
  }
}
```

Semantic interpretation:

- In GUI terms, a project window represents a wet-lab container context.
- A container may map to multiple candidate sequences/fragments.
- Explicit container objects are first-class state (`container_state`) and are
  indexed from sequence ids via `seq_to_latest_container`.
- Containers now also record `declared_contents_exclusive`:
  - `true` (default): the declared members are intended to be the full known
    contents of a clean vial/tube
  - `false`: the declared members are measured/known constituents of a more
    complex sample, and additional unlisted molecules may also be present
- Arrangements stay the semantic experiment-order layer.
- Racks are the linked physical placement layer and may host one or more
  arrangements without changing arrangement identity.

### Rack placement entities

- `RackProfileKind`
  - built-in physical carriers:
    - `small_tube_4x6`
    - `plate_6`
    - `plate_96`
    - `plate_384`
  - persisted custom snapshots use:
    - `custom`
- `RackProfileSnapshot`
  - persisted row/column/fill-direction/blocked-slot snapshot used by one saved rack
  - `fill_direction`
    - `row_major`
    - `column_major`
  - `blocked_coordinates[]`
    - normalized A1-style coordinate list
- `Rack`
  - one saved physical rack/plate draft
- `RackPlacementEntry`
  - one occupied A1-style coordinate on that rack
  - points back to:
    - `arrangement_id`
    - arrangement-local `order_index`
    - one `occupant`
- `RackOccupant`
  - `container`
  - `ladder_reference`

Rack-placement invariants:

- rack placement consumes arrangement order instead of duplicating experiment
  meaning in a second free-form list
- default placement is deterministic:
  - choose the smallest fitting built-in profile
  - fill row-major
  - use A1-style coordinates
- saved rack snapshots may then refine physical layout with:
  - `fill_direction = row_major|column_major`
  - `blocked_coordinates[]`
- A1-style row labels continue beyond `Z` as `AA`, `AB`, ...
- moving one sample or arrangement block is shift-neighbor by default; it
  preserves occupied order instead of creating arbitrary holes

### `gentle.rack_state.v1`

Purpose:

- provide one deterministic inspection payload for saved rack state
- keep GUI rack view and CLI/shell inspection on one shared state contract

Current shared entry point:

- `racks show RACK_ID`

Top-level structure:

- `schema`
- `rack`
- `placements[]`

Placement payload:

- `coordinate`
- `arrangement_id`
- `order_index`
- `role_label`
- `occupant`
  - `kind = container`
    - `container_id`
    - `container_name?`
    - `seq_id?`
  - `kind = ladder_reference`
    - `ladder_name`
  - `kind = empty`

### Operation

Current draft operations:

- `LoadFile { path, as_id? }`
- `CreateSequenceFromText { sequence_text, output_id?, name?, circular=false }`
  - creates a persistent synthetic project sequence from inline sequence text
  - whitespace is ignored and bases are stored upper-case
  - records ordinary lineage/container state so GUI, GUI Shell, Agent
    Assistant, CLI `op`, workflow, and MCP `op` callers can use the same
    created sequence id afterward
- `SaveFile { seq_id, path, format }`
- `RenderSequenceSvg { seq_id, mode, path }`
  - linear exports honor the current stored linear viewport in `display`
    (`linear_view_start_bp` / `linear_view_span_bp`) when that viewport is a
    proper subsequence crop
  - single-base `variation` features render as baseline markers in linear SVG
    output rather than as generic detached feature blocks
  - linear exports now also mark transcription starts/directions for
    strand-bearing `gene`/`mRNA`/`CDS`/`promoter` features and suppress
    unlabeled fallback coordinate text that would otherwise clutter
    figure-oriented exports
  - linear exports also prefer gene-style labels over accession-only transcript
    ids when possible and compact nearby repeated non-gene labels
  - direction-bearing `mRNA`/`promoter` bars render with arrowed ends, and the
    linear TSS cue uses a short hooked arrow so direction survives
    figure-oriented contexts
  - circular exports now use a transparent canvas and render single-base
    `variation` features as explicit radial markers on the DNA ring
  - circular exports also mark transcription starts for strand-bearing
    `gene`/`mRNA`/`CDS`/`promoter` features with a short arrow shaft plus
    direction arrowhead
  - circular exports also use a slightly larger ring and larger label fonts so
    figure-oriented construct maps stay readable when embedded in docs
  - `repeat_region` / `mobile_element` features carrying RepeatMasker/UCSC
    `rmsk`-style qualifiers (`repName`, `repClass`, `repFamily`, or
    `rmsk_*` / `repeat_*` / `rpt_*` aliases) use shared repeat subtype labels
    and colors across GUI maps, feature-tree grouping/filtering, and SVG
    export
  - SVG feature blocks, feature labels, restriction labels, and selected
    directional/variation markers carry stable `data-gentle-role` attributes
    plus `data-gentle-feature-kind` when applicable, so visual-regression tests
    and downstream tooling can inspect semantic glyphs without relying on
    raster screenshots
- `RestrictionSiteExpertView` is also the shared restriction-site detail
  contract for map hover/popover UI:
  - it exposes enzyme grouping, selected enzyme, site count, recognition span,
    explicit top/bottom cut positions, recognition/site sequence, complement,
    cut geometry, optional enzyme note, and REBASE URL
  - `tooltip_lines[]` carries concise presentation-ready lines with
    caret-marked top and bottom cut positions so GUI hover text, shared shell
    JSON, MCP tools, and future inspectors do not each invent separate
    cut-geometry presentation rules
- `RenderDotplotSvg { seq_id, dotplot_id, path, flex_track_id?, display_density_threshold?, display_intensity_gain?, overlay_x_axis_mode? }`
- `RenderFeatureExpertSvg { seq_id, target, path }`
  - shared renderer contract across GUI/CLI/JS/Lua for TFBS/restriction/splicing/isoform expert exports
  - splicing SVG includes explicit junction-support counts, frequency-encoded transcript-vs-exon matrix coloring, predicted exon->exon transition matrix support coloring, exon `len%3` (genomic-length modulo 3) cues, and CDS flank phase edge coloring (`0/1/2`) when transcript `cds_ranges_1based` are available
- `SummarizeTfbsScoreTracks { target, motifs, score_kind, clip_negative, path? }`
  - non-mutating continuous motif-score export for Promoter design and headless
    ClawBio/OpenClaw-style inspection
  - accepts the shared `SequenceScanTarget` operand so the same summary path
    works for stored `seq_id` spans and inline ASCII DNA
  - returns schema `gentle.tfbs_score_tracks.v1`
  - reports per-position forward and reverse score tracks for each requested
    TF/PSSM token across one shared span
  - `score_kind` accepts `llr_bits`, `llr_quantile`,
    `llr_background_quantile`, `llr_background_tail_log10`,
    `true_log_odds_bits`, `true_log_odds_quantile`,
    `true_log_odds_background_quantile`, or
    `true_log_odds_background_tail_log10`
  - `clip_negative=true` clamps negative scores to `0.0` for the bit-based
    kinds, which is the intended display mode for promoter-design plots when
    only positive support should be shown
  - the background-normalized kinds convert raw motif scores into modeled
    percentile or `-log10(tail)` views against a quantized IID random-DNA
    window model and zero out everything below the `0.95`
    random-background quantile
  - each returned track also carries a deterministic background-normalization
    block that now combines a larger deterministic random sample with
    theoretical motif-score bounds and modeled tail calibration, so the score
    tracks no longer flatten onto the old finite-sample `1.0` ceiling
  - the same report now also carries one pairwise correlation sidecar:
    - exact raw Pearson correlation over the displayed per-position signal
    - smoothed Pearson correlation over the same signal after centered boxcar
      smoothing (`25 bp`)
    - exact raw Spearman correlation over the same displayed signal
    - smoothed Spearman correlation over the same centered-boxcar signal
    - and a signed primary-peak offset so “are the strongest windows actually
      co-localized?” is explicit instead of guessed from the traces alone
  - the report’s TSS markers now also fall back to promoter-slice provenance
    when the selected sequence was derived through
    `ExtractGenomePromoterSlice`/`extract-promoter` with
    `annotation_scope=none`, so exported promoter plots can still carry one
    explicit transcription-start marker even without imported feature rows
  - `path` writes the same structured JSON report to disk for reuse outside the
    current adapter session
- `RenderTfbsScoreTracksSvg { target, motifs, score_kind, clip_negative, path }`
  - shared stacked SVG renderer for the same TFBS score-track payload
  - accepts the same `SequenceScanTarget` operand as
    `SummarizeTfbsScoreTracks`
  - writes a deterministic figure suitable for GUI/CLI/agent/README parity
  - reuses the `gentle.tfbs_score_tracks.v1` summary contract internally and
    also returns that report in `OpResult.tfbs_score_tracks`
  - the figure includes the same score-family-aware normalization context used
    by the JSON payload, surfaced either as `p99 / Δp99 / bg+` or as
    `theory max / peak q / -log10 tail`
  - stacked track labels now render as `TF name (JASPAR id)` when both are
    known, which avoids ambiguity for motifs whose names are also valid IUPAC
    strings (for example `MYC`)
  - TSS markers now render as one vertical dashed line with one top-level
    kinked arrow/label per position so the whole stacked figure reads as one
    shared DNA span instead of unrelated per-row ornaments
- `ScanTfbsHits { target, motifs, min_llr_bits?, min_llr_quantile?, per_tf_thresholds?, max_hits?, path? }`
  - non-mutating thresholded TFBS/JASPAR hit scan over either stored `seq_id`
    input or inline ASCII DNA through `SequenceScanTarget`
  - returns schema `gentle.tfbs_hit_scan.v1`
  - intended as the shared low-latency inspection path for CLI/shell/agents and
    future GUI selection-first actions before a caller explicitly chooses to
    materialize `TFBS` features
  - `path` writes the same structured JSON report to disk for reuse outside the
    current adapter session
- `RenderIsoformArchitectureSvg { seq_id, panel_id, expression_tsv_path?, path }`
- `RenderRnaStructureSvg { seq_id, path }`
- `RenderLineageSvg { path }`
- `RenderPoolGelSvg { inputs, path, ladders?, container_ids?, arrangement_id?, conditions?, render_options? }`
- `RenderProteinGelSvg { report_id, path, ladders? }`
- `RenderProteinGelReportsSvg { report_ids[], path, ladders? }`
- `RenderProteaseDigestGelSvg { seq_id?, report_id?, transcript_id?, proteases[], path, min_length_aa?, ladders? }`
- `RenderProtein2dGelSvg { report_id, path, ladders? }`
- `CreateArrangementSerial { container_ids, arrangement_id?, name?, ladders? }`
- `SetArrangementLadders { arrangement_id, ladders? }`
- `SetContainerDeclaredContentsExclusive { container_id, exclusive }`
- `CreateRackFromArrangement { arrangement_id, rack_id?, name?, profile? }`
- `PlaceArrangementOnRack { arrangement_id, rack_id }`
- `MoveRackPlacement { rack_id, from_coordinate, to_coordinate, move_block? }`
- `MoveRackSamples { rack_id, from_coordinates[], to_coordinate }`
- `MoveRackArrangementBlocks { rack_id, arrangement_ids[], to_coordinate }`
- `SetRackProfile { rack_id, profile }`
- `ApplyRackTemplate { rack_id, template }`
- `SetRackFillDirection { rack_id, fill_direction }`
- `SetRackProfileCustom { rack_id, rows, columns }`
- `SetRackBlockedCoordinates { rack_id, blocked_coordinates }`
- `ExportRackLabelsSvg { rack_id, path, arrangement_id?, preset }`
- `ExportRackFabricationSvg { rack_id, path, template }`
- `ExportRackIsometricSvg { rack_id, path, template }`
- `ExportRackHeroSvg { rack_id, path, template }`
- `ExportRackOpenScad { rack_id, path, template }`
- `ExportRackCarrierLabelsSvg { rack_id, path, arrangement_id?, template, preset }`
- `ExportRackSimulationJson { rack_id, path, template }`
- `RenderProtocolCartoonSvg { protocol, path }`
- `RenderProtocolCartoonTemplateSvg { template_path, path }`
- `ValidateProtocolCartoonTemplate { template_path }`
- `RenderProtocolCartoonTemplateWithBindingsSvg { template_path, bindings_path, path }`
- `ExportProtocolCartoonTemplateJson { protocol, path }`
- `ExportDnaLadders { path, name_filter? }`
- `ExportRnaLadders { path, name_filter? }`
- `ExportPool { inputs, path, pool_id?, human_id? }`
- `ExportPoolCollection { collection_subject, path, pool_id?, human_id? }`
- `ExportProcessRunBundle { path, run_id? }`
- `ExportLabAssistantInstructions { path, run_id?, title?, audience?, format? }`
- `Digest { input, enzymes, output_prefix? }`
- `Ligation { inputs, circularize_if_possible, protocol, output_id?, output_prefix?, unique? }`
- `MergeContainers { inputs, output_prefix? }`
- `Pcr { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrAdvanced { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrMutagenesis { template, forward_primer, reverse_primer, mutations, output_id?, unique?, require_all_mutations? }`
- `DesignPrimerPairs { ... }` (implemented baseline)
- `ExportPrimerDesignReport { report_id, path }`
- `DesignTerminalExonRtPrimerPool { request }` (implemented; ordered,
  transcript-aware fixed-adapter RT-primer-pool design)
- `ExportTerminalExonRtPrimerPoolReport { report_id, path }`
- `PcrOverlapExtensionMutagenesis { ... }` (implemented baseline; insertion/deletion/replacement overlap-extension flow)
- `DesignQpcrAssays { ... }` (implemented baseline; forward/reverse/probe)
- `TestCdnaPcr { seq_id, source_feature_id, forward_primer, reverse_primer, transcript_id?, transcript_order?, transcript_map_coordinate_mode?, min_amplicon_bp?, max_amplicon_bp?, max_mismatches?, require_3prime_exact_bases?, path?, svg_path?, materialize_products?, product_output_prefix?, product_gel_svg_path?, product_gel_ladders? }` (implemented baseline; transcript-derived cDNA assay test with opt-in product materialization/gel)
- `TestCdnaQpcr { seq_id, source_feature_id, forward_primer, reverse_primer, probe, transcript_id?, transcript_order?, transcript_map_coordinate_mode?, min_amplicon_bp?, max_amplicon_bp?, max_mismatches?, require_3prime_exact_bases?, path?, svg_path?, materialize_products?, product_output_prefix?, product_gel_svg_path?, product_gel_ladders? }` (implemented baseline; transcript-derived cDNA assay test with internal probe and opt-in product materialization/gel)
- `BuildTranscriptQpcrPanel { seq_id, source_feature_id, shared_qpcr_report_id, path? }` (implemented baseline; shared qPCR components plus transcript-characteristic forward-primer table)
- `TestCdnaQpcrFasta { cdna_fasta_paths[], forward_primer, reverse_primer, probe, transcript_id?, min_amplicon_bp?, max_amplicon_bp?, max_mismatches?, require_3prime_exact_bases?, path?, svg_path? }` (implemented baseline; cDNA/ncRNA FASTA/FASTA.gz screen with internal probe)
- `ComputeDotplot { seq_id, reference_seq_id?, span_start_0based?, span_end_0based?, reference_span_start_0based?, reference_span_end_0based?, mode, word_size, step_bp, max_mismatches?, tile_bp?, store_as? }` (implemented baseline, self + pairwise)
- `ComputeFlexibilityTrack { seq_id, span_start_0based?, span_end_0based?, model, bin_bp, smoothing_bp?, store_as? }` (implemented baseline)
- `DeriveSplicingReferences { seq_id, span_start_0based, span_end_0based, seed_feature_id?, scope?, output_prefix? }` (implemented baseline; emits derived DNA window + mRNA isoforms + exon-reference sequence)
- `AlignSequences { query?, target?, query_seq_id?, target_seq_id?, query_span_start_0based?, query_span_end_0based?, target_span_start_0based?, target_span_end_0based?, mode?, match_score?, mismatch_score?, gap_open?, gap_extend? }` (implemented baseline; `query`/`target` use `SequenceScanTarget` and can be stored `seq_id` or inline ASCII; legacy `*_seq_id` + span fields remain accepted; returns structured pairwise local/global report in `OpResult.sequence_alignment`)
- `ImportSequencingTrace { path, trace_id?, seq_id? }` (implemented baseline; imports one ABI/AB1 or SCF evidence file into the shared sequencing-trace store without mutating construct sequences)
- `ListSequencingTraces { seq_id? }`
- `ShowSequencingTrace { trace_id }`
- `ConfirmConstructReads { expected_seq_id, baseline_seq_id?, read_seq_ids?, trace_ids?, targets?, alignment_mode?, match_score?, mismatch_score?, gap_open?, gap_extend?, min_identity_fraction?, min_target_coverage_fraction?, allow_reverse_complement?, report_id? }` (implemented baseline; accepts already-loaded read sequences and/or imported sequencing traces as evidence inputs into one shared confirmation report, with optional baseline context for intended-edit vs reversion classification)
- `ReadAcquireStatus { manifest_path, cache_dir, work_dir }`
- `ReadAcquirePrepare { manifest_path, cache_dir, work_dir, analysis_format, read_layout, threads?, max_size?, min_free_gb?, drop_intermediate_fastq, continue_on_error }`
- `ReadAcquireInspect { sra_accession, cache_dir, work_dir }`
- `ReadAcquireCancel { sra_accession, cache_dir, work_dir }`
  - shared SRA-backed read acquisition uses external `prefetch`,
    `vdb-validate`, and `fasterq-dump`, but reports through
    `gentle.read_acquisition_report.v1`.
  - manifest TSV rows require `sample_id` and `sra_accession`; optional fields
    are `sample_name`, `assay_kind`, `read_layout`, `analysis_format`, and
    `note`.
  - per-run lifecycle uses
    `resource_key = "read_acquisition:<SRA_ACCESSION>"` and
    `lifecycle_status = missing|running|ready|failed|cancelled|stale`.
  - running activity JSON includes a cooperative `cancel_path`, phase/item,
    produced-byte estimate, and `monitored_free_bytes` /
    `minimum_free_bytes` when `min_free_gb` is configured.
  - activity JSON phases are `prefetch`, `validate_sra`, `dump_fastq`,
    `convert_fasta`, `verify_output`, and `complete`.
  - `ReadAcquireCancel` writes the run cancel marker; the active prepare loop
    terminates the external child process group and records
    `lifecycle_status = cancelled`.
  - `ReadAcquirePrepare` rechecks available disk while external SRA Toolkit
    phases run and fails early if free space drops below `min_free_gb`.
  - final `.sra` and prepared FASTA/FASTQ outputs are never deleted
    automatically; `drop_intermediate_fastq` only applies to FASTQ files that
    were converted into FASTA outputs.
- `InterpretRnaReads { seq_id, seed_feature_id, profile, input_path, input_format, scope, origin_mode?, target_gene_ids?, roi_seed_capture_enabled?, seed_filter, align_config, report_id?, report_mode?, checkpoint_path?, checkpoint_every_reads?, resume_from_checkpoint? }` (Nanopore cDNA phase-1 seed-filter pass; `multi_gene_sparse` expands local transcript-template indexing, while ROI capture remains planned)
- `AlignRnaReadReport { report_id, selection, align_config_override?, selected_record_indices? }` (Nanopore cDNA phase-2 retained-hit alignment pass; updates mapping/MSA/abundance report fields and re-ranks retained hits by alignment-aware retention rank)
- `PreflightRnaReadIsoforms { seq_id, seed_feature_id, scope?, seed_filter?, optimize_parameters?, positive_transcript_fasta_paths?, control_transcript_fasta_paths?, max_control_match_probability? }` (non-mutating target-transcript seed representation preflight; repeated positive FASTAs are hard must-pass transcript variants, repeated negative control FASTAs are grouped by inferred gene/symbol, and optimizer candidates are rejected when any positive transcript fails or any control group exceeds the configured match-probability ceiling)
- `ListRnaReadReports { seq_id? }`
- `ShowRnaReadReport { report_id }`
- `ExportRnaReadReport { report_id, path }`
- `ExportRnaReadHitsFasta { report_id, path, selection, selected_record_indices?, subset_spec? }`
- `ExportRnaReadSampleSheet { path, seq_id?, report_ids?, gene_ids?, complete_rule?, append? }`
- `ExportRnaReadTargetQuality { report_id, path, gene_ids, complete_rule? }`
- `ExportRnaReadExonPathsTsv { report_id, path, selection, selected_record_indices?, subset_spec? }`
- `ExportRnaReadExonAbundanceTsv { report_id, path, selection, selected_record_indices?, subset_spec? }`
- `ExportRnaReadScoreDensitySvg { report_id, path, scale, variant }`
- `ExportRnaReadAlignmentsTsv { report_id, path, selection, limit?, selected_record_indices?, subset_spec? }`
- `ExportRnaReadIsoformTriageTsv { report_id, path, selection, limit?, selected_record_indices?, subset_spec?, min_identity_fraction?, min_query_coverage_fraction?, min_confirmed_transition_fraction?, max_secondary_mappings? }`
- `ExportRnaReadAlignmentDotplotSvg { report_id, path, selection, max_points }`
- `ExtractRegion { input, from, to, output_id? }`
- `PrepareGenome { genome_id, catalog_path?, cache_dir?, timeout_seconds? }`
  - prepared-resource inspection invokes `blastdbcmd` and exposes
    `gentle.blast_database_inspection.v1` for the BLAST component independently
    of the overall preparation activity marker
  - the report records `index_kind = genomic_dna|transcriptome_cdna`, source
    genome/assembly/release, masking declaration, database prefix and version,
    sequence count, total bases, `blastdbcmd` executable/version, deterministic
    index-content fingerprint and algorithm, validation status, warnings, and
    compatible GENtle operations
  - `genomes status` / `helpers status` return the inspection under
    `components.blast_database`; a valid component remains reusable when an
    older whole-prepare activity record is `stale`
  - new reference-genome manifests default to `genomic_dna`; transcriptome/cDNA
    databases use the distinct `transcriptome_cdna` identity so their
    specificity interpretation cannot be conflated with genomic carryover
- Sequence-only BLAST resources use the shared-shell routes
  `genomes prepare-blast-resource` (alias `import-blast-resource`) and
  `genomes inspect-blast-resource` rather than `PrepareGenome`; an existing
  valid database whose catalog metadata needs correction uses the separate,
  non-destructive `genomes adopt-blast-resource` route:
  - the catalog entry declares `blast_index_kind = transcriptome_cdna`,
    `sequence_local` or `sequence_remote`, and optional `reference_name`,
    `reference_release`, and `blast_masking`
  - preparation returns `gentle.prepare_blast_sequence_resource.v1` and writes
    an engine-owned `gentle.blast_sequence_resource_manifest.v1`; callers do
    not author the prepared manifest
  - the manifest binds the source/digest, materialized FASTA, BLAST prefix and
    executable, index-content fingerprint plus its full-or-edge-sampled
    algorithm, index kind, masking, reference identity/release, and install
    time. It deliberately has no genome gene/transcript indexes
  - adoption requires an expected database fingerprint, validates the existing
    index with `blastdbcmd`, records the prepared FASTA checksum (and marks it
    verified only when an optional expected SHA-1 was supplied and matched),
    regenerates only its subject-annotation sidecar, writes
    `gentle.blast_resource_adoption.v1`, and never invokes `makeblastdb` or
    deletes index files. By default the catalog-source link is explicitly
    `attested`; `--verify-source-identities` upgrades it to `verified` only
    after every FASTA accession/length matches `blastdbcmd -entry all`
  - once a manifest records a database fingerprint, neither a changed sequence
    source nor changed index content at the same prefix can be accepted through
    ordinary preparation. Reviewed external replacement uses adoption;
    GENtle-driven replacement requires the explicit `--force-rebuild` flag
  - preparation also writes
    `gentle.blast_subject_annotation_index.v1`, a compact sidecar keyed by
    subject/transcript id. Ensembl FASTA deflines are parsed once for optional
    gene identity, symbol, biotypes, and description; GenBank records retain
    their versioned accession plus `/gene` and `/db_xref` values. Unknown
    identity stays explicitly unavailable. The sidecar is bound to both the
    prepared sequence checksum and its own checksum
  - inspection reuses `gentle.blast_database_inspection.v1`, including
    `blastdbcmd` validation, sequence count, content fingerprint, and compatible
    operations, plus subject-annotation status, record/annotated counts, source
    counts, path, and fingerprint
  - prepared sequence-only resources resolve through the same exhaustive
    `blastn-short` executor used by primer specificity. They do not fall back
    to a similarly named genome install
- `ExtractGenomeRegion { genome_id, chromosome, start_1based, end_1based, output_id?, annotation_scope?, max_annotation_features?, include_genomic_annotation?, catalog_path?, cache_dir? }`
  - `annotation_scope` accepts `none|core|full` and defaults to `core` when omitted.
  - `max_annotation_features` is an optional safety cap (0 or omitted = unlimited for explicit requests).
  - legacy `include_genomic_annotation` is still accepted (`true` -> `core`, `false` -> `none`) for compatibility.
  - operation results include `genome_annotation_projection` telemetry (requested/effective scope, feature counts, fallback metadata).
  - for helper genome IDs containing `pUC18`/`pUC19`, the engine applies a deterministic fallback MCS `misc_feature` annotation when source annotation does not already include an MCS feature and exactly one canonical MCS motif is found.
  - source-derived and fallback MCS features expose `mcs_expected_sites` with REBASE-normalized enzyme names when recognizable.
- `ExtractGenomeGene { genome_id, gene_query, occurrence?, output_id?, extract_mode?, promoter_upstream_bp?, annotation_scope?, max_annotation_features?, include_genomic_annotation?, catalog_path?, cache_dir? }`
  - `annotation_scope` accepts `none|core|full` and defaults to `core` when omitted.
  - `max_annotation_features` is an optional safety cap (0 or omitted = unlimited for explicit requests).
  - legacy `include_genomic_annotation` is still accepted (`true` -> `core`, `false` -> `none`) for compatibility.
  - operation results include `genome_annotation_projection` telemetry (requested/effective scope, feature counts, fallback metadata).
  - for helper genome IDs containing `pUC18`/`pUC19`, the same deterministic MCS fallback annotation behavior applies when an MCS feature is missing; non-unique motif matches are warned and skipped.
- `ExtractGenomePromoterSlice { genome_id, gene_query, occurrence?, transcript_id?, output_id?, upstream_bp?, downstream_bp?, annotation_scope?, max_annotation_features?, include_genomic_annotation?, catalog_path?, cache_dir? }`
  - derives one unclipped promoter slice directly from transcript TSS geometry instead of requiring a separate gene extraction + TSS recovery + region extraction chain.
  - when `transcript_id` is omitted, the engine deterministically chooses the outermost 5' transcript for the matched gene and warns when multiple transcript candidates exist.
  - `upstream_bp` and `downstream_bp` default to the shared promoter-window baseline (`1000` upstream, `200` into the transcribed area).
  - `annotation_scope` accepts `none|core|full` and defaults to `core` when omitted.
  - `max_annotation_features` is an optional safety cap (0 or omitted = unlimited for explicit requests).
  - legacy `include_genomic_annotation` is still accepted (`true` -> `core`, `false` -> `none`) for compatibility.
  - provenance rows record `gene_query`, `occurrence`, `transcript_id`,
    `tss_1based`, and both promoter flank lengths so GUI/CLI/agents can audit
    how the slice was derived and later recover an explicit TSS marker even
    when feature projection was intentionally disabled.
- `ExtendGenomeAnchor { seq_id, side, length_bp, output_id?, catalog_path?, cache_dir?, prepared_genome_id? }`
- `VerifyGenomeAnchor { seq_id, catalog_path?, cache_dir?, prepared_genome_id? }`
- `ProjectMicroarrayTrack { seq_id, manifest_path, contrasts, level, min_abs_logfc?, max_adj_p?, max_features?, clear_existing }`
- `ProjectGenomeInterval { source_genome_id, target_genome_id, projection_path, chrom, start_1based, end_1based, strand? }`
- `ListCutRunDatasets { filter?, catalog_path? }`
- `ShowCutRunDatasetStatus { dataset_id, catalog_path?, cache_dir? }`
- `PrepareCutRunDataset { dataset_id, catalog_path?, cache_dir? }`
- `ProjectCutRunDataset { seq_id, dataset_id, include_peaks?, include_signal?, clear_existing?, catalog_path?, cache_dir? }`
- `InterpretCutRunReads { seq_id, input_r1_path?, input_r2_path?, dataset_id?, catalog_path?, cache_dir?, input_format?, read_layout?, roi_flank_bp?, seed_filter?, align_config?, deduplicate_fragments?, report_id?, checkpoint_path?, checkpoint_every_reads? }`
- `ListCutRunReadReports { seq_id? }`
- `ShowCutRunReadReport { report_id }`
- `ExportCutRunReadCoverage { report_id, path, kind? }`
- `InspectCutRunRegulatorySupport { seq_id, dataset_ids, read_report_ids, catalog_path?, cache_dir?, promoter_search_start_0based?, promoter_search_end_0based_exclusive?, neighbor_window_bp, species_filters, path? }`
  - V1 is processed-evidence-first and currently reuses the shared anchored
    `ImportGenomeBedTrack` / `ImportGenomeBigWigTrack` projection behavior.
  - prepared CUT&RUN datasets now also expose one shared lifecycle contract:
    - `resource_key = "cutrun_dataset:<DATASET_ID>"`
    - `lifecycle_status = missing|running|ready|failed|cancelled|stale`
    - `current_activity` carries the persisted lease/heartbeat marker when a
      dataset prepare is active or when the last attempt ended as
      `failed|cancelled|stale`
  - duplicate `PrepareCutRunDataset` callers now reuse one active dataset
    prepare lease instead of materializing the same dataset twice in parallel.
  - V2 is ROI-first and interprets either ad hoc `FASTA`/`FASTQ` inputs or
    prepared catalog-linked raw reads against one selected genome-anchored
    region plus deterministic flanks.
  - when `roi_flank_bp = 0` and the anchor span exactly matches the imported
    sequence length, V2 maps directly against that anchored sequence and does
    not require the corresponding whole reference genome to be prepared;
    nonzero flanks still require the prepared reference sequence.
  - paired-end interpretation is first-class: mates are paired by normalized
    read id, concordant pairs emit fragment spans, and orphan/single-ended
    observations are retained as explicit report rows instead of being dropped.
  - when `dataset_id` is provided, the engine resolves `reads_r1` / `reads_r2`
    from the prepared CUT&RUN manifest and infers the raw-read format from the
    prepared file names.
  - CUT&RUN catalog entries may declare `reads_sra_accession`; preparation then
    delegates acquisition/conversion to the shared read-acquisition SRA path
    before writing the prepared CUT&RUN manifest. The built-in Rostock p73
    `E-MTAB-15709` shard uses this route, so full read downloads are deliberate
    `cutrun prepare` actions rather than routine status/list operations.
  - `ExportCutRunReadCoverage` writes TSV summaries for one saved read report:
    `coverage`, `cut_sites`, or `fragments`.
  - `InspectCutRunRegulatorySupport` is the first shared V3 reasoning surface:
    - it accepts one anchored `seq_id` plus repeated prepared `dataset_ids`
      and/or saved `read_report_ids`
    - optional `catalog_path` and `cache_dir` are used for every prepared
      dataset lookup and are echoed in the portable report, so an inspection
      can replay the same non-default V1 cache used by prepare/project
    - strong support windows can be derived from saved V2 read reports alone,
      from prepared signal-only evidence, or from prepared peak evidence
    - theoretical TFBS rows keep the legacy `confirmed_tfbs_rows` /
      `unconfirmed_tfbs_rows` split for existing consumers; each row also
      carries additive `support_status` with
      `confirmed|nearby|absent|motif_poor` so UIs can distinguish overlapping
      strong occupancy, nearby occupancy, no nearby occupancy, and weak/no
      motif sequence contexts without changing the old row vectors; rows also
      include additive `support_distance_bp` when nearest support is available
    - `confirmed` means a strong support window overlaps the theoretical TFBS;
      `nearby` means strong support is within the neighbor threshold but does
      not overlap; `absent` means no such support window was found; and
      `motif_poor` means the support window is weak for the requested motif or
      the motif cannot be resolved locally
    - motif-absent strong windows are reported separately and classified as
      `context_supported_by_other_motifs` or `motif_poor_supported`
    - if the assayed target factor cannot be resolved to a local motif, strong
      support windows are still reported through context-only motif reasoning
      instead of being dropped
    - motif-context scans intentionally use a high-confidence
      `motif_context_min_llr_quantile` of `0.95` by default, and resolved
      target motifs are excluded from the "other motifs" context list
    - recurring motif summaries are emitted for motifs found inside supported
      windows and in their neighbor flanks
  - default catalog: `assets/cutrun.json`; default prepared-cache root:
    `data/cutrun`; environment override: `GENTLE_CUTRUN_CACHE_DIR`.
  - `ProjectCutRunDataset` requires a genome-anchored sequence and rejects
    incompatible `supported_reference_genome_ids`.
  - status/projection payload schemas:
    `gentle.cutrun_dataset_list.v1`,
    `gentle.cutrun_prepared_manifest.v1`,
    `gentle.cutrun_dataset_status.v1`,
    `gentle.cutrun_dataset_projection.v1`,
    `gentle.cutrun_read_report.v1`,
    `gentle.cutrun_read_reports.v1`,
    `gentle.cutrun_read_coverage_export.v1`,
    `gentle.cutrun_regulatory_support.v1`,
    `gentle.gene_set_cutrun_regulatory_support.v1`.
  - `cutrun gene-set-regulatory-support` consumes a resolved gene-set promoter
    cohort or enough local inputs to build one, scores only prepared datasets
    and saved read reports, and keeps `evaluated` and `unevaluated` member
    states separate so set-level fractions use evaluated members only. The
    optional `--relationship manual|co-regulated|anti-co-regulated` flag echoes
    the declared cohort expectation and emits non-blocking `relationship_flags`
    and warnings from promoter occupancy over evaluated members only, without
    gating the report.

Microarray track projection notes:

- `ProjectMicroarrayTrack` projects prepared Clariom D/microarray contrast rows
  onto an existing genome-anchored sequence.
- raw operation example:
  `{"ProjectMicroarrayTrack":{"seq_id":"grch38_tp73","manifest_path":"data/publication_resources/rostock_p73_clariomd_e_mtab_14704/analysis/clariomd_probe_level/clariomd_microarray_track_manifest.json","contrasts":["AdTAp73alpha-AdGFP","AdTAp73beta-AdGFP"],"level":"probeset","min_abs_logfc":0.5,"max_adj_p":0.05,"max_features":5000,"clear_existing":true}}`
- manifest schema: `gentle.microarray_track_manifest.v1`
- projection report schema: `gentle.microarray_projection_report.v1`
- explicit interval-map projection report schema:
  `gentle.genome_coordinate_projection_report.v1`
- probe/probeset-region preflight schema: `gentle.probe_region_plan.v1`
- probe/probeset-region helper-output inspection schema:
  `gentle.probe_region_output_inspection.v1`
- probe/probeset-region backend execution schema:
  `gentle.probe_region_backend_run.v1`
- probe/probeset-region evidence interpretation schema:
  `gentle.probe_region_evidence_interpretation.v2`
- the manifest records dataset id, platform, normalization method, contrast
  order, coordinate system, supported `genome_id` aliases, and per-contrast TSV
  paths.
- `arrays probe-regions` accepts explicit CEL paths or a publication-resource
  dataset id plus gene, locus, transcript-cluster, or probeset selectors. The
  stage-one report is read-only: for publication-resource datasets such as
  `E-MTAB-14704`, it resolves declared local CEL paths and an available SDRF
  metadata file from the resource catalog while reporting any missing payloads
  without downloading them. It records filesystem status, CEL
  size/mtime-derived cache keys, parsed metadata previews, default condition
  contrasts, annotation/library readiness, explicit output/cache path status,
  backend-candidate readiness, normalized platform hints, planned outputs, and
  local dependency checks. Repeatable `--r-library-path PATH` values select
  additional R library roots for both preflight and the generated helper
  command. When no value is supplied, an existing workspace `.r-lib` is
  retained as a backward-compatible default; the resolved path is made explicit
  in the normalized request. The report separates
  `request.r_library_paths` from `r_library_paths_checked`, which also includes
  R's effective `.libPaths()`. Command and package probes are bounded to avoid a
  hung `Rscript` blocking preflight; `timed_out`, `probe_failed`, `missing`, and
  `unchecked` remain distinct dependency states. The `r_oligo` backend candidate includes the
  `scripts/probe_regions_oligo.R` helper path and, for explicit RMA/CEL
  requests, a suggested command. Legacy 3' IVT/CDF registry entries expose an
  `r_affy_cdf` candidate with the `scripts/probe_regions_affy.R` helper path;
  readiness depends on local R/`affy`, `limma`, `Biobase`, CDF, CEL, and
  annotation resources. Package probes use `system.file()` and
  `packageVersion()` without attaching namespaces, while the helper remains
  authoritative for executable runtime loading. The
  `affymetrix_power_tools` candidate
  recognizes user-supplied APT library directories/files containing at least
  PGF and CLF files, includes an optional MPS/meta-probesets file when present,
  and reports an explicit `apt-probeset-summarize -a rma-sketch ...` command
  for direct user execution. GENtle still does not run R, APT, package
  installation, downloads, or CEL summarization implicitly.
  When `--output DIR` is supplied, the same report is persisted as
  `DIR/plan.json` so agents and developers can inspect the exact backend
  readiness, dependency versions, cache keys, suggested commands, and planned
  outputs before any external tool is run.
- For `Clariom_D_Human`, the probe-region preflight also reports the local
  Thermo Fisher support ZIP paths in
  `annotation_source.vendor_support_files[]`. Both concise canonical ZIP names
  and browser-preserved `TFS-Assets_LSG_Support-Files_...` download names are
  accepted. These ZIPs are optional local inputs for probe/probeset coordinate
  and annotation development; they are not themselves a
  `gentle.microarray_track_manifest.v1` manifest. DNA-viewer array projection
  continues to consume a prepared manifest plus per-contrast TSVs whose
  coordinate build has been verified.
- `arrays run-probe-region-backend PLAN.json --allow-external-execution`
  (or `arrays run-probe-region-backend --plan PLAN.json --allow-external-execution`)
  consumes the persisted `gentle.probe_region_plan.v1` artifact and returns
  `gentle.probe_region_backend_run.v1`. Execution is explicit: the command
  refuses to run without the allow flag, refuses plans whose recorded preflight
  or selected backend readiness failed, and never downloads or installs CEL
  files, vendor resources, R packages, or APT. On success it captures
  stdout/stderr/exit status, validates the four-file helper-output contract,
  and rewrites `provenance.json` as
  `gentle.probe_region_backend_provenance.v1` with the actual command,
  dependency-version probes, input/output fingerprints, warnings, and errors.
- `arrays inspect-probe-region-output OUTPUT_DIR` consumes a completed explicit
  helper output directory without mutating project state. It validates
  `region_intensity_chrom_order.csv`, optional `sample_table.tsv`,
  `normalized_feature_matrix_manifest.json`, and `provenance.json`, then emits
  row/feature/transcript-cluster counts, sample/condition/logFC columns,
  chromosome and gene previews, bounded preview rows, coordinate/build
  declarations, target levels, artifact paths, provenance hints, warnings,
  and required-column errors. The same report separates `usable` from
  `projection_ready`: output can be readable while genome-anchored projection
  remains blocked until the helper declares compatible coordinate-system and
  genome-build metadata.
- `arrays import-apt-probe-region-output SUMMARY.tsv ANNOTATION.csv OUTPUT_DIR`
  converts an explicit APT summary table plus an explicit annotation/NetAffx
  coordinate table into the same helper-output directory contract. This route
  writes `region_intensity_chrom_order.csv`,
  `normalized_feature_matrix_manifest.json`, and `provenance.json`, then
  returns `gentle.probe_region_apt_import_report.v1` with an embedded
  `gentle.probe_region_output_inspection.v1`. It requires annotation columns
  for probeset/region id, chromosome, start, and stop, and accepts optional
  strand, transcript-cluster, probe-count, and gene-symbol columns. Optional
  `--metadata`, `--condition-column`, and `--sample-column` inputs match APT
  sample columns to conditions and append `mean_log2_*`, `sd_log2_*`, and
  default `log2FC_*` columns for native plotting/projection. When the
  annotation table also declares `probe_id` plus explicit PM probe
  `probe_start`/`probe_stop` coordinate columns, GENtle writes
  `probe_intensity_chrom_order.csv` with probe coordinates, x/y feature
  positions when present, parent probeset ids, and sample/derived value
  columns. Supplying `--probe-intensity TSV` plus optional
  `--probe-id-column NAME` reads true PM probe-level rows from that explicit
  matrix, recomputes condition summaries/logFC against the matrix's sample
  columns when metadata is present, and marks rows as `probe_level_input`.
  Without a probe matrix, the compatibility fallback still writes parent
  probeset-summary values marked as `parent_probeset_summary`. The import
  report and provenance include `probe_intensity_path`,
  `probe_intensity_source`, `probe_intensity_sample_columns`, and
  `missing_probe_intensity_count`. The provenance also records the declared
  backend identity, external-execution policy, local Rscript/APT version
  probes, input path/size/mtime/cheap-SHA1 fingerprints, and the replayable
  GENtle import command. It does not run APT or infer coordinates from vendor
  files by itself.
- `arrays render-probe-region-output-svg OUTPUT_DIR OUTPUT.svg` consumes the
  same completed helper directory and writes
  `gentle.probe_region_output_svg_export.v1` provenance for a deterministic
  SVG plot. The upper panel renders `mean_log2_*` condition intensity tracks;
  the lower panel renders `log2FC_*` tracks. It is still read-only with respect
  to project state and does not run R/APT.
- `arrays render-probe-region-evidence-svg REPORT.json OUTPUT.svg` consumes an
  exported `gentle.probe_region_evidence_interpretation.v2` report and writes a
  deterministic `gentle.probe_region_evidence_svg_export.v1` SVG summary. It
  renders one lane per transcript row, only the overlapped exon ranges and
  junction spans carried by the report, plus parent probeset spans and PM probe
  intervals grouped by parent feature id. This is a review-only
  transcript/exon-geometry constraint visualization: it does not synthesize
  non-overlapped exons, reconstruct full gene models, infer isoform support,
  assess probe sequence specificity, or decide multi-hit status.
  The GUI Splicing Expert consumes the same `junction_spans[]` rows to surface
  junction-targeting array probes as review-only evidence next to transcript
  geometry.
- `arrays project-probe-region-output SEQ_ID OUTPUT_DIR` projects selected or
  all `log2FC_*` helper-output rows into the existing genome-anchored array
  feature model (`gentle_array_*` qualifiers and
  `gentle.microarray_projection_report.v1`). Default `--level probe_region`
  reads `region_intensity_chrom_order.csv`; `--level pm_probe` reads
  `probe_intensity_chrom_order.csv` and only projects rows explicitly marked
  `probe_level_input`, preserving parent probeset ids and intensity-source
  qualifiers on the generated features. The route accepts direct coordinate
  compatibility when the helper output's declared `coordinate_system` or
  `genome_build` matches the target sequence anchor genome id, or an explicit
  `coordinate_projections[]`/`projection_maps[]` interval-map entry whose
  source matches the helper coordinate/build and whose target matches the
  anchor genome id.
- `arrays interpret-probe-region-evidence SEQ_ID` emits
  `gentle.probe_region_evidence_interpretation.v2` by comparing already
  projected `probe_region_output` array features with transcript/exon geometry
  on the target sequence. Optional `--gene LABEL`, `--level all|probe_region|pm_probe`,
  `--min-abs-logfc N`, `--threshold-source TEXT`, `--policy-sha256 SHA256`,
  and `--path FILE` filter/export the read-only report.
  The report records mapping status, overlap transcript ids, ambiguity tags,
  the projected array `platform` on each evidence row when declared,
  per-transcript compatible/constraining counts, and structured
  `transcript_mappings[]` rows with exon ordinals, exon ranges, junction
  spans, overlap base counts, conservative geometry scores, and score-basis
  guardrails. The report-level `coordinate_frame`, `coordinate_system`, and
  `coordinate_chromosome` fields define the shared plotting frame for evidence
  rows, exon ranges, and junction spans; each mapping also preserves
  `local_exon_ranges_1based` for sequence-local auditability. Transcript rows
  include a review-only `review_status` label for unique, shared,
  constraining, or absent geometry; the report explicitly does not infer
  isoform support, probe uniqueness, or biological validation.
  `threshold_source`, `policy_sha256`, and
  `interpretation_request_sha256` content-bind an effective effect gate when
  supplied. Threshold provenance without a threshold is rejected. Existing
  reports without these fields remain readable but may not borrow a threshold
  from a later study plan.
- Nearest-array navigation is a composition of shared contracts rather than a
  separate GUI operation. Query projected intervals with `features query
  SEQ_ID --qual-contains gentle_track_source=Array --nearest-to POSITION
  --limit 1 --include-qualifiers`; then use this interpretation operation when
  transcript/exon/junction geometry is required. Distance alone never implies
  transcript association.
- manifests may also include `coordinate_projections[]` entries with
  `source_genome_id`, `target_genome_id`, `method`, and `path`. These paths
  point at tab-delimited interval maps for build-to-build projection, currently
  intended for audited GRCh37/hg19 to GRCh38/hg38 array-track display.
- per-row TSV fields include `chrom`, `start_1based`, `end_1based`, `strand`,
  `feature_id`, `transcript_cluster_id`, `exon_id`, `probe_type`, `logFC`,
  `AveExpr`, `P.Value`, `adj.P.Val`, and optional junction/gene metadata.
- projection is strict: the target sequence must have `GenomeSequenceAnchor`
  provenance and either the manifest coordinate system/supported aliases must
  match the anchor `genome_id`, or a declared coordinate projection must map
  the manifest build into that anchor build; non-projectable/unverified array
  builds are rejected.
- GENtle streams only rows overlapping the anchored interval and materializes
  them as ordinary `track` features with `gentle_track_source=Array`,
  `gentle_array_dataset`, `gentle_array_platform`, `gentle_array_contrast`,
  `logFC`, `adj_P_Val`, `AveExpr`, `feature_id`,
  `transcript_cluster_id`, and `exon_id` qualifiers.
- projected rows retain both native and display coordinates in qualifiers:
  `gentle_array_native_*` records the original manifest build, while
  `genomic_*` and `chromosome/start_1based/end_1based` record the displayed
  anchor-build interval.

microRNA target-site scan notes:

- `gentle.mirna_target_scan.v1` reports exact canonical microRNA seed-site
  candidates over annotated transcript partitions. It is sequence evidence,
  not automatic biological validation.
- `gentle.mirna_seed_explanation.v1` reports the resolved query microRNA and
  canonical target motif table used by `mirna explain-seed`.
- `gentle.mirna_catalog_record.v1` reports catalog-backed mature-sequence,
  accession, alias, and source metadata used by `mirna catalog-show`.
- The shared shell commands are:
  - `mirna explain-seed MIRNA [--mature-sequence SEQ] [--format json]`
  - `mirna catalog-show MIRNA`
  - `mirna scan-target MIRNA TARGET [--mature-sequence SEQ] [--transcript-filter TEXT] [--regions CSV] [--seed-classes CSV] [--boundary-flank-bp N] [--species-note TEXT] [--evidence-note TEXT] [--format json]`
- Built-in v1 catalog coverage includes `hsa-miR-96-5p` / `MIMAT0000095` /
  NCBI Gene `407053`. Direct mature sequence input is also accepted and, when
  it exactly matches a built-in mature sequence, resolves back to that catalog
  record.
- Canonical target motifs are derived from the mature sequence:
  `8mer`, `7mer-m8`, `7mer-A1`, and `6mer`. For `hsa-miR-96-5p`, the table
  includes `GTGCCAAA`, `GTGCCAA`, `TGCCAAA`, and `TGCCAA`.
- Seed-class counts are intentionally overlapping in v1: an exact 8mer site may
  also appear as a 7mer or 6mer candidate when those classes are requested.
  Reports preserve every requested exact-class match instead of collapsing to
  only the most-specific class.
- Region aliases expand deterministically:
  `3utr`, `exon`, `intron`, and `boundary` become 3-prime UTR,
  coding/noncoding exons, introns, and both splice-boundary directions.
- Boundary scans use a 25 bp default flank on each side of the exon/intron
  junction. This is a compact donor/acceptor neighborhood for first-pass
  splicing-impact review; callers can override it with `--boundary-flank-bp`.
- Reports include grouped hits with transcript id, parent feature, sequence
  local 0-based half-open coordinates, 1-based display coordinates, strand,
  matched sequence, +/-20 bp region context, `evidence_tags`, notes, and
  summary counts for every requested region/seed-class pair.
- `evidence_tags` are typed string values. The v1 set includes
  `exact_seed_candidate`, `orthologous_experimental_context`, and the reserved
  `directly_validated_human_site` tag, which GENtle does not emit unless a
  future direct-human-validation source supplies it.
- For reverse-strand transcripts, `matched_sequence` and
  `region_context_sequence` are reported in transcript/mRNA orientation, while
  local and genomic coordinate fields remain in the displayed sequence's
  coordinate system. This keeps the biological pairing readable without hiding
  where the hit lands in the forward-coordinate viewer.
- The `orthologous_experimental_context` tag can record rat Tp73 PMID 37099528
  context for the hsa-miR-96-5p / TP73 example. GENtle does not emit
  `directly_validated_human_site` unless direct human validation is supplied by
  a future evidence source.

Catalog-backed reference/helper discovery notes:

- shared shell/CLI discovery commands `genomes list` and `helpers list` accept
  optional `--catalog PATH` and `--filter TEXT`
- `PATH` may point to either one JSON catalog file or a directory of top-level
  `*.json` fragments; directory fragments are merged deterministically by sorted
  filename and duplicate entry ids fail fast
- when `catalog_path` is omitted, the engine now resolves a deterministic
  discovery chain in this order:
  - built-in catalog file plus optional built-in fragment directory
  - system overlay file/directory under `/etc/gentle/catalogs/`
  - user overlay file/directory under `$XDG_CONFIG_HOME/gentle/catalogs/` or
    `~/.config/gentle/catalogs/`
  - project overlay file/directory under `PROJECT_ROOT/.gentle/catalogs/`
- the root locations for built-in/system/project discovery may be overridden in
  controlled environments via `GENTLE_ASSET_ROOT`,
  `GENTLE_SYSTEM_CONFIG_ROOT`, and `GENTLE_PROJECT_ROOT`
- adapters that need to preserve "use default discovery" intent through
  persisted operation/provenance records may emit
  `gentle://catalog/reference/default` or
  `gentle://catalog/helper/default`
- list results now include both the stable entry id array and richer `entries`
  metadata rows so frontends, agents, and ClawBio integrations can search and
  display the same catalog facts without re-encoding them
- helper/reference catalog entries may now carry typed discovery metadata such
  as `summary`, `aliases`, `tags`, `search_terms`, `species`, `helper_kind`,
  `host_system`, `sequence_availability`, `redistribution_status`,
  `biological_safety_note`, `usable_as_empty_backbone`, `procurement`, and
  optional structured `semantics`
- metadata-only helper/vector candidates may omit sequence and annotation
  sources when `sequence_availability`, `redistribution_status`, and
  `usable_as_empty_backbone` are present; prepare/status paths still require a
  real source before sequence indexing can proceed
- helper semantics vocabulary overlays are discovered separately from helper
  sequence catalogs so teams can extend normalized meaning without editing
  GENtle source:
  - built-in `assets/helper_semantics_vocabulary.json` plus optional
    `assets/helper_semantics_vocabulary.d/*.json`
  - system/user/project overlays at
    `catalogs/helper_semantics_vocabulary.json` and
    `catalogs/helper_semantics_vocabulary.d/*.json` under the same roots used
    for catalog discovery
  - vocabulary fragments use schema
    `gentle.helper_semantics_vocabulary.v1` and `terms[]` records with
    `axis`, canonical `value`, optional `label`, optional `description`,
    `aliases[]`, and optional `routine_hints[]`
- shared shell/CLI inspection command
  `helpers vocabulary list [--vocabulary PATH] [--filter TEXT]` returns the
  resolved vocabulary terms with deterministic source/provenance fields; MCP,
  JS, and Lua expose the same list through dedicated helpers so agents can
  inspect vocabulary facts before relying on helper interpretation output
- shared shell/CLI validation command
  `helpers vocabulary doctor [--vocabulary PATH] [--routine-catalog PATH]`
  returns `gentle.helper_semantics_vocabulary_doctor.v1` with source ordering,
  parsed fragment SHA-1 digests, duplicate canonical-term warnings, alias
  collision errors, malformed routine-hint diagnostics including missing
  `source_terms` targets, and unknown routine-family warnings marked as local
  extensions rather than authoritative catalog failures
- helper-list/status routes may now also expose an engine-owned normalized
  `interpretation` record derived from those helper fields:
  - `helper_id`
  - `description`, `summary`, `aliases`
  - `helper_kinds`, `host_systems`
  - `offered_functions`, `constraints`
  - `procurement_channels`, `local_variant_unpublished`
  - ontology-friendly `normalized_terms[]` with `axis`, normalized `value`,
    optional source-spelling `label`, derivation `source`, and optional
    vocabulary enrichment fields `vocabulary_label`,
    `vocabulary_description`, `vocabulary_source`
  - direct `routine_hints[]` with routine `family`, deterministic `rationale`,
    and `source_terms`, so planning/routine-ranking/MCP/ClawBio consumers do
    not need to re-derive first-order helper compatibility from prose
  - deterministic `components[]` and `relationships[]`
- that metadata is intended to stay compatible with the emerging
  reasoning/constraint engine and with later ontology-backed helper/vector
  inspection routes
- shared shell/CLI inspection command `helpers doctor-catalog [--catalog PATH]`
  returns `gentle.helper_vector_catalog_doctor.v1` with deterministic
  structured-field issues only; it does not infer biology, provenance, or
  commercial status from prose
- shared shell/CLI inspection command
  `helpers show-card [--catalog PATH] [--filter TEXT|--name TEXT]` returns
  `gentle.helper_vector_card.v1` cards as pure projections of helper catalog
  fields and structured helper semantics components, so GUI/CLI/agent surfaces
  can inspect vector candidates without re-parsing descriptions

Sequencing-trace evidence notes:

- raw traces are stored separately from `SequencingConfirmationReport`
  payloads; importing a trace does not run confirmation and does not mutate any
  sequence entry
- `ImportSequencingTrace` currently auto-detects:
  - ABIF/AB1 via `ABIF` magic bytes
  - SCF via `.scf` magic bytes
- stored `SequencingTraceRecord` payloads preserve:
  - file-supplied called bases
  - called-base confidence arrays when available
  - peak locations when available
  - raw per-channel intensity arrays when available
  - compact per-channel trace-length summaries
  - optional clip window metadata when present in the source file
  - optional sample/run/machine metadata when present in the source file
- trace-aware confirmation now reuses the same `SequencingConfirmationReport`
  model:
  - `ConfirmConstructReads` accepts `trace_ids` in addition to `read_seq_ids`
  - `ConfirmConstructReads` accepts optional `baseline_seq_id` so the expected
    construct remains primary truth while baseline context can distinguish
    intended edits from reference reversions
  - per-evidence rows expose evidence kind plus optional `trace_id`
  - target support/contradiction ids may now refer to imported trace ids when
    traces provide the relevant evidence
  - report payloads now include:
    - `baseline_seq_id?`
    - per-target `expected_bases?` / `baseline_bases?` for expected-edit loci
    - `variants[]` rows with observed allele, evidence id, confidence summary,
      peak center, and classification:
      `expected_match|intended_edit_confirmed|reference_reversion|unexpected_difference|low_confidence_or_ambiguous|insufficient_evidence`
  - persisted confirmation reports now project as lineage analysis artifacts in
    both the GUI lineage workspace and shared `RenderLineageSvg` export:
    nodes are keyed by `report_id`, attach to the expected construct plus
    optional baseline/reference sequence, and reopen the sequencing-confirmation
    specialist on that stored report in GUI adapters
- `SuggestSequencingPrimers { expected_seq_id, primer_seq_ids[], confirmation_report_id?, min_3prime_anneal_bp, predicted_read_length_bp }`
  - non-mutating helper for sequencing-confirmation review and primer coverage
    planning
  - `primer_seq_ids[]` may be empty when `confirmation_report_id` is present:
    that mode proposes fresh sequencing primers for unresolved loci using the
    expected construct plus the saved report context
  - returns `SequencingPrimerOverlayReport` with per-hit orientation, anneal
    span, predicted read span, optional coverage annotations against a
    persisted sequencing-confirmation report, per-problem guidance rows naming
    the best existing primer hit for unresolved targets or variant loci, and
    `proposals[]` rows for fresh primer candidates when no good existing hit is
    available
- `ImportBlastHitsTrack { seq_id, hits[], track_name?, clear_existing?, blast_provenance? }`
  - optional `blast_provenance` payload preserves invocation context
    (`genome_id`, `query_label`, `query_length`, `max_hits`, `task`,
    `blastn_executable`, `blast_db_prefix`, raw `command[]`, `command_line`,
    `catalog_path?`, `cache_dir?`, `options_override_json?`,
    `effective_options_json?`) for sequence-history/audit views.
- `SelectCandidate { input, criterion, output_id? }`
- `ImportIsoformPanel { seq_id, panel_path, panel_id?, strict }`
- `ImportUniprotSwissProt { path, entry_id? }`
- `FetchUniprotSwissProt { query, entry_id? }`
  - `query` is a UniProtKB/Swiss-Prot accession or entry name, for example
    `P04637` or `P53_HUMAN`; `entry_id` is the local GENtle metadata key used
    to store the fetched entry.
- `ImportUniprotEntrySequence { entry_id, output_id? }`
  - imports one first-class protein sequence plus projected UniProt feature
    annotations into regular project sequence state.
- `FetchEnsemblProtein { query, entry_id? }`
  - fetches one Ensembl transcript/protein-backed protein entry from Ensembl
    REST and persists it in `gentle.ensembl_protein_entries.v1`.
  - `query` is an Ensembl protein/translation stable id or accepted Ensembl
    lookup query; it is not a UniProt accession. `entry_id` is the local GENtle
    metadata key used to store the fetched entry.
- `ImportEnsemblProteinSequence { entry_id, output_id? }`
  - imports one stored Ensembl protein entry as a first-class protein sequence
    with imported Ensembl protein-feature annotations.
- `FetchEnsemblGene { query, species?, assembly?, flank_5prime_bp?, flank_3prime_bp?, entry_id? }`
  - fetches one Ensembl gene entry from Ensembl REST and persists it in
    `gentle.ensembl_gene_entries.v1`.
  - for human genes, `query` may be an HGNC-approved symbol such as `FUS` or
    `TP53`, or a stable Ensembl gene id such as `ENSG00000089280`; `species`
    should use Ensembl species names such as `homo_sapiens`; `assembly`, when
    supplied, must match the assembly reported by the lookup response;
    `flank_5prime_bp` and `flank_3prime_bp` request strand-aware genomic context
    around the gene; `entry_id` is the local GENtle metadata key used to store
    the fetched entry.
  - HGNC IDs such as `HGNC:11998` are nomenclature-record identifiers, not
    Ensembl, GenBank, or UniProt accessions; resolve them to an approved symbol
    or linked external accession before calling a database-specific fetch.
- `FetchEnsemblRegion { species, chromosome, start_1based, end_1based, strand?, output_id?, coord_system_version? }`
  - fetches one arbitrary Ensembl REST genomic region/ROI directly as a
    first-class DNA sequence without requiring a prepared local reference.
  - `species` is an Ensembl species name such as `homo_sapiens`; `chromosome`,
    `start_1based`, and `end_1based` are assembly coordinates for that species;
    `output_id` names the local GENtle sequence created from the fetched slice.
  - `strand` defaults to `+`; `-` requests the reverse-strand sequence from
    Ensembl and records `anchor_strand = "-"` in genome-extraction provenance.
  - the imported sequence gets a top-level `source` feature with organism,
    chromosome, genomic bounds, strand, REST source URL, and
    `synthetic_origin=ensembl_region_fetch`.
  - provenance is recorded under the same
    `gentle.provenance.genome_extractions[]` surface used by prepared
    reference extraction, with `sequence_source_type=ensembl_rest_region` and
    no local catalog/cache requirement.
- `ImportEnsemblGeneSequence { entry_id, output_id? }`
  - imports one stored Ensembl gene entry as a first-class DNA sequence with a
    top-level imported `gene` feature and Ensembl provenance qualifiers
  - when the stored Ensembl lookup was fetched with expanded transcript
    content, the import also materializes deterministic `mRNA`, `exon`, and
    `CDS` features with Ensembl transcript/protein identifiers, display names,
    biotype flags, canonical/Gencode-primary flags, and translation lengths.
    Negative-strand Ensembl genes are mapped into Ensembl's gene-oriented
    sequence coordinates so same-strand transcript/CDS features derive proteins
    without an extra reverse-complement step.
- `FetchGenBankAccession { accession, as_id? }`
  - `accession` is an NCBI GenBank nucleotide accession; `as_id` names the local
    GENtle sequence created from the fetched record.
- `FetchDbSnpRegion { rs_id, genome_id, flank_bp?, output_id?, annotation_scope?, max_annotation_features?, catalog_path?, cache_dir? }`
  - `rs_id` is a dbSNP identifier such as `rs9923231`; `genome_id` is a
    prepared GENtle genome catalog id; `output_id` names the local sequence
    extracted around that variant.
- `DeriveTranscriptSequences { seq_id, feature_ids[], scope?, output_prefix? }`
- `PlanExonSkippedIsoform { seq_id, transcript_feature_id, criteria[], plan_id? }`
- `MaterializeExonSkippedIsoform { plan_id, selected_candidate_ids[], output_prefix? }`
- `DeriveProteinSequences { seq_id, feature_ids[], feature_query?, scope?, output_prefix?, report_id? }`
  - this operation is self-sufficient and transcript-first: it does not depend
    on UniProt or any other external protein evidence source to decide what
    protein products exist
  - `protein` and `peptide` molecule labels are treated as one first-class
    protein family for shared import/export handling
  - also persists one `gentle.protein_derivation_report.v1` artifact keyed by
    stable `report_id` with stored `seq_id`, derived protein sequence ids, and
    `op_id` / `run_id` provenance for lineage/reopen paths
- `ReverseTranslateProteinSequence { seq_id, output_id?, speed_profile?, speed_mark?, translation_table?, target_anneal_tm_c?, anneal_window_bp? }`
- `ProjectUniprotToGenome { seq_id, entry_id, projection_id?, transcript_id? }`
  - persists one `gentle.uniprot_genome_projection.v1` artifact with stable
    `projection_id`, upstream `seq_id`/`entry_id`, and stored `op_id` /
    `run_id` provenance for lineage/reopen paths
- `GenerateCandidateSet { set_name, seq_id, length_bp, step_bp, feature_kinds[], feature_label_regex?, max_distance_bp?, feature_geometry_mode?, feature_boundary_mode?, feature_strand_relation?, limit? }`
- `GenerateCandidateSetBetweenAnchors { set_name, seq_id, anchor_a, anchor_b, length_bp, step_bp, limit? }`
- `DeleteCandidateSet { set_name }`
- `UpsertGuideSet { guide_set_id, guides[] }`
- `DeleteGuideSet { guide_set_id }`
- `FilterGuidesPractical { guide_set_id, config?, output_guide_set_id? }`
- `GenerateGuideOligos { guide_set_id, template_id, apply_5prime_g_extension?, output_oligo_set_id?, passed_only? }`
- `ExportGuideOligos { guide_set_id, oligo_set_id?, format: csv_table|plate_csv|fasta, path, plate_format? }`
- `ExportGuideProtocolText { guide_set_id, oligo_set_id?, path, include_qc_checklist? }`
- `ScoreCandidateSetExpression { set_name, metric, expression }`
- `ScoreCandidateSetDistance { set_name, metric, feature_kinds[], feature_label_regex?, feature_geometry_mode?, feature_boundary_mode?, feature_strand_relation? }`
- `FilterCandidateSet { input_set, output_set, metric, min?, max?, min_quantile?, max_quantile? }`
- `CandidateSetOp { op: union|intersect|subtract, left_set, right_set, output_set }`
- `ScoreCandidateSetWeightedObjective { set_name, metric, objectives[], normalize_metrics? }`
- `TopKCandidateSet { input_set, output_set, metric, k, direction?, tie_break? }`
- `ParetoFrontierCandidateSet { input_set, output_set, objectives[], max_candidates?, tie_break? }`
- `UpsertWorkflowMacroTemplate { name, description?, details_url?, parameters[], script }`
- `DeleteWorkflowMacroTemplate { name }`
- `UpsertCandidateMacroTemplate { name, description?, details_url?, parameters[], script }`
- `DeleteCandidateMacroTemplate { name }`
- `FilterByMolecularWeight { inputs, min_bp, max_bp, error, unique, output_prefix? }`
- `FilterByDesignConstraints { inputs, gc_min?, gc_max?, max_homopolymer_run?, reject_ambiguous_bases?, avoid_u6_terminator_tttt?, forbidden_motifs?, unique, output_prefix? }`
- `Reverse { input, output_id? }`
- `Complement { input, output_id? }`
- `ReverseComplement { input, output_id? }`
- `Branch { input, output_id? }`
- `SetDisplayVisibility { target, visible }`
- `SetLinearViewport { start_bp, span_bp }`
- `AnnotatePromoterWindows { input, gene_label?, transcript_id?, upstream_bp=1000, downstream_bp=200, collapse_mode=transcript|gene }`
  - derives strand-aware promoter windows from transcript TSS geometry
  - writes them back as ordinary `promoter` features with explicit generated
    qualifiers (`generated_by`, `promoter_source`, `gene`, `transcript_id`,
    `transcript_count`, `transcript_ids`, `upstream_bp`, `downstream_bp`)
  - exact duplicate promoter spans that differ only by downstream splice
    variation are collapsed at write-back time; the shared feature label keeps
    one promoter symbol and annotates how many transcripts contributed
  - generated promoter windows render distinctly from imported promoter
    features
- `SummarizeVariantPromoterContext { input, variant_label_or_id?, gene_label?, transcript_id?, promoter_upstream_bp=1000, promoter_downstream_bp=200, tfbs_focus_half_window_bp=100, path? }`
  - emits portable record schema `gentle.variant_promoter_context.v1`
  - reports chosen gene/transcript, promoter overlap, signed TSS distance,
    transcript ambiguity status, overlapping evidence rows, effect tags, and
    suggested assay ids
  - reuses TFBS summary logic for the `variant ± half_window` neighborhood when
    TFBS features are already present
- `SummarizeAlternativePromoterComparison { input, gene_label?, transcript_id?, promoter_upstream_bp=1000, promoter_downstream_bp=200, path? }`
  - emits portable record schema `gentle.alternative_promoter_comparison.v1`
  - groups transcript-derived promoter windows by exact DNA span so alternative
    promoter usage can be compared without stacking redundant downstream splice
    variation
  - each grouped row retains transcript-count / transcript-id provenance plus a
    representative transcript/TSS for GUI retargeting back into Promoter design
  - warnings make the transcript-level to DNA-level collapse explicit when
    several transcript TSS interpretations reduce to one genomic promoter span
- `SuggestPromoterReporterFragments { input, variant_label_or_id?, gene_label?, transcript_id?, retain_downstream_from_tss_bp=200, retain_upstream_beyond_variant_bp=500, max_candidates=5, path? }`
  - emits portable record schema `gentle.promoter_reporter_candidates.v1`
  - ranks transcript-aware, strand-aware promoter fragment candidates and marks
    one deterministic default recommendation
- `ListReporterCatalog { catalog_path?, filter?, limit?, path? }`
  - emits `gentle.reporter_catalog_report.v1`
  - validates the local reporter catalog, quarantines rows with missing
    sequence/provenance/license/safety gates, and adds deterministic sequence
    annotations such as length, GC, checksum, CDS sanity, and motif hits
  - reporter records prefer `sequence_sha256`; `sequence_sha1` remains a
    legacy compatibility field for older local catalogs and is not recomputed
    in-process
- `RecommendReporters { constraints, catalog_path?, limit?, path? }`
  - emits `gentle.reporter_recommendation.v1`
  - runs offline from the local catalog: hard constraints eliminate candidates
    first, then deterministic soft scores rank the remaining reporter records
  - rejected candidates remain in the report with machine-readable reasons
  - includes `biological_intent` so agents can distinguish promoter-reporter,
    luciferase, fusion, live, and endpoint reporter selection without parsing
    prose
- `ExportReporterCorpus { catalog_path?, path, format=json|jsonl }`
  - emits `gentle.reporter_corpus_export.v1`
  - writes the annotated, provenance-bearing reporter corpus for local
    retrieval or local AI training/evaluation prep; GENtle does not train a
    model in this V1
- `PlanReporterConstructHandoff { candidate_set_path, candidate_id?, catalog_path?, reporter_constraints?, reporter_backbone_seq_id?, reporter_backbone_load_path?, reference_fragment_seq_id?, alternate_fragment_seq_id?, output_prefix?, path? }`
  - emits `gentle.reporter_construct_handoff.v1`
  - consumes a saved `gentle.promoter_reporter_candidates.v1` JSON report
    rather than re-running promoter-fragment selection
  - includes `biological_intent =
    allele_paired_promoter_luciferase_reporter_handoff` for the V1 synthetic
    biology bridge
  - V1 supports the exact macro template
    `allele_paired_promoter_luciferase_reporter` and defaults reporter
    recommendation to luciferase-class candidates
  - reports typed macro-port readiness plus reporter-backbone resolution, and
    emits explicit follow-up commands for manual extraction, allele
    materialization, backbone loading, macro import, validation, and macro run
  - does not create constructs, fetch live registries, optimize codons, or make
    wet-lab claims
- `MaterializeVariantAllele { input, variant_label_or_id?, allele=reference|alternate, output_id? }`
  - phase-1 scope is single-nucleotide variants only
  - rejects indels, multi-allelic variants, or variants without explicit
    ref/alt qualifiers
  - preserves the variant feature on the derived output while marking the
    materialized allele
- `SetTopology { seq_id, circular }`
- `RecomputeFeatures { seq_id }`
- `SetParameter { name, value }` (purely in-silico project parameter change)

Isoform-panel operation semantics (current):

- `ImportIsoformPanel` loads curated panel resources with schema
  `gentle.isoform_panel_resource.v1` and binds them to one sequence context.
- `strict=true` enforces hard failure when panel transcript mapping fails;
  `strict=false` records warnings and keeps partial mappings.
- `RenderIsoformArchitectureSvg` emits a deterministic multi-mode architecture
  SVG derived from the same expert payload used by GUI/shell inspection.
  - the top panel remains coordinate-true, preserving genomic exon positions
    and introns
  - when CDS ranges are available for mapped transcripts, the top panel uses
    dual coding (faint full transcript exons + colored CDS blocks)
  - the lower panel now adds a compressed transcript/product coupling view:
    exon-chain transcript geometry, shared CDS colors, and isoform-local
    protein axes
  - shared CDS/exon colors are assigned by genomic exon family/position rather
    than by within-lane segment order, so overlapping versions of the same
    locus exon retain one color across different transcripts and proteoforms
  - lower protein rails no longer force one shared amino-acid axis; instead
    each isoform gets a local product axis while still staying linked back to
    the genomic panel through the colored CDS segments above
  - optional `expression_tsv_path` adds a right-side heatmap keyed by
    `isoform_id` and aligned to the rendered isoform rows; when omitted, the
    SVG remains the existing architecture-only export
- Isoform expression TSV v1 is tab-delimited with required header columns
  `isoform_id`, `sample_label`, and `value`.
  - values must be finite and non-negative
  - duplicate `isoform_id` + `sample_label` cells are rejected
  - rows for unknown isoforms are ignored with warnings in
    `expression_matrix.warnings`
  - the rendered matrix uses `Vec<Option<f64>>` values aligned to
    `sample_labels`, with absent cells rendered as missing
- `gentle.isoform_panel_resource.v1` supports optional protein reference-span
  hints per isoform:
  - `reference_start_aa` (1-based inclusive)
  - `reference_end_aa` (1-based inclusive)
  - when present, protein lanes clip domains within this span and project them
    onto the isoform-local product axis (useful for TP53 N-terminus/C-terminus
    class overlays and truncated proteoforms).
- `gentle.isoform_panel_resource.v1` also supports panel-level transcript
  geometry mode:
  - `transcript_geometry_mode: exon|cds` (default `exon`)
  - `cds` renders top-panel lanes from transcript CDS segments when available,
    falling back to exon geometry per transcript if CDS metadata is missing.
- `gentle.isoform_panel_resource.v1` can now carry structured local curation
  metadata at panel and isoform level:
  - `curation.source_kind`: `public_database`, `lab_curated`,
    `literature_curated`, `vendor_curated`, or `mixed`
  - `curation.source_label`, `curation.evidence[]`,
    `curation.validation_tags[]`, `curation.public_database_status`, and
    `curation.notes[]`
  - validation reports expose the panel `curation_source_kind` and
    `curated_isoform_count`, and per-isoform validation rows expose
    `curation_source_kind` plus `validation_tags[]`, so downstream
    GUI/CLI/ClawBio routes can distinguish local knowledge from public
    accession anchors instead of flattening both into one free-text source
    string.
- `gentle.isoform_panel_resource.v1` may also carry panel-level curation
  sidecars:
  - `evidence[]`: upstream or local records that support the panel
    (`evidence_id`, `source_type`, `accession`, optional `url`,
    `sequence_path`, `sequence_length_bp`, `sequence_sha256`, CDS span, and
    retrieval date)
  - `evaluations[]`: stored comparison or trust-assessment records with
    source evidence ids, isoform ids, row-level status/summary text, and small
    key/value metrics
  - validation reports include `evidence_count`, `evaluation_count`, and
    `evaluation_row_count`, and warn about dangling evidence/isoform
    references so curated resources can remain auditable across GUI, CLI, and
    headless workflows.

`LoadFile` import detection semantics (current):

- deterministic probe order: `SnapGene -> GenBank -> EMBL -> FASTA -> XML`
- SnapGene `.dna` files are supported for sequence/topology/features/basic
  notes import through the shared loader
- XML scope: `GBSet/GBSeq` and `INSDSet/INSDSeq` are supported
- other unsupported XML dialects return explicit schema/dialect diagnostics

`ExtendGenomeAnchor` side semantics:

- `side` accepts `five_prime` or `three_prime`.
- Direction is contextual to anchor strand.
- On anchor strand `-`, `five_prime` increases physical genomic position.
- If the anchor genome id is not prepared exactly, the engine can auto-resolve
  to one compatible prepared assembly-family entry (for example `GRCh38.p14`
  -> `Human GRCh38 Ensembl 116`).
- If multiple compatible prepared entries exist, extension fails with a
  deterministic options list so caller/GUI can choose explicitly.
- `prepared_genome_id` can be passed explicitly to force a specific prepared
  cache and bypass compatibility auto-selection.

`VerifyGenomeAnchor` semantics:

- Re-checks one anchored sequence against the selected prepared genome cache at
  recorded coordinates/strand.
- Writes one new provenance entry with `operation = VerifyGenomeAnchor` and
  `anchor_verified = true|false`.
- Returns an in-place state change (`changed_seq_ids`) for the same sequence id
  so GUI/CLI can refresh verification badges/status lines deterministically.

Local `SequenceAnchor` semantics (distinct from genome provenance anchoring):

- `SequenceAnchor` currently supports:
  - `Position { zero_based }`
  - `FeatureBoundary { feature_kind?, feature_label?, boundary, occurrence? }`
- `boundary` accepts `Start`, `End`, or `Middle`.
- This anchor model resolves in-sequence positions and is used for
  in-silico extraction/scoring workflows (`ExtractAnchoredRegion`,
  `GenerateCandidateSetBetweenAnchors`).

Adapter utility contracts (current, non-engine operations):

For narrative/operator guidance on when to use CLI, MCP, Agent Assistant, or an
external coding agent runtime, see:

- `docs/tutorial/01-01_agent_interfaces.md`

- `help [COMMAND ...] [--format text|json|markdown] [--interface ...]`
  - backed by the shared protocol `capability_registry` for built-in help
    rendering; during the migration, `docs/glossary.json` remains the
    compile-time seed for glossary command descriptors
  - `--format text` renders human-readable help
  - `--format json` renders machine-readable help catalog/topic payload; each
    command row includes the matching shared protocol `capability` descriptor,
    so clients can read syntax plus mutation/surfacing/schema metadata from one
    response
  - `--format markdown` renders documentation-ready markdown
  - `--interface` accepts:
    `all|cli-direct|cli-shell|gui-shell|gui-menu|js|lua|mcp`
    (`mcp` currently aliases to shared shell command docs)
  - glossary `interfaces` use `gui-shell` for GUI shell reachability and
    `gui-menu` for first-class GUI menu or command-palette surfacing.
- shared-shell isoform panel routes:
  - `panels import-isoform SEQ_ID PANEL_PATH [--panel-id ID] [--strict]`
  - `panels inspect-isoform SEQ_ID PANEL_ID`
  - `panels render-isoform-svg SEQ_ID PANEL_ID OUTPUT.svg [--expression-tsv PATH]`
  - `panels validate-isoform PANEL_PATH [--panel-id ID]`
- gene isoform evidence inspection reuses the shared feature-expert routes:
  - `inspect-feature-expert SEQ_ID isoform-evidence PANEL_ID [--annotation-release LABEL] [--rna-read-report-id ID]... [--probe-evidence PATH]... [--cdna-est-resource PATH]... [--expression-tsv PATH] [--occupancy-track NAME]... [--qpcr-report-id ID]...`
  - `render-feature-expert-svg SEQ_ID isoform-evidence PANEL_ID [same evidence options] OUTPUT.svg`
  - the report schema is `gentle.gene_isoform_evidence.v2`; inspection is a
    pure read and never creates qPCR assays or changes the sequence
  - legacy `gentle.gene_isoform_evidence.v1` payloads remain readable. Their
    new per-measurement and recommendation fields default empty and should be
    regenerated before contrast-aware interpretation
  - `transcripts[]` keeps biological `exon_family_ids_5_to_3` separate from
    `exon_family_ids_genomic_ascending`, which is essential for minus-strand
    genes
  - assembly-local geometry ids are stable coordinate labels:
    `EXF:{assembly}:{start}-{end}:{strand}` and
    `JCT:{assembly}:{low}-{high}:{strand}`; junction donor/acceptor fields are
    reported separately in transcript orientation
  - `exon_families[]` and `junctions[]` retain annotation-model counts and four
    independent components: specificity, dataset-relative abundance,
    contrast-specific responsiveness, and assayability. Missing evidence is
    `unknown` or `not_evaluated`, never numeric zero, and no aggregate utility
    score is inferred. Every abundance/responsiveness/assayability component
    retains `measurements[]` by evidence id, condition, value, and unit.
    Multiple measurements are never reduced by maximum magnitude; incompatible
    units remain separate and produce a warning instead of a numeric summary
  - `recommendation` is a deterministic rule-based tier (`assay_ready`,
    `evidence_prioritized`, `annotation_candidate`, or `not_evaluated`) plus
    recommended use and evidence ids. It is triage guidance, not a weighted
    biological score
  - `transcript_metrics[]` records exact annotation-derived protein identity
    SHA-256 and predicted unmodified molecular weight when a complete,
    unambiguous translated CDS is available. Ambiguous residues make mass
    unavailable rather than receiving an average residue mass; predicted mass
    is not evidence of protein expression or gel separation
  - `evidence_items[]` records typed source, method, provenance, target ids,
    family ids, and one of `observed`, `candidate`, `constraint_only`,
    `not_evaluated`, or `unknown`. Array-probe overlap is always
    `constraint_only`; it does not establish isoform support
  - human-facing adapters use the shared labels `observed evidence`,
    `candidate association`, `design constraint`, `not evaluated`, and
    `unresolved evidence` for those stable wire values. The labels are an
    interpretation aid only; the schema spellings and evidence status are not
    rewritten by the GUI
  - optional `occupancy_track_names[]` selects already projected BED/BigWig
    track names (`*` selects all overlapping projected tracks). The additive
    `occupancy_lanes[]` records source-specific local/genomic intervals and
    score ranges, while `occupancy_shared_abs_max_score` supplies the common
    SVG scale. Occupancy rows remain locus-level evidence and do not alter the
    abundance, responsiveness, specificity, or assayability components
  - persisted RNA-read and qPCR reports are resolved by report id. Probe
    interpretation JSON, expression TSV, and cDNA/EST resources remain
    explicit file inputs because those report stores do not yet exist
  - RNA/probe reports for another sequence, qPCR reports for another template,
    and cDNA/EST records with conflicting assembly/chromosome/strand provenance
    are warned about and not attached
  - cDNA/EST files use `gentle.cdna_est_evidence_resource.v1` with a
    `resource_id`, source, optional annotation release, and `records[]` carrying
    kind (`cdna`, `est`, `curated_transcript`, or `other`), accession/coordinates,
    exon-family or junction ids, support count, optional alignment fractions,
    and notes
  - provenance includes source ids plus file paths and SHA-256 digests. A
    report composes sequence evidence; it does not claim that an isoform is
    biologically validated
- gene-locus evidence composition reuses that ledger without changing it:
  - `inspect-feature-expert SEQ_ID gene-locus-evidence PANEL_ID [isoform-evidence options] [--probe-effect-table PATH]... [--probe-effect-contrast TOKEN]... [--probe-effect-coordinate-system ID] [--upstream-bp N] [--downstream-bp N] [--occupancy-layout JSON_OR_@FILE | --occupancy-track NAME ...] [--motif TOKEN]... [--score-kind KIND] [--motif-threshold N] [--motif-top-hits N] [--allow-negative]`
  - `render-feature-expert-svg SEQ_ID gene-locus-evidence PANEL_ID [same options] OUTPUT.svg`
  - the pure-read result schema is `gentle.gene_locus_evidence_display.v1`.
    It embeds the `gentle.gene_isoform_evidence.v2` ledger and adds
    `transcript_metrics[]`, annotation-backed `codon_markers[]`, optional
    `probe_effect_overlays[]`, grouped `occupancy_groups[]`, continuous
    `motif_tracks[]`, deduplicated junction `assay_overlays[]`, a combined
    `provenance[]` inventory, strand-aware locus/flank coordinates, and warnings
  - probe-effect tables are tab-separated and retain PSR intervals and JUC
    junction-edge geometry as distinct classes. Abundance columns follow the
    `log2_mean_*` convention and differential-activity columns follow
    `log2_*_minus_*`; repeat `--probe-effect-contrast` to select columns by id,
    display label, or source-column spelling. The explicit
    `--probe-effect-coordinate-system` must match the open sequence's genome
    anchor before rows are projected. The report carries the original feature
    ids, coordinates, contrast ids, values, PM-probe counts, row provenance,
    and separate abundance and differential display scales. These values are
    visualization evidence, not primer-binding proof, significance estimates,
    or direct isoform support
  - occupancy layout files use `gentle.gene_locus_occupancy_layout.v1`.
    Each group declares `group_id`, label, `scale_mode`
    (`shared_group`, `shared_all`, `independent`, or `fixed`), optional fixed
    scale/comparability rationale, and lanes with exact projected track name,
    optional display/condition labels, and role (`experimental`, `gfp_control`,
    `input_control`, `igg_control`, `positive_control`, `negative_control`, or
    `other`)
  - GENtle never infers cell line, condition, lane role, or cross-group
    comparability from a filename. `shared_all` without an explicit
    `cross_group_scale_justification` is retained but warned about
  - motif tracks reuse the existing TFBS scorer and active local JASPAR
    registry. A threshold controls labeled top hits only; the machine-readable
    continuous vectors remain in the report. Motif/occupancy colocation is not
    promoted to binding or regulation
  - transcript/CDS metrics reuse transcript derivation. Start/stop glyphs are
    emitted only from annotated CDS translations; noncoding/no-CDS transcripts
    are not mislabeled as incomplete CDS
  - qPCR assays are display-deduplicated by stable junction id while retaining
    every contributing assay id, transcript family, and transcript id. No assay
    is designed or persisted by this read-only route
  - the GUI `Locus figure` tab is a thin client of this same report and shared
    SVG renderer. The deterministic offline example and synthetic fixture are
    `docs/examples/workflows/patz1_gene_locus_evidence_offline.json` and
    `test_files/fixtures/gene_locus_evidence/patz1_offline_composer/`
- shared-shell UniProt routes:
  - `uniprot fetch QUERY [--entry-id ID]`
    - `QUERY` is a UniProtKB/Swiss-Prot accession or entry name, for example
      `P04637` or `P53_HUMAN`; `--entry-id` names the local GENtle metadata
      entry.
  - `uniprot import-swissprot PATH [--entry-id ID]`
    - `gentle.uniprot_entry.v1` preserves the parsed `CC   -!- ALTERNATIVE
      PRODUCTS` declaration (named count, names/synonyms, IsoIds, displayed or
      variant sequence status, and parse warnings) beside the raw record.
      `gentle.uniprot_projection_audit.v1` copies that complete inventory into
      `protein_isoform_inventory`; entries without named alternatives emit one
      explicit canonical target. The inventory is distinct from Ensembl link,
      current transcript geometry, and CDS/peptide concordance rows.
  - `uniprot list`
  - `uniprot show ENTRY_ID`
  - `uniprot map ENTRY_ID SEQ_ID [--projection-id ID] [--transcript ID]`
  - `uniprot projection-list [--seq SEQ_ID]`
  - `uniprot projection-show PROJECTION_ID`
  - `uniprot feature-coding-dna PROJECTION_ID FEATURE_QUERY [--transcript ID] [--mode genomic_as_encoded|translation_speed_optimized|both] [--speed-profile human|mouse|yeast|ecoli]`
  - `uniprot resolve-ensembl-links PROJECTION_ID [--transcript ID]`
  - `uniprot transcript-accounting PROJECTION_ID [--transcript ID]`
  - `uniprot compare-ensembl-exons PROJECTION_ID [--transcript ID] [--ensembl-entry ID]`
  - `uniprot compare-ensembl-peptide PROJECTION_ID [--transcript ID] [--ensembl-entry ID]`
  - `uniprot build-linked-transcript-inventory REQUEST_JSON_OR_@FILE`
    emits `gentle.uniprot_linked_transcript_inventory.v1`. Its request binds
    `inventory_id`, `assembly`, `annotation_release`, `source_resource_id`,
    `transcript_fasta_paths[]`, explicit `links[]` (`entry_id`, `isoform_id`,
    versioned `transcript_id`, optional locus/reference accession), and
    `output_path`. Source files and resolved mature cDNAs are SHA-256 bound.
    Missing and ambiguous sequences remain typed records without a cDNA digest;
    they cannot authorize a coverage claim. Exact cDNA identity is deliberately
    separate from locus identity and therefore does not establish genomic
    specificity.
  - `uniprot audit-projection PROJECTION_ID [--transcript ID] [--ensembl-entry ID] [--report-id ID]`
  - `uniprot audit-parity PROJECTION_ID [--transcript ID] [--ensembl-entry ID] [--report-id ID]`
  - `uniprot audit-list|show|export ...`
  - `uniprot audit-parity-list|show|export ...`
- shared-shell Ensembl protein routes:
  - `ensembl-protein fetch QUERY [--entry-id ID]`
    - `QUERY` is an Ensembl protein/translation stable id or accepted Ensembl
      lookup query; `--entry-id` names the local GENtle metadata entry. This is
      not the same identifier space as UniProt accessions.
  - `ensembl-protein list`
  - `ensembl-protein show ENTRY_ID`
  - `ensembl-protein import-sequence ENTRY_ID [--output-id ID]`
- shared-shell protease catalog routes:
  - `proteases list [--filter TEXT] [--output PATH]`
  - `proteases show QUERY [--output PATH]`
  - `proteases digest SEQ_ID PROTEASE[,PROTEASE...] [--output-prefix PREFIX] [--min-length-aa N] [--predict-only]`
  - `proteases digest-gel-svg SEQ_ID PROTEASE[,PROTEASE...] OUTPUT.svg [--min-length-aa N] [--ladder NAME] [--ladders CSV]`
  - semantics:
    - exposes the built-in protease catalog as deterministic JSON
    - search spans names, aliases, cleavage-pattern notation, specificity, and
      curated biotech-facing application tags
    - `ProteaseDigestProteinSequence` is the canonical engine operation for
      cleavage prediction/materialization on first-class protein or peptide
      sequences
    - `RenderProteaseDigestGelSvg` is the graphical companion for digest
      reports and keeps the same source-neutral report payload while producing
      an SVG suitable for ClawBio PNG rasterization
    - digest reports use `gentle.protease_digest_report.v1` and retain
      source-protein, transcript, derivation-mode, and translation-table
      provenance when present on transcript-derived proteins
    - peptide products are materialized as first-class `peptide` sequences
      unless the shell caller passes `--predict-only`
- shared-shell Ensembl gene routes:
  - `ensembl-gene fetch QUERY [--species NAME] [--assembly NAME] [--flank-bp N|--flank-5p-bp N --flank-3p-bp N] [--entry-id ID]`
    - `QUERY` is an approved gene symbol for the requested species (for human,
      an HGNC-approved symbol such as `FUS` or `TP53`) or a stable Ensembl gene
      id; `--entry-id` names the local GENtle metadata entry. `--assembly`
      verifies the resolved assembly rather than silently accepting a different
      build. `--flank-bp` applies the same strand-aware flank at both gene ends;
      the two directional options allow asymmetric context.
    - HGNC IDs such as `HGNC:11998` are stable nomenclature-record identifiers,
      not GenBank, UniProt, or Ensembl accessions. Resolve them to an approved
      symbol or linked database accession before using a database-specific fetch
      route.
  - `ensembl-gene list`
  - `ensembl-gene show ENTRY_ID`
  - `ensembl-gene import-sequence ENTRY_ID [--output-id ID]`
- shared-shell Ensembl region route:
  - `ensembl-region fetch SPECIES CHR START END [--strand +|-] [--output-id ID] [--coord-system-version VERSION]`
    - `SPECIES` is an Ensembl species name such as `homo_sapiens`; `CHR`,
      `START`, and `END` are assembly coordinates; `--output-id` names the
      local GENtle sequence created from the fetched interval.
  - equivalent compact form:
    `ensembl-region fetch SPECIES CHR:START..END[:STRAND] [--output-id ID]`
- shared-shell inline sequence creation route:
  - `sequence create --sequence-text DNA [--output-id ID] [--name TEXT] [--topology linear|circular]`
  - this is the mutating companion to read-only inline scans such as
    `features restriction-scan --sequence-text ...`; it creates a reusable
    project sequence before later commands operate on the same sequence id
- shared feature-expert route now also accepts transcript-first protein
  comparison with optional stored external evidence plus persisted UniProt
  projections as direct targets:
  - `inspect-feature-expert SEQ_ID protein-comparison [--transcript TRANSCRIPT_ID] [--ensembl-entry ENTRY_ID] [--feature-key KEY]... [--feature-key-not KEY]...`
  - `render-feature-expert-svg SEQ_ID protein-comparison [--transcript TRANSCRIPT_ID] [--ensembl-entry ENTRY_ID] [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
  - semantics:
    - derive transcript CDS/protein products directly from the current sequence
      state
    - do not require any stored UniProt projection or other external protein
      evidence
    - optional `--transcript` narrows the compare window to one transcript id
      while preserving the same source-neutral payload shape used by
      external-opinion-backed protein experts
    - optional `--ensembl-entry` resolves one persisted
      `gentle.ensembl_protein_entries.v1` record and layers its protein
      sequence/features onto the same transcript-first compare payload as an
      external protein opinion
    - optional `--feature-key` / `--feature-key-not` filters apply equally to
      transcript-first derived views with external protein opinions, so noisy
      imported Ensembl feature classes can be trimmed without changing the
      product geometry
  - `inspect-feature-expert SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]...`
  - `render-feature-expert-svg SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
  - semantics:
    - resolve one persisted `gentle.uniprot_genome_projections.v1` record
    - build a shared `IsoformArchitectureView`
    - transcript-native CDS/protein derivation is authoritative for transcript
      and product geometry
    - the persisted UniProt projection is consumed as one optional external
      protein opinion layered onto that transcript-native product model
    - each transcript/protein lane pair now clips to that transcript's mapped
      amino-acid coverage and projected feature spans, so truncated isoforms do
      not inherit full-length reference domains they do not encode
    - rendered expert SVGs now preserve both:
      - true genomic exon/intron position in the top panel
      - transcript/product geometry in a lower panel whose transcript columns
        are shared by genomic exon family while the protein rails stay
        isoform-local
      so skipped-exon proteoform differences become readable without throwing
      away locus context
    - if transcript features lack explicit `cds_ranges_1based`, the shared
      projection path now prefers compatible `CDS` features before falling back
      to whole-exon spans
    - `FeatureExpertTarget::UniprotProjection` now carries a shared
      `ProteinFeatureFilter { include_feature_keys[], exclude_feature_keys[] }`
    - default filtering hides UniProt `CONFLICT` features unless the caller
      explicitly re-includes them
    - topology/membrane-style UniProt features (`SIGNAL`, `TRANSIT`,
      `TOPO_DOM`, `TRANSMEM`, `INTRAMEM`) render in a dedicated lower band
      beneath the protein rail instead of competing with the standard domain
      label overlay
    - per-transcript compare status is now source-neutral:
      - `derived_only`
      - `consistent_with_external_opinion`
      - `low_confidence_external_opinion`
      - `no_transcript_cds`
      - `external_opinion_only`
    - inside the engine, external protein evidence now enters the Protein
      Expert through one source-neutral adapter boundary before view assembly,
      so future sources can reuse the same transcript-first comparison payload
      instead of forking the renderer/GUI contract
    - this comparison model is intentionally not UniProt-specific:
      Ensembl protein entries now populate the same fields through the same
      adapter boundary, and future providers should continue to reuse that
      transcript-first contract rather than replacing transcript-native
      translation
    - the same persisted projection now also appears in GUI lineage as one
      analysis artifact node linked from the source sequence and reopenable
      through the protein expert
- Protein residue genomic coordinate query semantics:
  - shared engine operation:
    `QueryProteinResidueGenomicCoordinates { seq_id, transcript_id?, residue_start_1based, residue_end_1based }`
  - shared-shell/CLI shorthand:
    `transcripts residue-genomic-coordinates SEQ_ID RESIDUE_START [RESIDUE_END]`
  - both routes resolve to one
    `gentle.protein_residue_genomic_coordinates.v1` report
  - optional `--transcript ID` narrows the query by transcript id, label, or
    `n-N` feature id
  - when executed as an engine operation, the report is exposed in
    `OpResult.protein_residue_genomic_coordinates`
  - each match reports the amino-acid residue, coding-orientation codon,
    per-base 1-based genomic positions, exon ordinals, reverse-strand coding
    order, and whether the codon crosses an exon junction
  - `transcript_feature_id` follows the same zero-based source-feature indexing
    used by adjacent transcript/projection reports so GUI/CLI reopen paths can
    reuse it directly
  - the query is transcript-native: it uses the same CDS resolution and
    translation-table selection as derived transcript/protein products, and it
    does not require UniProt or Ensembl protein evidence
- UniProt feature-coding DNA query semantics:
  - resolves one persisted `gentle.uniprot_genome_projection.v1` record
  - matches `FEATURE_QUERY` case-insensitively against mapped UniProt feature
    key/note text
  - returns one `gentle.uniprot_feature_coding_dna_query.v1` report
  - each transcript match includes:
    - `genomic_coding_dna`: spliced coding-strand DNA exactly as encoded in the
      current genome sequence
    - `translation_speed_optimized_dna`: optional preferred-codon alternative
      using the selected or inferred `TranslationSpeedProfile`
    - `exon_spans[]` and `exon_pairs[]` so GUI/CLI can report the exon or exon
      pair carrying the feature
  - exon ordinals follow transcript order, not raw genomic left-to-right
    position; reverse-strand transcript exon 1 is the transcript 5' exon
- shared-shell GenBank route:
  - `genbank fetch ACCESSION [--as-id ID]`
    - `ACCESSION` is an NCBI GenBank nucleotide accession; `--as-id` names the
      local GENtle sequence id created from the fetched record.
- shared-shell dbSNP route:
  - `dbsnp fetch RS_ID GENOME_ID [--flank-bp N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--catalog PATH] [--cache-dir PATH]`
    - `RS_ID` is a dbSNP identifier such as `rs9923231`; `GENOME_ID` is a
      prepared GENtle genome catalog id; `--output-id` names the local sequence
      extracted around the variant.
- shared-shell protocol-cartoon routes:
  - `protocol-cartoon list`
  - `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
  - `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
  - `protocol-cartoon template-validate TEMPLATE.json`
  - `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
  - `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
  - command surface is intentionally canonical: protocol-cartoon routes do not
    expose extra alias names

- Python adapter wrapper (`integrations/python/gentle_py`):
  - thin subprocess-based wrapper over `gentle_cli`
  - deterministic methods:
    - `capabilities()`
    - `state_summary()`
    - `op(operation)`
    - `workflow(workflow|workflow_path)`
    - `shell(line, expect_json=False)`
    - `render_dotplot_svg(seq_id, dotplot_id, output_svg, ...)`
  - raises structured `GentleCliError` with:
    - `code` (best-effort extracted stable code token)
    - `command`, `exit_code`, `stdout`, `stderr`
  - executable resolution order:
    1. constructor `cli_cmd`
    2. `GENTLE_CLI_CMD`
    3. `gentle_cli` on `PATH`
    4. repository fallback `cargo run --quiet --bin gentle_cli --`

- `gentle_mcp` (stdio MCP adapter, expanded UI-intent parity baseline)
  - MCP role:
    - request/response transport for tool execution (`tools/call`)
    - standardized capability discovery/negotiation (`tools/list`,
      `capabilities`, `help`)
  - current tools:
    - `capabilities`
    - `state_summary`
    - `runtime_status`
    - `restriction_site_detail` (shared restriction-site expert detail record)
    - `exon_skip_plan` (shared `transcripts exon-skip-plan` contract)
    - `exon_skip_materialize` (shared `transcripts exon-skip-materialize`
      contract; requires explicit `confirm=true`)
    - `op` (apply one `Operation`; requires explicit `confirm=true`)
    - `workflow` (apply one `Workflow`; requires explicit `confirm=true`)
    - `help`
    - `reference_catalog_entries` (shared `genomes list` catalog contract)
    - `helper_catalog_entries` (shared `helpers list` catalog contract)
    - `host_profile_catalog_entries` (shared `hosts list` catalog contract)
    - `ensembl_installable_genomes` (shared Ensembl discovery contract for currently installable candidates)
    - `construct_reasoning_graphs` (shared `construct-reasoning list-graphs` inspection contract)
    - `construct_reasoning_graph` (shared `construct-reasoning show-graph` inspection contract)
    - `helper_interpretation` (direct helper-construct interpretation lookup)
    - `ui_intents` (shared `ui intents` catalog)
    - `ui_intent` (shared `ui open|focus|close ...` resolution path)
    - `ui_prepared_genomes` (shared `ui prepared-genomes ...` query path)
    - `ui_latest_prepared` (shared `ui latest-prepared ...` query path)
  - successful mutating calls (`op`, `workflow`) persist state to the resolved
    `state_path`
  - UI-intent tools route through the shared shell parser/executor
    (`parse_shell_tokens` + `execute_shell_command_with_options`) and are
    required to remain non-mutating (`state_changed = false`)
  - tool handlers are adapter wrappers over existing deterministic engine/shell
    contracts (no MCP-only biology logic branch)
  - stdio framing/validation hardening:
    - `Content-Length` is required, duplicate headers are rejected
    - maximum accepted frame size is `8 MiB`
    - parsed JSON nesting depth is capped at `96`
    - `tools/call` params are strict (`name`, optional `arguments` only)
    - `tools/call.arguments` must be a JSON object

MCP query/introspection tool contracts (current):

- `runtime_status`
  - arguments:
    - optional `state_path`
  - behavior:
    - returns `gentle.runtime_status.v1`, with process-local live `frames[]`
      for the MCP process plus observed `activities[]` from existing
      persisted/project activity ledgers
    - mirrors shared shell `introspect runtime`; it does not write a
      runtime-status file

- `agent_systems`
  - arguments:
    - optional: `catalog_path`
  - behavior:
    - returns the same structured payload shape as shared shell
      `agents list [--catalog ...]`

- `agent_preflight`
  - arguments:
    - required: `system_id`
    - optional: `catalog_path`, `live`, `base_url`, `model`, `timeout_secs`,
      `connect_timeout_secs`, `read_timeout_secs`, `max_retries`,
      `max_response_bytes`
  - behavior:
    - returns the same structured payload shape as shared shell
      `agents preflight ...`
    - `live=true` routes through shared shell as `agents preflight --live`
      and includes the optional `live_probe` object

- `agent_models`
  - arguments:
    - required: `system_id`
    - optional: `catalog_path`, `base_url`
  - behavior:
    - returns the same structured payload shape as shared shell
      `agents discover-models ...`

- `agent_plan`
  - arguments:
    - required: `system_id`, `prompt`
    - optional: `state_path`, `catalog_path`, `base_url`, `model`,
      `timeout_secs`, `connect_timeout_secs`, `read_timeout_secs`,
      `max_retries`, `max_response_bytes`, `include_state_summary`,
      `max_candidates`, `allow_mutating_candidates`
  - behavior:
    - returns the same structured payload shape as shared shell
      `agents plan ...`

- `agent_execute_plan`
  - arguments:
    - required: `candidate_id`
    - one of: `plan` or `plan_path`
    - optional: `state_path`, `confirm`
  - behavior:
    - returns the same structured payload shape as shared shell
      `agents execute-plan ...`

- `restriction_site_detail`
  - arguments:
    - required: `seq_id`, `cut_pos_1based`
    - optional: `state_path`, `enzyme`, `recognition_start_1based`,
      `recognition_end_1based`
  - behavior:
    - returns the same structured payload shape as shared shell
      `inspect-feature-expert SEQ_ID restriction CUT_POS_1BASED ...`
    - remains non-mutating and rejects any shared-shell path that reports
      `state_changed=true`
  - result:
    - `kind = "restriction_site"`
    - `data` is the shared `RestrictionSiteExpertView`, including
      `tooltip_lines[]` for GUI-style hover/popover summaries

- `reference_catalog_entries`
  - arguments:
    - optional: `catalog_path`, `filter`
  - behavior:
    - returns the same structured payload shape as shared shell
      `genomes list [--catalog ...] [--filter ...]`
  - result:
    - `catalog_path`, `filter`, `genome_count`, `genomes[]`, `entries[]`

- `helper_catalog_entries`
  - arguments:
    - optional: `catalog_path`, `filter`
  - behavior:
    - returns the same structured payload shape as shared shell
      `helpers list [--catalog ...] [--filter ...]`
  - result:
    - `catalog_path`, `filter`, `genome_count`, `genomes[]`, `entries[]`
    - helper `entries[]` may carry normalized `interpretation` records

- `helper_semantics_vocabulary`
  - arguments:
    - optional: `vocabulary_path`, `filter`
  - behavior:
    - returns the same structured payload shape as shared shell
      `helpers vocabulary list [--vocabulary ...] [--filter ...]`
  - result:
    - `vocabulary_path`, `filter`, `term_count`, `terms[]`
    - each term carries `axis`, canonical `value`, optional `label`,
      optional `description`, `aliases[]`, `routine_hints[]`, and `source`

- `host_profile_catalog_entries`
  - arguments:
    - optional: `catalog_path`, `filter`
  - behavior:
    - returns the same structured payload shape as shared shell
      `hosts list [--catalog ...] [--filter ...]`
  - result:
    - `catalog_path`, `filter`, `profile_count`, `profile_ids[]`, `entries[]`

- `construct_reasoning_graphs`
  - arguments:
    - optional: `state_path`, `seq_id`
  - behavior:
    - returns the same structured payload shape as shared shell
      `construct-reasoning list-graphs [SEQ_ID]`
  - result:
    - graph summary list payload with shared summary rows for stored graphs

- `construct_reasoning_graph`
  - arguments:
    - required: `graph_id`
    - optional: `state_path`
  - behavior:
    - returns the same structured payload shape as shared shell
      `construct-reasoning show-graph GRAPH_ID`
  - result:
    - full portable graph payload plus the compact shared summary block used by
      CLI/GUI adapter-facing inspection surfaces

- `construct_reasoning_set_annotation_status`
  - arguments:
    - required: `graph_id`, `annotation_id`, `editable_status`
    - optional: `state_path`
  - behavior:
    - runs the same mutating shared shell command as
      `construct-reasoning set-annotation-status GRAPH_ID ANNOTATION_ID STATUS`
    - persists the updated project state when the command changes it
  - result:
    - updated graph payload
    - updated `annotation_candidate`
    - same compact shared summary block used by `construct_reasoning_graph`

- `construct_reasoning_write_annotation`
  - arguments:
    - required: `graph_id`, `annotation_id`
    - optional: `state_path`
  - behavior:
    - runs the same mutating shared shell command as
      `construct-reasoning write-annotation GRAPH_ID ANNOTATION_ID`
    - materializes one accepted or locked generated annotation candidate as an
      ordinary sequence feature when it is eligible
    - persists the updated project state when the command changes it
  - result:
    - refreshed graph payload
    - refreshed or already-backed `annotation_candidate`
    - `writeback` report in `gentle.annotation_candidate_writeback.v1` shape
    - same compact shared summary block used by `construct_reasoning_graph`

- `ensembl_installable_genomes`
  - arguments:
    - optional: `collection`, `filter`
  - behavior:
    - returns the same Ensembl discovery report used by GUI/CLI/JS/Lua for
      answering which genomes currently look installable because both FASTA and
      GTF species-directory listings are present
  - result:
    - `collection_filter`, `availability_basis`,
      `collection_latest_releases{}`, `candidates[]`, `warnings[]`

Shell/engine quick-install contracts:

- `genomes install-ensembl SPECIES_DIR ...`
- `helpers install-ensembl SPECIES_DIR ...`
  - behavior:
    - resolve current Ensembl FASTA/GTF files for one species directory
    - derive or accept an explicit target `genome_id`
    - write a real catalog entry before preparation starts
    - choose between two safe write modes:
      - `full_catalog`: update one writable JSON catalog file in place or as a standalone copy
      - `overlay_entry`: write only the new entry into an overlay fragment/file so default discovery does not duplicate built-in ids
    - run the existing prepare pipeline after the catalog write succeeds
  - result:
    - `preview.collection`, `preview.species_dir`, `preview.display_name`
    - `preview.file_stem`, `preview.release`
    - `preview.genome_id`
    - `preview.output_catalog_path`
    - `preview.catalog_write_mode`
    - `preview.catalog_entry_action`
    - `preview.sequence_remote`, `preview.annotations_remote`
    - `prepare_report`

- `helper_interpretation`
  - arguments:
    - required: `helper_id` (id or alias)
    - optional: `catalog_path`
  - behavior:
    - resolves one helper entry through the shared catalog/alias lookup logic
      and returns the normalized helper-construct interpretation if semantics
      are available
  - result:
    - `query`
    - `catalog_path`
    - `interpretation` (`null` when the entry exists but carries no structured
      helper semantics)

- `ui_intents`
  - arguments:
    - `state_path?` (optional; accepted for interface symmetry)
  - behavior:
    - executes shared shell command: `ui intents`
  - result:
    - structured payload schema: `gentle.ui_intents.v1`
    - includes stable `targets`, `target_details`, `target_metadata`,
      `commands`, and deterministic notes
    - `target_details[]` now carries per-target discoverability metadata:
      `title`, `detail`, `keywords`, `menu_path`, supported `actions`, and
      stable `optional_arguments`
    - `target_details[].arguments[]` is additive richer metadata; each entry
      carries `name`, `required`, and `detail` while `optional_arguments`
      remains available for legacy consumers
    - `target_metadata[]` is retained as a compatibility alias for older
      command-catalog clients that consumed only title/detail/keyword/action
      records

- `ui_intent`
  - arguments:
    - required: `action` (`open|focus|close`), `target`
    - optional: `state_path`, `seq_id`, `item_id`, `section`, `genome_id`,
      `helpers`, `catalog_path`, `cache_dir`, `filter`, `species`, `latest`
  - current stable targets:
    - `sequence-window` (requires `seq_id`; controls a loaded DNA sequence
      viewer without mutating the sequence record)
    - `recent-project` (requires the opaque GUI-host `item_id`; open only)
    - `tutorial-project` (requires the catalog `item_id`/chapter id; open only)
    - `configuration` (optional `section`: `external-applications`,
      `agent-systems`, `microarrays`, `graphics`, or `language`)
    - `prepared-references`
    - `prepare-reference-genome`
    - `retrieve-genome-sequence`
    - `blast-genome-sequence`
    - `import-genome-track`
    - `pcr-design`
    - `sequencing-confirmation`
    - `agent-assistant`
    - `prepare-helper-genome`
    - `retrieve-helper-sequence`
    - `blast-helper-sequence`
  - behavior:
    - executes shared shell command:
      - `ui open TARGET ...`, `ui focus TARGET ...`,
        `ui open recent-project ITEM_ID`,
        `ui open tutorial-project CHAPTER_ID`,
        `ui open|focus configuration [SECTION]`,
        `ui open sequence-window SEQ_ID`, `ui focus sequence-window SEQ_ID`,
        `ui close TARGET`, or `ui close sequence-window SEQ_ID`
    - for `target = prepared-references`, optional query flags can resolve
      `selected_genome_id` deterministically through the same helper path used
      by shared shell/CLI
    - parser guardrails are preserved:
      - query flags (`--helpers`, `--catalog`, `--cache-dir`, `--filter`,
        `--species`, `--latest`) are rejected for non-`prepared-references`
        targets
  - result:
    - generic structured payload schema: `gentle.ui_intent.v1`
    - specialized schemas: `gentle.ui_recent_project_intent.v1`,
      `gentle.ui_tutorial_project_intent.v1`, and
      `gentle.ui_configuration_intent.v1`
    - fields include `ui_intent`, `selected_genome_id`, optional
      `prepared_query`, `applied=false`, and deterministic `message`

- `ui_prepared_genomes`
  - arguments:
    - optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`, `filter`,
      `species`, `latest`
  - behavior:
    - executes shared shell command: `ui prepared-genomes ...`
  - result:
    - structured payload schema: `gentle.ui_prepared_genomes.v1`
    - includes `prepared_count`, sorted `genomes[]`, and `selected_genome_id`

- `ui_latest_prepared`
  - arguments:
    - required: `species`
    - optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`
  - behavior:
    - executes shared shell command: `ui latest-prepared SPECIES ...`
  - result:
    - structured payload schema: `gentle.ui_latest_prepared.v1`
    - includes `selected_genome_id` and nested `prepared_query` payload

MCP UI-intent JSON-RPC example (abbreviated):

```json
{
  "jsonrpc": "2.0",
  "id": 7,
  "method": "tools/call",
  "params": {
    "name": "ui_intent",
    "arguments": {
      "action": "open",
      "target": "prepared-references",
      "catalog_path": "assets/genomes.json",
      "species": "human",
      "latest": true
    }
  }
}
```

Result envelope shape:

```json
{
  "jsonrpc": "2.0",
  "id": 7,
  "result": {
    "isError": false,
    "structuredContent": {
      "schema": "gentle.ui_intent.v1",
      "selected_genome_id": "Human GRCh38 Ensembl 116",
      "applied": false
    }
  }
}
```

Adapter-equivalence guarantee for UI-intent tools:

- deterministic parity tests compare MCP UI-intent tool outputs with direct
  shared shell `ui ...` command outputs for:
  - intent catalog (`ui_intents`)
  - prepared query (`ui_prepared_genomes`)
  - latest helper (`ui_latest_prepared`)
  - open/focus/close intent resolution (`ui_intent`)

- `macros run/instance-list/instance-show/template-list/template-show/template-put/template-delete/template-import/template-run`
  - shared-shell macro adapter family for full operation/workflow scripting
  - template persistence is backed by engine operations
    `UpsertWorkflowMacroTemplate`/`DeleteWorkflowMacroTemplate`
  - `template-put` supports optional typed port contracts:
    - `--input-port PORT_ID:KIND[:one|many][:required|optional][:description]`
    - `--output-port PORT_ID:KIND[:one|many][:required|optional][:description]`
  - `template-import PATH` accepts:
    - one pack JSON file (`gentle.cloning_patterns.v1`)
    - one single-template JSON file (`gentle.cloning_pattern_template.v1`)
    - one directory tree (recursive `*.json` import; files must use one of the
      schemas above)
  - imports are transactional; if one template fails validation, no imported
    template changes are kept
  - expanded scripts can execute `op ...` and `workflow ...` statements and
    optionally roll back via `--transactional`
  - `template-run` supports non-mutating preflight mode via `--validate-only`
  - template-run responses now include a preflight payload
    (`gentle.macro_template_preflight.v1`) with warnings/errors and typed
    input/output port validation rows (`contract_source` indicates whether
    checks came from template metadata or routine catalog)
  - preflight includes cross-port semantic checks (alias/collision checks,
    input sequence/container consistency, and sequence-anchor semantics when
    sequence context is unambiguous)
  - routine-family semantic checks are now supported:
    - Gibson routines validate adjacent fragment overlap compatibility against
      configured overlap length before execution
    - Restriction routines validate enzyme-name resolution, duplicate-enzyme
      misuse, enzyme-site presence across bound input sequences, and common
      digest parameter sanity (`left_fragment`/`right_fragment`,
      `extract_from`/`extract_to`)
  - mutating `macros run` / `macros template-run` executions always persist one
    lineage macro-instance record (`ok`/`failed`/`cancelled`)
  - successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in error messages
  - `macros instance-list` and `macros instance-show` expose persisted lineage
    macro-instance records as first-class introspection contracts

- `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT] [--seq-id SEQ_ID]`
  - shared-shell/CLI routine catalog discovery surface
  - default catalog path: `assets/cloning_routines.json`
  - typed catalog schema: `gentle.cloning_routines.v1`
  - response schema: `gentle.cloning_routines_list.v1`
  - filters are case-insensitive; query performs substring match across
    routine id/title/family/status/template/tags/summary plus explainability
    metadata fields
  - when `--seq-id` is supplied, planning estimates also consume the active
    construct-reasoning graph for that sequence so variant-derived assay
    suggestions can bias routine-family ranking deterministically
- `routines explain ROUTINE_ID [--catalog PATH] [--seq-id SEQ_ID]`
  - shared-shell/CLI routine explainability surface
  - response schema: `gentle.cloning_routine_explain.v1`
  - returns one routine definition plus normalized explanation payload
    (purpose/mechanism/requires/contraindications/disambiguation/failure modes)
    and resolved confusing alternatives
  - when `--seq-id` is supplied, the response also includes the
    sequence-aware `routine_preference_context` plus a planning estimate for
    that routine, so explain-stage inspection stays aligned with list/compare
    ranking for the same construct-reasoning graph
- `routines compare ROUTINE_A ROUTINE_B [--catalog PATH] [--seq-id SEQ_ID]`
  - shared-shell/CLI deterministic routine comparison surface
  - response schema: `gentle.cloning_routine_compare.v1`
  - returns both routine definitions plus comparison payload:
    shared/unique tags, cross-reference status, aligned difference-matrix rows,
    and merged disambiguation questions
  - includes planning-aware estimate rows in comparison payload:
    - `estimated_time_hours`
    - `estimated_cost`
    - `local_fit_score`
    - `composite_meta_score`
  - `--seq-id` applies the same active-sequence construct-reasoning context
    used by `routines list`, so compare-stage estimates stay aligned with
    variant-driven assay preferences

- Planning meta-layer contracts (shared shell/CLI, engine-owned):
  - profile schema: `gentle.planning_profile.v1`
  - objective schema: `gentle.planning_objective.v1`
  - estimate schema: `gentle.planning_estimate.v1`
  - suggestion schema: `gentle.planning_suggestion.v1`
  - sync-status schema: `gentle.planning_sync_status.v1`
  - cloning-consultation schema:
    `gentle.planning_cloning_consultation.v1`
  - protein-expression handoff schema:
    `gentle.protein_expression_handoff.v1`
  - merge precedence for effective profile:
    - `global_profile -> confirmed_agent_overlay -> project_override`
  - purchasing latency heuristic in v1:
    - each missing required material class adds default
      `procurement_business_days_default` (default `10`) to estimate
      (Monday-Friday business-day model; no holiday calendar yet)
    - business-day delays are converted to `estimated_time_hours` with a
      deterministic weekend-aware factor (`24h * 7/5` per business day)
  - schema compatibility rule:
    - profile/objective payloads with mismatched schema ids are rejected
      (`InvalidInput`) instead of silently coerced
- `planning consult cloning [--seq-id SEQ_ID] [--objective JSON_OR_@FILE] [--profile-scope effective] [--format json|text]`
  - read-only consultation route for cloning-strategy and vector-choice
    planning
  - consumes the effective planning profile, current or supplied planning
    objective, host-profile catalog, helper/vector catalog, and the existing
    routine-planning estimate logic
  - `PlanningObjective.biological_intent` is optional; V1 recognizes
    `protein_expression_max_yield` and phrase-like inputs such as "give me the
    maximal amount of protein" as a read-only high-yield protein-expression
    planning intent
  - ranks one best strategy candidate for each of the 11 catalogued routine
    families (`restriction`, `gibson`, `sequence`, `pcr`, `crispr`,
    `golden_gate`, `gateway`, `topo`, `ta_gc`, `infusion`,
    `nebuilder_hifi`)
  - emits ranked `strategy_candidates[]`, `vector_candidates[]`,
    `missing_questions[]`, `local_constraints[]`, and
    `suggested_next_actions[]`
  - for `protein_expression_max_yield`, the report keeps the same cloning
    strategy/vector ranking path but adds explicit missing questions for total
    versus soluble/active/purified/secreted yield, expression chassis,
    folding/cofactor requirements, toxicity/induction tolerance, and
    scale/purification endpoint
  - v1 deliberately ranks helper/vector candidates only from structured catalog
    fields and leaves marker, promoter/expression, and MCS/site constraints as
    explicit `missing_questions[]`
  - `--seq-id` is accepted for traceability, but v1 does not yet consume
    construct-candidate graphs when ranking strategies or helper vectors
- `planning protein-expression-handoff [--seq-id SEQ_ID] [--objective JSON_OR_@FILE] [--requirements JSON_OR_@FILE] [--profile-scope effective] [--format json|text]`
  - read-only handoff route for high-yield protein-expression requests
  - emits `gentle.protein_expression_handoff.v1`
  - records `biological_intent = protein_expression_max_yield` for phrase-like
    objectives such as "give me the maximal amount of protein"
  - fields include `product_definition`, `product_readiness`,
    optional `requirements`, `sequence_context`, `cds_assessment`, `tag_assessment`,
    `host_chassis_candidates[]`, `vector_route_candidates[]`,
    `missing_questions[]`, `service_handoff_candidates[]`, `warnings[]`, and
    `suggested_next_actions[]`
  - `product_readiness` reports whether the product context is usable for
    expression review, whether translation was possible, any blockers, and the
    human-review gate that still applies before construct or service action
  - `sequence_context` summarizes loaded `--seq-id` context without mutation:
    sequence name/type, nucleotide or protein length, feature count, GC percent
    and bounded GC range, ambiguous bases, and relevant CDS/protein/tag
    annotation rows
  - `cds_assessment` prefers annotated CDS features when present and otherwise
    labels whole-sequence ORF/CDS checks as fallback; it reports nucleotide
    length, inferred protein length, start/stop sanity, internal stops,
    ambiguous codons, translation table, and explicit warnings
  - `tag_assessment` summarizes annotated affinity/solubility/epitope/signal
    tag context and keeps tag preference, position, cleavage, and retention
    policy as explicit review inputs
  - `--requirements` accepts exactly
    `gentle.protein_expression_requirements.v1`; unknown fields and a missing or
    different schema are rejected rather than silently upgraded
  - the requirements record groups reviewed answers by biological topic:
    `yield_goal`, `chassis`, `localization`, `folding`,
    `toxicity_induction`, `tag_policy`, `scale_purification`, and
    `outsourcing`, with optional reviewer/provenance notes
  - topic presence means that topic was explicitly reviewed; this allows an
    empty `folding` list to mean "reviewed, no special requirement known"
    rather than "not answered"
  - validated reviewed topics remove only their matching
    `missing_questions[]` rows; partial records leave all unrelated questions
    visible and do not alter host/vector ranking
  - `outsourcing.allowed = false` withholds provider preflight and quote next
    actions while preserving in-house cloning and product-review actions; the
    retained provider provenance row is marked
    `withheld_by_outsourcing_requirement` and carries no shell command
  - the bundled
    `docs/examples/external_services/geneart_protein_expression_request.json`
    remains the template/provenance source for the first V1 service scaffold
  - without a provider-ready product, the service candidate is explicitly
    `draft_example_review_required` and no provider command appears in either
    text or structured suggested actions
  - a provider-ready selected DNA or protein sequence produces a
    `product_draft_review_required` request preview whose `source_target`
    references the project `seq_id` or `protein_seq_id`; its preflight and,
    where readiness permits, quote actions serialize that same request rather
    than operating on the tutorial protein
  - project-bound request commands do not inline sequence letters; exact CDS
    boundaries, protein target, expression requirements, and outsourcing
    permission remain review gates before external handoff
  - no provider network call, order, optimization, or construct mutation is
    made
  - if no usable CDS/protein context is found, `missing_questions[]` asks for
    coding sequence, ORF, CDS annotation, or target-protein boundaries rather
    than choosing an expression route
  - if a CDS/protein context is inferable, `missing_questions[]` shifts toward
    expression-specific unknowns: total vs soluble/active/purified/secreted
    yield, purification endpoint, tag preference, host/chassis,
    toxicity/induction, PTMs/cofactors, secretion/localization, scale, and
    delivery endpoint
  - `product_definition.readiness.status` also drives
    `suggested_next_actions[]`: `whole_sequence_cds_candidate` and
    `annotated_cds_review_required` suggest GeneArt preflight, cloning
    consultation, and review-gated quote-packet preparation;
    `protein_sequence_review_required` suggests reverse-translation or
    provider protein-target handoff review; `needs_cds_boundary` suggests
    inspecting or marking CDS/ORF boundaries before expression planning
  - the optional text report renders provider commands from those same
    structured actions, so it cannot advertise preflight or quote work that
    JSON readiness withheld
  - this remains a review-gated handoff, not an autonomous construct designer:
    GENtle does not codon-optimize, mutate sequences, create constructs, query
    live providers, submit orders, or promise wet-lab yield
  - GUI follow-up: a dedicated Synthetic Biology / protein-expression handoff
    inspector should display existing `gentle.protein_expression_handoff.v1`
    reports later; this protocol slice intentionally keeps GUI exposure to the
    shared shell/CLI contract rather than adding a new dashboard
- `planning profile show [--scope global|project_override|confirmed_agent_overlay|effective]`
  - inspect one planning profile scope or merged effective profile
- `planning profile set JSON_OR_@FILE [--scope global|project_override|confirmed_agent_overlay]`
  - set/replace selected planning profile scope
- `planning profile clear [--scope global|project_override|confirmed_agent_overlay]`
  - clear selected planning profile scope
- `planning objective show`
  - inspect current planning objective
  - objective may now also carry:
    - optional `helper_profile_id`
    - optional `preferred_routine_families[]`
  - when present, routine-ranking routes synthesize helper-aware
    `routine_preference_context` and apply a transparent family-alignment bonus
  - the same synthesized context is now also reused by routine-decision traces
    and engine-owned macro-template suggestions so planner-facing adapters do
    not need their own helper-specific heuristics
- `planning objective set JSON_OR_@FILE`
  - set/replace planning objective
- `planning objective clear`
  - clear planning objective (engine defaults apply)
- `planning suggestions list [--status pending|accepted|rejected]`
  - list pending/resolved planning sync suggestions
- `planning suggestions accept SUGGESTION_ID`
  - accept suggestion and apply patch into confirmed overlay/objective
- `planning suggestions reject SUGGESTION_ID [--reason TEXT]`
  - reject suggestion with optional reason
- `planning sync status`
  - inspect planning sync lifecycle metadata
- `planning sync pull JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - register inbound advisory suggestion as pending
- `planning sync push JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - register outbound advisory suggestion as pending
  - payload for `planning sync pull|push`:
    - optional `profile_patch` (`gentle.planning_profile.v1`)
    - optional `objective_patch` (`gentle.planning_objective.v1`)
    - optional `message`
  - activation policy remains explicit user action (`accept`/`reject`);
    no auto-apply in v1

- `screenshot-window OUTPUT.png`
  - currently disabled by security policy
  - returns deterministic disabled message from shared shell/CLI/GUI command
    paths
  - kept as a reserved artifact-export adapter contract; approval of the narrow
    GUI Agent-help attachment flow does not enable headless capture

Agent Assistant image attachment extension:

- GUI Agent help requests may add `x_attachments[]` to
  `gentle.agent_request.v1`; ordinary text-only requests omit the field.
- Each entry uses schema `gentle.agent_attachment.v1` and carries `id`,
  `kind=image`, `file_name`, `mime_type`, absolute temporary `path`,
  `byte_len`, `sha256`, optional `source_window_title`, `capture_backend`, and
  pixel dimensions.
- Attachments are accepted only after explicit GUI capture and preview, are
  limited to four images / 20 MiB each, and are rejected when the selected
  `AgentSystemSpec.supports_image_attachments` value is false.
- External JSON-stdio adapters receive the validated typed attachment record.
  The bundled Codex and Pi bridges translate it to their native image arguments;
  native OpenAI, Anthropic, Mistral, and OpenAI-compatible transports encode
  the image as provider-native multimodal content.
- Remote request inspection redacts the local temporary path. Persisted
  `AgentConversationTurn.attachments[]` stores only path-free provenance and
  never stores image bytes.
- This is a GUI-host capability, not a capture operation callable by an agent,
  CLI, MCP, JS, or Lua adapter. Those adapters can consume explicitly supplied
  image-capable agent requests in future, but cannot initiate screen capture.
- A validated `AgentResponse.screenshot_request` may ask the GUI to show a
  consent card. GENtle resolves targets only from its current open-window model,
  excludes Agent Assistant, rejects ambiguous/closed targets, and binds one
  ephemeral approval to the originating turn, system, project generation, and
  selected viewport. Provider/project/target changes, conversation clearing,
  another request, decline, timeout, and replay fail closed. Approval performs
  one egui viewport capture only; the unchanged `x_attachments[]` preview and
  `Ask Agent` boundary remain the sole transport path. Native ScreenCaptureKit
  capture remains available only through the user's direct context-menu action.

- `agents list [--catalog PATH]`
  - Lists configured agent systems from catalog JSON.
  - Default catalog: `assets/agent_systems.json`.

- `agents preflight SYSTEM_ID [--live] [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N]`
  - Returns read-only transport/runtime metadata as `gentle.agent_preflight.v1`.
  - Default behavior is config-only for CLI/MCP compatibility.
  - `--live` adds a non-generating setup probe for `native_openai`,
    `native_anthropic`, `native_mistral`, `native_openai_compat`, and
    `pi_local_stdio`.
  - `live_probe` fields:
    - `probe_kind`: `model_discovery | command_shape`
    - `enabled`
    - `attempted_endpoints`
    - `selected_endpoint`
    - `reachable`
    - `auth_ok`
    - `model_list_ok`
    - `selected_model_seen`
    - `status_class`:
      `ok | missing_key | auth_failed | quota_or_billing | model_missing | endpoint_unreachable | unsupported_transport | provider_error`
    - `message`
    - `provider_error_code`
  - Endpoint/local-command policy:
    - OpenAI: `GET /models` with bearer auth
    - Anthropic: `GET /models` with `x-api-key` and `anthropic-version`
    - Mistral: `GET /models` with bearer auth
    - OpenAI-compatible: `/models`, then `/v1/models` fallback when the base
      URL is not already `/v1`
    - Pi Local: invoke the selected Pi executable with GENtle's intended
      `--print --no-session --no-tools --no-extensions --no-skills
      --no-prompt-templates --no-context-files --no-approve --system-prompt ...`
      command shape plus any selected `--model`, then append `--help` so the
      local CLI parses the flags without sending a prompt to a model
      (`probe_kind = command_shape`; authentication/model-list booleans are not
      claimed by this probe)
    - no chat/completion/responses request is made
    - quota/billing is reported only when that provider error appears during
      the model-list probe

- `agents discover-models SYSTEM_ID [--catalog PATH] [--base-url URL]`
  - Returns discovered model ids as `gentle.agent_models.v1`.
  - Native HTTP systems use model-list endpoints. `codex_local_stdio` instead
    reads visible ids from the local Codex `models_cache.json`, excludes hidden
    entries, and performs no provider request or credential-file read.
  - `pi_local_stdio` invokes `pi --list-models`, returns provider-qualified ids,
    and reports a deterministic setup error when Pi has no authenticated model.

- `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
  - Invokes one configured agent system via catalog transport.
  - `--base-url` applies a per-request runtime base URL override for native
    transports (`native_openai`, `native_anthropic`, `native_mistral`,
    `native_openai_compat`).
  - `--model` applies a per-request runtime model override for native
    transports (`native_openai`, `native_anthropic`, `native_mistral`,
    `native_openai_compat`), `codex_local_stdio`, and `pi_local_stdio`. The
    stdio bridges map the shared `GENTLE_AGENT_MODEL` override to their
    respective `--model MODEL` argument.
  - `--timeout-secs` applies a per-request timeout override for stdio/native
    transports (maps to `GENTLE_AGENT_TIMEOUT_SECS`).
  - `--connect-timeout-secs` applies a per-request HTTP connect timeout override
    for native transports (maps to `GENTLE_AGENT_CONNECT_TIMEOUT_SECS`).
  - `--read-timeout-secs` applies a per-request read timeout override for
    stdio/native transports (maps to `GENTLE_AGENT_READ_TIMEOUT_SECS`).
  - `--max-retries` applies a per-request transient retry budget override
    (maps to `GENTLE_AGENT_MAX_RETRIES`; `0` disables retries).
  - `--max-response-bytes` applies a per-request response body/output cap
    override (maps to `GENTLE_AGENT_MAX_RESPONSE_BYTES`).
  - `--no-state-summary` suppresses project context injection.
  - Suggested-command execution is per-suggestion only (no global always-execute).

- `agents plan SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--max-candidates N] [--no-state-summary] [--no-mutating-candidates]`
  - Accepts prose but returns typed `gentle.agent_plan_result.v1` candidates.
  - `shell` candidates compile to shared shell commands.
  - `op`, `workflow`, and `ui_intent` are first-class plan payload kinds on the
    stored JSON boundary.

- `agents execute-plan PLAN_JSON_OR_@FILE --candidate-id ID [--confirm]`
  - Executes one stored planner candidate and returns
    `gentle.agent_execution_result.v1`.
  - Candidate execution never silently re-plans.
  - Nested `agents ask` / `agents plan` / `agents execute-plan` shell payloads
    are rejected.

- JavaScript and Lua wrappers expose the same planner boundary through thin
  shared-shell-backed helpers:
  - `plan_agent_system(...)`
  - `execute_agent_plan(...)`

Machine-facing planner schemas:

- `gentle.agent_plan_request.v1`
  - fields:
    - `schema`
    - `system_id`
    - `prompt`
    - optional `state_summary`
    - optional `max_candidates`
    - optional `allow_mutating_candidates`
- `gentle.agent_plan_result.v1`
  - fields:
    - `schema`
    - `assistant_message`
    - `questions[]`
    - `candidates[]`
  - candidate fields:
    - `candidate_id`
    - `title`
    - `rationale`
    - `kind = shell | op | workflow | ui_intent`
    - `mutating`
    - `requires_confirmation`
    - `execution_mode = ask | auto`
    - exactly one payload field:
      - `shell_command`
      - `operation`
      - `workflow`
      - `ui_intent`
- `gentle.agent_execution_result.v1`
  - fields:
    - `schema`
    - `candidate_id`
    - `title`
    - `kind`
    - `mutating`
    - `requires_confirmation`
    - `confirmed`
    - `state_changed`
    - `output`

Agent bridge catalog schema (`gentle.agent_systems.v1`):

```json
{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "openai_gpt5_stdio",
      "label": "OpenAI GPT-5 (stdio bridge)",
      "description": "Optional human-readable description",
      "transport": "external_json_stdio",
      "command": ["openai-agent-bridge", "--model", "gpt-5"],
      "env": {},
      "working_dir": null
    },
    {
      "id": "openai_gpt5_native",
      "label": "OpenAI GPT-5 (native HTTP)",
      "transport": "native_openai",
      "model": "gpt-5",
      "base_url": "https://api.openai.com/v1",
      "env": {}
    },
    {
      "id": "anthropic_claude_sonnet_native",
      "label": "Claude Sonnet (native Anthropic HTTP)",
      "transport": "native_anthropic",
      "model": "claude-sonnet-4-6",
      "base_url": "https://api.anthropic.com/v1",
      "env": {}
    },
    {
      "id": "mistral_large_native",
      "label": "Mistral Large (native Mistral HTTP)",
      "transport": "native_mistral",
      "model": "mistral-large-latest",
      "base_url": "https://api.mistral.ai/v1",
      "env": {}
    }
  ]
}
```

Transport notes:

- `builtin_echo`: offline/demo transport.
- `external_json_stdio`: requires local bridge executable from `command[0]`.
  The default catalog includes `codex_local_stdio`, which runs
  `scripts/codex-agent-bridge` as one such bridge. That bridge delegates
  authentication to the local Codex CLI/App login, runs Codex in a conservative
  empty working directory, and emits strict `gentle.agent_response.v1` JSON for
  the existing GENtle validator. It resolves the executable from `CODEX_BIN`,
  then `codex` on `PATH`, then known local installs including the current macOS
  ChatGPT bundle's `Contents/Resources/codex` CLI and `Contents/MacOS/ChatGPT`
  app executable, followed by the legacy standalone Codex app path. Its model
  selector reads the CLI's non-secret local model metadata cache; leaving
  `Codex default` selected omits an explicit `--model` argument.
  The default catalog also includes `pi_local_stdio`, which runs
  `scripts/pi-agent-bridge`. That bridge delegates authentication to Pi,
  launches one ephemeral print request in an empty temporary directory, and
  disables tools, sessions, project context files, extensions, skills, and
  prompt templates. It resolves Pi from `PI_BIN`, `PATH`, or common local
  install paths. `pi --list-models` supplies provider-qualified model ids;
  leaving `Pi default` selected omits an explicit model argument.
- `native_openai`: built-in OpenAI HTTP adapter; requires `OPENAI_API_KEY`
  (environment or system-level `env` override in catalog entry).
- `native_anthropic`: built-in Anthropic Claude HTTP adapter; requires
  `ANTHROPIC_API_KEY` (environment or system-level `env` override in catalog
  entry).
- `native_mistral`: built-in Mistral HTTP adapter; requires
  `MISTRAL_API_KEY` (environment or system-level `env` override in catalog
  entry).
- `native_openai_compat`: built-in OpenAI-compatible local HTTP adapter
  (`/chat/completions`), intended for local services such as Jan/Msty/Ollama
  when they expose an OpenAI-compatible endpoint. API key is optional.
  `msty_mlx_local_compat_template` points at `http://localhost:11973/v1` for
  Msty/MLX Knife model servers; the older Msty gateway template remains on
  `http://localhost:11964`.
- `GENTLE_AGENT_BASE_URL` (or CLI `--base-url`) overrides catalog `base_url`
  per request for `native_openai`, `native_anthropic`, `native_mistral`, and
  `native_openai_compat`.
- `GENTLE_AGENT_MODEL` (or CLI `--model`) overrides catalog `model` per request
  for `native_openai`, `native_anthropic`, `native_mistral`, and
  `native_openai_compat`, `codex_local_stdio`, and `pi_local_stdio`.
- `GENTLE_AGENT_TIMEOUT_SECS` (or CLI `--timeout-secs`) overrides request
  timeout per attempt for agent transports.
- `GENTLE_AGENT_CONNECT_TIMEOUT_SECS` (or CLI `--connect-timeout-secs`)
  overrides HTTP connect timeout for native transports.
- `GENTLE_AGENT_READ_TIMEOUT_SECS` (or CLI `--read-timeout-secs`) overrides
  read timeout for stdio/native transports.
- `GENTLE_AGENT_MAX_RETRIES` (or CLI `--max-retries`) overrides transient retry
  count (`0` disables retries).
- `GENTLE_AGENT_MAX_RESPONSE_BYTES` (or CLI `--max-response-bytes`) overrides
  response-size cap per attempt (stdout/stderr or HTTP body).
- `native_openai_compat` requires a concrete model name; value `unspecified`
  is treated as missing and the request is rejected until a model is provided.
- `native_openai_compat` does not silently switch host/port; it uses catalog
  `base_url` or explicit `GENTLE_AGENT_BASE_URL`.

Agent request payload schema (`gentle.agent_request.v1`):

```json
{
  "schema": "gentle.agent_request.v1",
  "system_id": "openai_gpt5_stdio",
  "prompt": "User request text",
  "sent_at_unix_ms": 1768860000000,
  "state_summary": {},
  "x_introspection": {
    "schema": "gentle.agent_introspection_context.v1",
    "fact_expression_schema": "gentle.fact_expression.v1",
    "fact_graph_schema": "gentle.project_fact_graph.v1",
    "projection_scope": "engine_project_graph_without_external_evidence",
    "fact_limit": 128,
    "total_fact_count": 2,
    "included_fact_count": 2,
    "omitted_fact_count": 0,
    "truncated": false,
    "fact_type_counts": {
      "sequence.exists": 1,
      "sequence.kind": 1
    },
    "facts": [
      {
        "fact": "sequence.exists",
        "domain": "project",
        "subject": { "kind": "sequence", "id": "tp73" }
      },
      {
        "fact": "sequence.kind",
        "domain": "project",
        "subject": { "kind": "sequence", "id": "tp73" },
        "value": "dna"
      }
    ],
    "retrieval_routes": [
      { "purpose": "fact definitions and current state grouped by domain", "command": "introspect facts [--domain project|view|host|config] [--seq-id SEQ_ID]" },
      { "purpose": "complete current project fact graph", "command": "facts graph" },
      { "purpose": "evaluate one fact expression", "command": "facts eval FACT_EXPR_JSON_OR_@FILE" },
      { "purpose": "check capability readiness", "command": "introspect readiness [CAPABILITY_ID] [--arg NAME=VALUE ...] [--seq-id SEQ_ID]" },
      { "purpose": "verify declared effects after execution", "command": "introspect verify-effects CAPABILITY_ID [--arg NAME=VALUE ...] [--seq-id SEQ_ID]" }
    ],
    "notes": [
      "Missing open-world facts are unknown, not false; absence requires explicit proof evidence."
    ]
  },
  "x_conversation": {
    "schema": "gentle.agent_conversation.v1",
    "turns": [
      {
        "user_message": "Retrieve TP73.",
        "response": {
          "schema": "gentle.agent_response.v1",
          "assistant_message": "Which species should I use?",
          "questions": ["Which species should I use?"],
          "suggested_commands": []
        },
        "system_id": "codex_local_stdio",
        "system_label": "Codex Local",
        "completed_at_unix_ms": 1768859999000
      }
    ]
  },
  "x_local_references": {
    "schema": "gentle.agent_local_reference_context.v1",
    "catalog_entry_count": 2,
    "installed_reference_count": 1,
    "included_reference_count": 1,
    "omitted_reference_count": 0,
    "truncated": false,
    "references": [
      {
        "genome_id": "Human GRCh38 Ensembl 116",
        "species": "homo_sapiens",
        "description": "Human GRCh38 Ensembl 116",
        "aliases": ["GRCh38"],
        "tags": ["human", "grch38"],
        "sequence_source_type": "remote",
        "annotation_source_type": "remote",
        "sequence_ready": true,
        "annotation_ready": true,
        "fasta_index_ready": true,
        "gene_index_ready": true,
        "transcript_index_ready": true,
        "gene_extraction_ready": true
      }
    ],
    "retrieval_routes": [
      { "purpose": "list/filter reference catalog entries without downloading", "command": "genomes list [--filter TEXT]" },
      { "purpose": "inspect one reference install and its reusable components", "command": "genomes status GENOME_ID" }
    ],
    "warnings": []
  },
  "x_gui_context": {
    "schema": "gentle.agent_gui_context.v1",
    "host_available": true,
    "recent_project_count": 1,
    "recent_projects": [
      {
        "item_id": "recent-91a8f94b0dc2778a",
        "display_label": "tp73_project.gentle.json (projects)",
        "file_name": "tp73_project.gentle.json",
        "parent_label": "projects",
        "list_position": 1,
        "exists": true,
        "byte_count": 45210,
        "modified_at_unix_ms": 1787740200000,
        "current_project": false,
        "open_command": "ui open recent-project recent-91a8f94b0dc2778a"
      }
    ],
    "tutorial_project_count": 1,
    "included_tutorial_project_count": 1,
    "omitted_tutorial_project_count": 0,
    "tutorial_projects_truncated": false,
    "tutorial_projects": [
      {
        "chapter_id": "simple_pcr_selection_gui",
        "decimal_id": "04.01",
        "display_label": "04.01 Simple PCR selection",
        "title": "Simple PCR selection",
        "summary": "Open a worked PCR-selection project.",
        "group": "Primers, PCR & qPCR",
        "tier": "core",
        "example_id": "simple_pcr_selection_gui",
        "online": false,
        "review_status": "reviewed",
        "review_stale": false,
        "open_command": "ui open tutorial-project simple_pcr_selection_gui"
      }
    ],
    "configuration_sections": [
      {
        "section_id": "agent-systems",
        "title": "Agent Systems",
        "detail": "Configure agent providers, endpoints, credentials, models, and runtime limits.",
        "open_command": "ui open configuration agent-systems"
      }
    ],
    "warnings": []
  },
  "x_local_documents": {
    "schema": "gentle.agent_local_documents.v1",
    "max_document_count": 4,
    "max_document_bytes": 131072,
    "max_total_bytes": 262144,
    "linked_markdown_depth": 1,
    "total_included_byte_count": 39,
    "documents": [
      {
        "requested_path": "/checkout/docs/roadmap.md",
        "resolved_path": "/checkout/docs/roadmap.md",
        "source": "explicit_prompt",
        "media_type": "text/markdown",
        "original_byte_count": 39,
        "included_byte_count": 39,
        "truncated": false,
        "sha256": "sha256:7ca1345cb504834b64563bf72aa1a325756388e9e321d04d6bd997b4e4e41ff3",
        "content": "# Roadmap\nRun the GUI smoke checklist.\n"
      }
    ],
    "warnings": []
  }
}
```

`x_introspection` is an optional, backward-compatible request extension. GENtle
adds it together with `state_summary` when project-context injection is
enabled. It is generated from the same engine snapshot and contains at most 128
concrete facts. `fact_type_counts` is computed over the complete projection;
`total_fact_count`, `included_fact_count`, `omitted_fact_count`, and
`truncated` make bounded output explicit. The abbreviated example above shows
only the non-zero count entries relevant to its two facts; current GENtle
requests include every registered fact type, including zero counts.

`projection_scope = engine_project_graph_without_external_evidence` means the
snapshot does not claim request-specific evidence or GUI-host availability.
Those can be projected explicitly through the listed `introspect facts` route.

The extension is grounding, not authority. A model must not treat an omitted
fact or a missing open-world fact as false. When the bounded projection is
insufficient, `retrieval_routes` identifies the deterministic read-only command
to propose. Suggested-command readiness remains advisory, while the shared
parser and engine retain hard execution validation and confirmation policy.
Providers do not need filesystem access to use this context: the same payload
is sent to native HTTP transports and external stdio adapters such as Codex
Local.

`x_conversation` is an optional, backward-compatible extension. GENtle owns
the transcript because transports such as Codex Local may start a fresh,
ephemeral model process for every request. The current `prompt` is the newest
user message; `x_conversation.turns` contains successful earlier exchanges in
chronological order. GENtle stores at most 50 turns with the project and sends
at most the 12 most recent turns on a request. Provider credentials and API
keys are not part of this object. An adapter should reuse explicit facts from
the transcript (for example a species or sequence id) unless the current
prompt changes them, rather than asking for the same value again.
Prompt and response text are stored verbatim; provider credential fields are
excluded, but callers must not place secrets inside the prompt itself.

`x_local_references` is included on every request, independently of optional
project-state injection. It is a bounded inventory of catalog entries with an
existing prepared manifest, not a list of everything that could be downloaded.
The context contains stable catalog identity and component readiness but no
local paths. Building it performs no network request and no arbitrary filesystem
scan. Agents should prefer a compatible `gene_extraction_ready` reference and
compose `genomes extract-gene` followed by strand-aware `genomes extend-anchor`
commands, then use `ui open sequence-window` when the user asked to see the
result. They must not claim that catalog-only entries are installed. Direct
Ensembl retrieval remains a confirmation-gated fallback when no compatible
prepared reference is present.

`x_gui_context` is an optional, backward-compatible extension attached by the
live GUI Agent Assistant even when project-state injection is disabled. It
mirrors the current recent-project menu, generated executable tutorial catalog,
and global Configuration sections. Recent rows expose a deterministic opaque
`item_id`, filename and parent label, 1-based menu position, existence,
size/time metadata, current-project status, and an exact `open_command`.
Absolute recent-project paths are not sent. Only the live GUI host can resolve
`ui open recent-project ITEM_ID`, and it rechecks both the current list and file
existence before opening. Stale or unknown ids fail closed.

Tutorial rows carry chapter/example identity, title/summary, group/tier,
online status, review metadata, and exact chapter-opening commands. Total,
included, omitted, and truncation fields prevent a provider from interpreting
a bounded list as complete. Configuration rows use the shared five-section
vocabulary. Their commands only navigate to the existing tab; credentials,
executable paths, and other global values continue through GENtle's visible
Apply/Cancel model. Headless/CLI requests omit `x_gui_context` because they have
no live private GUI host catalog, while the same commands remain parser-valid
intent records for CLI and MCP.

`x_local_documents` is an optional, backward-compatible extension generated
when the current prompt explicitly contains an absolute path to a supported
UTF-8 text/document file (`.md`, `.txt`, `.json`, `.toml`, `.yaml`, `.csv`,
`.tsv`, and related text forms). GENtle copies bounded content into the request;
the provider never receives a filesystem handle. Limits are four documents,
128 KiB per document, and 256 KiB total. For an explicitly named Markdown
document, GENtle may include one level of prioritized local guide/runbook links
that resolve beneath the same directory tree. Each included document records
its requested/resolved path, source, media type, byte counts, truncation flag,
and SHA-256 digest. Read failures are retained in `warnings` so an agent does
not ask the user to paste a document that was actually included or silently
pretend that an unreadable document was reviewed.

This feature is an explicit disclosure boundary: the copied content is sent to
the selected local or cloud provider. It does not grant Pi, Codex, or an HTTP
provider ambient filesystem access, and document text cannot override GENtle's
response schema, parser validation, confirmation gates, or nested-agent ban.
Sequence/GenBank/FASTA files are not attached through this route; use GENtle's
normal import commands and project state for biological records.

Documentation context for agent systems:

- GENtle agent systems should be given the command documentation bundle before
  they propose commands: `docs/glossary.json`, `docs/cli.md` (especially the
  operand metavariable conventions), `docs/protocol.md`,
  `docs/ai_prompt_contract.md`, and, for biology/domain grounding,
  `docs/ai_cloning_primer.md`, `docs/ai_task_playbooks.md`,
  `docs/examples/ai_cloning_examples.md`, and optionally
  `docs/ai_glossary_extensions.json`.
- The glossary is a command index and parser-facing syntax contract; it is not
  a complete natural-language ontology. Placeholder operands such as `QUERY`,
  `ID`, `SEQ_ID`, `GENOME_ID`, `ENTRY_ID`, `PATH`, and `OUTPUT.svg` must be
  interpreted using the CLI manual and route-specific documentation.
- If an inner helper has not been given the relevant docs, or if a placeholder
  remains ambiguous after reading them, it should ask a clarifying question
  rather than invent a command, identifier, species alias, coordinate system, or
  local filesystem path.

Agent response payload schema (`gentle.agent_response.v1`):

```json
{
  "schema": "gentle.agent_response.v1",
  "assistant_message": "Text response",
  "questions": ["Optional follow-up question"],
  "suggested_commands": [
    {
      "title": "Optional short label",
      "preconditions": ["Optional requirement before the command can work"],
      "precondition_expr": {
        "all": [
          { "fact": "sequence.exists", "id": "demo_seq" }
        ]
      },
      "expected_outcomes": ["Optional observable effect expected if the command succeeds"],
      "expected_effects": [
        { "fact": "report.exists", "kind": "restriction_scan", "subject": "demo_seq" }
      ],
      "rationale": "Optional reason",
      "command": "state-summary",
      "execution": "ask"
    }
  ],
  "screenshot_request": {
    "id": "inspect-visible-map-1",
    "reason": "Inspect the visible feature lanes and disabled controls."
  }
}
```

`screenshot_request` is optional and additive. It may contain only `id` and
`reason`; the id is non-empty, ASCII-stable, and at most 128 characters, while
the trimmed human-readable reason is non-empty and at most 512 characters. One
object is the v1 maximum, and it counts as response content. A request does not
name a target and does not authorize or initiate capture. Arrays, paths,
coordinates, native window ids, capture commands, approval fields, unknown
fields, and excessive values fail response validation.

Native HTTP transports parse stochastic LLM text after provider extraction. If
that text is already a JSON object with `assistant_message`, `questions`, or
`suggested_commands`, or `screenshot_request`, GENtle may repair a
missing/non-string `schema` field to
`gentle.agent_response.v1` before validation. Native HTTP transports also
unwrap a single top-level Markdown `json` code fence before validation, because
local chat models often add that wrapper despite a JSON-only instruction.
`external_json_stdio` adapters receive the same narrowly bounded fence
unwrapping, but otherwise remain strict: they must emit the `schema` string
themselves, and prose around a fence or JSON object is rejected.

Agent command-scope declaration:

- The in-app Agent Assistant currently accepts only GENtle shared-shell
  commands that correspond to GENtle GUI/shared-engine capabilities.
- Local slash aliases are a small GENtle-owned compatibility layer, described
  by `gentle_protocol::shell_alias_registry()`. They are not operating-system
  shell commands and are not an OpenClaw gateway.
- In the GUI Agent Assistant prompt box, one-line local control commands are
  intercepted before contacting the model: `/...` slash commands plus bare
  `help`, `state-summary`, and `capabilities`. Mistyped slash commands fail
  locally instead of being reinterpreted as natural-language prompts.
- Bare absolute paths such as `/path/to/file` are not Ollama-style attachments
  in GENtle. Use `/open file PATH`, `/import file PATH`, or the GUI import path
  for local sequence files.
- Accepted read-only aliases:
  - `/help [TOPIC]`
  - `/list`
  - `/history` (canonical command: `history status`)
  - GUI presentation is host-specific but result-preserving: `/help` opens the
    built-in Shell Commands help tab, while `/list` renders the returned
    `EngineStateSummary` as a compact project overview. The latest structured
    local-command output remains copyable as JSON. Object-row context menus
    instantiate commands that are parsed and executed through the same shared
    shell; they do not introduce GUI-only engine operations. The main project
    overview remains the complete table/graph for lineage analyses and other
    project objects. Headless adapters continue to receive the same unchanged
    shell result payloads.
- Accepted session-history aliases:
  - `/undo` (canonical command: `history undo`)
  - `/redo` (canonical command: `history redo`)
  - Undo and redo are session-local state transitions. Agent suggestions for
    either transition require an explicit user action; GENtle rejects their
    automatic execution even when a provider labels the suggestion `auto`.
  - In the GUI, these aliases use the same guarded history handlers as
    `Edit -> Undo` / `Redo`: they are disabled while background jobs are active
    and refresh dependent windows and presentation caches after success.
- Accepted GUI/file aliases:
  - `/open` and `/import` request the GUI sequence-file picker. In headless
    shell execution they return a structured UI-intent payload explaining that
    a GUI host is required.
  - `ui open sequence-window SEQ_ID`,
    `ui focus sequence-window SEQ_ID`, and
    `ui close sequence-window SEQ_ID` control the viewer for an already-loaded
    sequence record. These are non-mutating GUI intents and do not refetch,
    reimport, or delete sequence data.
  - `/open sequence-window SEQ_ID` is a GENtle-local slash alias for
    `ui open sequence-window SEQ_ID`.
  - `ui close TARGET` requests closing a catalogued GUI tool/dialog target
    without mutating project data.
  - `/close sequence-window SEQ_ID` requests closing an open DNA sequence
    viewer for `SEQ_ID`. It is a non-mutating GUI intent and does not remove
    the loaded sequence record from project state.
  - `ui selection sequence-window SEQ_ID --range START..END` requests setting
    the DNA sequence-window selection in 0-based, end-exclusive coordinates.
    Headless shell execution returns `gentle.ui_sequence_selection_intent.v1`;
    the GUI Agent Assistant applies the selection when the loaded sequence is
    available.
  - `ui selection sequence-window SEQ_ID` requests inspecting the current GUI
    selection for that sequence window.
  - `display show TARGET`, `display hide TARGET`, and
    `display visibility TARGET on|off` update shared engine display state for
    feature/display targets such as `features`, `gene-features`,
    `cds-features`, `mrna-features`, `repeat-features`, `array-features`,
    `tfbs`, `restriction-enzymes`, `gc-contents`, `open-reading-frames`, and
    `methylation-sites`.
  - `/open file PATH [--id ID]` and `/import file PATH [--id ID]` load an exact
    user-provided sequence file through `LoadFile`.
- Accepted explicit sequence/fetch aliases:
  - `/paste sequence --sequence-text DNA [--id ID]`
  - `/features restriction-scan SEQ_ID [--enzyme NAME]`
    - `SEQ_ID` is an already-loaded GENtle sequence id; optional `--enzyme`
      narrows the restriction-site scan to one enzyme and may be repeated by
      using the full shared-shell `features restriction-scan` form.
  - `/fetch genbank ACCESSION [--id ID]`
    - `ACCESSION` is an NCBI GenBank nucleotide accession; `--id` names the
      local GENtle sequence id.
  - `/fetch ncbi ACCESSION [--id ID]` (synonym for GenBank accession fetch)
  - `/fetch uniprot QUERY [--id ID]`
    - `QUERY` is a UniProtKB/Swiss-Prot accession or entry name; `--id` names
      the local GENtle metadata entry.
  - `/fetch ensembl QUERY [--species NAME] [--assembly NAME] [--flank-bp N|--flank-5p-bp N --flank-3p-bp N] [--id ID] [--no-open]`
  - `/fetch ensembl-gene QUERY [--species NAME] [--assembly NAME] [--flank-bp N|--flank-5p-bp N --flank-3p-bp N] [--id ID] [--no-open]`
    - `QUERY` is an approved gene symbol for the requested species or a stable
      Ensembl gene id; for human symbols this means HGNC-approved symbols such
      as `FUS` or `TP53`, while HGNC IDs such as `HGNC:11998` need resolution
      before database-specific fetching.
    - `--assembly` fails closed if Ensembl resolves the gene on another or an
      unreported assembly. Flanks are measured from the strand-aware gene ends;
      `--flank-bp` sets both directions equally.
    - In the GUI Agent Assistant, `--no-open` suppresses the automatic DNA
      sequence viewer open after the fetched gene is imported into project
      state; headless shell fetches are unaffected.
  - `/fetch ensembl-protein QUERY [--id ID]`
    - `QUERY` is an Ensembl protein/translation stable id or accepted Ensembl
      lookup query, not a UniProt accession.
  - `/fetch ensembl-region SPECIES CHR START END [--strand +|-] [--id ID]`
    - `SPECIES` is an Ensembl species name such as `homo_sapiens`; coordinates
      are assembly positions; `--id` names the local GENtle sequence id.
  - `/fetch dbsnp RS_ID GENOME_ID [--id ID]`
    - `RS_ID` is a dbSNP identifier and `GENOME_ID` is a prepared GENtle genome
      catalog id.
- External fetch aliases are tagged as `CapabilityMutation::External` and
  should require explicit confirmation/network opt-in on agent surfaces.
- Unknown slash commands return a typed
  `gentle.shell_alias_rejection.v1` diagnostic with
  `supported_alternatives[]`.
- Slash commands such as `/grep`, `/find`, `/ls`, `/new`, and `/example` are
  intentionally rejected unless GENtle explicitly implements them later.
- OpenClaw-like filesystem, operating-system, and gateway commands are not part
  of this contract yet, though a future gateway layer may add them.
- When local file discovery is needed and the exact path is unknown, the agent
  should ask the user to locate the file by regular operating-system means
  rather than suggesting invented commands such as `fs.find` or
  `gentle.load_sequence`.
- On macOS, it is appropriate to suggest Finder search or Spotlight for that
  outside-GENtle discovery step.

Agent execution intent semantics:

- `chat`: explain/ask only, never executed as shell command.
- `ask`: executable suggestion requiring explicit user confirmation.
- `auto`: executable suggestion eligible for automatic execution only when
  caller enables `--allow-auto-exec`.
- Suggested-command `title` is the short user-facing intent shown in GUI
  tables. `preconditions[]` lists state requirements such as "sequence
  `demo_seq` exists"; it is advisory text, while the command parser still owns
  hard validity checks. `precondition_expr` is an optional machine-readable
  fact-graph expression. The GUI evaluates it advisory-only against the current
  project graph and displays readiness as `ready`, `blocked`, or `unknown`.
  `expected_outcomes[]` lists postcondition-like effects expected if the
  command succeeds, such as "a restriction-site report is available"; these are
  user-visible expectations to verify, not guarantees. `expected_effects[]` is
  an optional list of machine-readable fact-graph effects expected after
  success.
- Negative logic should be expressed as proof-backed positive facts when
  possible. For example, "no EcoRI site exists in `demo_seq`" should become a
  fact such as `restriction_site.absent` with an enzyme, subject, range, and
  basis report; it should not be inferred merely from the absence of a
  `restriction_site.present` fact.

Project fact graph:

GENtle's `state-summary` remains the compact project snapshot for humans,
adapters, and LLM context. The engine also exposes
`gentle.project_fact_graph.v1`, a deterministic typed projection for planning
and readiness checks. The current slice covers loaded sequence facts,
persisted reverse-translation, protein-derivation, primer/qPCR,
restriction-cloning handoff, sequencing-confirmation, CUT&RUN read, and
RNA-read interpretation report facts, persisted sequencing-trace evidence
records, curated isoform-panel bindings, persisted exon-skip selection plans,
plus explicit
restriction-site scan evidence:

- `sequence.exists`, `sequence.kind`, `sequence.length`, `sequence.circular`
  are closed-world facts over the loaded project state. `sequence.kind` uses
  normalized values `dna`, `rna`, or `protein`.
- `report.exists`, `restriction_site.present`, and `restriction_site.absent`
  are open-world facts and require a proof basis. Persisted primer-pair design
  reports project as `report.exists` with `value: "primer_design"`, persisted
  reverse-translation reports project with `value: "reverse_translation"`,
  persisted protein-derivation reports project with
  `value: "protein_derivation"`,
  persisted qPCR design reports project with `value: "qpcr_design"`,
  persisted restriction-cloning PCR handoff reports project with
  `value: "restriction_cloning_pcr_handoff"`, persisted
  sequencing-confirmation reports project with
  `value: "sequencing_confirmation"`, persisted CUT&RUN read reports project
  with `value: "cutrun_read"`, and persisted RNA-read interpretation reports
  project with `value: "rna_read"`. Each carries a hard-fact basis.
  Restriction-scan reports supplied through `--evidence` project as
  `value: "restriction_scan"`.
  Read-only report inspection commands use the same values for readiness, so a
  planner can prove whether a specific persisted report id is inspectable before
  issuing `show-*` commands.
- `sequencing_trace.exists` is a closed-world fact over imported sequencing
  trace evidence records stored in project metadata. `seq-trace list` and
  `ListSequencingTraces` are catalog-ready; `seq-trace show` and
  `ShowSequencingTrace` require `sequencing_trace.exists(TRACE_ID)`;
  `seq-trace import` and `ImportSequencingTrace` can verify the same fact when
  the caller supplied an explicit `TRACE_ID`.
- `isoform_panel.exists` and `isoform_panel.seq_id` are closed-world facts over
  imported curated isoform panels stored in project metadata.
  `isoform_panel.exists(PANEL_ID)` proves that the panel id is present, while
  `isoform_panel.seq_id(PANEL_ID) == SEQ_ID` proves that the panel was imported
  for the requested sequence context. `panels inspect-isoform`,
  `panels render-isoform-svg`, and `RenderIsoformArchitectureSvg` require both
  facts, so agents do not accidentally inspect an identically named panel bound
  to a different sequence.
- `exon_skip_plan.exists` is a closed-world fact over persisted exon-skip
  selection plans stored in project metadata. `transcripts exon-skip-plan`,
  `exon_skip_plan`, and `PlanExonSkippedIsoform` can verify this fact when the
  caller supplies a deterministic `PLAN_ID`; `transcripts
  exon-skip-materialize`, `exon_skip_materialize`, and
  `MaterializeExonSkippedIsoform` require it before consuming the plan.
- `introspect readiness` evaluates fact-annotated descriptors through their
  full `precondition_expr`, including `any` branches. This allows shared raw
  operation rows such as `ExportPrimerDesignReport` to express that either a
  `primer_design` or a `qpcr_design` report can satisfy readiness for the same
  operation payload.
- Raw primer/qPCR design operation rows `DesignPrimerPairs`,
  `DesignInsertionPrimerPairs`, and `DesignQpcrAssays` are fact-annotated with
  the same engine-owned semantics as their shell routes: the template sequence
  must exist, and a deterministic `REPORT_ID` can be verified as the generated
  `report.exists` effect.
- Raw cDNA assay-test operation rows `TestCdnaPcr` and `TestCdnaQpcr` are
  fact-annotated over `sequence.exists(SEQ_ID)` and declare conservative
  `may_on_success` effects because JSON/SVG artifacts, materialized assay
  products, and product-gel outputs depend on optional payload fields. Shell
  routes `primers test-cdna-pcr` and `primers test-cdna-qpcr` expose the same
  readiness contract. `TestCdnaQpcrFasta` and
  `primers test-cdna-qpcr-fasta` have no project-state precondition because
  their transcript catalogs are external FASTA/FASTA.gz paths validated during
  execution.
- `ProjectGenomeInterval` is fact-annotated as a no-project-precondition
  coordinate-projection operation over an external interval-map file; it returns
  a projection report but does not declare durable project facts.
- `ReverseTranslateProteinSequence` is fact-annotated with the same
  engine-owned semantics as `reverse-translate run`: the input sequence must
  exist and be protein-kind, and a deterministic `OUTPUT_ID` can be verified as
  the generated coding-DNA sequence.
- Raw transcript/protein/splicing derivation rows `DeriveTranscriptSequences`,
  `DeriveProteinSequences`, and `DeriveSplicingReferences` require
  `sequence.exists(SEQ_ID)` and declare conservative `may_on_success`
  sequence-creation effects because output ids are feature/prefix/uniqueness
  derived. `DeriveProteinSequences` also verifies
  `report.exists(REPORT_ID) == protein_derivation` when a deterministic
  report id is supplied.
- Protease catalog and digest routes are fact-annotated where existing facts
  can express readiness. `proteases list` and `proteases show` are
  no-precondition catalog reads. `proteases digest`,
  `ProteaseDigestProteinSequence`, and `proteases digest-gel-svg` require the
  input sequence to exist and have `sequence.kind == "protein"`. Digest
  materialization declares only a `may_on_success` sequence effect because
  peptide ids are prefix/index-derived. `RenderProteaseDigestGelSvg` can be
  ready from either an explicit protein sequence or a persisted
  `protein_derivation` report, and SVG render operations model output paths as
  `artifact.written` external handoffs. `RenderProteinGelSvg`,
  `RenderProteinGelReportsSvg`, and `RenderProtein2dGelSvg` require
  `report.exists(REPORT_ID) == protein_derivation`; grouped report rendering
  should repeat readiness checks for every report id in the payload.
- `MaterializeVariantAllele` is fact-annotated with the same engine-owned
  semantics as `variant materialize-allele`: the input sequence must exist, and
  a deterministic `OUTPUT_ID` can be verified as the generated allele sequence.
- `AlignSequences` is fact-annotated with the same readiness semantics as
  `align compute`: both query and target sequence ids must exist, and the
  operation declares no project-state effects because the alignment is returned
  in the operation result.
- Dotplot and flexibility analysis payloads project closed-world project facts:
  `dotplot.exists(DOTPLOT_ID)` and `flexibility_track.exists(TRACK_ID)`.
  `dotplot compute`, `ComputeDotplot`, `dotplot overlay-compute`,
  `ComputeDotplotOverlay`, `flex compute`, and `ComputeFlexibilityTrack`
  require loaded sequence inputs and can verify those payload facts when
  deterministic ids are supplied. `dotplot show`, `flex show`, and
  `RenderDotplotSvg` use the stored payload facts for readiness;
  `RenderDotplotSvg` models the SVG path as an `artifact.written` external
  handoff.
- Candidate-window optimization routes project a closed-world
  `candidate_set.exists(SET_NAME)` fact from persisted candidate-set metadata.
  Candidate generation requires `sequence.exists(SEQ_ID)` and verifies the
  named set. Show/metrics/score/delete routes require the set fact;
  filter/top-k/pareto/set-op routes require input set facts and verify the
  named output set. Delete rows verify a hard `not(candidate_set.exists)` effect
  after successful deletion, using the same closed-world fact evaluator as
  readiness.
- Guide-design routes project closed-world `guide_set.exists(GUIDE_SET_ID)`,
  `guide_filter_report.exists(GUIDE_SET_ID)`, and
  `guide_oligo_set.exists(OLIGO_SET_ID)` facts from persisted guide metadata.
  Upsert verifies the guide-set fact; practical filtering verifies the filter
  report and, when an output id is supplied, the filtered guide set; oligo
  generation verifies the oligo-set fact. Guide deletes verify a hard
  `not(guide_set.exists)` effect after successful deletion. Guide
  CSV/FASTA/protocol exports use `artifact.written` external handoff effects
  rather than local project-state mutations.
- Wet-lab container authoring routes project closed-world
  `container.exists(CONTAINER_ID)`, `arrangement.exists(ARRANGEMENT_ID)`, and
  `rack.exists(RACK_ID)` facts from persisted project state. Container
  exclusivity updates require and preserve the container fact. Serial
  arrangement creation can verify a deterministic arrangement id when supplied;
  rack creation can verify a deterministic rack id when supplied. Rack
  move/profile/template/block commands require and preserve `rack.exists`, while
  rack SVG/OpenSCAD/simulation exports model output paths as `artifact.written`
  external handoffs. Single-container transform operations (`DigestContainer`,
  `LigationContainer`, and `FilterContainerByMolecularWeight`) require
  `container.exists(CONTAINER_ID)` and declare `may_on_success` product effects
  because concrete product ids are derived from execution parameters rather than
  a single deterministic output binding.
- List-valued pool/container rows (`MergeContainers`, `MergeContainersById`,
  `Ligation`, `FilterByMolecularWeight`, `FilterByDesignConstraints`,
  `ExportPool`, `ExportPoolCollection`, `RenderPoolGelSvg`, and their
  shell/adapter aliases) use
  descriptor-side `foreach_arg` atoms for bound readiness. Sequence-list rows
  expand each supplied `INPUTS` id into `sequence.exists(id)` checks,
  container-list rows expand each supplied `CONTAINER_IDS` id into
  `container.exists(id)` checks, and mixed gel-render rows accept any one valid
  source mode: sequence ids, container ids, or `ARRANGEMENT_ID`. Artifact-
  producing pool rows still declare `artifact.written(OUTPUT_PATH)` external
  handoff effects.
- `screenshot-window` is represented in the catalog as a reserved host/view
  capability. It requires `ui.host_available == true` and a
  `host.tool_available(SCREENSHOT_CAPTURE) == true` fact, and remains blocked
  in normal builds/policies where screenshot capture is disabled.
- Macro-template routes project closed-world
  `workflow_macro_template.exists(TEMPLATE_NAME)`,
  `candidate_macro_template.exists(TEMPLATE_NAME)`, and
  `macro_instance.exists(MACRO_INSTANCE_ID)` facts from persisted project
  metadata and lineage state. Candidate/workflow template upsert verifies the
  corresponding template fact. Template show/delete/run routes require the
  corresponding template fact; template deletes verify a hard
  `not(..._macro_template.exists)` effect after successful deletion. Workflow
  template runs may record a macro-instance lineage row, and
  `macros instance-show` requires the `macro_instance.exists` fact.
- Protein/gene metadata routes project closed-world
  `uniprot_entry.exists(ENTRY_ID)`,
  `ensembl_gene_entry.exists(ENTRY_ID)`, and
  `ensembl_protein_entry.exists(ENTRY_ID)` facts from stored UniProt and
  Ensembl metadata entries. Fetch/import routes can verify explicit entry ids;
  show routes and metadata-backed sequence imports require the matching stored
  metadata fact. Metadata-backed sequence imports can also verify
  `sequence.exists(OUTPUT_ID)` when a deterministic output id was supplied.
- UniProt projection routes project closed-world
  `uniprot_projection.exists(PROJECTION_ID)` facts. `uniprot map` and
  `ProjectUniprotToGenome` require both `uniprot_entry.exists(ENTRY_ID)` and
  `sequence.exists(SEQ_ID)` and can verify a supplied deterministic
  `PROJECTION_ID`. `uniprot projection-list` is ready without project
  preconditions. Projection-specific show, feature-coding, Ensembl-link
  resolution, transcript-accounting, Ensembl-comparison, and audit-generation
  routes require the projection fact. Persisted UniProt projection audit and
  audit-parity reports still project as `report.exists(REPORT_ID)` facts, so
  their show/export routes can be readiness-checked.
- Reference/helper genome extraction routes are fact-aware for the project
  sequence they create, not for prepared-cache existence. `genomes|helpers
  extract-region`, `extract-gene`, and `extract-promoter`, the raw
  `ExtractGenome*` operation rows, and JS/Lua aliases can verify
  `sequence.exists(OUTPUT_ID)` when a deterministic output id was supplied.
  Catalog lookup, prepared-cache compatibility, annotation filters, and repeat
  materialization remain execution-time checks. Anchor extension and
  verification routes require `sequence.exists(SEQ_ID)`.
- Raw persisted-report operation rows mirror the shell report readiness model
  where the report kind is unambiguous. `ListSequencingConfirmationReports`,
  `ListCutRunReadReports`, and `ListRnaReadReports` are catalog-ready with no
  project preconditions. `ShowSequencingConfirmationReport`,
  `ShowCutRunReadReport`, and `ShowRnaReadReport` require the matching
  `report.exists` fact and declare no effects.
  `ExportSequencingConfirmationReport`, `ExportSequencingConfirmationSupportTsv`,
  `ExportCutRunReadCoverage`, and `ExportRnaReadReport` require the matching
  `report.exists` fact and model the output path as an `artifact.written`
  external handoff.
- `artifact.written` is a host-domain open-world fact for external file or
  handoff artifacts. It is used as an `external_handoff` effect for commands
  such as `render-svg`, `render-rna-svg`, `render-lineage-svg`, and
  `features export-bed`; `render-feature-expert-svg` follows the same SVG
  handoff pattern. The raw engine operation rows `ExportFeaturesBed`,
  `ExportSequenceContextBundle`, `RenderSequenceSvg`, `RenderRnaStructureSvg`,
  `RenderFeatureExpertSvg`, `RenderTfbsScoreTrackCorrelationSvg`, and
  `RenderLineageSvg` expose the same readiness and external-handoff model for
  registry-driven adapters.
  `InspectSequenceContextView` exposes the same sequence-readiness model
  without an artifact effect. Project readiness for those commands is still
  evaluated from their project inputs.
- `config.param` is a config-domain closed-world fact for engine-owned
  behavioral parameters, persisted display settings, metadata-backed BLAST
  option parameters, and accepted `set-param` compatibility aliases. The subject
  kind is `other`, the subject id is the canonical parameter name or accepted
  alias, and the fact value is the current JSON value. Alias facts intentionally
  mirror their canonical stored value so readiness/effect checks can bind the
  spelling used by a GUI, shell, MCP, or agent route.
- Read-only sequence inspectors such as `rna-info`, `features query`, and
  `features tfbs-summary` use `sequence.exists` readiness and do not project
  persistent report effects. The raw `SummarizeTfbsRegion` and
  `QueryProteinResidueGenomicCoordinates` operation rows expose the same
  sequence-readiness model. Variant promoter/reporter inspectors
  (`variant promoter-context`, `variant reporter-fragments`,
  `SummarizeVariantPromoterContext`, and
  `SuggestPromoterReporterFragments`) use the same sequence-readiness model
  and model optional JSON paths as `artifact.written` external handoffs rather
  than persistent project reports. `variant annotate-promoters`,
  `AnnotatePromoterWindows`, and `AnnotateTfbs` are sequence-gated mutating
  annotation rows without hard fact effects until feature freshness/count facts
  exist. `inspect-feature-expert` currently uses the same sequence-level
  readiness; target-specific facts can narrow this later.
- The known fact vocabulary is registered in engine protocol code and is also
  appended to the Agent Assistant system prompt, so prompt grounding and
  deterministic evaluation share one list of fact names.
- Unknown fact names evaluate to `unknown` rather than failing the payload.
- Boolean expressions use three-valued Kleene logic: `not(unknown)` remains
  `unknown`.
- Absence is never inferred from a missing `restriction_site.present` fact.
  Use a proof-backed `restriction_site.absent` fact from a covering
  zero-hit restriction scan.
- Open-world proof facts are built from persisted report metadata or explicit
  evidence bundles today. Restriction-site scan reports are not auto-persisted
  into project state, and a future evidence ledger such as
  `state.metadata["fact_evidence"]` should be added deliberately rather than as
  an implicit side effect.

Shared-shell routes:

- `ui selection sequence-window SEQ_ID [--range START..END|--start N --end N]`
  - Returns `gentle.ui_sequence_selection_intent.v1` in headless shell
    contexts; GUI hosts may apply it to open or pending DNA sequence windows.
- `display show TARGET`, `display hide TARGET`,
  `display visibility TARGET on|off`
  - Apply `SetDisplayVisibility` and return
    `gentle.display_visibility_result.v1`.
  - In introspection, `display` has no project preconditions and declares a
    `view.visible_tracks` `view_session` effect. The current display-layer
    booleans are also projected as a closed-world `view.visible_tracks` fact
    on the `host` UI subject so words-only clients can read the visible layer
    state back deterministically. The fact vocabulary also includes
    `view.selection` for `ui selection` effects; selection effects remain
    descriptors rather than verified hard postconditions.
- `SetLinearViewport { start_bp, span_bp }`
  - In introspection, the raw operation row declares no preconditions and
    verifies the closed-world `view.viewport` fact for the
    `linear_sequence` UI subject. The effect value is an object with
    `start_bp` and `span_bp` fields.
- `ui intents`
  - Returns `gentle.ui_intents.v1`.
  - In introspection, this is a no-precondition, fact-annotated catalog route
    for discovering GUI intent targets and command forms.
- `facts graph [--evidence SCAN.json ...]`
  - Returns `gentle.project_fact_graph.v1`.
  - `--evidence` accepts JSON emitted by `features restriction-scan ... --path`.
- `facts eval FACT_EXPR_JSON_OR_@FILE [--evidence SCAN.json ...]`
  - Returns `gentle.fact_evaluation.v1` with `truth`,
    `unmet_atoms[]`, and `unknown_atoms[]`.
  - Readiness mapping is advisory: `satisfied -> ready`,
    `unsatisfied -> blocked`, `unknown -> unknown`.
- `introspect facts [--domain project|view|host|config] [--seq-id SEQ_ID] [--evidence SCAN.json ...] [--ui-host true|false]`
  - Returns `gentle.introspection.v1` route `facts`.
  - Groups fact read-back by domain. Existing project facts keep their names
    (`sequence.*`, `report.*`, `restriction_site.*`) and carry additive
    `domain: "project"`.
  - `--seq-id` filters projected facts to a concrete sequence subject.
  - Always projects `ui.host_available` in the `view` domain. Headless shell
    contexts default to `false`; GUI-attached callers may pass
    `--ui-host true` when evaluating view-intent readiness.
  - Projects deterministic `host.tool_available` rows for configured agent
    systems from the same non-live availability checks used by
    `agents list`/`agents preflight`.
- `introspect capabilities [--kind operation|view_intent|host_config]`
  - Returns `gentle.introspection.v1` route `capabilities`.
  - Projects the shared protocol capability registry so glossary commands,
    engine operations, and MCP tools are discoverable through the same
    words-only route.
  - For glossary-backed command rows, each row's `registry` object carries the
    descriptor-backed `usage`, `interfaces`, and `aliases` fields alongside the
    lower-level input/output schema and adapter surfacing metadata.
  - Rows with validated fact preconditions/effects carry
    `annotation_status: "fact_annotated"`; registry-only rows remain visible as
    `annotation_status: "registry_only"` without invented readiness semantics.
    The payload labels this as
    `annotation_scope: "registry_with_fact_annotated_slice"`.
  - The top-level `introspect facts`, `introspect capabilities`,
    `introspect readiness`, `introspect verify-effects`, `introspect runtime`,
    and `introspect all`
    shell routes are also fact-annotated as read-only self-description/status
    routes with no project-state preconditions.
- `introspect runtime`
  - Returns `gentle.runtime_status.v1`, with process-local live `frames[]` and
    best-effort observed `activities[]` from existing persisted/project
    ledgers.
  - This route does not write a runtime-status file and it is not a historical
    job log; active frames disappear when the scoped work guard exits.
  - Shell command dispatch pushes a `shell_command` frame while the command is
    running. Genome preparation also pushes a nested `resource_prepare` frame
    and mirrors its current phase/progress fields.
  - GUI background jobs such as reference preparation, tutorial-project
    opening, track import, BLAST, dbSNP fetch, Agent Assistant requests, model
    discovery, and ClawBio runs own live frames while their task records are
    active.
  - Observed `activities[]` reuse existing ledgers only: genome prepare
    `.prepare_activity.json` markers from default reference/helper catalogs,
    CUT&RUN shared-asset activity markers from default discovery, and
    project-backed BLAST async jobs. Each activity is tagged with its source,
    scope, lifecycle status, and observation class such as `live`,
    `cross_process`, `completed`, `failed`, `cancelled`, `stale`, or
    `unknown`.
  - MCP exposes the same snapshot as the `runtime_status` tool for the MCP
    server process and accepts optional `state_path` so project-backed async
    registries can be inspected. A standalone `gentle_mcp` process cannot see
    unsaved GUI process-local frames unless a future explicit cross-process
    bridge is added; the snapshot says this in `observability_notes[]`.
  - On Unix, `gentle`, `gentle_cli`, and `gentle_mcp` install a SIGUSR1 waiter
    thread at startup. Sending SIGUSR1 prints the process-local runtime
    snapshot to STDERR; it does not create or update any file.
- `introspect readiness [CAPABILITY_ID] [--arg NAME=VALUE ...] [--seq-id SEQ_ID] [--readiness ready|blocked|unknown] [--evidence SCAN.json ...] [--ui-host true|false]`
  - Returns `gentle.introspection.v1` route `readiness`.
  - `--seq-id` is shorthand for `--arg SEQ_ID=...` and asks about
    sequence-subject readiness for fact-annotated descriptors.
  - `--readiness` filters returned rows after evaluation.
  - Without args, templated atoms such as `{"subject":{"arg":"SEQ_ID"}}`
    report `unknown` with reason `unbound argument`.
  - Self-description/status commands such as `help`, `capabilities`,
    `state-summary`, and `history status` have empty precondition expressions
    and are ready in catalog mode.
  - Fact-layer commands `facts graph` and `facts eval` are fact-annotated,
    no-mutation routes; `facts eval` validates its supplied expression at
    execution time rather than as a projected project-state precondition.
  - `set-param` and the lower-case adapter alias `set_parameter` are
    fact-annotated as `host_config` routes with no project preconditions. Their
    verified effect binds `PARAM_NAME` to a `config.param` subject and parses
    `PARAM_VALUE` as JSON before comparing it with the projected parameter
    value.
  - Direct sequence-derivation operation rows such as `Reverse`, `Complement`,
    `ReverseComplement`, `Branch`, and `ExtractRegion` are fact-annotated over
    raw operation payload fields. Their bound readiness requires
    `sequence.exists(INPUT_SEQ_ID)`, and hard-effect verification checks
    `sequence.exists(OUTPUT_ID)` when the caller supplies the deterministic
    output id that execution created.
  - Project/session utility rows are fact-annotated without existing project
    facts. `history undo`, `history redo`, `load-project`, and `load_project`
    declare `may_on_success` effects because concrete fact changes depend on
    the undo/redo stack or loaded project file. `save-project` and
    `save_project` model external project-state output as
    `artifact.written(OUTPUT_PATH)`. `load_dna` mirrors `LoadFile` for JS/Lua
    adapters and can verify `sequence.exists(OUTPUT_ID)` only when a stable
    adapter output id is bound.
  - External sequence creation routes are fact-annotated for shell and raw
    operation callers: `LoadFile`, `load_dna`, `genbank fetch`,
    `FetchGenBankAccession`, `ensembl-region fetch`, `FetchEnsemblRegion`,
    `dbsnp fetch`, `FetchDbSnpRegion`, `FetchUniprotLinkedGenBank`,
    `ImportUniprotEntrySequence`, `ensembl-gene import-sequence`,
    `ImportEnsemblGeneSequence`, `ensembl-protein import-sequence`, and
    `ImportEnsemblProteinSequence`. These routes have no project-state
    preconditions; effect verification checks `sequence.exists(OUTPUT_ID)` only
    when the caller binds the deterministic `as_id`/`output_id` used during
    execution.
  - Raw core sequence operation rows are fact-annotated where their state
    contracts are discrete. `SaveFile` requires `sequence.exists(SEQ_ID)` and
    models the output path as an `artifact.written` external handoff. `Digest`
    and the glossary `digest` alias require
    `sequence.exists(INPUT_SEQ_ID)` but declare no hard effect because fragment
    ids are prefix/index-derived. `ExtractAnchoredRegion` requires the same
    input sequence fact but declares no hard effect because candidate ids are
    output-prefix/rank-derived. `SelectCandidate` requires
    `sequence.exists(INPUT_SEQ_ID)` and verifies `sequence.exists(OUTPUT_ID)`
    when execution used a deterministic selected-candidate id. `Pcr`,
    `PcrAdvanced`, and `PcrMutagenesis` require
    `sequence.exists(TEMPLATE_SEQ_ID)` and verify `sequence.exists(OUTPUT_ID)`
    when execution used a deterministic product id.
    `PcrOverlapExtensionMutagenesis` currently models template readiness only
    because candidate ids are prefix/rank-derived.
  - The `render-dotplot-svg` shell alias is fact-annotated with the same
    sequence/dotplot readiness contract and external SVG handoff effect as the
    raw `RenderDotplotSvg` operation row.
  - `SetTopology` is fact-annotated over raw operation payload fields. Its
    bound readiness requires `sequence.exists(SEQ_ID)`, and hard-effect
    verification checks the closed-world `sequence.circular(SEQ_ID)` fact
    against the JSON boolean `CIRCULAR` binding.
  - `RecomputeFeatures` is fact-annotated as a readiness-only operation over
    raw operation payload fields. Its bound readiness requires
    `sequence.exists(SEQ_ID)`, but it declares no hard fact effect until
    computed-feature freshness is represented in the project fact graph.
  - `SetLinearViewport` is fact-annotated over raw operation payload fields.
    It has no project preconditions and verifies the closed-world
    `view.viewport` value for the `linear_sequence` UI subject using nested
    argument bindings for `START_BP` and `SPAN_BP`.
  - Persisted-report export commands such as `reverse-translate export-report`,
    `primers export-report`, `primers export-qpcr-report`, and
    `primers export-restriction-cloning-handoff`, plus
    `seq-confirm export-report`, `seq-confirm export-support-tsv`,
    `cutrun export-coverage`, `rna-reads export-report`, and specialized
    RNA-read artifact exports (`rna-reads export-hits-fasta`,
    `rna-reads export-target-quality`, `rna-reads export-paths-tsv`,
    `rna-reads export-abundance-tsv`, `rna-reads export-score-density-svg`,
    `rna-reads export-alignments-tsv`,
    `rna-reads export-isoform-triage-tsv`, and
    `rna-reads export-alignment-dotplot-svg`) use `report.exists` readiness and
    model the output path as an `artifact.written` `external_handoff` effect.
  - Sequencing-confirmation report commands project
    `report.exists == sequencing_confirmation`; `seq-confirm list-reports` is
    catalog-ready with no project preconditions. `seq-confirm run` and the raw
    `ConfirmConstructReads` row require `sequence.exists(EXPECTED_SEQ_ID)` and
    verify `report.exists(REPORT_ID) == sequencing_confirmation` when an
    explicit report id was supplied. Optional repeated read sequences and trace
    ids are exposed as readable evidence inputs; list-valued readiness binding
    remains a later introspection refinement.
  - Sequencing-primer overlay commands `seq-primer suggest` and
    `SuggestSequencingPrimers` require `sequence.exists(EXPECTED_SEQ_ID)`;
    optional primer sequence ids and confirmation reports are readable inputs
    rather than hard preconditions.
  - Sequencing-trace evidence commands project
    `sequencing_trace.exists(TRACE_ID)`; trace listing is catalog-ready, trace
    show requires the trace fact, and explicit-id trace import can verify the
    same fact after execution.
  - Report-list routes such as `cutrun list-read-reports` and
    `rna-reads list-reports` are fact-annotated as no-precondition catalog
    reads. `cutrun show-read-report`, `cutrun export-coverage`,
    `rna-reads show-report`, `rna-reads export-report`, and the specialized
    RNA-read artifact exports use the projected `cutrun_read` / `rna_read`
    report facts for bound readiness.
    `rna-reads align-report` / `AlignRnaReadReport` also use the projected
    `rna_read` report fact and verify that the report remains present after the
    mutating alignment pass. `rna-reads materialize-hits` /
    `MaterializeRnaReadHitSequences` use `rna_read` report readiness but do not
    yet declare a hard sequence-creation effect because created ids depend on
    the selected hits. `rna-reads preflight-isoforms` /
    `PreflightRnaReadIsoforms` use `sequence.exists(SEQ_ID)` readiness and do
    not mutate project state.
    `rna-reads show-alignment`, `rna-reads show-alignments`,
    `rna-reads summarize-gene-support`, `rna-reads inspect-gene-support`,
    `rna-reads inspect-alignments`, and `rna-reads inspect-concatemers` use the
    same `rna_read` report fact for bound readiness; the shell routes with an
    optional output path model it as an `artifact.written` external handoff.
    `SummarizeRnaReadGeneSupport` and `InspectRnaReadGeneSupport` use the same
    `rna_read` report fact and model their optional JSON path as an external
    artifact handoff.
  - Primer helper readbacks are fact-annotated where existing facts can express
    readiness: `primers preflight` is ready without project state,
    `primers seed-from-feature`, `primers seed-from-splicing`, and
    `primers restriction-cloning-vector-suggestions` require
    `sequence.exists(SEQ_ID)`, and
    `primers seed-restriction-cloning-handoff` requires both
    `report.exists(PRIMER_REPORT_ID) == primer_design` and
    `sequence.exists(DESTINATION_VECTOR_SEQ_ID)`. cDNA assay test routes remain
    registry-only until their conditional product-materialization effects are
    modeled.
  - No-project local catalog/report operations such as
    `SummarizeJasparEntries`, `BenchmarkJasparRegistry`, `ListJasparCatalog`,
    `ResolveTfQueries`, `ListReporterCatalog`, and `RecommendReporters` are
    fact-annotated as catalog-ready and model optional JSON `path` outputs as
    external artifact handoffs.
  - Reporter shell routes `reporters list` and `reporters recommend` mirror
    the raw reporter catalog operations as catalog-ready optional JSON artifact
    handoffs. `reporters export-corpus` and raw `ExportReporterCorpus` are
    also catalog-ready and model their required JSON/JSONL output path as an
    `artifact.written` external handoff.
  - Service status/provider catalog routes (`services status`,
    `services providers list`, `services providers doctor`,
    `services delivery-route`, `services project-preflight`,
    `services project-quote`, `services handoff`, and `services guide`) are
    fact-annotated as ready without project state. Doctor, quote, and handoff
    output paths are modeled as external artifact handoffs. These routes
    prepare local classification, validation, guide, or handoff records only;
    they do not submit vendor orders. `services route-project-source` remains
    separate because its preconditions depend on the selected project object
    kind.
  - Planning read-only routes (`planning consult cloning`,
    `planning protein-expression-handoff`, `planning profile show`,
    `planning objective show`, `planning suggestions list`, and
    `planning sync status`) are fact-annotated as ready without project state
    and declare no side effects. Planning profile/objective set/clear,
    suggestion accept/reject, and sync pull/push are also ready without project
    facts, but declare only a non-verifiable `may_on_success` planning-state
    effect because profile/objective/suggestion facts are not projected yet.
  - Shell-level resource/catalog inspection routes such as
    `resources summarize-jaspar`, `resources status`,
    `resources sync-rebase`, `resources sync-jaspar`,
    `resources sync-ucsc-rmsk`, `resources import-gene-list-cache`,
    `resources import-ontology-assignment-cache`,
    `resources import-co-regulated-cache`, `resources install-ucsc-rmsk`,
    `resources prepare-ucsc-rmsk-index`,
    `resources sync-jaspar-remote-metadata`, `resources benchmark-jaspar`,
    `resources suggest-ucsc-rmsk-index`, `resources list-jaspar`,
    `resources inspect-jaspar`, `resources resolve-tf-query`,
    `resources list-publication-datasets`,
    `resources status-publication-dataset`, `genomes validate-catalog`,
    `helpers validate-catalog`, `helpers vocabulary list`, and
    `helpers vocabulary doctor` are fact-annotated as ready without project
    state. Optional JSON output paths are external artifact handoffs; catalog,
    motif, dataset, and file path validation remains an execution-time concern.
  - Local cache/resource mutation routes (`cache clear` and
    `resources prepare-publication-dataset`) are ready without project facts,
    but declare only a non-verifiable `may_on_success` effect because concrete
    filesystem/cache changes depend on supplied paths, download choices, and
    execution-time validation.
  - Host/helper/protease/microRNA catalog helper routes such as `hosts list`,
    `list_host_profile_catalog_entries`, `host_profile_catalog_entries`,
    `list_helper_catalog_entries`, `helper_catalog_entries`,
    `helper_semantics_vocabulary`, `helper_interpretation`, `proteases list`,
    `proteases show`, `mirna explain-seed`, and `mirna catalog-show` are
    fact-annotated as ready without project state. Protease list/show JSON
    output paths are modeled as optional external artifact handoffs; catalog
    lookup validation remains an execution-time concern.
  - Reference/helper genome cache inspection routes (`genomes status`,
    `helpers status`, `genomes genes`, `helpers genes`) and JS/Lua helper
    adapters (`list_reference_genomes`, `list_reference_catalog_entries`,
    `is_reference_genome_prepared`, `list_reference_genome_genes`) are
    fact-annotated as ready without project state. This is catalog-mode
    readiness only: concrete genome ids, prepared caches, catalog paths, regex
    filters, and biotype filters remain execution-time validation concerns
    until dedicated reference-resource facts are projected.
  - Built-in ladder catalog routes (`ladders list`, `inspect_dna_ladders`,
    `inspect_rna_ladders`, `list_dna_ladders`, `list_rna_ladders`) are
    fact-annotated as catalog-ready without project state. Ladder export routes
    (`ladders export`, `export_dna_ladders`, `export_rna_ladders`,
    `ExportDnaLadders`, `ExportRnaLadders`) additionally model their JSON
    output path as an `artifact.written` external handoff.
  - Agent-system catalog routes (`agents list`, `agent_systems`,
    `list_agent_systems`) are fact-annotated as catalog-ready without project
    state because they only enumerate configured systems.
  - `agents preflight`, `agents discover-models`, `agent_preflight`, and
    `agent_models` are fact-annotated host-config routes that require
    `host.tool_available(SYSTEM_ID) == true`. Preflight may emit an
    external-handoff report; model discovery declares no project effects.
    Asking and planning (`agents ask`, `ask_agent_system`, `agents plan`, and
    `agent_plan`) use the same host-tool availability precondition. Plan
    execution (`agents execute-plan`, `agent_execute_plan`) is payload-ready
    without project facts because readiness depends on the supplied plan JSON
    and candidate id; concrete project effects are command-dependent and are
    advertised only as `may_on_success`.
  - Read-acquisition routes (`reads acquire status`, `reads acquire prepare`,
    `reads acquire inspect`, `reads acquire cancel`, and raw `ReadAcquire*`
    rows) are fact-annotated as payload/path-ready without project facts.
    Prepare and cancel declare `artifact.written(WORK_DIR)` as an external
    handoff because they update external acquisition state. Introspection
    readiness does not validate manifest contents, SRA accessions, tools, or
    filesystem paths; those checks remain execution-time behavior.
  - Adapter parity aliases `state_summary`, `reference_catalog_entries`,
    `ui_intents`, `ui_prepared_genomes`, and `ui_latest_prepared`, plus shell
    routes `ui prepared-genomes` and `ui latest-prepared`, are fact-annotated
    with the same no-project readiness as their shared contracts
    (`state-summary`, reference catalog listing, and deterministic UI
    catalog/prepared-genome query routes).
  - Generic GUI intent requests (`ui open`, `ui focus`, `ui close`, and
    `ui_intent`) are fact-annotated as view intents that require
    `ui.host_available == true`. A headless adapter may still receive the
    structured intent payload, but readiness remains blocked until an attached
    GUI host can apply it.
  - Protocol-cartoon catalog/render/template rows are fact-annotated as ready
    without project state. Render/export rows model SVG or JSON outputs as
    `artifact.written` external handoffs; template input path validation remains
    an execution-time file concern rather than a project fact.
  - Prepared-cache inspection, CUT&RUN dataset catalog/status inspection, and
    array helper inspection rows are ready without project state. External
    file/directory/catalog validation remains an execution-time concern;
    `arrays render-probe-region-output-svg` additionally models its SVG output
    as an `artifact.written` external handoff.
  - Genome-track import and array projection rows (`tracks import-bed`,
    `tracks import-bigwig`, `tracks import-vcf`, their raw/MCP operation rows,
    `genomes|helpers blast-track`, `arrays project-microarray-track`,
    `ProjectMicroarrayTrack`, and `arrays project-probe-region-output`)
    require `sequence.exists` for the loaded target sequence. They currently
    remain readiness-only because feature freshness/track-update facts are not
    projected yet.
  - Tracked genome signal-file subscription rows
    (`tracks tracked list|add|remove|clear|apply`) have no project-fact
    preconditions. List is read-only; add/remove/clear/apply declare only a
    non-verifiable `may_on_success` project-state effect because imported
    features depend on current anchors and external file validation.
  - Sequence-scan reporting/rendering rows (`FindRestrictionSites`,
    `features tfbs-score-tracks-svg`, `RenderTfbsScoreTracksSvg`,
    `SummarizeTfbsScoreTracks`, `features tfbs-track-similarity`, and
    `SummarizeTfbsTrackSimilarity`) require `sequence.exists(SEQ_ID)` when the
    scan target is a loaded project sequence. Inline sequence text and
    motif/enzyme validation remain execution-time checks; JSON/SVG output paths
    are modeled as `artifact.written` external handoffs where present.
  - Catalog/list routes for candidate sets, guide sets, workflow macros,
    candidate macro templates, and routine catalogs are ready without project
    state. Routes that inspect a named persisted set/template remain
    registry-only until those object-existence facts are projected.
  - Construct-reasoning graph list routes (`construct-reasoning list-graphs`,
    `construct_reasoning_graphs`) are ready without project state; optional
    `SEQ_ID` is a filter, not a readiness precondition. Named graph
    inspection/action routes require
    `construct_reasoning_graph.exists(GRAPH_ID)`. Status/writeback actions
    verify that the graph remains present, and graph export models the JSON
    output as an `artifact.written` external handoff.
  - Persisted analysis-payload list routes (`dotplot list`, `flex list`) are
    ready without project state; optional `SEQ_ID` is a filter, not a readiness
    precondition. Show/render-by-id routes remain registry-only until
    dotplot/flex-track existence facts are projected.
  - Local metadata/catalog routes (`genomes list`, `helpers list`,
    `ensembl-gene list`, `ensembl-protein list`, and
    `gene-groups list|show|resolve|doctor`) are ready without project state.
    Catalog/group/token validation remains an execution-time concern, and
    gene-group JSON outputs are modeled as optional `artifact.written` external
    handoffs.
  - Ensembl discovery/catalog-maintenance routes (`genomes ensembl-available`,
    `helpers ensembl-available`, `ensembl_installable_genomes`,
    `list_ensembl_installable_genomes`, `genomes preview-ensembl-specs`,
    `helpers preview-ensembl-specs`, `genomes update-ensembl-specs`, and
    `helpers update-ensembl-specs`) are ready without project state. Update
    rows model the external catalog output as
    `artifact.written(OUTPUT_CATALOG_PATH)`, while concrete catalog paths and
    Ensembl/network availability remain execution-time concerns.
  - With args, the engine instantiates the atoms and evaluates them through the
    same fact evaluator used by `facts eval`.
- `introspect verify-effects CAPABILITY_ID [--arg NAME=VALUE ...] [--seq-id SEQ_ID] [--evidence SCAN.json ...] [--ui-host true|false]`
  - Returns `gentle.introspection.v1` route `verify-effects`.
  - Evaluates only effects with `effect_kind: "must_on_success"` for the named
    fact-annotated capability.
  - `verified=true` means all hard postconditions are satisfied against the
    current fact graph plus supplied evidence. `view_session`,
    `external_handoff`, and `may_on_success` effects remain descriptor metadata
    and are not asserted as local postconditions.
- `introspect all [--evidence SCAN.json ...] [--ui-host true|false]`
  - Returns facts, registry-backed capability descriptors, and catalog-mode
    readiness in one aggregate `gentle.introspection.v1` payload.

Initial fact-expression shape:

```json
{
  "all": [
    { "fact": "sequence.exists", "id": "demo_seq" },
    { "fact": "sequence.kind", "id": "demo_seq", "equals": "dna" }
  ]
}
```

Initial absence-proof effect shape:

```json
{
  "fact": "restriction_site.absent",
  "subject": { "kind": "sequence", "id": "demo_seq" },
  "enzyme": "EcoRI",
  "range": { "kind": "whole_sequence" },
  "basis": {
    "report_id": "restriction_scan_report_id",
    "report_kind": "restriction_scan",
    "evidence_class": "hard_fact"
  }
}
```

Agent schema/compatibility policy:

- `schema` is mandatory for catalog/request/response JSON objects.
- Supported major versions (current): `gentle.agent_systems.v1`,
  `gentle.agent_request.v1`, `gentle.agent_response.v1`.
- Future incompatible major versions (for example `.v2`) are rejected with a
  deterministic schema-unsupported error.
- Response validation is strict for canonical fields:
  - top-level allowed: `schema`, `assistant_message`, `questions`,
    `suggested_commands` plus extension keys prefixed with `x_` or `x-`
  - `suggested_commands[]` allowed: `title`, `preconditions`,
    `precondition_expr`, `expected_outcomes`, `expected_effects`, `rationale`,
    `command`, `execution` plus extension keys prefixed with `x_` or `x-`
  - unsupported canonical fields (for example `commands`, `mode`) are rejected

Execution safety rules:

- There is no global always-execute mode.
- Execution is per suggestion:
  - explicit run (`--execute-index`, `--execute-all`, GUI row `Run`)
  - optional auto-run only for `execution = auto` + `--allow-auto-exec`
- Recursive `agents ask`, `agents plan`, and `agents execute-plan` execution
  from suggested commands is blocked.

Failure-handling policy for external adapters:

- Adapter invocations use bounded retry with exponential backoff for transient
  failures.
- OpenAI `429` with `insufficient_quota` is treated as non-transient (no retry)
  and returned with the original API error body plus billing/usage guidance.
- Missing/unreachable adapter binaries fail gracefully with deterministic
  adapter-unavailable errors.
- CLI/shell errors are stable and prefixed for scripting, e.g.:
  - `AGENT_INVALID_INPUT`
  - `AGENT_SCHEMA_VALIDATION`
  - `AGENT_SCHEMA_UNSUPPORTED`
  - `AGENT_ADAPTER_UNAVAILABLE`
  - `AGENT_ADAPTER_TRANSIENT`
  - `AGENT_ADAPTER_FAILED`
  - `AGENT_RESPONSE_PARSE`
  - `AGENT_RESPONSE_VALIDATION`

ClawBio/OpenClaw integration scaffold schemas:

- integration path:
  `integrations/clawbio/skills/gentle-cloning/`
- included helper launchers:
  - `gentle_local_checkout_cli.sh` for local editable GENtle checkouts
  - `gentle_apptainer_cli.sh` for Apptainer/Singularity-backed `:cli` images
- wrapper request schema: `gentle.clawbio_skill_request.v1`
  - `mode`: `skill-info|capabilities|state-summary|shell|op|workflow|construct-reasoning-list-inspections|construct-reasoning-run-inspection|gene-protein-2d-gel|exon-skip-plan|exon-skip-materialize|agent-plan|agent-execute-plan|raw`
  - optional: `state_path`, `timeout_secs`
  - optional: `expected_artifacts[]`
    - wrapper-declared output files to copy into the ClawBio output bundle
      after command execution
    - relative paths resolve from the actual execution working directory
    - declared `.svg` paths remain the internal render/provenance source, but
      the wrapper now also rasterizes them into messenger-ready PNG artifacts
      in the output bundle at fixed deterministic scale `2.0`
    - text-bearing SVG rasterization requires at least one font face visible
      to `resvg`; `svg-png` reports `font_face_count` and fails early if an
      SVG contains text but no fonts are available. Headless deployments can
      install system fonts or set `GENTLE_SVG_FONT_FILE` /
      `GENTLE_SVG_FONT_DIR`.
    - `svg-pdf` uses the same rasterization path and embeds the rendered image
      into a single-page PDF for documentation/handoff bundles.
  - mode-specific:
    - `skill-info`: reports ClawBio skill/catalog metadata without invoking
      `gentle_cli`
    - `shell`: `shell_line`
    - `op`: `operation` (JSON object/string)
    - `workflow`: `workflow` or `workflow_path`
      - relative `workflow_path` resolves via current working directory, then
        `GENTLE_REPO_ROOT`, then the local GENtle repo containing the scaffold
        when discoverable
    - `construct-reasoning-list-inspections`: `graph_id`, with optional
      `fact_id`, `annotation_id`, `candidate_id`, `evidence_id`, `seq_id`,
      `action_kind`, and `summary_id`; the wrapper builds
      `construct-reasoning list-inspection-actions ...` and returns the same
      structured list payload as CLI shell callers
    - `construct-reasoning-run-inspection`: `graph_id` plus `action_id`, with
      optional dotplot tuning/render fields `word_size`, `step_bp`,
      `max_mismatches`, `tile_bp`, `dotplot_id`, and `render_svg_path`; the
      wrapper builds `construct-reasoning run-inspection-action ...`
    - `gene-protein-2d-gel`: `gene_symbol`, optional `species` (default
      `homo_sapiens`), optional `source` (currently `ensembl`), and optional
      `ladders[]`; the wrapper builds a deterministic workflow that fetches
      the expanded Ensembl gene record, imports transcript/exon/CDS features,
      derives protein-coding mRNA products, renders a 2D pI-vs-kDa SVG, and
      promotes the declared SVG into the PNG-first artifact bundle
    - `exon-skip-plan`: `seq_id`, `transcript_feature_id`, optional
      `skip_candidate_ids[]`, `skip_intervals_1based[]`,
      `overlap_intervals_1based[]`, `length_mod3_values[]`,
      `coding_mod3_values[]`, `coding_contexts[]`,
      `cds_phase_entry_kinds[]`, `feature_query`, and `plan_id`; stored
      candidate rows report transcript position, exon support frequency,
      flanking intron lengths, full-exon and coding-only modulo-3 frame cues,
      optional CDS boundary phases, stable `cds_phase_entry_kind`, and CDS
      phase warnings when a CDS-overlapping whole-exon skip changes coding
      length modulo 3
    - `exon-skip-materialize`: `plan_id`, optional `candidate_ids[]`,
      optional `output_prefix`, required `confirm=true`, and optional
      `return_items[]` (`genbank`, `cdna_fasta`, `amino_acid_sequence`,
      `amino_acid_fasta`) so ClawBio can declare exactly what it wants sent
      back after the reviewed plan is executed
    - `agent-plan`: `system_id`, `prompt`, optional planner/runtime overrides
      such as `catalog_path`, `base_url`, `model`, `max_candidates`,
      `include_state_summary`, and `allow_mutating_candidates`
    - `agent-execute-plan`: `plan` or `plan_path`, `candidate_id`, optional
      `confirm`
    - `raw`: `raw_args[]`
- wrapper result schema: `gentle.clawbio_skill_result.v1`
  - `status`: `ok|command_failed|timeout|failed|degraded_demo|incomplete|verification_failed`
  - includes resolver details, executed command, exit code, stdout/stderr, and
    generated artifact paths
  - `stdout_json` is populated when the wrapped `gentle_cli` stdout parses as
    JSON
  - `chat_summary_lines[]` is populated when `stdout_json.schema` is
    `gentle.sequence_context_view.v1`, so ClawBio/OpenClaw can relay the
    compact sequence-context summary before attaching larger SVG/BED artifacts
  - when no domain-specific summary is available but a command completed with
    parseable output, the wrapper synthesizes a short execution summary so
    chat renderers still show the command and output shape instead of only
    report section headings
  - `artifacts.collected[]` may enumerate declared output files copied into the
    wrapper bundle with `declared_path`, `bundle_path`, `source_path`,
    `copied_path`, and optional `derived_from`
  - figure-oriented `preferred_artifacts[]` now points at PNG bundle outputs
    rather than SVG paths
    - single-figure runs promote one rasterized `generated/...png`
    - multi-figure runs promote one best-first
      `generated/clawbio_storyboard.png`
    - original SVGs may still remain in the bundle as provenance/supporting
      artifacts, but only the best-first SVG is rasterized immediately so
      legacy media collectors do not auto-send additional PNGs
    - additional figures are represented by request-first `continue_artifact`
      actions whose nested request narrows `expected_artifacts[]`
  - `artifact_summary` may carry
    `gentle.clawbio_artifact_bundle_summary.v1` with:
    - `best_first_artifact`
    - `preferred_artifact_count`
    - `displayable_artifact_count`
    - `collected_artifact_count`
    - `continuation_action_count`
    - short `summary_lines[]` suitable for chat/report previews
  - requests may opt into `claim_attribution_mode = strict`; the wrapper then
    emits `gentle.clawbio_claim_ledger.v1` in both
    `result.json.claim_ledger` and `claim_ledger.json`
    - optional request `input_claims[]` contains non-empty caller-supplied
      statements; strict reports render each with `[input]`, and
      `claim_ledger.input_claims[]` records its deterministic claim id plus the
      exact `/request/input_claims/N` evidence pointer and source-value hash
    - display prefixes distinguish `[gentle]` executable-derived statements,
      `[clawbio]` orchestration/presentation, `[input]` caller assumptions,
      named `[external:TOOL]` results, and `[unverified]` prose
    - `[gentle]` attribution requires an exact JSON pointer or source scope,
      source-value/scope SHA-256, and named deterministic projection; prose is
      downgraded when this binding is unavailable
    - reserved prefixes in caller text are escaped and cannot self-assign
      `[gentle]` or `[external:TOOL]` authority
    - `processing_steps[]` records the request digest and exact GENtle command;
      `claims[]`, `input_claims[]`, `warning_claims[]`, and `artifacts[]`
      retain source tool, evidence pointer, and scientific-content authority
    - raw `stdout_json` remains unchanged; strict attribution is an additive
      audit/presentation layer and cannot upgrade any GENtle evidence state
  - every run carries `execution_manifest` with schema
    `gentle.clawbio_execution_manifest.v1` and writes the same payload to
    `reproducibility/execution_manifest.json`
    - `request_binding` records the normalized request SHA-256, raw request-file
      SHA-256, and content-bound input files
    - `state_binding.before|after` records the explicit state-file path,
      presence, size, and SHA-256
    - `wrapper` binds the generic skill script/catalog/descriptor bytes;
      `runtime` records resolver details, mandatory delegated-run GENtle
      version preflight, and hashes for resolvable launcher/image/binary files
    - `execution_steps[]` separates process execution from
      `native_result.reported_status_bindings[]`; a completed computation may
      retain a native biological `fail` without becoming an execution failure
    - `native_result` binds raw stdout/stderr, parsed JSON, and stable native
      report/run/operation identifiers; `artifacts[]` binds emitted files
    - the receipt is topic-neutral and contains no discussion, participant,
      permission, vote, freshness, withdrawal, or interpretation state
  - optional request `input_bindings[]` accepts `binding_id`, `path`, optional
    `role`/`media_type`, and optional expected SHA-256; missing or mismatched
    files fail before GENtle execution, while mutation or removal during a run
    produces `verification_failed` and retains both pre-run and post-run
    observations in the receipt
  - optional request `delegation` uses
    `gentle.clawbio_skill_delegation.v1`
    - the wrapper verifies the co-deployed source skill catalog, descriptor,
      version, intent, plan step, hashes, and exact `delegate_contract`
    - optional `resolved_slots` are recorded, while the normalized resolved
      wrapper request is always bound independently
    - routing provenance records the selected route but does not claim that
      natural-language dispatch is deterministic
    - confirmation-gated delegated routes return `approval_required` without
      running the scientific command and write
      `gentle.clawbio_execution_proposal.v1`; its canonical
      `approval_basis` binds the normalized request, exact command,
      source and generic-wrapper descriptors/contracts, runtime, state,
      scientific inputs, resolved paths, selected material/reference spaces,
      material tool/cache environment, and pinned primer backend
      - the included local Cargo launcher/fallback is reduced to the built,
        hash-bound `gentle_cli` executable before proposal creation; an OCI
        runtime must use an immutable image digest rather than a mutable tag
    - `gentle.clawbio_approved_execution_request.v1` references the stored
      proposal and carries a `gentle.clawbio_execution_approval.v1` assertion
      for its exact digest. The wrapper reloads the request from that proposal
      and rejects descriptor, command, runtime, material-environment, state,
      input, path, or backend drift before execution
    - approver identity and authorization belong to the caller/OpenClaw
      control plane. GENtle verifies the assertion/digest binding but does not
      claim to authenticate the approver
    - direct structured requests without delegation remain backward
      compatible; read-only delegated list/show/preflight routes may remain
      automatic
  - optional post-run discussion-moderation handoff:
    - `discussion_moderation_handoff.py` consumes an explicit
      `gentle.clawbio_skill_result.v1`, its matching
      `gentle.clawbio_execution_manifest.v1`, and caller context schema
      `gentle.clawbio_discussion_moderation_handoff_context.v1`
    - it emits `gentle.clawbio_discussion_moderation_handoff.v1` plus a
      content-addressed `gentle.clawbio_native_result_artifact.v1`
    - the adapter's stable ID, semantic version, and script SHA-256 are recorded
      in both the request environment and provider receipt
    - every caller input ref must resolve exactly once to a stable
      `request_binding.content_bound_inputs[]` row with the same before/after
      size and SHA-256; inferred, fuzzy, or path-only joins are rejected
    - `analysis_run_intake.request_digest` uses discussion-moderation's
      canonical digest fields and excludes paths, timestamps, raw request-file
      serialization, and receipt identity so relocated exact replay remains
      idempotent
    - the packet carries a computed typed-evidence object and exact provider
      hashes but does not open or edit a ledger; the consumer must validate
      permissions/current refs, record output evidence first, recompute the
      digest, and then record the analysis run through its owning commands
    - provider execution outcome, evidence-envelope validation, and native
      biological statuses stay distinct. A completed scientific `fail` remains
      a computed failure result, while failed execution has no output evidence
      and incomplete or `not_run` evidence is never promoted to success
  - browser/OpenClaw inline image display remains a later ClawBio-side task;
    this phase is limited to PNG-first bundle production inside `gentle_rs`
  - wrapper `agent-plan` / `agent-execute-plan` modes intentionally share the
    typed GENtle planner boundary instead of routing machine consumers through
    the chat-oriented `agents ask` UX
  - wrapper-level action envelopes in `result.json.suggested_actions[]`,
    `result.json.preferred_demo_actions[]`, and
    `result.json.blocked_actions[].action` are self-describing enough for
    ClawBio/OpenClaw to retain them across one chat session and execute a
    user-confirmed follow-up without re-inferring the originating skill:
    - required fields: `action_id`, `label`, `kind`, `skill_alias`, and
      `requires_confirmation`
    - optional fields: `rationale`, `timeout_secs`, `shell_line`,
      `expected_artifacts`, `resource_key`, `lifecycle_status`, and
      skill-opaque nested `request`
    - `skill_alias` is currently `gentle-cloning` for this wrapper and is the
      routing key a generic OpenClaw action-retention layer should use
    - `action_id` is deterministic for a given producer/label within a release
      and is suitable for selecting one action from the current result/session;
      it is not a global permanent identifier unless a specific action kind
      documents stronger stability
    - confirmation invariant: when a nested action `request` carries
      `confirm=true`, the wrapper-level envelope must also set
      `requires_confirmation=true`; the inverse is not required because some
      chat-level confirmations protect long-running, external, or setup-like
      shell routes that do not need a nested GENtle `confirm` flag
- service handoff payload: `gentle.service_handoff.v1`
  - produced by `services handoff [--scope NAME] [--output PATH]`
  - embeds the normal `gentle.service_readiness.v1` status as
    `service_readiness`
  - adds chat-gateway decision fields:
    - `readiness[]` with `resource_key`, `display_name`, `resource_kind`,
      `prepared`, `lifecycle_status`, and short status text
    - `suggested_actions[]` with deterministic GENtle shell commands that are
      safe to offer next
    - `running_actions[]` with refresh/status commands for already-active
      shared prepares
    - `blocked_actions[]` for useful but not immediately executable setup,
      such as ATtRACT sync before a local ZIP path is known
    - `preferred_demo_actions[]` for low-friction demo commands that can be
      shown by ClawBio without inventing route logic
    - `status_overview` with lifecycle counts, an overall
      `ready|setup_needed|setup_running|attention_needed` state, and one
      recommended next action for compact chat handoff
    - `environment_hints[]` for deployment variables such as
      `GENTLE_CLI_CMD`, `GENTLE_REPO_ROOT`, `GENTLE_REFERENCE_CACHE_DIR`,
      `GENTLE_HELPER_CACHE_DIR`, and `GENTLE_CUTRUN_CACHE_DIR`
  - the ClawBio wrapper normalizes `stdout_json.suggested_actions[]`,
    `stdout_json.preferred_demo_actions[]`, and
    `stdout_json.blocked_actions[]` into top-level `result.json` fields,
    adding nested `gentle.clawbio_skill_request.v1` payloads where an action is
    executable so conversational confirmation/execution does not need to parse
    the raw GENtle payload
  - `skill-info` now also reports `ui_intent_support` metadata describing how
    the wrapper surfaces shared UI-intent discovery through `capabilities`
  - `capabilities` now performs a best-effort auxiliary `ui intents` probe
    through the same resolved runtime
    - on success, `result.json.ui_intent_catalog` carries the shared
      `gentle.ui_intents.v1` payload and `suggested_actions[]` gains
      `kind = ui_intent` entries with structured `ui_intent` metadata plus
      executable `ui open TARGET` shell requests
    - on failure or older runtimes, `result.json.ui_intent_catalog_error`
      captures the non-fatal probe failure while the main `capabilities`
      request still succeeds
- Telegram guide payload: `gentle.telegram_guide.v1`
  - produced by
    `services guide --channel telegram [--section SECTION] [--gene SYMBOL]`
  - intended for bench-user chat orientation, while
    `services handoff --scope clawbio` remains the operator/setup route
  - fields:
    - `summary_lines[]` with a compact introduction and optional gene prompt
    - `readiness_summary_lines[]` copied from the current service-readiness
      view
    - `menu_sections[]` for the available guide chapters
    - `suggested_actions[]` as lightweight navigation links
    - `blocked_actions[]` copied from GENtle-owned setup/readiness state
  - guide navigation uses normal action records with
    `kind = guide_section`, `requires_confirmation = false`, and a `shell_line`
    that invokes another guide section such as
    `services guide --channel telegram --section tfbs --gene TP73`
  - when `--gene` is omitted, gene-aware sections use deterministic defaults
    such as `TERT`/`TP73` for promoter-TFBS, `TP73` for gene context, and
    `TP73`/`TP53` for isoform demos
- reproducibility outputs:
  - `report.md`
  - `result.json`
  - `reproducibility/commands.sh`
  - `reproducibility/environment.yml`
  - `reproducibility/execution_manifest.json`
  - `reproducibility/checksums.sha256`
- included bootstrap example requests:
  - `request_services_status.json`
  - `request_services_telegram_guide.json`
  - `request_services_handoff.json`
  - `request_genomes_list_human.json`
  - `request_genomes_status_grch38.json`
  - `request_genomes_prepare_grch38.json`
  - `request_helpers_status_puc19.json`
  - `request_helpers_prepare_puc19.json`
- included follow-on example requests:
  - `request_genomes_extract_gene_tp53.json`
  - `request_helpers_blast_puc19_short.json`
  - `request_workflow_vkorc1_planning.json`
  - `request_protocol_cartoon_gibson_svg.json`
    - declares `expected_artifacts[]` so the generated SVG is copied into the
      output bundle under `generated/...` and rasterized into the PNG-first
      ClawBio contract
  - `request_cdna_pcr_test_demo_direct.json`,
    `request_cdna_qpcr_taqman_test_demo_direct.json`, and
    `request_workflow_cdna_pcr_qpcr_assay_test_offline.json`
    - declare cDNA assay transcript-map SVG artifacts generated by
      `--svg` / `svg_path`, so ClawBio can show where PCR/qPCR products are
      functional across the tested transcript set
  - `request_workflow_cdna_pcr_qpcr_product_gel_nonspecific_offline.json`,
    `request_cdna_pcr_products_gel_demo_direct.json`, and
    `request_cdna_qpcr_taqman_products_gel_demo_direct.json`
    - request product materialization / `product_gel_svg_path`, so detected
      cDNA PCR/qPCR products are promoted into GENtle sequences, grouped into
      one vial/container, and exported as product-gel SVG artifacts for
      ClawBio-first display

Planned operation refinements:

- `MergeContainers { inputs, output_prefix? }`
  - Explicitly models wet-lab mixing of multiple tubes/pools.
- Protocol-based ligation:
  - `Ligation { input_container, protocol, output_container?, ... }`
  - `protocol` determines allowed end joins.
  - Initial protocol values:
    - `sticky`
    - `blunt`
  - Future protocol values may include established ligation workflows
    represented as named presets.

Current parameter support:

- `max_fragments_per_container` (default `80000`)
  - limits digest fragment output per operation
  - also serves as ligation product-count limit guard
- `require_verified_genome_anchor_for_extension` (default `false`)
  - when `true`, `ExtendGenomeAnchor` requires anchor provenance with
    `anchor_verified=true`
  - anchors with `anchor_verified=false` or missing verification status are
    rejected in strict mode
  - alias parameters accepted: `strict_genome_anchor_verification`,
    `strict_anchor_verification`
- `genome_anchor_prepared_fallback_policy` (default `single_compatible`)
  - controls how `ExtendGenomeAnchor` / `VerifyGenomeAnchor` resolve anchor
    genome ids when exact prepared cache id is not present.
  - accepted values:
    - `off` (no compatibility fallback; must match exact prepared id)
    - `single_compatible` (auto-fallback only when one compatible prepared
      cache exists)
    - `always_explicit` (never auto-fallback; require explicit selection even
      when only one compatible prepared cache exists)
  - alias parameters accepted: `genome_anchor_fallback_mode`,
    `genome_anchor_prepared_mode`
- primer-design backend controls:
  - `primer_design_backend` (default `auto`)
    - accepted values: `auto`, `internal`, `primer3`
    - `auto` tries Primer3 and falls back deterministically to internal scoring
      with explicit warning + fallback reason in report metadata
  - `primer3_executable` (default `"primer3_core"`)
    - executable path/name used when backend is `primer3` or `auto`
    - alias parameters accepted: `primer3_backend_executable`, `primer3_path`
- `feature_details_font_size` (default `9.0`, range `8.0..24.0`)
  - controls GUI font size for the feature tree entries and feature range details
- `regulatory_feature_max_view_span_bp` (default `50000`, range `>= 0`)
  - hides regulatory feature overlays in linear view when current view span
    exceeds this threshold (`0` disables regulatory overlays)
- `gc_content_bin_size_bp` (default `100`, range `>= 1`)
  - controls GC-content aggregation bin size for linear/circular rendering and
    SVG export
- Linear DNA-letter routing parameters:
  - `linear_show_sequence_bases` (default `true`)
    - master visibility switch for DNA base letters in the linear map
    - aliases accepted by `SetParameter`:
      - `linear_show_sequence_letters`
      - `show_linear_sequence_bases`
      - `show_linear_sequence_letters`
    - when `false`, adaptive/forced modes still report routing diagnostics, but
      the active base-letter mode is `OFF`
  - `linear_sequence_letter_layout_mode` (default `AutoAdaptive`)
    - supported canonical modes:
      - `auto|adaptive|auto_adaptive`
      - `standard|standard_linear`
      - `helical|continuous_helical`
      - `condensed_10_row|condensed`
    - auto mode uses deterministic viewport-density tiers:
      - `<= 1.5x`: standard
      - `<= 2x`: helical (if compressed letters enabled)
      - `<= 10x`: condensed-10 (if compressed letters enabled)
      - `> 10x`: `OFF`
  - `linear_sequence_helical_letters_enabled` (default `true`)
    - applies to auto mode only (allows/disallows compressed auto tiers)
  - `linear_sequence_helical_phase_offset_bp` (range `0..9`)
    - seam offset used by helical/condensed row mapping
  - reverse/helical strand geometry controls:
    - `linear_show_double_strand_bases` / `linear_show_reverse_strand_bases`
      (bool alias pair; controls reverse-strand letter visibility)
    - `linear_helical_parallel_strands` (default `true`)
      - `true`: forward/reverse helical slant stays parallel
      - `false`: forward/reverse helical slant is mirrored (cross-over look)
    - `reverse_strand_visual_opacity` (range `0.2..1.0`, default `0.55`)
      - shared reverse-strand emphasis in linear map and sequence panel
- Legacy linear-letter threshold knobs are compatibility-only and return
  deterministic deprecated no-op messages (no routing effect):
  - `linear_sequence_base_text_max_view_span_bp`
  - `linear_sequence_helical_max_view_span_bp`
  - `linear_sequence_condensed_max_view_span_bp`
- VCF display filter parameters (shared GUI/SVG state):
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
  - `vcf_display_required_info_keys` (CSV string or string array)
- TFBS display filter parameters (shared GUI/SVG state):
  - `show_tfbs`
  - `tfbs_display_use_llr_bits`
  - `tfbs_display_min_llr_bits`
  - `tfbs_display_use_llr_quantile`
  - `tfbs_display_min_llr_quantile`
  - `tfbs_display_use_true_log_odds_bits`
  - `tfbs_display_min_true_log_odds_bits`
  - `tfbs_display_use_true_log_odds_quantile`
  - `tfbs_display_min_true_log_odds_quantile`
- Repeat display parameter (shared GUI/SVG state):
  - `show_repeat_features`
  - equivalent visibility target: `SetDisplayVisibility { target: "RepeatFeatures", visible }`
- Microarray-array display parameter (shared GUI/SVG state):
  - `show_array_features`
  - equivalent visibility target: `SetDisplayVisibility { target: "ArrayFeatures", visible }`
- Restriction-enzyme display parameters (shared GUI/SVG state):
  - `show_restriction_enzymes`
  - `show_restriction_enzyme_sites` (bool alias)
  - `restriction_enzyme_display_mode`
  - `restriction_display_mode` (string alias)
    - supported values:
      - `preferred_only`
      - `preferred_and_unique`
      - `unique_only`
      - `all_in_view`
  - `preferred_restriction_enzymes`
  - `preferred_restriction_enzymes_csv` (CSV alias)
  - `restriction_preferred_enzymes` (CSV/string-array alias)
    - accepts either a CSV string or a string array
- BLAST options-layer parameters:
  - `blast_options_override` (JSON object or `null`)
    - project-level BLAST option layer merged before per-command request JSON
    - supports the same keys as request JSON (`task`, `max_hits`, `thresholds`)
  - `blast_options_defaults_path` (string path or `null`)
    - optional defaults-file path used ahead of project/request layers
    - if unset, engine falls back to `assets/blast_defaults.json`

Current ligation protocol behavior:

- `protocol` is mandatory.
- If `protocol = Blunt`, ligation enumerates ordered input pairs with blunt-end
  compatibility checks.
- If `protocol = Sticky`, ligation enumerates ordered input pairs with sticky-end
  overhang compatibility checks.
- `unique = true` requires exactly one product.

`FilterByMolecularWeight` semantics:

- Applies a bp-range filter across provided input sequence ids.
- Effective accepted range is expanded by `error`:
  - `effective_min = floor(min_bp * (1 - error))`
  - `effective_max = ceil(max_bp * (1 + error))`
- `unique = true` requires exactly one match, otherwise the operation fails.

`FilterByDesignConstraints` semantics:

- Applies practical design-constraint filters across provided input sequence ids.
- Optional GC bounds:
  - `gc_min` and/or `gc_max` (fractional range `0.0..1.0`)
  - when both are provided, `gc_min <= gc_max` is required
- Optional homopolymer cap:
  - `max_homopolymer_run >= 1`
  - rejects candidates with a longer A/C/G/T run
- `reject_ambiguous_bases` (default `true`):
  - rejects sequences containing non-ACGT letters
- `avoid_u6_terminator_tttt` (default `true`):
  - rejects sequences containing `TTTT`
- Optional `forbidden_motifs`:
  - IUPAC motifs; reject when motif appears on either strand
- `unique = true` requires exactly one match, otherwise the operation fails.

Guide-design semantics:

- Guide sets persist in `ProjectState.metadata["guide_design"]`
  (`schema = gentle.guide_design.v1`) and include:
  - guide sets
  - practical-filter reports
  - oligo sets
  - audit log entries for guide operations/exports
- `UpsertGuideSet`:
  - normalizes guide fields and validates required properties
  - sorts by rank (then guide id) and rejects duplicate `guide_id` within one set
- `FilterGuidesPractical`:
  - applies deterministic practical filters over one guide set
  - supports GC bounds, global/per-base homopolymer limits, ambiguous-base
    rejection, U6 `TTTT` avoidance, dinucleotide repeat cap, forbidden motifs,
    and required 5' base checks
  - can emit a passed-only output guide set (`output_guide_set_id`)
  - always persists a structured per-guide report with reasons/warnings/metrics
- `GenerateGuideOligos`:
  - generates forward/reverse oligos using a named template
  - supports optional 5' G extension and passed-only mode
  - persists generated oligo records in named oligo sets
- `ExportGuideOligos`:
  - exports an oligo set as `csv_table`, `plate_csv` (96/384), or `fasta`
  - records export actions in the guide-design audit log
- `ExportGuideProtocolText`:
  - exports a deterministic human-readable protocol text artifact
  - optional QC checklist can be included/excluded

Candidate-set semantics:

- `GenerateCandidateSet` creates a persisted candidate window set over one source
  sequence and computes baseline metrics for each candidate.
- `GenerateCandidateSetBetweenAnchors` creates a persisted candidate window set
  constrained to the in-sequence interval between two local anchors.
- `ScoreCandidateSetExpression` computes a derived metric from an arithmetic
  expression over existing metrics.
- `ScoreCandidateSetDistance` computes feature-distance metrics against filtered
  feature targets.
- `FilterCandidateSet` keeps/drops candidates by absolute bounds and/or quantile
  bounds for a named metric.
- `CandidateSetOp` supports set algebra (`union`, `intersect`, `subtract`) over
  candidate identity (`seq_id`, `start_0based`, `end_0based`).
- `ScoreCandidateSetWeightedObjective` computes one metric from weighted
  objective terms (`maximize`/`minimize` per term, optional normalization).
- `TopKCandidateSet` selects an explicit top-k subset for one metric with a
  deterministic tie-break policy.
- `ParetoFrontierCandidateSet` keeps non-dominated candidates for multiple
  objectives (`maximize`/`minimize` per objective), with optional tie-break
  truncation.
- Workflow macro templates are persisted in project metadata:
  - `UpsertWorkflowMacroTemplate` stores/replaces named templates
  - `DeleteWorkflowMacroTemplate` removes templates
  - each template now carries `template_schema`
    (`gentle.cloning_macro_template.v1`) so cloning-operation macro intent is
    explicit at engine level
  - optional `details_url` can link to external protocol/reference material
  - optional typed `input_ports`/`output_ports` can be persisted directly in
    template metadata (same port shape as routine catalog ports)
  - template expansion/binding is exposed through adapter command surfaces
    (`macros template-*`, including `macros template-import PATH`)
  - expanded scripts run through shared shell execution (`macros run`) and can
    orchestrate full cloning operations via `op ...` or `workflow ...` payloads
  - shipped starter assets:
    - legacy pack: `assets/cloning_patterns.json` (`gentle.cloning_patterns.v1`)
    - hierarchical catalog: `assets/cloning_patterns_catalog/**/*.json`
      (`gentle.cloning_pattern_template.v1`, one template per file)
    - Gibson baseline template:
      `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`
    - Restriction baseline template:
      `assets/cloning_patterns_catalog/restriction/digest_ligation/digest_ligate_extract_sticky.json`
- Typed cloning-routine catalog baseline:
  - manifest: `assets/cloning_routines.json`
  - schema: `gentle.cloning_routines.v1`
  - typed routine metadata fields include routine family/status/tags, linked
    template name/path, and typed input/output port declarations
  - includes Gibson + restriction family baselines:
    - `gibson.two_fragment_overlap_preview`
    - `restriction.digest_ligate_extract_sticky`
  - adapter discovery surface:
    `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT] [--seq-id SEQ_ID]`
  - explainability and comparison surfaces:
    - `routines explain ROUTINE_ID [--catalog PATH] [--seq-id SEQ_ID]`
    - `routines compare ROUTINE_A ROUTINE_B [--catalog PATH] [--seq-id SEQ_ID]`
- Macro-instance lineage baseline:
  - mutating `macros run` / `macros template-run` append one
    `LineageMacroInstance` record in project lineage state for success and
    failure pathways
  - records include deterministic `macro_instance_id`, optional
    `routine_id/template_name`, typed bound inputs/outputs, emitted `op_id`s,
    status, and optional `status_message`
  - lineage graph + lineage SVG consume these records as macro box nodes with
    explicit input/output edges where sequence/container references resolve
- Candidate macro templates are persisted in project metadata:
  - `UpsertCandidateMacroTemplate` stores/replaces named templates
  - `DeleteCandidateMacroTemplate` removes templates
  - optional `details_url` can link to external protocol/reference material
  - template expansion/binding is exposed through adapter command surfaces
    (`candidates template-*`)
- Between-anchor generation augments baseline metrics with anchor-aware fields
  (`distance_to_anchor_a_bp`, `distance_to_anchor_b_bp`,
  `distance_to_nearest_anchor_bp`, interval span metadata).

Feature-distance geometry controls (candidate generation and distance scoring):

- `feature_geometry_mode` (optional, default `feature_span`):
  - `feature_span`: one interval per feature using whole-feature bounds
  - `feature_parts`: one interval per explicit location part (ignores intronic
    gaps for multipart features)
  - `feature_boundaries`: boundary points of explicit location parts
- `feature_boundary_mode` (optional, default `any`):
  - `any`, `five_prime`, `three_prime`, `start`, `end`
  - only meaningful when `feature_geometry_mode = feature_boundaries`
- `feature_strand_relation` (optional, default `any`):
  - `any`, `same`, `opposite`
  - current engine interpretation is sequence-forward relative
    (`same = '+'`, `opposite = '-'` feature strand)
- Directed-boundary interpretation:
  - on `+` strand: `five_prime = start`, `three_prime = end`
  - on `-` strand: `five_prime = end`, `three_prime = start`
  - for unknown strand, `five_prime`/`three_prime` conservatively include both
    boundaries.

`RenderPoolGelSvg` semantics:

- Accepts explicit `inputs` (sequence ids) and an output `path`; all explicit
  inputs become co-migrating members of one sample lane.
- Optional `container_ids` renders one lane per referenced stored container.
- Optional `arrangement_id` renders one lane per stored serial-arrangement lane
  in declared `lane_container_ids` order and carries matching role labels.
  `ArrangementMode::Plate` is rejected dynamically because plate placement is
  not a serial gel lane layout.
- Optional `conditions` carries one shared gel-run profile for the whole render:
  - `agarose_percent`
  - `buffer_model` (`tae` / `tbe`)
  - `topology_aware`
- Computes pool migration from sequence bp length plus one deterministic
  heuristic condition model:
  - agarose/buffer reshape the shared migration curve
  - topology-aware mode distinguishes:
    - generic `circular`
    - `supercoiled`
    - `relaxed circular`
    - `nicked/open circular`
    - `linear`
  - when sequence names/definitions/comments carry explicit circular-form hints,
    those richer forms are used; otherwise circular sequences fall back to the
    generic `circular` model
- Computes sample-band brightness from estimated DNA mass proxy, not only
  multiplicity.
- Applies one deterministic co-migration grouping threshold so nearby fragments
  that would plausibly collapse into one readout are exported as one merged
  band annotation instead of several visually indistinguishable bars.
- Chooses one or two ladders to span pool range:
  - from explicit `ladders` list when provided
  - otherwise from saved arrangement ladders when `arrangement_id` is present
    and the arrangement stores a ladder choice
  - otherwise from built-in ladder catalog (auto mode)
- Renders ladder lanes plus pooled band lanes as SVG artifact.
- Optional `render_options` is presentation-only and defaults to:
  - `lane_label_layout: auto`: preserve horizontal short names, wrap isolated
    long names, and angle adjacent/difficult long names; accepted explicit
    values are `horizontal`, `wrapped`, `staggered`, and `angled`
  - `band_label_layout: auto`: retain an in-gel band annotation only when its
    estimated monospace width fits before the next lane or gel edge; `panel`
    keeps detail text solely in the fragment table and `inline` requests the
    legacy always-inline behavior
  - `isoform_marker_mode: auto`: detect stable Ensembl and RefSeq transcript
    accessions embedded in product identifiers, normalize away accession version
    suffixes, and assign each identity a deterministic color, relative marker
    position, and `O`/`I` binary legend code; `off` suppresses this layer
  - merged bands retain one marker per represented isoform; marker position and
    color repeat across lanes, while the binary code is repeated in the fragment
    table as a color-independent textual join key
  - wrapped/angled labels may increase only the SVG canvas height; lane, band,
    migration, and fragment-table biology remain unchanged
- The SVG also includes a compact fragment table for non-ladder lanes:
  - observed apparent size
  - actual bp
  - topology-form hint
  - estimated mass proxy
  - merged-band annotation when several fragments co-migrate
- Shared-shell `render-pool-gel-svg` output additionally returns
  `gel_band_rows[]` and `gel_summary_lines[]`, so text-first frontends such as
  Telegram can explain lane/band sizes without parsing SVG text or depending on
  host font rasterization.
  - source fragment labels
- When serial lanes carry interpretable roles such as `vector`, `insert_*`, and
  `product`, the right-hand detail panel also adds `Comparison hints` for:
  - insert vs fine ladder sizing
  - vector vs product size-shift reading
  - simple product-minus-vector vs summed-insert consistency checks
  - role badges under the lane labels so generic container names still read as
    `VECTOR`, `INSERT`, or `PRODUCT`
- When merged bands are present, the right-hand detail panel also adds a short
  `Merged-band notes` block that explains observed band position vs the
  underlying actual source-size span.

`RenderProteinGelSvg` semantics:

- Accepts one persisted `ProteinDerivationReport` id plus an output `path`.
- Reuses the first-class protein sequence entries already created by
  `DeriveProteinSequences`; it does not re-derive proteins from DNA at render
  time.
- Computes molecular weight from amino-acid sequence content and plots bands
  on a log-scaled kDa axis.
- Chooses one or two protein ladders deterministically, or honors the explicit
  `ladders` override when provided.
- Renders one lane per derived protein plus the selected ladder lane(s), using
  the derived protein `product` qualifier as the visible label when available
  and keeping transcript/protein accessions in the detail text. A provenance
  panel records the source report and selection summary.
- This is the canonical protein-gel route for transcript-native demos such as
  the TP73 isoform workflow; `RenderPoolGelSvg` remains DNA/bp-based.

`RenderProteinGelReportsSvg` semantics:

- Accepts two or more persisted `ProteinDerivationReport` ids plus an output
  SVG `path`; one report is also accepted for deterministic scripting parity.
- Reuses the same first-class protein sequence entries as `RenderProteinGelSvg`
  and never re-derives proteins at render time.
- Renders one gel column per report/gene and one band per derived isoform
  product inside that column, with ladder lanes placed at the side(s) by the
  shared protein-ladder resolver.
- This is the canonical multi-gene 1D protein-gel route for ClawBio guide
  figures such as PATZ1/TP73/TP53/TP63/SP1/BACH2 isoform molecular-weight
  comparisons.

`RenderProteaseDigestGelSvg` semantics:

- Accepts either one first-class protein/peptide `seq_id`, or one
  `ProteinDerivationReport` `report_id` plus an optional `transcript_id` for
  selecting the derived protein row, one or more protease names or aliases, and
  an output SVG `path`.
- The `report_id` path is preferred inside workflows that derive proteins
  immediately before rendering, because it follows the actual sequence id
  created in the current state even when prior runs forced a uniqueness suffix.
- Uses the same built-in protease catalog and cleavage matcher as
  `ProteaseDigestProteinSequence`, but runs in prediction-only mode and does
  not materialize peptide sequences.
- Returns the same `gentle.protease_digest_report.v1` in `OpResult` while also
  rendering one protein-gel lane per retained peptide product.
- `min_length_aa` filters short peptides before rendering; `ladders` can
  request deterministic protein ladder overlays.
- This is the ClawBio/OpenClaw-friendly graphical insight path for protease
  digestion because declared SVG artifacts are rasterized by the wrapper into
  PNG-first messenger outputs.

`RenderProtein2dGelSvg` semantics:

- Accepts one persisted `ProteinDerivationReport` id plus an output `path`.
- Reuses the same first-class protein sequence entries created by
  `DeriveProteinSequences`; it does not re-derive proteins at render time.
- Estimates isoelectric point from amino-acid composition using the shared
  deterministic Henderson-Hasselbalch model in `AminoAcids`.
- Places one protein spot per derived protein on a 2D gel with pI on the
  X axis and log-scaled molecular weight on the Y axis.
- Uses the same deterministic ladder-selection logic for Y-axis references as
  the 1D protein-gel renderer, or honors the explicit `ladders` override.
- Renders plot-internal pI and molecular-weight axis labels plus ladder-derived
  kDa tick labels so messenger PNG previews remain interpretable without
  relying on surrounding report text.
- Renders derived protein `product` labels when available, keeps transcript /
  protein accessions in the spot details, labels spots directly in the plot,
  and records the source report plus transcript-selection summary in the
  provenance panel.
- This is the canonical 2D protein-gel route for transcript-native spot-map
  demos and for ClawBio's parameterized Ensembl `gene-protein-2d-gel` request
  mode; `RenderProteinGelSvg` remains the 1D lane-based route and
  `RenderPoolGelSvg` remains DNA/bp-based.

`CreateArrangementSerial` semantics:

- Persists an ordered serial-lane setup over stored containers.
- Optional `ladders` can store one symmetric ladder or one left/right ladder
  pair for later gel preview/export reuse.
- Also materializes one default physical rack draft:
  - choose the smallest built-in rack/plate profile that fits the arrangement
    payload plus ladder-reference positions
  - place the arrangement block row-major in that rack
  - link the arrangement back to the resulting `default_rack_id`

`SetArrangementLadders` semantics:

- Mutates an existing serial arrangement in place.
- `ladders = null` clears back to shared engine auto ladder selection.
- one ladder name means the same ladder is used on both sides during
  arrangement-based gel preview/export.
- two ladder names mean explicit left/right ladder selection.

`SetContainerDeclaredContentsExclusive` semantics:

- Mutates one persisted container’s `declared_contents_exclusive` flag.
- `exclusive = true` means the listed members are intended to be the full known
  contents of the vial/tube.
- `exclusive = false` means the listed members are known/measured constituents
  of a more complex sample that may also contain unlisted molecules.

`CreateRackFromArrangement` semantics:

- Creates one new physical rack/plate draft from one stored arrangement.
- Optional `profile` overrides the default smallest-fitting profile choice.
- If `profile` is omitted, the engine chooses in this order:
  - `small_tube_4x6`
  - `plate_96`
  - `plate_384`
- `plate_6` is available for explicit cell-culture layouts, but it is not the
  implicit default for small cloning/bench arrangements.
- Placement is row-major and preserves arrangement order.
- Ladder-bearing arrangements reserve left/right ladder-reference positions in
  the same contiguous block.

`PlaceArrangementOnRack` semantics:

- Places one arrangement onto an existing rack as one contiguous block at the
  next free region in fill order.
- Existing rack occupants stay in order; the appended arrangement does not
  reorder earlier blocks.
- Shared racks are therefore possible without losing arrangement identity.

`MoveRackPlacement` semantics:

- Moves one occupied rack coordinate within one saved rack.
- `move_block=false` means move one sample within its arrangement block and
  shift neighboring occupied positions to preserve order.
- `move_block=true` means move the whole arrangement block and shift later
  occupied blocks in fill order.
- The operation is order-preserving by design; it does not treat arbitrary
  holes as the primary editing model.

`MoveRackSamples` semantics:

- Moves two or more selected samples together within one saved rack.
- `from_coordinates[]` must all resolve to occupied positions from the same
  arrangement block.
- The shared engine normalizes the selected samples against the rack's current
  occupied order; it preserves that rack order even if the request lists the
  source coordinates differently.
- The selected samples move as one contiguous combined group within that same
  arrangement block.
- Neighboring occupied positions shift in fill order to keep the block
  contiguous.
- This is the shared engine contract behind rack-editor sample multi-select
  moves.

`MoveRackArrangementBlocks` semantics:

- Moves two or more selected arrangement blocks together within one saved rack.
- `arrangement_ids` are normalized against the rack's current occupied order;
  the shared engine preserves the rack-ordering of selected blocks even if the
  request lists them in another order.
- The selected blocks move as one contiguous combined group.
- Later occupied blocks shift in fill order to keep the rack contiguous.
- This is the shared engine contract behind rack-editor multi-select moves.

`SetRackProfile` semantics:

- Reprojects one saved rack onto another built-in profile.
- Existing arrangement order is preserved while coordinates are reflowed under
  the target profile geometry.
- Existing fill direction is preserved.
- Existing blocked coordinates are preserved when still in-bounds for the new
  geometry; out-of-bounds blocked coordinates are dropped deterministically.

`ApplyRackTemplate` semantics:

- Applies one engine-owned quick-authoring template on top of an existing rack
  snapshot.
- Built-in templates:
  - `bench_rows`
    - `fill_direction = row_major`
    - `blocked_coordinates = []`
  - `plate_columns`
    - `fill_direction = column_major`
    - `blocked_coordinates = []`
  - `plate_edge_avoidance`
    - `fill_direction = column_major`
    - `blocked_coordinates = outer perimeter of the current profile`
- Existing arrangement order is preserved while occupied coordinates are
  reflowed onto the resulting available slots.
- `plate_edge_avoidance` requires at least a `3 x 3` profile so an interior
  region remains after blocking the perimeter.

`SetRackFillDirection` semantics:

- Reprojects one saved rack onto the same geometry with a different fill order.
- Supported values:
  - `row_major`
  - `column_major`
- Existing arrangement order is preserved while occupied coordinates are
  reassigned under the new fill order.

`SetRackProfileCustom` semantics:

- Reprojects one saved rack onto one custom A1-style geometry.
- `rows` and `columns` are persisted directly in the rack profile snapshot.
- Existing fill direction is preserved.
- Existing blocked coordinates are preserved when still in-bounds for the new
  geometry; out-of-bounds blocked coordinates are dropped deterministically.
- Existing arrangement order is preserved while coordinates are reflowed under
  the custom geometry.
- A1-style row labels continue beyond `Z` as `AA`, `AB`, ...

`SetRackBlockedCoordinates` semantics:

- Persists one normalized blocked/reserved coordinate set on the rack profile.
- Blocked coordinates are excluded from placement capacity and fill-order
  reflow.
- Existing arrangement order is preserved while occupied coordinates are
  reassigned onto the remaining available positions.
- Duplicate blocked coordinates are removed deterministically.

`ExportRackLabelsSvg` semantics:

- Writes one deterministic SVG label sheet for a saved rack.
- Optional `arrangement_id` restricts output to labels belonging to one
  arrangement block on that rack.
- `preset` is engine-owned and defaults to `compact_cards`.
- Built-in presets:
  - `compact_cards`
  - `print_a4`
  - `wide_cards`
- The sheet now carries one scope header:
  - rack title
  - arrangement scope when requested
  - preset id
- Individual cards prioritize printable lane identity over verbose provenance:
  - rack id + coordinate
  - role
  - ladder name or sequence id / pool summary
  - bp length/topology when sequence-backed
  - container id as a final compact provenance hint when space remains
- Verbose arrangement/origin text is intentionally not repeated on every card;
  it lives in the shared sheet header instead so compact-card exports remain
  printable.

`ExportRackFabricationSvg` semantics:

- Writes one deterministic top-view fabrication/planning SVG for a saved rack.
- Uses one engine-owned physical carrier template layered on the current saved
  rack snapshot rather than inventing a second placement model.
- Current built-in physical templates:
  - `storage_pcr_tube_rack`
  - `pipetting_pcr_tube_rack`
  - `cell_culture_plate`
- The export consumes:
  - rack geometry (`rows`, `columns`, blocked coordinates)
  - saved rack occupancy and arrangement ids for visual planning markers
- Intended downstream uses:
  - fabrication sketching
  - bench planning
  - future simulation adapters

`ExportRackIsometricSvg` semantics:

- Writes one deterministic pseudo-3D/isometric SVG for a saved rack.
- Uses the same engine-owned physical carrier template family as
  `ExportRackFabricationSvg` and `ExportRackOpenScad`.
- Consumes the same linked rack snapshot:
  - rows / columns / blocked coordinates
  - occupied placements
  - arrangement ids and colors
- Intended downstream uses:
  - README / documentation hero figures
  - bench-facing communication assets
  - presentation-ready rack review without leaving the shared engine path
- For the cell-culture template, occupied wells render as flat culture-well
  fills rather than PCR tubes.

`ExportRackHeroSvg` semantics:

- Writes one deterministic README-facing hero SVG for a saved rack.
- Supports the built-in `storage_pcr_tube_rack`, `pipetting_pcr_tube_rack`,
  and `cell_culture_plate` physical templates.
- Uses the same linked rack snapshot as the technical rack exports, but renders
  a presentation-specific top-down rack/plate:
  - a rack/plate outline with a clipped upper-left orientation corner
  - row and column labels
  - tight circular empty-slot rims and floors/caps
  - saved-arrangement labels and subtle rings for occupied slots
  - a bottom caption strip
- This route intentionally remains separate from `ExportRackOpenScad`; OpenSCAD
  stays the CAD/printing source, while hero SVG stays dependency-free and
  deterministic for documentation.

`ExportRackOpenScad` semantics:

- Writes one deterministic parameterized OpenSCAD source file for a saved rack.
- Uses the same engine-owned physical carrier template family as
  `ExportRackFabricationSvg`.
- Current OpenSCAD export intentionally favors geometry over embedded text:
  - tube openings
  - outer carrier body
  - front label-strip recess
- Printable labels/front strips remain separate shared exports rather than
  being baked permanently into the 3D geometry in this first baseline.

`ExportRackCarrierLabelsSvg` semantics:

- Writes one deterministic carrier-matched SVG sheet for a saved rack.
- Uses the same engine-owned physical carrier template family as
  `ExportRackFabricationSvg` and `ExportRackOpenScad`.
- Supports optional arrangement scoping:
  - whole-rack export when `arrangement_id` is omitted
  - one arrangement/module export when `arrangement_id` is provided
- Built-in presets:
  - `front_strip_and_cards`
  - `front_strip_only`
  - `module_cards_only`
- Current baseline can emit:
  - one front-strip label sized from the selected physical template
  - one module card per arrangement in scope

`ExportRackSimulationJson` semantics:

- Writes one deterministic machine-readable JSON export for downstream
  simulation adapters.
- Uses the same engine-owned physical carrier template family as the physical
  SVG/OpenSCAD exports.
- Current baseline includes:
  - rack/profile metadata
  - selected physical template geometry
  - arrangement block summaries
  - one slot record per physical coordinate with:
    - row/column/coordinate
    - fill ordinal
    - blocked status
    - physical center in mm
    - occupant/arrangement metadata when occupied

`RenderDotplotSvg` semantics:

- Inputs:
  - `seq_id` (owner sequence id for the stored dotplot payload)
  - `dotplot_id` (stored payload id from `ComputeDotplot` or `ComputeDotplotOverlay`)
  - `path` (output SVG)
  - optional `flex_track_id` (adds flexibility panel in same SVG)
  - optional `display_density_threshold` and `display_intensity_gain` (display tuning)
  - optional `overlay_x_axis_mode` for overlay payloads:
    `percent_length | left_aligned_bp | right_aligned_bp | shared_exon_anchor | query_anchor_bp`
  - optional `overlay_anchor_exon` for anchored overlay exports
    (`{ start_1based, end_1based }`)
- Ownership checks:
  - dotplot payload must belong to `seq_id`
  - optional flexibility track must also belong to `seq_id`
- Output:
  - deterministic SVG dotplot artifact; operation is non-mutating
  - pairwise, annotated self, and overlay payloads render a genome-context side
    rail beside the reference/genomic y-axis when `reference_annotation`
    intervals are present; supported intervals include exon and materialized
    RepeatMasker/UCSC repeat context
  - genome-context rail SVG elements carry stable `data-gentle-role`
    attributes (`genome-context-track`, `genome-context-interval`,
    `genome-context-strand`, `genome-context-intron-guide`) plus
    `data-gentle-feature-kind` / `data-gentle-lane` markers for machine
    inspection
  - overlay payloads render all stored `query_series` with a legend;
    `overlay_x_axis_mode` chooses whether transcript queries are shown as
    normalized percent length, as left/right aligned base-pair coordinates, as
    shared-exon-anchored coordinates, or as explicit per-query-anchor
    coordinates.
  - `shared_exon_anchor` filters the overlay down to series that contain the
    selected reference exon and shifts each supporting transcript so the exon
    start lands at the maximal observed transcript-local start for that exon.
  - `query_anchor_bp` filters the overlay down to series carrying stored
    `query_anchor_0based` values and shifts each supporting query so the anchor
    base lands at the maximal observed anchor position across the drawn series.
  - flexibility panel is suppressed there because the overlay does not expose
    one shared query-axis coordinate system.

`DeriveTranscriptSequences` semantics:

- Inputs:
  - `seq_id`
  - optional `feature_ids[]`
  - optional splicing `scope`
  - optional `output_prefix`
- Behavior:
  - derives one spliced transcript/cDNA sequence per admitted
    `mRNA`/`transcript` feature (or per transcript admitted by the selected
    splicing scope).
  - preserves transcript provenance on the derived sequence through synthetic
    local `mRNA` and `exon` features.
  - when coding context is available, also derives a synthetic local `CDS`
    feature and attached protein-translation qualifiers on the derived
    transcript sequence.
  - CDS/protein derivation now resolves from, in order:
    - explicit transcript `cds_ranges_1based`
    - matching source `CDS` features that fit within the transcript exons
    - optional `/codon_start` or `phase`
    - optional `/transl_table`
    - source/transcript/CDS `organism` and `organelle` context
  - translation-table resolution is deterministic:
    - explicit `/transl_table` on CDS, transcript, or source wins
    - plastid/chloroplast-like organelles default to NCBI table `11`
    - vertebrate mitochondrial context (currently human and mouse) defaults to
      NCBI table `2`
    - common invertebrate mitochondrial context (currently fruit-fly and
      nematode matching) defaults to NCBI table `5`
    - yeast mitochondrial context defaults to NCBI table `3`
    - *Escherichia coli* source context defaults to NCBI table `11`
    - unresolved mitochondrial context without explicit `/transl_table` still
      falls back to table `1` and emits an explicit warning
  - the derived transcript still keeps translation qualifiers locally, and
    `DeriveProteinSequences` can then materialize the corresponding peptides as
    first-class protein sequence entries.
- Derived feature qualifiers:
  - derived `mRNA` may now include:
    - `cds_ranges_1based`
    - `protein_length_aa`
    - `translation_table`
    - `translation_table_label`
    - `translation_table_source`
    - `translation_context_organism`
    - `translation_context_organelle`
    - `translation_speed_profile_hint`
    - `translation_speed_profile_source`
    - `translation_speed_reference_species`
    - `derived_protein_translation`
  - derived synthetic `CDS` may now include:
    - `translation`
    - `codon_start`
    - `transl_table`
    - `translation_table_label`
    - `translation_table_source`
    - `protein_length_aa`
    - `terminal_stop_trimmed`
    - `translation_speed_profile_hint`
    - `translation_speed_profile_source`
    - `translation_speed_reference_species`
    - zero or more `translation_warning` qualifiers
- Translation-speed preparation:
  - transcript derivation now records a normalized
    `translation_speed_profile_hint` when the source organism resolves to one
    of the initial target species:
    - `human`
    - `mouse`
    - `yeast`
    - `ecoli`
  - provenance is now explicit through:
    - `translation_speed_profile_source`
    - `translation_speed_reference_species`
  - the bundled `mouse` speed profile now uses a dedicated mouse reference
    column (`Mus musculus domesticus`) instead of the earlier rat proxy
- Output:
  - additive sequence creation through regular `OpResult.created_seq_ids`
  - deterministic messages/warnings about CDS absence, translation-table
    fallback, partial codons, ambiguous codons, or internal stops.

`PlanExonSkippedIsoform` / `MaterializeExonSkippedIsoform` semantics:

- `PlanExonSkippedIsoform` is the phase-1 selection surface.
  - Inputs: `seq_id`, `transcript_feature_id`, typed `criteria[]`, optional
    `plan_id`.
  - Output: `gentle.exon_skip_selection_plan.v1`.
  - Stores the plan under `ProjectState.metadata["exon_skip_selection_plans"]`.
  - Candidate exons expose stable IDs (`exon_1`, `exon_2`, ...), genomic
    coordinates, support/constitutive status, selection sources, matched
    feature IDs, and CDS/frame warnings.
  - Criteria support manual candidate IDs, explicit intervals, current map
    selection intervals, feature-overlap queries, and reserved reasoning-source
    candidate IDs.
- `MaterializeExonSkippedIsoform` is phase 2 and consumes a stored plan.
  - Inputs: `plan_id`, optional `selected_candidate_ids[]`, optional
    `output_prefix`, optional `return_kinds[]`.
    - supported return kinds:
      `genbank|cdna_fasta|amino_acid_sequence|amino_acid_fasta`
  - Output: `gentle.exon_skip_materialization.v1`.
  - Does not re-plan; stale source transcript coordinates are rejected.
  - Rejects empty/all-skipped isoforms.
  - Creates both a retained-exon cDNA/mRNA sequence and a genomic annotation
    product carrying the exon-skipped transcript model.
  - Generated feature qualifiers record the plan ID, source sequence/feature,
    skipped exon candidate IDs, and synthetic origin.
  - When `return_kinds[]` is supplied, `return_payloads[]` carries the
    requested machine-facing handoff text. This is the preferred ClawBio/OpenClaw
    route for saying whether the caller wants the adjusted GenBank entry, the
    retained-exon cDNA FASTA, the amino-acid sequence, or amino-acid FASTA.

`DeriveProteinSequences` semantics:

- Inputs:
  - `seq_id`
  - optional `feature_ids[]`
  - optional splicing `scope`
  - optional `output_prefix`
- Behavior:
  - derives one first-class protein sequence per selected/admitted transcript
  - transcript-native translation is the authoritative product derivation path;
    it does not consult UniProt or other external protein sources
  - `protein` and `peptide` molecule labels are treated as the same first-class
    protein family in shared FASTA/import/export handling
  - uses annotated CDS translation when available
  - if CDS annotation is absent, falls back deterministically to:
    - an inferred ATG-start ORF on the derived transcript
    - otherwise the longest stop-free reading-frame segment
  - emits one full-span local `Protein` feature on the derived peptide with:
    - transcript/source provenance
    - derivation mode (`annotated_cds`, `inferred_orf`, `heuristic_longest_frame`)
    - translation-table context
    - organism/organelle context when available
    - optional `translation_speed_profile_hint`
    - optional `translation_speed_profile_source`
    - optional `translation_speed_reference_species`
- Output:
  - additive sequence creation through regular `OpResult.created_seq_ids`
  - persists one transcript-native protein-derivation report with stable
    `report_id`, created protein sequence ids, and stored `op_id` / `run_id`
    provenance for lineage/reopen flows
  - deterministic warnings when CDS annotation is missing or heuristic
    inference had to be used
  - when `feature_query` is supplied, the engine deterministically admits only
    the matched transcript features before translation and stores that query on
    the persisted report so later protein-gel exports can reopen the same
    selection without guessing

`ReverseTranslateProteinSequence` semantics:

- Inputs:
  - `seq_id` (must resolve to a protein sequence)
  - optional `output_id`
  - optional `speed_profile`:
    - `human`
    - `mouse`
    - `yeast`
    - `ecoli`
  - optional `speed_mark`:
    - `fast`
    - `slow`
  - optional `translation_table`
  - optional `target_anneal_tm_c`
  - optional `anneal_window_bp`
- Behavior:
  - generates one synthetic coding DNA sequence for the selected protein
  - codon choice is deterministic and translation-table-aware
  - translation-table resolution is deterministic in this order:
    - explicit request
    - existing protein feature translation qualifiers
    - organism/organelle context inferred from attached features
  - `speed_mark=fast` biases toward preferred codons for the selected/bundled
    species profile
  - `speed_mark=slow` biases away from the preferred codon when synonymous
    choices exist
  - optional `target_anneal_tm_c` applies a lightweight local suffix Tm
    heuristic over `anneal_window_bp` windows to mildly steer codon choice; it
    is advisory rather than a full sequence optimizer
- Output:
  - additive sequence creation through regular `OpResult.created_seq_ids`
  - persisted reverse-translation report metadata under
    `reverse_translation_reports`:
    - stable `report_id`
    - `op_id` / `run_id`
    - source protein sequence id
    - created coding-sequence id
    - effective translation table / label / source
    - optional organism/organelle context
    - resolved speed profile / source / reference species
    - optional speed mark and annealing-Tm steering inputs
    - codon-choice diagnostics:
      - `preferred_synonymous_choice_count`
      - `alternative_synonymous_choice_count`
      - `fallback_unknown_codon_count`
      - `gc_fraction`
      - `realized_anneal_tm_c`
    - warnings captured during codon selection / qualifier synthesis
  - `OpResult.reverse_translation_report` now carries the same resolved
    provenance in one portable record:
    - created coding-sequence id
    - protein length / coding length
    - translation table / label / source
    - optional organism/organelle context
    - requested/resolved translation-speed profile
    - speed-profile source / reference species
    - optional speed mark
    - optional annealing-target heuristic settings
    - warnings
  - shared shell / direct CLI inspection route:
    - `reverse-translate run PROTEIN_SEQ_ID ...`
    - `reverse-translate list-reports [PROTEIN_SEQ_ID]`
    - `reverse-translate show-report REPORT_ID`
    - `reverse-translate export-report REPORT_ID OUTPUT.json`
  - one full-span synthetic local `CDS` feature with:
    - protein provenance
    - translation table/label
    - translation table source
    - optional translation context organism/organelle qualifiers
    - optional speed profile/source/mark qualifiers
    - optional speed reference-species qualifier
    - optional annealing-target qualifiers
    - zero or more `reverse_translation_warning` qualifiers

`RenderProtocolCartoonSvg` semantics:

- Inputs:
  - `protocol` (built-in ids now include `gibson.two_fragment`,
    `gibson.single_insert_dual_junction`, `pcr.assay.pair`,
    `pcr.assay.pair.no_product`, `pcr.assay.pair.with_tail`, and
    `pcr.assay.qpcr`)
  - `path` (output SVG)
- Behavior:
  - renders a deterministic protocol-cartoon strip through one engine route,
    independent of GUI/CLI entry point.
  - emits canonical conceptual step order for the requested protocol as an
    ordered event-sequence model.
  - template representation baseline is now available in engine internals:
    - schema id: `gentle.protocol_cartoon_template.v1`
    - sparse template rows (event/molecule/feature) are resolved with
      deterministic defaults into render-ready specs.
    - built-in protocol families should now be composed from shared internal
      figure building blocks (feature spans, strand-specific tails, linear
      molecule rows, event rows) rather than ad-hoc per-protocol struct
      literals; this keeps future PCR/Gibson growth on one composition model
  - internal model used by renderer:
    - event -> molecules -> feature fragments
    - molecule topology supports `linear|circular`
    - linear molecules may carry end styles
      (`NotShown`, `Continuation`, `Blunt`, or
      `Sticky { polarity: FivePrime|ThreePrime, nt }`)
    - feature fragments can optionally render different top-strand and
      bottom-strand colors and lengths, plus strand-specific nicks after a
      segment boundary; this is useful for annealed overlaps, exonuclease
      chew-back cartoons with single-stranded tails, and polymerase-filled
      intermediates that still require ligase
  - malformed protocol cartoon specs fail validation and render deterministic
    invalid-spec SVG diagnostics instead of panicking.
- Output:
  - deterministic SVG artifact; operation is non-mutating.

`RenderProtocolCartoonTemplateSvg` semantics:

- Inputs:
  - `template_path` (JSON file path, schema `gentle.protocol_cartoon_template.v1`)
  - `path` (output SVG)
- Behavior:
  - reads template JSON from disk and parses it deterministically.
  - resolves sparse event/molecule/feature rows using deterministic defaults
    (action/caption/topology/end styles/feature length/palette).
  - validates resolved cartoon semantics before rendering.
- Output:
  - deterministic SVG artifact; operation is non-mutating.

`ValidateProtocolCartoonTemplate` semantics:

- Inputs:
  - `template_path` (JSON file path, schema `gentle.protocol_cartoon_template.v1`)
- Behavior:
  - reads and parses template JSON deterministically.
  - resolves sparse defaults and validates resolved cartoon semantics.
  - emits validation diagnostics through operation result messages; no SVG is
    written.
- Output:
  - non-mutating validation result suitable for pre-render checks in CLI/GUI
    flows.

`RenderProtocolCartoonTemplateWithBindingsSvg` semantics:

- Inputs:
  - `template_path` (JSON file path, schema `gentle.protocol_cartoon_template.v1`)
  - `bindings_path` (JSON file path, schema
    `gentle.protocol_cartoon_template_bindings.v1`)
  - `path` (output SVG)
- Behavior:
  - loads template and binding payloads.
  - applies deterministic ID-targeted overrides (defaults, event, molecule,
    feature) and then resolves the bound template.
  - validates resolved semantics before SVG rendering.
- Output:
  - deterministic SVG artifact; operation is non-mutating.

`ExportProtocolCartoonTemplateJson` semantics:

- Inputs:
  - `protocol` (built-in protocol cartoon id, for example `gibson.two_fragment`
    or `pcr.assay.pair` / `pcr.assay.pair.with_tail` / `pcr.assay.qpcr`)
  - `path` (output JSON file)
- Behavior:
  - materializes the canonical built-in template
    (`gentle.protocol_cartoon_template.v1`) for the requested protocol.
  - writes deterministic pretty JSON suitable for user editing/tweaking.
- Output:
  - deterministic JSON artifact; operation is non-mutating.

`ExportProcessRunBundle` semantics:

- Exports a deterministic JSON run bundle artifact (`gentle.process_run_bundle.v1`)
  for reproducibility/audit.
- Inputs:
  - `path` (required): output JSON file
  - `run_id` (optional): when set, only operation-log rows for that `run_id`
    are exported; when omitted, all operation-log rows are exported.
- Payload sections:
  - `inputs`:
    - per-operation extracted input references
      (`sequence_ids`, `container_ids`, `arrangement_ids`, candidate/guide sets,
      genome ids, file inputs)
    - aggregated referenced ids and inferred `root_sequence_ids`
  - `parameter_overrides`:
    - chronological `SetParameter` overrides with `op_id`, `record_index`,
      parameter `name`, and exact JSON `value`
  - `decision_traces`:
    - optional routine-assistant trace rows
      (`gentle.routine_decision_trace.v1`) captured in project metadata and
      exported for routine-selection reproducibility (`trace_id`, selected
      routine/alternatives, disambiguation questions/answers, binding snapshot,
      helper-aware `routine_preference_context`, candidate planning-score
      snapshots, suggested macro templates, ordered `preflight_history`,
      canonical `preflight_snapshot`, execution outcome, export events)
  - `operation_log`:
    - selected immutable operation records (`run_id`, operation payload, result)
  - `outputs`:
    - created/changed sequence ids
    - final sequence summaries for affected ids
    - container/arrangement ids created by selected operations
    - file artifact paths produced by selected operations
  - `parameter_snapshot`:
    - full current engine parameter snapshot at export time.
  - `construct_reasoning`:
    - portable construct-reasoning payload for ClawBio/OpenClaw and other
      agent-facing consumers
    - `seq_ids_considered`: deterministic union of referenced plus
      created/changed sequence ids from the exported run slice
    - `summaries`: compact per-sequence reasoning summaries
      (objective, fact/decision coverage, host/helper/medium context,
      interpreted growth signals, supported selection rules, and warning lines)
    - `graphs`: the selected stored `gentle.construct_reasoning_graph.v1`
      payloads themselves for full offline inspection
- Failure modes:
  - empty `path` => `InvalidInput`
  - unknown filtered `run_id` (no selected rows) => `NotFound`

`ExportLabAssistantInstructions` semantics:

- Exports a deterministic bench handoff
  (`gentle.lab_assistant_instructions.v2`) for non-IT lab assistants from the
  current operation journal and sequence/container state.
- Supported output formats are `markdown`, `odt`, and `docx`. When `format` is
  omitted, GENtle infers the format from `path`; unknown extensions fall back to
  `markdown`.
- `odt` and `docx` reports automatically embed a generated lineage overview
  graphic when SVG-to-PNG rasterization succeeds. If rasterization fails, the
  report is still written and `warning_lines[]` records the reason.
- Inputs:
  - `path` (required): output report path
  - `run_id` (optional): when set, only operation-log rows for that `run_id`
    are summarized; when omitted, all operation-log rows are summarized
  - `title` (optional): heading for the handoff; otherwise inferred from the
    first recognized cloning operation or falls back to "GENtle cloning handoff"
  - `audience` (optional): human-facing audience label, defaulting to
    "Lab assistant"
  - `format` (optional): one of `markdown`, `odt`, or `docx`
- Returned `OpResult.lab_assistant_instructions` payload:
  - `output_format`: realized output format
  - `material_rows[]`: sequence, container, and arrangement IDs with display
    names, source role, length/topology where applicable, members, and notes
  - `step_sections[]`: "Before starting", "Design-derived bench sequence",
    and "After the bench work" sections
  - `embedded_visuals[]`: graphical overview rows embedded in editable document
    formats, including visual id, label, format, source, and pixel dimensions
  - `checkpoint_lines[]`, `safety_lines[]`, and `record_keeping_lines[]`
  - `warning_lines[]` when no recorded operations or no cloning-specific
    operations were available
- The export names design intent and verification checkpoints. It deliberately
  does not invent volumes, temperatures, incubation times, kit choices, or
  biosafety approvals that are not present in GENtle state; local SOPs and
  supervisor-approved conditions remain authoritative.
- Failure modes:
  - empty `path` => `InvalidInput`
  - unknown filtered `run_id` (no selected rows) => `NotFound`

Construct reasoning graph foundation (implemented first slice):

- Shared portable records now exist for:
  - `gentle.construct_objective.v1`
  - `gentle.design_evidence.v1`
  - `gentle.design_fact.v1`
  - `gentle.design_decision_node.v1`
  - `gentle.construct_candidate.v1`
  - `gentle.construct_reasoning_inspection_action.v1`
  - `gentle.construct_reasoning_graph.v1`
  - `gentle.construct_reasoning_store.v1`
  - `gentle.host_profile_catalog.v1`
- Current engine-backed scope:
  - project metadata key: `construct_reasoning`
  - objective upsert/store round-trip
  - construct-candidate schema now also reserves additive
    `protein_to_dna_handoff` detail for ranked DNA handoff suggestions,
    including:
    - strategy (`transcript_native_reuse`, `feature_coding_dna`,
      `reverse_translated_synthetic`)
    - source protein sequence id plus source artifact references
    - optional transcript/projection/Ensembl/feature-query context
    - amino-acid coverage summary
    - preserved / relaxed constraints
    - translation-table / speed-profile / speed-mark summary
    - provenance-oriented score, warnings, and next-step recommendations
  - construct-objective schema now reserves additive host/helper context
    fields:
    - `propagation_host_profile_id`
    - `expression_host_profile_id`
    - `host_route[]`
    - `medium_conditions[]`
    - `helper_profile_id`
    - `required_host_traits[]`
    - `forbidden_host_traits[]`
    - optional `intended_tasks[]` for explicit workflow intent using the same
      `pcr|nanopore_sequencing|read_mapping|cloning_stability|construct_maintenance`
      vocabulary as task severity. A missing field is legacy/unspecified and
      permits conservative positive inference; an explicit empty list means
      none of those tasks applies.
  - design-evidence schema now also reserves additive non-sequence context
    fields:
    - `scope`
    - `host_profile_id`
    - `host_route_step_id`
    - `helper_profile_id`
    - `medium_condition_id`
    - optional typed `repeat_annotation` with source kind, repeat name, class,
      family, and source references. Legacy `repeat_*=` note parsing remains a
      compatibility fallback for older graphs rather than the primary model.
  - deterministic read-only graph build from:
    - construct-objective context such as selected propagation/expression host
      profiles, helper profile, host-route steps, medium conditions, and
      required/forbidden host traits
    - existing sequence facts: restriction sites plus sequence-feature spans
      such as exon/CDS/gene/transcript/UTR/promoter/TFBS/variant when present
  - deterministic hard-rule fact/decision population for:
    - propagation-host context
    - expression-host context
    - host-transition status
    - host-route restriction/methylation review from:
      - explicit route-step trait text such as `hsdR- M+`, `dam+`, `dcm+`,
        `hsdR+`, or `MDRS+`
      - deterministic sequence motif tallies for Dam, Dcm, and EcoK-like
        target sites
    - adapter/linker restriction-capture review from:
      - explicit construct-objective capture plans that record:
        - capture restriction enzyme
        - minimal vs MCS-like adapter style
        - whether blunt insert ends are required
      - optional insert-methylation protection intent
      - optional extra retrieval-site enzymes on longer adapters
      - deterministic internal-site tallies for every named adapter motif
        (capture site plus extra retrieval-site enzymes) on the insert
        sequence
      - review flags when:
        - the capture enzyme cannot be resolved from the active catalog
        - any named adapter motif already occurs on the insert without
          protection
        - all named adapter motifs already occur on the insert, so the adapter
          may contribute little retrieval discrimination
        - insert methylation is requested and therefore still requires
          enzyme-specific methylation-sensitivity review before assuming only
          newly ligated adapter sites will cut
    - growth/condition context from structured medium-condition interpretation
      (for example nutrient omission, antibiotic selection agent, heat shock,
      and temperature signals)
    - helper/MCS context
    - variant-effect context derived from overlap of mapped variant markers
      against promoter/enhancer/TFBS/CDS/exon/UTR/splice evidence already in
      the graph
      - per-variant summaries now keep transcript-level ambiguity explicit via
        `transcript_context_status` and `transcript_effect_summaries[]`
        instead of flattening multi-transcript consequences away
    - variant-assay context that maps the same deterministic overlap rules onto
      first assay-family suggestions such as:
      - promoter/regulatory reporter follow-up
      - allele-paired coding-expression comparison
      - minigene splice follow-up
      - UTR reporter / translation comparison
    - selection/complementation context built from engine-owned
      selection/complementation rules, currently seeded with:
      - the proline-rescue baseline (`proA`/`proB`-style annotated construct
        features plus proline-style medium conditions)
      - helper-backed selectable-marker context from normalized helper-profile
        interpretation (for example helper semantics such as `AmpR` /
        `ampicillin_selection`)
  - active-sequence graph refresh helper that reuses the existing graph/objective
    identity when rebuilding deterministic evidence after sequence changes
  - JSON export helper for one stored graph
  - protein-to-DNA handoff reasoning helper (`BuildProteinToDnaHandoffReasoning`)
    that persists inside the same construct-reasoning graph/store instead of
    creating a separate report family:
    - transcript-native CDS reuse when the selected protein is backed by a
      stored project-native derivation report
    - mapped feature-coding DNA when a stored UniProt projection plus feature
      query can recover coding DNA for the selected protein context
    - reverse-translated synthetic fallback when no stronger DNA provenance is
      available or intentionally selected
    - ranking metadata now surfaces in graph summaries through:
      - `contains_protein_to_dna_handoff`
      - `protein_to_dna_handoff_candidate_count`
      - `protein_to_dna_source_protein_seq_ids[]`
  - host-profile catalog loading/list projection from the shared starter JSON
    catalog (`assets/host_profiles.json`) with filter matching across ids,
    aliases, genotype/phenotype tags, notes, and source notes
  - graph-backed annotation-candidate layer:
    - `gentle.annotation_candidate.v1`
    - one portable record per annotation-grade/supporting sequence interval
    - carries support links from one visible span back to the facts/decisions
      it contributes to
    - editable `draft|accepted|rejected|locked` status now persists on the
      candidate itself so reviewed annotation state survives graph refresh
  - graph-backed annotation-summary layer:
    - `gentle.annotation_candidate_summary.v1`
    - one portable collapsed summary per overlapping annotation-candidate
      cluster/signature on a sequence
    - carries grouped annotation ids, transcript-context/effect-tag summaries,
      support labels, and review-status summary so genomic views can stay calm
      without hiding the reasoning lineage
    - intended as the shared summary surface for future similarity-analysis
      annotations as well, instead of creating a second adapter-local collapse
      path
  - annotation-candidate write-back report:
    - `gentle.annotation_candidate_writeback.v1`
    - reports whether one accepted/locked generated candidate created a new
      ordinary sequence feature or was already backed by one
  - graph-level `inspection_actions[]`:
    - portable recommended inspection actions currently used for
      repeat/similarity dotplot handoffs
    - deterministic `action_id` derived from graph/source/mode/focus
    - dotplot mode (`self_forward` or `self_reverse_complement`) chosen from
      structured evidence roles, tags, and orientation signals rather than
      display labels
    - `rationale` carries the source fact/annotation/summary reason for why
      this action is recommended, so adapters do not need to reconstruct that
      explanation from GUI labels
    - `driving_evidence_ids[]` names the specific evidence rows behind the
      action, while source fact/annotation/summary id lists let adapters attach
      the same action to the right inspector row
    - `repeat_family_provenances[]` preserves every overlapping/nested curated
      repeat family with typed name/class/family, evidence ids, source refs,
      confidence, and agreement strength (`provenance_only`,
      `curated_annotation`, `class`, or `family`). The singular
      `repeat_family_provenance` remains a compatibility projection of the
      first row. Standalone Alu-like heuristic rows remain soft hypotheses
      without curated provenance.
    - actions are persisted reasoning-snapshot outputs. Reading or normalizing
      a graph does not silently regenerate non-empty `inspection_actions[]`;
      refresh is the explicit operation that recomputes them.
  - graph-level optional `input_fingerprint`:
    - SHA-256 identities for the complete source sequence/feature snapshot and
      normalized objective plus an explicit reasoning rule-set version
    - optional `source_artifact_kind`, `source_artifact_id`, and
      `source_artifact_sha256` bind a graph to the exact computational report
      it explains; primer-selection graphs currently use
      `source_artifact_kind = primer_design_report`
    - live readers report `current`, `stale`, or `unknown` freshness with
      reasons; older graphs without a fingerprint remain readable as `unknown`
    - changed sequence/features, changed objective, changed rule version, or a
      changed/missing recognized source artifact makes a fingerprinted graph
      stale until it is explicitly refreshed
  - fact-level `task_severities[]`:
    - compact rule-based task interpretation for repeat/similarity facts
    - each row separates intrinsic evidence concern from objective priority:
      `base_severity`/`base_score` describe the evidence itself, while
      `severity`/`score` are the effective objective-aware result;
      `objective_adjustment` records the transparent signed difference
    - `applicability` is `unknown|applicable|not_applicable`, and
      `applicability_basis` is
      `unspecified|explicit_objective_tasks|legacy_objective_inference`
    - initial task vocabulary is `pcr`, `nanopore_sequencing`,
      `read_mapping`, `cloning_stability`, and `construct_maintenance`
    - initial severity vocabulary is `none`, `low`, `medium`, `high`
    - engine-generated scores are transparent 0.0-1.0 rule-derived scalars:
      `0.0` maps to `none`, `(0.0, 0.25)` to `low`, `[0.25, 0.60)` to
      `medium`, and `>= 0.60` to `high`; objective adjustments are explicit and
      the row rationale says why. Explicitly non-applicable tasks preserve
      their intrinsic score while setting effective priority to zero.
    - legacy free-text inference only establishes positive applicability and
      never infers `not_applicable`; explicit `intended_tasks[]` is authoritative
    - task severity layers on top of the generic fact and evidence rows; it
      does not create additional annotation candidates or map overlays
- Current GUI-backed scope:
  - sequence-window `Reasoning` display toggle
  - read-only linear DNA-map overlay that auto-refreshes from the engine-owned
    graph and now paints only annotation-grade/supporting sequence spans
    directly on-sequence:
    - generated or inferred annotation-like spans (for example promoter windows)
    - cDNA-confirmed exon/splice/CDS/transcript spans
    - sequence intervals explicitly referenced by stored facts or decisions
    - raw low-level evidence such as unrestricted restriction-site spam stays
      in the graph/inspector and is no longer drawn by default just because it
      exists as `sequence_span` evidence
  - GUI-side hover/selection inspection for one evidence span at a time
    - selected spans now surface their annotation-candidate source kind,
      supporting fact/decision labels, and collapsed annotation-summary context
      when available
  - side-panel construct-reasoning inspector section for non-sequence facts and
    decision steps (host/helper/host-route restriction-methylation/medium/
    growth/selection context) without pretending they are coordinate spans
  - side-panel `Annotation summaries` section for graph-backed collapsed
    annotation rows that keep overlapping candidates readable on long loci
  - side-panel `Annotation candidates` section for graph-backed candidate
    annotations that now support shared-engine accept/reject/draft curation
    instead of only read-only inspection
    - accepted or locked generated candidates now also expose a shared-engine
      `Write Back` action that materializes them as ordinary sequence features
  - shared shell/CLI route:
    - `construct-reasoning list-inspection-actions GRAPH_ID [--fact-id ID] [--annotation-id ID] [--candidate-id ID] [--evidence-id ID] [--seq-id ID] [--action-kind KIND] [--summary-id ID]`
    - lists the same portable `inspection_actions[]` objects that the GUI
      attaches to fact, annotation, evidence, sequence, and summary rows
    - `construct-reasoning run-inspection-action GRAPH_ID ACTION_ID [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID] [--render-svg OUTPUT.svg]`
    - resolves one `action_id`, computes its dotplot through `ComputeDotplot`,
      returns the resolved `compute_parameters`, and optionally writes SVG
      evidence through `RenderDotplotSvg`
    - graph show/list/action responses include live snapshot freshness;
      execution rejects a fingerprinted stale graph and asks the caller to
      refresh instead of applying old action coordinates to new inputs
    - `construct-reasoning set-annotation-status GRAPH_ID ANNOTATION_ID draft|accepted|rejected|locked`
    - updates the persisted graph in place and returns the updated candidate
      plus the same compact summary block exposed by
      `construct-reasoning show-graph`
    - `construct-reasoning write-annotation GRAPH_ID ANNOTATION_ID`
    - writes back one accepted or locked generated candidate as a feature,
      refreshes the graph, and returns the `gentle.annotation_candidate_writeback.v1`
      report alongside the refreshed candidate and summary block
  - Planning-window `Host Profile Browser` backed by the same shared catalog so
    host/strain traits can be inspected without reparsing raw JSON
  - GUI-only role/class visibility filters layered on top of the same shared
    engine-owned graph payload (no adapter-local graph recomputation)
  - ClawBio/OpenClaw-facing run-bundle export integration:
    - deterministic per-sequence summary rows for concise agent consumption
      including additive variant effect tags and suggested assay-family ids
    - embedded stored reasoning graphs for full offline inspection/replay
  - planning/routine handoff now also consumes the same graph:
    - `routine_preference_context` records can carry
      `construct_reasoning_seq_id`,
      `construct_strategy_derived_preferred_routine_families`,
      `variant_effect_tags`,
      `variant_suggested_assay_ids`, and
      `variant_derived_preferred_routine_families`
    - `routines list` / `routines compare` planning estimates and GUI Routine
      Assistant traces can therefore explain when adapter/linker capture or
      variant-derived assay context boosted one routine family over another
  - sequence-complexity / similarity / mobile-element reasoning baseline:
    - generated `repeat_region` spans for:
      - low-complexity windows
      - homopolymer / tandem-repeat segments
      - direct-repeat clusters
      - inverted-repeat clusters
    - generated `mobile_element` spans for cautious `Alu-like SINE candidate`
      calls
    - fact/decision layer now includes:
      - `sequence_complexity_context`
      - `repeat_architecture_context`
      - `mobile_element_context`
      - `similarity_operational_risk_context`
      - `pcr_operational_risk_context`
      - `nanopore_operational_risk_context`
      - `mapping_operational_risk_context`
      - `cloning_stability_context`
    - materialized RepeatMasker/UCSC `rmsk`-style repeat annotations that
      overlap internal repeat/mobile-element signals are kept as separate
      reliable evidence rows and summarized into `curated_repeat_support[]`,
      including repeat name/class/family and supporting evidence ids
    - repeat/similarity dotplot handoff buttons consume graph-level
      `inspection_actions[]`; adapters no longer need to infer those actions
      from labels or GUI-only context matching
- Current evidence-class rules:
  - restriction sites => `hard_fact`
  - dbSNP / VCF-generated variant markers => `hard_fact`
  - exon/splice annotations with explicit cDNA-style qualifier hints =>
    `hard_fact`
  - imported/derived sequence annotations => `reliable_annotation`
  - TFBS-style annotations => `soft_hypothesis`
  - materialized RepeatMasker/UCSC `rmsk`-style repeat/mobile-element
    annotations => `reliable_annotation`
  - generic mobile-element annotations without curated repeat-family qualifiers
    => `soft_hypothesis`
- Current evidence-scope behavior:
  - graph builds now emit both:
    - `sequence_span` evidence for mapped restriction/annotation features
    - non-sequence construct-objective context evidence
      (`host_profile`, `host_transition`, `medium_condition`,
      `helper_profile`, `whole_construct`) when the objective carries those
      fields
  - GUI DNA overlays intentionally keep rendering only `sequence_span`
    evidence; non-sequence evidence stays in the portable graph payload rather
    than being faked as coordinate spans
- Not in this slice yet:
  - construct-candidate ranking
  - curated host/helper profile catalog loading and biological compatibility
    scoring against those catalogs
  - editable reasoning/decision GUI surfaces beyond the current read-only span
    overlay plus inspector summary

Protocol-cartoon family growth direction (planned):

- Gibson specialist work now validates the abstraction-first protocol-cartoon
  strategy:
  - one protocol family should be expressed as a canonical
    `gentle.protocol_cartoon_template.v1` template plus deterministic bindings
  - the renderer should grow a shared collection of reusable figure building
    blocks that protocol families compose, instead of embedding protocol-
    specific drawing code in each built-in cartoon
  - protocol growth in count/shape (for example multi-fragment Gibson) should
    prefer repeated events, repeated molecules, and binding-level overrides
    rather than renderer-specific special cases
- Implemented baseline:
  - built-in Gibson cartoons now compose from shared internal building blocks
    for duplex spans, strand-specific tails, linear molecule rows, and event
    rows
  - this is intentionally still mechanism-first: Gibson cartoons describe
    fragment flow and achieved homology relationships, not full primer objects
    or low-level PCR parameter details
- PCR-assay cartoons should follow the same rule. The shared renderer should
  remain chemistry-agnostic and continue to render only ordered events,
  molecules, and features; PCR-specific meaning belongs in template structure
  and bindings, not in new PCR-only drawing primitives.
- PCR cartoon purpose:
  - explain assay intent and artifact flow, not every thermocycler sub-step
  - stay aligned with engine-owned operations, reports, and lineage-visible
    artifacts
  - keep lower-level primer sequences/thermocycler details in adjacent textual
    reports or inspectors rather than the cartoon itself
- Canonical PCR assay scene vocabulary should stay stable across modalities:
  - source template/context event
  - target/ROI event (selected span, feature-derived span, or queued region)
  - assay setup event (forward/reverse pair and optional probe or staged
    inner/outer sets may be named, but do not need literal primer glyphs)
  - amplification event
  - product/artifact event (amplicon, extracted copy, report, export, or
    explicit no-accepted-pairs outcome)
- Implemented PCR baseline in the `pcr.assay.*` family:
  - `pcr.assay.pair`: base strip with one selected ROI, one assay-setup lane,
    one amplification step, one amplicon/report outcome, and explicit forward/
    reverse primer glyphs with 5'/3' orientation
  - `pcr.assay.pair.no_product`: same family with an explicit report-only
    terminal state when no accepted primer pair yields a product
  - `pcr.assay.pair.with_tail`: insertion-first strip with requested extension
    sequences + insertion anchors, anchor-adjacent primer windows, and carried-
    through inserted terminal tails in the final amplicon
  - `pcr.oe.substitution`: six-step overlap-extension substitution strip with
    primer set `a`..`f`, first-step product haplotypes (`AB`/`CD`/`EF`),
    strand-specific anneal-gap geometry, and polymerase fill
- Implemented qPCR baseline in the same `pcr.assay.*` family:
  - `pcr.assay.qpcr`: same base strip enriched with one internal probe window
    plus explicit forward/reverse primer glyphs, a retained probe marker, and
    one quantitative readout terminal state
- Planned PCR modality adaptation should continue through the same
  `pcr.assay.*` protocol-cartoon family:
  - nested PCR: same family with two amplification stages (outer -> inner)
    instead of one, reusing the same event vocabulary
  - inverse PCR: same family with circular-template bindings and outward-facing
    primer semantics
  - batch/multiplex/tiling: repeated assay groups or repeated output lanes in
    bindings, not new renderer semantics per assay count
  - empty/failure outcomes: report/artifact nodes can render without product
    nodes when no accepted assay is produced
- Recommended rollout order:
  - extend the shipped PCR/qPCR baseline to queued batch PCR without changing
    renderer semantics
  - add nested, inverse, long-range, and multiplex variants as further
    template/binding expansions
- Naming/design rule:
  - do not introduce one built-in protocol id per assay count or minor UI view
  - prefer one stable protocol family with bindings that carry assay modality,
    stage count, molecule presence, and repeated-lane structure
  - keep generated explanatory strips exportable through the existing
    `protocol-cartoon ...` routes

RNA secondary-structure semantics:

- Inspection API:
  - `GentleEngine::inspect_rna_structure(seq_id)`
  - Runs `RNAfold --noPS` on the normalized RNA sequence and returns a structured text report (`stdout`/`stderr`, dot-bracket `structure`, optional `mfe_kcal_per_mol`, and command metadata).
- Export operation:
  - `RenderRnaStructureSvg { seq_id, path }`
  - Runs `RNAfold --noPS`, writes a temporary two-line `sequence` + dot-bracket input file, then runs `rnapkin -o <path> <input>`.
  - Expects SVG/PNG output according to the extension accepted by `rnapkin`; current operation naming remains SVG-oriented for existing callers.
- Input constraints:
  - accepted only for single-stranded RNA (`molecule_type` `RNA` or `ssRNA`)
  - empty sequence is rejected
- Runtime dependency:
  - external `RNAfold` and `rnapkin` executables are required for rendering
  - executable path resolution order for folding:
    1. env var `GENTLE_RNAFOLD_BIN`
    2. fallback executable name `RNAfold` in `PATH`
  - executable path resolution order for drawing:
    1. env var `GENTLE_RNAPKIN_BIN`
    2. fallback executable name `rnapkin` in `PATH`

DNA ladder catalog semantics:

- Inspection API:
  - `GentleEngine::inspect_dna_ladders(name_filter?)`
  - Returns structured ladder metadata:
    - `schema` (`gentle.dna_ladders.v1`)
    - `ladder_count`
    - `ladders[]` (`name`, `loading_hint`, `min_bp`, `max_bp`, `band_count`, `bands`)
- Export operation:
  - `ExportDnaLadders { path, name_filter? }`
  - Writes the same structured payload to JSON at `path`.
  - Optional `name_filter` applies case-insensitive name matching before export.

RNA ladder catalog semantics:

- Inspection API:
  - `GentleEngine::inspect_rna_ladders(name_filter?)`
  - Returns structured ladder metadata:
    - `schema` (`gentle.rna_ladders.v1`)
    - `ladder_count`
    - `ladders[]` (`name`, `loading_hint`, `min_nt`, `max_nt`, `band_count`, `bands`)
- Export operation:
  - `ExportRnaLadders { path, name_filter? }`
  - Writes the same structured payload to JSON at `path`.
  - Optional `name_filter` applies case-insensitive name matching before export.

Historical screenshot artifact contract (still disabled):

- Guardrail:
  - command is currently rejected by security policy even when
    `--allow-screenshots` is provided.
- Command surface:
  - direct CLI: `gentle_cli screenshot-window OUTPUT.png`
  - shared shell (CLI and GUI shell panel): `screenshot-window OUTPUT.png`
- Scope and safety:
  - captures only the active/topmost GENtle window
  - window lookup is native AppKit in-process (no AppleScript automation path)
  - command is primarily intended for GUI shell contexts with an active window
  - rejects full-desktop capture and non-GENtle targets
  - rejects request if no eligible active GENtle window is available
  - its historical backend was macOS `screencapture`; it is not reused by the
    GUI Agent-help flow, whose optional native path uses ScreenCaptureKit
- Output:
  - writes an image file at caller-provided `OUTPUT` path (custom filename
    supported)
  - recommended default image format is inferred from extension (e.g. `.png`)
- Result payload shape:

```json
{
  "schema": "gentle.screenshot.v1",
  "path": "docs/images/gui-main.png",
  "window_title": "GENtle - pGEX-3X",
  "captured_at_unix_ms": 1768860000000,
  "pixel_width": 1680,
  "pixel_height": 1020,
  "backend": "macos.screencapture"
}
```

Operation progress/cancellation semantics:

- `apply_with_progress` and workflow progress callbacks receive
  `OperationProgress` updates.
- Callback return value:
  - `true`: continue
  - `false`: request cancellation
- Current event families:
  - `PrimerDesign`
  - `Tfbs`
  - `GenomePrepare`
  - `GenomeTrackImport`
  - `DbSnpFetch`
  - `ReadAcquisition`
  - `RnaReadInterpret`
- Current cancellation support:
  - internal pair-PCR/qPCR design now emits staged `PrimerDesign` snapshots
    (`candidate_enumeration`, `pair_evaluation`, `probe_evaluation`,
    `complete`) with candidate/evaluation counters, but does not yet expose
    cooperative cancellation
  - genome preparation supports cooperative cancellation plus optional
    `timeout_seconds` timeboxing and reports deterministic cancellation/timeout
    outcomes.
  - read acquisition supports cooperative cancellation through the progress
    callback and explicit cancel-marker requests (`ReadAcquireCancel` /
    `reads acquire cancel`); while an external SRA Toolkit process is running,
    GENtle emits `ReadAcquisition` snapshots, updates the activity JSON, checks
    available disk against `min_free_gb`, and terminates the child process group
    on cancellation or disk-threshold failure.
  - shared reference/helper prepare is now duplicate-safe at the semantic
    target level:
    - one target may actively prepare at a time
    - duplicate callers reuse the active run and observe `lifecycle_status =
      "running"` instead of starting redundant work
    - stale heartbeat markers are converted into `lifecycle_status = "stale"`
      and can then be retried safely
    - on Unix, a `running` marker with a dead `owner_pid` is also converted to
      `stale`; other platforms retain the existing lock/heartbeat behavior
  - genome-track imports support cooperative cancellation and return partial
    import warnings.
  - dbSNP fetch currently emits staged progress events (`validate_input`,
    `inspect_prepared_genome`, `contact_server`, `wait_response`,
    `parse_response`, `resolve_placement`, `extract_region`,
    `attach_variant_marker`) but does not yet expose cooperative cancellation.
  - primer/qPCR design now emits machine-readable progress snapshots for
    candidate-generation and search stages:
    - `forward_candidates`
    - `reverse_candidates`
    - `pair_search`
    - `pair_search_complete`
    - qPCR-only: `probe_candidates`, `assay_search`,
      `assay_search_complete`
    - Primer3-backed routes emit `primer3_run`, zero or more live
      `primer3_search` rows, and optional `fallback_to_internal`
    Each `PrimerDesignProgress` row carries backend request/use,
    ROI bounds, candidate counts, pair-evaluation limits, accepted
    pair/assay counts, and a `done` flag so agent/CLI consumers can detect
    stalled or overly broad searches without parsing GUI state.
    - a `primer3_search` row additionally carries optional
      `primer3_progress = {record, completed, bound}`. Within one input record,
      `completed` is nondecreasing and `bound` is nonincreasing according to
      Primer3's bounded-progress contract. `completed / bound` is suitable for
      a progress indicator, but the counters measure candidate evaluations and
      must not be described as elapsed time or an ETA.
    - Search records created by splitting a broad interval overlap by at most
      `max_primer_length - 1` bases so a primer crossing the split is not lost.
      Consequently, plan-level candidate and search-space totals are sums of
      work performed per record, not counts of unique primer placements; an
      overlap may be evaluated in both records. Removing that repeated work
      requires an algorithmic de-duplication design, not an accounting-only
      subtraction.
    - GENtle detects support through `primer3_core --help` and requests this
      stream only when `--progress` is advertised. On Unix it may send
      `SIGUSR1` after an already-confirmed progress stream becomes quiet,
      prompting another bounded status row. Help-probed legacy binaries run
      directly without the flag; if advertised support is rejected at runtime,
      GENtle retries once in legacy mode. Legacy runs retain coarse stages but
      cannot supply `primer3_progress` counters.
    - returning `false` from the operation progress callback cancels a live
      Primer3 child. In `auto` mode, cancellation is terminal and does not
      trigger an internal-backend fallback.
  - RNA-read interpretation uses cooperative callback checks while emitting
    periodic progress snapshots (including seed-confirmation histogram bins).

`Pcr` semantics (current):

- Exact primer matching on linear templates.
- Enumerates all valid amplicons formed by forward-primer matches and downstream
  reverse-primer binding matches.
- `unique = true` requires exactly one amplicon; otherwise fails.
- `output_id` may only be used when exactly one amplicon is produced.

`PcrAdvanced` semantics:

- Primer spec fields:
  - `sequence` (full primer, 5'->3')
  - `anneal_len` (3' suffix length used for template binding)
  - `max_mismatches` (allowed mismatches within anneal part)
  - `require_3prime_exact_bases` (hard exact-match requirement at primer 3' end)
  - `library_mode` (`Enumerate` or `Sample`) for degenerate/IUPAC primers
  - `max_variants` cap for primer-library expansion
  - `sample_seed` deterministic seed when `library_mode = Sample`
- Supports 5' tails and mismatch-mediated mutagenesis.
- Supports degenerate/randomized synthetic primers via IUPAC codes.
- Product is constructed from:
  - full forward primer sequence
  - template interior between forward and reverse anneal windows
  - reverse-complement of full reverse primer sequence

`PcrMutagenesis` semantics:

- Builds on `PcrAdvanced` primer behavior.
- Accepts explicit SNP intents:
  - `zero_based_position`
  - `reference`
  - `alternate`
- Validates reference bases against the template.
- Filters amplicons to those that introduce requested SNPs.
- `require_all_mutations` (default `true`) controls whether all or at least one
  mutation must be introduced.

`DesignTerminalExonRtPrimerPool` contract (implemented baseline):

- Purpose:
  - design one ordered pool of sequence-specific reverse-transcription oligos
    that share a fixed 5-prime adapter
  - resolve each target through GENtle's annotated Splicing Expert and mature-
    transcript model rather than through adapter-local sequence slicing
  - search exact-length candidate windows near the transcript-oriented start
    of each selected terminal exon
  - persist sufficient geometry, sequence, ranking, interaction, and operation
    provenance for GUI, CLI, MCP, JavaScript, Lua, and Python consumers
- Operation payload:

```json
{
  "DesignTerminalExonRtPrimerPool": {
    "request": {
      "schema": "gentle.terminal_exon_rt_primer_pool_request.v1",
      "fixed_adapter_5prime": "ACTTGCCTGTCGCTCTATCTTC",
      "variable_length_bp": 22,
      "terminal_exon_search_window_bp": 250,
      "max_candidates_per_target": 10,
      "targets": [
        {
          "seq_id": "tp73_locus",
          "source_feature_id": 0,
          "transcript_id": "NM_005427.4",
          "label": "TP73"
        },
        {
          "seq_id": "patz1_locus",
          "source_feature_id": 0,
          "transcript_id": "NM_032052.4",
          "label": "PATZ1"
        }
      ],
      "report_id": "tp73_patz1_terminal_exon_rt_pool"
    }
  }
}
```

The example identifiers illustrate the shape of the request. Callers must use
the exact sequence, feature, and transcript identifiers present in the loaded
GENtle project.

- Request semantics:
  - `targets[]` order identifies the requested pool members and is preserved
    in the report. GENtle first retains a bounded candidate list per target,
    then chooses one candidate per target jointly so an avoidable interaction
    with an earlier target does not become fixed merely because of caller order
  - `source_feature_id` identifies a member of the intended Splicing Expert
    group. When that group contains more than one transcript,
    `transcript_id` is required; GENtle does not guess an isoform
  - `variable_length_bp` is exact. Each candidate must fit wholly inside the
    terminal exon
  - `terminal_exon_search_window_bp` is measured from the terminal exon's
    transcript-oriented start, not from the numerically lower genomic
    coordinate. Reverse-strand transcripts therefore use the same biological
    convention without special adapter behavior
  - `fixed_adapter_5prime` is normalized and must contain only `A`, `C`, `G`,
    and `T`. The full oligo is exactly `fixed_adapter_5prime` followed by the
    selected sequence-specific segment
- Candidate/report semantics:
  - `target_segment_5_to_3` is the selected mature-transcript segment;
    `variable_primer_5_to_3` is its reverse complement; and
    `full_oligo_5_to_3` is the fixed adapter plus that variable primer
  - transcript and source intervals are 0-based, half-open. Source intervals
    remain increasing coordinates even when the transcript is on strand `-`
  - per-target candidate retention uses a deterministic lower-is-better
    lexicographic rank over:
    variable 3-prime complementarity to the adapter, adapter/variable maximum
    complementarity, complete-oligo self 3-prime complementarity, maximum
    prior-pool 3-prime complementarity, complete-oligo self complementarity,
    maximum prior-pool complementarity, distance from the terminal-exon start,
    and finally the variable sequence. The final pool is selected jointly from
    those retained lists by a deterministic bounded search. It protects
    complete-oligo and variable-segment self-complementarity and homopolymer
    quality first, then minimizes the worst complete-oligo and variable-only
    mutual 3-prime runs, general mutual complementary runs, summed complete-
    oligo 3-prime runs, candidate ranks, terminal-exon distance, and ordered
    variable sequences. `pool_selection_policy` and
    `pool_selection_states_evaluated` make that calculation inspectable
  - `variable_tm_c` is descriptive only. Tm does not participate in this
    ranking and the report states that policy explicitly
  - `selected_pool_interactions[]` records variable-only and complete-oligo
    pairwise complementarity for every selected pair
  - schemas are `gentle.terminal_exon_rt_primer_pool_request.v1` and
    `gentle.terminal_exon_rt_primer_pool.v1`; reports are persisted in the
    existing primer-design report store and participate in normal operation
    undo/history
  - by default, no whole-transcriptome or whole-genome off-target screen is
    implied. When the request includes `genomic_specificity`, GENtle executes
    exhaustive `blastn-short` searches for every selected variable segment
    against the named prepared `genomic_dna` index. The additive
    `genomic_specificity` report binds the index kind, assembly/release,
    masking, prefix, sequence/base counts, BLAST database version, inspection
    tool, content fingerprint, every perfect full-length hit, and the raw HSP
    count. GENtle separately compares that exact match with the selected
    candidate's genome-anchored source interval. `passed` means every selected
    segment has exactly one perfect full-length genomic match and that match is
    the intended anchored locus; missing anchors produce `review_required`,
    while a coordinate mismatch fails readiness. This is genomic uniqueness
    evidence, not experimental validation and not whole-transcriptome
    specificity

`ExportTerminalExonRtPrimerPoolReport` writes one persisted
`gentle.terminal_exon_rt_primer_pool.v1` record to the requested JSON path.

`DesignPrimerPairs` contract (implemented baseline):

- Purpose:
  - propose ranked forward/reverse primer pairs for one linear template under
    explicit constraints
  - provide deterministic, machine-readable reports that can be consumed by
    GUI/CLI/scripting/agents
  - bridge PCR primer roles from concept to explicit retrieval/design
    constraints and, for accepted candidates, to in-silico oligo
    instantiations
- Primer lifecycle semantics:
  - choosing a PCR creates primer roles before concrete sequences are known
  - choosing a region/amplicon turns those roles into constraints over a
    substance to design, retrieve, order, or select from local stock
  - choosing a candidate primer satisfies the constraint with a concrete oligo
    design, but does not imply wet-lab availability
  - persisted primer reports and materialized primer sequence artifacts are
    therefore design/provenance records; future inventory/procurement records
    should link to them rather than replacing them
  - multiple chosen primer/probe instantiations may later be collected into one
    reviewed oligo order form, analogous to an assembly grouping its
    participating parts; the order form groups line items for procurement but
    does not itself prove that any physical oligo has arrived
- Operation payload:

```json
{
  "DesignPrimerPairs": {
    "template": "seq_id",
    "roi_start_0based": 1000,
    "roi_end_0based": 1600,
    "forward": {
      "min_length": 20,
      "max_length": 30,
      "location_0based": null,
      "start_0based": null,
      "end_0based": null,
      "min_tm_c": 55.0,
      "max_tm_c": 68.0,
      "min_gc_fraction": 0.35,
      "max_gc_fraction": 0.70,
      "max_anneal_hits": 1,
      "non_annealing_5prime_tail": null,
      "fixed_5prime": null,
      "fixed_3prime": null,
      "required_motifs": [],
      "forbidden_motifs": [],
      "locked_positions": []
    },
    "reverse": {
      "min_length": 20,
      "max_length": 30,
      "location_0based": null,
      "start_0based": null,
      "end_0based": null,
      "min_tm_c": 55.0,
      "max_tm_c": 68.0,
      "min_gc_fraction": 0.35,
      "max_gc_fraction": 0.70,
      "max_anneal_hits": 1,
      "non_annealing_5prime_tail": null,
      "fixed_5prime": null,
      "fixed_3prime": null,
      "required_motifs": [],
      "forbidden_motifs": [],
      "locked_positions": []
    },
    "pair_constraints": {
      "require_roi_flanking": false,
      "required_amplicon_motifs": [],
      "forbidden_amplicon_motifs": [],
      "fixed_amplicon_start_0based": null,
      "fixed_amplicon_end_0based_exclusive": null,
      "rejected_near_miss_limit": null
    },
    "min_amplicon_bp": 120,
    "max_amplicon_bp": 1200,
    "max_tm_delta_c": 2.0,
    "max_pairs": 200,
    "report_id": "tp73_roi_primers_v1"
  }
}
```

- `max_tm_delta_c`, `max_pairs`, `report_id`, and `pair_constraints` are optional in current
  implementation:
  - `max_tm_delta_c` default: `2.0`
  - `max_pairs` default: `200`
  - `report_id` default: auto-generated deterministic-safe id stem
  - `pair_constraints` default:
    `{"require_roi_flanking":false,"required_amplicon_motifs":[],"forbidden_amplicon_motifs":[],"fixed_amplicon_start_0based":null,"fixed_amplicon_end_0based_exclusive":null,"rejected_near_miss_limit":null}`
  - `pair_constraints.rejected_near_miss_limit` retains a bounded,
    deterministic subset of evaluated rejections:
    - omitted/`null`: `20`
    - `0`: disabled
    - maximum: `100`
    - `DesignPrimerPairs` and `DesignInsertionPrimerPairs` retain pair-level
      rows
    - `DesignQpcrAssays` retains evaluated pair/probe assay-level rows; probe
      candidate-generation failures remain aggregate-only
- Side constraints (`forward`, `reverse`, and qPCR `probe`) accept optional
  sequence-level filters:
  - `non_annealing_5prime_tail` (added to the final oligo but excluded from
    anneal Tm/GC/hit calculations)
  - `fixed_5prime`, `fixed_3prime`
  - `required_motifs[]`, `forbidden_motifs[]`
  - `locked_positions[]` entries (`offset_0based`, single IUPAC `base`)
  - an ordinary primer can be fixed exactly while GENtle designs its partner
    by setting that side's `min_length` and `max_length` to the complete
    canonical primer length and `fixed_5prime` to the complete 5'-to-3'
    sequence. The internal backend enumerates only that sequence; the Primer3
    backend additionally emits `SEQUENCE_PRIMER` or
    `SEQUENCE_PRIMER_REVCOMP`. Shorter `fixed_5prime` values remain prefix
    constraints and are not misrepresented as exact supplied primers.
- Built-in primer-ranking heuristics (internal and Primer3 pair-ranking stage):
  - preferred primer length window: `20..30 bp` (outside window is penalized)
  - 3' GC clamp preference (`G/C` at terminal 3' base)
  - secondary-structure risk penalty (homopolymer and self-complementary runs)
  - primer-dimer risk penalty (global and 3'-anchored complementary runs)
  - both backends rank retained pairs with the same GENtle additive model
    `gentle_primer_pair_rank_v1`; Primer3 proposes candidates but does not
    replace the report's GENtle score

- Report schema:
  - `gentle.primer_design_report.v1`
  - deterministic ordering by score then tie-break fields
  - backend metadata block:
    - `backend.requested` (`auto|internal|primer3`)
    - `backend.used` (`internal|primer3`)
    - optional `backend.fallback_reason`
    - optional `backend.primer3_executable`
    - optional `backend.primer3_version`
  - each pair includes:
    - forward/reverse primer sequence and genomic binding window
    - per-primer diagnostics:
      - `length_bp`
      - `anneal_length_bp`
      - `non_annealing_5prime_tail_bp`
      - `three_prime_base`
      - `three_prime_gc_clamp`
      - `longest_homopolymer_run_bp`
      - `self_complementary_run_bp`
    - estimated `tm_c` and `gc_fraction` for annealing segment only
    - anneal-hit counts per side
    - amplicon start/end/length
    - pair dimer diagnostics:
      - `primer_pair_complementary_run_bp`
      - `primer_pair_3prime_complementary_run_bp`
    - rule-pass flags and aggregate score
    - `score_terms[]`, whose `raw_value * weight = contribution` rows sum to
      the existing score; stable terms cover baseline, Tm delta, amplicon
      length fit, extra anneal hits, primer-length preference, homopolymer and
      self-complementarity excess, pair/global and pair/3'-end
      complementarity excess, and 3' GC-clamp balance
  - report-level score interpretation:
    - `score_decomposition_status`, `score_decomposition_reason`
    - `score_model = gentle_primer_pair_rank_v1`
    - `score_direction = higher_is_better`
    - status is `pass` when retained pairs carry exact terms and `not_run`
      when no pair was retained
    - no residual or backend-specific contribution is fabricated
  - optional rejection summary buckets (for explainability):
    - out-of-window
    - GC/Tm out of bounds
    - non-unique anneal
    - primer sequence-constraint failure
    - pair constraint failure
    - amplicon-size or ROI-coverage failure
    - pair-evaluation-limit skipped count
  - bounded selection provenance:
    - `rejected_near_misses[]` contains only evaluated pair-shaped candidates,
      ordered by score (scored before unscored) and stable coordinates/sequence
      tie-breaks
    - each row carries the shared `PrimerDesignRejectionReason` census
      vocabulary in `reasons[]`, exact `failed_checks[]`, coordinates,
      sequences, and score when available
    - `near_miss_capture` records status, scope, requested/effective limit,
      eligible/retained/omitted rows, and deterministic comparison work
    - internal capture is `incomplete` if the pair-evaluation ceiling was
      reached; single-primer rejection buckets remain aggregate-only
    - Primer3 capture is always `incomplete` when enabled because only
      Primer3-returned pairs rejected by GENtle post-filters are inspectable;
      Primer3-internal and single-primer rejections remain aggregate-only
  - construct-reasoning projection:
    - `construct_reasoning_graph_id` links the report to one persisted
      `gentle.construct_reasoning_graph.v1`
    - selected pairs become weighted-rule decision nodes; bounded rejected
      pair intervals become `ContextEvidence` with report/op/run provenance
    - the graph fingerprint binds the exact primer-design report content, so
      report drift makes graph freshness stale
  - region-level biological exclusion sources are not consulted by this
    selector slice. In particular, no repeat, common-variant, or
    paralogue-shared interval evidence is emitted, and absence must not be read
    as checked and clear. Homopolymer diagnostics remain candidate score terms,
    not a merged excluded-region track.
    New reports record `excluded_region_analysis_status: not_run` and an
    explanatory `excluded_region_analysis_reason`; legacy reports omit both.
  - generic construct-graph refresh never rewrites this report-bound graph; a
    stale selection graph is refreshed by rerunning `DesignPrimerPairs`.
  - existing `ConstructReasoningOverlay::from_graph()` projects the evaluated
    near-miss evidence rows into the linear DNA map; no separate
    primer-exclusion track or GUI-local scoring is used.
  - mutating artifact materialization per accepted pair:
    - one forward-primer sequence (`..._fwd`)
    - one reverse-primer sequence (`..._rev`)
    - one predicted amplicon sequence (`..._amplicon`) built from:
      full forward primer + template interior + reverse-primer reverse-complement
      (including non-annealing 5' tails)
    - one per-pair pool container containing all three artifacts
  - optional insertion-anchored context:
    - `insertion_context` (present when `DesignInsertionPrimerPairs` is used)
    - records requested forward/reverse anchor positions, extension sequences,
      window constraints, shift budget, and per-pair compensation rows
      (`forward_anchor_shift_bp`, `reverse_anchor_shift_bp`,
      compensation segments, and compensated primer/tail strings)

`ExportPrimerDesignReport` contract:

- Purpose:
  - write a persisted `gentle.primer_design_report.v1` record to JSON so
    wrappers such as ClawBio can attach the ranked primer-pair report as a
    deterministic artifact.
- Operation payload:
  - `report_id`: existing primer-design report id.
  - `path`: JSON output path.
- This is the preferred workflow-owned export step after `DesignPrimerPairs`
  when a chat adapter needs a shareable report rather than only stdout.

`AssessPrimerPairSpecificity` contract (implemented local BLAST baseline):

- Purpose:
  - confirm one saved or explicit primer pair against a prepared local
    reference-genome BLAST database
  - provide Primer-BLAST-like specificity evidence without remote NCBI
    submission or GUI-only state
- Operation payload:

```json
{
  "AssessPrimerPairSpecificity": {
    "primer_report_id": "tp73_roi_primers_v1",
    "pair_rank": 1,
    "target_genome_id": "GRCh38.p14",
    "policy": {
      "specificity_check": "report_only",
      "max_target_amplicon_bp": 4000,
      "readiness_max_target_amplicon_bp": 1000,
      "readiness_max_target_amplicon_bp_source": "explicit_override",
      "readiness_max_target_amplicon_bp_reason": "The caller explicitly set the ordinary assay/readiness ceiling to 1000 bp.",
      "report_detail_mode": "compact",
      "min_primer_coverage_fraction": 0.8,
      "max_3prime_mismatches": 0,
      "three_prime_window_bp": 5,
      "min_total_mismatches_to_unintended_target": 2,
      "allow_same_gene_splice_variants": false,
      "max_hits_per_primer": 500
    },
    "path": "tp73_pair1_specificity.json"
  }
}
```

- Explicit primers are accepted instead of `primer_report_id`:
  - `forward_primer`: forward oligo/annealing sequence
  - `reverse_primer`: reverse oligo/annealing sequence
- For saved `PrimerDesignReport` pairs, the recorded
  `non_annealing_5prime_tail_bp` is removed before BLAST. Full oligos and tails
  remain in report provenance.
- Target scope in v4:
  - prepared genomic-DNA or transcriptome-cDNA BLAST indexes
  - runs through GENtle's existing BLAST index/preflight machinery
  - local `blastn-short` is used with short-query settings; GENtle validates
    the database sequence count through `blastdbcmd` and sets
    `-max_target_seqs` above that count, so auxiliary contigs are not silently
    excluded
  - `max_hits_per_primer` is a post-search review threshold and never becomes
    the BLAST subject cap
- Report schema:
  - `gentle.primer_specificity_report.v4` with
    `gentle.primer_specificity_policy.v2`
  - both `AssessPrimerPairSpecificity` and
    `ImportPrimerPairSpecificityHandoff` persist the report in project
    metadata. Its stable, content-derived `report_id` binds the primer pair,
    intended-target model, database content, and effective policy
  - follows the shared computational-artifact vocabulary:
    `report_id`, `op_id`, `run_id`, `generated_at_unix_ms`,
    `primary_seq_id`, `related_seq_ids`, `external_inputs`,
    `request_summary`, `effective_settings_summary`, `reopen_hint`, and
    `export_kinds`
  - `external_inputs` binds prepared-BLAST database identity/release and
    content fingerprint, the genome catalog, and any cited primer-design
    report or externally executed specificity handoff
  - `design_provenance` cites the persisted primer-design report and selected
    pair when resolvable. Explicit raw/commercial/literature primer strings do
    not acquire a reconstructed design rationale; their design provenance is
    `not_run`
  - a cited design pair and the normalized full oligos actually assessed are
    fingerprinted independently with
    `sha256_canonical_forward_reverse_full_sequence_json_v1`. Matching digests
    yield `pass`; content drift yields `fail`; an unresolved citation remains
    `incomplete`
  - `characterization_dimensions[]` gives independent
    `pass|fail|incomplete|not_run` states for design provenance, oligo-pair QC,
    genomic specificity, transcriptome specificity, intended isoform coverage,
    junction/genomic-carryover interpretation, search completeness,
    known-variant screening, and repeat/low-complexity screening. A cited
    design report is not promoted into a freshly rerun oligo-QC pass, and no
    variant/repeat clearance is claimed when those analyses were not run
  - includes BLAST binary preflight, per-primer BLAST invocation provenance,
    BLAST database content identity, input primers, policy, aggregated
    warnings, primer hits, candidate amplicons, and a summary badge
  - primer hits report identity, coverage, aligned-HSP mismatches, unaligned
    query bases, their sum as `effective_mismatches`, exact 3' terminal mismatch
    count where prepared subject sequence can be fetched, strand, and 1-based
    subject coordinates; wrapped BLAST identifiers such as
    `gb|KI270750.1|` are normalized against the prepared FASTA index while the
    raw identifier remains available
  - BLAST remains the exhaustive locus finder, while
    `policy.full_alignment.mode = disabled|best_effort|required` controls a
    second full-query semiglobal alignment against a fetched subject window.
    Newly constructed requests default to `best_effort`; legacy serialized
    policies lacking this additive field retain the former `disabled` behavior.
    The stored alignment records the
    complete query/subject strings, relation line, CIGAR, indels, internal and
    3-prime mismatches, binding coordinates, scoring parameters, and subject
    window hash. `required` makes a missing full alignment an incomplete
    search; no local HSP may stand in for unavailable end-to-end evidence
  - new primer-specificity handoffs request BLAST `qlen`, `qseq`, and `sseq`.
    A validated HSP spanning the complete primer is converted directly into
    the same persisted full-alignment evidence, avoiding one reference lookup
    and dynamic-programming alignment per HSP. Partial or legacy 13-column
    HSPs retain the subject-window semiglobal path. The stored `algorithm`
    distinguishes these methods; only GENtle's normal completeness, intended
    product, allowlist, and finalization gates can produce a native pass
  - prepared transcriptome resources add optional `subject_annotation` to
    primer hits and candidate amplicons. It includes the full and normalized
    transcript stable id, full and normalized gene stable id, `gene_symbol`,
    biotypes, description, and explicit annotation source when those values
    were present in the prepared reference. GenBank subjects instead retain
    their versioned accession, `/gene`, and `/db_xref` values. Missing values
    remain absent and use `annotation_source = unavailable`
  - BLAST tabular output remains compact and does not request `stitle`;
    descriptive identity is joined from the prepared sidecar instead of
    repeating a potentially long title on every HSP
  - `query_coverage_fraction` is computed independently for every HSP from its
    inclusive `qstart..qend` span divided by primer length; BLAST `qcovs`
    cannot promote a short HSP because that field is aggregated per subject.
    Unaligned query bases count as effective mismatches when candidate products
    are ranked and screened
    and therefore a short otherwise mismatch-free HSP cannot look like a full
    clean primer match
  - candidate amplicons include forward/reverse products plus Primer-BLAST-style
    forward/forward and reverse/reverse warning products. Genomic left/right
    ordering is derived from subject strand and coordinates rather than primer
    role, so minus-strand targets are paired without swapping assay roles
  - accepted hits are indexed by subject, strand/orientation, and coordinate
    before product pairing. `compaction.pairing_candidate_comparison_count`
    records the bounded work performed instead of relying on wall-clock
    measurements
  - `intended_target.model = transcript_set` supports assays intended to
    amplify multiple transcripts. `expected_products[]` records each intended
    subject and optional exact product range in its `genomic_dna` or
    `transcriptome_cdna` target space; a transcriptome pass requires exactly
    one compatible product for every declared transcript
  - intended transcript matching ignores only an optional final numeric
    version suffix on an Ensembl transcript stable id. For example,
    `ENST00000465287` matches observed `ENST00000465287.1`, while the full
    observed id remains in hit/product provenance. Arbitrary dotted ids,
    Ensembl gene ids, and non-Ensembl accessions remain exact-match
  - gene ids and symbols are descriptive annotation only. Transcript subject
    ids remain authoritative for intended-product matching; another transcript
    of the same annotated gene is still unintended unless it was declared
  - intended genomic products are resolved only from an explicit prepared-FASTA
    subject and computed genomic binding interval; cDNA amplicon length is
    retained as provenance but is never used to identify a genomic product
  - a junction-spanning RT-PCR primer may legitimately have no contiguous
    intended genomic product; in a genomic-DNA database that case is evaluated
    as carryover/off-target evidence rather than forced into an invented target
  - unintended compatible products fail when their combined effective
    mismatches remain below `min_total_mismatches_to_unintended_target`
  - `reviewed_off_target_allowlist[]` can retain an exact, human-reviewed
    product without relabelling it as intended. Each allowance binds target
    space, subject, optional exact coordinates, reviewer, reason, and an
    evidence SHA-256. It never converts an incomplete search into a pass and
    is deliberately narrower than same-gene or same-family blanket waivers
  - `intended_target`, `genomic_specificity`, and
    `transcriptome_specificity` keep genomic-DNA carryover/off-target evidence
    separate from whole-transcriptome/cDNA cross-amplification. Only the
    assessment matching the inspected BLAST index kind is populated by one run;
    the other remains `not_run`
  - `intended_target.expected_products[]` records explicit expected subjects,
    optional one-based ranges or product lengths, transcript ids, and evidence
    source per BLAST target space. `genomic_target_geometry_known` distinguishes
    resolved annotation geometry from a missing geometry; GENtle does not infer
    a genomic product from cDNA length. If a transcript product has neither
    portable coordinates nor a known length, exactly one compatible product on
    that named transcript is required
  - each target-space assessment records expected and observed intended-product
    counts. Intended isoform coverage is derived from those counts independently
    of the transcriptome-specificity verdict, so an off-target product can fail
    specificity without erasing evidence that all intended isoforms were covered
  - `search_completeness` is an enforceable result state containing the
    validated database sequence count, required `-max_target_seqs`, minimum
    observed command limit, command count, and explanatory reason. If
    completeness is not proven, `summary.status = incomplete` and
    `specificity_pass = false` regardless of observed hits
  - `policy.readiness_max_target_amplicon_bp` is the ordinary assay gate.
    `policy.max_target_amplicon_bp` remains the broader exploratory ceiling.
    Compatible products between them carry `long_product_warning = true` and
    are summarized separately without becoming ordinary in-window specificity
    failures. `readiness_max_target_amplicon_bp_source` records
    `legacy_policy_max`, `explicit_override`, or
    `transcript_assay_panel_allowed_range`, while
    `readiness_max_target_amplicon_bp_reason` records the human-readable
    derivation (for example, the panel's declared 70-1,000 bp range)
  - newly created reports default to compact detail:
    `compaction` records raw and retained hit/product counts and compact mode
    omits rejected HSPs plus nonviable products from inline JSON. Full mode is
    an explicit lossless/debug option. Imported handoffs record raw TSV path,
    size, and SHA-256 under `raw_detail_artifacts`
  - compatibility: a v3 report without `compaction`,
    `within_readiness_amplicon_range`, or `raw_detail_artifacts` deserializes
    as full detail, treats its historical products as inside the one legacy
    ceiling, and leaves raw-artifact citations empty. A v1 policy missing
    v2 fields likewise keeps the old single-ceiling/full-detail meaning
  - `primers list-reports` lists specificity summaries alongside
    primer-design summaries. `primers show-report REPORT_ID` and
    `primers export-report REPORT_ID OUTPUT.json` accept either persisted
    report kind, so headless adapters can reopen or export the artifact
    without rerunning BLAST
  - project fact projection exposes persisted specificity artifacts as
    `report.exists = "primer_specificity"`. Reports linked to a project
    sequence are also projected as GUI lineage analysis nodes; PCR Designer
    reopens the specificity summary and its design citation
  - `primers specificity-alignment-html REPORT_ID OUTPUT.html` is a pure
    projection of the stored v4 JSON. It renders every retained full-primer
    alignment and its provenance without rerunning BLAST or changing state
- Additive policy vocabulary reserved for design-time parity:
  - `specificity_check = none|report_only|require_pass`
  - `specificity_target_genome_id`
  - `allow_same_gene_splice_variants`
  - mask booleans: `avoid_known_variants`, `avoid_rmsk_repeats`,
    `avoid_low_complexity`
  - transcript-aware policy enums:
    `exon_junction_policy = no_preference|must_span|must_not_span` and
    `intron_separation_policy = no_preference|must_separate_by_intron`
  - the current implementation stores/report-routes these knobs, while hard
    design-time rejection by genome-wide specificity and binding-site masks
    remains a follow-up.

External BLAST handoff for wrapper-owned execution:

- The handoff lifecycle is also available through the shared engine contract:
  - `PreparePrimerPairSpecificityHandoff` returns
    `OpResult.primer_specificity_handoff` and writes the deterministic command
    bundle without launching `blastn`;
  - `ImportPrimerPairSpecificityHandoff` returns
    `OpResult.primer_specificity_report` after reading the completed declared
    TSVs and can optionally write that report through `path`.
  - These operation payloads can be submitted unchanged through CLI `op` or
    `workflow`, MCP `op`, JavaScript `apply_operation`, and Lua
    `apply_operation`; the convenience shell commands below use the same
    operations.
- `primers specificity-plan` resolves the same saved or explicit primer pair,
  prepared genome, BLAST database, and effective policy without launching
  `blastn`.
- It writes a deterministic `gentle.primer_specificity_handoff.v1` bundle with:
  - one query FASTA and one structured command for each primer;
  - authoritative `program` plus `args[]` fields (the rendered `command_line`
    is convenience text only);
  - explicit `-out` paths, accepted exit code `0`, database prefix plus current
    content fingerprint/index kind, resolved genome/catalog/cache provenance,
    a subject cap above the validated database sequence count,
    `search_completeness`, intended-target geometry, policy, an
    `all_commands_success` completion policy, and an import command.
- An outer scheduler owns process lifetime and completion. It should run both
  commands, require a successful exit code, and only then call
  `primers specificity-import HANDOFF.json`.
- Replanning the same deterministic handoff removes its previously declared
  output TSVs before returning, so a fresh plan cannot silently reuse stale
  BLAST results.
- Import requires the declared output files but does not infer process state
  from their size: an empty BLAST TSV is a valid no-hit result. It parses the
  completed outputs and applies the same hit, amplicon, and pass/fail logic as
  the inline `AssessPrimerPairSpecificity` path.
- Import and whole-panel finalization re-inspect the database and refuse a
  handoff when the content fingerprint or index kind changed at the same
  prefix. Legacy handoffs without a proven `-max_target_seqs` remain readable,
  but their reports are `incomplete` and cannot become accepted evidence.
- GENtle never executes a `command_line` read from a handoff. Adapters should
  dispatch the stored `program` and `args[]` directly.

Whole-panel external specificity acceptance:

- `primers transcript-assay-specificity-plan PANEL_REPORT_ID --target-genome
  GENOME_ID --output-dir DIR` emits
  `gentle.transcript_assay_panel_specificity_handoff.v1`. It binds the current
  transcript-assay panel digest, selected assay ids/ranks and annealing
  sequences, `gentle.primer_specificity_policy.v2`, prepared-genome identity,
  BLAST database prefix/content fingerprint/index kind/options, nested handoff
  schemas, and structured commands.
- The adjacent
  `gentle.transcript_assay_panel_specificity_execution_manifest.v1` template is
  process evidence, not a biological decision. For every declared command the
  scheduler returns `command_id`, `assay_id`, `exit_code`, `output_path`,
  `output_size_bytes`, and `output_sha256`. A completed empty output has size
  zero and the SHA-256 of empty bytes; it is distinct from a missing output or
  an absent/non-success exit code.
- The scheduler invokes `primers transcript-assay-specificity-finalize
  HANDOFF.json EXECUTION_MANIFEST_JSON_OR_@FILE` even when one process fails.
  GENtle validates command coverage, uniqueness, exit status, byte identity,
  current panel identity, primer/policy/database provenance, and then applies
  the same specificity interpretation as the inline route.
- Finalization returns
  `gentle.transcript_assay_panel_specificity_acceptance.v2` with exactly one of
  `pass`, `specificity_fail`, `not_assessed`, or `incomplete`.
  `specificity_fail` means all process evidence was complete but at least one
  assay failed GENtle's biological policy. `not_assessed` means all commands
  and output identities were complete but intended-target geometry was
  unavailable for at least one assay. `incomplete` covers
  failed/missing/duplicate execution, stale panel or primer state, altered
  handoffs, and provenance mismatch. Separate assay-id arrays identify
  biological failures, not-assessed results, incomplete interpretations, and
  failed external processes.
- Only `pass` sets `accepted = true`. A complete panel result (`pass`,
  `specificity_fail`, or `not_assessed`) atomically retains all per-assay
  assessments on the persisted panel, keyed by assay and BLAST target kind.
  An `incomplete` result remains available as the returned or exported
  acceptance evidence but attaches neither a partial acceptance nor partial
  assay assessments to the panel. This keeps failed/incomplete execution
  distinct from `not_run` without making a half-validated panel look current.
  An optional `--path` writes the same acceptance object returned by the shell
  command.
- `primers transcript-assay-specificity-redesign REQUEST_JSON_OR_@FILE`
  consumes a content-bound `specificity_fail` acceptance plus the exact
  original `DesignTranscriptAssayPanel` operation. It retains passing assays,
  excludes rejected primer footprints, and either proposes replacement pairs
  preserving the same biological target binding or emits a narrowly scoped
  infeasibility statement. Every substitute records the rejected pair,
  failure report/acceptance, causes, changed decision, backend, attempt, and
  `replaces_assay_id`; it remains a candidate until it passes the same cDNA
  and genomic gates. The route works on a detached clone and does not mutate
  the approved source state.
- NCBI e-PCR is not part of this contract. A future provider-neutral
  Primer-BLAST evidence importer may supplement, but must not weaken, the
  reproducible prepared-genome local-BLAST gate.

Simple PCR constraint handoff:

- ClawBio should treat a simple PCR request as four explicit constraints before
  calling `DesignPrimerPairs`:
  - template sequence id or retrieval workflow,
  - core ROI that must be inside the amplicon,
  - allowed forward/reverse flank windows,
  - min/max amplicon and primer/Tm/GC limits.
- For the smallest safe default, first extract a compact context around the
  core ROI, then express the flank windows in that extracted template. This
  keeps internal search bounded and avoids accidental whole-locus primer scans.
- Canonical offline example:
  `docs/examples/workflows/simple_pcr_primer_design_offline.json`.
- The example also exports `artifacts/simple_pcr_demo_primers.protocol.svg`
  with `RenderProtocolCartoonSvg { protocol: "pcr.assay.pair", ... }`, so chat
  adapters can show the user how the core ROI, primer windows, primers, and
  final amplicon relate before exposing the JSON primer report.

`DesignInsertionPrimerPairs` contract (implemented MVP):

- Purpose:
  - insertion-first wrapper around deterministic pair-primer design when the
    user already knows insert extensions and requested insertion anchors.
- Operation payload shape:

```json
{
  "DesignInsertionPrimerPairs": {
    "template": "seq_id",
    "insertion": {
      "requested_forward_3prime_end_0based_exclusive": 620,
      "requested_reverse_3prime_start_0based": 700,
      "forward_extension_5prime": "GAATTC",
      "reverse_extension_5prime": "CTCGAG",
      "forward_window_start_0based": 560,
      "forward_window_end_0based_exclusive": 650,
      "reverse_window_start_0based": 660,
      "reverse_window_end_0based_exclusive": 760,
      "max_anchor_shift_bp": 12
    },
    "forward": {
      "min_length": 20,
      "max_length": 30
    },
    "reverse": {
      "min_length": 20,
      "max_length": 30
    },
    "pair_constraints": {
      "require_roi_flanking": false
    },
    "min_amplicon_bp": 120,
    "max_amplicon_bp": 1200,
    "max_tm_delta_c": 2.0,
    "max_pairs": 200,
    "report_id": "tp73_insert_v1"
  }
}
```

- MVP behavior:
  - the insertion block is normalized first (IUPAC extension validation +
    anchor/window bounds checks)
  - forward/reverse primer windows are enforced from insertion windows
  - forward/reverse non-annealing tails are set from insertion extensions
  - primer design backend selection remains identical to `DesignPrimerPairs`
    (`auto|internal|primer3`)
  - resulting report is the same primer-report schema with populated
    `insertion_context` rows for shift/compensation inspection
  - no dedicated GUI form yet; operation is available through `op`/workflow
    payloads.

`PrepareRestrictionCloningPcrHandoff` contract (implemented v1):

- Purpose:
  - take one persisted `DesignPrimerPairs` result pair and turn it into a
    cloning-aware handoff with restriction-site tails matched against a chosen
    destination vector
  - keep core primer proposal/ranking unchanged while creating new extended
    primer artifacts and one restriction-ready amplicon artifact for downstream
    digest/ligation staging
- Operation payload shape:

```json
{
  "PrepareRestrictionCloningPcrHandoff": {
    "template": "seq_id",
    "primer_report_id": "tp73_pairs_v1",
    "pair_index": 0,
    "destination_vector_seq_id": "pgl3_mcs",
    "mode": "directed_pair",
    "forward_enzyme": "EcoRI",
    "reverse_enzyme": "HindIII",
    "forward_leader_5prime": "GC",
    "reverse_leader_5prime": "AT"
  }
}
```

- Baseline behavior:
  - validates that `primer_report_id` belongs to `template` and that
    `pair_index` exists in the persisted report
  - derives:
    - extended forward primer =
      `forward_leader_5prime + forward restriction site + original forward primer`
    - extended reverse primer =
      `reverse_leader_5prime + reverse restriction site + original reverse primer`
    - predicted tailed amplicon from the full extended primer sequences
  - preserves annealing Tm/GC/hit semantics from the original annealing segment
    while recomputing full-oligo secondary-structure and pair-dimer heuristics
    as advisory diagnostics
  - blocking compatibility checks:
    - vector site absent or non-unique
    - tailed amplicon site counts imply internal collisions instead of only
      terminally added restriction sites
    - directed-pair order disagrees with vector MCS order
      (`mcs_expected_sites`) or, if absent, unique-cut order by cut position
  - successful runs materialize graph-visible artifacts:
    - one extended forward primer sequence
    - one extended reverse primer sequence
    - one predicted tailed amplicon sequence
    - one per-handoff pool container
  - successful runs also persist structured downstream hints:
    - one staged `Workflow` containing the runnable `PcrAdvanced` + insert
      digest + vector digest steps
    - one `PcrAdvanced` operation payload using the full tailed oligos with
      preserved `anneal_len`
    - one insert `Digest` payload
    - one vector `Digest` payload
    - one ligation JSON snippet placeholder
- Report schema:
  - `gentle.restriction_cloning_pcr_handoff.v1`
  - key fields include:
    - `template`, `primer_report_id`, `pair_index`, `pair_rank`
    - `destination_vector_seq_id`
    - `mode`, selected enzymes, optional leaders
    - original and extended primer records
    - created artifact ids
    - tailed amplicon length plus 5'/3' sequence previews
    - extended pair dimer diagnostics
    - `compatibility` summary with vector-site counts, insert-site counts,
      cut positions, blocking errors, and warnings
    - `workflow_hints` with suggested downstream operation payloads, including
      one staged pre-ligation workflow and the remaining ligation placeholder

`PcrOverlapExtensionMutagenesis` contract (implemented baseline):

- Purpose:
  - deterministic overlap-extension insertion/deletion/replacement mutagenesis
    planning + staged product materialization in the main operation graph.
- Operation payload shape:

```json
{
  "PcrOverlapExtensionMutagenesis": {
    "template": "seq_id",
    "edit_start_0based": 620,
    "edit_end_0based_exclusive": 640,
    "insert_sequence": "GGTACC",
    "constraints": {
      "overlap_bp": 24,
      "outer_forward": {
        "min_length": 20,
        "max_length": 30
      },
      "outer_reverse": {
        "min_length": 20,
        "max_length": 30
      },
      "inner_forward": {
        "min_length": 18,
        "max_length": 28
      },
      "inner_reverse": {
        "min_length": 18,
        "max_length": 28
      }
    },
    "output_prefix": "tp73_oe_mut"
  }
}
```

- Baseline behavior:
  - `edit_start_0based..edit_end_0based_exclusive` defines the replaced region
    on the original template.
    - insertion: `edit_start == edit_end` and `insert_sequence` non-empty
    - deletion: `insert_sequence` empty and `edit_end > edit_start`
    - replacement: both deletion and insertion are non-empty
  - inner primers are chosen upstream/downstream of the edit and receive dynamic
    5' overlap tails derived from the mutant sequence so stage-1 products share
    one explicit overlap segment (minimum `overlap_bp`).
  - outer primers amplify both stage-1 fragments and the stage-2 final mutant
    amplicon.
  - operation materializes graph-visible artifacts:
    - primers: `..._outer_fwd`, `..._outer_rev`, `..._inner_fwd`, `..._inner_rev`
    - stage-1 products: `..._stage1_left`, `..._stage1_right`
    - final stage-2 mutant: `..._mutant`
    - three per-stage pool containers (left, right, final)
  - operation warnings include deterministic candidate-search limit notices when
    the combinatorial search budget is exhausted.
  - insertion/replacement runs now also emit
    `OpResult.protocol_cartoon_preview` for built-in protocol
    `pcr.oe.substitution`, including deterministic
    `flank_bp`/`overlap_bp`/`insert_bp` geometry and bound template overrides
    (`gentle.protocol_cartoon_template_bindings.v1`) for adapter rendering.

`DesignQpcrAssays` contract (implemented baseline):

- Purpose:
  - propose ranked qPCR assays with three oligos (forward primer, reverse
    primer, internal probe) for one linear template.
- Operation payload shape:
  - same core fields as `DesignPrimerPairs` plus:
    - `probe` (`PrimerDesignSideConstraint`)
    - `max_probe_tm_delta_c` (probe Tm distance to mean primer Tm)
    - `max_assays` (result cap)
    - optional `transcript_targeting`
      - `source_feature_id`
      - `mode = shared_gene | distinguish_transcript`
      - `transcript_id` required for `distinguish_transcript`
      - optional
        `specificity_evidence = junction_only | unique_exon_or_chain | either_prefer_junction`
  - `pair_constraints` is supported identically to `DesignPrimerPairs` and
    applies to forward/reverse pair proposal before probe selection.
- Current baseline behavior:
  - forward/reverse pair generation follows the same backend routing as
    `DesignPrimerPairs` (`auto|internal|primer3` for pair proposal).
  - probe selection is deterministic, constrained to amplicon interior, and
    reuses the same side sequence-constraint fields (`fixed_5prime`,
    `fixed_3prime`, motifs, locked positions).
  - probe Tm gating is enforced via `max_probe_tm_delta_c`; ranking now
    prefers probe-based qPCR layouts where the probe Tm is about 7.5 C above
    the forward/reverse primer mean when such candidates are available.
  - when `transcript_targeting` is present, qPCR design first derives
    transcript-local exon/junction ROIs from the selected splicing group on
    `TargetGroupTargetStrand` and only then reuses the normal qPCR backend
    logic on those transcript templates.
  - `shared_gene` mode prefers assays whose amplicon context is shared across
    all transcripts in the selected gene/group and records an explicit fallback
    summary if only a largest supported subset is assayable.
  - `distinguish_transcript` mode supports three evidence policies:
    - `junction_only`
      - requires at least one primer spanning a transcript-unique exon-exon
        junction
    - `unique_exon_or_chain`
      - allows transcript discrimination through a primer placed on a
        transcript-unique exon or exon-chain context
    - `either_prefer_junction`
      - prefers a junction-spanning distinguishing primer when available
      - otherwise falls back to a transcript-unique exon/exon-chain assay
  - transcript-distinguishing modes fail explicitly when the requested evidence
    type cannot be satisfied.
- Report schema:
  - `gentle.qpcr_design_report.v1`
  - includes ranked `assays[]` with forward/reverse/probe oligos, amplicon
    window, and rule flags.
  - each retained assay includes exact additive `score_terms[]` under
    `score_model = gentle_qpcr_assay_rank_v1` and
    `score_direction = higher_is_better`:
    - all inherited primer-pair ranking terms
    - `probe_tm_offset_from_preferred`, the absolute distance from the
      preferred probe-minus-primer-mean Tm offset
    - `probe_amplicon_midpoint_distance`, the absolute probe/amplicon midpoint
      distance in bases
    - zero-weight observational terms for probe self-complementarity,
      probe/primer complementarity, and 3'-anchored probe/primer
      complementarity; v1 records these diagnostics without changing the
      established score
  - `score_decomposition_status` is `pass` when a retained assay has exact
    terms and `not_run` when no assay was retained.
  - bounded selection provenance includes:
    - `rejected_near_misses[]` for evaluated pair/probe combinations rejected
      by probe placement or Tm checks
    - `near_miss_capture` with status, scope, requested/effective limit, and
      eligible/retained/omitted deterministic work counts
    - a parallel `QpcrDesignRejectionReason` vocabulary with
      `count_for_reason` reconciliation against the nested primer/probe census
    - `incomplete` status for Primer3-hidden pair rejection space, internal
      pair-evaluation truncation, and transcript-local rejections that cannot
      be projected to source coordinates
  - `construct_reasoning_graph_id` links the report to a
    report-content-fingerprinted graph. Retained assays become weighted-rule
    decisions and bounded coordinate-bearing rejected assays become
    non-verdict `ContextEvidence`; report drift makes the graph stale.
  - when transcript-aware targeting is active, persisted reports also include:
    - report-level `transcript_targeting`
    - report-level `transcript_targeting_result`
      - chosen mode/transcript
      - requested transcript label when available
      - transcript count considered
      - selected support count/fraction
      - realized specificity evidence when transcript-distinguishing mode was
        requested
      - whether shared-mode fallback was used
      - deterministic warnings
    - per-assay `transcript_context`
      - design transcript id/label
      - support count/fraction and supported transcript ids
      - mapped source exon ranges for forward/reverse/probe/amplicon
      - covered junction labels
      - whether each oligo spans a junction
      - genomic-DNA carryover risk based on junction oligos/probes, the mapped
        genomic equivalent of the cDNA amplicon, and the configured max amplicon
        window
      - realized specificity evidence
      - whether the assay satisfies the requested targeting intent
  - includes `best_assay_probe_placement` and `best_assay_summary` so
    shell/CLI/GUI reopen flows can inspect one compact persisted explanation of
    the current top retained assay without re-deriving it locally.
  - includes qPCR rejection summary with pair-level and probe-level counters.
  - region-level repeat, variant, and paralogue exclusion evidence remains
    absent because those sources are not consulted during qPCR selection;
    absence is `not_run`, never an implied clear result. New reports record
    `excluded_region_analysis_status: not_run` plus an explanatory
    `excluded_region_analysis_reason`; legacy reports omit both.

`TestCdnaPcr` / `TestCdnaQpcr` / `TestCdnaQpcrFasta` contract (implemented baseline):

- Purpose:
  - test already-chosen PCR/qPCR oligos against transcript-derived cDNA
    templates for one annotated splicing group.
  - use the same transcript-native cDNA derivation path as transcript-aware
    qPCR design, so assay testing is independent of external protein or
    transcript evidence services.
  - optionally screen supplied cDNA/ncRNA FASTA or FASTA.gz transcript catalogs
    with the same primer/probe/product rules, for example complete Ensembl
    `cdna.all.fa.gz` plus `ncrna.fa.gz` files.
- Operation payload shape:
  - `seq_id`
  - `source_feature_id`
    - existing zero-based GENtle feature id used to resolve the transcript/gene
      splicing group
    - used only by `TestCdnaPcr` and `TestCdnaQpcr`
  - `cdna_fasta_paths[]`
    - one or more transcript FASTA/FASTA.gz files for `TestCdnaQpcrFasta`
    - records are streamed directly from disk and may include large Ensembl
      cDNA/ncRNA catalogs
  - `forward_primer`
  - `reverse_primer`
  - `probe` for `TestCdnaQpcr` and `TestCdnaQpcrFasta`
  - optional `transcript_id`
    - limits testing to one derived transcript template when supplied
  - optional `transcript_order`
    - transcript-derived tests only; controls report/map row order
  - optional `transcript_map_coordinate_mode`
    - transcript-derived tests only; `cdna` keeps the existing transcript-local
      cDNA axes, while `genomic_aligned` draws the map on the shared
      source/genomic axis
  - optional `min_amplicon_bp` / `max_amplicon_bp`
  - optional `max_mismatches`
    - defaults to `0`
  - optional `require_3prime_exact_bases`
    - enforced for forward/reverse primer hits
  - optional `path`
    - writes the same report JSON returned in the operation result
  - optional `svg_path`
    - writes the report-owned transcript-map SVG to a separate file for
      ClawBio/OpenClaw artifact collection or direct CLI review
  - optional `materialize_products`
    - default `false`; when true, detected transcript-derived cDNA products are
      created as linear GENtle sequence entries and grouped into one
      singleton/pool product container
  - optional `product_output_prefix`
    - deterministic prefix for materialized product sequence ids
  - optional `product_gel_svg_path`
    - writes a pool-gel SVG for the materialized product container; supplying
      this path implies product materialization even when
      `materialize_products` is omitted
  - optional `product_gel_ladders`
    - ladder names passed to the shared DNA pool-gel renderer
- Current baseline behavior:
  - templates are derived from the selected splicing group on
    `TargetGroupTargetStrand`.
  - forward primer hits are searched on the transcript-derived cDNA strand.
  - reverse primer hits are searched as reverse-complement binding sites.
  - qPCR probes are accepted on either cDNA orientation, but must fall inside
    the primer-bounded amplicon interior.
  - matching is IUPAC-aware and exact by default, with optional mismatch
    tolerance and deterministic 3' exactness gating for primers.
  - UTR/noncoding transcript bases are not guessed from genomic context; only
    derived cDNA templates from the selected transcript features are tested.
  - FASTA screens preserve Ensembl-style record ids from the first header token,
    accept `.gz` input, and can take multiple files in one run.
  - FASTA screens count every scanned record in `transcript_count` but, unless
    `transcript_id` is supplied, report only detected transcript rows to avoid
    huge reports for complete transcript catalogs.
  - transcript-derived PCR/qPCR tests remain non-mutating by default. When
    product materialization is requested, each detected cDNA product becomes a
    derived linear DNA sequence, all products from the assay are grouped into
    one product container, and multiple nonspecific products are represented as
    one pool lane in the existing `RenderPoolGelSvg` renderer.
  - product materialization is repeat-safe. Product rows carry deterministic
    assay signatures, repeat runs reuse matching sequence ids and matching
    product containers, and structured output separates
    `created_product_seq_ids[]` from `reused_product_seq_ids[]`.
  - product-gel output includes text-first companion data:
    `gel_band_rows[]` for machine display and `gel_summary_lines[]` for
    compact Telegram/terminal narration of lane/band sizes and merged
    non-specific products.
  - if no products are detected, materialization does not fail the assay; it
    creates no product sequence/container and reports why no product gel was
    written.
  - every report includes `transcript_map` with schema
    `gentle.cdna_assay_transcript_map.v1`; the embedded SVG draws each shown
    transcript on its own cDNA coordinate axis, with amplicons, symbolic
    forward/reverse primer hits, probe hits, group/source exon labels (`E1`,
    `E2`, ...) ordered by transcript strand, source-exon identity
    colors/patterns, and exon-junction ticks overlaid where products are
    functional; requested forward/reverse primer and probe sequences are shown
    once in the legend, forward primer glyphs are one-sided above the cDNA
    axis, reverse primer glyphs are one-sided below it, probe glyphs are
    one-sided on the detected probe-binding strand, crowded transcript sets
    auto-wrap into columns, and SVG tooltips retain the transcript-local exon
    ordinal when it differs from the group label.
  - `transcript_map.column_count` and `transcript_map.rows_per_column` describe
    the selected SVG layout so adapters can size previews without parsing the
    SVG.
  - `transcript_order` is included at report level and in `transcript_map`.
    It defaults to `transcript_id`; transcript-derived tests may request
    `genomic_first_exon`, `genomic_last_exon`, or `antisense_first_exon`.
  - `transcript_map_coordinate_mode` is included at report level, and
    `transcript_map.coordinate_mode` repeats the realized map axis. The default
    is `cdna`; `genomic_aligned` SVGs use shared source/genomic coordinates for
    exon blocks, primer/probe glyphs, and amplicon source spans.
- Report schema:
  - `gentle.cdna_assay_test_report.v1`
  - report-level fields include assay kind, source sequence/feature, group
    label, strand, requested transcript id, primers/probe, construct-length
    summary, genomic-DNA carryover risk summary, oligo QC,
    mismatch/size settings, transcript/product counts, overall status, and
    warnings.
  - `construct_lengths` records forward/reverse/probe oligo lengths plus the
    distinct detected amplicon construct lengths, including min/max when any
    product is found.
  - `oligo_qc` uses schema `gentle.oligo_qc_report.v1` and performs a
    deterministic exact reverse-complement run screen over supplied oligos:
    per-oligo self-complementary and 3'-self runs, pairwise oligo
    complementarity, and primer-side 3' extension-risk warnings. This is a
    portable first-pass QC layer, not a thermodynamic replacement for Primer3 /
    wet-lab assay validation.
  - `oligo_qc.method_reference` records that the vocabulary follows Primer3's
    public `SELF_ANY` / `SELF_END` and `COMPL_ANY` / `COMPL_END` distinction,
    while the implementation remains independent GENtle Rust code with no
    vendored or translated Primer3 source.
  - report-level fields include assay kind, template source kind, source
    paths, source sequence/feature, group label, strand, requested transcript
    id, primers/probe, mismatch/size settings, transcript/product counts,
    overall status, and warnings.
  - `pair_id` is the canonical ordered forward/reverse physical-oligo identity
    shared with transcript-panel handoffs. `assay_test_id` additionally binds
    that pair to the template source, optional probe, transcript filter, assay
    limits, and display settings, so repeated tests of one pair remain distinct
    and auditable. Older reports omit both fields and remain readable.
  - `transcript_map` includes `artifact_id`,
    `media_type = image/svg+xml`, canvas dimensions, row/omission counts,
    product count, a compact summary, and the SVG text. FASTA reports use the
    same detected-row rule for the figure unless `transcript_id` requests one
    specific record.
  - per-transcript rows include transcript feature id, transcript id/label,
    optional source path, strand, cDNA length, status, transcript-local exon
    segments when source mapping is available, primer/probe hits, and products.
    Transcript-derived rows prefer transcript/product descriptions over bare
    gene symbols for `transcript_label`, so graphical maps can distinguish
    isoform biology instead of repeating the group label.
  - primer/probe hit rows include local cDNA coordinates, binding sequence,
    binding orientation, mismatch count, mapped source ranges, covered junction
    labels, and whether the hit spans a junction.
  - product rows include local cDNA amplicon coordinates/length, hit indices,
    optional probe hit indices, mapped source ranges, covered junction labels,
    whether the product spans a junction, the genomic-equivalent span/length
    when source mapping is available, and a per-product genomic-DNA carryover
    risk rationale. `product_sequence_sha256` identifies the exact mature-cDNA
    template span in transcript 5-prime-to-3-prime orientation, including both
    binding regions and excluding non-annealing tails or primer-induced
    substitutions. A shared digest proves exact predicted product-sequence
    identity; different digests do not imply that products resolve on a gel.
- Optional materialization schema:
  - `gentle.cdna_assay_product_materialization.v1`
  - fields include assay kind, source sequence/feature, group label, detected
    product count, all product sequence ids, `created_product_seq_ids[]`,
    `reused_product_seq_ids[]`, per-product transcript/amplicon/probe/carryover
    rows, product container id/kind plus `container_created`, optional product
    gel SVG path, `gel_band_rows[]`, `gel_summary_lines[]`, output prefix,
    `idempotent_reuse`, and warnings.

PrimerBank lookup and cDNA continuation (implemented):

- `SearchPrimerBank { request, source_html_path?, path? }` is the shared,
  read-only engine operation used by CLI and Shell adapters. It returns the
  typed report in `OpResult.primerbank_search_report`; `source_html_path`
  selects reproducible offline parsing instead of a live request and `path`
  optionally exports the same full report.
- `gentle.primerbank_search.v1` is a non-mutating typed projection of one
  PrimerBank HTML search response. `query` records the submitted query,
  query-kind, and species selector; `source_url`/`source_kind` distinguish a
  live lookup from parsing saved HTML; `usage_policy_url`, citations, and
  warnings retain external-source context.
- `species_check` records the requested species, observed labels, matched,
  mismatched, and unresolved row counts, plus one typed status: `matched`,
  `mismatch`, `unresolved`, or `not_requested`. Every `genes[]` row carries
  its own `species_match_status`. GENtle evaluates the returned record because
  PrimerBank's exact-ID lookup does not enforce its submitted species filter.
- `genes[]` retains NCBI Gene, GenBank, protein, species, coding-DNA length,
  and description fields where the response provides them. Every
  `primer_pairs[]` row retains the PrimerBank id, amplicon length, detail URL,
  and separate forward/reverse records with sequence, length, Tm, and raw plus
  normalized coordinates.
- PrimerBank primer locations are recorded as 1-based inclusive positions on
  PrimerBank's coding-sequence record. They are not silently projected onto a
  GENtle genomic sequence or contemporary transcript annotation. Reverse
  rows retain their descending raw start/end and also expose an ascending
  interval for display.
- `validation_status = not_assessed_by_gentle` is intentional. Catalog
  presence, a compatible cDNA product, and experimental validation remain
  separate statements. GENtle performs individual lookups only; it does not
  bundle or mirror the PrimerBank database.
- `gentle.primerbank_cdna_test.v1` joins one exact catalog pair and its source
  provenance to the unchanged output from GENtle's existing transcript-aware
  `TestCdnaPcr` route. It tests compatibility with current project transcript
  models but does not establish whole-genome specificity; use the existing
  specificity handoff/import lifecycle for that independent check. `--path`
  writes this complete wrapper; `--svg` writes the nested cDNA transcript map.
  It also records `expected_species`, `primerbank_species`, and
  `species_match_status`. When the selected sequence or transcript carries an
  `organism` annotation, `target_sequence_species` and
  `target_sequence_species_match_status` independently compare that project
  object with the explicit expected species. A known target-sequence mismatch
  or PrimerBank-record mismatch stops before PCR testing; a missing target
  organism remains visible as `unresolved` with a warning because the required
  explicit species can guide the run but cannot independently validate the
  sequence annotation. The PrimerBank record itself must always be `matched`.

Primer variant screening (implemented):

- `ScreenPrimerVariants { request, path?, evidence_dir? }` is a read-only,
  engine-owned operation. It screens every declared physical pair against one
  local VCF/VCF.gz in a single streaming pass; it does not download variant
  resources or infer genomic coordinates from historical catalog positions.
- `gentle.primer_variant_screen_request.v1` requires a candidate assembly and
  explicit current-assembly binding geometry. Each forward, reverse, or probe
  oligo has one or more 1-based closed genomic segments with strand, oligo
  offsets, and the expected increasing-reference sequence. Multiple segments
  represent junction-spanning binding without inventing an intronic span.
  Optional amplicon segments make non-oligo overlap visible.
- The variant source is either a direct VCF path plus provenance or a local
  `gentle.primer_variant_resource_manifest.v1`. A manifest pins source name,
  release, population, retrieval label, assembly, optional AF INFO key,
  optional INFO keys to retain as uninterpreted annotations, optional content
  SHA-256, and contig aliases. Relative VCF paths resolve next to the manifest.
- The operation verifies request/source/VCF assembly declarations, resolves
  contigs conservatively, and checks each overlapping VCF REF allele against
  the declared reference sequence. Incompatible assembly or reference evidence
  yields `incompatible_reference`; it is never treated as a clear screen.
- Overlap rows distinguish `primer`, `probe`, and `amplicon_only`, report
  strand-aware oligo positions and distance from the 3-prime end, and classify
  SNVs, MNVs, insertions, deletions, and complex alleles. Indels and complex
  alleles remain conservative overlap evidence without haplotype realignment.
- Missing or malformed allele frequency is `null`, never zero. With a maximum
  allowed frequency, unknown AF remains relevant. Probe overlaps follow the
  requested `relevant`, `report_only`, or `ignore` policy; amplicon-only rows
  remain descriptive.
- Setting `degenerate_rescue_minimum_frequency` enables an a posteriori
  mixed-base rescue screen. Only verified, passing, simple primer SNVs with a
  known allele frequency at or above that threshold qualify. Filtered calls,
  missing frequency, indels/MNVs, probe/amplicon-only variants, and
  incompatible REF evidence remain visible with an explicit ineligibility
  reason. Changes in the configured critical 3-prime window additionally
  require `allow_critical_3prime_degenerate_rescue=true`.
- `gentle.primer_variant_degenerate_rescue.v1` records the original and
  adjusted forward/reverse sequences, strand-aware primer-oriented alleles,
  per-position IUPAC codes, contributing variant provenance, synthesis-mixture
  complexity, and a new sequence-derived `pair_id`. The IUPAC sequence denotes
  a mixed oligo synthesis, not a genotype or one heterozygous molecule, and the
  new physical pair requires fresh specificity, thermodynamic, and experimental
  validation. Retained VCF annotations such as locally supplied splice
  consequences are context only; GENtle does not infer a splicing effect.
- Adjusted IUPAC primers can be screened again against the same assembly-aware
  geometry. The cDNA PCR/qPCR matcher evaluates IUPAC symbols as nucleotide
  sets (the semantic equivalent of regular-expression character classes), so
  each represented allele can match without expanding the oligo into every
  sequence combination.
- `gentle.primer_variant_screen.v1` wraps deterministic, source-fingerprinted
  `gentle.primer_variant_evidence.v1` reports. Identical forward/reverse
  sequences share one physical `pair_id` only when their binding geometry is
  identical; all candidate-source rows remain attached. `evidence_dir` writes
  one directly consumable evidence JSON per pair for
  `BuildExperimentalAssayHandoff`.

Primer-design shell command family (implemented):

- Shared-shell family:
  - `primers design REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers design-qpcr REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers design-group-target REQUEST_JSON_OR_@FILE [--path OUTPUT.json] [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers primerbank search QUERY [--by gene-symbol|gene-id|genbank|protein|primerbank-id|keyword] [--species human|mouse|all] [--html SAVED.html] [--path OUTPUT.json]`
  - `primers primerbank show PRIMERBANK_ID [--species human|mouse|all] [--html SAVED.html] [--path OUTPUT.json]`
  - `primers primerbank test-cdna SEQ_ID FEATURE_ID PRIMERBANK_ID --species human|mouse [--html SAVED.html] [--transcript-id ID] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--path OUTPUT.json] [--svg OUTPUT.svg]`
  - `primers import-external-pairs INPUT.json|tsv SEQ_ID FEATURE_ID [--format auto|json|tsv] [--report-id ID] [--transcript-id ID] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--specificity-target-genome GENOME_ID] [--specificity-catalog PATH] [--specificity-cache-dir DIR] [--artifact-output-dir DIR] [--materialize-products] [--product-gel-ladder NAME ...] [--path OUTPUT.json]`
  - `primers screen-variants REQUEST_JSON_OR_@FILE [--path OUTPUT.json] [--evidence-dir DIR]`
  - `primers test-cdna-pcr SEQ_ID FEATURE_ID --forward SEQ --reverse SEQ [--transcript-id ID] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--path OUTPUT.json] [--svg OUTPUT.svg] [--materialize-products] [--product-output-prefix PREFIX] [--product-gel-svg OUTPUT.svg] [--product-gel-ladder NAME ...]`
  - `primers test-cdna-qpcr SEQ_ID FEATURE_ID --forward SEQ --reverse SEQ --probe SEQ [--transcript-id ID] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--path OUTPUT.json] [--svg OUTPUT.svg] [--materialize-products] [--product-output-prefix PREFIX] [--product-gel-svg OUTPUT.svg] [--product-gel-ladder NAME ...]`
  - `primers transcript-qpcr-panel SEQ_ID FEATURE_ID SHARED_QPCR_REPORT_ID [--path OUTPUT.json]`
  - `primers design-transcript-assay-panel OPERATION_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers design-transcript-assay-panel SEQ_ID FEATURE_ID [--assay-kind endpoint-rt-pcr|sybr-qpcr|taqman-qpcr] [--cdna-synthesis oligo-dt|random-hexamers|gene-specific|mixed] [--objective pan-transcript|one-per-class|minimal-discrimination-panel|isoform-end-matrix] [--coverage-policy require-all|best-effort] [--coverage-universe JSON_OR_@FILE] [--assay-tier routine-common-region-screen|isoform-discrimination|long-range-structure-discovery] [--preferred-min-amplicon-bp N --preferred-max-amplicon-bp N] [--junctions JSON_OR_@FILE] [--junction-evidence PATH ...] [--junction-evidence-priority required|preferred] [--min-3prime-junction-overlap-bp N] [--min-5prime-junction-overlap-bp N] [--annotation-release TEXT] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-assays-per-class N] [--max-mismatches N] [--require-3prime-exact-bases N] [--oligo-dt-5prime-risk-threshold-bp N] [--report-id ID] [--path OUTPUT.json] [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers experimental-handoff PANEL_REPORT_ID [--policy JSON_OR_@FILE] [--variant-evidence PATH ...] [--order-form-id ID] [--path OUTPUT.json] [--order-table OUTPUT.tsv] [--gel-svg OUTPUT.svg] [--gel-ladder NAME ...]`
  - `primers test-cdna-qpcr-fasta CDNA_FASTA[.gz] [CDNA_FASTA[.gz] ...] --forward SEQ --reverse SEQ --probe SEQ [--transcript-id ID] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--path OUTPUT.json] [--svg OUTPUT.svg]`
  - `primers preflight [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers prepare-restriction-cloning REQUEST_JSON_OR_@FILE`
  - `primers seed-restriction-cloning-handoff PRIMER_REPORT_ID VECTOR_SEQ_ID [--pair-rank N] [--mode single_site|directed_pair] [--forward-enzyme NAME] [--reverse-enzyme NAME] [--forward-leader SEQ] [--reverse-leader SEQ]`
  - `primers restriction-cloning-vector-suggestions SEQ_ID`
  - `primers list-restriction-cloning-handoffs`
  - `primers show-restriction-cloning-handoff REPORT_ID`
  - `primers export-restriction-cloning-handoff REPORT_ID OUTPUT.json`
  - `primers seed-from-feature SEQ_ID FEATURE_ID`
  - `primers seed-from-splicing SEQ_ID FEATURE_ID`
  - `primers seed-qpcr-from-feature SEQ_ID FEATURE_ID`
  - `primers seed-qpcr-from-splicing SEQ_ID FEATURE_ID [--mode shared_gene|distinguish_transcript] [--transcript-id ID] [--specificity-evidence junction_only|unique_exon_or_chain|either_prefer_junction]`
  - `primers list-reports`
  - `primers show-report REPORT_ID`
  - `primers export-report REPORT_ID OUTPUT.json`
  - `primers list-qpcr-reports`
  - `primers show-qpcr-report REPORT_ID`
  - `primers export-qpcr-report REPORT_ID OUTPUT.json`
  - `primers list-transcript-assay-panels`
  - `primers show-transcript-assay-panel REPORT_ID`
  - `primers export-transcript-assay-panel REPORT_ID OUTPUT.json`
  - `primers list-transcript-assay-fallbacks`
  - `primers show-transcript-assay-fallback EXECUTION_ID`
  - `primers export-transcript-assay-fallback EXECUTION_ID OUTPUT.json`
  - `primers oligo-order create REQUEST_JSON_OR_@FILE`
  - `primers oligo-order from-primer-report REPORT_ID --pair-rank N[,N...] [--form-id ID] [--scale TEXT] [--purification TEXT] [--modification TEXT ...]`
  - `primers oligo-order from-qpcr-report REPORT_ID --assay-rank N[,N...] [--include-probe true|false] [--form-id ID] [--scale TEXT] [--purification TEXT] [--modification TEXT ...]`
  - `primers oligo-order list`
  - `primers oligo-order show FORM_ID`
  - `primers oligo-order export FORM_ID OUTPUT.json`
  - `primers oligo-order route FORM_ID`
  - `primers oligo-order quote FORM_ID [--provider metabion|geneart] [--service-kind KIND] [--output-dir DIR]`
  - `primers oligo-order review-dedup FORM_ID [--reviewer NAME] [--duplicate-action keep-separate] [--note TEXT]`
- `primers design` expects an operation payload whose root variant is
  `{"DesignPrimerPairs": {...}}`.
- `primers design-qpcr` accepts either:
  - an operation payload whose root variant is `{"DesignQpcrAssays": {...}}`
  - a full `gentle.qpcr_seed_request.v1` payload carrying one runnable
    `operation.DesignQpcrAssays`

Related-sequence group-target primer design:

- `DesignPrimerGroupTarget` and `primers design-group-target` implement the
  local, deterministic counterpart of Primer-BLAST's related-sequence group
  target mode. The request supplies at least two loaded linear sequence ids,
  side/pair constraints, product limits, a strict `require_all` or explicit
  `best_effort` coverage policy, and bounded alignment/search budgets.
- Unless explicitly selected, the longest member is the representative with a
  stable sequence-id tie break. GENtle globally aligns every other member to
  that representative, records the alignments and member sequence SHA-256s,
  and derives exact conserved representative intervals. A member insertion
  splits an interval: individually matching bases on either side are never
  treated as one contiguous primer footprint.
- Optional `target_start_0based` and
  `target_end_0based_exclusive` constrain a core that the pair must flank.
  Without a target, GENtle evaluates bounded left/right combinations of the
  conserved intervals. Long intervals are split into overlapping records so
  every possible primer footprint remains represented.
- Candidate counts use the same local side filters and pair geometry as the
  internal designer. Alignment cells, record count, per-record pair work, and
  total pair work are checked before Primer3. An over-budget request fails as
  `group_target_alignment_space_too_broad` or `search_space_too_broad`; raising
  a limit is an explicit request change.
- Each retained candidate is tested against every supplied member through the
  shared cDNA product scanner. `require_all` fails unless one pair yields
  exactly one in-range product on every member. `best_effort` never hides
  uncovered or multi-product members; it returns `completion_status = partial`
  and the complete pair-by-member matrix.
- Output schema `gentle.primer_group_target_design.v1` binds the normalized
  request digest, representative and member sequence digests, alignments,
  exact common intervals, bounded search records, exact Primer3 Boulder input
  when Primer3 is used, accepted pairs, member products, backend provenance,
  and warnings. It is returned in
  `OpResult.primer_group_target_design`; `--path` writes the same report.
- Placement intervals remain exact across the supplied group even when
  mismatch-tolerant product evaluation is explicitly requested. This avoids
  implying that a tolerated post-design mismatch was a conserved primer site.
- The operation consumes caller-provided sequences only. It does not retrieve
  related records, infer orthology/paralogy, or replace final genomic and
  whole-cDNA specificity assessment.

Primer-pair alternative frontier:

- New `PrimerDesignReport` and `TranscriptAssayPanelReport` projections carry
  an additive `pareto_frontier`. It compares only candidates that already pass
  hard design constraints and retains non-dominated tradeoffs across the
  dimensions available to that design: product coverage/ambiguity,
  annotation-common-region evidence, practicality tier, existing candidate
  score, Tm delta, and pair/3-prime complementarity.
- The existing score and declared panel objective still select the primary
  assay. The frontier explains credible alternatives; it does not rescue a
  failed pair or replace specificity, variant, repeat, or wet-lab review.
- Frontier comparison is deterministically capped at 2,000 accepted
  candidates and output at 25 alternatives. A capped report uses
  `status = bounded_non_dominated_accepted_candidate_projection`, records both
  accepted and evaluated counts, and sets `truncated = true`.
- Older report JSON without these additive fields remains readable and yields
  an empty frontier.
- External primer-pair import contracts:
  - JSON input schema: `gentle.external_primer_pair_batch.v1`. TSV uses the
    columns `source_kind`, `provider`, `catalogue_id`, `source_url`,
    `claimed_accession`, `aliases`, `forward_sequence_5_to_3`,
    `reverse_sequence_5_to_3`, `claimed_target`, `validation_claims`, and
    `annotations_json`
  - `source_kind` is `external`, `commercial_catalogue`, `literature`, or
    `laboratory`. Provider is required for every row; commercial-catalogue rows
    additionally require `catalogue_id`
  - sequence normalization removes whitespace and copied position digits,
    uppercases bases, maps RNA `U` to DNA `T`, and then rejects any
    non-IUPAC character with row, role, and input-character position
  - sequence-derived primer ids depend only on the normalized oligo sequence.
    The oriented pair id depends only on normalized forward plus reverse sequence.
    Duplicate pair rows are evaluated once but retain all source provenance,
    aliases, claimed accessions, validation claims, source-file SHA-256, and
    deterministic source-record ids. The report also carries a normalized,
    row-order-independent batch SHA-256; automatic report ids bind that digest
    to the cDNA/specificity/materialization settings so different evaluations
    cannot silently overwrite one another
  - output schema: `gentle.external_primer_pair_import_report.v1`, also returned
    in `OpResult.external_primer_pair_import_report`. Every unique pair carries
    GENtle-computed length, Tm, Tm method/assumptions, GC fraction/percent,
    delta-Tm, 3-prime clamp, homopolymer/self-complementarity metrics, shared
    oligo/inter-oligo QC, cDNA transcript products and map, genomic-DNA carryover
    assessment, specificity state, optional product materialization/gel, and all
    retained source rows
  - `tm_c` is primer melting temperature from GENtle's shared primer Tm model,
    not a PCR annealing temperature. No cycling condition is inferred
  - provider targeting and validation text is provenance only. It never changes
    cDNA coverage, carryover, QC, or specificity results, and
    `vendor_claims_used_as_biological_evidence` is always false
  - specificity is explicitly `not_run` without
    `--specificity-target-genome`; vendor claims cannot turn that state into a
    pass. When specificity is requested, GENtle derives the intended
    transcript set and product ranges from its computed cDNA assay, then
    projects genomic target geometry only when the project sequence has a
    provenance-bearing genome anchor. A complete search can therefore pass on
    computed geometry, but never on provider claims. Product materialization
    and gel rendering remain opt-in because they create first-class project
    products
- `primers test-cdna-pcr` and `primers test-cdna-qpcr` are non-mutating assay
  checks over transcript-derived cDNA templates and return
  `gentle.cdna_assay_test_report.v1`; `--path` persists that same report and
  `--svg` writes the embedded transcript-map SVG. `--transcript-order` controls
  transcript-map/report row order for transcript-derived tests only.
  `--map-coordinate-mode` selects the map axis: `cdna` keeps each transcript on
  its own cDNA-length axis, while `genomic_aligned` draws exon blocks,
  primer/probe glyphs, and amplicon source spans on the shared source/genomic
  coordinate range so common primer loci line up across isoforms. The returned
  report carries `transcript_map_coordinate_mode`, and the embedded
  `transcript_map` carries `coordinate_mode`; in genomic-aligned SVGs dashed
  product bridges mark spliced products whose source ranges skip genomic
  sequence while amplicon labels remain cDNA lengths. Shell output includes
  `preferred_artifacts[]` when `--svg` is supplied so ClawBio can promote the
  same map through its PNG-first artifact bundle.
- Supplying `--materialize-products` or `--product-gel-svg` to
  `primers test-cdna-pcr` / `primers test-cdna-qpcr` switches the command from
  report-only to state-changing: detected products are materialized as derived
  DNA sequences, grouped into a product container, and optionally rendered as a
  gel. Shell output then includes `materialization` plus gel-first
  `preferred_artifacts[]` so ClawBio can show nonspecific products as bands
  before the transcript map.
- `primers transcript-qpcr-panel` is a non-mutating compatibility helper for
  fixed-component transcript-targeted qPCR panels. It consumes a stored shared-gene qPCR report,
  reuses that assay's shared reverse primer and probe, and emits
  `gentle.transcript_qpcr_panel.v1` with shared oligo records plus one
  transcript row per admitted cDNA template. Rows prefer a transcript-specific
  exon-junction forward primer when possible, allow a single-exon/exon-chain
  characteristic forward primer when that is the valid/efficient evidence, and
  emit deterministic `not_found` rows when no single forward primer can
  distinguish the transcript while retaining the shared reverse/probe product.
  Byte-identical mature cDNAs instead receive
  `not_distinguishable_between_members`; this statement is never inferred for
  merely similar or difficult-to-prime transcripts. New multi-transcript work
  should prefer `primers design-transcript-assay-panel`.
  All source positions are machine-readable as local 0-based/exclusive ranges
  plus 1-based/inclusive display ranges; genomic 1-based/inclusive ranges and
  reference strand labels are included when a genome anchor is available.
- Transcript assay panel schema:
  - `gentle.transcript_assay_panel.v2`
  - `gentle.transcript_assay_panel_feasibility.v1` is the read-only,
    pre-Primer3 endpoint-matrix assessment returned by
    `primers inspect-transcript-assay-feasibility`. It binds the exact
    `DesignTranscriptAssayPanel` operation SHA-256 and carries stable
    first/terminal class and reaction ids, transcript-local end windows,
    effective primer/product limits, proven structural blockers, and the count
    of reactions for which Primer3 is still warranted. Workload fields retain
    the warranted total/max template bp and candidate-pair request upper bound
    so a scheduler can apply an explicit timeout policy without GENtle claiming
    an exact duration. It proves only certain
    annotation/geometry negatives; `primer_search_required` is not a promise
    that sequence composition or thermodynamics will yield a pair. The MCP
    `transcript_assay_feasibility` tool delegates to this same shared shell and
    engine path
  - feasibility and completed panel reports may embed
    `gentle.transcript_assay_primer_search_plan.v1`. This sequence-aware plan
    is produced before Primer3 and binds the exact operation digest, effective
    `search_policy`, transcript-local allowed intervals, rejected starts and
    reasons, exon/junction context, quality-filtered candidate estimates, a
    conservative placement upper bound, and the exact bounded Primer3 fields.
    The default policy offers no primer side more than 96 bp in one record (or
    the caller's larger explicit minimum primer length for legacy requests) and
    caps per-record, per-target, and operation-wide estimated pair work. A
    required target with no admissible interval is `no_admissible_regions`; a
    declared budget that cannot contain the search is
    `search_space_too_broad`. Both stop before Primer3. Increasing a budget is
    an explicit request change and never changes `coverage_policy`. Required
    endpoint/junction targets receive budget before GENtle-generated optional
    anchors; optional anchors may be omitted with a warning, but final assay
    coverage is still evaluated normally
  - `gentle.transcript_assay_cdna_similarity_map.v1` is an optional,
    content-bound planning input referenced by
    `search_policy.cdna_similarity_map = {path, expected_sha256}`. It records
    transcript-local 0-based/exclusive intervals supported by a whole-cDNA
    similarity search, their intended-family/paralog/other-subject
    classification, source subject ids, map and target-resource ids, database
    content fingerprint, and exact BLAST program/task/version/options
    provenance. Every interval also binds the SHA-256 of the mature-cDNA
    template whose coordinates it uses. Active transcript ids with a changed
    template digest or out-of-bounds interval fail before Primer3; unrelated
    transcript rows may coexist in one multi-transcript map
  - `BuildTranscriptAssayCdnaSimilarityMap` and
    `primers build-transcript-assay-cdna-similarity-map` are the engine-owned
    producer. They query each byte-distinct mature-cDNA class against an
    already prepared, fingerprinted transcriptome/cDNA BLAST resource, retain
    exact BLAST invocation provenance, and classify only identities established
    by annotation or caller-supplied paralog ids. Similarity alone never creates
    a paralog claim. A design may request the producer with
    `search_policy.build_cdna_similarity_map`; the generated report is embedded
    and bound before Primer3. The read-only feasibility route does not execute
    BLAST and says so explicitly
  - the v1 map disposition is deliberately advisory: `informative` leaves the
    search order unchanged and `deprioritize` adds overlap evidence to bounded
    records. Within each biological target, records with less advisory overlap
    run first, followed by deterministic candidate-count and record-id
    tie-breakers. The map does not emit `SEQUENCE_EXCLUDED_REGION`, suppress a
    primer window, claim specificity, or waive complete-cDNA/genomic BLAST.
    Search plans echo the map id, file digest, database fingerprint, affected
    interval ids, and overlap bp. Callers may bind a prebuilt map or request
    the engine-owned producer; these are mutually exclusive
  - `search_policy.interval_evidence[]` unifies caller-supplied variation,
    repeat/low-complexity, and similarity intervals. `exclude` rows become
    deterministic Primer3 `SEQUENCE_EXCLUDED_REGION` fields;
    `deprioritize` rows affect bounded-record ordering only; `informative`
    rows remain visible without changing selection. The optional project
    projection maps current known-variant and RepeatMasker features into each
    mature-cDNA template with strand-aware, content-checked coordinates
  - `search_policy.primer3_chemistry` records either executable defaults, an
    explicit Primer3-2.x computational baseline, or reviewed custom overrides.
    Overrides are restricted to documented salt, dNTP, DNA, DMSO, formamide,
    and Tm-model tags and are copied into every emitted Boulder record. A
    profile is calculation provenance, not a claim about the laboratory mix
  - `search_policy.max_junction_single_side_match_bp` optionally caps either
    contiguous side of a realized junction-spanning primer. It complements
    the existing 3-prime/5-prime minimum overlaps and is enforced identically
    by preflight candidate counting, internal design, and Primer3 post-filtering
  - `search_policy.min_intervening_intron_bp` optionally requires a retained
    pair to flank an annotated intron of at least that length. A genuinely
    junction-spanning primer also satisfies the anti-genomic-carryover design
    intent because it has no contiguous genomic binding site. The rule is a
    design constraint, not a genomic specificity pass, and is enforced in the
    same preflight/internal/Primer3-post-filter paths
  - `DesignTranscriptAssayPanel.search_policy` is additive and optional. Old
    requests use the documented default. During a progress-capable Primer3
    run, progress rows retain native bounded work counters and add the GENtle
    search target/record ids, record ordinal/count, deterministic pair-work
    bounds, elapsed time, and a separately labelled uncertain linear runtime
    projection. After ten minutes, a projection above five hours is marked
    suspicious; after thirty minutes, a projection above ten hours stops only
    that bounded record. Five hours without usable progress also stops the
    record. Other precomputed records may continue, but strict final coverage
    still fails when the stopped record leaves a required class/reaction
    uncovered; that strict failure carries the structured runtime-reduction
    receipts, including the stopped child's output byte counts and hashes.
    Malformed future `PRIMER3_PROGRESS` rows are retained in stderr
    provenance and surfaced in the Primer3 explanation rather than silently
    treated as valid counters
  - new reports persist `source_genome_anchor` at design time, including the
    source sequence, genome/reference id, chromosome, one-based inclusive
    interval, strand, and verification state. Whole-panel genomic specificity
    therefore remains reproducible after live project provenance changes
  - old v2 reports without `source_genome_anchor` remain readable. Finalization
    may use current project provenance as a compatibility fallback; if no
    trustworthy target geometry is available, a complete external search is
    reported as `not_assessed`, never as pass or biological failure
  - the shell command's operation-JSON form accepts the externally tagged
    `DesignTranscriptAssayPanel` operation unchanged. It therefore exposes the
    complete engine request, including side/pair constraints and optional
    probe fields that are intentionally not repeated as dozens of convenience
    flags. The same object is accepted by workflows, MCP `op`, JavaScript
    `apply_operation`, and Lua `apply_operation`
  - `DesignTranscriptAssayPanel` derives mature cDNA templates from the shared
    splicing engine, groups only byte-identical sequences into exact
    equivalence classes, generates mode-appropriate primer candidates from one
    representative per class, and confirms candidates with the mismatch-aware
    cDNA assay path.
  - `coverage_universe` is a third, independent request axis. It does not
    replace `objective` or `coverage_policy`: it selects the mandatory
    biological targets, the objective chooses how to assay them, and the
    policy decides whether incomplete coverage is accepted. Omission uses
    `all_annotated_cdna_classes` and is not serialized, preserving historical
    operation bytes and approval/checkpoint digests. Other kinds are
    `explicit_transcripts` and `uniprot_supported_isoforms`.
  - UniProt-supported coverage is derived from one or more content-bound
    `gentle.uniprot_projection_audit.v1` sources. New audits carry the complete
    parsed UniProt `ALTERNATIVE PRODUCTS` inventory. Each mandatory target is
    `(UniProt entry id, named isoform id)`, not an Ensembl transcript row. An
    entry without named alternative products contributes one canonical protein
    target. Every mapped Ensembl/current-annotation transcript remains grouped
    under its protein target, while GENtle deterministically chooses one
    representative assay template by audit status and transcript id. This
    prevents one protein isoform from becoming many mandatory Primer3 searches.
  - Add a content-bound `gentle.uniprot_linked_transcript_inventory.v1` source
    to make the linked-transcript denominator authoritative across primary,
    patch, and alternate loci. Resolution verifies the inventory's canonical
    content hash and requires its assembly and annotation release to match the
    assay request before Primer3 runs. Every linked record remains in the
    denominator. A record outside the loaded locus may inherit cDNA coverage
    only when its exact mature-cDNA SHA-256 matches an assessed local template;
    that identity does not assess the record's distinct genomic locus for
    carryover specificity. Missing, ambiguous, stale, or unmatched mandatory
    records stop before Primer3 rather than silently shrinking the denominator.
    Audit-only legacy requests remain executable, but their linked-transcript
    coverage summary is unavailable rather than authoritative.
  - Successful inventory resolution now persists one
    `linked_transcript_resolutions[]` row per linked record. Each row keeps the
    exact versioned source transcript id and mature-cDNA SHA-256, records the
    digest-matched local transcript ids, and joins that digest to exactly one
    assessed panel equivalence group. Local-locus identity compares normalized
    Ensembl stable ids while retaining exact versions in provenance; digest-
    identical records with another stable id remain genomically unassessed.
  - Normalization records each audit SHA-256 plus the exact protein-target
    inventory in `required_uniprot_isoforms`; execution requires those hashes
    and rejects changed content. An old audit without the inventory must be
    regenerated. A named isoform with no current transcript mapping is
    `unresolved` and stops before Primer3 even under `best_effort`. Reports
    retain every mapped transcript and exact-cDNA class, the representative
    template, coverage assay ids, and whether the selected assays distinguish
    each protein target from the others. Shared/indistinguishable coverage is
    reported but is acceptable when every mandatory protein target is covered.
  - UniProt identity, Ensembl cross-reference mapping, current GenBank/Ensembl
    transcript geometry, and CDS/peptide audit status remain separate evidence
    layers because those annotation sources need not specify identical coding
    models. PCR measures mapped mature cDNA and does not establish protein-
    isoform expression or repair a CDS disagreement. Use a separate
    all-annotation plan when transcript-model coverage itself is the objective.
  - This additive axis deliberately preserves alternative plans. A user may
    keep separate named all-annotation, UniProt-prioritized, and explicit
    transcript plans, each with any existing objective, rather than replacing
    one plan with a supposedly universal answer.
  - `assay_kind` is `endpoint_rt_pcr`, `sybr_qpcr`, or `taqman_qpcr`.
    Requests and older reports that omit it retain the original
    `taqman_qpcr` forward/reverse/probe behavior. Endpoint and SYBR records
    always carry the canonical primer pair and leave `probe` plus the legacy
    TaqMan `assay` object absent.
  - objectives are `pan_transcript`, `one_per_class`,
    `minimal_discrimination_panel`, `maximally_informative_panel`, and
    endpoint-only `isoform_end_matrix`.
    The endpoint objective derives first-exon/first-junction and terminal-exon/
    last-junction classes, retains only combinations supported by an annotated
    mature transcript, and permits one physical primer pair to reference more
    than one supported end reaction.
  - `maximally_informative_panel` requires
    `informative_selection = gentle.transcript_assay_informative_selection_policy.v1`.
    Its whole-panel `max_assays` budget is applied after required evidence
    assays, then a deterministic bounded-greedy selector maximizes equal-weight
    target coverage, observable target-pair separation, and existing assay
    preference in that order. Separation may be presence/absence or predicted
    single-product size, but a size distinction counts only at or above
    `minimum_resolvable_size_difference_bp`. The report carries
    `bounded_greedy_v1`, never a global-optimum claim, plus recomputable
    incremental/leave-one-out counts, redundant alternatives, and unresolved
    target pairs. Long-range candidates and zero-increment assays remain
    excluded unless the explicit policy allows them.
  - With `uniprot_supported_isoforms`, the existing objective names operate on
    mandatory protein targets: pan coverage reaches every target, one-per-class
    selects coverage per protein target, discrimination tries to separate
    target detection patterns, and endpoint reactions derive from one recorded
    representative per target. Exact-cDNA classes and every mapped transcript
    remain in the report. Consequently `completion_status=complete` can coexist
    with informational uncovered transcript classes; target completion is
    stated separately in `uncovered_coverage_target_ids` and each target's
    `coverage_status`. Failure to distinguish two covered targets is recorded,
    not silently promoted to an uncovered target.
  - `assay_tier` is an independent experimental-purpose axis:
    `routine_common_region_screen`, `isoform_discrimination`, or
    `long_range_structure_discovery`. It does not replace the selection
    objective. The routine tier requires `pan_transcript` and confirms a
    common region only from the intersection of transcript annotation source
    intervals; product detection and Clariom intensity cannot create that
    structural claim.
  - `practicality_policy` records a preferred product range inside the allowed
    `min_amplicon_bp..max_amplicon_bp` range. Objective-specific biological
    coverage is considered first, preferred routine length second, and the
    existing primer candidate score last. A selected product below the
    preferred minimum is `allowed_nonpreferred`; only a product above the
    preferred maximum is `long_range_fallback`. No universal preferred cutoff
    is invented when the caller does not supply one.
  - each pair summary carries a concise `selection_explanation` and at most
    five deterministic `considered_alternatives[]`. PSR and JUC inputs remain
    separate `selection_evidence[].evidence_kind` rows: JUC can constrain a
    junction target, while overlapping PSR evidence is contextual support and
    never proof that an exon is common.
  - endpoint mode defaults to `200..10000` bp and refuses a configured ceiling
    above 10,000 bp. Its `end_classes[]`, `end_reactions[]`, and
    `band_size_matrix[]` make differently sized transcript products explicit.
    Endpoint-gel band intensity is rough or semi-quantitative, not a
    quantitative transcript-abundance measurement. Endpoint searches default
    to four requested candidate pairs per annotated end reaction when
    `max_assays_per_class` is omitted; an explicit value remains authoritative
  - SYBR mode defaults to short products and never fabricates an internal
    probe. `short_sybr_junction_assays[]` is the primer-only subset whose
    selected forward or reverse primer satisfies a requested junction overlap.
  - `junctions[]` accepts transcript-local boundaries, genomic inclusive intron
    spans, or one-based adjacent exon ordinals. `junction_evidence_paths[]`
    consumes `gentle.probe_region_evidence_interpretation.v2` reports and turns
    each JUC row with junction geometry into an explicit audited junction target;
    raw array activity remains design evidence, not isoform validation.
  - requested junctions bypass the six-anchor automatic-search cap. Each
    request receives a `junction_evaluations[]` row with resolved transcripts,
    local positions, selected assay ids, and a reason when it cannot be
    assayed. `required` targets participate in strict completion; `preferred`
    targets are reported but non-blocking.
  - Primer3 requests use `SEQUENCE_OVERLAP_JUNCTION_LIST`,
    `PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION`, and
    `PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION`. Junction-driven requests omit
    `SEQUENCE_TARGET`, because Primer3 requires pairs to flank a target and
    that would conflict with requiring a primer to overlap the same boundary.
    The internal backend enforces the same overlap geometry before ranking,
    and both paths are independently checked when the report records
    `forward`, `reverse`, or `neither`.
  - coverage policy defaults to `require_all`. An unsatisfied strict request is
    an error that enumerates uncovered equivalence classes and persists no
    report. `best_effort` must be selected explicitly and returns
    `completion_status = partial`, warnings, uncovered class ids, and unresolved
    class pairs where applicable. For endpoint matrices, hard annotation/window
    blockers reject `require_all` before Primer3. Explicit `best_effort` skips
    those reactions while retaining `structurally_impossible` status and
    blocker details; a feasible search that finds no acceptable pair instead
    records `primer_search_exhausted`. Cancellation/timeout remains an
    execution error and is not reclassified as either scientific outcome
  - `DesignTranscriptAssayPanelWithFallback` is the approval-bound exception
    to asking again after a strict run. It executes the nested `require_all`
    operation unchanged and derives a separate `best_effort` /
    `maximally_informative_panel` operation only from the machine-typed
    `transcript_assay_coverage_infeasible` result after bounded search. Invalid
    input, unresolved UniProt mappings, absent evidence or backend, search-
    budget refusal, specificity failure, cancellation, digest drift, and
    internal errors never trigger it. Legacy requests and
    `fallback_submission.mode=never` retain the original fail-closed behavior.
    A persisted `gentle.transcript_assay_fallback_execution.v1` binds the
    strict operation digest/failure identity, policy digest/schema, search
    policy digest, engine revision, fallback operation digest, and exact field
    diff. For ordinary repository builds, the engine and selection-audit
    revision is the package version plus the full Git `HEAD` commit; source
    archives without Git retain the package version alone. This identifies the
    committed source baseline and is not a clean-working-tree attestation. The
    strict outcome is not overwritten. Any generated fallback panel
    is forcibly `completion_status=partial`, keeps all uncovered targets
    visible, and carries one deterministic combined virtual gel or the typed
    `no_predicted_products` reason. This does not imply specificity, validation,
    UniProt completeness, experimental success, or order readiness.
  - the detection matrix records `no_product`, `single_product`, or
    `multiple_products` for every selected assay and transcript row. Transcript
    interpretation is typed as `specific`, `shared_family`, `no_product`, or
    `not_distinguishable_between_members`.
  - every new `selected_assays[]` row carries
    `primer_pair_summary` (`gentle.primer_pair_summary.v2`), an additive
    communication projection assembled from the canonical pair, detection
    matrix, junction, specificity-followup, and backend records. It repeats the
    assay id, design transcript, forward/reverse sequence (explicitly
    5-prime-to-3-prime), oligo and annealing lengths, `tm_c`, GC fraction and
    unrounded percent, binding positions, canonical designed-amplicon
    coordinates/length, pair `tm_delta_c`, predicted transcript products/sizes,
    concise oligo-QC status/reasons, junction matches,
    `whole_genome_specificity_status`, GENtle package version, requested/used
    backend, and optional Primer3 version. `length_nt` must equal the returned
    sequence length, and `tm_delta_c` is copied from and checked against the two
    canonical melting temperatures; report consumers must not recompute Tm.
    The QC block interprets the pair's stored rule flags and metrics; summary
    generation does not rerun sequence or thermodynamic analysis. Compatible
    older reports are enriched when read/exported. If a legacy payload omitted
    boolean pair-rule flags, serde defaults them to `false`; GENtle cannot
    safely infer whether those rules failed or were never recorded, so a rerun
    is required before treating that enriched QC block as a historical result.
  - `design_amplicon_start_0based`,
    `design_amplicon_end_0based_exclusive`, and
    `design_amplicon_length_bp` are copied verbatim from the selected canonical
    pair on `design_transcript_id`. They describe the amplicon used during pair
    design even if the later cross-transcript matrix calls no product.
    `predicted_amplicon_lengths_bp` is a separate sorted, deduplicated union of
    product lengths in that matrix and may therefore be empty.
  - `gc_percent` is the unrounded convenience projection
    `gc_fraction * 100.0`. Primer binding coordinates are strand-agnostic,
    zero-based half-open footprints on the design-transcript cDNA; for the
    reverse primer they do not encode the 5-prime-to-3-prime order of
    `sequence_5_to_3`.
  - summary v2 keeps machine identity separate from human naming:
    - `assay_id` remains the sequence-derived pair/probe identity;
      `forward.primer_id` and `reverse.primer_id` are stable hashes of the
      respective primer sequences.
    - `display_label` values are annotation-context labels such as
      `GENE_E2_F`, `GENE_E2|E3_R`, and `GENE_E2F-E6R`. They use exon ordinals
      from `exon_numbering_reference_transcript_id` in transcript 5-prime-to-
      3-prime order and may change when the annotation reference changes.
      They never replace the immutable ids.
    - `aliases[]` preserves optional literature/lab names when a caller has
      supplied them. GENtle does not derive or fabricate aliases for de-novo
      designs. `origin` distinguishes `de_novo`, `legacy_literature`, and
      `legacy_lab`; `selection_role` independently permits `anchor` or
      `companion` without conflating either with origin.
    - `satisfied_design_objective` and structured `selection_reasons[]` rows
      explain the panel-level reason for choosing the pair. Each reason has a
      typed `code`, readable `message`, and `related_ids[]`; callers do not
      need to parse one composite provenance string. The separate booleans
      `forward|reverse.primer_spans_junction`, `amplicon_spans_junction`, and
      `selected_because_of_junction_evidence` must not be treated as synonyms.
  - `selection_evidence[]` is a repeatable structured evidence table. A
    Clariom/JUC-derived row records `probe_region_influenced`, its
    required/preferred status, platform, PSR/JUC feature id, genomic region,
    exon junction, contrast, available statistic/intensity source, and source
    schema/path/SHA-256. It never records `probe_sequence_reused` unless an
    actual probe sequence has been supplied and matched. The current
    transcript-panel adapter consumes projected region geometry and therefore
    emits only `probe_region_influenced`.
    - Differential-junction eligibility is narrower than junction geometry or
      `required` priority. `absolute_effect_threshold` and
      `differential_eligibility` state whether the upstream interpretation
      report supplied a contrast, numeric effect, and explicit absolute-effect
      threshold that the JUC row passed. Without all three, eligibility is
      `not_assessed`, `missing_contrast`, or `missing_measurement`; GENtle does
      not infer differential transcript contribution merely because a Clariom
      junction was supplied.
    - Each selected pair also reports its rank-ordered
      `incremental_binary_distinction_count`, leave-one-out
      `exclusive_binary_distinction_count`, newly separated exact-cDNA class
      pairs, and assays with an identical binary detection signature. A
      threshold-qualified required JUC can therefore be marked
      `retained_despite_zero_marginal_discrimination` while its selection
      reason names the contrast, descriptive effect, threshold, and separate
      validation role. This is not a statistical-significance or isoform-
      abundance claim.
    - `disposition` separates preferred geometry, preferred differential
      evidence, required validation obligations, below-threshold rows, and
      incomplete provenance. Preferred evidence can influence ranking but
      never masquerades as an obligation. Observation, transcript-projection,
      and matched-target counts are reported separately.
    - `selection_audit_status` defaults to `not_computed_legacy`. Current
      designs and stored reports with canonical detection matrices receive a
      deterministic `binary_detection_leave_one_out_v1` recomputation,
      generator revision, and stable coverage-target identity when available;
      absent historical values are never presented as computed zeroes. The
      pair-level `selection_operation_sha256` binds the audit to the normalized
      panel-design request; an empty legacy digest remains missing provenance.
  - `selection_provenance_status = de_novo_no_external_selection_evidence`
    makes an evidence-free design explicit; empty aliases or evidence rows are
    not retroactive evidence. Reports predating v2 are refreshed with
    `legacy_report_selection_provenance_unavailable` and
    `geometry_not_persisted_in_legacy_report_re_run_required`, preserving
    sequence/Tm/product data while requiring a rerun for trustworthy labels or
    selection provenance. For such legacy reports, requested primer-junction
    overlap booleans can be restored from canonical `junction_matches[]`, but a
    remaining `false` `primer_spans_junction` or `amplicon_spans_junction`
    value is not authoritative for unrequested junctions because complete exon
    geometry was not persisted; rerun the design before interpreting absence.
    Legacy oligo `origin` remains `unknown` rather than being inferred as
    `de_novo`.
  - `tm_c` is primer melting temperature, never a recommended PCR annealing
    temperature. This summary has no annealing-temperature field because the
    panel request does not supply a complete chemistry/polymerase model.
    `genomic_carryover_status = not_evaluated` is likewise explicit: requested
    junction overlap is reported separately and does not establish a genomic-
    template or whole-genome specificity pass. Exact per-transcript hit/exon/
    carryover geometry remains available from `gentle.cdna_assay_test_report.v1`.
    Blank legacy genomic-confirmation values are normalized to `not_run`.
  - `provenance.gentle_version` identifies the GENtle binary that generated the
    communication projection. It is not claimed to be the version that
    originally designed a legacy pair, because older reports did not persist
    that value.
  - each matrix cell also carries `oligo_dt_5prime_reach`. For oligo-dT cDNA
    and each predicted product, `required_cdna_reach_from_3prime_end_bp`
    measures from the annotated mature-transcript 3-prime end to the product's
    most 5-prime required base. This is sequence geometry, not a measured
    reverse-transcription length. Without
    `oligo_dt_5prime_risk_threshold_bp`, GENtle reports the distance as
    `distance_reported_unthresholded` and makes no categorical risk call. With
    an explicit positive threshold, product-bearing cells are classified as
    `within_configured_threshold` or `elevated_5prime_risk`; the former does
    not guarantee complete reverse transcription. Cells lacking a required
    oligo target are `structural_target_absent`, while cells with targets but
    without compatible product geometry are `indeterminate`. Non-oligo-dT
    requests remain `not_applicable`.
  - exact-match filtering is a negative-only optimization: GENtle uses it only
    when mismatch tolerance is zero and all binding sequences are unambiguous
    DNA. Every retained candidate is evaluated by the full cDNA assay path;
    mismatch-tolerant requests are never discarded for lacking an exact hit.
  - `order_ready_primers[]` contains deterministic forward/reverse rows (plus a
    probe only for TaqMan), while `specificity_followups[]` records that the
    local cDNA matrix ran and provides the existing prepared-genome BLAST route
    for genomic confirmation. It does not claim that BLAST/e-PCR confirmation
    has already run.
  - `provenance` records the annotation release, admitted transcript ids,
    primer backend, and SHA-256 identified Clariom/JUC inputs.
  - when `cdna_synthesis = oligo_dt`, endpoint reports warn that incomplete
    reverse transcription can underrepresent long 5-prime regions, that a
    missing long band is not proof of isoform absence, and that short end-
    specific assays or 5-prime RACE may be required.
  - reports are stored in the existing `gentle.primer_design_reports.v1`
    project metadata store without changing that store schema. The individual
    report carries the v2 schema above and is available through list/show/export
    shell routes.
- Gene isoform-assay study planning:
  - `primers plan-gene-isoform-study REQUEST_JSON_OR_@FILE` consumes
    `gentle.gene_isoform_assay_study_plan_request.v1` and, during planning,
    returns `gentle.gene_isoform_assay_study_plan.v1` without running primer
    design
  - `--normalize-only --normalized-request OUTPUT.json` resolves every policy
    default, sorts declared inputs, hashes the isoform-evidence report, hashes
    each optional JUC/PSR input and prior plan, and preserves observations as
    `user_supplied_unvalidated`; it also records the evidence report id/schema,
    policy digest, optional evidence report ids/schemas, and prior plan id.
    In normalize-only mode stdout is this normalized request object itself,
    not an envelope; callers may persist it verbatim. This complete normalized
    record is the first approval basis
  - the plan reports transcript count, byte-identical mature-cDNA classes,
    annotation-backed informative regions, separately labelled array support,
    abundance, exact-label responsiveness in annotation-informative regions,
    assayability, missing
    evidence, transparent decision factors, the automatic recommendation, and
    any explicit override plus reason. Raw Clariom activity is descriptive and
    is not promoted to formal differential expression
  - the request may declare one `coverage_universe`; normalization resolves
    and content-binds it before the first approval, the plan echoes a complete
    `coverage_resolution`, and every emitted assay operation carries the same
    universe. This changes only the mandatory target set. It neither chooses
    an objective nor weakens coverage. Compare alternative universes by keeping
    separate named plans rather than silently rewriting an approved plan
  - the request may also declare
    `fallback_submission.mode=preapproved_informative_partial`. Normalization
    binds the complete versioned fallback policy into `request_sha256`; only
    eligible strict one-per-class/minimal-discrimination steps are wrapped,
    and the exact wrapper bytes are then bound into `operation_batch_sha256`
    and `approved_workflow_sha256`. A blank fallback report id is derived
    deterministically from the strict report id before those second-stage
    digests are computed. Changing the assay budget, size-resolution threshold,
    or any other fallback field invalidates both approval bases.
  - selected profiles are `routine_common_region_screen`,
    `targeted_junction_validation`, `isoform_discrimination`,
    `comprehensive_isoform_dossier`, and `long_range_structure_discovery`.
    Missing evidence remains unknown and cannot by itself reduce study depth
  - `planned_operations[]` stores each complete
    `DesignTranscriptAssayPanel` payload and digest. `operation_batch_sha256`
    binds their exact order; `approved_workflow_sha256` binds the exact
    canonical workflow bytes written by `--workflow OUTPUT.json`. Endpoint
    planned operations additionally echo the geometry-only feasibility report.
    Its `operation_sha256` fingerprints the canonical computational request
    with the output-only `path` normalized away, so relocating an export does
    not change feasibility provenance. It is not an approval digest; the
    planned-operation, batch, and workflow digests continue to bind the exact
    payload including its destination path;
    `policy.endpoint_coverage_policy` defaults to `require_all`, and an
    explicitly approved `best_effort` value is carried unchanged into the
    emitted endpoint operation
  - `primers execute-gene-isoform-study-workflow PLAN_JSON_OR_@FILE
    WORKFLOW_JSON_OR_@FILE` verifies both digests before executing the already
    parsed workflow. A missing legacy workflow digest or any byte/operation
    mismatch fails before the first operation runs. Success returns
    `gentle.gene_isoform_assay_study_workflow_execution.v1` with both verified
    digests and the ordinary workflow results
  - multi-study second-stage approval is additive to that single-study route:
    `primers compose-gene-isoform-study-workflow-batch REQUEST_JSON_OR_@FILE`
    consumes `gentle.gene_isoform_assay_study_workflow_batch_request.v1`, whose
    ordered `entries[]` name one plan path and workflow path each. It verifies
    every existing plan/workflow binding and returns
    `gentle.gene_isoform_assay_study_workflow_batch.v1` with canonical paths,
    exact file hashes, workflow run ids, per-workflow operation hashes/counts,
    one combined ordered-operation hash, and a self-binding batch-basis hash.
    It is a pure read and executes no assay operation. The combined digest
    hashes the compact JSON serialization of all workflow `ops` concatenated in entry order; the
    batch-basis digest hashes compact JSON of the complete batch with
    `batch_basis_sha256` cleared
  - `primers execute-gene-isoform-study-workflow-batch BATCH_JSON_OR_@FILE
    [--checkpoint-dir DIR] [--reuse-proposal PROPOSAL_JSON_OR_@FILE
    --approve-reuse-sha256 SHA256] [--on-gene-failure abort|continue]`
    rechecks the batch-basis hash and every referenced file, plan identity,
    workflow-byte hash, operation hash/count, and combined ordered digest
    before executing the first workflow. It also evaluates every exact endpoint
    operation after those digest checks but before the first workflow: a
    structurally impossible strict operation rejects the complete batch with
    its machine-readable feasibility report, starts no Primer3 process, and
    applies no operation. Success returns
    `gentle.gene_isoform_assay_study_workflow_batch_execution.v1`. This keeps a
    multi-gene second stage inside one state-bound command instead of approving
    several commands against a shared state that the first command then
    changes. Execution occurs in a detached engine. SIGUSR1/runtime status
    names the current gene, plan, workflow and operation ordinals plus honest
    completed, failed, and remaining workflow counts
  - `--on-gene-failure` selects the unit of atomicity. Under the default
    `abort` the live project is replaced only after the complete batch
    succeeds. Under `continue` each gene is atomic instead of the batch: a
    failing gene is restored to an in-memory boundary snapshot taken when that
    gene started, so its earlier operations are never committed, and execution
    proceeds with the remaining precomputed genes. The report then carries
    `execution_atomicity: "per_gene"`, `atomic_detached_execution: false`,
    `batch_complete: false`, honest `completed_workflow_count` /
    `failed_workflow_count` / `skipped_operation_count`, a `gene_failures[]`
    entry per skipped gene (ordinals, `reason_code`, `detail`,
    `rolled_back_operation_count`), and a `status` of `completed` or `failed`
    on every `executions[]` entry. A batch in which no gene completed returns
    an error and commits nothing, exactly as `abort` does
  - before a `continue` run starts, GENtle rejects any approved operation that
    declares an external output path. Filesystem side effects cannot
    participate in the in-memory per-gene rollback; callers must omit those
    paths or use `abort`. After rollback, `last_checkpoint_manifest` is reset
    to the checkpoint present at the failed gene's boundary. Continue-mode
    checkpoints are emitted only after a whole gene succeeds, so the run does
    not add a checkpoint advertising a discarded partial-gene state. Abort
    mode continues to emit per-operation recovery checkpoints
  - `continue` relaxes execution atomicity only. Approval and verification
    failures — batch-basis, plan, workflow, operation digests, endpoint
    feasibility — still refuse the complete batch under either policy, and the
    requested coverage policy is untouched: a gene failing `require-all`
    coverage is skipped, never downgraded
  - with `--checkpoint-dir`, every successful approved operation writes a
    `gentle.gene_isoform_assay_study_checkpoint.v1` manifest and a hash-bound
    detached-engine snapshot. A later batch failure leaves the live state
    unchanged but retains those private checkpoint artifacts and exact result
    records
  - a checkpoint describes an exact contiguous ordered prefix. Under
    `continue` the first skipped gene puts a hole in the completed set, so
    checkpoint writing freezes at the last contiguous prefix and the report
    records `reuse_checkpoint.frozen`,
    `reuse_checkpoint.frozen_at_completed_operation_count`, and
    `reason_code: "gene_failure_broke_exact_prefix"` rather than advancing a
    prefix that no longer exists. Because such a run also commits the genes
    that succeeded, the baseline it was proposed against is gone and that
    frozen proposal is correctly refused afterwards; recovery is a new batch
    containing only the failed genes, composed against the committed baseline.
    For the same reason `continue` refuses to start from a reuse prefix that
    ends inside a gene rather than on a gene boundary, since no in-memory
    boundary exists for that gene
  - `primers inspect-gene-isoform-study-reuse BATCH_JSON_OR_@FILE
    --checkpoint-dir DIR [--path OUTPUT.json]` is read-only. It selects only the
    longest checkpoint whose baseline state, GENtle package/build/executable
    hash, and complete ordered operation prefix match the target batch, and
    emits `gentle.gene_isoform_assay_study_reuse_proposal.v1`. Appending a gene
    can therefore reuse an unchanged prefix; changing GENtle or any prefix
    payload cannot
  - reuse requires both the stored proposal and its exact
    `proposal_sha256`. Execution revalidates the proposal, manifest and engine
    bytes before importing the checkpoint into a detached run. Approval means
    “transfer this exact prior computation”; it does not validate assay biology
    and cannot waive an identity mismatch. Without those two flags, GENtle
    starts fresh even when compatible checkpoints exist
  - prior plans and unvalidated observations support an explicit next
    iteration. Retained assay ids are recorded but do not silently alter the
    automatic evidence recommendation
- Readiness-bound order bridge:
  - `primers oligo-order from-experimental-handoff HANDOFF.json
    --expected-sha256 SHA256` accepts only
    `gentle.experimental_assay_handoff.v1` rows whose matching card is
    `order_ready`, blocker-free, and bound to the embedded named readiness
    policy
  - every generated order line retains the handoff digest, readiness-policy id
    and digest, complete readiness row, card id, and assay/pair/oligo ids. It
    does not accept an unbound broad readiness label and never submits an order.
    Publication and order-sheet export recheck those identities and hashes
    against the included handoff instead of trusting copied status text
- Gene transcript-assay routine composition:
  - `primers compose-gene-assay-routine REQUEST_JSON_OR_@FILE [--path OUTPUT.json]`
    produces `gentle.gene_transcript_assay_routine.v1`
  - the request names an exported `gentle.gene_isoform_evidence.v1|v2` path,
    optional expected SHA-256, and persisted transcript-assay panel report ids.
    The input may be either a bare evidence report or the tagged
    `FeatureExpertView::IsoformEvidence` JSON emitted by the shared expert
    inspection path; the digest always covers the original input bytes
  - composition is a pure read: it does not rerun primer design, cDNA product
    tests, BLAST, or e-PCR and does not mutate the project
  - the report carries the isoform-evidence digest, current canonical panel
    digests, existing specificity acceptance status, uncovered transcript
    classes, endpoint/junction/common-control roles, order-table rows, and a
    conservative experimental sequence. Missing specificity is represented by
    absent status plus `specificity_accepted = false`, never as pass
  - a supplied expected isoform-evidence digest is checked before composition;
    specificity acceptance remains bound to the existing panel digest and is
    rejected as stale when those values differ
  - antibody/epitope compatibility and scalar evidence weighting are not
    inferred in v1. Additive antibody evidence requires a separately
    provenance-backed contract rather than guessing from protein mass or probe
    intensity
- Experimental assay handoff schemas:
  - `gentle.experimental_assay_handoff.v1` is a deterministic, read-only
    per-panel package. `BuildExperimentalAssayHandoff` consumes one persisted
    `gentle.transcript_assay_panel.v2`, runs the shared cDNA assay test once for
    every selected pair, and emits one `gentle.experimental_assay_card.v1` per
    pair plus a compact order/readiness table. It does not redesign primers,
    submit an order, or record a wet-lab result.
  - additive `coverage_summary` separates the broader annotation denominator
    (transcript records and exact distinct mature-cDNA sequences) from the
    UniProt-supported denominator (protein targets, linked transcript records,
    and distinct linked mature-cDNA sequences). It carries uncovered ids,
    annotation release, source report ids/digests, and plain summary lines.
    Legacy handoffs deserialize without this field; renderers must show the
    denominator as unavailable rather than displaying zero or inferring pass.
    Patch-locus records remain separate transcript records and may share one
    exact mature-cDNA sequence denominator.
    `assessed_*` and `unassessed_*` fields keep annotation members outside a
    narrowed panel distinct from members that were assessed and found
    uncovered. Optional group-level totals on `coverage_resolution` preserve
    the complete pre-narrowing exact-cDNA denominator; absence in a legacy
    report means unavailable, not zero.
    Linked-record numerators are evaluated independently from protein-target
    coverage: a record is covered only when its own digest-linked equivalence
    group is covered, and a distinct linked cDNA is covered only when the group
    assigned to that digest is covered. Legacy panel payloads without the
    record-level join report linked coverage as unavailable instead of
    promoting every record through a covered protein target.
  - additive `qa_aggregate_summary` derives typed totals once from the
    handoff's assay cards and readiness rows. It reports assessed/pass/fail/
    incomplete/not-evaluated counts separately for transcriptome specificity,
    genomic carryover, and critical oligo QC, plus order-ready totals and exact
    blocker/assay ids. Presentation adapters render these engine-owned totals
    rather than recomputing readiness.
  - canonical oligo identity removes all whitespace, uppercases, and maps
    `U` to `T`; full SHA-256-derived `oligo_id` and ordered `pair_id` values are
    authoritative for joins. Existing `primer_id` and `assay_id` fields remain
    backward-compatible report/display identities. `tube_id` is a short human
    label and is not the machine join key.
  - each card embeds the exact policy schema/version and every gate outcome.
    The shipped v1 default requires critical oligo QC, annotation provenance,
    a `genomic_carryover` pass, and a separate `transcriptome_specificity`
    pass. RT-PCR and qPCR cards therefore cannot become order-ready from a
    genomic search alone. The automatically generated cDNA assay test and
    optional `gentle.primer_variant_evidence.v1` are surfaced but absence alone
    is non-blocking; an evaluated failure remains a blocker. A linked order
    form retains its existing duplicate-review gate.
  - readiness states are `candidate`, `specificity_checked`, and `order_ready`;
    `wet_lab_validated` is reserved for a separate human-authored observation,
    never inferred from sequence evidence. Endpoint-gel abundance is described
    only as ordinal/semi-quantitative.
  - exact product-sequence classes are grouped by the cDNA product digest.
    Gel resolvability is evaluated separately and only when a named
    `GelRunConditions` policy is supplied. Distinct sequence digests never
    constitute a claim of visible band separation.
  - optional `virtual_gel` uses
    `gentle.experimental_assay_virtual_gel.v1`. It provenance-binds the source
    package/panel digest, render options, selected ladders, gel conditions,
    product count, predicted-empty card ids, and the emitted SVG SHA-256. Every
    assay card retains one sample lane so primer-pair columns are comparable.
    The SVG carries visible `[gentle]` attribution and machine-readable
    provenance attributes. Default visualization conditions are explicitly not
    a gel-resolution pass.
  - optional order-form formulations preserve modifications, scale, and
    purification as procurement identities. Tm, GC, secondary-structure, and
    specificity facts remain calculations for the unmodified annealing
    sequence and are not silently reinterpreted for modified chemistry. The
    package records the linked form id and a digest of its serialized content;
    card identity likewise binds the complete linked variant-evidence report.
  - the optional variant-evidence input is provenance-bound by pair id,
    reference assembly, source/release, population, retrieval time, and content
    digest. GENtle does not infer population frequencies from primer sequence
    alone.
- `primers preflight` returns `gentle.primer3_preflight.v1` with the requested
  backend plus configured-executable token, default-fallback marker, effective
  executable, resolved path, working directory, and reachability/version/error
  diagnostics. Its additive `progress_supported` field is `true` when
  `--help` advertises `--progress`, `false` when readable help omits it, and
  `null` when the capability could not be determined.
- a bare configured name such as the default `primer3_core` selects among every
  `PATH` match rather than the first: `--progress` support wins first, because
  bounded search reporting and SIGUSR1 status depend on it; then the highest
  version, compared component by component so `2.10.0` outranks `2.6.1`; then
  the latest filesystem creation time, falling back to the modification time
  where no birth time is recorded. A configured value naming a path is
  authoritative and used unchanged. The additive `candidates[]` records each
  considered build's path, `progress_supported`, `version`, `created_unix_ms`
  and `selected` flag, and `selection_reason` states why the winner won.
  `executable` remains the requested name; the chosen build is `resolved_path`.
  The selection is cached per candidate set and file identity so every
  operation in one run binds the same executable, which checkpoint reuse
  depends on.
- `primers prepare-restriction-cloning` expects an operation payload whose root
  variant is `{"PrepareRestrictionCloningPcrHandoff": {...}}`.
- `primers seed-restriction-cloning-handoff` is non-mutating and returns a
  validated `PrepareRestrictionCloningPcrHandoff` payload plus the
  recommendation context used to select or validate the enzyme pair.
- `primers restriction-cloning-vector-suggestions` is non-mutating and returns
  the same MCS-first / unique-cutter suggestion ordering the GUI PCR Designer
  uses for the selected destination vector.
- restriction-cloning saved-report helpers mirror the existing primer/qPCR
  lifecycle:
  - list stored handoff summaries
  - inspect one persisted handoff report by id
  - export one persisted handoff report to JSON
- `primers seed-from-feature` and `primers seed-from-splicing` are
  non-mutating helper commands that resolve an ROI and emit seeded operation
  payloads for both pair-PCR and qPCR design.
- `primers seed-qpcr-from-feature` and `primers seed-qpcr-from-splicing` are
  qPCR-only non-mutating helper commands that resolve an ROI and emit one
  seeded `DesignQpcrAssays` payload plus the built-in qPCR protocol-cartoon id
  (`pcr.assay.qpcr`) for shell/CLI/ClawBio promotion.
- qPCR-only seed helpers default to probe-based/TaqMan-like constraints:
  - primer length 18-24 bp
  - primer Tm 63-73 C with pair delta <= 3 C
  - probe length 20-30 bp
  - probe Tm 71-80 C, intended to land roughly 5-10 C above the primer mean
  - amplicon range 80-200 bp
- `primers seed-qpcr-from-splicing` additionally supports transcript-aware
  seeding:
  - `--mode shared_gene|distinguish_transcript`
  - `--transcript-id ID` when distinguishing one transcript from competing
    isoforms in the same gene/group
  - `--specificity-evidence junction_only|unique_exon_or_chain|either_prefer_junction`
    when `--mode distinguish_transcript`
- those qPCR-only seed helpers now also emit a deterministic `rationale`
  payload so agent callers can explain why the ROI was chosen and reuse
  GENtle's expected default assay limits without reverse-engineering them from
  the operation body.
- the qPCR designer can now use the same seeded report geometry to drive an
  in-window preview summary, so the shared shell payload and the GUI stay on
  one deterministic qPCR setup story.
- `primers design` and `primers show-report` additionally include
  `simple_pcr_pairs`, a derived helper array that summarizes each accepted pair
  in simple-PCR terms:
  - amplicon coordinates/length
  - left/right distance from the core ROI
  - left/right overlap into the core ROI
  - `left_to_core_label` / `right_to_core_label`
  - `flanks_core_cleanly`
  - `tm_delta_c` and score
- `primers oligo-order` persists first-class order forms inside
  `ProjectState.metadata["primer_design_reports"].oligo_order_forms` without
  changing the surrounding store schema (`gentle.primer_design_reports.v1`).
  Individual forms carry schema `gentle.oligo_order_form.v1`.
  - `create` accepts a generic `OligoOrderFormCreateRequest` JSON payload with
    `form_id?`, `target_label?`, `source_note?`, and `line_items[]`.
  - `from-primer-report` expands each requested pair rank into forward and
    reverse line items with source report, rank, role, template, operation/run,
    and source-coordinate provenance when available.
  - `from-qpcr-report` expands each requested assay rank into forward, reverse,
    and probe line items by default; `--include-probe false` suppresses the
    probe line.
  - line ids are deterministic hashes of form id, source kind/report/rank/role,
    and the normalized procurement tuple. They are not UUIDs and do not use
    time or randomness.
  - creation never merges duplicate line items. `duplicate_groups[]` records
    exact procurement duplicates keyed by `sequence_5_to_3 + modifications +
    scale + purification`; `sequence_reuse_groups[]` records same-sequence
    reuse with different procurement tuples. When duplicate groups exist,
    `duplicate_review.status` starts as `review_required`.
  - `review-dedup` marks duplicate review complete with the currently supported
    action `keep_separate`; it does not remove or merge any line item.
  - `route` converts the reviewed or unreviewed form into a provider-neutral
    `ExternalServiceDeliveryRouteRequest` whose `source_target.kind` is
    `oligo_order_form`, carrying `form_id`, `line_items[]`, duplicate-review
    metadata, duplicate groups, and sequence-reuse groups, then returns the
    standard `gentle.external_service_delivery_route.v1` recommendation.
  - `quote` builds a concrete `gentle.external_service_request.v1` from the
    same order-form source target and delegates to the existing external-service
    quote-handoff path. It refuses forms whose duplicate review is still
    `review_required`; run `primers oligo-order review-dedup FORM_ID` first.
    When provider or service kind is omitted, GENtle uses the delivery-route
    recommendation. `--output-dir` writes the same JSON/CSV/Markdown bundle as
    `services project-quote`.
  - oligo-order quote bundles preserve neutral line-item fields such as
    `line_id`, `line_no`, `name`, `role`, `sequence_5_to_3`, `length_nt`,
    `modifications`, `scale`, `purification`, report/rank provenance, and
    duplicate/reuse group hints for manual vendor-template transfer.
- Response schemas:
  - `gentle.primer_seed_request.v1`
  - `gentle.qpcr_seed_request.v1`
  - `gentle.primer_design_report.v1`
  - `gentle.primer_design_report_list.v1`
  - `gentle.terminal_exon_rt_primer_pool_request.v1`
  - `gentle.terminal_exon_rt_primer_pool.v1`
  - `gentle.qpcr_design_report.v1`
  - `gentle.qpcr_design_report_list.v1`
  - `gentle.oligo_order_form.v1`
  - `gentle.oligo_order_form_list.v1`
  - `gentle.oligo_order_form_create_result.v1`
  - `gentle.oligo_order_form_export.v1`
  - `gentle.oligo_order_duplicate_review.v1`
  - `gentle.restriction_cloning_pcr_handoff.v1`
  - `gentle.restriction_cloning_pcr_handoff_seed.v1`
  - `gentle.restriction_cloning_vector_enzyme_suggestions.v1`
- `gentle.oligo_order_form.v1` payload fields:
  - `form_id`, `target_label`, optional `source_note`
  - `created_at_unix_ms`, `updated_at_unix_ms`
  - `line_items[]`
    - `line_id`, `line_no`, `name`, `role`, `sequence_5_to_3`, `length_nt`
    - `modifications[]`, `scale`, `purification`, optional `notes`
    - `provenance`
      - `source_kind`, `report_id`, `report_schema`, `template`
      - optional `op_id`, `run_id`, `pair_rank`, `assay_rank`
      - `role`, `source_coordinates_0based[]`
  - `duplicate_groups[]`
    - `group_id`, `line_ids[]`, `sequence_5_to_3`, `modifications[]`,
      `scale`, `purification`
  - `sequence_reuse_groups[]`
    - `group_id`, `line_ids[]`, `sequence_5_to_3`,
      `procurement_tuple_count`
  - `duplicate_review`
    - `status` (`not_required`, `review_required`, or `reviewed`)
    - `default_action` (`keep_separate` in Phase A)
    - optional `reviewer`, `reviewed_at_unix_ms`, `note`
- `gentle.primer_seed_request.v1` payload fields:
  - `template`
  - `source` (`kind=feature|splicing`, `feature_id`, and splicing metadata when available)
  - `roi_start_0based`
  - `roi_end_0based_exclusive`
  - `operations.design_primer_pairs` (`{"DesignPrimerPairs": ...}`)
  - `operations.design_qpcr_assays` (`{"DesignQpcrAssays": ...}`)
- `gentle.qpcr_seed_request.v1` payload fields:
  - `template`
  - `source` (`kind=feature|splicing`, `feature_id`, and splicing metadata when available)
  - `roi_start_0based`
  - `roi_end_0based_exclusive`
  - `rationale`
    - `summary`
    - `why_this_roi`
    - `recommended_defaults`
      - `min_amplicon_bp`
      - `max_amplicon_bp`
      - `max_tm_delta_c`
      - `max_probe_tm_delta_c`
      - `primer_tm_c`
      - `probe_tm_c`
      - `probe_tm_offset_target_c`
      - `max_assays`
  - `operation` (`{"DesignQpcrAssays": ...}`)
    - splicing-seeded qPCR requests may include
      `transcript_targeting.source_feature_id`,
      `transcript_targeting.mode`, optional
      `transcript_targeting.transcript_id`, and optional
      `transcript_targeting.specificity_evidence`
  - `protocol_cartoon`
    - `protocol` (`pcr.assay.qpcr`)
    - `summary`
    - `default_output_svg`
- `gentle.restriction_cloning_vector_enzyme_suggestions.v1` payload fields:
  - `suggestions.seq_id`
  - `suggestions.selected_mcs[]` (preferred MCS-annotated cutters that are
    currently unique/usable)
  - `suggestions.other_unique[]` (other unique cutters on the vector)
  - `suggestions.missing_mcs[]` (annotated MCS cutters that were named but are
    not currently uniquely usable on the vector)
  - `suggestions.recommended_single_site[]`
    (`enzyme`, `cut_position_0based`)
  - `suggestions.recommended_directed_pairs[]`
    (`order_source`, `forward_enzyme`, `reverse_enzyme`,
    `forward_cut_position_0based`, `reverse_cut_position_0based`)
- `gentle.restriction_cloning_pcr_handoff_seed.v1` payload fields:
  - `primer_report_id`
  - `template`
  - `destination_vector_seq_id`
  - `pair_index`, `pair_rank`
  - `selected_pair`
  - `selected_pair_core_geometry`
  - `mode`
  - `forward_enzyme`, `reverse_enzyme`
  - `forward_leader_5prime`, `reverse_leader_5prime`
  - `selection_source`
  - optional `suggestion_order_source`
  - `vector_suggestions`
  - `operation` (`{"PrepareRestrictionCloningPcrHandoff": ...}`)

Feature formula shell contract (implemented):

- Shared-shell command:
  - `features formula SEQ_ID EXPR`
  - `EXPR` may be a single coordinate term (`=gene.tss`) or a range
    expression (`=gene.upstream(1000) .. gene.tss`)
- Execution semantics:
  - non-mutating engine inspection over one sequence
  - uses the same engine-owned feature-coordinate resolver as GUI selection
    and primer/qPCR ROI fields
  - supports occurrence selectors, label selectors, signed offsets, `tss`, and
    strand-aware `upstream(N)`
  - range expressions normalize endpoint order after resolution, so
    reverse-strand upstream/TSS formulas yield a local
    `start_0based..end_0based_exclusive` interval
- Response schema:
  - `gentle.feature_coordinate_formula_resolution.v1`
  - coordinate results include `coordinate_0based`
  - range results include `range.start_0based`,
    `range.end_0based_exclusive`, and `range.length_bp`

Feature-location edit contract (implemented):

- Schemas:
  - `gentle.feature_location_edit.v1` remains the exact simple-location
    request/report shape.
  - `gentle.feature_location_edit.v2` is emitted when `segment_index` selects
    one child in a supported flat compound. The same Rust contract carries
    optional v2 fields; simple serialization omits them and retains the v1 key
    set.
- Shared operations:
  - `PreviewFeatureLocationEdit` validates an edit and emits the complete
    before/after report without changing project state or creating an undo
    checkpoint.
  - `EditFeatureLocation` applies the same request only when
    `expected_feature_fingerprint_sha256` matches the complete current feature
    record (kind, exact location structure, and ordered qualifiers).
  - fingerprints use
    `sha256_gb_io_feature_serde_json_ordered_qualifiers_v1` and the repository
    `sha256:`-prefixed digest representation.
- Shared-shell command:
  - `features edit-location SEQ_ID FEATURE_INDEX [--segment-index INDEX] --start-1based N --end-1based-inclusive M [--dry-run] [--expected-feature-fingerprint-sha256 SHA] [--path OUT.json]`
  - run `--dry-run` first, then pass its
    `before_feature_fingerprint_sha256` when applying.
- Scope is intentionally strict:
  - v1 edits exact `Range` and `Complement(Range)` locations.
  - v2 edits one existing exact child of flat `Join(Range...)`,
    `Order(Range...)`, `Complement(Join(...))`, or
    `Complement(Order(...))`.
  - v2 clones the location AST and replaces only the selected child; container
    kind, complement wrapper, child count/order, qualifiers, and every
    unedited node are preserved.
  - nested, fuzzy, between-base, external, bond, one-of, gap, non-monotonic,
    and circular cross-origin compounds are rejected rather than normalized.
- Strand is preserved. Reports expose local half-open coordinates, 1-based
  inclusive coordinates, and strand-aware 5-prime/3-prime positions.
  `related_features[]` lists annotations sharing the old start or end boundary
  for human review; those annotations are never changed automatically.
  Segment reports additionally expose `target_scope=segment`,
  `compound_context` (operator, stored direction/index, biological segment
  number, and count), typed `compound_validation_warnings[]`, and
  `related_segment_boundaries[]` with both interval-boundary roles.
- Invalid, empty, negative, inverted, out-of-range, and unsupported edits are
  errors. Segment overlap, departure from the compound's established stored
  direction, and CDS coding-length deltas not divisible by three are
  non-blocking review warnings; GENtle does not rewrite `/codon_start`.
- The fingerprint is an optimistic lock on the complete current feature, not a
  signature of the requested coordinates. GUI changes to feature, segment, or
  coordinates invalidate the cached preview before Apply is enabled.

Feature-record curation contract (implemented):

- Schema:
  - `gentle.feature_record_curation.v2`
  - v1 is the Create/Delete-only predecessor; v2 adds strict Split/Merge
    request and outcome variants without changing v1 field meanings
- Shared operations:
  - `PreviewFeatureRecordCuration`
  - `ApplyFeatureRecordCuration`
  - both accept a tagged `FeatureRecordCurationRequest` with
    `operation_kind=create|delete|split|merge`.
- Shared-shell commands:
  - `features create SEQ_ID --kind KIND --start-1based N --end-1based-inclusive M [--strand forward|reverse] [--qualifier KEY[=VALUE] ...] [--dry-run] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
  - `features delete SEQ_ID FEATURE_INDEX [--dry-run] [--expected-feature-fingerprint-sha256 SHA] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
  - `features split SEQ_ID FEATURE_INDEX --split-before-1based N [--dry-run] [--expected-feature-fingerprint-sha256 SHA] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
  - `features merge SEQ_ID FIRST_FEATURE_INDEX SECOND_FEATURE_INDEX [--dry-run] [--expected-first-feature-fingerprint-sha256 SHA] [--expected-second-feature-fingerprint-sha256 SHA] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
  - run `--dry-run` first. Create apply requires the returned annotation-state
    fingerprint; Delete apply requires both that fingerprint and the deleted
    feature fingerprint; Split has the same two-lock rule; Merge requires the
    annotation-state fingerprint plus both source-feature fingerprints.
- Create appends one exact `Range` or `Complement(Range)` to the ordered feature
  table. Preview deliberately does not promise the eventual feature index.
  Qualifiers are an ordered list of `{key, value}` records: duplicate keys,
  interleaving, empty string values, and valueless qualifiers remain distinct.
- Delete removes one complete existing feature record and accepts any location
  shape already represented by `gb-io`, including complex/fuzzy locations.
  The report carries a lossless serialized location, human display location,
  exact ordered qualifiers, and the number of later feature indices that shift
  down by one.
- Split accepts only exact simple `Range` or `Complement(Range)` locations and
  an internal boundary. It replaces table index `i` with genomic-left and
  genomic-right records at `[i, i+1]`, retaining the original kind, strand, and
  exact ordered qualifier vector on both outputs. The report names both output
  indices and the number of later indices shifted up.
- Merge accepts only two distinct, exactly touching simple records with equal
  kind, strand, and exact ordered qualifier vectors. It spans their combined
  genomic interval, retains the lower source table index, removes the higher
  index, and reports the resulting index shift. Gaps, overlaps, compound or
  fuzzy locations, strand differences, and qualifier conflicts are rejected;
  GENtle does not infer a metadata reconciliation.
- Annotation-state fingerprints use
  `sha256_sequence_id_length_topology_ordered_gb_io_features_serde_json_v1`.
  They cover sequence id, sequence length, topology, and every complete feature
  in stored order, and remain stable across project save/load. They are
  optimistic-concurrency locks, not biological signatures.
- `review_candidates[]` reports geometric overlap and exact shared recognized
  INSDC identifiers (`locus_tag`, `gene`, `protein_id`, `transcript_id`) as
  informational evidence only. It never interprets overlap as an error or a
  shared identifier as a dependency, and it never edits another annotation.
- All applies are ordinary full-checkpoint mutations with undo/redo.
  Nested-location creation, dependency propagation, qualifier reconciliation,
  and automatic transcript repair remain outside v2.

Feature-query shell contract (implemented):

- Shared-shell command:
  - `features query SEQ_ID [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--nearest-to POSITION] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
- Execution semantics:
  - non-mutating engine inspection over one sequence’s feature table
  - deterministic filter pipeline:
    kind include/exclude, optional range relation (`overlap|within|contains`),
    strand filter, label contains/regex, qualifier filters, and length bounds
  - deterministic ordering by requested sort key with stable tie-breaks +
    pagination (`offset`/`limit`)
  - optional `nearest_to_0based` / `--nearest-to POSITION` replaces the primary
    ordering with nearest-base interval distance; equal distances are ordered
    by start, end, and feature id. The point must lie within the local sequence.
- Response schema:
  - `gentle.sequence_feature_query_result.v1`
  - fields include:
    - `seq_id`, `sequence_length_bp`, `total_feature_count`
    - `matched_count`, `returned_count`, `offset`, `limit`, optional
      `nearest_to_0based`
    - normalized `query`
    - `rows[]` with `feature_id`, `kind`, `start_0based`,
      `end_0based_exclusive`, `length_bp`, optional `distance_to_query_bp`,
      `strand`, `label`, `labels[]`, and optional qualifier maps when requested
      (`--include-qualifiers`)

Feature BED export contract (implemented):

- Shared-shell command:
  - `features export-bed SEQ_ID OUTPUT.bed [--coordinate-mode auto|local|genomic] [--include-restriction-sites] [--restriction-enzyme NAME] [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--nearest-to POSITION] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
- Raw/shared operation:
  - `{"ExportFeaturesBed":{"query":{"seq_id":"tp53_region","kind_in":["gene","mRNA","exon","CDS"]},"path":"artifacts/tp53_region.features.bed","coordinate_mode":"auto","include_restriction_sites":false,"restriction_enzymes":[]}}`
- Execution semantics:
  - non-mutating export built on the same feature-query filter contract used by
    `features query`
  - when `--limit` / `query.limit` is omitted, the exporter writes all matching
    rows instead of defaulting to the paged query window
  - `coordinate_mode=auto` prefers genomic BED coordinates whenever a feature
    carries `chromosome`, `genomic_start_1based`, and `genomic_end_1based`;
    otherwise the row falls back to local `SEQ_ID` coordinates
  - `include_restriction_sites=true` appends deterministic REBASE-derived
    `restriction_site` rows, filtered by the same range/strand/label/qualifier
    options; `restriction_enzymes[]` narrows those rows to selected enzymes
  - TFBS/JASPAR annotations remain ordinary feature rows, so `kind_in=["TFBS"]`
    exports the current binding-site table after `AnnotateTfbs`

Restriction-site scan contract (implemented):

- Shared-shell commands:
  - `features restriction-scan SEQ_ID [--range START..END|--start N --end N] [--enzyme NAME] [--max-sites-per-enzyme N] [--no-cut-geometry] [--path FILE.json]`
  - `features restriction-scan --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] [--range START..END|--start N --end N] [--enzyme NAME] [--max-sites-per-enzyme N] [--no-cut-geometry] [--path FILE.json]`
- Raw/shared operation:
  - `{"FindRestrictionSites":{"target":{"kind":"seq_id","seq_id":"tp53_region","span_start_0based":700,"span_end_0based_exclusive":1200},"enzymes":["EcoRI","SmaI"],"include_cut_geometry":true}}`
- Execution semantics:
  - non-mutating scan backed by the same shared REBASE enzyme catalog used by
    other restriction-aware flows
  - accepts either stored `seq_id` targets or inline ASCII DNA text through one
    contract
  - when `enzymes[]` is omitted, the shared preferred restriction-enzyme list
    is used, with fallback to the default preferred set
  - the report includes both local scan coordinates and source-sequence
    coordinates so selection/visible-span calls stay easy to interpret
- TFBS/JASPAR scan contract (implemented):
  - Shared-shell commands:
    - `features tfbs-scan SEQ_ID --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--min-llr-bits VALUE] [--min-llr-quantile VALUE] [--per-tf-min-llr-bits TF=VALUE] [--per-tf-min-llr-quantile TF=VALUE] [--max-hits N] [--path FILE.json]`
    - `features tfbs-scan --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--min-llr-bits VALUE] [--min-llr-quantile VALUE] [--per-tf-min-llr-bits TF=VALUE] [--per-tf-min-llr-quantile TF=VALUE] [--max-hits N] [--path FILE.json]`
  - Raw/shared operation:
    - `{"ScanTfbsHits":{"target":{"kind":"seq_id","seq_id":"tp53_region","span_start_0based":700,"span_end_0based_exclusive":1200},"motifs":["SP1","TP73"],"min_llr_quantile":0.95,"max_hits":50}}`
  - Execution semantics:
    - non-mutating PWM scan backed by the same shared local motif registry and
      IUPAC scoring helpers used by `AnnotateTfbs`,
      `SummarizeTfbsScoreTracks`, and JASPAR presentation/expert routes
    - accepts either stored `seq_id` targets or inline ASCII DNA text through
      one contract
    - returns both local scan coordinates and source-sequence coordinates so
      selection/visible-span calls remain interpretable without materializing a
      new sequence record
    - full-span circular targets also mark origin-crossing rows with
      `wraps_origin=true`; those rows report modulo-local end coordinates while
      keeping the full wrapped `matched_sequence`
    - `max_hits` caps the total retained row count across the whole scan in
      deterministic motif order and sets `truncated_at_max_hits=true` when
      later hits are skipped

Repeat-environment cohort contract (implemented baseline):

- Shared-shell commands:
  - `features repeat-query GENOME_ID --rmsk PATH [--rep-class CLASS] [--rep-family FAMILY] [--rep-name NAME] [--alias ALIAS] [--chromosome CHR] [--range START..END] [--limit N] [--path FILE.json]`
  - `features repeat-overlaps SEQ_ID --index RMSK_INTERVAL_INDEX.json [--range START..END] [--limit N] [--path FILE.json]`
  - `features materialize-repeats SEQ_ID --index RMSK_INTERVAL_INDEX.json [--max-features N] [--append] [--path FILE.json]`
  - `features repeat-cohort GENOME_ID --rmsk PATH [--rep-class CLASS] [--rep-family FAMILY] [--rep-name NAME] [--alias ALIAS] [--geometry repeat_midpoint|transcript_5utr_start|pol2_promoter_upstream|cds_stop_context] [--upstream-bp N] [--downstream-bp N] [--limit N] [--catalog PATH] [--cache DIR] [--path FILE.json]`
  - `features window-cohort-tfbs COHORT_JSON --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--cache DIR] [--path FILE.json]`
  - `genomes extract-region|extract-gene|extract-promoter ... [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features]`
- Raw/shared operations:
  - `{"QueryRepeatAnnotations":{"genome_id":"Human GRCh38 Ensembl 116","rmsk_path":"data/rmsk.txt.gz","filter":{"normalized_alias":"LINE/L1"},"limit":100}}`
  - `{"QueryRepeatOverlaps":{"seq_id":"grch38_tp53","rmsk_index_path":"data/resources/ucsc.rmsk.hg38.interval-index.json","start_0based":0,"end_0based_exclusive":5000,"limit":100}}`
  - `{"MaterializeRepeatFeatures":{"seq_id":"grch38_tp53","rmsk_index_path":"data/resources/ucsc.rmsk.hg38.interval-index.json","max_features":1000,"clear_existing":true}}`
  - `{"BuildRepeatEnvironmentCohort":{"genome_id":"Human GRCh38 Ensembl 116","rmsk_path":"data/rmsk.txt.gz","filter":{"normalized_alias":"LINE/L1"},"geometry_mode":"pol2_promoter_upstream","upstream_bp":2000,"downstream_bp":2000,"limit":50}}`
  - `{"SummarizeWindowCohortTfbs":{"cohort":{...},"motifs":["SP1","TP73"],"score_kind":"llr_background_tail_log10","clip_negative":true}}`
- Execution semantics:
  - `QueryRepeatAnnotations` parses UCSC `rmsk.txt` / `rmsk.txt.gz` rows and
    preserves deterministic repeat identity fields: `repName`, `repClass`,
    `repFamily`, strand, chromosome, 1-based genomic span, and optional
    score/divergence fields.
  - Repeat filters accept raw class/family/name values plus normalized aliases
    such as `LINE/L1`, `SINE/Alu`, and `LTR/ERV`; malformed rows are counted
    and skipped instead of aborting the whole deterministic query.
  - `QueryRepeatOverlaps` requires a sequence with genome-extraction anchor
    provenance. It projects prepared-index intervals into local sequence
    coordinates, clips optional local query ranges, flips local strand on
    reverse anchors, resolves common chromosome aliases, and returns
    `gentle.sequence_repeat_overlap.v1`.
  - `gentle.sequence_repeat_overlap.v1` includes both row-level overlaps and
    neutral metrics: `query_length_bp`, `total_overlap_bp`,
    `covered_query_bp`, `coverage_fraction`, `nearest_repeat_distance_bp`, and
    class/family/alias summary rows.
  - `MaterializeRepeatFeatures` uses the same projection service and writes
    ordinary `repeat_region` features with `gentle_generated=ucsc_rmsk`,
    `rmsk_*`, `repeat_*`, and original `repName`/`repClass`/`repFamily`
    qualifiers plus score/divergence fields. By default it clears prior
    generated UCSC-rmsk features before
    writing; `--append` keeps existing generated rows and skips duplicates by
    stable annotation id. The report schema is
    `gentle.repeat_feature_materialization.v1`.
  - Reference extraction commands with `--rmsk-index` apply
    `MaterializeRepeatFeatures` immediately to the first created anchored
    sequence and return a sibling `repeat_materialization` result object.
  - `BuildRepeatEnvironmentCohort` creates one row per selected repeat locus,
    then projects transcript/gene context when a prepared genome catalog is
    available. The report stores all geometry windows, the selected geometry,
    and explicit missing-reason strings for unavailable 5'UTR/TSS/CDS-stop
    contexts.
  - Supported geometry modes are `repeat_midpoint`,
    `transcript_5utr_start`, `pol2_promoter_upstream`, and
    `cds_stop_context`. The default flank on each side is 2000 bp.
  - `SummarizeWindowCohortTfbs` consumes the stored cohort window records and
    scores the selected windows with the same TFBS score-track machinery used
    by promoter and inline-DNA analyses; minus-strand windows are
    reverse-complemented before scoring so cohort rows stay strand-oriented.
  - RNA-read evidence fields are part of the cohort row shape but remain
    annotation/ranking extension points in this baseline; repeat selection is
    never hard-filtered by RNA evidence.
- Response/report schemas:
  - `gentle.repeat_annotation_query.v1`
  - `gentle.repeat_environment_cohort.v1`
  - `gentle.window_cohort_tfbs.v1`
- File format:
  - BED6 core columns:
    `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`
  - deterministic extra columns:
    `kind`, `row_id`, `coordinate_source`, `qualifiers_json`
- Response/report schema:
  - `gentle.sequence_feature_bed_export.v1`
  - fields include:
    - `seq_id`, `path`, `coordinate_mode`
    - `matched_sequence_feature_count`, `matched_restriction_site_count`,
      `matched_row_count`
    - `exportable_row_count`, `exported_row_count`
    - `local_coordinate_row_count`, `genomic_coordinate_row_count`
    - `skipped_missing_genomic_coordinates`
    - `bed_columns[]`

Sequence-context view summary contract (implemented baseline):

- Raw/shared operation:
  - `{"InspectSequenceContextView":{"seq_id":"rs9923231_16_31093368_31099368","mode":"linear","viewport_start_0based":2400,"viewport_end_0based_exclusive":3501,"include_visible_classes":["gene","mrna","variation","tfbs"],"coordinate_mode":"genomic","limit":20}}`
- Execution semantics:
  - non-mutating compact summary for one DNA-sequence-view context
  - reuses the shared feature-query pipeline instead of inventing a second
    frontend-only inspection path
  - defaults to the current stored linear viewport when no explicit viewport is
    supplied and the engine display state already carries one
  - defaults to the currently enabled visible feature classes when
    `include_visible_classes[]` is omitted
  - surfaces a short chat-friendly/text-friendly summary alongside a compact
    feature table so callers like ClawBio can answer first and attach larger
    SVG/BED artifacts second
- Response/report schema:
  - `gentle.sequence_context_view.v1`
  - fields include:
    - `seq_id`, `sequence_length_bp`
    - `mode`, `coordinate_mode`
    - `viewport_start_0based`, `viewport_end_0based_exclusive`,
      `viewport_span_bp`
    - optional `genome_anchor` plus optional genomic viewport bounds when the
      sequence is genome-anchored
    - `visible_classes[]` with class ids, feature-kind coverage, counts, and
      prominent labels
    - `matched_feature_count`, `returned_feature_count`, `limit`
    - `rows[]` with local bounds, label metadata, and optional genomic feature
      coordinates when qualifiers or anchor geometry make them available
    - `summary_lines[]` for compact relay into chat/report layers

Sequence-context bundle export contract (implemented baseline):

- Raw/shared operation:
  - `{"ExportSequenceContextBundle":{"seq_id":"rs9923231_16_31093368_31099368","mode":"linear","viewport_start_0based":2400,"viewport_end_0based_exclusive":3501,"coordinate_mode":"genomic","include_feature_bed":true,"include_text_summary":true,"include_restriction_sites":false,"restriction_enzymes":[],"output_dir":"artifacts/rs9923231.context_bundle"}}`
- Execution semantics:
  - non-mutating artifact export for one DNA-sequence-view context
  - writes one deterministic output directory containing:
    - `context.svg`
    - `context_summary.json`
    - optional `context_summary.txt`
    - optional `context_features.bed`
    - `bundle.json` manifest
  - reuses the same shared `InspectSequenceContextView` report for the compact
    summary component instead of inventing a second export-only summary format
  - reuses the existing `ExportFeaturesBed` path for the BED companion when
    `include_feature_bed=true`
  - defaults to current display-visible classes for the BED feature query, so
    the coordinate-bearing table tracks the same viewer focus as the summary
- Response/report schema:
  - `gentle.sequence_context_bundle.v1`
  - fields include:
    - `seq_id`, `output_dir`
    - `svg_path`, `summary_json_path`, optional `summary_text_path`,
      optional `feature_bed_path`, `bundle_json_path`
    - `include_text_summary`, `include_feature_bed`,
      `include_restriction_sites`, `restriction_enzymes[]`
    - `artifacts[]` with deterministic presentation metadata:
      `artifact_id`, `path`, `media_type`, `artifact_kind`, `caption`,
      `recommended_use`, `presentation_rank`, `is_best_first_artifact`
    - `best_first_artifact_id`, `best_first_artifact_path`
    - embedded `sequence_context_view`
    - optional embedded `feature_bed_export`

Dotplot + flexibility operation contract (implemented baseline):

- Dotplot operation:
  - `ComputeDotplot { seq_id, reference_seq_id?, span_start_0based?, span_end_0based?, reference_span_start_0based?, reference_span_end_0based?, mode, word_size, step_bp, max_mismatches?, tile_bp?, store_as?, inspection_provenance? }`
  - `ComputeDotplotOverlay { owner_seq_id, reference_seq_id, reference_span_start_0based?, reference_span_end_0based?, queries[], word_size, step_bp, max_mismatches?, tile_bp?, store_as? }`
  - `mode`: `self_forward | self_reverse_complement | pair_forward | pair_reverse_complement`
  - pair modes require `reference_seq_id` and use the optional
    `reference_span_start_0based` / `reference_span_end_0based` for the
    y/reference axis.
  - `ComputeDotplotOverlay` is reference-centered and requires at least one
    query spec; each query uses `pair_forward` or `pair_reverse_complement`
    against the same reference span
  - stores payload schema `gentle.dotplot_view.v3`
  - payload includes:
    - `owner_seq_id`
    - optional `op_id` / `run_id` identifying the operation that created the
      stored payload
    - optional `inspection_provenance`
      (`gentle.dotplot_inspection_provenance_citation.v1`) when
      construct reasoning requested the dotplot; it binds graph/action
      identity and input fingerprint, source facts/annotations/candidates/
      summaries/evidence, rationale, and the resolved dotplot request
    - citation status is `pass`, `fail`, or `unknown`: request or stored-source
      mismatches and stale graph snapshots fail verification, while legacy
      graph snapshots whose freshness cannot be established remain `unknown`
    - absence of `inspection_provenance` means manual or otherwise unknown
      provenance and is distinct from a present citation with `status=fail`
    - shared reference span + seed parameters
    - sparse match points (`points[]`)
    - per-query-bin reference-distribution boxplot summary
      (`boxplot_bin_count`, `boxplot_bins[]` with
      `min/q1/median/q3/max + hit_count`)
    - `query_series[]` for multi-query overlays
      - each series may also carry optional `query_anchor_0based` +
        `query_anchor_label` values for curated cross-family/domain-anchored
        rendering with `query_anchor_bp`
    - optional `reference_annotation` genome-context side rail for pairwise,
      annotated self, and overlay dotplot payloads
      - intervals carry `kind`, `label`, optional `strand`, `lane`,
        optional `color_rgb`, and optional human-readable `detail`
      - exon intervals are not filtered to one selected gene; sense and
        antisense exons remain separate so opposite-strand context survives
      - materialized RepeatMasker/UCSC `rmsk` repeat annotations
        (`repeat_region` / `mobile_element` features with repeat qualifiers)
        are represented as repeat intervals
      - dotplot computation/rendering does not perform live repeat-index
        queries; repeat context must already be present as sequence features
    - optional `query_series[].transcript_feature_id` when overlays originate
      from locus transcript lanes
    - optional `overlay_anchor_exons[]` carrying precomputed shared-exon anchor
      metadata for `shared_exon_anchor` rendering/export
  - guardrails:
    - `word_size >= 1`
    - `step_bp >= 1`
    - query/reference spans must satisfy `0 <= start < end <= sequence_len`
    - pair-evaluation safety limit is enforced for latency protection
    - point count is capped with deterministic truncation warning
- Flexibility operation:
  - `ComputeFlexibilityTrack { seq_id, span_start_0based?, span_end_0based?, model, bin_bp, smoothing_bp?, store_as? }`
  - `model`: `at_richness | at_skew`
  - stores payload schema `gentle.flexibility_track.v1`
  - guardrails:
    - `bin_bp >= 1`
    - same span validation contract as dotplot
    - optional smoothing uses deterministic moving-average bins
- Metadata persistence:
  - metadata key: `dotplot_analysis`
  - store schema: `gentle.dotplot_analysis_store.v1`
  - both dotplots and flexibility tracks are persisted under this key
  - `DotplotViewSummary` carries the same optional operation/run identity and
    inspection citation as the full view for listing and audit without loading
    sparse match points
- Shared-shell command family:
  - `dotplot compute SEQ_ID [--reference-seq REF_SEQ_ID] [--start N] [--end N] [--ref-start N] [--ref-end N] [--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
  - `dotplot overlay-compute OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - convenience wrapper over `ComputeDotplotOverlay`
    - if `--reference-seq` is omitted, adapters default the shared reference to
      `OWNER_SEQ_ID`
    - each `--query-spec` must deserialize into one `DotplotOverlayQuerySpec`
  - `dotplot list [SEQ_ID]`
  - `dotplot show DOTPLOT_ID`
  - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]`
  - `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]` (alias)
  - `flex compute SEQ_ID [--start N] [--end N] [--model at_richness|at_skew] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]`
  - `flex list [SEQ_ID]`
  - `flex show TRACK_ID`

Splicing-reference derivation + pairwise alignment operation contract (implemented baseline):

- Splicing-reference derivation operation:
  - `DeriveSplicingReferences { seq_id, span_start_0based, span_end_0based, seed_feature_id?, scope?, output_prefix? }`
  - emits multiple derived sequence outputs from one genomic window:
    - DNA window (`..._dna`)
    - one mRNA sequence per transcript lane (`..._mrna_*`, transcript orientation, `T->U`)
    - exon-consecutive artificial reference sequence (`..._exon_reference`)
  - shared splicing expert payloads now also carry splice-boundary motif
    annotations per intron boundary:
    - donor/acceptor dinucleotide (`motif_2bp`)
    - paired intron signature (`paired_motif_signature`, for example `GT-AG`)
    - paired boundary position (`partner_position_1based`)
    - high-level motif class / annotation
      (`gt_ag_major_canonical`, `gc_ag_major_noncanonical`,
      `at_ac_minor_u12_like`, `at_ag_minor_u12_like`,
      `other_noncanonical`)
  - shared splicing expert payloads also carry conservative intron-signal
    heuristics per intron:
    - donor/acceptor positions and intron length
    - best branchpoint-like adenine in the usual 18-40 nt acceptor-proximal
      window (`branchpoint_position_1based`, `branchpoint_motif`,
      `branchpoint_score`, `branchpoint_annotation`)
    - best acceptor-proximal polypyrimidine-rich tract
      (`polypyrimidine_start_1based`, `polypyrimidine_end_1based`,
      `polypyrimidine_fraction`, `polypyrimidine_annotation`)
    - these are intentionally framed as heuristics, not a splice predictor
  - if `seed_feature_id` is omitted, engine selects one overlapping mRNA feature deterministically from the requested span
  - default `scope`: `target_group_target_strand`
- ATtRACT splice-aware evidence inspection:
  - normalized runtime motif snapshot schema:
    - `gentle.attract_motifs.v1`
    - records one deterministic ATtRACT motif row per normalized
      `ATtRACT_db.txt` entry, including:
      - `pwm_row_count`
      - `consensus_only_row_count`
      - optional `snapshot_fingerprint` on the snapshot wrapper so runtime and
        session-loaded resources can be identified by content, not just by item
        count
      - `entry_id`
      - `matrix_id`
      - `gene_name`
      - `organism`
      - `motif_iupac`
      - optional `pfm` rows (`a`, `c`, `g`, `t`) when a `pwm.txt` block could
        be mapped to the row’s `matrix_id`
      - optional provenance fields such as `experiment`, `family`, `domain`,
        `pubmed_id`, and `quality_score`
      - `model_kind` (`consensus_iupac` or `pwm_counts`)
      - `pwm_present` so downstream consumers can distinguish “PWM existed in
        the archive” from “this particular row still fell back to
        consensus/IUPAC”
  - inspection payload schema:
    - `gentle.attract_splicing_evidence.v1`
    - `AttractSplicingEvidenceSettings`
      - `scope`
      - `transcript_strand_only`
      - `boundary_flank_bp`
      - optional `requested_organism`
      - `allow_species_fallback`
      - `minimum_quality_score`
      - `minimum_match_quantile`
      - `pwm_mapping_policy`:
        - `strict_same_length` (default)
        - `windowed_submatrix`
      - `compare_alternate_policy`
    - `AttractSpeciesMatchMode`:
      - `exact_organism`
      - `fallback_all_compatible`
    - `AttractRegionClass`:
      - `exon_body`
      - `donor_flank`
      - `acceptor_flank`
      - `intron_body`
    - summary rows:
      - grouped by factor / organism / matrix id
      - include strongest match score, its score kind/quantile, motif-quality
        maximum, hit count, region-class counts, PWM provenance counts
        (`exact_length_pwm_hits`, `windowed_pwm_hits`, `consensus_only_hits`),
        and supporting transcript ids
    - hit rows:
      - transcript identity/strand
      - factor / organism / matrix id / motif
      - region class
      - genomic coordinates
      - local coordinates within the scanned exon/intron window
      - matched sequence
      - `match_score`
      - `match_score_kind`
      - optional `match_score_quantile`
      - motif `quality_score`
      - exact-species flag
      - PWM provenance:
        - `pwm_mapping_status`
        - `mapping_policy_used`
        - optional `pfm_subwindow_start_1based`
        - optional `pfm_subwindow_end_1based`
    - provenance:
      - active ATtRACT source label
      - active resource item count
      - active resource PWM-backed row count
      - active resource consensus-only row count
      - optional `active_resource_fingerprint`
      - optional `alternate_policy_summary`
      - requested/resolved organism
      - species-match mode
      - scan warnings
    - aggregate hit counters:
      - `pwm_scored_hit_count`
      - `exact_length_pwm_hit_count`
      - `windowed_pwm_hit_count`
      - `consensus_hit_count`
  - current v1 scan semantics:
    - transcript-strand aware by default
    - scans exon bodies and intron bodies from the selected splicing group
    - promotes hits to `donor_flank` / `acceptor_flank` when they fall within
      `boundary_flank_bp` of the corresponding splice boundary
    - exact-organism match is preferred; broader fallback is explicit and
      recorded in the payload
    - normalized consensus/IUPAC matching remains the conservative candidate
      gate
    - under `strict_same_length`, only rows whose motif length exactly matches
      the linked PWM width use PWM `llr_bits` ranking
    - under `windowed_submatrix`, rows whose linked PWM is longer than the
      motif may reuse one unique consensus-compatible PWM subwindow; those hits
      are labeled as `llr_bits_windowed` with the chosen subwindow range
    - rows without a mapped PWM block continue to use deterministic
      consensus/IUPAC exact matching only
    - interpretation guidance:

Splicing-expert composition note:

- the GUI now also combines `SplicingExpertView.intron_signals` with an
  already-cached `gentle.attract_splicing_evidence.v1` payload to produce one
  intron-centered regulatory interpretation table:
  - branchpoint-like / polypyrimidine heuristics stay visible as conservative
    splice-mechanistic context
  - cached RBP hits are summarized per intron as donor-flank, acceptor-flank,
    and intron-body support counts
  - no extra biology is recomputed in the GUI; the view is a deterministic
    composition of existing shared payloads
      - `strict_same_length` is the canonical conservative mode
      - `windowed_submatrix` is an engine-supported GENtle heuristic, not an
        ATtRACT-published mapping rule
      - `compare_alternate_policy` is intended to help reviewers quantify how
        much that heuristic changes the evidence before switching the active
        mapping mode
  - contract boundary:
    - this is an RBP/splicing interpretation payload, not the generic
      TFBS/PSSM score-track payload
    - future PWM/PSSM-backed ATtRACT scoring should reuse shared motif math
      under the engine without collapsing these higher-level payloads into one
      another
- Pairwise alignment operation:
  - `AlignSequences { query?, target?, query_seq_id?, target_seq_id?, query_span_start_0based?, query_span_end_0based?, target_span_start_0based?, target_span_end_0based?, mode?, match_score?, mismatch_score?, gap_open?, gap_extend? }`
  - preferred operands: `query` and `target` as `SequenceScanTarget` (`SeqId` or `InlineSequence`)
  - legacy operands: `query_seq_id`/`target_seq_id` plus optional legacy span fields remain accepted for existing workflows
  - `mode`: `global | local` (default `global`)
  - scoring defaults: `match=2`, `mismatch=-3`, `gap_open=-5`, `gap_extend=-1`
  - returns structured payload `sequence_alignment` with spans, score, coverage, identity, and CIGAR-like compact operations string
  - non-mutating operation (no sequence/container state mutation)
- Shared-shell command family:
  - `splicing-refs derive SEQ_ID START_0BASED END_0BASED [--seed-feature-id N] [--scope all_overlapping_any_strand|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
  - `align compute QUERY_SEQ_ID TARGET_SEQ_ID [--query-start N] [--query-end N] [--target-start N] [--target-end N] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N]`
  - `align compute --query-sequence-text DNA --target-sequence-text DNA [--query-id-hint TEXT] [--target-id-hint TEXT] [--query-range START..END] [--target-range START..END] [--mode global|local]`
  - `attract inspect-splicing SEQ_ID FEATURE_ID [--scope ...] [--organism NAME] [--flank-bp N] [--min-score X] [--min-match-quantile Q] [--pwm-mapping strict_same_length|windowed_submatrix] [--compare-policies] [--all-transcripts] [--no-fallback]`

RNA-read interpretation contract (Nanopore cDNA phase-1 baseline):

- Operations:
  - `InterpretRnaReads { seq_id, seed_feature_id, profile, input_path, input_format, scope, origin_mode?, target_gene_ids?, roi_seed_capture_enabled?, seed_filter, align_config, report_id?, report_mode?, checkpoint_path?, checkpoint_every_reads?, resume_from_checkpoint? }`
  - `AlignRnaReadReport { report_id, selection, align_config_override?, selected_record_indices? }`
  - `PreflightRnaReadIsoforms { seq_id, seed_feature_id, scope?, seed_filter?, optimize_parameters?, positive_transcript_fasta_paths?, control_transcript_fasta_paths?, max_control_match_probability? }`
  - `SummarizeRnaReadGeneSupport { report_id, gene_ids, selected_record_indices?, complete_rule?, path? }`
  - `InspectRnaReadGeneSupport { report_id, gene_ids, selected_record_indices?, complete_rule?, cohort_filter?, path? }`
  - `RunRnaReadBatchMap { manifest_path, seq_id, seed_feature_id, gene_ids[], out_dir, profile?, input_format?, scope?, origin_mode?, target_gene_ids?, roi_seed_capture_enabled?, seed_filter?, align_config?, report_mode?, align_selection?, complete_rule?, concatemer_settings?, concatemer_limit?, continue_on_error? }`
  - persisted RNA-read reports are regular computational artifacts:
    - stored reports now carry `report_id`, `op_id`, and `run_id`
    - `ListRnaReadReports` surfaces the same `op_id` / `run_id` linkage in
      summary rows
    - GUI lineage projects persisted reports as analysis nodes and can reopen
      the `RNA-read Mapping` workspace directly on the selected report
  - `seed_feature_id` may reference an `mRNA`, `transcript`, `ncRNA`,
    `misc_RNA`, `exon`, `gene`, or `CDS` feature; transcript-template
    admission then follows the selected splicing-scope rules around that seed.
  - implemented profile: `nanopore_cdna_v1`
  - implemented input format: `fasta` (`.fa/.fasta`, optional `.fa.gz/.fasta.gz`; `.sra` must be converted externally in phase-1)
  - default seed/filter constants:
    - `kmer_len=10`
    - `seed_stride_bp=1`
    - `min_seed_hit_fraction=0.30` (bootstrap default; future SNR calibration track can override policy)
    - `min_weighted_seed_hit_fraction=0.05`
    - `min_unique_matched_kmers=12`
    - `min_chain_consistency_fraction=0.40`
    - `max_median_transcript_gap=4.0`
    - `min_confirmed_exon_transitions=1`
    - `min_transition_support_fraction=0.05`
    - weighted-hit definition:
      - `weighted_hit_fraction = sum(1 / occurrence_count(seed_bits)) / tested_kmers`
      - `occurrence_count` is measured inside the active scoped seed index
    - seed pass gate:
      - `raw_hit_fraction >= min_seed_hit_fraction`
      - `AND weighted_hit_fraction >= min_weighted_seed_hit_fraction`
      - `AND unique_matched_kmers >= min(min_unique_matched_kmers, tested_kmers)`
      - `AND chain_consistency_fraction >= min_chain_consistency_fraction`
      - `AND median_transcript_gap <= max_median_transcript_gap`
      - `AND confirmed_transitions >= min_confirmed_exon_transitions`
      - `AND confirmed_transition_fraction >= min_transition_support_fraction`
  - phase-1 seed-span behavior:
    - full-read hashing is always used for every read
    - seed-start density is controlled by `seed_stride_bp`
    - default density is one start per base (`seed_stride_bp=1`)
  - sparse-origin behavior:
    - `origin_mode` accepts `single_gene|multi_gene_sparse` (default
      `single_gene`)
    - `target_gene_ids[]` and `roi_seed_capture_enabled` are persisted in the
      report payload for deterministic follow-up runs
    - `multi_gene_sparse` expands local transcript-template indexing with
      transcripts matched from `target_gene_ids[]`
    - `roi_seed_capture_enabled=true` is currently a deterministic no-op with
      explicit warning in report `warnings[]` until the ROI capture layer is
      implemented
  - report compaction and resume behavior:
    - `report_mode=full` keeps retained top hits exactly as ranked
    - `report_mode=seed_passed_only` keeps a smaller retained subset for later
      inspection/alignment:
      - retained hits that passed the composite seed gate
      - retained hits at or above raw `min hit`
      - counters still remain based on the full stream
    - `checkpoint_path` + `checkpoint_every_reads` writes deterministic JSON
      snapshots (`gentle.rna_read_interpret_checkpoint.v1`) during streaming
    - `resume_from_checkpoint=true` resumes from the checkpoint snapshot and
      fast-forwards already-processed records deterministically
  - phase-2 alignment behavior:
    - `AlignRnaReadReport` loads a persisted report and reprocesses a selected
      retained subset (`all|seed_passed|aligned`)
    - phase-2 progress events now emit once per selected retained row
      (`update_every_reads=1`) so adapters can show visible row-by-row advance
    - GUI/shared-shell default selection is `seed_passed`
    - optional `selected_record_indices[]` (0-based stored `record_index`)
      overrides the selection preset and aligns only the explicit subset
    - `selection=all` remains available when you deliberately want the broader
      rescued-retained working set to receive round-2 similarity/coverage
      scores
    - `selection=aligned` means rerun phase 2 only on retained rows that
      already have a stored mapping from an earlier phase-2 pass
    - if `selection=seed_passed` matches no retained hits and no explicit
      record indices were supplied, the engine falls back to retained rows at
      or above raw `min hit`, and if that is still empty, to the highest
      phase-1 score retained row
    - aligner configuration uses `align_config_override` when supplied,
      otherwise the report-stored `align_config`
    - mapping backend uses `bio::alignment::pairwise::banded` with
      `align_band_bp` as band width (`w`) and transcript-seed `kmer_len` as
      seed length (`k`), plus deterministic dense fallback when the banded
      solver yields no mapping
    - phase-2 pairwise alignment evaluates both query orientations for every
      retained row (stored query plus reverse complement) and keeps the
      best-scoring deterministic candidate, preferring semiglobal over local
      and non-reversed over reversed only as later tie-breakers
    - selected retained rows are pairwise aligned regardless of whether their
      recomputed composite seed-pass flag remains true; the seed-pass result is
      still recomputed and stored independently for later inspection
  - native batch-map behavior:
    - manifest v1 is TSV and requires `sample_id` plus either `input_path` or
      `sra_accession`
    - `input_path` rows currently accept FASTA only (`.fa`, `.fasta`,
      `.fa.gz`, `.fasta.gz`) and run the same engine `InterpretRnaReads` plus
      `AlignRnaReadReport` paths as manual single-sample workflows
    - rows with only `sra_accession` are not fetched by default; the batch
      report marks them `needs_preparation` and writes
      `sra_preparation_plan.tsv` plus `sra_preparation_commands.sh`
    - with `prepare_sra=true`, GENtle runs `ReadAcquirePrepare` first, stores
      the nested `read_acquisition_report`, and maps from the prepared FASTA
      path
    - default batch settings are `report_mode=full`, `align_selection=all`,
      `complete_rule=near`, `max_secondary_mappings=5`, and
      `continue_on_error=true`
    - if `origin_mode` is omitted, one requested gene uses `single_gene` and
      multiple requested/target genes use `multi_gene_sparse`
    - `RunRnaReadBatchMap` writes one bundle under `out_dir`:
      `batch_report.json`, `batch_summary.tsv`,
      `gene_screen_summary.tsv`, `sample_sheet.tsv`, `isoform_support.tsv`,
      `concatemer_partner_summary.tsv`, optional SRA preparation/read-acquisition
      files, and per-sample gene-support / concatemer JSON
    - `gene_screen_summary.tsv` uses schema
      `gentle.rna_read_gene_screen_summary.v1`; it normalizes the native
      batch-map output to the same conservative seed-passed support table used
      by the pancreas gene-screen shell helper, including all-read q90/q95/q99
      length context and seed-passed support-read length summaries
    - `batch_summary.tsv` is the sample-level dashboard substrate; it includes
      status/error/warning fields, read counts/fractions, origin-class counts,
      requested/matched/missing genes, accepted-target and other-gene metrics,
      complete-near/strict/exact counts, concatemer suspicion counts, and JSON
      artifact paths
    - `isoform_support.tsv` aggregates accepted target reads by resolved
      transcript/isoform, including fragment/complete counts, mean length,
      mean identity/coverage, and exon/exon-pair/direct-transition support
      JSON columns
    - `concatemer_partner_summary.tsv` aggregates conservative suspicion
      evidence for recurring partner genes/transcripts; target-plus-partner
      rows remain evidence, not proof of a chimera
    - updated report fields include:
      - per-hit mapping fields (`best_mapping`, `secondary_mappings`)
      - per-hit `msa_eligible` and `msa_eligibility_reason`
      - aggregate `read_count_aligned` and `retained_count_msa_eligible`
      - refreshed seed/path diagnostics
        (`transition_support_rows`, `isoform_support_rows`)
      - refreshed mapped support rows
        (`exon_support_frequencies`, `junction_support_frequencies`,
        `mapped_isoform_support_rows`)
      - `mapped_isoform_support_rows[]` also carries optional
        `dominant_triage_bin` and `triage_bin_counts` fields so adapters can
        show conservative categorical isoform-read support while still
        accepting older reports that lack those fields
      - mapped exon/junction support is derived from aligned transcript-template
        offsets first and falls back to coarse genomic-span overlap only for
        legacy mappings that do not carry template offsets
      - deterministic retained-hit re-ranking by alignment-aware retention rank
    - `export-isoform-triage-tsv` is a conservative read-level triage export:
      it bins retained reads as `known_isoform_confirmed`,
      `known_isoform_ambiguous`, `gene_supported_no_isoform_call`, or
      `off_target_or_bad_seed` from already-computed seed, origin,
      transcript-template alignment, transition, and secondary-mapping fields.
      It deliberately does not call novel isoforms and should not be used as an
      expression heatmap source.
  - alignment inspection behavior:
    - `rna-reads inspect-alignments` accepts coarse `selection` plus a
      structured subset contract:
      - `effect_filter = all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only`
      - `sort_key = rank|identity|coverage|score`
      - `search = free-text match over read ids, transcript ids/labels,
        effect labels, and `#record_index` labels`
      - `selected_record_indices[]` provides the explicit subset for
        `selected_only`
      - `score_density_variant = all_scored|composite_seed_gate|retained_replay_current_controls`
      - optional `score_density_seed_filter_override` carries the current
        seed-gate controls when an adapter requests retained-only replay under
        current controls
      - `score_bin_index` + `score_bin_count` provide a formal
        score-density-bin subset for reproducible histogram-driven inspection
        within that chosen histogram population
    - inspection payload now includes:
      - `aligned_count`: aligned rows admitted by coarse `selection`
      - `subset_match_count`: aligned rows matching the structured subset
        before `limit`
      - `row_count`: returned rows after `limit`
      - `subset_spec`: normalized structured subset object echoed back in the
        response for deterministic replay
    - row `rank` remains the original alignment-aware retention rank even when
      subset sorting reorders the returned rows
  - fragment/concatemer suspicion audit behavior:
    - `InspectRnaReadConcatemers` is a non-mutating saved-report audit aimed at
      diagnosing whether retained cDNA reads look more like fragment fusions
      than coherent full transcript molecules
    - output schema:
      - `gentle.rna_read_concatemer_inspection.v1`
    - settings:
      - `internal_homopolymer_min_bp`
      - `end_margin_bp`
      - `max_primary_query_coverage_fraction`
      - `min_secondary_identity_fraction`
      - `max_secondary_query_overlap_fraction`
      - `adapter_fasta_path`
      - `adapter_min_match_bp`
      - `fragment_min_bp`
      - `fragment_max_parts`
      - `fragment_min_identity_fraction`
      - `fragment_min_query_coverage_fraction`
      - `transcript_fasta_paths[]`
      - `transcript_index_paths[]`
    - current evidence signals:
      - low primary `query_coverage_fraction`
      - internal poly(A) run away from both read ends
      - internal poly(T) run away from both read ends
      - disjoint secondary mappings with limited query overlap
      - phase-1 local-block / partial-origin classification
      - optional internal adapter-like matches sourced from an external FASTA
      - optional iterative fragment decomposition over the admitted
        transcript-template set, with per-read distinct-gene/group counts
    - row payload now also includes:
      - `internal_adapter_hit_count`
      - `top_internal_adapter_label`
      - `top_internal_adapter_match_bp`
      - `fragment_origin_gene_count`
      - `fragment_origin_gene_ids[]`
      - `adapter_hits[]`
      - `fragment_origins[]`
    - top-level payload now also includes:
      - `internal_adapter_match_count`
      - `multi_gene_fragment_count`
    - current severity semantics are intentionally conservative:
      - `strong` requires at least one disjoint secondary mapping plus another
        signal
      - `moderate` requires at least two non-background signals
      - `weak` indicates only one current signal
    - adapters should treat the payload as a suspicion-ranking surface rather
      than as a definitive chimera caller
    - if the source report was aligned with `max_secondary_mappings = 0`, the
      payload includes a warning that the disjoint-secondary branch of the
      audit was unavailable
  - reusable external transcript-catalog index behavior:
    - schema:
      - `gentle.rna_read_transcript_catalog_index.v1`
    - built by the shared-shell command:
      - `rna-reads build-transcript-index OUTPUT.json --transcript-fasta PATH [--transcript-fasta PATH ...] [--kmer-len N]`
    - payload fields:
      - `schema`
      - `generated_at_unix_ms`
      - `seed_kmer_len`
      - `source_paths[]`
      - `transcript_count`
      - `gene_count`
      - `warnings[]`
      - `templates[]`
    - each template row preserves:
      - `transcript_id`
      - `transcript_label`
      - `gene_id`
      - `strand`
      - `sequence`
      - precomputed `kmer_positions`
    - `InspectRnaReadConcatemers` can reuse one or more prepared indexes via
      `transcript_index_paths[]`
    - current scope:
      - this is a reusable prepared transcript-catalog layer for repeated
        concatemer audits
      - it is not yet a heavier fully joined global inverted seed index across
        all external catalogs
  - on-demand pairwise-alignment detail behavior:
    - the engine can reconstruct the exact phase-2 read-vs-transcript-template
      alignment for one retained row from the saved report plus admitted
      transcript-template set
    - detail payload schema:
      - `gentle.rna_read_alignment_detail.v1`
    - payload includes:
      - selected retained row id (`record_index`, `header_id`)
      - transcript/template target identity
      - phase-2 `alignment_mode`
      - alignment backend (`banded` or `dense_fallback`)
      - aligned query/template spans, full template length, score, identity,
        query coverage, transcript-template coverage, and CIGAR
      - aligned `query / relation / target` text rows for manual inspection of
        low-complexity or partial confirmations
  - exact-subset export behavior:
    - `ExportRnaReadHitsFasta`, `ExportRnaReadExonPathsTsv`,
      `ExportRnaReadExonAbundanceTsv`, `ExportRnaReadAlignmentsTsv`, and
      `ExportRnaReadIsoformTriageTsv` accept optional
      `selected_record_indices[]`
    - when present, the explicit 0-based stored `record_index` subset
      overrides the coarse `selection` preset
    - these exports also accept optional `subset_spec`, a human-readable formal
      description such as `filter=... | sort=... | search=...`; when provided,
      the exported artifact records both the explicit `record_index` subset and
      the subset definition that produced it
    - intended for exporting the exact contributor reads surfaced by mapped
      `Audit` actions in the GUI
  - target-gene cohort summary behavior:
    - `SummarizeRnaReadGeneSupport` is non-mutating and consumes one persisted
      aligned RNA-read report
    - required `gene_ids[]` are normalized/deduplicated and matched
      case-insensitively against the same splicing group-label logic already
      used for transcript grouping
    - output schema:
      - `gentle.rna_read_gene_support_summary.v1`
    - base cohort:
      - retained rows with `best_mapping` present
      - optionally intersected with explicit `selected_record_indices[]`
    - accepted target cohort:
      - base-cohort rows whose `best_mapping.transcript_feature_id` resolves to
        one of the requested matched genes/groups
    - complete/fragment split:
      - `complete_rule = near|strict|exact` controls which accepted rows land
        in the `complete` cohort
      - `fragment` is the remaining accepted-target cohort
      - summary still reports nested `complete_strict_count` and
        `complete_exact_count` regardless of the chosen `complete_rule`
    - support attribution is derived from phase-2 mapped support, not phase-1
      `exon_path`
    - per-cohort output blocks:
      - `all_target`
      - `fragments`
      - `complete`
    - each block includes:
      - `read_count`
      - `exon_support[]`
      - `exon_pair_support[]`
      - `direct_transition_support[]`
    - top-level distribution sidecars also summarize target-fragment quality:
      - `evaluated_read_lengths`
        - read-length distribution over the aligned/evaluable base cohort used
          for gene-support classification
      - `accepted_target_read_lengths`
        - total read-length distribution for accepted target-positive reads
      - `accepted_target_fragment_lengths`
        - aligned query-span distribution for the accepted target fragments
      - `accepted_target_query_coverage`
        - fraction-of-read-covered distribution for those accepted target
          fragments (`query_coverage_fraction`)
      - each length summary carries:
        - `sample_count`
        - `mean_length_bp`
        - `min_length_bp`
        - `q25_length_bp`
        - `median_length_bp`
        - `q75_length_bp`
        - `max_length_bp`
        - `p95_length_bp`
        - exact `length_counts[]`
      - each fraction summary carries:
        - `sample_count`
        - `mean_fraction`
        - `median_fraction`
        - `p95_fraction`
        - exact `bin_counts[]`
    - row semantics:
      - `exon_support[]`: each exon counted at most once per read
      - `exon_pair_support[]`: every ordered exon_i -> exon_j pair observed in
        the mapped exon order once per read, including skipped pairs like
        `1->3`
      - `direct_transition_support[]`: neighboring exon steps only, so
        `1->2` is counted but skipped pairs like `1->3` are not
      - all support fractions are normalized by the enclosing cohort size
      - exon and pair rows carry deterministic gene-level exon ordinals plus
        genomic coordinates for auditability
    - when `path` / shell `--output` is provided, the exact same JSON payload
      returned to the caller is also written to disk
    - subordinate-artifact provenance fields:
      - `generated_at_unix_ms`
      - summary op `op_id` / `run_id`
      - source RNA-read report `generated_at_unix_ms` / `op_id` / `run_id`
  - target-gene cohort audit behavior:
    - `InspectRnaReadGeneSupport` is non-mutating and shares the same
      requested-gene matching, selected-record restriction, accepted-target
      logic, and `complete_rule` classification used by
      `SummarizeRnaReadGeneSupport`
    - output schema:
      - `gentle.rna_read_gene_support_audit.v1`
    - evaluation universe:
      - all selected saved-report rows after `selected_record_indices[]`
        filtering, including unaligned retained rows
    - grouped top-level subset handles:
      - `accepted_target_record_indices[]`
      - `fragment_record_indices[]`
      - `complete_record_indices[]`
      - `complete_strict_record_indices[]`
      - `complete_exact_record_indices[]`
    - row status values:
      - `unaligned`
      - `aligned_other_gene`
      - `accepted_fragment`
      - `accepted_complete`
    - row payload includes:
      - `record_index`, `header_id`
      - resolved `gene_id` when available
      - aligned transcript identity (`transcript_feature_id`,
        `transcript_id`, `transcript_label`)
      - machine-readable `status_reason`
      - `full_length_exact`, `full_length_near`, `full_length_strict`, and
        derived `full_length_class`
      - `mapped_exon_ordinals[]`
      - ordered `exon_pairs[]`
      - ordered `direct_transition_pairs[]`
      - phase-2 `score`, `identity_fraction`, `query_coverage_fraction`
      - `passed_seed_filter` as provenance only
    - `cohort_filter = all|accepted|fragment|complete|rejected` limits the
      returned `rows[]` set without changing the grouped top-level subset
      arrays
    - when `path` / shell `--output` is provided, the exact same JSON payload
      returned to the caller is also written to disk
    - subordinate-artifact provenance fields:
      - `generated_at_unix_ms`
      - audit op `op_id` / `run_id`
      - source RNA-read report `generated_at_unix_ms` / `op_id` / `run_id`
- Report persistence:
  - report schema: `gentle.rna_read_report.v1`
  - metadata store schema: `gentle.rna_read_reports.v1`
  - metadata key: `rna_read_reports`
  - `rna-reads list-reports` summary rows include sparse-origin request
    provenance:
    - `origin_mode`
    - `input_orientation_mode` / `input_orientation_label`
      (`cdna_oriented` / `cDNA-oriented` or `direct_rna` / `direct-RNA`,
      derived from `seed_filter.cdna_poly_t_flip_enabled`)
    - `target_gene_count`
    - `roi_seed_capture_enabled`
  - report payload now includes per-report:
    - `exon_support_frequencies[]`
    - `junction_support_frequencies[]`
    - `score_density_bins[]` (`all_scored` phase-1 histogram)
    - `seed_pass_score_density_bins[]` (`composite_seed_gate` histogram)
    - exact read-length histograms (`length_bp -> count`) for
      deterministic subset auditing:
      - `read_length_counts_all`
      - `read_length_counts_seed_passed`
      - `read_length_counts_aligned`
      - `read_length_counts_full_length_exact`
      - `read_length_counts_full_length_near`
      - `read_length_counts_full_length_strict`
      - checkpoint snapshots mirror these vectors so resume/restart keeps
        histogram accumulation deterministic
    - storage/streaming controls:
      - `report_mode` (`full` or `seed_passed_only`)
      - `checkpoint_path` / `checkpoint_every_reads`
      - `resumed_from_checkpoint`
    - request provenance fields:
      - `origin_mode`
      - `target_gene_ids[]`
      - `roi_seed_capture_enabled`
    - `origin_class_counts` (running/final deterministic class tallies)
  - per-hit payload now includes:
    - `origin_class`
    - `origin_reason`
    - `origin_confidence`
    - `strand_confidence`
    - `origin_candidates[]` (selected/plus/minus/seed-chain candidate hints)
    - `best_mapping.alignment_mode` (`semiglobal` preferred, with deterministic
      local fallback when quality is better)
    - `best_mapping.query_reverse_complemented` (whether phase-2 had to
      reverse-complement the stored read to fit the chosen transcript-template
      mapping)
  - alignment inspection payload schema:
    - `gentle.rna_read_alignment_inspection.v1`
    - produced by non-mutating shared-shell inspection command
      `rna-reads inspect-alignments`
    - each row now carries:
      - phase-1 transcript-assignment fields
        (`phase1_primary_transcript_id`, `seed_chain_transcript_id`,
        `exon_path_transcript_id`, `transcript_exon_path`, `exon_path`,
        `exon_transitions_confirmed/total`, `selected_strand`,
        `reverse_complement_applied`)
      - phase-2 best-mapping fields
        (`transcript_id`, `transcript_label`, `strand`, `alignment_mode`,
        `target_start_1based`, `target_end_1based`, `target_length_bp`,
        `identity_fraction`, `query_coverage_fraction`,
        `target_coverage_fraction`, `score`, `secondary_mapping_count`)
      - full-length classification flags (derived deterministically from
        transcript-template coverage and current alignment threshold):
        - `full_length_exact` (`100%` template coverage)
        - `full_length_near` (`>=95%` template coverage)
        - `full_length_strict`
          (`near` + both template ends within `15 bp` + identity above
          active alignment threshold)
      - deterministic comparison field `alignment_effect`
        (`confirmed_assignment`, `reassigned_transcript`,
        `aligned_without_phase1_assignment`)
      - mapped-support attribution arrays for the best mapping
        (`mapped_exon_support[]`, `mapped_junction_support[]`)
    - top-level payload now also carries:
      - `aligned_count`
      - `subset_match_count`
      - `row_count`
      - `limit`
      - normalized `subset_spec`
        (`effect_filter`, `sort_key`, `search`,
        `selected_record_indices[]`, `score_density_variant`,
        `score_bin_index`, `score_bin_count`)
- Sample-sheet export:
  - operation: `ExportRnaReadSampleSheet { path, seq_id?, report_ids?, gene_ids?, complete_rule?, append? }`
  - export schema: `gentle.rna_read_sample_sheet_export.v1`
  - output: TSV with run/read metrics, sparse-origin request provenance
    (`report_mode`, `input_orientation_mode`, `input_orientation_label`,
    `origin_mode`, `target_gene_count`, `target_gene_ids_json`,
    `roi_seed_capture_enabled`), JSON-serialized exon/junction frequency
    columns, and `origin_class_counts_json` for cohort-level downstream
    analysis.
  - additional per-report columns include `mean_read_length_bp`; when one or
    more `gene_ids[]` are requested the same row also carries
    accepted-target counts/fractions, fragment vs complete counts,
    `gene_support_mean_assigned_read_length_bp`, and JSON-serialized
    exon / exon-pair / direct-transition support tables for the requested gene
    cohort.
- Target-quality comparison export:
  - operation:
    `ExportRnaReadTargetQuality { report_id, path, gene_ids, complete_rule? }`
  - comparison bundle schema:
    `gentle.rna_read_target_quality_comparison_bundle.v1`
  - operation result schema:
    `gentle.rna_read_target_quality_export.v1`
  - entry semantics:
    - each exported entry stores:
      - `gentle_version`
      - source report/profile/scope/origin metadata
      - requested/matched/missing gene ids
      - the all-read length distribution used for context
      - the full shared `RnaReadGeneSupportSummary` target-quality payload
    - entries are keyed deterministically so repeated exports with the same
      source report + target + settings replace their prior slot instead of
      duplicating it
  - path behavior:
    - `.json` writes the comparison bundle directly
    - if a `.json` file already exists and is:
      - an existing comparison bundle: merge/append into it
      - a legacy single `RnaReadGeneSupportSummary`: wrap it into a new bundle
        and append the new entry
    - `.svg` writes a rendered comparison figure plus a same-basename
      `.bundle.json` sidecar bundle
    - if the requested SVG path already exists without a reusable GENtle
      sidecar bundle, the export preserves that file and writes a sibling
      `*_compare.svg` plus sidecar instead of overwriting
- Shared-shell command family:
  - `reads acquire status MANIFEST.tsv --cache-dir DIR --work-dir DIR`
  - `reads acquire prepare MANIFEST.tsv --cache-dir DIR --work-dir DIR [--analysis-format fasta|fastq] [--read-layout single_end|paired_end|split_spot] [--threads N] [--max-size SIZE] [--min-free-gb N] [--drop-intermediate-fastq] [--continue-on-error]`
  - `reads acquire inspect RUN_ACCESSION --cache-dir DIR --work-dir DIR`
  - `reads acquire cancel RUN_ACCESSION --cache-dir DIR --work-dir DIR`
  - `rna-reads interpret SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile nanopore_cdna_v1] [--format fasta] [--scope all_overlapping_any_strand|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--seed-stride-bp N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
  - `rna-reads batch-map MANIFEST.tsv --seq-id SEQ_ID --seed-feature-id FEATURE_ID --gene GENE_ID [--gene GENE_ID ...] --out-dir OUT [--target-gene GENE_ID]... [--origin-mode single_gene|multi_gene_sparse] [--report-mode full|seed_passed_only] [--align-selection all|seed_passed|aligned] [--complete-rule near|strict|exact] [--max-secondary-mappings N] [--continue-on-error|--fail-fast] [--prepare-sra] [--read-cache-dir DIR] [--read-work-dir DIR] [--drop-intermediate-fastq] [--transcript-fasta PATH]... [--transcript-index PATH]...`
  - `rna-reads align-report REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
  - `rna-reads preflight-isoforms SEQ_ID FEATURE_ID [--scope all_overlapping_any_strand|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--positive-transcript-fasta PATH ...] [--must-pass-transcript-fasta PATH ...] [--control-transcript-fasta PATH ...] [--optimize-parameters] [--max-control-match-probability F] [seed filter options]`
  - `rna-reads list-reports [SEQ_ID]`
  - `rna-reads show-report REPORT_ID`
  - `rna-reads show-alignment REPORT_ID RECORD_INDEX`
  - `rna-reads show-alignments REPORT_ID [--gene GENE_ID ... --cohort all|accepted|fragment|complete|rejected --complete-rule near|strict|exact | --record-indices i,j,k] [--limit N] [--output PATH]`
  - `rna-reads summarize-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--output PATH]`
  - `rna-reads inspect-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--cohort all|accepted|fragment|complete|rejected] [--output PATH]`
  - `rna-reads inspect-alignments REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--effect-filter all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only] [--sort rank|identity|coverage|score] [--search TEXT] [--record-indices i,j,k] [--score-bin-variant all_scored|composite_seed_gate] [--score-bin-index N] [--score-bin-count M]`
  - `rna-reads inspect-concatemers REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--internal-homopolymer-min-bp N] [--end-margin-bp N] [--max-primary-query-cov F] [--min-secondary-identity F] [--max-secondary-query-overlap F] [--adapter-fasta PATH] [--adapter-min-match-bp N] [--fragment-min-bp N] [--fragment-max-parts N] [--fragment-min-identity F] [--fragment-min-query-cov F] [--transcript-fasta PATH]... [--transcript-index PATH]...`
  - `rna-reads build-transcript-index OUTPUT.json [--kmer-len N] --transcript-fasta PATH [--transcript-fasta PATH ...]`
  - `rna-reads allele-hash-screen --gene GENE --transcript-fasta PATH (--variant-table PATH | --vcf PATH --transcript-map PATH [--vcf-sample SAMPLE]) [--from-rna-report REPORT_ID] [--read-file PATH ...] [--read-pair R1,R2 ...] [--salmon-unmapped-names PATH] [--salmon-mappings-sam PATH] [--read-id-allowlist PATH] [--kmer-len N] [--min-unique-kmer-hits N] [--min-informative-reads N] [--balanced-band-lo F] [--balanced-band-hi F] [--max-inline-read-calls N] --out OUT_DIR`
  - `rna-reads materialize-hits REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--output-prefix PREFIX]`
  - `rna-reads export-report REPORT_ID OUTPUT.json`
  - `rna-reads export-hits-fasta REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-sample-sheet OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--gene GENE_ID]... [--complete-rule near|strict|exact] [--append]`
  - `rna-reads export-target-quality REPORT_ID OUTPUT.{svg|json} --gene GENE_ID [--gene GENE_ID ...] [--complete-rule near|strict|exact]`
  - `rna-reads export-paths-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-abundance-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-dexseq-annotation-gff REPORT_ID OUTPUT.gff`
  - `rna-reads export-dexseq-counts-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads verify-dexseq REPORT_ID OUTPUT.gff OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT] [--r-library-path PATH ...]`
  - `rna-reads export-score-density-svg REPORT_ID OUTPUT.svg [--scale linear|log] [--variant all_scored|composite_seed_gate]`
  - `rna-reads export-alignments-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-alignment-dotplot-svg REPORT_ID OUTPUT.svg [--selection all|seed_passed|aligned] [--max-points N]`
  - shell output convenience fields:
    - `rna-reads list-reports` includes `summary_rows[]` with concise
      human-readable provenance lines (`mode`, `origin`, target count,
      ROI-capture flag, read counters)
    - `rna-reads show-report` includes `summary` with the same provenance
      framing for one report
    - `rna-reads show-alignment` returns a `gentle.rna_read_alignment_display.v1`
      wrapper with:
      - `report_id`
      - exact saved `record_index`
      - `alignment` = the portable `RnaReadAlignmentDisplay` detail record used
        by the GUI `Show alignment` pane
    - `rna-reads show-alignments` returns a
      `gentle.rna_read_alignment_display_batch.v1` wrapper with:
      - `report_id`, `seq_id`, and `seed_feature_id`
      - `selection_mode = record_indices|gene_cohort`
      - requested/matched/missing gene ids for cohort mode
      - `cohort_filter`, `complete_rule`, `limit`, and final
        `selected_record_indices`
      - `entries[]` containing `record_index`, `header_id`, and the same
        `RnaReadAlignmentDisplay` payload as repeated `show-alignment`
      - `skipped_records[]` and `skipped_record_indices[]` for reads without
        `best_mapping`
    - `rna-reads summarize-gene-support` returns the full
      `gentle.rna_read_gene_support_summary.v1` payload directly, including
      `requested_gene_ids`, `matched_gene_ids`, `missing_gene_ids`,
      selected-record echo fields, target-fragment quality distributions,
      per-cohort support tables, and both the summary op plus source-report
      provenance ids/timestamps
    - `rna-reads inspect-gene-support` returns the full
      `gentle.rna_read_gene_support_audit.v1` payload directly, including
      grouped cohort record-index arrays plus row-level `status`,
      `status_reason`, full-length fields, mapped exon/junction audit data,
      and both the audit op plus source-report provenance ids/timestamps
    - `rna-reads batch-map` returns the full
      `gentle.rna_read_batch_map_report.v1` payload and writes the durable
      bundle files listed above; successful FASTA rows also persist ordinary
      RNA-read reports in engine state under deterministic `rna_batch_*`
      report ids unless a manifest `report_id` column overrides them
    - `rna-reads inspect-alignments` returns aligned rows ranked by
      alignment-aware retention score (mapping + seed metrics), plus a
      structured `subset_spec` payload (`effect_filter`, `sort_key`, `search`,
      `selected_record_indices`, `score_density_variant`, `score_bin_index`,
      `score_bin_count`) and
      `subset_match_count`
    - `rna-reads inspect-concatemers` returns the full
      `gentle.rna_read_concatemer_inspection.v1` payload directly, including
      machine-readable signal counts, the normalized thresholds used for the
      audit, ranked suspicious rows, per-read adapter hits and fragment-origin
      explanations, repeated non-primary partner-gene/transcript summaries
      across the suspicious cohort, explicit `selected_record_indices` /
      `subset_match_count` provenance for exact subset reruns, and warning
      text when secondary-mapping evidence was unavailable
    - repeated `--transcript-fasta` arguments append multiple transcript
      catalogs to one audit run; the serialized settings now preserve those as
      `transcript_fasta_paths[]`, while legacy single-path payloads still
      deserialize via `transcript_fasta_path`
    - `rna-reads preflight-isoforms` returns the full
      `gentle.rna_read_isoform_preflight.v1` payload directly, including
      target transcript seed-pass rows, external must-pass positive transcript
      rows, separate grouped negative-control summaries
      (`TP53`, `TP63`, or inferred symbols from supplied FASTA headers/paths),
      the selected `recommended_seed_filter`,
      `threshold_recommendation` with control-margin diagnostics plus
      paste-ready `seed_filter_cli_fragment` / `interpret_command_fragment`,
      and a reproducible
      `recommended_command_fragment` that preserves every
      `--positive-transcript-fasta` / `--control-transcript-fasta` path
    - `fragment_max_parts = 0` is now a supported "signals only" mode for the
      concatemer audit:
      - it disables fragment decomposition entirely
      - adapter/homopolymer/disjoint-secondary/phase-1 signals still run
      - intended use: first-pass triage before a narrower transcriptome-backed
        rerun on selected `record_index` values
    - `rna-reads build-transcript-index` returns the full
      `gentle.rna_read_transcript_catalog_index.v1` payload directly and also
      writes the same JSON to the requested output path
    - `rna-reads allele-hash-screen` returns the full
      `gentle.rna_allele_hash_screen.v3` payload directly and writes the same
      JSON, a read-call TSV, and reference/hap1/hap2 transcript FASTA files
      under the requested output directory. The deterministic path accepts
      an explicit transcript-coordinate variant TSV with columns
      `transcript_id`, `cdna_pos_1based`, `ref`, `alt`, and `genotype`;
      optional columns are `variant_id`/`id` and `phase_set`. Phased
      genotypes (`0|1`, `1|0`) materialize hap1/hap2 transcripts, while
      unphased slash genotypes are labelled `unphased_allele_level_only`
      rather than inventing phase. Distinct phase sets remain independent
      `phase_blocks`; global haplotype FASTA output is marked `inferred=false`
      unless one coherent phased block supports it. The alternative VCF path
      accepts a local VCF plus a transcript-map TSV with `chrom`,
      `genomic_pos_1based`, `transcript_id`, and `cdna_pos_1based` columns
      (optional `strand`). Only explicit PASS biallelic SNVs are projected;
      `--vcf-sample` is required when the VCF contains multiple samples, and
      every projected reference allele must match the transcript sequence.
      `--from-rna-report REPORT_ID` resolves the same accepted target-gene
      cohort as `rna-reads inspect-gene-support`: retained rows must have a
      `best_mapping` assigned to the requested gene/group, while seed-pass
      state is provenance rather than an acceptance gate. Their stored
      sequences are screened directly, so no extracted FASTA/FASTQ is needed
      for that source. `--salmon-unmapped-names` and
      `--salmon-mappings-sam` reuse the target-rescue ID parsers to select
      Salmon-unassigned and target-transcript-mapped IDs from explicit
      `--read-file`/`--read-pair` sequence inputs. Salmon ID files alone are
      rejected because they do not provide every selected read sequence.
      Repeated `--read-pair R1,R2` inputs are streamed in lockstep and counted
      as fragments, with an invalid-base boundary preventing cross-mate k-mers.
      Report fields include `schema`, `gene`, `phase_mode`, `params`,
      physical read/fragment/evidence-observation counts, `phase_blocks`,
      `gene_representation`, optional `rna_report_expression_support`,
      `output_files`, `haplotype_fastas`, `transcript_summaries`,
      `variant_summaries`, `classification_counts`, `source_provenance[]`,
      `reads[]`, and `warnings[]`. Each read call carries
      `source_origins[]`; provenance counts can overlap when an explicit read
      is also selected into a Salmon cohort. Reports using the v1 shape still
      deserialize with empty source-provenance fields. `reads[]` is capped by
      `--max-inline-read-calls` (default
      10,000), while the streamed TSV and aggregate counts remain complete.
      Read calls classify evidence as `hap1`, `hap2`, `alternate`,
      `reference_only`, `ambiguous`, `uninformative`, or `off_target` and carry
      matched transcript ids, supporting variant ids, and local phase-block
      calls. V3 attaches an optional `representation` assessment to every
      variant, transcript, and phase-block summary and one top-level gene
      assessment. The assessment records the verdict, hap1 fraction,
      haplotype-informative count, expression-weight basis, applied thresholds,
      coverage/phase/depth caveats, and an optional advisory two-sided binomial
      value under p=0.5. Defaults are 10 informative observations and an
      inclusive balanced band of 0.40 through 0.60. Unphased units never become
      haplotype calls; mixed roll-ups exclude them with an explicit caveat.
      Roll-ups across multiple phased blocks retain the requested aggregate
      verdict but add `block_local_phase_labels_aggregated`, because block-local
      labels do not reconstruct a linked gene-wide haplotype.
      When sourced from an RNA report, retained target-gene/transcript support
      is attached as expression context while the operative weight remains the
      screen's own allele-informative depth; support below the configured
      informative-depth floor adds a non-blocking caveat. V1/v2 payloads
      deserialize with no invented assessment. The report is sequence evidence
      only and does not claim biological allelic imbalance, significance, or
      causation.
    - `rna-reads materialize-hits` returns a
      `gentle.rna_read_hit_materialization.v1` wrapper with:
      - the mutating `result` (`OpResult`) from
        `MaterializeRnaReadHitSequences`
      - `created_sequence_count`
      - `created_sequences[]` with per-sequence `seq_id`, `name`, and
        `length_bp`
    - `rna-reads export-target-quality` returns the full
      `gentle.rna_read_target_quality_export.v1` payload directly, including
      `requested_path`, `written_path`, optional `bundle_path`, `entry_count`,
      and append/reuse flags so adapters can tell whether the export extended
      an existing comparison artifact or started a fresh one
- Alignment-detail inspection:
  - shared shell:
    `rna-reads show-alignment REPORT_ID RECORD_INDEX`
  - output:
    - exact saved-report pairwise alignment detail for one retained read
    - non-mutating; intended for cluster/headless review of the same
      read-vs-template evidence the GUI exposes under `Show alignment`
- Alignment-TSV export:
  - operation:
    `ExportRnaReadAlignmentsTsv { report_id, path, selection, limit?, selected_record_indices?, subset_spec? }`
  - export schema: `gentle.rna_read_alignment_tsv_export.v1`
  - output: ranked alignment rows as TSV with:
    - leading `#` metadata lines for report provenance (`selection`, `limit`,
      `selected_record_indices`, `subset_spec`, `profile`, `scope`, `origin_mode`)
    - seed-screen sampling/gating context (`k`, `seed_stride_bp`,
      overlap/order-density wording, seed thresholds)
    - alignment config summary (`min_identity_fraction`,
      `max_secondary_mappings`)
    - phase-1 transcript/path diagnostics
    - phase-2 mapping metrics
    - `alignment_effect`
    - compact mapped exon/junction attribution columns
    - optional top-`N` truncation via `limit`
- Score-density SVG export:
  - `rna-reads export-score-density-svg` writes the same report summary used by
    the GUI plus seed-screen provenance in the SVG header:
    - `variant = all_scored|composite_seed_gate|retained_replay_current_controls`
    - `profile`, `report_mode`, `scope`, `origin_mode`
    - seed-filter summary with `k`, `seed_stride_bp`, thresholds, and
      overlap/order-density wording
    - optional `replay_seed_filter` summary when the export uses retained-only
      replay under current controls
    - whether bins were stored in the report or derived from retained hits
- Alignment-dotplot export:
  - operation:
    `ExportRnaReadAlignmentDotplotSvg { report_id, path, selection, max_points }`
  - export schema: `gentle.rna_read_alignment_dotplot_svg_export.v1`
  - output: SVG scatter of query coverage vs identity for aligned hits with
    score-colored points and report-config threshold guide.
- Read-sequence materialization:
  - shared shell:
    `rna-reads materialize-hits REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--output-prefix PREFIX]`
  - operation:
    `MaterializeRnaReadHitSequences { report_id, selection, selected_record_indices?, output_prefix? }`
  - output:
    - creates one ordinary project sequence per selected retained RNA-read hit
    - exact `selected_record_indices` takes precedence over coarse `selection`
    - intended for downstream dotplots/manual inspection of saved-report
      outliers without re-reading the FASTA input
- `rna-reads export-hits-fasta` header extensions:
  - optional `selected_record_indices[]` overrides the coarse selection preset
    for exact saved-report subset export
  - optional `subset_spec` records the formal subset definition that produced
    that explicit `record_index` subset
  - `exon_path_tx=<transcript_id|none>`
  - `transcript_exon_path=<ordinal_path|none>` uses transcript-local exon
    numbering and is the primary human-facing path; `:` marks hash-confirmed
    adjacent exon transitions and `-` marks unconfirmed adjacency
  - `exon_path=<ordinal_path|none>` is the DEXSeq-style genomic projection into
    disjoint exonic-part bins, numbered by ascending genomic coordinate; `_`
    marks adjacent bins inside the same transcript-local exon (no splice),
    while `:` / `-` retain confirmed / unconfirmed splice junction semantics
  - schema note: these global projection ordinals are renumbered relative to
    pre-exonic-parts releases; `transcript_exon_path` is unchanged and remains
    the primary human-facing field
  - `exon_transitions=<confirmed>/<total>`
  - `rc_applied=<true|false>` (automatic cDNA poly-T reverse-complement
    normalization marker)
  - `origin_class=<...>` plus `origin_conf=<...>` and `strand_conf=<...>`
- `rna-reads export-exon-paths-tsv` and `rna-reads export-exon-abundance-tsv`
  now begin with the same `#` report/seed-screen provenance block used by the
  alignment TSV export, minus alignment-only fields; optional `subset_spec`
  records the formal subset definition alongside `selected_record_indices`.
  Abundance exports count per-bin support from `exon_path` and exclude `_`
  intra-exon adjacencies from junction rows.
- DEXSeq interoperability uses one persisted join table on
  `RnaReadInterpretationReport.exonic_part_bins[]`:
  - each bin stores `global_ordinal`, aggregate `gene_id`, per-aggregate-gene
    `exonic_part_number`, one-based inclusive coordinates, strand, supporting
    transcript IDs, and the existing partition's `constitutive` flag
  - `global_ordinal` is exactly the integer encoded in `hit.exon_path`; part
    numbers restart at `1` for each aggregate gene and increase by genomic
    coordinate even on minus-strand loci
  - aggregate gene IDs are the sorted `+`-joined set of genes whose
    transcripts support the bin; unresolved genes use a stable
    `GENtle_unassigned_<seq_id>` identifier
  - reports written before this additive field deserialize normally with an
    empty table, but both DEXSeq exporters reject them with a request to re-run
    alignment
- `ExportRnaReadDexseqAnnotationGff` /
  `rna-reads export-dexseq-annotation-gff` write a selection-independent
  flattened GFF with `aggregate_gene` and `exonic_part` rows. Exonic-part
  attributes are `transcripts`, zero-padded `exonic_part_number`, and
  `gene_id`. The result schema is
  `gentle.rna_read_dexseq_annotation_gff_export.v1`.
- `ExportRnaReadDexseqCountsTsv` /
  `rna-reads export-dexseq-counts-tsv` write the corresponding strict
  two-column HTSeq-style count file. Ordinary row IDs are
  `<gene_id>:E<zero-padded-three-digit-part-number>` and are generated from the
  same persisted table as the GFF. A selected retained read contributes once
  per touched ordinal after per-read deduplication. No metadata comments or
  extra columns are written; selection, explicit record indices, and
  `subset_spec` are carried by the JSON result schema
  `gentle.rna_read_dexseq_counts_tsv_export.v1`.
  Both export result schemas report the total written `row_count` alongside
  aggregate-gene and exonic-part counts.
  - `_empty` counts selected rows with empty `exon_path`
  - `_lowaqual` counts selected rows that failed the seed gate
  - `_notaligned` counts selected rows without `best_mapping`
  - `_ambiguous` and `_ambiguous_readpair` are currently `0` because those
    HTSeq categories are not retained by GENtle
  - special-row diagnostics are independently observed and may overlap
  - use the per-sample count files and shared flattened GFF with
    `DEXSeqDataSetFromHTSeq(countfiles, sampleData, design, flattenedfile)`
- `VerifyRnaReadDexseqExports` / `rna-reads verify-dexseq` call both existing
  exporters, then preflight `Rscript` and the Bioconductor `DEXSeq` package
  without installing anything. When both are present, the bounded
  `scripts/rna_read_dexseq_verify.R` helper constructs a real
  `DEXSeqDataSetFromHTSeq()` from the two files. Result schema
  `gentle.rna_read_dexseq_verification.v1` carries both export summaries,
  dependency rows, requested/effective R library paths, the reproducible
  verifier command, `verifier_status`, and an optional `DEXSEQ_OK` stdout
  summary. Missing dependencies and timeouts remain structured inspection
  outcomes rather than being mistaken for an accepted DEXSeq pair.
- cDNA/direct-RNA normalization controls in `seed_filter`:
  - `cdna_poly_t_flip_enabled` (default `true`)
  - `poly_t_prefix_min_bp` (default `18`): minimum T support used by the
    tolerant 5' poly-T-head detector (minor interruptions in the head are
    accepted)
  - report-list/detail/export surfaces derive stable presentation markers from
    `cdna_poly_t_flip_enabled`: `input_orientation_mode=cdna_oriented` with
    label `cDNA-oriented`, or `input_orientation_mode=direct_rna` with label
    `direct-RNA`
- Scope/strand semantics for `InterpretRnaReads`:
  - `all_overlapping_any_strand`: all overlapping transcripts on any strand,
    including antisense/opposite-strand genes relative to the selected target
    gene/group
  - `target_group_any_strand`: target-group transcripts only, any annotated
    strand allowed
  - `all_overlapping_target_strand`: all overlapping transcripts on the
    selected target gene/group's annotated strand only
  - `target_group_target_strand`: target-group transcripts on the selected
    target gene/group's annotated strand only
  - scoring note:
    any-strand modes score against the union of admitted target-gene-strand and
    antisense/opposite-strand templates; target-gene-strand modes exclude
    antisense/opposite-strand templates.
  - seed-index note:
    indexed seeds include annotated exon-body and exon-exon transition k-mers
    for admitted transcripts.

Async BLAST shell contract (agent/MCP-ready baseline):

- Shared-shell families (both `genomes` and `helpers` scopes):
  - `blast-start GENOME_ID QUERY_SEQUENCE ...`
  - `blast-status JOB_ID [--with-report]`
  - `blast-cancel JOB_ID`
  - `blast-list`
- Deterministic job payload schemas:
  - `gentle.blast_async_start.v1`
  - `gentle.blast_async_status.v1`
  - `gentle.blast_async_cancel.v1`
  - `gentle.blast_async_list.v1`
- External-binary preflight payload:
  - `blast-start` responses now include `binary_preflight` with schema
    `gentle.blast_external_binary_preflight.v1`.
  - payload includes deterministic `blastn` and `makeblastdb` probe rows with:
    `found`, `version`, `executable`, and resolved `path` diagnostics.
  - equivalent preflight payload is also emitted by synchronous shared-shell
    routes `prepare`, `blast`, and `blast-track`.
- Job status contract:
  - `job_id` stable per process
  - non-terminal states: `queued | running`
  - terminal states: `completed | failed | cancelled`
  - scheduler metadata:
    - `max_concurrent_jobs`
    - `running_jobs`
    - `queued_jobs`
    - `queue_position` (present while state is `queued`)
  - optional final `report` on `blast-status --with-report`
- Durability/restart semantics:
  - BLAST async status snapshots are persisted in project metadata as
    `blast_async_jobs` (`gentle.blast_async_job_store.v1`).
  - On restart/reload, recovered jobs that were previously non-terminal but no
    longer have an active worker context are normalized deterministically:
    - `cancel_requested=true` -> `cancelled`
    - otherwise -> `failed` with explicit restart/reload interruption reason.
  - `blast-start`, `blast-status`, `blast-cancel`, and `blast-list` may mark
    shell state as changed when they persist updated async job snapshots.
- Introspection:
  - Read-only async BLAST status/list rows (`genomes blast-status`,
    `helpers blast-status`, `genomes blast-list`, `helpers blast-list`,
    `blast_async_status`, `blast_async_list`) are fact-annotated with no
    project-state preconditions.
  - These rows may refresh async-job metadata while polling/listing, but they
    declare no hard biological project effects; job ids and optional terminal
    reports remain execution-time validation concerns.
  - Async BLAST start/cancel and synchronous BLAST execution remain
    registry-only until their external-binary, prepared-index, cancellation,
    and report/materialization semantics are modeled explicitly.
- Scheduler policy:
  - async BLAST jobs are executed by a bounded FIFO scheduler (queue + worker slots)
  - default concurrency uses host CPU parallelism
  - optional override via environment variable
    `GENTLE_BLAST_ASYNC_MAX_CONCURRENT` (clamped to `1..256`)
- `gentle_mcp` exposes equivalent tool routes:
  - `blast_async_start`
  - `blast_async_status`
  - `blast_async_cancel`
  - `blast_async_list`

### Workflow

```json
{
  "run_id": "string",
  "ops": ["Operation", "Operation", "..."]
}
```

Notes:

- Splicing Expert `Nanopore cDNA interpretation` uses this same workflow shape
  when you click `Copy Workflow JSON`.
- `Prepare Workflow Op` in the same panel writes `run_id`/`ops` into the GUI
  workflow runner so the exact `InterpretRnaReads` payload can be rerun through
  the generic workflow path.

### OpResult

```json
{
  "op_id": "op-1",
  "created_seq_ids": ["..."],
  "changed_seq_ids": ["..."],
  "warnings": ["..."],
  "messages": ["..."],
  "genome_annotation_projection": null,
  "sequence_alignment": null
}
```

### Error

```json
{
  "code": "InvalidInput|NotFound|Unsupported|Io|Internal",
  "message": "human-readable explanation"
}
```

## State model

`gentle_cli` persists engine state in JSON (`.gentle_state.json` by default).

This supports:

- resumable multi-step workflows
- external inspection
- reproducibility and audit trails

## Recommended AI-agent flow

1. Query `capabilities`
2. Import or initialize state
3. Apply one operation at a time, checking warnings/errors
4. Save/export artifacts
5. Optionally export final state for handoff

## Planned next additions

- richer sequence-editing and annotation operation set
- ligation protocol presets with sticky/blunt compatibility derivation
- render/view model endpoint for frontend-independent graphical representation
- schema publication for strict client-side validation
- CRISPR guide-design next phase:
  - off-target search/ranking contracts
  - on-target efficacy model integration hooks
  - guide-design macro/template expansion into deterministic `Workflow` JSON
  - see draft: `docs/rna_guides_spec.md`
