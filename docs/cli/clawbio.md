# ClawBio and OpenClaw

<!-- manual-counterparts:start -->
Related GUI reference: [GENtle Shell (GUI)](../gui/gentle-shell-gui.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## ClawBio/OpenClaw skill scaffold

Short skill descriptions select workflows; they do not certify them. Shared
`introspect capabilities` includes optional capability-maturity discovery, and
ClawBio `skill-info`/`intents` advertise the generated sidecar without launching
GENtle. `capabilities` keeps its existing meaning. Missing or stale metadata
remains unassessed and does not disable normal use. See the
[evidence and compatibility contract](../capability_maturity.md).

GENtle ships a copy-ready ClawBio skill scaffold at:

- `integrations/clawbio/skills/gentle-cloning/`
- ready-to-paste catalog object:
  `integrations/clawbio/skills/gentle-cloning/catalog_entry.json`
- included local-checkout launcher:
  `integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh`
- included Apptainer/Singularity helper launcher:
  `integrations/clawbio/skills/gentle-cloning/gentle_apptainer_cli.sh`

Purpose:

- expose deterministic GENtle CLI routes as a ClawBio skill,
- position GENtle as a bridge from biological observations, including
  patient/cohort signals, to sequence-grounded mechanistic follow-up and
  wet-lab planning,
- make clear that prepared Ensembl/reference assets and BLAST-capable indices
  are reusable local infrastructure rather than GENtle-only artifacts,
- keep reproducibility artifacts (`report.md`, `result.json`,
  `reproducibility/*`) for each run.

Logical capability split inside that one runtime alias:

- runtime and resource readiness
- genomic context
- TFBS analysis
- restriction analysis
- splicing expert
- isoform architecture
- protein isoform gel / 2D-gel rendering
- experimental follow-up

This is a ClawBio-facing routing/documentation split, not separate
executables. It exists so wrappers and tutorials can call out the real shared
GENtle command surfaces already present behind the single `gentle-cloning`
skill.

Recommended OpenClaw framing for broad discovery questions:

- `observation -> mechanistic hypothesis -> GENtle sequence/context analysis -> wet-lab validation plan`
- GENtle does not prove causality by itself; it translates a prioritized signal
  into sequence-grounded inspection and reproducible follow-up artifacts.

Recommended first-time route for ClawBio/BioClaw users:

```bash
export GENTLE_REPO_ROOT=/home/clawbio/GENtle
export GENTLE_CLI_CMD=/home/clawbio/ClawBio/skills/gentle-cloning/gentle_local_checkout_cli.sh

cd /home/clawbio/ClawBio
python clawbio.py run gentle-cloning --demo
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_intents_runtime.json --output /tmp/gentle_intents_runtime
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_runtime_version.json --output /tmp/gentle_runtime_version
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_list_human.json --output /tmp/gentle_list_human
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_list_gst.json --output /tmp/gentle_list_helpers
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_hosts_list_deor.json --output /tmp/gentle_list_hosts
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_ensembl_available_human.json --output /tmp/gentle_ensembl_human
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_install_ensembl_mouse.json --output /tmp/gentle_install_mouse
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_shell_state_summary.json --output /tmp/gentle_state_summary
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_services_telegram_guide.json --output /tmp/gentle_telegram_guide
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_status_grch38.json --output /tmp/gentle_status_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_prepare_grch38.json --output /tmp/gentle_prepare_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_blast_grch38_short.json --output /tmp/gentle_grch38_blast
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_prepare_puc19.json --output /tmp/gentle_prepare_puc19
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genbank_fetch_pbr322.json --output /tmp/gentle_fetch_pbr322
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_dbsnp_fetch_rs9923231.json --output /tmp/gentle_fetch_rs9923231
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_extract_gene_tp53.json --output /tmp/gentle_extract_tp53
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_export_bed_grch38_tp53_gene_models.json --output /tmp/gentle_tp53_bed
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_blast_puc19_short.json --output /tmp/gentle_puc19_blast
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_vkorc1_planning.json --output /tmp/gentle_vkorc1_planning
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_svg_pgex_fasta_circular.json --output /tmp/gentle_pgex_map
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_export_bed_pgex_fasta_tfbs_restriction.json --output /tmp/gentle_pgex_bed
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_svg_pgex_fasta_linear_tfbs.json --output /tmp/gentle_pgex_tfbs_map
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_tfbs_summary_pgex_fasta.json --output /tmp/gentle_pgex_tfbs_summary
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_inspect_feature_expert_pgex_fasta_tfbs.json --output /tmp/gentle_pgex_tfbs_expert
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_feature_expert_pgex_fasta_tfbs_svg.json --output /tmp/gentle_pgex_tfbs_expert_svg
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_svg_pgex_fasta_linear_restriction.json --output /tmp/gentle_pgex_restriction_map
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_inspect_feature_expert_pgex_fasta_restriction_ecori.json --output /tmp/gentle_pgex_restriction_expert
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json --output /tmp/gentle_pgex_restriction_expert_svg
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_tp53_isoform_architecture_online.json --output /tmp/gentle_tp53_isoform_workflow
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_inspect_feature_expert_tp53_isoform.json --output /tmp/gentle_tp53_isoform_text
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_feature_expert_tp53_isoform_svg.json --output /tmp/gentle_tp53_isoform_expert
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_tp53_splicing_expert_svg.json --output /tmp/gentle_tp53_splicing_expert
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_inspect_feature_expert_tp53_splicing.json --output /tmp/gentle_tp53_splicing_text
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_primers_preflight_auto.json --output /tmp/gentle_primers_preflight
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_seed_primers_tp53_splicing.json --output /tmp/gentle_tp53_primer_seed
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_seed_qpcr_tp53_feature.json --output /tmp/gentle_tp53_qpcr_feature_seed
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_seed_qpcr_tp53_splicing.json --output /tmp/gentle_tp53_qpcr_seed
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_seed_qpcr_tp53_splicing_specific_junction.json --output /tmp/gentle_tp53_qpcr_specific_seed
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_design_qpcr_taqman_tp53_operation.json --output /tmp/gentle_tp53_taqman_design
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_cdna_qpcr_taqman_test_demo_direct.json --output /tmp/gentle_cdna_taqman_test
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_cdna_pcr_qpcr_product_gel_nonspecific_offline.json --output /tmp/gentle_cdna_product_gel_workflow
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_cdna_pcr_products_gel_demo_direct.json --output /tmp/gentle_cdna_pcr_products_gel
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_cdna_qpcr_taqman_products_gel_demo_direct.json --output /tmp/gentle_cdna_qpcr_products_gel
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_primer_reports_list.json --output /tmp/gentle_primer_reports
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_primer_report_show_demo.json --output /tmp/gentle_primer_report_show
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_primer_report_export_demo.json --output /tmp/gentle_primer_report_export
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_qpcr_reports_list.json --output /tmp/gentle_qpcr_reports
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_qpcr_report_show_demo.json --output /tmp/gentle_qpcr_report_show
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_qpcr_report_export_demo.json --output /tmp/gentle_qpcr_report_export
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_p53_family_query_anchor_dotplot.json --output /tmp/gentle_p53_family_anchor
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_tp73_variant1_trypsin_digest_gel.json --output /tmp/gentle_tp73_trypsin_digest_gel
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_protocol_cartoon_gibson_svg.json --output /tmp/gentle_gibson_graphics
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_protocol_cartoon_qpcr_svg.json --output /tmp/gentle_qpcr_graphics
```

Notes:

- `request_intents_runtime.json` returns the installed wrapper's
  `gentle.clawbio_skill_intents_runtime.v1` payload, including the
  `INTENTS.json` descriptor hash and route/request-mode summary; ClawBio should
  use it for explicit descriptor refresh or drift checks rather than as a
  startup dependency.
- examples carrying `state_path: ".gentle_state.json"` expect a project state
  file in the working directory
- `request_render_svg_pgex_fasta_circular.json` is a common follow-on graphics
  route after `request_workflow_file.json`, which loads `pgex_fasta` into that
  state
- `request_export_bed_pgex_fasta_tfbs_restriction.json` is a matching
  follow-on tabular route on that same `pgex_fasta` state; it annotates TFBS,
  exports TFBS rows, and appends selected restriction-site rows into one BED
  bundle artifact
- `request_render_svg_pgex_fasta_linear_tfbs.json` and
  `request_render_svg_pgex_fasta_linear_restriction.json` are matching
  follow-on DNA-window graphics routes on that same `pgex_fasta` state
- `request_tfbs_summary_pgex_fasta.json`,
  `request_inspect_feature_expert_pgex_fasta_tfbs.json`, and
  `request_render_feature_expert_pgex_fasta_tfbs_svg.json` are follow-on TFBS
  routes after `request_render_svg_pgex_fasta_linear_tfbs.json`
- `request_primers_preflight_auto.json`,
  `request_seed_primers_tp53_splicing.json`,
  `request_seed_qpcr_tp53_feature.json`,
  `request_seed_qpcr_tp53_splicing.json`,
  `request_seed_qpcr_tp53_splicing_specific_junction.json`,
  `request_design_qpcr_taqman_tp53_operation.json`, and
  `request_cdna_qpcr_taqman_test_demo_direct.json` are typed ClawBio routes
  over the shared `primers ...` shell family; they cover backend readiness,
  PCR/qPCR seeding, transcript-aware TaqMan seeding, direct qPCR design, and
  direct cDNA qPCR assay testing without requiring callers to hand-write shell
  strings. The direct cDNA assay examples also declare transcript-map SVG
  artifacts so ClawBio can attach a PNG view of where the assay is functional.
- `request_workflow_cdna_pcr_qpcr_product_gel_nonspecific_offline.json`
  demonstrates the opt-in product-materialization path for transcript-derived
  cDNA PCR/qPCR assay tests: detected products become GENtle sequence entries,
  are grouped into one vial/container, and are rendered as a pool-gel lane so
  non-specific products appear as multiple bands. The direct product-gel
  examples (`request_cdna_pcr_products_gel_demo_direct.json` and
  `request_cdna_qpcr_taqman_products_gel_demo_direct.json`) are follow-on
  routes for a state that already contains the demo locus.
- `request_primer_reports_list.json`, `request_primer_report_show_demo.json`,
  `request_primer_report_export_demo.json`, `request_qpcr_reports_list.json`,
  `request_qpcr_report_show_demo.json`, and
  `request_qpcr_report_export_demo.json` cover the matching saved-report
  lifecycle from ClawBio
- `request_inspect_feature_expert_pgex_fasta_restriction_ecori.json` and
  `request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json` are
  follow-on restriction expert routes after `request_workflow_file.json`
- `request_export_bed_grch38_tp53_gene_models.json` is a follow-on annotation
  export after `request_genomes_extract_gene_tp53.json`
- `request_genomes_blast_grch38_short.json` is a follow-on search route after
  `request_genomes_prepare_grch38.json`; it exercises the shared
  reference-genome BLAST path against the prepared GRCh38 catalog entry
- `request_render_feature_expert_tp53_isoform_svg.json` is a follow-on expert
  route after `request_workflow_tp53_isoform_architecture_online.json` (or an
  equivalent prior isoform-panel import)
- `request_workflow_tp53_splicing_expert_svg.json` replays one deterministic
  offline splicing-expert workflow from the bundled
  `docs/figures/tp53_ensembl116_panel_source.gb` source asset and collects the
  rendered expert SVG into the output bundle
- graphics contract for ClawBio/OpenClaw:
  - graphics requests still declare SVG engine outputs in `expected_artifacts[]`
  - the wrapper rasterizes the best-first SVG into one deterministic PNG bundle
    artifact at fixed scale `2.0`
  - `result.json.preferred_artifacts[]` now points at the PNG bundle paths for
    figures, including `generated/clawbio_storyboard.png` when one run yields
    multiple related graphics
  - `result.json.artifact_summary` summarizes the best-first artifact,
    displayable/collected artifact counts, and continuation actions for
    one-image-per-reply chat gateways
  - remaining SVG figures stay in the bundle as provenance/supporting artifacts
    and are exposed through request-first `continue_artifact` actions rather
    than immediate extra PNG files
  - browser/OpenClaw inline image rendering remains a later ClawBio-side
    attachment/UI task
- capabilities/interaction contract for ClawBio/OpenClaw:
  - `request_capabilities.json` still runs the normal `gentle_cli capabilities`
    route
  - the wrapper now also performs a best-effort shared `ui intents` probe and,
    when available, lifts that payload into `result.json.ui_intent_catalog`
  - those capabilities replies may now also carry `suggested_actions[]` with
    `kind = ui_intent`, executable `ui open TARGET` shell lines, and attached
    `ui_intent` metadata so operators can be handed into current GUI surfaces
    without ClawBio hard-coding the target list

`gentle_local_checkout_cli.sh` defaults these paths when unset:

- `CARGO_TARGET_DIR=$GENTLE_REPO_ROOT/target`
- `GENTLE_REFERENCE_CACHE_DIR=$GENTLE_REPO_ROOT/data/genomes`
- `GENTLE_HELPER_CACHE_DIR=$GENTLE_REPO_ROOT/data/helper_genomes`

That keeps both builds and prepared caches inside the local GENtle checkout, so
first-run ClawBio preparation does not depend on any external/shared worktree
target layout.

Relative `workflow_path` values inside the wrapper resolve in this order:

1. current working directory
2. `GENTLE_REPO_ROOT` when set
3. the local GENtle repo containing the scaffold, when discoverable

When a resolved workflow lives inside a GENtle repo, the wrapper also executes
from that repo root so repo-relative assets referenced by the workflow keep
working after the scaffold is copied into a separate ClawBio checkout.

Quick scaffold usage (inside the scaffold directory):

```bash
python gentle_cloning.py --demo --output /tmp/gentle_clawbio_demo
python gentle_cloning.py --skill-info --output /tmp/gentle_clawbio_skill_info
python gentle_cloning.py --input examples/request_capabilities.json --output /tmp/gentle_clawbio_caps
```

When copied into a ClawBio checkout (`skills/gentle-cloning/`):

```bash
python clawbio.py run gentle-cloning --demo
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_file.json --output /tmp/gentle_clawbio_run
```

Command resolution order used by the wrapper:

1. `--gentle-cli` explicit command string
2. `GENTLE_CLI_CMD` environment variable
3. `gentle_cli` on `PATH`
4. repository-local fallback: `cargo run --quiet --bin gentle_cli --`

Typical Apptainer route for ClawBio/OpenClaw:

```bash
apptainer pull gentle.sif docker://ghcr.io/smoe/gentle_rs:cli
export GENTLE_CLI_CMD='skills/gentle-cloning/gentle_apptainer_cli.sh /absolute/path/to/gentle.sif'
```

Included first-run bootstrap request examples:

- `request_runtime_version.json`
- `request_version_installed.json`
- `request_genomes_list_human.json`
- `request_helpers_list_gst.json`
- `request_hosts_list_deor.json`
- `request_genomes_ensembl_available_human.json`
- `request_genomes_install_ensembl_mouse.json`
- `request_shell_state_summary.json`
- `request_services_telegram_guide.json`
- `request_services_telegram_guide_{readiness,gene_context,tfbs,inline_dna,cloning,isoforms,follow_up}.json`
- `request_genomes_status_grch38.json`
- `request_genomes_prepare_grch38.json`
- `request_helpers_status_puc19.json`
- `request_helpers_prepare_puc19.json`

Included follow-on analysis/planning/graphics examples:

- `request_genbank_fetch_pbr322.json`
- `request_dbsnp_fetch_rs9923231.json`
- `request_genomes_extract_gene_tp53.json`
- `request_export_bed_grch38_tp53_gene_models.json`
  - follow-on route after `request_genomes_extract_gene_tp53.json`
  - exports the extracted TP53 gene/mRNA/exon/CDS table to one BED6+4 artifact
- `request_genomes_blast_grch38_short.json`
  - follow-on route after `request_genomes_prepare_grch38.json`
  - exercises the shared `genomes blast ...` route against the prepared GRCh38
    Ensembl 116 reference catalog entry
- `request_helpers_blast_puc19_short.json`
- `request_render_svg_pgex_fasta_circular.json`
  - expects a state containing `pgex_fasta`, for example after
    `request_workflow_file.json`
- `request_export_bed_pgex_fasta_tfbs_restriction.json`
  - same `pgex_fasta` follow-on route, but exports TFBS/JASPAR hits plus
    selected restriction-site rows into one BED artifact
- `request_render_svg_pgex_fasta_linear_tfbs.json`
  - same `pgex_fasta` follow-on route, but with explicit JASPAR/TFBS display
  filtering before linear SVG export
- `request_tfbs_summary_pgex_fasta.json`
  - same `pgex_fasta` follow-on route, but emits grouped TFBS summary text
- `request_inspect_feature_expert_pgex_fasta_tfbs.json`
  - same `pgex_fasta` follow-on route, but opens one generated TFBS feature in
    textual expert form
- `request_render_feature_expert_pgex_fasta_tfbs_svg.json`
  - same `pgex_fasta` follow-on route, but renders one generated TFBS feature
    to expert SVG
- `request_render_svg_pgex_fasta_linear_restriction.json`
  - same `pgex_fasta` follow-on route, but with explicit restriction display
    settings before linear SVG export
- `request_inspect_feature_expert_pgex_fasta_restriction_ecori.json`
  - same `pgex_fasta` follow-on route, but inspects the EcoRI cleavage context
    in textual expert form
- `request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json`
  - same `pgex_fasta` follow-on route, but renders the EcoRI cleavage context
    to expert SVG
- `request_workflow_vkorc1_planning.json`
  - strongest graphical variant-follow-up answer in the current scaffold
  - the wrapper now promotes one best-first
    `generated/clawbio_storyboard.png` messenger artifact while keeping the
    underlying SVG figures as provenance/supporting outputs
- `request_workflow_tp53_isoform_architecture_online.json`
  - runs the canonical TP53 isoform workflow example and collects the rendered
    architecture SVG into the ClawBio bundle
- `request_inspect_feature_expert_tp53_isoform.json`
  - follow-on text companion after `request_workflow_tp53_isoform_architecture_online.json`
- `request_render_feature_expert_tp53_isoform_svg.json`
  - renders the same TP53 isoform architecture through the shared
    `render-feature-expert-svg ... isoform ...` expert route
- `request_workflow_tp53_splicing_expert_svg.json`
  - replays a deterministic offline splicing-expert workflow from the bundled
    TP53 Ensembl 116 panel-source GenBank asset and collects the rendered SVG
- `request_inspect_feature_expert_tp53_splicing.json`
  - follow-on text companion after `request_workflow_tp53_splicing_expert_svg.json`
- `request_workflow_p53_family_query_anchor_dotplot.json`
  - replays the anchored p53-family comparison with TP73 on the shared
    reference axis and TP63 plus TP53 aligned by the conserved motif
    `CATGTGTAACAG`
- `request_workflow_tp73_variant1_trypsin_digest_gel.json`
  - derives TP73 transcript variant 1 protein from the bundled offline TP73
    GenBank asset, renders a trypsin peptide digest gel SVG, and lets the
    ClawBio wrapper expose the PNG raster as the chat-facing preferred artifact
- `request_protocol_cartoon_gibson_svg.json`
  - declares `expected_artifacts[]` so the generated SVG is copied into the
    ClawBio output bundle under `generated/...` and rasterized into the
    PNG-first messenger contract
- shipped BED-export request examples now cover both common follow-on
  surfaces:
  - shell mode:
    `request_export_bed_grch38_tp53_gene_models.json`
  - workflow/op mode:
    `request_export_bed_pgex_fasta_tfbs_restriction.json`
  - both ride the shared routes directly:
    - `features export-bed SEQ_ID OUTPUT.bed [--coordinate-mode auto|local|genomic] [--include-restriction-sites] [--restriction-enzyme NAME ...] [feature-query filters ...]`
    - `ExportFeaturesBed { query, path, coordinate_mode, include_restriction_sites, restriction_enzymes[] }`
    - `InspectSequenceContextView { seq_id, mode?, viewport_start_0based?, viewport_end_0based_exclusive?, include_visible_classes[], coordinate_mode?, limit? }`
      for one compact chat-facing summary of the currently interesting genomic
      environment before attaching a larger SVG/BED bundle
    - `ExportSequenceContextBundle { seq_id, mode?, viewport_start_0based?, viewport_end_0based_exclusive?, coordinate_mode?, include_feature_bed?, include_text_summary?, include_restriction_sites?, restriction_enzymes[], output_dir }`
      - the resulting bundle manifest now carries ranked artifact metadata plus
        one explicit `best_first_artifact_*` selection so chat/report layers can
        choose a first figure deterministically
      for one deterministic DNA-window export directory containing SVG +
      summary JSON/text + optional BED companion artifacts

ClawBio troubleshooting note:

- If `python clawbio.py run gentle-cloning ...` reports `Unknown skill 'gentle-cloning'`,
  the copied scaffold is newer than the runtime `clawbio.py` registry on that
  machine.
- Regenerating `skills/catalog.json` alone is not enough when `clawbio.py`
  still lacks the `gentle-cloning` entry.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Agent Bridge](agents.md) | [Next: JavaScript Shell](javascript.md)
<!-- manual-page-nav:end -->
