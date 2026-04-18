# gentle-cloning skill scaffold

This folder is a ClawBio/OpenClaw-ready skill scaffold for GENtle.

It is intended to position GENtle as a deterministic bridge from biological
observations, including patient/cohort-derived leads, to sequence-grounded
mechanistic follow-up and wet-lab planning. It also makes clear that prepared
Ensembl/reference assets and BLAST-capable indices are reusable local
infrastructure, not GENtle-only byproducts.

## Logical capability split

The runtime alias is still `gentle-cloning`, but the intended ClawBio-facing
surface is now six explicit sub-capabilities:

- genomic context
- TFBS analysis
- restriction analysis
- splicing expert
- isoform architecture
- variant follow-up

This is a documentation/routing split, not six separate executables. The goal
is to make it obvious that GENtle already exposes dedicated shared command
surfaces for those tasks instead of one vague "cloning" box.

## Files

- `SKILL.md`: skill metadata + routing/instructions
- `gentle_cloning.py`: wrapper executable
- `gentle_local_checkout_cli.sh`: local-checkout launcher for `cargo run --bin gentle_cli --`
- `gentle_apptainer_cli.sh`: Apptainer/Singularity launcher for `gentle_cli`
- `catalog_entry.json`: ready-to-paste object for ClawBio `skills/catalog.json`
- `examples/*.json`: request payload examples, including bootstrap, extract/BLAST, planning, and graphics flows
- `tests/test_gentle_cloning.py`: minimal wrapper tests

## Positioning for OpenClaw

When OpenClaw answers broad questions such as "How does GENtle help me?", the
intended framing is:

- GENtle helps move from an observation to a reproducible follow-up.
- For patient/cohort signals, that means:
  `statistical observation -> mechanistic hypothesis -> sequence/context analysis -> wet-lab validation plan`.
- GENtle does not prove causality by itself.
- GENtle can bootstrap reusable local Ensembl/reference assets, including
  BLAST-capable indices, that may also be useful to external bioinformatics
  tools. Its added value is deterministic preparation, cataloging, provenance,
  and reuse in the same workflow.

### Ensembl availability answers

For user questions like "Can GENtle get data from Ensembl?", the intended
ClawBio/OpenClaw answer pattern is:

- do not stop at "no direct remote database access"
- immediately say that GENtle can prepare and reuse a local
  Ensembl-backed reference copy
- inspect the current local state first with:
  - `services status`
  - `genomes status "Human GRCh38 Ensembl 116"`
  - `resources status`
  - `genomes ensembl-available --collection vertebrates --filter human`
- if the reference is not prepared yet, say so explicitly and offer the local
  prepare path rather than presenting the request as a dead end

This keeps the answer truthful while still surfacing the practical path
forward.

### Investigating stale Ensembl answers

If ClawBio still replies only with "GENtle cannot access remote databases
directly" after you updated the skill files, the most likely problem is that
the chat layer did not actually invoke the refreshed `gentle-cloning` skill.

The shortest investigation path is:

1. verify the copied remote skill bundle really contains:
   - `examples/request_services_status.json`
   - the `Ensembl availability answers` section in `SKILL.md`
2. run the skill directly from the ClawBio checkout:
   - `python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_services_status.json --output /tmp/gentle_services_status`
   - `python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_status_grch38.json --output /tmp/gentle_status_grch38`
   - the wrapper now resolves relative `--input` paths not only from the
     current cwd, but also from the copied skill layout itself, so these
     example paths keep working even when ClawBio launches the skill from a
     different working directory
3. inspect whether those direct runs succeed and produce `report.md` /
   `result.json`
   - failed runs now also surface the failing command, execution cwd, exit
     code, and a short stderr/stdout preview directly in `result.json.error`
     plus a dedicated `result.json.failure_summary` block
4. if the direct runs succeed but the chat answer is still the old one, the
   issue is not GENtle capability; it is ClawBio routing/prompting/caching
5. in that case, restart the ClawBio process that owns the chat session and
   re-test with a wording that strongly invites the skill, for example:
   - `Please use the GENtle cloning skill to tell me whether Ensembl-backed human reference data is available or can be prepared locally.`

If a direct run now fails with a message like `Unknown shell command
'services'`, treat that as a GENtle version mismatch first:

- the copied skill scaffold is newer than the installed `gentle_cli` binary on
  PATH
- update the installed binary, or point `GENTLE_CLI_CMD` at
  `gentle_local_checkout_cli.sh` for an updated checkout
- then re-run the same example request before debugging anything deeper

This investigation separates three cases cleanly:

- GENtle skill not invoked
- GENtle skill invoked but status command failed
- GENtle skill invoked and worked, but ClawBio summarized it poorly

## Recommended first-time route: local GENtle checkout

This is the intended newcomer path for ClawBio/BioClaw users who want to keep
GENtle editable in place instead of starting from a container.

```bash
cd /home/clawbio
git clone https://github.com/smoe/gentle_rs.git GENtle

export GENTLE_REPO_ROOT=/home/clawbio/GENtle
export GENTLE_CLI_CMD=/home/clawbio/ClawBio/skills/gentle-cloning/gentle_local_checkout_cli.sh

cd /home/clawbio/ClawBio
python clawbio.py run gentle-cloning --demo
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_list_human.json --output /tmp/gentle_list_human
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_list_gst.json --output /tmp/gentle_list_helpers
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_hosts_list_deor.json --output /tmp/gentle_list_hosts
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_ensembl_available_human.json --output /tmp/gentle_ensembl_human
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_install_ensembl_mouse.json --output /tmp/gentle_install_mouse
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_shell_state_summary.json --output /tmp/gentle_state_summary
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_services_status.json --output /tmp/gentle_services_status
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_status_grch38.json --output /tmp/gentle_status_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_resources_status.json --output /tmp/gentle_resources_status
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_prepare_grch38.json --output /tmp/gentle_prepare_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_blast_grch38_short.json --output /tmp/gentle_grch38_blast
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_prepare_puc19.json --output /tmp/gentle_prepare_puc19
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genbank_fetch_pbr322.json --output /tmp/gentle_fetch_pbr322
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_dbsnp_fetch_rs9923231.json --output /tmp/gentle_fetch_rs9923231
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_inspect_sequence_context_rs9923231_vkorc1.json --output /tmp/gentle_rs9923231_context_summary
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_export_sequence_context_bundle_rs9923231_vkorc1.json --output /tmp/gentle_rs9923231_context_bundle
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_svg_rs9923231_vkorc1_linear.json --output /tmp/gentle_rs9923231_context_svg
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_vkorc1_context_svg_auto_prepare.json --output /tmp/gentle_rs9923231_context_demo
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_export_bed_rs9923231_vkorc1_context_features.json --output /tmp/gentle_rs9923231_context_bed
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
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_p53_family_query_anchor_dotplot.json --output /tmp/gentle_p53_family_anchor
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_protocol_cartoon_gibson_svg.json --output /tmp/gentle_gibson_graphics
```

Notes:

- examples carrying `state_path: ".gentle_state.json"` expect a project state
  file in the working directory
- `request_render_svg_pgex_fasta_circular.json` is a common follow-on graphics
  route after `request_workflow_file.json`, which loads `pgex_fasta` into that
  state
- `request_render_svg_rs9923231_vkorc1_linear.json` is a matching genomic-context
  follow-on route after `request_dbsnp_fetch_rs9923231.json`; it renders the
  fetched VKORC1 / rs9923231 locus as a linear DNA-window SVG
- `request_workflow_vkorc1_context_svg_auto_prepare.json` is the lowest-hanging
  graphical demo route:
  it first ensures `Human GRCh38 Ensembl 116` is prepared locally, then fetches
  `rs9923231` and exports one linear genomic-context SVG into the wrapper
  bundle
- `request_inspect_sequence_context_rs9923231_vkorc1.json` is the compact
  summary follow-on after `request_dbsnp_fetch_rs9923231.json`; it returns the
  shared `InspectSequenceContextView` report so chat layers can answer with one
  concise viewport/context summary before attaching larger figures
  - the wrapper now promotes that report into `result.json.stdout_json` and
    `result.json.chat_summary_lines[]` so OpenClaw/ClawBio can relay the
    summary directly
- `request_export_sequence_context_bundle_rs9923231_vkorc1.json` is the
  matching one-shot bundle route after `request_dbsnp_fetch_rs9923231.json`;
  it writes the SVG, summary JSON/text, BED companion, and bundle manifest into
  one deterministic artifact directory for the same viewport
  - the wrapper also promotes the nested compact summary into
    `result.json.chat_summary_lines[]` when the exported bundle contains the
    shared `sequence_context_view` payload
  - and it now promotes any bundle-owned ranked artifact metadata into
    `result.json.preferred_artifacts[]`, so ClawBio can choose the best first
    figure without ad hoc path guessing
- `request_export_bed_rs9923231_vkorc1_context_features.json` is the matching
  coordinate export after `request_dbsnp_fetch_rs9923231.json`; it writes the
  fetched locus' gene/mRNA/variation rows with genomic coordinates into one BED
  artifact
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
  routes after `request_render_svg_pgex_fasta_linear_tfbs.json`; they expose
  grouped TFBS text summary, TFBS expert text, and TFBS expert SVG from the
  same annotated `pgex_fasta` state
- `request_inspect_feature_expert_pgex_fasta_restriction_ecori.json` and
  `request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json` are
  follow-on restriction expert routes after `request_workflow_file.json`; they
  inspect and render the EcoRI cleavage context on the loaded pGEX sequence
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

`gentle_local_checkout_cli.sh` defaults these paths when they are unset:

- `CARGO_TARGET_DIR=$GENTLE_REPO_ROOT/target`
- `GENTLE_REFERENCE_CACHE_DIR=$GENTLE_REPO_ROOT/data/genomes`
- `GENTLE_HELPER_CACHE_DIR=$GENTLE_REPO_ROOT/data/helper_genomes`

That keeps builds and prepared caches inside the local GENtle checkout, which
avoids the shared-target write issues that can show up in more locked-down
ClawBio deployments.

## Container smoke test

```bash
export GENTLE_CLI_CMD='docker run --rm -i -v "$PWD":/work -w /work ghcr.io/smoe/gentle_rs:cli'
python gentle_cloning.py --demo --output /tmp/gentle_clawbio_demo
```

Apptainer / Singularity alternative:

```bash
apptainer pull gentle.sif docker://ghcr.io/smoe/gentle_rs:cli
export GENTLE_CLI_CMD="$PWD/gentle_apptainer_cli.sh $PWD/gentle.sif"
python gentle_cloning.py --demo --output /tmp/gentle_clawbio_demo
```

This launcher uses `apptainer exec ... gentle_cli ...`. That is still the
intended ClawBio/OpenClaw route, even though the headless `:cli` image also
supports quick smoke tests such as `apptainer run gentle.sif capabilities`.

Alternative runtimes:

- local GENtle checkout via `gentle_local_checkout_cli.sh`
- local `gentle_cli` on `PATH`
- `--gentle-cli "<command>"`
- repository-local `cargo run --quiet --bin gentle_cli --`

## Request schema

- `schema`: `gentle.clawbio_skill_request.v1`
- `mode`: `capabilities|state-summary|shell|op|workflow|raw`
- optional: `state_path`, `timeout_secs`, `ensure_reference_prepared`

Mode-specific fields:

- `shell`: `shell_line`
- `op`: `operation`
- `workflow`: `workflow` or `workflow_path`
- `raw`: `raw_args[]`

`ensure_reference_prepared` is an opt-in wrapper preflight:

- required: `genome_id`
- optional: `catalog_path`, `cache_dir`, `status_timeout_secs`,
  `prepare_timeout_secs`

When present, the wrapper runs `genomes status ...` first and only runs
`genomes prepare ...` when the requested reference is not yet prepared. The
before/after status payloads and exact preflight commands are then written into
`report.md`, `result.json`, and `reproducibility/commands.sh`.

Relative `workflow_path` values are resolved in this order:

1. current working directory
2. `GENTLE_REPO_ROOT` when set
3. the local GENtle repo containing the copied scaffold, when discoverable

When a resolved workflow lives inside a GENtle repo, the wrapper also executes
from that repo root so repo-relative assets referenced by the workflow keep
working after the scaffold is copied into a separate ClawBio checkout.

Included first-run bootstrap requests:

- `examples/request_genomes_list_human.json`
- `examples/request_helpers_list_gst.json`
- `examples/request_hosts_list_deor.json`
- `examples/request_genomes_ensembl_available_human.json`
- `examples/request_genomes_install_ensembl_mouse.json`
- `examples/request_shell_state_summary.json`
- `examples/request_services_status.json`
- `examples/request_genomes_status_grch38.json`
- `examples/request_resources_status.json`
- `examples/request_genomes_prepare_grch38.json`
- `examples/request_helpers_status_puc19.json`
- `examples/request_helpers_prepare_puc19.json`

Included follow-on analysis/planning/graphics requests:

- `examples/request_genbank_fetch_pbr322.json`
- `examples/request_dbsnp_fetch_rs9923231.json`
- `examples/request_inspect_sequence_context_rs9923231_vkorc1.json`
  - follow-on route after `examples/request_dbsnp_fetch_rs9923231.json`
  - returns one compact viewport-aware context summary for chat/report layers
  - prefer `result.json.chat_summary_lines[]` as the first reply, then attach
    SVG/BED follow-ons as needed
- `examples/request_export_sequence_context_bundle_rs9923231_vkorc1.json`
  - follow-on route after `examples/request_dbsnp_fetch_rs9923231.json`
  - exports one deterministic directory containing SVG + summary JSON/text +
    BED + manifest for the same genomic-context window
- `examples/request_render_svg_rs9923231_vkorc1_linear.json`
  - follow-on route after `examples/request_dbsnp_fetch_rs9923231.json`
  - renders the fetched VKORC1 / rs9923231 locus as a linear genomic-context SVG
- `examples/request_workflow_vkorc1_context_svg_auto_prepare.json`
  - lowest-hanging graphical demo for remote ClawBio/OpenClaw installs
  - auto-checks/prepares `Human GRCh38 Ensembl 116`, fetches `rs9923231`, and
    exports a compact linear genomic-context SVG into the wrapper bundle
- `examples/request_resources_status.json`
  - reports which integrated external resource snapshots are active right now
    (`REBASE`, `JASPAR`) and records `ATtRACT` explicitly as
    not-yet-integrated, including its current published ZIP download URL
- `examples/request_services_status.json`
  - reports one combined readiness view across canonical references, helper
    backbones, and active external resource snapshots so chat/report layers can
    answer "what can this GENtle instance work with right now?" from one
    deterministic artifact
  - when a prepare/reindex run is in flight, the same report can also say
    whether GENtle is currently downloading, indexing, cancelled, or failed
    instead of only reporting static prepared/not-prepared state
- `examples/request_export_bed_rs9923231_vkorc1_context_features.json`
  - follow-on route after `examples/request_dbsnp_fetch_rs9923231.json`
  - exports the fetched locus' gene/mRNA/variation rows with genomic coordinates
- `examples/request_genomes_extract_gene_tp53.json`
- `examples/request_export_bed_grch38_tp53_gene_models.json`
  - follow-on route after `examples/request_genomes_extract_gene_tp53.json`
  - exports the extracted TP53 gene/mRNA/exon/CDS table to one BED6+4 artifact
- `examples/request_genomes_blast_grch38_short.json`
  - follow-on route after `examples/request_genomes_prepare_grch38.json`
  - exercises the shared `genomes blast ...` route against the prepared GRCh38
    Ensembl 116 reference catalog entry
- `examples/request_helpers_blast_puc19_short.json`
- `examples/request_render_svg_pgex_fasta_circular.json`
  - expects a state containing `pgex_fasta`, for example after running
    `examples/request_workflow_file.json`
- `examples/request_export_bed_pgex_fasta_tfbs_restriction.json`
  - same `pgex_fasta` follow-on route, but exports TFBS/JASPAR hits plus
    selected restriction-site rows into one BED artifact
- `examples/request_render_svg_pgex_fasta_linear_tfbs.json`
  - same `pgex_fasta` follow-on route, but with explicit JASPAR/TFBS display
  filtering before linear SVG export
- `examples/request_tfbs_summary_pgex_fasta.json`
  - same `pgex_fasta` follow-on route, but emits grouped TFBS summary text for
    a defined focus/context window
- `examples/request_inspect_feature_expert_pgex_fasta_tfbs.json`
  - same `pgex_fasta` follow-on route, but opens one generated TFBS feature in
    textual expert form
- `examples/request_render_feature_expert_pgex_fasta_tfbs_svg.json`
  - same `pgex_fasta` follow-on route, but renders one generated TFBS feature
    to expert SVG
- `examples/request_render_svg_pgex_fasta_linear_restriction.json`
  - same `pgex_fasta` follow-on route, but with explicit restriction display
    settings before linear SVG export
- `examples/request_inspect_feature_expert_pgex_fasta_restriction_ecori.json`
  - same `pgex_fasta` follow-on route, but inspects the EcoRI cleavage context
    in textual expert form
- `examples/request_render_feature_expert_pgex_fasta_restriction_ecori_svg.json`
  - same `pgex_fasta` follow-on route, but renders the EcoRI cleavage context
    to expert SVG
- `examples/request_workflow_vkorc1_planning.json`
- `examples/request_workflow_tp53_isoform_architecture_online.json`
  - replays the canonical TP53 isoform workflow example and collects the
    rendered architecture SVG into the ClawBio bundle
- `examples/request_inspect_feature_expert_tp53_isoform.json`
  - follow-on text companion after
    `examples/request_workflow_tp53_isoform_architecture_online.json`
- `examples/request_render_feature_expert_tp53_isoform_svg.json`
  - renders the same TP53 isoform architecture through the shared
    `render-feature-expert-svg ... isoform ...` expert route
- `examples/request_workflow_tp53_splicing_expert_svg.json`
  - replays a deterministic offline splicing-expert workflow from the bundled
  TP53 Ensembl 116 panel-source GenBank asset and collects the rendered SVG
- `examples/request_inspect_feature_expert_tp53_splicing.json`
  - follow-on text companion after
    `examples/request_workflow_tp53_splicing_expert_svg.json`
- `examples/request_workflow_p53_family_query_anchor_dotplot.json`
  - replays the cross-family anchored dotplot with TP73 as the shared
    reference axis and TP63 plus TP53 aligned by the conserved motif
    `CATGTGTAACAG`
- `examples/request_protocol_cartoon_gibson_svg.json`
  - uses `expected_artifacts[]` so the generated SVG is copied into the
  wrapper output bundle under `generated/...`
- shipped BED-export examples now cover both common follow-on surfaces:
  - shell mode:
    `request_export_bed_grch38_tp53_gene_models.json`
  - workflow/op mode:
    `request_export_bed_pgex_fasta_tfbs_restriction.json`
  - both ride the shared routes directly:
    - `features export-bed SEQ_ID OUTPUT.bed [--coordinate-mode auto|local|genomic] [--include-restriction-sites] [--restriction-enzyme NAME ...] [feature-query filters ...]`
    - `ExportFeaturesBed { query, path, coordinate_mode, include_restriction_sites, restriction_enzymes[] }`
  - together they cover genome annotation, TFBS/JASPAR rows, and optional
    deterministic REBASE restriction-site rows through one BED6+4 export
    contract

## Troubleshooting

- If `python clawbio.py run gentle-cloning ...` says `Unknown skill 'gentle-cloning'`, your
  ClawBio runtime registry is older than the scaffold.
- Check `python clawbio.py list | grep gentle-cloning`.
- Rebuilding `skills/catalog.json` alone is not enough when the runtime
  `clawbio.py` registry does not yet carry the `gentle-cloning` entry.
- Update the ClawBio checkout or add the current `gentle-cloning` runtime
  registration block before retrying.
