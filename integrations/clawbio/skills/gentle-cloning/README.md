# gentle-cloning skill scaffold

This folder is a ClawBio/OpenClaw-ready skill scaffold for GENtle.

## Files

- `SKILL.md`: skill metadata + routing/instructions
- `gentle_cloning.py`: wrapper executable
- `gentle_local_checkout_cli.sh`: local-checkout launcher for `cargo run --bin gentle_cli --`
- `gentle_apptainer_cli.sh`: Apptainer/Singularity launcher for `gentle_cli`
- `catalog_entry.json`: ready-to-paste object for ClawBio `skills/catalog.json`
- `examples/*.json`: request payload examples, including bootstrap, extract/BLAST, planning, and graphics flows
- `tests/test_gentle_cloning.py`: minimal wrapper tests

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
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_status_grch38.json --output /tmp/gentle_status_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_prepare_grch38.json --output /tmp/gentle_prepare_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_prepare_puc19.json --output /tmp/gentle_prepare_puc19
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genbank_fetch_pbr322.json --output /tmp/gentle_fetch_pbr322
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_dbsnp_fetch_rs9923231.json --output /tmp/gentle_fetch_rs9923231
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_extract_gene_tp53.json --output /tmp/gentle_extract_tp53
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_blast_puc19_short.json --output /tmp/gentle_puc19_blast
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_vkorc1_planning.json --output /tmp/gentle_vkorc1_planning
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_svg_pgex_fasta_circular.json --output /tmp/gentle_pgex_map
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_svg_pgex_fasta_linear_tfbs.json --output /tmp/gentle_pgex_tfbs_map
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_svg_pgex_fasta_linear_restriction.json --output /tmp/gentle_pgex_restriction_map
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_tp53_isoform_architecture_online.json --output /tmp/gentle_tp53_isoform_workflow
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_render_feature_expert_tp53_isoform_svg.json --output /tmp/gentle_tp53_isoform_expert
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_protocol_cartoon_gibson_svg.json --output /tmp/gentle_gibson_graphics
```

Notes:

- examples carrying `state_path: ".gentle_state.json"` expect a project state
  file in the working directory
- `request_render_svg_pgex_fasta_circular.json` is a common follow-on graphics
  route after `request_workflow_file.json`, which loads `pgex_fasta` into that
  state
- `request_render_svg_pgex_fasta_linear_tfbs.json` and
  `request_render_svg_pgex_fasta_linear_restriction.json` are matching
  follow-on DNA-window graphics routes on that same `pgex_fasta` state
- `request_render_feature_expert_tp53_isoform_svg.json` is a follow-on expert
  route after `request_workflow_tp53_isoform_architecture_online.json` (or an
  equivalent prior isoform-panel import)

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
- optional: `state_path`, `timeout_secs`

Mode-specific fields:

- `shell`: `shell_line`
- `op`: `operation`
- `workflow`: `workflow` or `workflow_path`
- `raw`: `raw_args[]`

Relative `workflow_path` values are resolved in this order:

1. current working directory
2. `GENTLE_REPO_ROOT` when set
3. the local GENtle repo containing the copied scaffold, when discoverable

Included first-run bootstrap requests:

- `examples/request_genomes_list_human.json`
- `examples/request_helpers_list_gst.json`
- `examples/request_hosts_list_deor.json`
- `examples/request_genomes_ensembl_available_human.json`
- `examples/request_genomes_install_ensembl_mouse.json`
- `examples/request_shell_state_summary.json`
- `examples/request_genomes_status_grch38.json`
- `examples/request_genomes_prepare_grch38.json`
- `examples/request_helpers_status_puc19.json`
- `examples/request_helpers_prepare_puc19.json`

Included follow-on analysis/planning/graphics requests:

- `examples/request_genbank_fetch_pbr322.json`
- `examples/request_dbsnp_fetch_rs9923231.json`
- `examples/request_genomes_extract_gene_tp53.json`
- `examples/request_helpers_blast_puc19_short.json`
- `examples/request_render_svg_pgex_fasta_circular.json`
  - expects a state containing `pgex_fasta`, for example after running
    `examples/request_workflow_file.json`
- `examples/request_render_svg_pgex_fasta_linear_tfbs.json`
  - same `pgex_fasta` follow-on route, but with explicit JASPAR/TFBS display
    filtering before linear SVG export
- `examples/request_render_svg_pgex_fasta_linear_restriction.json`
  - same `pgex_fasta` follow-on route, but with explicit restriction display
    settings before linear SVG export
- `examples/request_workflow_vkorc1_planning.json`
- `examples/request_workflow_tp53_isoform_architecture_online.json`
  - replays the canonical TP53 isoform workflow example and collects the
    rendered architecture SVG into the ClawBio bundle
- `examples/request_render_feature_expert_tp53_isoform_svg.json`
  - renders the same TP53 isoform architecture through the shared
    `render-feature-expert-svg ... isoform ...` expert route
- `examples/request_protocol_cartoon_gibson_svg.json`
  - uses `expected_artifacts[]` so the generated SVG is copied into the
    wrapper output bundle under `generated/...`
- No canned BED-export request is shipped yet, but ClawBio can already call
  the shared routes directly:
  - shell mode:
    `features export-bed SEQ_ID OUTPUT.bed [--coordinate-mode auto|local|genomic] [--include-restriction-sites] [--restriction-enzyme NAME ...] [feature-query filters ...]`
  - raw/op mode:
    `ExportFeaturesBed { query, path, coordinate_mode, include_restriction_sites, restriction_enzymes[] }`
  - this covers genome annotation, TFBS/JASPAR features, and optional
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
