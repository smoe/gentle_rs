# GENtle CLI Manual

This page documents command-line entry points for GENtle.

## Overview

GENtle currently provides six binaries:

- `gentle`: graphical desktop app (GUI)
- `gentle_cli`: JSON operation/workflow CLI for automation and AI tools
- `gentle_js`: interactive JavaScript shell (optional build target; requires
  feature `js-interface`)
- `gentle_lua`: interactive Lua shell (optional build target; requires feature
  `lua-interface`)
- `gentle_examples_docs`: generates adapter snippets and tutorial artifacts from canonical protocol examples
- `gentle_mcp`: MCP stdio server (guarded mutating + UI-intent parity baseline;
  includes standardized capability discovery via `tools/list`,
  `capabilities`, and `help`)

In addition, the GUI includes an embedded `Shell` panel that uses the same
shared shell parser/executor as `gentle_cli shell`.

Structured command glossary:

- `docs/glossary.json` is the machine-readable command glossary used for
  per-command help rendering (`help ...`) and catalog export
  (`help --format json|markdown`).
- Shell-only dry-run catalog preview routes are available through the shared
  shell, for example:
  - `gentle_cli shell 'genomes preview-ensembl-specs --catalog assets/genomes.json'`
  - `gentle_cli shell 'helpers preview-ensembl-specs --catalog assets/helper_genomes.json'`

Structured workflow examples:

- canonical source files: `docs/examples/workflows/*.json`
- schema: `gentle.workflow_example.v1`
- includes GUI-first parity skeleton:
  - `docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json`
  - merged canonical tutorial:
    `docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`
- on-demand snippet generation (CLI/shell/JS/Lua):
  - `cargo run --bin gentle_examples_docs -- generate`
- validation only:
  - `cargo run --bin gentle_examples_docs -- --check`
- tutorial generation:
  - `cargo run --bin gentle_examples_docs -- tutorial-generate`
- tutorial drift check:
  - `cargo run --bin gentle_examples_docs -- tutorial-check`
- test gating:
  - `always` examples run in default tests
  - `online` examples run only with `GENTLE_TEST_ONLINE=1`
  - set `GENTLE_SKIP_REMOTE_TESTS=1` to force-skip all remote-resource tests
    (takes precedence over `GENTLE_TEST_ONLINE`)
- `skip` examples are syntax-checked only

Architecture invariant: all adapters/frontends above route cloning/business
behavior through the same shared engine.

Catalog path note:

- `--catalog PATH` accepts either one JSON file or one directory of `*.json`
  catalog fragments.
- Directory fragments are merged in sorted filename order; duplicate entry ids
  currently fail fast instead of silently overriding each other.
- When `--catalog` is omitted, GENtle now uses a deterministic discovery chain
  instead of one hard-coded file:
  - built-in `assets/*.json` plus optional built-in `assets/*.d/`
  - system overlays under `/etc/gentle/catalogs/`
  - user overlays under `$XDG_CONFIG_HOME/gentle/catalogs/` or
    `~/.config/gentle/catalogs/`
  - project overlays under `PROJECT_ROOT/.gentle/catalogs/`
- Built-in/system/project roots can be overridden in controlled environments
  via `GENTLE_ASSET_ROOT`, `GENTLE_SYSTEM_CONFIG_ROOT`, and
  `GENTLE_PROJECT_ROOT`.
- `genomes list` and `helpers list` now also accept `--filter TEXT` and search
  ids plus catalog metadata such as aliases, tags, search terms, and the new
  helper semantic/procurement fields.
- Helper catalogs may now carry richer metadata (`summary`, `aliases`, `tags`,
  `search_terms`, `helper_kind`, `host_system`, `procurement`, `semantics`)
  while still using the same preparation/indexing pipeline as before.

## ClawBio/OpenClaw skill scaffold

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

- genomic context
- TFBS analysis
- restriction analysis
- splicing expert
- isoform architecture
- variant follow-up

This is a ClawBio-facing routing/documentation split, not six separate
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
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_list_human.json --output /tmp/gentle_list_human
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_list_gst.json --output /tmp/gentle_list_helpers
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_hosts_list_deor.json --output /tmp/gentle_list_hosts
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_ensembl_available_human.json --output /tmp/gentle_ensembl_human
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_install_ensembl_mouse.json --output /tmp/gentle_install_mouse
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_shell_state_summary.json --output /tmp/gentle_state_summary
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
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_seed_qpcr_tp53_splicing.json --output /tmp/gentle_tp53_qpcr_seed
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_p53_family_query_anchor_dotplot.json --output /tmp/gentle_p53_family_anchor
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_protocol_cartoon_gibson_svg.json --output /tmp/gentle_gibson_graphics
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_protocol_cartoon_qpcr_svg.json --output /tmp/gentle_qpcr_graphics
```

Notes:

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
- `request_seed_qpcr_tp53_splicing.json` is a follow-on shell route on the
  TP53 splicing example state; it emits the non-mutating
  `gentle.qpcr_seed_request.v1` payload for the saved splicing group
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

- `request_genomes_list_human.json`
- `request_helpers_list_gst.json`
- `request_hosts_list_deor.json`
- `request_genomes_ensembl_available_human.json`
- `request_genomes_install_ensembl_mouse.json`
- `request_shell_state_summary.json`
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
- `request_protocol_cartoon_gibson_svg.json`
  - declares `expected_artifacts[]` so the generated SVG is copied into the
    ClawBio output bundle under `generated/...`
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

## Python module wrapper (`gentle_py`)

GENtle also ships a thin Python wrapper over `gentle_cli`:

- module path: `integrations/python/gentle_py/`
- quick docs: `integrations/python/README.md`

Quick usage:

```bash
python3 - <<'PY'
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path("integrations/python").resolve()))
from gentle_py import GentleClient

client = GentleClient(state_path=".gentle_state.json")
print(json.dumps(client.capabilities(), indent=2))
PY
```

Deterministic API surface:

- `capabilities()`
- `state_summary()`
- `op(operation)`
- `workflow(workflow|workflow_path=...)`
- `shell(line, expect_json=False)`
- `render_dotplot_svg(seq_id, dotplot_id, output_svg, ...)`

CLI resolution order:

1. constructor `cli_cmd`
2. `GENTLE_CLI_CMD`
3. `gentle_cli` on `PATH`
4. repository fallback `cargo run --quiet --bin gentle_cli --`

Resource update capability status:

- `gentle_cli`: supported (`resources sync-rebase`, `resources sync-jaspar`, `resources sync-jaspar-remote-metadata`, `resources summarize-jaspar`, `resources benchmark-jaspar`, `resources list-jaspar`, `resources inspect-jaspar`)
- `gentle_js`: supported (`sync_rebase`, `sync_jaspar`)
- `gentle_lua`: supported (`sync_rebase`, `sync_jaspar`)

Reference genome capability status:

- `gentle_cli`: supported via shared engine operations (`PrepareGenome`, `ExtractGenomeRegion`, `ExtractGenomeGene`) and shell-level `genomes/helpers blast`
- `gentle_js`: supported via dedicated helpers (`list_reference_genomes`, `list_reference_catalog_entries`, `list_helper_catalog_entries`, `list_host_profile_catalog_entries`, `list_ensembl_installable_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`
- `gentle_lua`: supported via dedicated helpers (`list_reference_genomes`, `list_reference_catalog_entries`, `list_helper_catalog_entries`, `list_host_profile_catalog_entries`, `list_ensembl_installable_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`

Construct-reasoning inspection capability status:

- `gentle_cli`: supported via shared shell/direct commands (`construct-reasoning list-graphs`, `construct-reasoning show-graph`, `construct-reasoning set-annotation-status`, `construct-reasoning write-annotation`, `construct-reasoning export-graph`)
- `gentle_js`: supported via dedicated shared-shell-backed helpers (`list_construct_reasoning_graphs`, `show_construct_reasoning_graph`, `set_construct_reasoning_annotation_status`, `write_back_construct_reasoning_annotation`)
- `gentle_lua`: supported via dedicated shared-shell-backed helpers (`list_construct_reasoning_graphs`, `show_construct_reasoning_graph`, `set_construct_reasoning_annotation_status`, `write_back_construct_reasoning_annotation`)
- `gentle_mcp`: supported via thin tool wrappers over the same shared shell contracts (`construct_reasoning_graphs`, `construct_reasoning_graph`, `construct_reasoning_set_annotation_status`, `construct_reasoning_write_annotation`)

Agent-assistant capability status:

- `gentle_cli`: supported via shared-shell command family (`agents list`, `agents ask`) and direct forwarding (`gentle_cli agents ...`)
- `gentle_js`: supported via helper wrappers (`list_agent_systems`, `ask_agent_system`) over shared shell execution
- `gentle_lua`: supported via helper wrappers (`list_agent_systems`, `ask_agent_system`) over shared shell execution
- GUI: supported via standalone `Agent Assistant` viewport using the same bridge and command executor

Candidate-set capability status:

- `gentle_cli`: supported as first-class `candidates ...` commands and shared-shell `candidates ...` commands, backed by shared engine operations (`GenerateCandidateSet`, `GenerateCandidateSetBetweenAnchors`, `DeleteCandidateSet`, `ScoreCandidateSetExpression`, `ScoreCandidateSetDistance`, `ScoreCandidateSetWeightedObjective`, `TopKCandidateSet`, `ParetoFrontierCandidateSet`, `FilterCandidateSet`, `CandidateSetOp`, `UpsertCandidateMacroTemplate`, `DeleteCandidateMacroTemplate`)
- `gentle_js`: supported via `apply_operation` with the same candidate-set operations
- `gentle_lua`: supported via `apply_operation` with the same candidate-set operations

Guide-design capability status:

- `gentle_cli`: supported as first-class `guides ...` commands and shared-shell `guides ...` commands, backed by shared engine operations (`UpsertGuideSet`, `DeleteGuideSet`, `FilterGuidesPractical`, `GenerateGuideOligos`, `ExportGuideOligos`, `ExportGuideProtocolText`)
- `gentle_js`: supported via `apply_operation` with the same guide-design operations
- `gentle_lua`: supported via `apply_operation` with the same guide-design operations

Primer-design report capability status:

- `gentle_cli`: supported via shared-shell `primers ...` commands and direct forwarding (`gentle_cli primers ...`), backed by `DesignPrimerPairs`, `DesignQpcrAssays`, and the post-design cloning handoff operation `PrepareRestrictionCloningPcrHandoff`, plus non-mutating ROI seed helpers (`primers seed-from-feature`, `primers seed-from-splicing`), a non-mutating restriction-cloning handoff request seeder (`primers seed-restriction-cloning-handoff`), vector suggestion helpers (`primers restriction-cloning-vector-suggestions`), and persisted report inspect/export helpers for primer, qPCR, and restriction-cloning handoff reports. `--progress` now also streams `progress primers ...` lines for shared-shell primer/qPCR design commands, not only raw `op` / `workflow` JSON execution.
- `gentle_js`: supported via `apply_operation` (`DesignPrimerPairs`, `DesignQpcrAssays`, `PrepareRestrictionCloningPcrHandoff`) plus shared-shell execution for report listing/show/export
- `gentle_lua`: supported via `apply_operation` (`DesignPrimerPairs`, `DesignQpcrAssays`, `PrepareRestrictionCloningPcrHandoff`) plus shared-shell execution for report listing/show/export

Dotplot/flexibility capability status:

- `gentle_cli`: supported via shared-shell/direct commands:
  - `dotplot compute|overlay-compute|list|show|render-svg`
  - `render-dotplot-svg`
  - `transcripts derive`
  - `flex compute|list|show`
  - `splicing-refs derive`
  - `align compute`
  backed by `ComputeDotplot`, `ComputeDotplotOverlay`, `DeriveTranscriptSequences`,
  `ComputeFlexibilityTrack`, `DeriveSplicingReferences`,
  `AlignSequences`, and `RenderDotplotSvg`.
  - `dotplot compute` supports self and pairwise modes via
    `--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement`
    with optional `--reference-seq`, `--ref-start`, and `--ref-end`.
  - `dotplot overlay-compute` accepts repeated `--query-spec JSON_OR_@FILE`
    payloads for shared-reference multi-query overlays; if `--reference-seq` is
    omitted, the owner sequence id is reused as the shared reference.
  - `dotplot render-svg` / `render-dotplot-svg` accept
    `--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor`
    so overlay exports can switch between normalized transcript length,
    left-aligned bp, right-aligned bp, and shared-exon-anchor layouts
    without recomputing the dotplot.
  - `transcripts derive` supports:
    - full-sequence derivation (all `mRNA`/`transcript` features)
    - selected feature derivation (`--feature-id`)
    - splicing-scope constrained derivation from one seed feature
      (`--scope ...` with exactly one `--feature-id`)
    - additive synthetic local `CDS` + translated-protein qualifiers when the
      admitted transcript has resolvable CDS context and translation-table
      metadata
- `gentle_js`: baseline support via `apply_operation` (`ComputeDotplot`,
  `DeriveTranscriptSequences`, `ComputeFlexibilityTrack`,
  `DeriveSplicingReferences`, `AlignSequences`) plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.
- `gentle_lua`: baseline support via `apply_operation` (`ComputeDotplot`,
  `DeriveTranscriptSequences`, `ComputeFlexibilityTrack`,
  `DeriveSplicingReferences`, `AlignSequences`) plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.
- `gentle_py`: baseline support via `op(...)` plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.

RNA-read interpretation capability status (Nanopore cDNA phase-1):

- `gentle_cli`: supported via shared-shell/direct commands:
  - `rna-reads interpret`
  - `rna-reads align-report`
  - `rna-reads list-reports`
  - `rna-reads show-report`
  - `rna-reads summarize-gene-support`
  - `rna-reads inspect-gene-support`
  - `rna-reads inspect-alignments`
  - `rna-reads export-report`
  - `rna-reads export-hits-fasta`
  - `rna-reads export-sample-sheet`
  - `rna-reads export-paths-tsv`
  - `rna-reads export-abundance-tsv`
  - `rna-reads export-score-density-svg`
  - `rna-reads export-alignments-tsv`
  - `rna-reads export-alignment-dotplot-svg`
  backed by `InterpretRnaReads`, `AlignRnaReadReport`,
  `ListRnaReadReports`, `ShowRnaReadReport`,
  `SummarizeRnaReadGeneSupport`, `InspectRnaReadGeneSupport`,
  `ExportRnaReadReport`,
  `ExportRnaReadHitsFasta`, `ExportRnaReadSampleSheet`,
  `ExportRnaReadExonPathsTsv`, `ExportRnaReadExonAbundanceTsv`, `ExportRnaReadScoreDensitySvg`,
  `ExportRnaReadAlignmentsTsv`, and `ExportRnaReadAlignmentDotplotSvg`.
  Input supports FASTA plus gzipped FASTA (`.fa/.fasta` and `.fa.gz/.fasta.gz`).
  Concatenated gzip members are accepted for gzipped FASTA input as well.
  Progress output includes periodic `progress rna-reads ...` lines during
  `apply_with_progress` runs.
  Phase split:
  - `interpret`: seed-filter pass (Nanopore phase-1)
    - `align-report`: retained-hit alignment pass (phase-2) that updates
      mapping fields, MSA-eligibility counters, exon-transition rows, and
      exon/junction abundance frequencies in the persisted report; retained hits
      are re-ranked by alignment-aware retention rank after alignment.
    - supports explicit row filtering via
      `--record-indices i,j,k` (0-based stored `record_index` values);
      when provided, this overrides `--selection`.
    - default adapter behavior now uses `selection=all` so rescued retained
      rows also receive round-2 similarity/coverage scores.
  - `summarize-gene-support`: non-mutating target-gene cohort summary over one
    saved aligned RNA-read report.
    - `--gene GENE` is required and repeatable; matches are case-insensitive
      against the same transcript group-label logic used by the splicing view.
    - optional `--record-indices i,j,k` restricts the base aligned cohort to an
      exact saved-report subset before target-gene filtering.
    - `--complete-rule near|strict|exact` controls the fragment-vs-complete
      cohort split; strict and exact counts are still reported explicitly in the
      payload.
    - output is JSON on stdout by default; `--output PATH` writes the exact same
      payload to disk.
    - payload tables include:
      - per-cohort exon support
      - ordered exon-pair support (including skipped pairs like `1->3`)
      - neighboring direct-transition support only (for example `1->2`)
  - `inspect-gene-support`: non-mutating row-level audit of the same
    target-gene cohort logic used by `summarize-gene-support`.
    - `--gene GENE` is required and repeatable, with the same
      case-insensitive group-label matching.
    - optional `--record-indices i,j,k` restricts the evaluation universe
      before cohort classification.
    - `--complete-rule near|strict|exact` controls whether accepted rows land
      in `accepted_fragment` vs `accepted_complete`.
    - `--cohort all|accepted|fragment|complete|rejected` filters the returned
      row list without changing the grouped top-level record-index arrays.
    - output is JSON on stdout by default; `--output PATH` writes the exact
      same payload to disk.
    - rows include resolved gene/transcript identity when available, status +
      machine-readable `status_reason`, full-length flags/class, mapped exon
      ordinals, ordered exon pairs, neighboring direct transitions, phase-2
      score/identity/coverage, and `passed_seed_filter`.
  - `inspect-alignments`: non-mutating ranked alignment inspection
    over persisted report hits, with optional structured subset controls:
    `--effect-filter`, `--sort`, `--search`, and `--record-indices`.
  - `export-alignments-tsv`: non-mutating ranked alignment-row TSV export for
    downstream table-based triage.
  - `export-alignment-dotplot-svg`: non-mutating dotplot-like scatter export
    (coverage vs identity, score-colored points) for aligned report hits.
- `gentle_js`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_lua`: baseline support via `apply_operation` for the same operation
  family.
- GUI workflow parity:
  - Splicing Expert `Nanopore cDNA interpretation` now exposes
    `Prepare Workflow Op` and `Copy Workflow JSON` so the exact
    `InterpretRnaReads` payload can be executed via regular workflow routes.
  - Example:
    - `gentle_cli workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`
    - or paste copied JSON directly:
      `gentle_cli workflow '{"run_id":"workflow_rna_reads_tp73_ncbi_cdna_srr32957124","ops":[{"InterpretRnaReads":{"seq_id":"tp73.ncbi","seed_feature_id":0,"profile":"nanopore_cdna_v1","input_path":"reads.fa.gz","input_format":"fasta","scope":"all_overlapping_both_strands","seed_filter":{"kmer_len":10,"seed_stride_bp":1,"min_seed_hit_fraction":0.30,"min_weighted_seed_hit_fraction":0.05,"min_unique_matched_kmers":12,"max_median_transcript_gap":4.0,"min_chain_consistency_fraction":0.40,"min_confirmed_exon_transitions":1,"min_transition_support_fraction":0.05,"cdna_poly_t_flip_enabled":true,"poly_t_prefix_min_bp":18},"align_config":{"band_width_bp":24,"min_identity_fraction":0.60,"max_secondary_mappings":0},"report_id":"cdna_reads"}}]}'`

Sequencing-confirmation capability status (called-read plus trace-aware engine baseline):

- `gentle_cli`: supported via shared-shell commands plus the current GUI
  specialist:
  - `seq-trace import PATH [--trace-id ID] [--seq-id ID]`
  - `seq-trace list [SEQ_ID]`
  - `seq-trace show TRACE_ID`
  - `seq-confirm run EXPECTED_SEQ_ID [--baseline BASELINE_SEQ_ID] --reads ID[,ID...] [--trace-id ID]... [--trace-ids ID[,ID...]]`
  - `seq-confirm list-reports [EXPECTED_SEQ_ID]`
  - `seq-confirm show-report REPORT_ID`
  - `seq-confirm export-report REPORT_ID OUTPUT.json`
  - `seq-confirm export-support-tsv REPORT_ID OUTPUT.tsv`
  - `seq-primer suggest EXPECTED_SEQ_ID [--primers ID[,ID...]] [--confirmation-report REPORT_ID] [--min-3prime-anneal-bp N] [--predicted-read-length-bp N]`
  raw-trace intake is backed by `ImportSequencingTrace`,
  `ListSequencingTraces`, and `ShowSequencingTrace`.
  backed by `ConfirmConstructReads`, `ListSequencingConfirmationReports`,
  `ShowSequencingConfirmationReport`, `ExportSequencingConfirmationReport`, and
  `ExportSequencingConfirmationSupportTsv`. Sequencing-primer overlays are
  backed by `SuggestSequencingPrimers`.
  Current scope:
  - expected construct plus already-loaded read sequences
  - imported ABI/AB1/SCF traces can now participate directly in
    `seq-confirm run` without first materializing a project sequence
  - default full-span confirmation when no explicit targets are supplied
  - explicit junction targets via `--junction LEFT_END_0BASED`
    and `--junction-flank N`
  - forward and reverse-complement read evaluation through the same shared
    engine report contract
  - raw ABI/AB1/SCF intake stores trace evidence separately from
    confirmation reports:
    - preserves file-supplied called bases
    - preserves called-base confidence arrays when available
    - preserves peak locations when available
    - preserves raw per-channel chromatogram curves and clip windows when
      available from the source trace
    - keeps one imported trace store that can be listed/shown without running
      confirmation
  - trace-aware confirmation keeps one shared report model:
    - sequence-backed and trace-backed evidence rows land in the same
      confirmation report
    - trace-backed rows expose their `trace_id` and evidence kind
    - target support/contradiction ids can now reflect imported trace ids
    - optional `--baseline` context lets confirmation classify one locus as:
      intended edit, reference reversion, unexpected difference, or
      low-confidence/insufficient evidence
    - `seq-confirm show-report` now exposes the persisted `baseline_seq_id`
      plus per-variant rows for expected edits inferred from the
      baseline-vs-expected diff
  - sequencing-primer suggestion now has two modes through the same
    `SuggestSequencingPrimers` operation:
    - overlay existing primer ids onto the expected construct
    - or, with only `--confirmation-report`, emit fresh primer proposals for
      unresolved loci that still lack a good existing hit
    - when existing primer hits do exist but sit outside the preferred
      sequencing window, the report now keeps both the best existing hit and
      the top fresh proposal
  - GUI parity is now available through `Patterns -> Sequencing Confirmation...`
    and command palette `Sequencing Confirmation`, backed by the same
    `ConfirmConstructReads` / report-store path:
    - optional baseline sequence input is available in the run form
    - the review pane now includes a variant list and chromatogram curve view
      for trace-backed loci
    - sequencing-primer overlays can now be suggested in the same specialist
      from already-loaded primer sequence ids
    - the overlay report can optionally annotate which saved-report targets and
      variant loci fall inside each predicted primer-derived read span
    - unresolved targets and variant loci now also get one shared
      recommendation row naming the best existing primer hit to clarify them
    - older traces without stored curve arrays stay usable for confirmation,
      but the GUI asks you to re-import them before chromatogram review
  Not yet included in this phase:
  - GUI-side raw trace import
  - chromatogram editing or base re-calling
  - full whole-trace browsing beyond the current variant-focused window
  - lineage/artifact projection for confirmation reports
- `gentle_js`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_lua`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_py`: baseline support via `op(...)` for the same operation family.

Feature-expert SVG parity status:

- `gentle_cli`: supported via direct commands and shared-shell command path
  (`inspect-feature-expert`, `render-feature-expert-svg`)
- `gentle_js`: supported via `render_feature_expert_svg(...)` wrapper over
  shared engine operation `RenderFeatureExpertSvg`
- `gentle_lua`: supported via `render_feature_expert_svg(...)` wrapper over
  shared engine operation `RenderFeatureExpertSvg`
- splicing expert SVG now includes:
  - explicit junction transition support labels/table
  - frequency-encoded transcript-vs-exon matrix cell coloring
  - predicted exon-to-exon transition matrix with frequency-encoded cell
    coloring
  - CDS flank phase coloring (`0/1/2`) on exon left/right edges when
    transcript `cds_ranges_1based` are available
  - exon-length modulo-3 (`len%3`) cues on exon headers (heuristic frame cue)
  - shared output semantics across CLI/JS/Lua because all routes call the same
    engine renderer

UniProt mapping capability status:

- shared shell (`gentle_cli shell`, GUI shell): supported via `uniprot ...` commands
  - `uniprot fetch QUERY [--entry-id ID]`
  - `uniprot import-swissprot PATH [--entry-id ID]`
  - `uniprot list` / `uniprot show ENTRY_ID`
  - `uniprot map ENTRY_ID SEQ_ID [--projection-id ID] [--transcript ID]`
  - `uniprot projection-list [--seq SEQ_ID]`
  - `uniprot projection-show PROJECTION_ID`
  - `uniprot feature-coding-dna PROJECTION_ID FEATURE_QUERY [--transcript ID] [--mode genomic_as_encoded|translation_speed_optimized|both] [--speed-profile human|mouse|yeast|ecoli]`
  - `uniprot resolve-ensembl-links PROJECTION_ID [--transcript ID]`
  - `uniprot transcript-accounting PROJECTION_ID [--transcript ID]`
  - `uniprot compare-ensembl-exons PROJECTION_ID [--transcript ID] [--ensembl-entry ID]`
  - `uniprot compare-ensembl-peptide PROJECTION_ID [--transcript ID] [--ensembl-entry ID]`
  - `uniprot audit-projection PROJECTION_ID [--transcript ID] [--ensembl-entry ID] [--report-id ID]`
  - `uniprot audit-parity PROJECTION_ID [--transcript ID] [--ensembl-entry ID] [--report-id ID]`
  - `uniprot audit-list [--seq SEQ_ID]`
  - `uniprot audit-show REPORT_ID`
  - `uniprot audit-export REPORT_ID OUTPUT.json`
  - `uniprot audit-parity-list [--seq SEQ_ID]`
  - `uniprot audit-parity-show REPORT_ID`
  - `uniprot audit-parity-export REPORT_ID OUTPUT.json`
- shared feature-expert route: supported for stored UniProt genome projections
  via the same expert command family used by splicing/isoform inspection
  - `inspect-feature-expert SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]...`
  - `render-feature-expert-svg SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
  - default behavior hides UniProt `CONFLICT` annotations unless they are
    explicitly re-included with `--feature-key CONFLICT`
  - each transcript now gets its own projected protein-coverage rail, so
    missing or partial domains stay transcript-specific instead of every
    isoform inheriting the full reference-protein annotation set
  - when transcript features do not carry `cds_ranges_1based`, GENtle now
    prefers compatible `CDS` features before falling back to exon spans, which
    keeps RefSeq-style records closer to the actual encoded product
  - exported UniProt-projection SVGs now combine:
    - a coordinate-true genomic transcript/exon panel that preserves exon
      positions and introns
    - a shared genomic-exon-column transcript/product panel that reuses one
      exon-family color across aligned transcript segments and translated
      peptide segments
    - isoform-local protein axes, which make skipped-exon domain loss easier to
      read without discarding the genomic panel above
  - those CDS/exon colors are now keyed by genomic exon family/position, not by
    local segment order inside each transcript/product lane, so vertically
    matching locus exons keep the same color even when transcript boundaries
    differ slightly
  - membrane/topology-style features (`SIGNAL`, `TRANSIT`, `TOPO_DOM`,
    `TRANSMEM`, `INTRAMEM`) render in a dedicated lower band beneath the main
    protein rail so dense labels stay readable
- shared shell (`gentle_cli shell`, GUI shell): GenBank accession import
  - `genbank fetch ACCESSION [--as-id ID]`
- shared shell (`gentle_cli shell`, GUI shell): dbSNP-guided annotated region extraction
  - `dbsnp fetch RS_ID GENOME_ID [--flank-bp N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--catalog PATH] [--cache-dir PATH]`
- engine operations behind those commands:
  - `FetchUniprotSwissProt`, `ImportUniprotSwissProt`, `ProjectUniprotToGenome`,
    `FetchGenBankAccession`, `FetchDbSnpRegion`
- reverse-translation is now also available through dedicated shared-shell /
  direct CLI commands:
  - `reverse-translate run PROTEIN_SEQ_ID [--output-id ID] [--speed-profile human|mouse|yeast|ecoli] [--speed-mark fast|slow] [--translation-table N] [--target-anneal-tm-c F] [--anneal-window-bp N]`
  - `reverse-translate list-reports [PROTEIN_SEQ_ID]`
  - `reverse-translate show-report REPORT_ID`
  - `reverse-translate export-report REPORT_ID OUTPUT.json`
  backed by `ReverseTranslateProteinSequence` plus persisted
  reverse-translation report inspection/export routes.
- construct-reasoning protein-to-DNA handoff is now also available through the
  shared shell:
  - `construct-reasoning build-protein-dna-handoff SEQ_ID PROTEIN_SEQ_ID [--transcript TRANSCRIPT_ID] [--projection-id ID] [--ensembl-entry ID] [--feature-query TEXT] [--ranking-goal balanced_provenance|native_fidelity|expression_optimized] [--speed-profile human|mouse|yeast|ecoli] [--speed-mark fast|slow] [--translation-table N] [--target-anneal-tm-c F] [--anneal-window-bp N] [--objective-id ID] [--graph-id ID]`
  - `construct-reasoning list-graphs [SEQ_ID]`
  - `construct-reasoning show-graph GRAPH_ID`
  - `construct-reasoning set-annotation-status GRAPH_ID ANNOTATION_ID draft|accepted|rejected|locked`
  - `construct-reasoning write-annotation GRAPH_ID ANNOTATION_ID`
  - `construct-reasoning export-graph GRAPH_ID OUTPUT.json`
  backed by `BuildProteinToDnaHandoffReasoning` plus the existing persisted
  construct-reasoning graph store.
  - `construct-reasoning show-graph` now returns both the full portable graph
    and a compact CLI-facing `summary` block:
    - `summary_lines`
    - `warning_lines`
    - `fact_summaries`
    - current fact summaries now include adapter-capture review plus the new
      similarity-derived predictor rows for PCR/amplification,
      nanopore/direct-sequencing, repeat-driven mapping, and cloning stability
  - `construct-reasoning set-annotation-status` updates one persisted
    annotation candidate in place and returns:
    - the updated portable graph
    - the updated `annotation_candidate`
    - the same shared `summary` block exposed by `construct-reasoning show-graph`
  - `construct-reasoning write-annotation` materializes one accepted or locked
    generated annotation candidate as an ordinary sequence feature and returns:
    - the refreshed portable graph
    - the refreshed or already-backed `annotation_candidate`
    - a `writeback` report (`gentle.annotation_candidate_writeback.v1`) with
      `created`, `already_present`, and the resulting `feature_id`
    - the same shared `summary` block exposed by `construct-reasoning show-graph`
  - for adapter/linker restriction-capture facts, that summary now surfaces the
    same shared-engine review as the GUI inspector:
    - per-motif conflicts across capture plus retrieval sites
    - the stronger `all named adapter motifs already occur on the insert` case
    - `methylation may rescue this, but the enzyme-specific methylation
      knowledge base is still incomplete`
- additional protein-sequence operations are still available through the shared
  engine contract (`apply_operation` / workflow JSON):
  - `ImportUniprotEntrySequence`
  - `DeriveProteinSequences`
    - transcript-first and self-sufficient; UniProt is optional comparison
      evidence, not the source of truth for what peptide products exist

## Build and run

From the repository root:

```bash
cargo run --bin gentle
cargo run --features js-interface --bin gentle_js
cargo run --features lua-interface --bin gentle_lua
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_mcp -- --help
cargo run --bin gentle -- --version
cargo run --bin gentle_cli -- --version
```

For optimized builds:

```bash
cargo run --release --bin gentle
cargo run --release --features js-interface --bin gentle_js
cargo run --release --features lua-interface --bin gentle_lua
cargo run --release --bin gentle_cli -- capabilities
cargo run --release --bin gentle_examples_docs -- --check
cargo run --release --bin gentle_examples_docs -- tutorial-check
cargo run --release --bin gentle_mcp -- --help
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- --version
```

Important: with `cargo run`, arguments for GENtle binaries must come after `--`.
Example: `cargo run --bin gentle -- --version`.

JavaScript and Lua shells are compile-time optional so default `cargo check`
and `cargo build` do not pull in their heavier runtime dependencies. The Python
wrapper in `integrations/python/gentle_py` remains separate from Cargo
features.

Release-build policy:

- default local builds stay lean unless JS/Lua support is requested explicitly
- release packaging builds enable `--features script-interfaces`
- for a local optimized build that matches release behavior, use
  `cargo build --release --features script-interfaces`

## Protocol-first example source

Canonical workflow examples are adapter-neutral JSON files:

- source: `docs/examples/workflows/*.json`
- generated adapter snippets: `docs/examples/generated/*.md`
- tutorial runtime manifest: `docs/tutorial/manifest.json`
  (generated from `docs/tutorial/sources/`)
- committed generated tutorial output: `docs/tutorial/generated/`

Regenerate snippets on demand:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Validate example files and schema without writing output:

```bash
cargo run --bin gentle_examples_docs -- --check
```

Generate tutorial markdown + retained artifacts:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
```

Validate committed tutorial-generated output:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
```

## `gentle` (GUI launcher)

`gentle` starts the graphical application.

```bash
cargo run --bin gentle
cargo run --bin gentle -- path/to/project.gentle.json
cargo run --bin gentle -- --version
```

If a project path is provided, GENtle opens that project at startup.

Current behavior:

- opens GUI windows
- opens an empty project by default
- optionally opens a project path passed on startup
- supports project/sequence actions from the File menu

## `gentle_js` (JavaScript shell)

`gentle_js` starts an interactive JavaScript REPL backed by GENtle data structures.
This binary is available when Cargo feature `js-interface` is enabled.

```bash
cargo run --features js-interface --bin gentle_js
```

Exit methods:

- type `exit`
- `Ctrl+C`
- `Ctrl+D`

### JavaScript shell functions

1. `load_dna(path)`
   - Loads a DNA sequence from file.
2. `write_gb(seq, path)`
   - Writes a sequence object to a GenBank file.
3. `load_project(path)` / `save_project(state, path)`
   - Loads/saves GENtle project JSON.
4. `capabilities()`
   - Returns shared-engine capabilities.
5. `state_summary(state)`
   - Returns normalized sequence/container/display summary.
6. `inspect_dna_ladders(name_filter)`
   - Returns built-in ladder catalog as structured JSON.
   - `name_filter` is optional (`null`/`""` means all ladders).
7. `export_dna_ladders(path, name_filter)`
   - Writes ladder catalog JSON to `path`.
   - Optional `name_filter` limits exported ladders by case-insensitive name match.
8. `apply_operation(state, op)`
   - Applies one engine operation to a project state.
   - `op` may be a JS object or JSON string.
   - Returns `{ state, result }`.
9. `apply_workflow(state, workflow)`
   - Applies a workflow to a project state.
   - `workflow` may be a JS object or JSON string.
   - Returns `{ state, results }`.
10. `import_pool(state, input_pool_json, prefix)`
   - Imports a `.pool.gentle.json` artifact into `state` via shared adapter logic.
   - `prefix` is optional (`null`/`""` defaults to `pool`).
   - Returns `{ state, state_changed, output }`.
11. `sync_rebase(input, output, commercial_only)`
   - Parses REBASE/Bairoch input and writes a REBASE resource JSON snapshot.
   - `output` is optional (`null`/`""` uses default runtime resource path).
12. `sync_jaspar(input, output)`
   - Parses JASPAR PFM text and writes motif resource JSON snapshot.
   - `output` is optional (`null`/`""` uses default runtime resource path).
13. `list_reference_genomes(catalog_path)`
    - Lists genome IDs from the genome catalog.
    - `catalog_path` is optional (`null`/`""` uses default catalog).
14. `list_reference_catalog_entries(catalog_path, filter)`
    - Lists structured reference-catalog entries, including typed metadata.
    - `filter` is optional and matches the same search surface as `genomes list --filter`.
15. `list_helper_catalog_entries(catalog_path, filter)`
    - Lists structured helper-catalog entries, including optional normalized
      `interpretation` records.
    - `filter` is optional and matches the same search surface as `helpers list --filter`.
16. `list_host_profile_catalog_entries(catalog_path, filter)`
    - Lists structured host-profile catalog entries for construct-reasoning host/strain lookup.
    - `filter` is optional and matches ids, aliases, species, strain, genotype tags, phenotype tags, and notes.
17. `list_ensembl_installable_genomes(collection, filter)`
    - Lists Ensembl species directories that currently appear installable because both FASTA and GTF listings are present.
    - `collection` is optional (`all`, `vertebrates`, `metazoa`).
18. `list_construct_reasoning_graphs(state, seq_id)`
    - Lists stored construct-reasoning graphs plus compact shared-shell summary rows.
    - `seq_id` is optional and limits rows to one active sequence.
19. `show_construct_reasoning_graph(state, graph_id)`
    - Returns one stored construct-reasoning graph plus the same compact summary block exposed by `construct-reasoning show-graph`.
20. `set_construct_reasoning_annotation_status(state, graph_id, annotation_id, status)`
    - Updates one stored annotation-candidate review status and returns the updated state plus shared-shell mutation output.
21. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
21. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
22. `prepare_genome(state, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
23. `extract_genome_region(state, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
24. `extract_genome_gene(state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features, extract_mode, promoter_upstream_bp)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
25. `blast_reference_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json)`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
    - `options_json` is optional and may override any quick option (`max_hits`, `task`) plus threshold fields.
26. `blast_helper_genome(helper_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json)`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.
27. `set_parameter(state, name, value)`
    - Convenience wrapper for engine `SetParameter`.
28. `set_vcf_display_filter(state, options)`
    - Convenience wrapper that updates one or more VCF display-filter parameters
      via `SetParameter` in one call.
29. `list_agent_systems(catalog_path)`
    - Lists configured agent systems from catalog JSON.
    - `catalog_path` is optional (`null`/`""` uses `assets/agent_systems.json`).
30. `ask_agent_system(state, system_id, prompt, options)`
    - Invokes one configured agent system through shared shell execution.
    - Returns `{ state, state_changed, output }`.
    - `state` may be `null` to use an empty/default project state.
    - `options` supports:
      - `catalog_path`
      - `allow_auto_exec` / `execute_all`
      - `execute_indices` (1-based)
      - `include_state_summary` (default `true`)
25. `render_dotplot_svg(state, seq_id, dotplot_id, output_svg, options)`
    - Convenience wrapper around engine `RenderDotplotSvg`.
    - `options` supports:
      - `flex_track_id` (or `flexTrackId`)
      - `display_density_threshold`
      - `display_intensity_gain`

### JavaScript example

Use generated adapter snippets to stay synchronized with canonical workflow JSON:

- `docs/examples/generated/load_and_digest_pgex.md`
- `docs/examples/generated/load_branch_reverse_complement_pgex_fasta.md`
- `docs/examples/generated/guides_filter_and_generate_oligos.md`
- `docs/examples/generated/digest_ligation_extract_region_minimal.md`
- `docs/examples/generated/contribute_gentle_development_baseline.md`
- `docs/examples/generated/guides_export_csv_and_protocol.md`

## `gentle_mcp` (MCP stdio server)

`gentle_mcp` starts a Model Context Protocol server over stdio.
It provides both tool execution (`tools/call`) and capability
discovery/negotiation (`tools/list`, `capabilities`, `help`).

Current tools:

- `capabilities`
- `state_summary`
- `op` (apply one operation; requires `confirm=true`)
- `workflow` (apply one workflow; requires `confirm=true`)
- `help`
- `reference_catalog_entries` (structured reference catalog rows via the shared `genomes list` contract)
- `helper_catalog_entries` (structured helper catalog rows with normalized helper `interpretation` when available)
- `host_profile_catalog_entries` (structured host-profile catalog rows via the shared `hosts list` contract)
- `ensembl_installable_genomes` (shared Ensembl discovery report for currently installable candidates)
- `helper_interpretation` (direct helper-construct interpretation lookup by id or alias)
- `ui_intents` (discover deterministic UI-intent contracts)
- `ui_intent` (run deterministic `ui open|focus` intent resolution path)
- `ui_prepared_genomes` (run deterministic prepared-genome query path)
- `ui_latest_prepared` (resolve latest prepared genome for one species)
- `blast_async_start` (start async BLAST job through shared shell route)
- `blast_async_status` (poll async BLAST job status and optional report)
- `blast_async_cancel` (request cancellation for async BLAST job)
- `blast_async_list` (list known async BLAST jobs in scope)

Tool-parity rule:

- MCP UI-intent tools execute through the same shared parser/executor route as
  CLI shell commands (`ui ...`), and return the same structured payload
  contracts.
- No MCP-only biology/UI logic branches are allowed.

Data-boundary and trust model:

- `gentle_mcp` reads/writes project state via the resolved `state_path`
  (default `.gentle_state.json`).
- It does not attach to another running GUI process's in-memory engine state.
- Unsaved GUI-only edits are therefore not visible to MCP until they are
  persisted to disk and read from that same state path.
- If MCP is pointed at the same project state file used by active GUI/CLI
  workflows, that shared visibility is intentional.
- Treat connected MCP clients as trusted local operators; `confirm=true`
  expresses explicit mutating intent but is not an authentication boundary.

Run:

```bash
cargo run --bin gentle_mcp
cargo run --bin gentle_mcp -- --state path/to/project.gentle.json
```

Minimum MCP JSON-RPC flow:

1. Initialize:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}
```

2. List tools:

```json
{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}
```

3. Call one tool:

```json
{
  "jsonrpc":"2.0",
  "id":3,
  "method":"tools/call",
  "params":{
    "name":"ui_intent",
    "arguments":{
      "action":"open",
      "target":"prepared-references",
      "catalog_path":"assets/genomes.json",
      "species":"human",
      "latest":true
    }
  }
}
```

`ui_intent` arguments:

- required:
  - `action`: `open|focus`
  - `target`: one of:
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
- optional:
  - `state_path`
  - `genome_id`
  - `helpers`
  - `catalog_path`
  - `cache_dir`
  - `filter`
  - `species`
  - `latest`

`ui_prepared_genomes` arguments:

- optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`, `filter`,
  `species`, `latest`

`ui_latest_prepared` arguments:

- required: `species`
- optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`

`reference_catalog_entries` arguments:

- optional: `catalog_path`, `filter`

`helper_catalog_entries` arguments:

- optional: `catalog_path`, `filter`

`host_profile_catalog_entries` arguments:

- optional: `catalog_path`, `filter`

`ensembl_installable_genomes` arguments:

- optional: `collection`, `filter`

`helper_interpretation` arguments:

- required: `helper_id` (catalog id or alias)
- optional: `catalog_path`

`blast_async_start` arguments:

- required: `genome_id`, `query_sequence`
- optional: `state_path`, `helpers`, `max_hits`, `task`, `options_json`,
  `catalog_path`, `cache_dir`

`blast_async_status` arguments:

- required: `job_id`
- optional: `helpers`, `with_report`

`blast_async_cancel` arguments:

- required: `job_id`
- optional: `helpers`

`blast_async_list` arguments:

- optional: `helpers`

`ui_intents` arguments:

- optional: `state_path` (accepted for API symmetry)

MCP tool result envelope behavior:

- `result.isError = false`: success; inspect `result.structuredContent`.
- `result.isError = true`: deterministic tool-level error; inspect
  `result.content[0].text`.
- UI-intent structured payloads match shared shell schemas:
  - `gentle.ui_intents.v1`
  - `gentle.ui_intent.v1`
  - `gentle.ui_prepared_genomes.v1`
  - `gentle.ui_latest_prepared.v1`
- catalog-entry MCP tools match the shared shell list payload shape:
  - `catalog_path`
  - `filter`
  - `genome_count`
  - `genomes[]`
  - `entries[]`
- `helper_interpretation` returns:
  - `query`
  - `catalog_path`
  - `interpretation` (`null` when the helper exists but has no structured semantics)

Mutating tool safety:

- `op` and `workflow` reject execution unless `confirm=true` is explicitly set.
- UI-intent tools are currently non-mutating; if a routed command ever reports
  state mutation unexpectedly, MCP returns an explicit tool error.

## `gentle_lua` (Lua shell)

`gentle_lua` starts an interactive Lua REPL backed by GENtle data structures.
This binary is available when Cargo feature `lua-interface` is enabled.

```bash
cargo run --features lua-interface --bin gentle_lua
```

Exit methods:

- type `exit`
- `Ctrl+C`
- `Ctrl+D`

### Lua shell functions

1. `load_dna(filename)`
   - Loads a DNA sequence from file.
2. `write_gb(seq, filename)`
   - Writes a sequence object to a GenBank file.
3. `load_project(filename)` / `save_project(state, filename)`
   - Loads/saves GENtle project JSON.
4. `capabilities()`
   - Returns shared-engine capabilities.
5. `state_summary(project)`
   - Returns normalized sequence/container/display summary.
6. `inspect_dna_ladders([name_filter])`
   - Returns built-in ladder catalog as structured table.
7. `export_dna_ladders(output_json, [name_filter])`
   - Writes ladder catalog JSON to disk.
8. `apply_operation(project, op)`
   - Applies one engine operation; `op` can be Lua table or JSON string.
   - Returns table with `state` and `result`.
9. `apply_workflow(project, workflow)`
   - Applies workflow; `workflow` can be Lua table or JSON string.
   - Returns table with `state` and `results`.
10. `import_pool(project, input_pool_json, [prefix])`
   - Imports a `.pool.gentle.json` artifact into `project` via shared adapter logic.
   - Returns table with `state`, `state_changed`, and `output`.
11. `sync_rebase(input, output, commercial_only)`
   - Parses REBASE/Bairoch input and writes a REBASE resource JSON snapshot.
   - `output` and `commercial_only` are optional.
12. `sync_jaspar(input, output)`
   - Parses JASPAR PFM text and writes motif resource JSON snapshot.
   - `output` is optional.
13. `list_reference_genomes(catalog_path)`
    - Lists genome IDs from the genome catalog.
14. `list_reference_catalog_entries([catalog_path], [filter])`
    - Lists structured reference-catalog entries, including typed metadata.
15. `list_helper_catalog_entries([catalog_path], [filter])`
    - Lists structured helper-catalog entries, including optional normalized
      `interpretation` records.
16. `list_host_profile_catalog_entries([catalog_path], [filter])`
    - Lists structured host-profile catalog entries for construct-reasoning host/strain lookup.
17. `list_ensembl_installable_genomes([collection], [filter])`
    - Lists current Ensembl candidates where both FASTA and GTF species listings exist.
18. `list_construct_reasoning_graphs(project, [seq_id])`
    - Lists stored construct-reasoning graphs plus compact shared-shell summary rows.
19. `show_construct_reasoning_graph(project, graph_id)`
    - Returns one stored construct-reasoning graph plus the same compact summary block exposed by `construct-reasoning show-graph`.
20. `set_construct_reasoning_annotation_status(project, graph_id, annotation_id, status)`
    - Updates one stored annotation-candidate review status and returns the updated project state plus shared-shell mutation output.
21. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
21. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
22. `prepare_genome(project, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
23. `extract_genome_region(project, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
24. `extract_genome_gene(project, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features, extract_mode, promoter_upstream_bp)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
25. `blast_reference_genome(genome_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir], [options_json])`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
    - `options_json` is optional and may override any quick option (`max_hits`, `task`) plus threshold fields.
26. `blast_helper_genome(helper_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir], [options_json])`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.
27. `set_parameter(project, name, value)`
    - Convenience wrapper for engine `SetParameter`.
28. `set_vcf_display_filter(project, opts)`
    - Convenience wrapper that updates one or more VCF display-filter parameters
      via `SetParameter` in one call.
29. `list_agent_systems([catalog_path])`
    - Lists configured agent systems from catalog JSON.
    - `catalog_path` is optional (defaults to `assets/agent_systems.json`).
30. `ask_agent_system(project, system_id, prompt, [catalog_path], [allow_auto_exec], [execute_all], [execute_indices], [include_state_summary])`
    - Invokes one configured agent system through shared shell execution.
    - Returns table with `state`, `state_changed`, and `output`.
    - `project` may be `nil` to use an empty/default project state.
25. `render_dotplot_svg(project, seq_id, dotplot_id, output_svg, [flex_track_id], [display_density_threshold], [display_intensity_gain])`
    - Convenience wrapper around engine `RenderDotplotSvg`.

### Lua example

Use generated adapter snippets to stay synchronized with canonical workflow JSON:

- `docs/examples/generated/load_and_digest_pgex.md`
- `docs/examples/generated/load_branch_reverse_complement_pgex_fasta.md`
- `docs/examples/generated/guides_filter_and_generate_oligos.md`
- `docs/examples/generated/digest_ligation_extract_region_minimal.md`
- `docs/examples/generated/contribute_gentle_development_baseline.md`
- `docs/examples/generated/guides_export_csv_and_protocol.md`

## File format expectations

Current CLI workflows rely on sequence files supported by internal loaders:

- GenBank
- EMBL
- FASTA
  - default interpretation: synthetic blunt `dsDNA`
  - optional FASTA-header metadata tokens for synthetic oligos:
    - `molecule=ssdna` for single-stranded DNA
    - `molecule=rna` for RNA (input `T` is normalized to `U`)
    - `molecule=dsdna` plus optional overhangs:
      - `forward_5=...` (alias `f5=...`)
      - `forward_3=...` (alias `f3=...`)
      - `reverse_5=...` (alias `r5=...`)
      - `reverse_3=...` (alias `r3=...`)
- NCBI GenBank XML (`GBSet/GBSeq`)
  - currently rejected with explicit diagnostics for unsupported XML dialects
    (for example `INSDSet/INSDSeq`)

Example FASTA headers:

- `>oligo_ss molecule=ssdna`
- `>oligo_rna molecule=rna`
- `>oligo_ds molecule=dsdna f5=GATC r5=CTAG`

## `gentle_cli` (automation/agent CLI)

`gentle_cli` exposes operation-based control with JSON input and output.
It is intended for scripted workflows and AI agents.

State is persisted to `.gentle_state.json` by default (override with `--state PATH`).

### Commands

```bash
cargo run --bin gentle_cli -- --version
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_cli -- state-summary
cargo run --bin gentle_cli -- containers set-exclusive container-3 false
cargo run --bin gentle_cli -- help
cargo run --bin gentle_cli -- help candidates generate
cargo run --bin gentle_cli -- help --format json
cargo run --bin gentle_cli -- op '<operation-json>'
cargo run --bin gentle_cli -- op op.json
cargo run --bin gentle_cli -- workflow '<workflow-json>'
cargo run --bin gentle_cli -- workflow docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/rna_reads_interpret_cdna_tp73_template.json
cargo run --bin gentle_cli -- --progress op '<operation-json>'
cargo run --bin gentle_cli -- --progress-stdout workflow '<workflow-json>'
cargo run --bin gentle_cli -- export-state state.json
cargo run --bin gentle_cli -- import-state state.json
cargo run --bin gentle_cli -- save-project project.gentle.json
cargo run --bin gentle_cli -- load-project project.gentle.json
cargo run --bin gentle_cli -- render-svg pgex linear pgex.linear.svg
cargo run --bin gentle_cli -- render-svg pgex circular pgex.circular.svg
cargo run --bin gentle_cli -- render-dotplot-svg tp73_cdna dotplot_primary tp73.dotplot.svg
cargo run --bin gentle_cli -- render-dotplot-svg tp73_cdna dotplot_primary tp73.dotplot.svg --flex-track flex_primary --display-threshold 0.2 --intensity-gain 1.8
cargo run --bin gentle_cli -- render-dotplot-svg tp73_genomic_local tp73_multi_isoform_overlay_dotplot tp73.overlay.right.svg --overlay-x-axis right_aligned_bp --display-threshold 0.08 --intensity-gain 2.8
cargo run --bin gentle_cli -- render-dotplot-svg tp73_genomic_local tp73_multi_isoform_overlay_dotplot tp73.overlay.anchor.svg --overlay-x-axis shared_exon_anchor --overlay-anchor-exon 55034..55276 --display-threshold 0.08 --intensity-gain 2.8
cargo run --bin gentle_cli -- render-rna-svg rna_seq rna.secondary.svg
cargo run --bin gentle_cli -- rna-info rna_seq
cargo run --bin gentle_cli -- render-lineage-svg lineage.svg
cargo run --bin gentle_cli -- protocol-cartoon list
cargo run --bin gentle_cli -- protocol-cartoon render-svg gibson.two_fragment gibson.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.pair pcr_pair.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.pair.with_tail pcr_pair_with_tail.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.oe.substitution pcr_oe_substitution.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.qpcr qpcr.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-template-svg docs/examples/protocol_cartoon/demo_template.json demo.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon template-validate docs/examples/protocol_cartoon/demo_template.json
cargo run --bin gentle_cli -- protocol-cartoon render-with-bindings docs/examples/protocol_cartoon/demo_template.json docs/examples/protocol_cartoon/demo_bindings.json demo.bound.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon template-export gibson.two_fragment gibson.template.json
cargo run --bin gentle_cli -- gibson preview @docs/examples/plans/gibson_destination_first_single_insert.json --output gibson.preview.json
cargo run --bin gentle_cli -- gibson apply @docs/examples/plans/gibson_destination_first_single_insert.json
cargo run --bin gentle_cli -- shell 'help'
cargo run --bin gentle_cli -- shell 'state-summary'
cargo run --bin gentle_cli -- shell 'op @op.json'
cargo run --bin gentle_cli -- render-pool-gel-svg frag_1,frag_2 digest.gel.svg
cargo run --bin gentle_cli -- render-pool-gel-svg frag_1,frag_2 digest.gel.svg --ladders "NEB 100bp DNA Ladder,NEB 1kb DNA Ladder"
cargo run --bin gentle_cli -- render-pool-gel-svg - digest.gel.svg --containers container-3,container-8
cargo run --bin gentle_cli -- render-pool-gel-svg - digest.gel.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- render-gel-svg - digest.gel.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- arrange-serial container-3,container-8 --id arrangement-2 --name "Digest run A" --ladders "NEB 100bp DNA Ladder"
cargo run --bin gentle_cli -- racks create-from-arrangement arrangement-2 --rack-id rack-1 --name "Bench day A"
cargo run --bin gentle_cli -- racks place-arrangement arrangement-3 --rack rack-1
cargo run --bin gentle_cli -- racks move rack-1 --from A2 --to A4
cargo run --bin gentle_cli -- racks move rack-1 --from A1 --to B1 --block
cargo run --bin gentle_cli -- racks move-samples rack-1 --from A1 --from A3 --to A2
cargo run --bin gentle_cli -- racks move-blocks rack-1 --arrangement arrangement-2 --arrangement arrangement-3 --to B2
cargo run --bin gentle_cli -- racks show rack-1
cargo run --bin gentle_cli -- racks labels-svg rack-1 arrangement-2.labels.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- racks fabrication-svg rack-1 rack.fabrication.svg --template storage_pcr_tube_rack
cargo run --bin gentle_cli -- racks openscad rack-1 rack.scad --template pipetting_pcr_tube_rack
cargo run --bin gentle_cli -- racks set-profile rack-1 plate_96
cargo run --bin gentle_cli -- ladders list
cargo run --bin gentle_cli -- ladders list --filter NEB
cargo run --bin gentle_cli -- ladders list --molecule rna
cargo run --bin gentle_cli -- ladders export dna_ladders.snapshot.json
cargo run --bin gentle_cli -- ladders export dna_ladders.neb.json --filter NEB
cargo run --bin gentle_cli -- ladders export rna_ladders.snapshot.json --molecule rna
cargo run --bin gentle_cli -- export-pool frag_1,frag_2 digest.pool.gentle.json "BamHI+EcoRI digest pool"
cargo run --bin gentle_cli -- import-pool digest.pool.gentle.json imported
cargo run --bin gentle_cli -- services status
cargo run --bin gentle_cli -- resources status
cargo run --bin gentle_cli -- resources sync-rebase rebase.withrefm data/resources/rebase.enzymes.json --commercial-only
cargo run --bin gentle_cli -- resources sync-jaspar JASPAR2026_CORE_non-redundant_pfms_jaspar.txt data/resources/jaspar.motifs.json
cargo run --bin gentle_cli -- resources sync-jaspar-remote-metadata --filter TP --limit 50 data/resources/jaspar.remote_metadata.json
cargo run --bin gentle_cli -- resources summarize-jaspar --motif SP1 --motif REST --random-length 10000 --seed 123 --output jaspar.summary.json
cargo run --bin gentle_cli -- resources benchmark-jaspar --random-length 10000 --seed 123 --output data/resources/jaspar.registry_benchmark.json
cargo run --bin gentle_cli -- resources list-jaspar --filter TP --limit 50 --output jaspar.catalog.json
cargo run --bin gentle_cli -- resources inspect-jaspar SP1 --random-length 10000 --seed 123 --fetch-remote --output jaspar.expert.json
scripts/benchmark_jaspar_interest_catalog.sh run --motif SP1 --out-dir data/resources/jaspar_interest_catalog
scripts/benchmark_jaspar_interest_catalog.sh merge --out-dir data/resources/jaspar_interest_catalog
cargo run --bin gentle_cli -- agents list
cargo run --bin gentle_cli -- agents list --catalog assets/agent_systems.json
cargo run --bin gentle_cli -- agents ask builtin_echo --prompt "summarize current project state"
cargo run --bin gentle_cli -- agents ask builtin_echo --prompt "ask: Which sequence should I use?" --execute-index 1
cargo run --bin gentle_cli -- agents ask local_llama_compat --prompt "summarize project context" --base-url http://localhost:11964 --model deepseek-r1:8b
cargo run --bin gentle_cli -- op '{"PrepareGenome":{"genome_id":"ToyGenome","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeRegion":{"genome_id":"ToyGenome","chromosome":"chr1","start_1based":1001,"end_1based":1600,"output_id":"toy_chr1_1001_1600","annotation_scope":"core","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeGene":{"genome_id":"ToyGenome","gene_query":"MYGENE","occurrence":1,"output_id":"toy_mygene","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"DesignInsertionPrimerPairs":{"template":"toy_chr1_1001_1600","insertion":{"requested_forward_3prime_end_0based_exclusive":220,"requested_reverse_3prime_start_0based":300,"forward_extension_5prime":"GAATTC","reverse_extension_5prime":"CTCGAG","forward_window_start_0based":170,"forward_window_end_0based_exclusive":245,"reverse_window_start_0based":275,"reverse_window_end_0based_exclusive":360,"max_anchor_shift_bp":12},"forward":{"min_length":20,"max_length":30},"reverse":{"min_length":20,"max_length":30},"pair_constraints":{"require_roi_flanking":false},"min_amplicon_bp":120,"max_amplicon_bp":1200,"max_tm_delta_c":2.0,"max_pairs":50,"report_id":"insertion_demo_v1"}}'
cargo run --bin gentle_cli -- genomes list --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes list --catalog assets/helper_genomes.json
cargo run --bin gentle_cli -- genomes validate-catalog --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes update-ensembl-specs --catalog assets/genomes.json --output-catalog exports/genomes.updated.json
cargo run --bin gentle_cli -- genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP53$" --biotype protein_coding --limit 20
cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 3600
cargo run --bin gentle_cli -- genomes remove-prepared "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes remove-catalog-entry "Human GRCh38 Ensembl 113" --catalog exports/genomes.custom.json
cargo run --bin gentle_cli -- genomes blast "Human GRCh38 Ensembl 116" ACGTACGTACGT --task blastn-short --max-hits 10 --options-json '{"thresholds":{"min_identity_percent":97.0,"min_query_coverage_percent":80.0}}' --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-region "Human GRCh38 Ensembl 116" 1 1000000 1001500 --output-id grch38_chr1_slice --annotation-scope core --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- tracks import-bed grch38_tp53 data/chipseq/peaks.bed.gz --name H3K27ac --min-score 10 --clear-existing
cargo run --bin gentle_cli -- tracks import-bigwig grch38_tp53 data/chipseq/signal.bw --name ATAC --min-score 0.2 --clear-existing
cargo run --bin gentle_cli -- tracks import-vcf grch38_tp53 data/variants/sample.vcf.gz --name Variants --min-score 20 --clear-existing
cargo run --bin gentle_cli -- helpers list
cargo run --bin gentle_cli -- helpers validate-catalog
cargo run --bin gentle_cli -- helpers update-ensembl-specs --catalog assets/helper_genomes.json --output-catalog exports/helper_genomes.updated.json
cargo run --bin gentle_cli -- helpers status "Plasmid pUC19 (online)"
cargo run --bin gentle_cli -- helpers prepare "Plasmid pUC19 (online)" --cache-dir data/helper_genomes --timeout-secs 600
cargo run --bin gentle_cli -- helpers remove-prepared "Plasmid pUC19 (online)" --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- helpers genes "Plasmid pUC19 (online)" --filter bla --limit 20
cargo run --bin gentle_cli -- helpers blast "Plasmid pUC19 (online)" ACGTACGTACGT --task blastn-short --max-hits 10 --options-json '{"thresholds":{"min_identity_percent":95.0}}' --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- cache inspect --references --cache-dir data/genomes
cargo run --bin gentle_cli -- cache clear derived-indexes-only --references --cache-dir data/genomes --prepared-id "Human GRCh38 Ensembl 116"
cargo run --bin gentle_cli -- cache clear selected-prepared --both --cache-dir data/genomes --cache-dir data/helper_genomes --prepared-path data/helper_genomes/localproject
cargo run --bin gentle_cli -- cache clear all-prepared-in-cache --both --include-orphans
cargo run --bin gentle_cli -- candidates generate sgrnas chr1_window --length 20 --step 1 --feature-kind gene --max-distance 500 --limit 5000
cargo run --bin gentle_cli -- candidates score sgrnas gc_balance "100 * (gc_fraction - at_fraction)"
cargo run --bin gentle_cli -- candidates score-weighted sgrnas priority --term gc_fraction:0.7:max --term distance_to_seq_start_bp:0.3:min --normalize
cargo run --bin gentle_cli -- candidates top-k sgrnas sgrnas_top --metric priority --k 100 --direction max
cargo run --bin gentle_cli -- candidates pareto sgrnas sgrnas_front --objective gc_fraction:max --objective distance_to_seq_start_bp:min --max-candidates 200
cargo run --bin gentle_cli -- candidates filter sgrnas sgrnas_q95 --metric gc_balance --min-quantile 0.95
cargo run --bin gentle_cli -- candidates macro @candidate_flow.gsh
cargo run --bin gentle_cli -- candidates macro --transactional --file candidate_flow.gsh
cargo run --bin gentle_cli -- candidates template-put scan_tp53 --script 'generate ${set_name} ${seq_id} --length ${len} --step 1' --param set_name --param seq_id=grch38_tp53 --param len=20
cargo run --bin gentle_cli -- candidates template-run scan_tp53 --bind set_name=tp53_candidates --transactional
cargo run --bin gentle_cli -- guides list
cargo run --bin gentle_cli -- guides put tp73_guides --json '[{"guide_id":"g1","seq_id":"tp73","start_0based":100,"end_0based_exclusive":120,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17}]'
cargo run --bin gentle_cli -- guides filter tp73_guides --config '{"gc_min":0.3,"gc_max":0.7,"avoid_u6_terminator_tttt":true}' --output-set tp73_guides_pass
cargo run --bin gentle_cli -- guides oligos-generate tp73_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set tp73_lenti --passed-only
cargo run --bin gentle_cli -- guides oligos-export tp73_guides exports/tp73_guides.csv --format csv_table --oligo-set tp73_lenti
cargo run --bin gentle_cli -- guides protocol-export tp73_guides exports/tp73_guides.protocol.txt --oligo-set tp73_lenti
cargo run --bin gentle_cli -- shell 'macros template-list'
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns.json'
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns_catalog'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --transactional'
cargo run --bin gentle_cli -- shell 'macros run --transactional --file cloning_flow.gsh'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_pass_only true'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_required_info_keys ["AF","DP"]'
cargo run --bin gentle_cli -- shell 'set-param tfbs_display_min_llr_quantile 0.95'
cargo run --bin gentle_cli -- shell 'set-param show_restriction_enzymes true'
cargo run --bin gentle_cli -- shell 'set-param restriction_enzyme_display_mode "preferred_only"'
cargo run --bin gentle_cli -- shell 'set-param preferred_restriction_enzymes ["EcoRI","BamHI"]'
cargo run --bin gentle_cli -- shell 'panels import-isoform grch38_tp53 assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1'
cargo run --bin gentle_cli -- shell 'panels inspect-isoform grch38_tp53 tp53_isoforms_v1'
cargo run --bin gentle_cli -- shell 'panels render-isoform-svg grch38_tp53 tp53_isoforms_v1 exports/tp53_isoform_architecture.svg'
cargo run --bin gentle_cli -- shell 'panels validate-isoform assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1'
cargo run --bin gentle_cli -- inspect-feature-expert grch38_tp53 isoform tp53_isoforms_v1
cargo run --bin gentle_cli -- render-feature-expert-svg grch38_tp53 isoform tp53_isoforms_v1 exports/tp53_isoform_architecture.svg
```

You can pass JSON from a file with `@file.json` or a bare existing file path.
When loading from file path, an initial shebang line (`#!...`) is ignored so
executable script files can embed JSON payloads directly.
`workflow` accepts both raw workflow payloads (`{"run_id":"...","ops":[...]}`)
and wrapped protocol example payloads (`{"workflow":{...}}`).

Global CLI options:

- `--state PATH`: use a non-default project state file
- `--progress` or `--progress-stderr`: print live progress events to `stderr`
- `--progress-stdout`: print live progress events to `stdout`
- `--allow-screenshots`: currently rejected (screenshot bridge disabled by security policy)

Current progress events include TFBS annotation updates, genome-prepare
updates (download/index phases), genome-track import updates, RNA-read
interpretation/alignment updates, and primer/qPCR design updates.
Primer/qPCR progress lines expose backend choice/fallback plus candidate and
search counts, for example `progress primers ... stage=pair_search ...` and
`progress primers ... stage=assay_search_complete ...`.
When `--progress-stdout` is used, progress lines are emitted before the final JSON output.

`state-summary` output includes:

- sequences
- containers
  - each container row now includes `declared_contents_exclusive`
- display visibility flags

Project aliases:

- `save-project PATH` aliases `export-state PATH`
- `load-project PATH` aliases `import-state PATH`

Shared shell command:

- `shell '<command line>'`
  - Parses and executes commands through the same shared parser/executor used by
    the GUI `Shell` panel (`src/engine_shell.rs`).
  - Shell-supported commands:
    - `help`
    - `capabilities`
    - `state-summary`
    - `containers set-exclusive CONTAINER_ID true|false`
    - `load-project PATH`
    - `save-project PATH`
    - `render-svg SEQ_ID linear|circular OUTPUT.svg`
    - `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor] [--overlay-anchor-exon START..END]`
    - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor] [--overlay-anchor-exon START..END]`
    - `dotplot overlay-compute OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - `render-rna-svg SEQ_ID OUTPUT.svg`
    - `rna-info SEQ_ID`
    - `render-lineage-svg OUTPUT.svg`
    - `protocol-cartoon list`
    - `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
    - `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
    - `protocol-cartoon template-validate TEMPLATE.json`
    - `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
    - `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
    - `render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]`
    - `render-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]`
    - `arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]`
    - `arrange-set-ladders ARR_ID [--ladders NAME[,NAME]]`
    - `racks create-from-arrangement ARR_ID [--rack-id ID] [--name TEXT] [--profile small_tube_4x6|plate_96|plate_384]`
    - `racks place-arrangement ARR_ID --rack RACK_ID`
    - `racks move RACK_ID --from A1 --to B1 [--block]`
    - `racks move-blocks RACK_ID --arrangement ARR_ID [--arrangement ARR_ID ...] --to B1`
    - `racks show RACK_ID`
    - `racks labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--preset compact_cards|print_a4|wide_cards]`
    - `racks fabrication-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
    - `racks isometric-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
    - `racks openscad RACK_ID OUTPUT.scad [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
    - `racks set-profile RACK_ID small_tube_4x6|plate_96|plate_384`
    - `racks set-custom-profile RACK_ID ROWS COLUMNS`
    - `ladders list [--molecule dna|rna] [--filter TEXT]`
    - `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
    - `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
    - `export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]`
    - `import-pool INPUT.pool.gentle.json [PREFIX]`
    - `resources sync-rebase INPUT.withrefm_or_URL [OUTPUT.rebase.json] [--commercial-only]`
    - `resources sync-jaspar INPUT.jaspar_or_URL [OUTPUT.motifs.json]`
    - `resources sync-jaspar-remote-metadata [--motif TOKEN ...] [--motifs CSV] [--all] [--filter TOKEN] [--limit N] [--output OUTPUT.json]`
    - `resources summarize-jaspar [--motif TOKEN ...] [--motifs CSV] [--all] [--random-length N] [--seed N] [--output OUTPUT.json]`
    - `resources benchmark-jaspar [--random-length N] [--seed N] [--output OUTPUT.json]`
    - `resources list-jaspar [--filter TOKEN] [--limit N] [--fetch-remote] [--output OUTPUT.json]`
    - `resources inspect-jaspar MOTIF [--random-length N] [--seed N] [--fetch-remote] [--output OUTPUT.json]`
    - `agents list [--catalog PATH]`
    - `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
    - `genomes list [--catalog PATH]`
    - `genomes ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
    - `genomes install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
    - `genomes validate-catalog [--catalog PATH]`
    - `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
    - `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
    - `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast-start GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast-status JOB_ID [--with-report]`
    - `genomes blast-cancel JOB_ID`
    - `genomes blast-list`
    - `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `helpers list [--catalog PATH]`
    - `helpers ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
    - `helpers install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
    - `helpers validate-catalog [--catalog PATH]`
    - `hosts list [--catalog PATH] [--filter TEXT]`
    - `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
    - `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
    - `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `helpers blast-start HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `helpers blast-status JOB_ID [--with-report]`
    - `helpers blast-cancel JOB_ID`
    - `helpers blast-list`
    - `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-vcf SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
    - `macros instance-list`
    - `macros instance-show MACRO_INSTANCE_ID`
    - `macros template-list`
    - `macros template-show TEMPLATE_NAME`
    - `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port PORT_ID:KIND[:one|many][:required|optional][:description]]`
    - `macros template-delete TEMPLATE_NAME`
    - `macros template-import PATH`
    - `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]`
    - `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT] [--seq-id SEQ_ID]`
    - `routines explain ROUTINE_ID [--catalog PATH] [--seq-id SEQ_ID]`
    - `routines compare ROUTINE_A ROUTINE_B [--catalog PATH] [--seq-id SEQ_ID]`
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
    - `guides list`
    - `guides show GUIDE_SET_ID [--limit N] [--offset N]`
    - `guides put GUIDE_SET_ID (--json JSON|@FILE|--file PATH)`
    - `guides delete GUIDE_SET_ID`
    - `guides filter GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]`
    - `guides filter-show GUIDE_SET_ID`
    - `guides oligos-generate GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]`
    - `guides oligos-list [--guide-set GUIDE_SET_ID]`
    - `guides oligos-show OLIGO_SET_ID`
    - `guides oligos-export GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]`
    - `guides protocol-export GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]`
    - `features query SEQ_ID [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
    - `features export-bed SEQ_ID OUTPUT.bed [--coordinate-mode auto|local|genomic] [--include-restriction-sites] [--restriction-enzyme NAME] [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
    - `features tfbs-summary SEQ_ID --focus START..END [--context START..END] [--min-focus-count N] [--min-context-count N] [--limit N]`
    - `features tfbs-score-tracks-svg SEQ_ID OUTPUT.svg --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|true_log_odds_bits|true_log_odds_quantile] [--allow-negative]`
    - `features restriction-scan SEQ_ID [--range START..END|--start N --end N] [--enzyme NAME] [--max-sites-per-enzyme N] [--no-cut-geometry] [--path FILE.json]`
    - `features restriction-scan --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] [--range START..END|--start N --end N] [--enzyme NAME] [--max-sites-per-enzyme N] [--no-cut-geometry] [--path FILE.json]`
    - `variant annotate-promoters SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--upstream-bp N] [--downstream-bp N] [--collapse transcript|gene]`
    - `variant promoter-context SEQ_ID [--variant ID] [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--tfbs-focus-half-window-bp N] [--path FILE.json]`
    - `variant reporter-fragments SEQ_ID [--variant ID] [--gene-label LABEL] [--transcript-id ID] [--retain-downstream-from-tss-bp N] [--retain-upstream-beyond-variant-bp N] [--max-candidates N] [--path FILE.json]`
    - `variant materialize-allele SEQ_ID --allele reference|alternate [--variant ID] [--output-id ID]`
    - `primers design REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers design-qpcr REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers preflight [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers seed-from-feature SEQ_ID FEATURE_ID`
    - `primers seed-from-splicing SEQ_ID FEATURE_ID`
    - `primers seed-qpcr-from-feature SEQ_ID FEATURE_ID`
    - `primers seed-qpcr-from-splicing SEQ_ID FEATURE_ID`
    - `primers list-reports`
    - `primers show-report REPORT_ID`
      - includes `simple_pcr_pairs` with per-pair left/right distance from the
        core ROI, overlap flags, and flanking labels for quick CLI inspection
    - `primers export-report REPORT_ID OUTPUT.json`
    - `primers list-qpcr-reports`
    - `primers show-qpcr-report REPORT_ID`
      - persisted qPCR report output now includes `best_assay_summary` plus
        machine-readable `best_assay_probe_placement`
    - `primers export-qpcr-report REPORT_ID OUTPUT.json`
    - `dotplot compute SEQ_ID [--reference-seq REF_SEQ_ID] [--start N] [--end N] [--ref-start N] [--ref-end N] [--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - `dotplot overlay-compute OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - `dotplot list [SEQ_ID]`
    - `dotplot show DOTPLOT_ID`
    - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor] [--overlay-anchor-exon START..END]`
    - `transcripts derive SEQ_ID [--feature-id N ...] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
    - `flex compute SEQ_ID [--start N] [--end N] [--model at_richness|at_skew] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]`
    - `flex list [SEQ_ID]`
    - `flex show TRACK_ID`
    - `splicing-refs derive SEQ_ID START_0BASED END_0BASED [--seed-feature-id N] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
    - `align compute QUERY_SEQ_ID TARGET_SEQ_ID [--query-start N] [--query-end N] [--target-start N] [--target-end N] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N]`
    - `rna-reads interpret SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile nanopore_cdna_v1] [--format fasta] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--seed-stride-bp N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
    - `rna-reads align-report REPORT_ID [--selection all|seed_passed|aligned] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
    - `rna-reads list-reports [SEQ_ID]`
    - `rna-reads show-report REPORT_ID`
    - `rna-reads summarize-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--output PATH]`
    - `rna-reads inspect-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--cohort all|accepted|fragment|complete|rejected] [--output PATH]`
    - `rna-reads inspect-alignments REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--effect-filter all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only] [--sort rank|identity|coverage|score] [--search TEXT] [--record-indices i,j,k] [--score-bin-variant all_scored|composite_seed_gate] [--score-bin-index N] [--score-bin-count M]`
    - `rna-reads export-report REPORT_ID OUTPUT.json`
    - `rna-reads export-hits-fasta REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-sample-sheet OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--gene GENE_ID]... [--complete-rule near|strict|exact] [--append]`
    - `rna-reads export-paths-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-abundance-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-score-density-svg REPORT_ID OUTPUT.svg [--scale linear|log] [--variant all_scored|composite_seed_gate]`
    - `rna-reads export-alignments-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-alignment-dotplot-svg REPORT_ID OUTPUT.svg [--selection all|seed_passed|aligned] [--max-points N]`
    - `rna-reads align-report` re-ranks retained hits by alignment-aware
      retention rank after mapping refresh.
    - `rna-reads inspect-alignments` returns ranked aligned rows suitable for
      read-level inspection without mutating report payloads, and now echoes
      the formal structured subset spec (`effect_filter`, `sort_key`,
      `search`, `selected_record_indices`, `score_density_variant`,
      `score_bin_index`, `score_bin_count`) plus `subset_match_count` in the
      JSON payload so agent-driven inspection stays reproducible. Rows now also
      include full-length flags:
      - `full_length_exact`
      - `full_length_near`
      - `full_length_strict`
    - `rna-reads summarize-gene-support` returns a machine-readable
      `gentle.rna_read_gene_support_summary.v1` payload for one or more target
      genes from a saved aligned report:
      - base cohort = aligned retained rows, optionally narrowed by
        `--record-indices`
      - accepted target cohort = base rows whose best mapping resolves to one
        of the requested genes
      - cohort tables are emitted for `all_target`, `fragments`, and
        `complete`
      - `exon_pair_support` includes skipped ordered pairs like `1->3`
      - `direct_transition_support` includes neighboring exon steps only
      - `--output PATH` writes the same JSON returned on stdout
    - `rna-reads inspect-gene-support` returns a row-level
      `gentle.rna_read_gene_support_audit.v1` payload for the same saved-report
      target-gene evaluation:
      - `status` is one of `unaligned`, `aligned_other_gene`,
        `accepted_fragment`, or `accepted_complete`
      - grouped top-level arrays echo the exact accepted/fragment/complete
        record indices regardless of the chosen row `--cohort` filter
      - rows carry machine-readable `status_reason`, full-length flags/class,
        mapped exon ordinals, exon pairs, direct transitions, and phase-2
        alignment score/identity/coverage fields
    - `rna-reads export-alignments-tsv` writes the same ranked alignment rows
      in TSV form for downstream filtering/sorting; `--record-indices`
      exports an exact saved-report subset and overrides coarse `--selection`;
      `--subset-spec` records the formal subset definition that produced it.
      TSV rows now include:
      - `target_coverage_fraction`
      - `full_length_exact`
      - `full_length_near`
      - `full_length_strict`
      - `full_length_class` (`exact` / `strict_end` / `near` / `partial`)
    - `rna-reads export-paths-tsv` and `rna-reads export-abundance-tsv` now
      accept the same `--record-indices` exact-subset override plus optional
      `--subset-spec` provenance.
    - `rna-reads export-alignment-dotplot-svg` emits a dotplot-style alignment
      scatter (coverage vs identity) with score-based point coloring.
    - `rna-reads export-hits-fasta` headers include seed metrics and exon-path
      annotations:
      - `--record-indices i,j,k` exports the exact saved-report subset and
        overrides coarse `--selection`
      - `--subset-spec TEXT` records the formal subset definition in the FASTA
        headers
      - `exon_path_tx=<transcript_id|none>`
      - `exon_path=<ordinal_path|none>` where `:` marks hash-confirmed adjacent
        exon transitions and `-` marks unconfirmed adjacency
      - `exon_transitions=<confirmed>/<total>`
      - `rc_applied=<true|false>` indicating automatic cDNA poly-T
        reverse-complement normalization was applied
      - `origin_class=<...>` with `origin_conf` / `strand_conf` summary values
        from deterministic origin classification
    - cDNA/direct-RNA seed-normalization semantics:
      - default (`--cdna-poly-t-flip`): reads with a T-rich 5' head are
        reverse-complement normalized before hashing/scoring
        - head detection tolerates small interruptions/mismatches in the
          poly-T stretch near the 5' end
      - direct RNA (`--no-cdna-poly-t-flip`): reads are scored as provided
      - `--poly-t-prefix-min-bp` controls the minimum T support threshold used
        by the 5' head detector for automatic cDNA flip (default `18`)
    - phase-1 seed-span policy:
      - all reads are hashed across their full span (no 3-window sampling)
      - `--seed-stride-bp` is the real phase-1 density knob: `1` means one
        ordered k-mer start per base, while larger values trade sensitivity for
        speed
    - sparse-origin options:
      - `--origin-mode` accepts `single_gene` (default) and
        `multi_gene_sparse`
      - `--target-gene` (repeatable) and `--roi-seed-capture` are persisted in
        the report payload for deterministic follow-up runs
      - `multi_gene_sparse` expands template indexing with local-annotation
        matches from `--target-gene`
      - `--roi-seed-capture` remains planned and is surfaced as a deterministic
        warning until implemented
    - `rna-reads list-reports` rows expose sparse-origin request provenance:
      `origin_mode`, `target_gene_count`, and `roi_seed_capture_enabled`
      - shell/CLI JSON output also includes `summary_rows[]` for quick
        human-readable triage (`mode`, `origin`, target count, ROI-capture
        flag, and read counters)
    - `rna-reads show-report` output includes a `summary` field with the same
      compact provenance framing, plus aligned/full-length percentages and
      compact read-length histogram summaries for:
      - all reads
      - seed-passed reads
      - aligned reads
      - full-length exact / near / strict subsets
    - tutorial reference (TP53 basis + multi-gene sparse mapping):
      `docs/tutorial/generated/chapters/12_tp53_multi_gene_sparse_mapping_online.md`
    - batch target-gene cohort tutorial:
      `docs/tutorial/rna_read_batch_gene_support_cli.md`
    - `rna-reads export-sample-sheet` includes sparse-origin provenance columns
      (`report_mode`, `origin_mode`, `target_gene_count`,
      `target_gene_ids_json`, `roi_seed_capture_enabled`),
      `mean_read_length_bp`, and `origin_class_counts_json` alongside
      exon/junction frequency JSON fields
    - when `rna-reads export-sample-sheet` is given one or more `--gene`
      arguments, each report row also carries target-gene cohort metrics:
      - accepted-target counts/fractions
      - fragment vs complete counts (`near|strict|exact`)
      - `gene_support_mean_assigned_read_length_bp`
      - JSON-serialized exon support, ordered exon-pair co-presence, and
        neighboring direct-transition support
    - report compaction and checkpoint options:
      - `--report-mode full` (default): persist retained top hits as ranked
      - `--report-mode seed_passed_only`: persist only retained hits that
        passed the seed gate (stream counters are unchanged)
      - `--checkpoint-path PATH --checkpoint-every-reads N`: persist
        deterministic resume snapshots during streaming
      - `--resume-from-checkpoint`: load checkpoint snapshot and continue from
        the saved processed-read index
    - phase-1 seed-pass gate policy:
      - `pass = raw_hit_fraction >= min_seed_hit_fraction AND weighted_hit_fraction >= min_weighted_seed_hit_fraction AND unique_matched_kmers >= min(min_unique_matched_kmers, tested_kmers) AND chain_consistency_fraction >= min_chain_consistency_fraction AND median_transcript_gap <= max_median_transcript_gap AND confirmed_transitions >= min_confirmed_transitions AND confirmed_transition_fraction >= min_transition_support_fraction`
      - `chain_consistency_fraction` is the fraction of matched seed observations
        explained by one coherent transcript-offset chain
      - `weighted_hit_fraction` is occurrence-normalized inside the scoped seed
        index: `sum(1/occurrence_count(seed_bits))/tested_kmers`
      - defaults:
        - `--kmer-len 10`
        - `--min-seed-hit-fraction 0.30`
        - `--min-weighted-seed-hit-fraction 0.05`
        - `--min-unique-matched-kmers 12`
        - `--min-chain-consistency-fraction 0.40`
        - `--max-median-transcript-gap 4.0`
        - `--min-confirmed-transitions 1`
        - `--min-transition-support-fraction 0.05`
    - RNA-read scope semantics:
      - `all_overlapping_both_strands`: index all overlapping transcripts on
        both strands
      - `target_group_any_strand`: index only target-group transcripts, both
        strands allowed
      - `all_overlapping_target_strand`: index all overlapping transcripts on
        target strand only
      - `target_group_target_strand`: index only target-group transcripts on
        target strand
      - strand-scoring note:
        both-strand modes score against the union of admitted `+/-` templates;
        target-strand modes exclude opposite-strand templates from scoring.
      - seed-index note:
        indexed seeds include exon-body and exon-exon transition k-mers for
        admitted transcripts.
    - Primer request payload notes (`primers design` / `primers design-qpcr`):
      - `forward`/`reverse`/`probe` side constraints support optional
        sequence filters:
        `non_annealing_5prime_tail`, `fixed_5prime`, `fixed_3prime`,
        `required_motifs[]`, `forbidden_motifs[]`, `locked_positions[]`.
      - `non_annealing_5prime_tail` is added to the oligo sequence while
        anneal `tm_c`/`gc_fraction`/`anneal_hits` remain computed on the
        template-binding segment only.
      - `pair_constraints` is optional and supports:
        `require_roi_flanking`, amplicon motif filters, and fixed amplicon
        start/end coordinates.
    - Primer ROI seed helper notes (`primers seed-from-feature` / `primers seed-from-splicing`):
      - returns non-mutating schema `gentle.primer_seed_request.v1`
      - includes `template`, source metadata, `roi_start_0based`,
        `roi_end_0based_exclusive`
      - includes ready-to-run operations:
        `operations.design_primer_pairs` (`DesignPrimerPairs`) and
        `operations.design_qpcr_assays` (`DesignQpcrAssays`)
    - qPCR-only seed helper notes
      (`primers seed-qpcr-from-feature` /
      `primers seed-qpcr-from-splicing`):
      - returns non-mutating schema `gentle.qpcr_seed_request.v1`
      - includes `template`, source metadata, `roi_start_0based`,
        `roi_end_0based_exclusive`
      - includes `rationale.summary`, `rationale.why_this_roi`, and
        `rationale.recommended_defaults`
      - includes one ready-to-run `operation` (`DesignQpcrAssays`)
      - includes built-in protocol-cartoon metadata for `pcr.assay.qpcr`, so
        shell/CLI/ClawBio flows can promote the same qPCR strip without hard-
        coding that protocol id elsewhere
    - Restriction-cloning handoff notes (`primers prepare-restriction-cloning`):
      - expects an `Operation` payload whose root variant is
        `PrepareRestrictionCloningPcrHandoff`
      - creates extended forward/reverse primer artifacts plus one predicted
        tailed amplicon artifact and one per-handoff container
      - persists report schema
        `gentle.restriction_cloning_pcr_handoff.v1`
      - shell/direct output includes the saved handoff report in the command
        response `output`
    - Restriction-cloning handoff seed helper
      (`primers seed-restriction-cloning-handoff PRIMER_REPORT_ID VECTOR_SEQ_ID ...`):
      - non-mutating helper returning schema
        `gentle.restriction_cloning_pcr_handoff_seed.v1`
      - resolves one saved primer-pair report + pair rank into a ready-to-run
        `PrepareRestrictionCloningPcrHandoff` operation payload
      - if enzymes are omitted, uses the top recommended single-site or
        directed-pair suggestion for the destination vector
      - if directed-pair enzymes are provided explicitly, validates that their
        order matches the same MCS/unique-cut ordering used by vector
        suggestions
    - Restriction-cloning vector suggestion notes
      (`primers restriction-cloning-vector-suggestions SEQ_ID`):
      - non-mutating helper for command-line/agent reasoning against one
        destination vector
      - returns schema
        `gentle.restriction_cloning_vector_enzyme_suggestions.v1`
      - reports `selected_mcs`, `other_unique`, `missing_mcs`,
        `recommended_single_site`, and `recommended_directed_pairs`
      - mirrors the GUI preference order:
        annotated `mcs_expected_sites` first, then other unique cutters
    - Restriction-cloning saved-report helpers:
      - `primers list-restriction-cloning-handoffs`
      - `primers show-restriction-cloning-handoff REPORT_ID`
      - `primers export-restriction-cloning-handoff REPORT_ID OUTPUT.json`
    - Feature query helper notes (`features query`):
      - non-mutating structured result schema:
        `gentle.sequence_feature_query_result.v1`
      - deterministic filters over feature kind, range relation, strand, labels,
        qualifiers, and length
      - deterministic ordering (`feature_id|start|end|kind|length`) with
        `offset`/`limit` paging suitable for agent iteration
    - Feature BED export helper notes (`features export-bed`):
      - non-mutating structured result schema:
        `gentle.sequence_feature_bed_export.v1`
      - writes BED6 plus deterministic extra columns:
        `kind`, `row_id`, `coordinate_source`, `qualifiers_json`
      - exports all matching rows by default when `--limit` is omitted; use
        `--offset`/`--limit` only when you want a bounded subset
      - `--coordinate-mode auto` prefers genomic BED coordinates when features
        carry `chromosome` + `genomic_start_1based` +
        `genomic_end_1based`; otherwise it falls back to local `SEQ_ID`
        coordinates
      - `--include-restriction-sites` appends deterministic REBASE-derived
        `restriction_site` rows, and `--restriction-enzyme NAME` can be
        repeated to keep only selected enzymes
    - Restriction-site scan helper notes (`features restriction-scan`):
      - non-mutating structured result schema:
        `gentle.restriction_site_scan.v1`
      - accepts either stored `SEQ_ID` or inline ASCII DNA via
        `--sequence-text`
      - `--range` / `--start` / `--end` restrict the scan to one local span on
        the chosen operand
      - when no `--enzyme` is supplied, the shared preferred restriction-enzyme
        list is used; if that list is empty, GENtle falls back to the default
        preferred enzyme set
      - the report includes both local scan coordinates and source-sequence
        coordinates, plus optional cleavage geometry unless
        `--no-cut-geometry` is set
    - `panels import-isoform SEQ_ID PANEL_PATH [--panel-id ID] [--strict]`
    - `panels inspect-isoform SEQ_ID PANEL_ID`
    - `panels render-isoform-svg SEQ_ID PANEL_ID OUTPUT.svg`
    - `panels validate-isoform PANEL_PATH [--panel-id ID]`
    - `set-param NAME JSON_VALUE`
      - anchor-related examples:
        - `set-param require_verified_genome_anchor_for_extension true`
        - `set-param genome_anchor_prepared_fallback_policy "single_compatible|always_explicit|off"`
    - promoter-design TF score-track examples:
      - structured summary through the raw operation bridge:
        - `gentle_cli op '{"type":"SummarizeTfbsScoreTracks","seq_id":"tp73_context","motifs":["TP73","SP1","BACH2","PATZ1"],"start_0based":0,"end_0based_exclusive":4000,"score_kind":"llr_quantile","clip_negative":true}'`
      - shared SVG export through the shell/CLI route:
        - `gentle_cli shell 'features tfbs-score-tracks-svg tp73_context docs/figures/tp73_upstream_tfbs_score_tracks.svg --motif TP53 --motif TP63 --motif TP73 --motif PATZ1 --motif SP1 --motif BACH2 --motif REST --range 15564..16764 --score-kind llr_quantile --allow-negative'`
      - `--score-kind` lets you inspect either raw bit scores or empirical
        quantiles; the quantile modes are often the better first pass when you
        want to compare motif prominence within one local window
      - clip-negative mode remains the intended default for presentation-oriented
        promoter plots on the bit-based score kinds because it suppresses
        negative-only windows and leaves only positive motif support;
        `--allow-negative` is available when a figure should show the full
        continuous score landscape even when only a subset of factors crosses
        into positive support
    - `op <operation-json-or-@file>`
    - `workflow <workflow-json-or-@file>`
    - `screenshot-window OUTPUT.png` (currently disabled by security policy)
    - `help [COMMAND ...] [--format text|json|markdown] [--interface all|cli-direct|cli-shell|gui-shell|js|lua|mcp]`
  - Use single quotes around JSON payloads to preserve whitespace:
    - `gentle_cli shell 'workflow {"run_id":"r1","ops":[]}'`
  - Structured help export for automation:
    - `gentle_cli help --format json`
    - `gentle_cli help --format markdown`

Screenshot bridge:

- Current status: disabled by security policy.
- `screenshot-window` currently returns a deterministic disabled-policy message.
- `--allow-screenshots` is rejected by argument parsing.
- Manual documentation updates remain the current approach for image artifacts.

Isoform architecture panel workflow:

- canonical panel resource:
  - `assets/panels/tp53_isoforms_v1.json`
- operation route (JSON):
  - `ImportIsoformPanel` + `RenderIsoformArchitectureSvg`
- shared shell route:
  - `panels import-isoform ...`
  - `panels inspect-isoform ...`
  - `panels render-isoform-svg ...`
  - `panels validate-isoform ...`
- direct expert route:
  - `inspect-feature-expert SEQ_ID isoform PANEL_ID`
  - `render-feature-expert-svg SEQ_ID isoform PANEL_ID OUTPUT.svg`
  - same command family for splicing:
    - `inspect-feature-expert SEQ_ID splicing FEATURE_ID`
    - `render-feature-expert-svg SEQ_ID splicing FEATURE_ID OUTPUT.svg`
  - same command family for transcript-first protein comparison, with or
    without stored external protein evidence:
    - `inspect-feature-expert SEQ_ID protein-comparison [--transcript TRANSCRIPT_ID] [--ensembl-entry ENTRY_ID] [--feature-key KEY]... [--feature-key-not KEY]...`
    - `render-feature-expert-svg SEQ_ID protein-comparison [--transcript TRANSCRIPT_ID] [--ensembl-entry ENTRY_ID] [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
    - this route opens the same transcript-first Protein Expert payload used by
      the GUI `Open Derived Protein Expert` action
    - `--ensembl-entry` layers one stored Ensembl protein entry onto that same
      transcript-first compare surface as optional external evidence without
      redefining which transcript-derived products exist
    - `--feature-key` / `--feature-key-not` can trim imported Ensembl or other
      external-opinion feature classes when the comparison view gets crowded
  - same command family for persisted UniProt protein mappings:
    - `inspect-feature-expert SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]...`
    - `render-feature-expert-svg SEQ_ID uniprot-projection PROJECTION_ID [--feature-key KEY]... [--feature-key-not KEY]... OUTPUT.svg`
    - `CONFLICT` is hidden by default; explicit include/exclude keys let you
      trim noisy UniProt classes such as `VARIANT`, `STRAND`, or `HELIX`
  - direct Ensembl protein metadata routes:
    - `ensembl-protein fetch QUERY [--entry-id ID]`
    - `ensembl-protein list`
    - `ensembl-protein show ENTRY_ID`
    - `ensembl-protein import-sequence ENTRY_ID [--output-id ID]`
    - fetch/import keep Ensembl as optional external evidence plus an ordinary
      first-class protein import path; they do not replace transcript-native
      derivation as the authoritative product model
  - direct coding-DNA query for one projected UniProt feature:
    - `uniprot feature-coding-dna PROJECTION_ID FEATURE_QUERY [--transcript ID] [--mode genomic_as_encoded|translation_speed_optimized|both] [--speed-profile human|mouse|yeast|ecoli]`
    - returns one structured report per matching transcript feature span, including:
      - amino-acid interval
      - genomic coding DNA as encoded by the genome (spliced, coding-strand orientation)
      - optional preferred-codon translation-speed DNA
      - exon attribution and ordered exon pairs when the feature crosses a splice junction
  - reusable UniProt/Ensembl audit primitives now also exist for external AI orchestration:
    - `uniprot resolve-ensembl-links ...`
    - `uniprot transcript-accounting ...`
    - `uniprot compare-ensembl-exons ...`
    - `uniprot compare-ensembl-peptide ...`
  - the integrated audit persists a human-facing report plus a local unsent maintainer-email draft:
    - `uniprot audit-projection ...`
    - `uniprot audit-parity ...`
    - `uniprot audit-list|show|export ...`
    - `uniprot audit-parity-list|show|export ...`

Pool exchange commands:

- `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
  - Exports explicit sequence IDs (`IDS` is comma-separated) with topology and
    overhang fields into a versioned JSON pool artifact.
  - topology now preserves richer gel-oriented circular-form hints when GENtle
    can infer them (`supercoiled`, `relaxed_circular`, `nicked_circular`)
  - Adds `human_id` at pool level and per-member `human_id`.
- `import-pool INPUT.pool.gentle.json [PREFIX]`
  - Imports pool members into current state; generated IDs are prefixed.
  - richer topology strings such as `supercoiled` or `nicked_circular` are
    preserved as gel-topology hints on imported circular molecules.
- `export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]`
  - Calls engine operation `ExportProcessRunBundle`.
  - Writes a deterministic process/audit bundle JSON
    (`gentle.process_run_bundle.v1`) with:
    - extracted operation inputs
    - chronological parameter overrides
    - selected operation log records
    - output summaries (created/changed ids, exported artifact paths)
  - Omitting `--run-id` exports all operation-log rows.

Rendering export commands:

- `render-svg SEQ_ID linear|circular OUTPUT.svg`
  - Calls engine operation `RenderSequenceSvg`.
  - Linear exports honor the current stored linear viewport when one is set,
    so workflow-driven `SetLinearViewport` crops carry through into the SVG.
  - Single-base `variation` features render as DNA-baseline markers, which is
    useful for dbSNP-driven context figures.
  - Linear exports now also add transcription-start hooked arrows for
    strand-bearing `gene`/`mRNA`/`CDS`/`promoter` features and suppress
    unlabeled fallback coordinate labels that would otherwise clutter figure
    exports.
  - When transcript/gene features only carry accession-like ids, exports now
    prefer gene-style labels when available and suppress redundant nearby
    duplicates for cleaner communication figures.
  - Direction-bearing `mRNA`/`promoter` bars now render with arrowed ends, and
    the linear transcription-start cue uses a short hooked arrow rather than
    only a flat triangle marker.
  - Circular exports now use a transparent canvas, show single-base
    `variation` features as radial DNA-ring markers, and add
    transcription-start arrow shafts/arrowheads for strand-bearing
    `gene`/`mRNA`/`CDS`/`promoter` features. Circular figure exports also use
    a slightly larger ring and larger label fonts for readability.
- `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor] [--overlay-anchor-exon START..END]`
  - Calls engine operation `RenderDotplotSvg`.
  - `DOTPLOT_ID` must exist in stored dotplot payloads (`dotplot compute ...` / GUI compute).
  - `SEQ_ID` is the owner sequence id of the stored payload; for reference-centered
    overlay payloads this is the shared reference/locus DNA sequence, not one
    particular query isoform.
  - `--flex-track` optionally overlays one stored flexibility track in the same SVG.
    Overlay dotplots suppress that panel because each query series uses its own
    normalized x-axis.
  - `--display-threshold` and `--intensity-gain` apply the same density/contrast controls as GUI dotplot display.
  - Overlay payload exports render all stored query series with a legend and a
    merged-exon side track when reference exon annotation is available.
  - `shared_exon_anchor` requires `--overlay-anchor-exon START..END` using the
    stored reference-exon coordinates from the overlay payload.
  - `query_anchor_bp` uses explicit `query_anchor_0based` values already stored
    in the overlay payload. This is the shared route for curated cross-family
    or domain-anchor figures; it does not require `--overlay-anchor-exon`.
  - Alias: `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]`.
- `dotplot overlay-compute OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
  - Calls engine operation `ComputeDotplotOverlay`.
  - `OWNER_SEQ_ID` identifies the stored payload owner and is also reused as
    the shared reference when `--reference-seq` is omitted.
  - Each `--query-spec` must deserialize to one `DotplotOverlayQuerySpec`
    record, for example:
    `{"seq_id":"tp53_iso_1","label":"Isoform 1","mode":"pair_forward","color_rgb":[29,78,216]}`
  - Query specs can optionally narrow one isoform/query span via
    `span_start_0based` / `span_end_0based` and switch per-series orientation
    with `mode=pair_forward|pair_reverse_complement`.
  - Cross-family/domain-anchor overlays can also carry
    `query_anchor_0based` plus optional `query_anchor_label`, for example:
    `{"seq_id":"tp63_family_query","label":"TP63","query_anchor_0based":1044,"query_anchor_label":"shared core motif CATGTGTAACAG","mode":"pair_forward","color_rgb":[220,38,38]}`
  - `--ref-start` / `--ref-end` constrain the shared reference span on the
    y-axis for every stored series.
- `render-rna-svg SEQ_ID OUTPUT.svg`
  - Calls engine operation `RenderRnaStructureSvg`.
  - Accepts only single-stranded RNA (`molecule_type` of `RNA`/`ssRNA`).
  - Uses external `rnapkin` executable (set `GENTLE_RNAPKIN_BIN` to override executable path).
- `render-lineage-svg OUTPUT.svg`
  - Calls engine operation `RenderLineageSvg`.
- `protocol-cartoon list`
  - Lists deterministic protocol-cartoon IDs currently supported by the engine.
- `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonSvg`.
  - Built-in protocol IDs currently include `gibson.two_fragment`,
    `gibson.single_insert_dual_junction`, `pcr.assay.pair`,
    `pcr.assay.pair.no_product`, `pcr.assay.pair.with_tail`,
    `pcr.oe.substitution`, and `pcr.assay.qpcr`.
- `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonTemplateSvg`.
  - Loads template JSON schema `gentle.protocol_cartoon_template.v1` and resolves
    sparse defaults before rendering.
- `protocol-cartoon template-validate TEMPLATE.json`
  - Calls engine operation `ValidateProtocolCartoonTemplate`.
  - Validates template JSON schema `gentle.protocol_cartoon_template.v1` and resolved cartoon semantics without writing SVG output.
- `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonTemplateWithBindingsSvg`.
  - Applies binding JSON schema `gentle.protocol_cartoon_template_bindings.v1` to the template before resolving and rendering.
- `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
  - Calls engine operation `ExportProtocolCartoonTemplateJson`.
  - Exports canonical built-in template JSON (`gentle.protocol_cartoon_template.v1`)
    for the requested protocol ID.
  - Protocol-cartoon commands intentionally use one canonical namespace with no
    extra alias routes, which keeps scripted and AI-facing command selection
    unambiguous.
- `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
  - Calls the shared non-mutating Gibson preview derivation path.
  - Accepts the single-insert destination-first subset of
    `gentle.gibson_assembly_plan.v1`.
  - Optional plan input:
    `validation_policy.desired_unique_restriction_site_enzyme_name`
    requests one new unique REBASE site on a terminal overlap when the current
    defined-site single-insert plan can support it.
  - The plan payload is resolved against the current project state, so the
    referenced `seq_id` values must already exist in memory or in the loaded
    state file.
  - Returns one `gentle.gibson_assembly_preview.v1` payload with:
    - resolved opening/junction overlaps
    - Gibson primer suggestions (`5' overlap + 3' priming`)
    - validation/advisory findings
    - structured suggested design adjustments for the current priming-window
      blocker when the overlap side already resolves
    - optional structured requested-unique-site outcome
    - protocol-cartoon bindings for the factual Gibson strip
- `gibson apply PLAN_JSON_OR_@FILE`
  - Calls engine operation `ApplyGibsonAssemblyPlan`.
  - Consumes the same plan JSON accepted by `gibson preview`.
  - Creates new project sequences for the Gibson outputs:
    - left insert primer
    - right insert primer
    - assembled product
  - Also records one serial arrangement for gel-oriented review:
    - original vector lane
    - ordered insert lane(s)
    - assembled product lane
    - recommended DNA ladders flanking the samples on export
  - Transfers destination and insert features onto the assembled product
    deterministically.
  - Partially consumed destination annotations are trimmed or projected when a
    truthful product-side rewrite is available.
  - Surviving MCS annotations are revalidated against actual restriction sites
    on the assembled product, including new sites introduced by the insert.
  - Returns the standard mutating `OpResult`, including created sequence IDs,
    warnings, messages, and the operation ID used for provenance/reopen flows.
- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
  - Calls the shared non-mutating prepared-cache inspection path.
  - Returns one `gentle.prepared_cache_inspection.v1` payload with:
    - selected cache roots
    - manifest-backed prepared installs
    - orphaned remnants under those roots
    - per-entry artifact groups, byte totals, and file counts
  - If no `--cache-dir` is given, default roots are:
    - `--references`: `data/genomes`
    - `--helpers`: `data/helper_genomes`
    - `--both`: both roots above
- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache [--references|--helpers|--both] [--cache-dir PATH ...] [--prepared-id ID ...] [--prepared-path PATH ...] [--include-orphans]`
  - Calls the shared prepared-cache cleanup path.
  - Returns one `gentle.prepared_cache_cleanup.v1` payload with exact affected
    paths, removed byte totals, and per-entry cleanup results.
  - `blast-db-only` and `derived-indexes-only` apply only to manifest-backed
    prepared installs and require one or more `--prepared-id` or
    `--prepared-path` selectors.
  - `selected-prepared` removes only the selected prepared installs and may
    also remove orphaned remnants when `--include-orphans` is set.
  - `--prepared-path` is the precise selector when duplicate prepared ids exist
    across multiple selected cache roots.
  - `all-prepared-in-cache` clears all prepared installs under the selected
    roots and ignores `--prepared-id` / `--prepared-path`.
  - Cleanup is conservative by design:
    - only selected known cache roots are touched
    - catalog JSON and project state files remain unchanged
    - `.gentle_state.json`, MCP/runtime files, backdrop/runtime caches, and
      `target/` are out of scope
- `render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID] [--agarose-pct FLOAT] [--buffer tae|tbe] [--topology-aware true|false]`
  - Calls engine operation `RenderPoolGelSvg`.
  - Use `IDS` as a comma-separated sequence-id list, or pass `-`/`_` when using `--containers` or `--arrangement`.
  - `--containers` renders one lane per container ID.
  - `--arrangement` renders lanes from a stored serial arrangement.
  - Optional `--ladders` overrides both auto mode and any saved arrangement
    ladder pair.
  - If `--ladders` is omitted for `--arrangement`, engine uses the saved
    arrangement ladder pair when present, otherwise falls back to auto ladder
    selection.
  - If `--ladders` is omitted without `--arrangement`, engine auto-selects one
    or two ladders based on pool bp range.
  - `--agarose-pct` / `--buffer` / `--topology-aware` override the shared
    deterministic gel-run profile for that render.
  - Defaults are `1.0%`, `TAE`, topology-aware `true`.
  - when topology-aware mode is on, explicit sequence hints like
    `supercoiled`, `relaxed circular`, or `nicked/open circular` refine
    apparent migration beyond the old generic circular-vs-linear split.
- `render-gel-svg ...`
  - Alias for `render-pool-gel-svg ...` with identical semantics.
- `arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]`
  - Calls engine operation `CreateArrangementSerial`.
  - Persists a serial lane setup that can be reused by `render-pool-gel-svg --arrangement`.
- `arrange-set-ladders ARR_ID [--ladders NAME[,NAME]]`
  - Calls engine operation `SetArrangementLadders`.
  - Persists one or two ladder names on an existing serial arrangement.
  - Omit `--ladders` to clear back to shared engine auto ladder selection.
- `containers set-exclusive CONTAINER_ID true|false`
  - Calls engine operation `SetContainerDeclaredContentsExclusive`.
  - `true` restores a clean-vial / declared-only interpretation.
  - `false` marks the listed members as a known/measured subset of a more
    complex sample.
- `racks create-from-arrangement ARR_ID [--rack-id ID] [--name TEXT] [--profile small_tube_4x6|plate_96|plate_384]`
  - Calls engine operation `CreateRackFromArrangement`.
  - Creates one deterministic physical rack draft from one stored arrangement.
  - If `--profile` is omitted, engine chooses the smallest built-in profile
    that fits the arrangement payload plus ladder-reference slots.
- `racks place-arrangement ARR_ID --rack RACK_ID`
  - Calls engine operation `PlaceArrangementOnRack`.
  - Appends the arrangement as one contiguous block onto the next free region
    in the target rack.
- `racks move RACK_ID --from A1 --to B1 [--block]`
  - Calls engine operation `MoveRackPlacement`.
  - Without `--block`, moves one sample within its arrangement block using
    shift-neighbor semantics.
  - With `--block`, moves the whole arrangement block and shifts later
    occupied positions in fill order.
- `racks move-samples RACK_ID --from A1 [--from A2 ...] --to B1`
  - Calls engine operation `MoveRackSamples`.
  - Moves multiple selected samples together as one contiguous group within a
    single arrangement block.
  - Shared engine preserves the selected samples in current rack order even if
    the CLI request lists them in another order.
- `racks move-blocks RACK_ID --arrangement ARR_ID [--arrangement ARR_ID ...] --to B1`
  - Calls engine operation `MoveRackArrangementBlocks`.
  - Moves the selected arrangement blocks together as one contiguous group.
  - Shared engine preserves the selected blocks in current rack order even if
    the CLI request lists them in another order.
- `racks show RACK_ID`
  - Returns `gentle.rack_state.v1` JSON with the saved rack record and resolved
    occupied positions.
- `racks labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--preset compact_cards|print_a4|wide_cards]`
  - Calls engine operation `ExportRackLabelsSvg`.
  - Writes one deterministic label sheet for the whole rack or one selected
    arrangement block on that rack.
  - Built-in presets:
    - `compact_cards` (current compact two-column cards)
    - `print_a4` (denser print-oriented A4 sheet)
    - `wide_cards` (single-column wider cards)
- `racks fabrication-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
  - Calls engine operation `ExportRackFabricationSvg`.
  - Writes one top-view fabrication/planning SVG for the full saved rack.
  - Current built-in physical templates:
    - `storage_pcr_tube_rack`
    - `pipetting_pcr_tube_rack`
- `racks isometric-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
  - Calls engine operation `ExportRackIsometricSvg`.
  - Writes one pseudo-3D/isometric SVG for the full saved rack.
  - Intended for presentation-ready rack review and README/documentation
    figures while still using the same saved rack/template state as the other
    physical exports.
- `racks openscad RACK_ID OUTPUT.scad [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
  - Calls engine operation `ExportRackOpenScad`.
  - Writes one parameterized OpenSCAD source file for the full saved rack.
- `racks carrier-labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack] [--preset front_strip_and_cards|front_strip_only|module_cards_only]`
  - Calls engine operation `ExportRackCarrierLabelsSvg`.
  - Writes one carrier-matched front-strip plus module-label SVG sheet for the
    full rack or one selected arrangement block on that rack.
  - Built-in carrier-label presets:
    - `front_strip_and_cards`
    - `front_strip_only`
    - `module_cards_only`
- `racks simulation-json RACK_ID OUTPUT.json [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
  - Calls engine operation `ExportRackSimulationJson`.
  - Writes one machine-readable rack geometry/placement JSON export for
    downstream simulation adapters.
- `racks set-profile RACK_ID small_tube_4x6|plate_96|plate_384`
  - Calls engine operation `SetRackProfile`.
  - Reflows occupied positions onto another built-in rack/plate profile while
    preserving arrangement order.
- `racks apply-template RACK_ID bench_rows|plate_columns|plate_edge_avoidance`
  - Calls engine operation `ApplyRackTemplate`.
  - Applies one shared authoring shortcut on top of the existing rack snapshot.
  - Built-in templates:
    - `bench_rows`
      - row-major fill
      - clears blocked coordinates
    - `plate_columns`
      - column-major fill
      - clears blocked coordinates
    - `plate_edge_avoidance`
      - column-major fill
      - blocks the outer perimeter
- `racks set-fill-direction RACK_ID row_major|column_major`
  - Calls engine operation `SetRackFillDirection`.
  - Reflows occupied positions onto the same rack geometry under a different
    deterministic fill order.
- `racks set-custom-profile RACK_ID ROWS COLUMNS`
  - Calls engine operation `SetRackProfileCustom`.
  - Reflows occupied positions onto one custom A1-style rack geometry.
  - A1-style row labels continue beyond `Z` as `AA`, `AB`, ...
- `racks set-blocked RACK_ID COORD[,COORD...]...`
  - Calls engine operation `SetRackBlockedCoordinates`.
  - Reserves one normalized blocked-coordinate set on the rack and reflows
    occupied positions across the remaining available slots.
- `racks set-blocked RACK_ID --clear`
  - Clears all blocked rack coordinates.

RNA secondary-structure text command:

- `rna-info SEQ_ID`
  - Returns a JSON report from `rnapkin -v -p <sequence>`.
  - Includes `stdout`/`stderr` and command metadata.
  - Accepts only single-stranded RNA (`molecule_type` of `RNA`/`ssRNA`).
  - Uses external `rnapkin` executable (set `GENTLE_RNAPKIN_BIN` to override executable path).

DNA ladder catalog commands:

- `ladders list [--molecule dna|rna] [--filter TEXT]`
  - Returns engine-inspected ladder catalog (`schema`, `ladder_count`, `ladders`).
  - Default molecule is `dna`; use `--molecule rna` for RNA ladders.
  - Optional `--filter` keeps only ladders whose names contain `TEXT` (case-insensitive).
- `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
  - Calls engine operation `ExportDnaLadders` or `ExportRnaLadders`.
  - Writes structured ladder catalog JSON to disk.

Resource sync commands:

- `services status`
  - Reports one combined service-readiness snapshot for common GENtle-backed
    analysis flows.
  - Covers canonical reference genomes (currently `Human GRCh38 Ensembl 116`),
    helper backbones (currently `Plasmid pUC19 (online)`), and active external
    resource snapshots in one record.
  - When a prepare/reindex run is active, it also surfaces live dependency
    activity (`running`, `failed`, `cancelled`) plus the current prepare phase
    and progress hints instead of only saying `prepared`/`not_prepared`.
  - Includes compact `summary_lines` that chat/report layers such as
    ClawBio/OpenClaw can surface directly before deciding whether they should
    fetch, prepare, or render.
- `resources status`
  - Reports which integrated external resource snapshots are active right now.
  - Covers built-in/runtime status for `REBASE` and `JASPAR`, including item
    counts and whether GENtle is currently using a runtime override or the
    built-in snapshot.
  - JASPAR status now also reports the persisted remote-metadata snapshot used
    by the JASPAR catalog/expert cache:
    - `path`
    - `exists`
    - `valid`
    - `item_count`
    - `error`
  - Also reports normalized ATtRACT runtime status:
    - `known_external_only` = homepage/ZIP URL known, but no active normalized
      snapshot
    - `runtime_snapshot` = `data/resources/attract.motifs.json` is present and
      active
    - `session_override` = the current session has loaded an ATtRACT snapshot
      from a non-default path
  - ATtRACT status now also exposes deterministic provenance for the active and
    runtime snapshots:
    - `runtime_fingerprint`
    - `active_fingerprint`
  - Includes the current ATtRACT homepage plus published ZIP download URL:
    `https://attract.cnic.es/attract/static/ATtRACT.zip`.
- `resources sync-rebase INPUT.withrefm [OUTPUT.rebase.json] [--commercial-only]`
  - Parses REBASE/Bairoch-style records (`withrefm`) into GENtle restriction-enzyme JSON.
  - `INPUT` may be a local file path or an `https://...` URL.
  - Default output: `data/resources/rebase.enzymes.json`.
  - At runtime, GENtle auto-loads this file (if present) and overrides embedded restriction enzymes.
- `resources sync-jaspar INPUT.jaspar.txt [OUTPUT.motifs.json]`
  - Parses JASPAR PFM text into motif snapshot JSON with IUPAC consensus.
  - `INPUT` may be a local file path or an `https://...` URL.
  - Default output: `data/resources/jaspar.motifs.json`.
  - `ExtractAnchoredRegion` can resolve TF motifs by ID/name from this local registry.
  - GENtle ships with built-in motifs in `assets/jaspar.motifs.json` (currently generated from JASPAR 2026 CORE non-redundant); this command provides local updates/extensions.
- `resources summarize-jaspar [--motif TOKEN ...] [--motifs CSV] [--all] [--random-length N] [--seed N] [--output OUTPUT.json]`
  - Presents local JASPAR entries through one deterministic score-oriented report.
  - For each selected entry, GENtle derives:
    - a maximizing sequence,
    - a minimizing sequence,
    - and `llr_bits` / `true_log_odds_bits` score distributions on one deterministic pseudorandom DNA background.
  - Omitting `--motif` and `--motifs` summarizes the full active local JASPAR registry.
  - Default background length: `10000 bp`.
  - Default seed: one built-in deterministic seed so agent/CLI runs remain reproducible.
- `resources benchmark-jaspar [--random-length N] [--seed N] [--output OUTPUT.json]`
  - Benchmarks the full active local JASPAR registry through one deterministic
    shared background and writes an export-ready drift snapshot.
  - Default output: `data/resources/jaspar.registry_benchmark.json`.
  - Includes:
    - the full per-entry presentation rows used by `summarize-jaspar`
    - one score-family aggregate summary for `llr_bits`
    - one score-family aggregate summary for `true_log_odds_bits`
    - top motifs by maximum score and by positive random-background hit rate
  - This is the intended routine artifact for regression/drift checks across
    JASPAR updates or scoring changes.
- `resources list-jaspar [--filter TOKEN] [--limit N] [--fetch-remote] [--output OUTPUT.json]`
  - Lists local JASPAR entries as a portable catalog view.
  - Includes per-row:
    - motif id
    - optional TF name
    - consensus IUPAC sequence
    - motif length
    - and, when `--fetch-remote` is requested, a compact remote summary with
      collection/class/family/tax-group fields plus species-count previews for
      the returned subset
  - `--filter` matches motif id, TF name, and consensus.
  - `--limit` caps the returned subset before optional remote lookups.
  - When `--fetch-remote` is requested, GENtle first reuses any matching rows
    from `data/resources/jaspar.remote_metadata.json` and only refreshes
    missing/selected rows live when necessary.
- `resources sync-jaspar-remote-metadata [--motif TOKEN ...] [--motifs CSV] [--all] [--filter TOKEN] [--limit N] [--output OUTPUT.json]`
  - Refreshes and persists reusable JASPAR remote metadata for one selected
    subset of the local registry.
  - Default output: `data/resources/jaspar.remote_metadata.json`.
  - Stored rows include motif id/name, consensus, motif length, full remote
    metadata payload, and a compact summary for catalog views.
- `resources sync-attract INPUT.zip_or_URL [OUTPUT.attract.json]`
  - Normalizes the published ATtRACT ZIP into GENtle’s runtime motif snapshot.
  - `INPUT` may be a local ZIP path or an `https://...` URL.
  - Default output: `data/resources/attract.motifs.json`.
  - Current v1 normalization reads `ATtRACT_db.txt`, keeps deterministic
    consensus/IUPAC motif records, and when `pwm.txt` blocks can be mapped by
    `Matrix_id`, stores those PFM rows so splice-aware inspection can add
    PWM-backed ranking on top of the same candidate hits.
- `attract inspect-splicing SEQ_ID FEATURE_ID [--scope ...] [--organism NAME] [--flank-bp N] [--min-score X] [--min-match-quantile Q] [--pwm-mapping strict_same_length|windowed_submatrix] [--compare-policies] [--all-transcripts] [--no-fallback]`
  - Runs the shared engine-owned splice-aware ATtRACT inspection route over the
    selected Splicing Expert group.
  - Default behavior:
    - transcript-strand only
    - exact organism match first
      - optional fallback to all normalized ATtRACT entries when no exact match
      exists
      - region classification returned as `exon_body`, `donor_flank`,
      `acceptor_flank`, or `intron_body`
    - `--min-match-quantile` defaults to `0.99` and applies to PWM-backed rows;
      rows without a mapped PWM block keep exact consensus/IUPAC matching
  - Returns grouped RBP summary rows plus detailed hit rows from the same
    payload the GUI uses.
  - Payload provenance now also includes `active_resource_fingerprint` so CLI,
    GUI, and exported reports can tie the evidence back to the exact normalized
    ATtRACT snapshot content rather than only a motif count.
- `resources inspect-jaspar MOTIF [--random-length N] [--seed N] [--fetch-remote] [--output OUTPUT.json]`
  - Expands one local JASPAR entry into an expert-oriented report.
  - Includes:
    - the full count matrix,
    - a sequence-logo payload,
    - `llr_bits` and `true_log_odds_bits` score distributions with compact
      histogram bins,
    - and optional remote JASPAR metadata/species assignments when
      `--fetch-remote` is requested.
  - Default background length: `10000 bp`.
  - Default seed: one built-in deterministic seed so repeated expert
    inspection stays reproducible.
  - When `--fetch-remote` is requested, GENtle reuses any matching persisted
    remote metadata row first and only then falls back to a live refresh.

Agent bridge commands:

Conceptual/tutorial companion:

- `docs/agent_interfaces_tutorial.md` (role map + interface comparison:
  CLI/shared shell vs MCP vs Agent Assistant vs external coding agents).

- `agents list [--catalog PATH]`
  - Lists configured agent systems from catalog JSON.
  - Default catalog path: `assets/agent_systems.json`.
  - Includes availability status (`available`, `availability_reason`) so callers
    can skip systems that are not currently runnable.
- `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
  - Invokes one configured agent system and returns message/questions/suggested shell commands.
  - `--base-url` sets a per-request runtime endpoint override (maps to
    `GENTLE_AGENT_BASE_URL`) for native transports.
  - `--model` sets a per-request model override (maps to `GENTLE_AGENT_MODEL`)
    for native transports.
  - `--timeout-secs` sets per-request timeout override (maps to
    `GENTLE_AGENT_TIMEOUT_SECS`) for stdio/native transports.
  - `--connect-timeout-secs` sets per-request HTTP connect timeout (maps to
    `GENTLE_AGENT_CONNECT_TIMEOUT_SECS`) for native transports.
  - `--read-timeout-secs` sets per-request read timeout (maps to
    `GENTLE_AGENT_READ_TIMEOUT_SECS`) for stdio/native transports.
  - `--max-retries` sets per-request transient retry budget (maps to
    `GENTLE_AGENT_MAX_RETRIES`; `0` means no retries).
  - `--max-response-bytes` caps adapter response size per attempt (maps to
    `GENTLE_AGENT_MAX_RESPONSE_BYTES`).
  - for `native_openai_compat`, model must resolve to a concrete value (catalog
    model or `--model`); `unspecified` is rejected.
  - `--no-state-summary` disables project-context injection in the request.
  - External adapter responses must be strict `gentle.agent_response.v1` JSON.
  - Unknown canonical response fields are rejected (extensions must use `x_` or `x-` prefix).
  - Adapter calls retry transient failures with exponential backoff before returning an error.
  - Suggested commands are executed only when explicitly selected
    (`--execute-all`, `--execute-index`) or when `--allow-auto-exec` is enabled
    and the suggestion intent is `auto`.
  - Recursive `agents ask` execution from suggested commands is blocked by design.
  - Failures use deterministic error prefixes for scripting, e.g.
    `AGENT_ADAPTER_UNAVAILABLE`, `AGENT_ADAPTER_TRANSIENT`,
    `AGENT_RESPONSE_PARSE`, `AGENT_RESPONSE_VALIDATION`.
  - Supported transports:
    - `builtin_echo` (offline/demo)
    - `external_json_stdio` (local bridge executable)
    - `native_openai` (built-in OpenAI HTTP adapter, requires `OPENAI_API_KEY`)
    - `native_openai_compat` (built-in local OpenAI-compatible adapter for
      Jan/Msty/Ollama-style `/chat/completions` endpoints; key optional)
      - endpoint host/port come from catalog `base_url` (or `--base-url` if set);
        GENtle does not silently switch to a different host/port

Prompting users should send to agents:

- Prefer compact structured prompts so local models stay reliable.
- Use this template:

```text
Objective:
<one clear goal>

Context:
<sequence/genome/helper ids>

Inputs:
- anchors/coordinates:
- feature labels/kinds:

Constraints:
- length:
- GC range:
- motifs/sites:
- strand assumptions:

Output wanted:
- exact commands
- validation checklist

Execution policy:
chat-only | ask-before-run | allow-auto-exec
```

- For local models without cloning background, preload:
  - `docs/ai_cloning_primer.md`
  - `docs/ai_task_playbooks.md`
  - `docs/ai_prompt_contract.md`
  - `docs/examples/ai_cloning_examples.md`
  - optional: `docs/ai_glossary_extensions.json`

Genome convenience commands:

- `genomes list [--catalog PATH]`
  - Lists available genomes in the catalog.
- `genomes ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
  - Lists Ensembl species directories that currently look installable because both FASTA and GTF listings are present.
  - Returns a read-only discovery report with current listing URLs and latest release numbers seen per collection.
- `genomes install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
  - Resolves concrete current Ensembl FASTA/GTF URLs for one species directory, writes a real catalog entry, and immediately runs the normal prepare workflow.
  - When the active catalog is a writable single JSON file, GENtle updates that file in place.
  - When the active catalog comes from discovery/overlay mode, GENtle writes only the new entry into an overlay fragment instead of copying the whole merged catalog.
  - `--output-catalog PATH` can point to a catalog JSON file or to a directory where a fragment file should be created.
  - `--genome-id ID` overrides the default deterministic id derived from species + assembly + release.
- `genomes validate-catalog [--catalog PATH]`
  - Verifies catalog JSON schema/entry rules and that each entry resolves usable
    sequence/annotation source definitions.
- `genomes update-ensembl-specs [--catalog PATH] [--output-catalog PATH]`
  - Refreshes explicit pinned Ensembl URLs/specs for catalog rows that carry
    `ensembl_template` metadata.
  - Older pinned release rows stay in the catalog.
  - If a newer current species release exists, GENtle adds or refreshes the
    newest pinned row rather than overwriting the old one.
  - If the active catalog file is not writable, `--output-catalog PATH` is
    required so the updated catalog can be written as a copy.
- `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
  - Shows whether the genome cache is prepared/indexed.
  - Also reports `effective_cache_dir`, `requested_catalog_key`,
    `requested_family`, and `compatible_prepared_options` so "not prepared"
    states explain exactly which cache root was inspected and whether a
    same-family prepared install already exists.
  - When the selected cache root does not contain the requested genome, output
    also includes a ready-to-run `prepare_command` hint.
  - Also reports resolved source details: `sequence_source_type`,
    `annotation_source_type`, `sequence_source`, `annotation_source`.
  - Also reports optional mass metadata:
    `nucleotide_length_bp`, `molecular_mass_da`, `molecular_mass_source`.
- `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Lists indexed genes from prepared cache (paged by `--limit`/`--offset`).
  - `--filter` is a case-insensitive regular expression.
  - `--biotype` can be repeated to constrain matches to selected biotypes.
- `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
  - Runs engine `PrepareGenome`.
  - Output now includes `binary_preflight` (`gentle.blast_external_binary_preflight.v1`)
    with `makeblastdb`/`blastn` found/missing/version/path diagnostics captured
    before prepare work starts.
- `genomes remove-prepared GENOME_ID [--catalog PATH] [--cache-dir PATH]`
  - Deletes only the prepared install directory for one genome from the selected cache.
  - Catalog JSON is left unchanged.
- `genomes remove-catalog-entry GENOME_ID [--catalog PATH] [--output-catalog PATH]`
  - Deletes only one row from the selected genome catalog JSON.
  - Prepared cache files are left unchanged until `remove-prepared` is run explicitly.
  - If the active catalog file is not writable, `--output-catalog PATH` is
    required so the edited catalog can be written as a copy.
- `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Runs `blastn` against prepared genome cache/index.
  - `--task` defaults to `blastn-short`; accepted values: `blastn-short`, `blastn`.
  - `--options-json` / `--options-file` accept a JSON object that can override quick options and include thresholds (`max_evalue`, `min_identity_percent`, `min_query_coverage_percent`, `min_alignment_length_bp`, `min_bit_score`, `unique_best_hit`).
  - Output includes `binary_preflight` with explicit BLAST tool diagnostics.
- `genomes blast-start GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Starts one async BLAST job and returns a stable `job_id`.
  - Jobs run through a bounded FIFO scheduler and may return initial state
    `queued` or `running` depending on current slot availability.
  - Start payload includes `binary_preflight` with explicit BLAST tool diagnostics.
- `genomes blast-status JOB_ID [--with-report]`
  - Polls async BLAST job status; `--with-report` includes final report payload when available.
  - Status includes scheduler metadata (`max_concurrent_jobs`, `running_jobs`,
    `queued_jobs`, and `queue_position` while queued).
- `genomes blast-cancel JOB_ID`
  - Requests cooperative cancellation for one async BLAST job.
- `genomes blast-list`
  - Lists known async genome-BLAST jobs in the current process.
- `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeRegion`.
  - `--annotation-scope` selects projection policy (`none`, `core`, `full`).
  - Default scope is `core` when neither scope nor legacy flags are set.
  - `--max-annotation-features N` applies a deterministic safety cap
    (`0` disables the cap for explicit unrestricted transfer).
  - legacy `--include-genomic-annotation` maps to `--annotation-scope core`.
  - legacy `--no-include-genomic-annotation` maps to `--annotation-scope none`.
  - Result payload includes `genome_annotation_projection` telemetry.
- `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeGene`.
- `genomes extract-promoter GENOME_ID QUERY [--occurrence N] [--transcript-id ID] [--output-id ID] [--upstream-bp N] [--downstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomePromoterSlice`.
  - Derives one unclipped promoter slice directly from transcript TSS geometry.
  - When `--transcript-id` is omitted, GENtle deterministically uses the
    outermost 5' transcript among the matched gene’s transcript records and
    warns when multiple transcript candidates existed.
  - Defaults are `--upstream-bp 1000` and `--downstream-bp 200`.
- `genomes extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Runs engine `ExtendGenomeAnchor`.
  - Extends an already genome-anchored sequence in-silico on contextual `5'` or `3'`.
  - If exact anchor genome id is not prepared but one compatible assembly-family
    cache exists, extension auto-uses that cache and emits a warning.
  - If multiple compatible prepared caches exist, command fails and lists
    explicit options.
- `genomes verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Runs engine `VerifyGenomeAnchor`.
  - Re-checks one anchored sequence against prepared reference sequence at
    recorded coordinates and records `anchor_verified` in provenance.

Helper convenience commands:

- `helpers list [--catalog PATH]`
  - Same behavior as `genomes list`, but defaults to `assets/helper_genomes.json`.
  - Helper `entries[]` rows now also include an optional normalized
    `interpretation` record so ClawBio/agents/planners can consume one shared
    helper-meaning layer instead of reparsing raw catalog prose.
- `helpers ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
  - Same discovery report shape as `genomes ensembl-available`, but exposed under the helper-family command tree for contract symmetry across adapters.
- `helpers install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
  - Same quick-install contract as `genomes install-ensembl`, but targeting the helper-catalog/cache defaults.
- `helpers validate-catalog [--catalog PATH]`
  - Same behavior as `genomes validate-catalog`, with helper-catalog default.
- `helpers update-ensembl-specs [--catalog PATH] [--output-catalog PATH]`
  - Same behavior as `genomes update-ensembl-specs`, with helper-catalog default.
- `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes status`, with helper-catalog default
    (including length/mass metadata fields).

Host convenience commands:

- `hosts list [--catalog PATH] [--filter TEXT]`
  - Lists structured host/strain profile rows used by construct-reasoning host-fit logic.
  - Default catalog is the bundled starter catalog at `assets/host_profiles.json`.
  - Output now also includes the same optional normalized `interpretation`
    record for helper rows that carry structured helper semantics.
- `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Same behavior as `genomes genes`, with helper-catalog default.
- `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
  - Same behavior as `genomes prepare`, with helper-catalog default.
  - `--timeout-secs N`: optional prepare-job timebox.
- `helpers remove-prepared HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes remove-prepared`, with helper-catalog default.
- `helpers extract-promoter HELPER_ID QUERY [--occurrence N] [--transcript-id ID] [--output-id ID] [--upstream-bp N] [--downstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-promoter`, with helper-catalog default.
- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
  - Lists prepared installs and orphaned remnants under the selected cache
    roots with artifact-group byte totals.
  - Default roots are `data/genomes` for references and
    `data/helper_genomes` for helpers.
  - Runtime defaults can be redirected for sibling worktrees via:
    - `GENTLE_REFERENCE_CACHE_DIR`
    - `GENTLE_HELPER_CACHE_DIR`
- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache [--references|--helpers|--both] [--cache-dir PATH ...] [--prepared-id ID ...] [--prepared-path PATH ...] [--include-orphans]`
  - Conservatively deletes derived cache artifacts or prepared installs inside
    the selected roots only.
  - Partial modes keep cached FASTA/annotation sources and manifests intact so
    reindex-from-cached-files remains possible.
  - Catalog JSON, `.gentle_state.json`, MCP/runtime files, backdrop/runtime
    caches, and `target/` are not treated as cache.
- `racks create-from-arrangement ARR_ID [--rack-id ID] [--name TEXT] [--profile small_tube_4x6|plate_96|plate_384]`
  - Same shared rack-draft creation path used by GUI `Open Rack`.
- `racks place-arrangement ARR_ID --rack RACK_ID`
  - Same shared rack-append path used by GUI `Place on Existing Rack...`.
- `racks move RACK_ID --from A1 --to B1 [--block]`
  - Same shared rack move/reorder contract used by the rack editor.
- `racks move-samples RACK_ID --from A1 [--from A2 ...] --to B1`
  - Same shared multi-sample rack move contract used by rack-editor
    `Command`/`Ctrl`-selected sample groups.
- `racks move-blocks RACK_ID --arrangement ARR_ID [--arrangement ARR_ID ...] --to B1`
  - Same shared multi-block rack move contract used by rack-editor
    Command/Ctrl-selected arrangement chips.
- `racks show RACK_ID`
  - Emits structured rack placement JSON (`gentle.rack_state.v1`).
  - Rack profile payload includes current `fill_direction` and
    `blocked_coordinates`.
- `racks labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--preset compact_cards|print_a4|wide_cards]`
  - Exports deterministic rack or arrangement-scoped label SVG through the
    same engine path as GUI `Labels SVG`.
- `racks fabrication-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
  - Exports one deterministic top-view fabrication/planning SVG through the
    same engine path as GUI `Fabrication SVG...`.
- `racks isometric-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
  - Exports one deterministic pseudo-3D/isometric rack SVG through the same
    engine path as GUI `Isometric SVG...`.
- `racks openscad RACK_ID OUTPUT.scad [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack]`
  - Exports one deterministic parameterized OpenSCAD file through the same
    engine path as GUI `OpenSCAD...`.
- `racks set-profile RACK_ID small_tube_4x6|plate_96|plate_384`
  - Reflows one saved rack onto another built-in profile.
- `racks set-custom-profile RACK_ID ROWS COLUMNS`
  - Reflows one saved rack onto a custom A1-style geometry while preserving
    arrangement order.
- `helpers remove-catalog-entry HELPER_ID [--catalog PATH] [--output-catalog PATH]`
  - Same behavior as `genomes remove-catalog-entry`, with helper-catalog default.
- `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast`, with helper-catalog default.
- `helpers blast-start HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast-start`, with helper-catalog default.
- `helpers blast-status JOB_ID [--with-report]`
  - Same behavior as `genomes blast-status`, scoped to helper jobs.
- `helpers blast-cancel JOB_ID`
  - Same behavior as `genomes blast-cancel`, scoped to helper jobs.
- `helpers blast-list`
  - Same behavior as `genomes blast-list`, scoped to helper jobs.
- `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-region`, with helper-catalog default.
  - For helper IDs containing `pUC18`/`pUC19`, GENtle auto-attaches a
    canonical MCS `misc_feature` when no MCS annotation is present in source
    annotation and exactly one canonical MCS motif is found.
  - MCS features expose `mcs_expected_sites` with REBASE-normalized enzyme names
    when recognizable from source/fallback annotation text.
- `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-gene`, with helper-catalog default.
  - pUC18/pUC19 helper extractions apply the same automatic MCS fallback
    annotation behavior when applicable (non-unique motif matches are warned and
    skipped).
- `helpers extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Same behavior as `genomes extend-anchor`, with helper-catalog default.
- `helpers verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Same behavior as `genomes verify-anchor`, with helper-catalog default.

Workflow macro commands (`gentle_cli shell 'macros ...'`):

- `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
  - Executes semicolon/newline-separated shell statements.
  - Existing file paths are auto-loaded even without `@` (shebang-friendly).
  - Supports transactional rollback (`--transactional`) when any statement fails.
  - Designed for full cloning workflows through `op ...` and `workflow ...`
    statements (Digest/Ligation/PCR/ExtractRegion/container ops, etc.).
  - All runs persist a lineage macro-instance record:
    - success: status `ok`
    - failure: status `failed` (or `cancelled` when cancellation-like error text is detected)
  - Successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in the error message.
- `macros instance-list`
  - Lists recorded macro-instance lineage rows.
  - Response schema: `gentle.lineage_macro_instances.v1`.
- `macros instance-show MACRO_INSTANCE_ID`
  - Shows one recorded macro-instance lineage row.
  - Response schema: `gentle.lineage_macro_instance.v1`.
- `macros template-list`
  - Lists persisted workflow macro templates.
- `macros template-show TEMPLATE_NAME`
  - Shows one persisted workflow template definition.
- `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port PORT_ID:KIND[:one|many][:required|optional][:description]]`
  - Creates/updates a named workflow macro template in project metadata.
  - Placeholders in script use `${param_name}` and must be declared via `--param`.
  - Optional `--details-url URL` records external protocol/reference details for
    catalog display.
  - Optional typed port contracts (`--input-port` / `--output-port`) are stored
    in template metadata and used as first-class preflight contract source.
- `macros template-delete TEMPLATE_NAME`
  - Deletes one persisted workflow template.
- `macros template-import PATH`
  - Imports workflow macro templates from:
    - one pack JSON file (`gentle.cloning_patterns.v1`)
    - one single-template JSON file (`gentle.cloning_pattern_template.v1`)
    - one directory tree (recursive `*.json` import)
  - If one template fails validation, no imported template changes are kept.
- `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]`
  - Expands a named template with provided bindings/defaults, then executes it as
    a workflow macro script.
  - Executes typed preflight checks before mutation using template port
    contracts when present, otherwise routine catalog mapping.
  - Preflight now includes richer semantics:
    - cross-port alias/collision checks,
    - input sequence vs input container consistency checks,
    - sequence-anchor semantic checks against the bound input sequence when unambiguous,
    - Gibson-family overlap checks (adjacent suffix/prefix validation against
      configured overlap length on Gibson routines),
    - restriction-family digest checks (enzyme-name resolution, duplicate-enzyme
      misuse, recognition-site presence across bound input sequences, and
      digest parameter sanity for fragment indices/extract range).
  - `--validate-only` runs preflight only and never mutates state.
  - Mutating runs now always record a macro-instance lineage row:
    - success: `ok`
    - preflight/execute failure: `failed`/`cancelled`
  - Successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in the error message.

Typed routine catalog command (`gentle_cli routines ...` or `gentle_cli shell 'routines ...'`):

- `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT] [--seq-id SEQ_ID]`
  - Lists typed cloning routines from catalog JSON (`gentle.cloning_routines.v1`).
  - `--family`, `--status`, `--tag`: exact case-insensitive filters.
  - `--query`: case-insensitive substring match across id/title/family/status/template/tags/summary plus explainability metadata fields.
  - `--seq-id`: refreshes construct reasoning for that sequence and feeds any
    variant-derived assay preferences into the planning estimate / ranking.
  - Default catalog path: `assets/cloning_routines.json`.
- `routines explain ROUTINE_ID [--catalog PATH] [--seq-id SEQ_ID]`
  - Returns structured explainability payload for one routine.
  - Response schema: `gentle.cloning_routine_explain.v1`.
  - Includes resolved alternatives plus purpose/mechanism/requirements, contraindications, disambiguation questions, and failure modes.
  - `--seq-id` adds the same construct-reasoning planning context used by
    `routines list` / `routines compare`, including any variant-derived assay
    preferences and the resulting planning estimate for the explained routine.
- `routines compare ROUTINE_A ROUTINE_B [--catalog PATH] [--seq-id SEQ_ID]`
  - Returns deterministic side-by-side comparison payload for two routines.
  - Response schema: `gentle.cloning_routine_compare.v1`.
  - Includes shared/unique vocabulary tags, difference-matrix rows, and merged disambiguation questions.
  - Also includes planning-aware estimate rows:
    `estimated_time_hours`, `estimated_cost`, `local_fit_score`,
    `composite_meta_score`.
  - `--seq-id` applies the same sequence-aware construct-reasoning preference
    context used by `routines list`.

Planning meta-layer commands (`gentle_cli planning ...` or `gentle_cli shell 'planning ...'`):

- Planning schemas:
  - `gentle.planning_profile.v1`
  - `gentle.planning_objective.v1`
  - `gentle.planning_estimate.v1`
  - `gentle.planning_suggestion.v1`
  - `gentle.planning_sync_status.v1`
- Effective profile merge precedence:
  - `global_profile -> confirmed_agent_overlay -> project_override`
- Schema compatibility:
  - payloads declaring a non-matching planning schema id are rejected
    (not auto-upgraded silently).
- `planning profile show [--scope global|project_override|confirmed_agent_overlay|effective]`
  - Shows one profile scope or merged effective profile.
- `planning profile set JSON_OR_@FILE [--scope global|project_override|confirmed_agent_overlay]`
  - Sets/replaces one profile scope.
- `planning profile clear [--scope global|project_override|confirmed_agent_overlay]`
  - Clears one profile scope.
- `planning objective show`
  - Shows current objective weights/guardrails.
  - Objective payloads may also include optional `helper_profile_id` and
    `preferred_routine_families[]` fields for helper-aware routine ranking.
- `planning objective set JSON_OR_@FILE`
  - Sets/replaces current planning objective.
- `planning objective clear`
  - Clears objective and falls back to engine defaults.
- `planning suggestions list [--status pending|accepted|rejected]`
  - Lists planning sync suggestions.
- `planning suggestions accept SUGGESTION_ID`
  - Accepts a pending suggestion and applies its patches.
- `planning suggestions reject SUGGESTION_ID [--reason TEXT]`
  - Rejects a pending suggestion.
- `planning sync status`
  - Shows sync metadata and pending count.
- `planning sync pull JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - Registers a pending inbound advisory suggestion.
- `planning sync push JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - Registers a pending outbound advisory suggestion.
- Sync payload shape (`pull`/`push`):
  - optional `profile_patch` (`gentle.planning_profile.v1`)
  - optional `objective_patch` (`gentle.planning_objective.v1`)
  - optional `message`

Planning estimate rule (v1 purchasing simplification):

- Missing required local material adds procurement delay using
  `procurement_business_days_default` (default `10`) unless overridden on the
  item (`inventory.<class>.procurement_business_days`).
- Delay is applied once per missing required material class (deduplicated).
- Business-day model in v1: Monday-Friday only (no holiday calendar).
- Business-day delays are converted to `estimated_time_hours` with a
  deterministic weekend-aware factor (`24h * 7/5` per business day).
- When helper/routine preferences are present in the planning objective,
  `routines list` and `routines compare` emit a synthesized
  `routine_preference_context` and a transparent
  `routine_family_alignment_bonus` inside each planning estimate explanation.
- When `--seq-id` is present and the sequence has construct reasoning with
  variant assay suggestions, that same context also contributes:
  - `construct_reasoning_seq_id`
  - `variant_effect_tags`
  - `variant_suggested_assay_ids`
  - `variant_derived_preferred_routine_families`
  - `routine_family_alignment_sources += variant_derived`

Shipped starter assets:

- Legacy pack:
  - `assets/cloning_patterns.json` (`gentle.cloning_patterns.v1`)
- Hierarchical catalog (one template file per macro):
  - `assets/cloning_patterns_catalog/**/*.json`
  - each file schema: `gentle.cloning_pattern_template.v1`
  - folder hierarchy is used by GUI `Patterns` menu hierarchy
- Typed routine manifest:
  - `assets/cloning_routines.json`
  - schema: `gentle.cloning_routines.v1`
  - includes `gibson.two_fragment_overlap_preview`
    -> `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`
- Import commands:
  - `gentle_cli shell 'macros template-import assets/cloning_patterns.json'`
  - `gentle_cli shell 'macros template-import assets/cloning_patterns_catalog'`

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

Guide-design commands (`gentle_cli guides ...` and `gentle_cli shell 'guides ...'`):

- `guides list`
  - Lists persisted guide sets from guide-design metadata.
- `guides show GUIDE_SET_ID [--limit N] [--offset N]`
  - Pages guide rows from one guide set.
- `guides put GUIDE_SET_ID (--json JSON|@FILE|--file PATH)`
  - Creates/replaces one guide set from JSON array payload.
  - Input payload is `Vec<GuideCandidate>` objects.
- `guides delete GUIDE_SET_ID`
  - Deletes one guide set.
- `guides filter GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]`
  - Applies practical guide constraints and persists a filter report.
  - Optional `--output-set` writes passed guides into a new set.
- `guides filter-show GUIDE_SET_ID`
  - Returns persisted practical-filter report (reasons, warnings, metrics).
- `guides oligos-generate GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]`
  - Generates guide oligos from selected template.
  - Built-ins include `lenti_bsmbi_u6_default` and `plain_forward_reverse`.
- `guides oligos-list [--guide-set GUIDE_SET_ID]`
  - Lists persisted oligo sets (optionally filtered by source guide set).
- `guides oligos-show OLIGO_SET_ID`
  - Shows one oligo set with generated forward/reverse oligos.
- `guides oligos-export GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]`
  - Exports oligo rows for ordering or plate layouts.
- `guides protocol-export GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]`
  - Exports a human-readable wet-lab protocol text for generated oligos.

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

### Example operations

Load a file:

```json
{"LoadFile":{"path":"test_files/pGEX-3X.gb","as_id":"pgex"}}
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

Pairwise sequence alignment (global/local) with structured result payload:

```json
{"AlignSequences":{"query_seq_id":"tp73_refs_mrna_NM_001204186","target_seq_id":"tp73_refs_exon_reference","mode":"local","match_score":2,"mismatch_score":-3,"gap_open":-5,"gap_extend":-1}}
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
  --output-id vkorc1_rs9923231_promoter_alternate
```

These direct shell/CLI routes keep the promoter classification, fragment
suggestion, and allele materialization on the same shared engine contracts used
by operation JSON, workflows, and the updated tutorial path.

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

Notes:

- `PrepareGenome` is intended as a one-time setup step per genome and cache
  location.
- During prepare, GENtle also attempts to create a BLAST nucleotide index
  (`makeblastdb`) under the genome install directory.
- HTTP-based source downloads are resumable (Range requests with retry/backoff),
  and completed installs persist SHA-1 checksums for sequence/annotation files.
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

Save as GenBank:

```json
{"SaveFile":{"seq_id":"pgex_frag_1","path":"frag1.gb","format":"GenBank"}}
```

### Current limitations in the new operation layer

- `import-pool` is currently an adapter-level utility contract (CLI/GUI shared
  shell + JS/Lua wrappers) and not yet an engine operation.

## Error behavior

- REPLs print runtime errors without exiting automatically.
- Invalid file paths or unsupported content produce load/write errors.
- Some functions assume computed annotations are available after `load_dna`.
