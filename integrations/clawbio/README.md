# ClawBio integration scaffold

This directory contains a copy-ready skill scaffold for exposing GENtle through
the ClawBio/OpenClaw skill system.

Current scaffold:

- `skills/gentle-cloning/SKILL.md`
- `skills/gentle-cloning/gentle_cloning.py`
- `skills/gentle-cloning/gentle_apptainer_cli.sh`
- `skills/gentle-cloning/examples/*.json`
- `skills/gentle-cloning/tests/*`
- `skills/gentle-cloning/catalog_entry.json` (ready-to-paste
  `skills/catalog.json` entry)
- `experimental_followup_request_catalog.json` (machine-readable intent catalog
  for ClawBio planners, including optional decision contexts such as upstream
  variant-prioritization assumptions)
- `generate_experimental_followup_catalog_graph.py` (deterministic generator
  for the catalog-derived Mermaid graph)
- `experimental_followup_catalog_graph.mmd` (generated graph of the current
  catalog entries and links)
- `experimental_followup_graph.md` (annotated graph explaining how request
  origin, intent, evidence, artifacts, routine planning, and confirmation gates
  relate)
- `experimental_followup_requests.md` (requests for ClawBio developers who
  want prompts such as "determine the effect of this SNP", "what should we do
  with this differentially expressed gene", or "characterize this splice
  variant" to route through external evidence plus GENtle assay-planning
  artifacts)

Intended use:

1. Copy `integrations/clawbio/skills/gentle-cloning/` into a ClawBio checkout
   under `skills/gentle-cloning/`.
2. Set `GENTLE_CLI_CMD` to the recommended Docker-backed GENtle CLI route, for
   example:

   ```bash
   export GENTLE_CLI_CMD='docker run --rm -i -v "$PWD":/work -w /work ghcr.io/smoe/gentle_rs:cli'
   ```

   Alternatives:
   - use the included Apptainer/Singularity launcher:

     ```bash
     apptainer pull gentle.sif docker://ghcr.io/smoe/gentle_rs:cli
     export GENTLE_CLI_CMD="$PWD/skills/gentle-cloning/gentle_apptainer_cli.sh $PWD/gentle.sif"
     ```

     This launcher intentionally keeps using `apptainer exec ... gentle_cli ...`
     for deterministic wrapper execution, even though `apptainer run
     gentle.sif capabilities` also works on the headless `:cli` image.

   - install `gentle_cli` locally and rely on `PATH`
   - pass `--gentle-cli "<command>"`
   - use repository-local `cargo run --quiet --bin gentle_cli --`
3. Run:
   - `python clawbio.py run gentle-cloning --demo`
   - `python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_skill_info.json --output <dir>`
   - `python clawbio.py run gentle-cloning --input <request.json> --output <dir>`

   The demo now starts with one deterministic graphical GENtle export so the
   first ClawBio reply can show a figure. The broader `capabilities` command is
   presented afterward as the suggested next step.

Catalog registration:

1. Copy `integrations/clawbio/skills/gentle-cloning/catalog_entry.json`
   into your ClawBio checkout.
2. Add the object under `skills[]` in `skills/catalog.json` (or regenerate
   catalog via ClawBio's `scripts/generate_catalog.py` flow).

The wrapper emits the standard ClawBio-style bundle:

- `report.md`
- `result.json`
  - includes raw `stdout`/`stderr`
  - includes `stdout_json` when wrapped GENtle output is valid JSON
  - includes `chat_summary_lines[]` for compact `InspectSequenceContextView`
    replies
- `reproducibility/commands.sh`
- `reproducibility/environment.yml`
- `reproducibility/checksums.sha256`
