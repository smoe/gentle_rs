# ClawBio integration scaffold

This directory contains a deployable skill scaffold for exposing GENtle through
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

Deployment into a ClawBio checkout:

Choose one of these layouts:

- Local development: use a symbolic link so ClawBio always sees the live GENtle
  checkout.

  ```bash
  cd /path/to/ClawBio
  mkdir -p skills
  ln -sfn /path/to/gentle_rs/integrations/clawbio/skills/gentle-cloning \
    skills/gentle-cloning
  ```

  This is the preferred developer setup when GENtle and ClawBio are on the same
  machine/filesystem. If ClawBio caches skill metadata or prompt files, restart
  the ClawBio process after changing the GENtle checkout.

- Remote checkout or container staging: use `rsync` so executable bits and
  deletions stay in sync.

  ```bash
  rsync -a --delete \
    /path/to/gentle_rs/integrations/clawbio/skills/gentle-cloning/ \
    /path/to/ClawBio/skills/gentle-cloning/
  ```

- Shell-only transfer between two directories: use a tandem tar stream.

  ```bash
  mkdir -p /path/to/ClawBio/skills
  (cd /path/to/gentle_rs/integrations/clawbio/skills && tar -cf - gentle-cloning) \
    | (cd /path/to/ClawBio/skills && tar -xpf -)
  ```

- Simple one-off install: copy the directory into the ClawBio checkout under
  `skills/gentle-cloning/`.

After placing the skill directory:

1. Set `GENTLE_CLI_CMD` to the recommended Docker-backed GENtle CLI route, for
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
2. Run:
   - `python clawbio.py run gentle-cloning --demo`
   - `python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_skill_info.json --output <dir>`
   - `python clawbio.py run gentle-cloning --input <request.json> --output <dir>`

   The demo now starts with one deterministic graphical GENtle export so the
   first ClawBio reply can show a figure. The broader `capabilities` command is
   presented afterward as the suggested next step.

Catalog registration:

1. Register the object from `skills/gentle-cloning/catalog_entry.json`.
   - symlink setups can read it directly from the linked skill directory
   - `rsync`, tar, and copy setups place it under the ClawBio checkout
2. Add the object under `skills[]` in `skills/catalog.json` (or regenerate
   catalog via ClawBio's `scripts/generate_catalog.py` flow).

Recognition smoke tests:

1. Confirm the skill files are visible from the ClawBio checkout:

   ```bash
   ls -l skills/gentle-cloning/SKILL.md
   ls -l skills/gentle-cloning/examples/request_runtime_version.json
   ```

2. Confirm the ClawBio runtime registry knows the skill:

   ```bash
   python clawbio.py list | grep -F gentle-cloning
   ```

   If this fails, Roboterry cannot route to GENtle yet. The problem is skill
   registration/catalog loading, not GENtle execution.

3. Confirm direct skill execution works without the chat router:

   ```bash
   python clawbio.py run gentle-cloning \
     --input skills/gentle-cloning/examples/request_runtime_version.json \
     --output /tmp/gentle_runtime_version
   cat /tmp/gentle_runtime_version/result.json
   grep -R "GENtle ClawBio skill wrapper invoked" /tmp/gentle_runtime_version
   ```

   Expected result: `status` is `ok`, `request.mode` is `version`, and
   `chat_summary_lines[]` contains the installed GENtle runtime version.
   The `result.json` and `report.md` files also carry the invocation marker
   `GENtle ClawBio skill wrapper invoked`; this marker is safe to grep for in
   ClawBio/Roboterry subprocess logs.

4. Confirm Roboterry/chat routing invokes the skill, not only prose:

   ```text
   Please use the gentle-cloning skill to run the GENtle runtime version check.
   ```

   A good reply should mention that the `gentle-cloning` skill was run and
   should produce or reference a ClawBio output bundle. If direct execution
   works but Roboterry does not invoke the skill, the remaining issue is
   Roboterry routing, prompt policy, or process caching. Restart the Roboterry
   / ClawBio chat process after updating linked or copied skill files.

Roboterry failure interpretation:

- If steps 1-3 pass but Roboterry replies with prose such as
  "I don't have the capability to run GENtle", GENtle and the ClawBio skill are
  working. The chat layer did not dispatch to `clawbio.py run gentle-cloning`.
- In that case, inspect the Roboterry process/configuration rather than the
  GENtle wrapper:
  - confirm Roboterry is running from the same ClawBio checkout where
    `python clawbio.py list | grep -F gentle-cloning` succeeds
  - confirm Roboterry has been restarted after skill/catalog changes
  - confirm Roboterry has tool/skill execution enabled for the current chat
    channel or user
  - inspect Roboterry logs for the raw Telegram text, selected route/skill,
    and whether any `clawbio.py run gentle-cloning ...` subprocess was started
  - if no subprocess was started, fix Roboterry routing/prompt policy; if a
    subprocess was started, inspect that run bundle's `result.json`
  - if the subprocess stdout/stderr is logged, grep for
    `GENtle ClawBio skill wrapper invoked`

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
