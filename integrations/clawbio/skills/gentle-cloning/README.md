# gentle-cloning skill scaffold

This folder is a ClawBio/OpenClaw-ready skill scaffold for GENtle.

## Files

- `SKILL.md`: skill metadata + routing/instructions
- `gentle_cloning.py`: wrapper executable
- `gentle_local_checkout_cli.sh`: local-checkout launcher for `cargo run --bin gentle_cli --`
- `gentle_apptainer_cli.sh`: Apptainer/Singularity launcher for `gentle_cli`
- `catalog_entry.json`: ready-to-paste object for ClawBio `skills/catalog.json`
- `examples/*.json`: request payload examples, including first-run reference/helper preparation
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
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_status_grch38.json --output /tmp/gentle_status_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_genomes_prepare_grch38.json --output /tmp/gentle_prepare_grch38
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_helpers_prepare_puc19.json --output /tmp/gentle_prepare_puc19
```

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

Included first-run bootstrap requests:

- `examples/request_genomes_list_human.json`
- `examples/request_genomes_status_grch38.json`
- `examples/request_genomes_prepare_grch38.json`
- `examples/request_helpers_status_puc19.json`
- `examples/request_helpers_prepare_puc19.json`
