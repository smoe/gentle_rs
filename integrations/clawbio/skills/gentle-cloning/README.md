# gentle-cloning skill scaffold

This folder is a ClawBio/OpenClaw-ready skill scaffold for GENtle.

## Files

- `SKILL.md`: skill metadata + routing/instructions
- `gentle_cloning.py`: wrapper executable
- `gentle_apptainer_cli.sh`: Apptainer/Singularity launcher for `gentle_cli`
- `catalog_entry.json`: ready-to-paste object for ClawBio `skills/catalog.json`
- `examples/*.json`: request payload examples
- `tests/test_gentle_cloning.py`: minimal wrapper tests

## Recommended smoke test

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
