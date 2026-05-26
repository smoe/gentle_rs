# GENtle Local Agent Handoff

Audience: Codex, Claude, OpenClaw, and other coding/bench agents that need to
use an already-working local GENtle-through-ClawBio installation.

This is a small routing memory, not a new GENtle interface. Use it when an
agent needs to run GENtle locally and should not rediscover or reimplement the
ClawBio skill setup.

## Boundary

- GENtle owns deterministic cloning, sequence, resource, and artifact
  execution.
- ClawBio owns the skill runner, request bundle, output bundle, and
  chat-facing confirmation workflow.
- OpenClaw or an external coding agent may own the broader gateway/orchestration
  layer.
- A native OpenClaw/Codex/Claude wrapper should therefore delegate to the
  existing ClawBio `gentle-cloning` runner unless the user explicitly asks for
  low-level GENtle debugging.

Do not create a second GENtle command vocabulary in the wrapper. Do not
invent filesystem-search, shell, or OpenClaw-style gateway commands. Route
through the existing ClawBio request examples or typed request JSON.

## When to Use GENtle

Use the `gentle-cloning` skill for:

- cloning design and method comparison,
- plasmid/vector and helper-construct work,
- restriction, Gibson, Golden Gate, Gateway, PCR, qPCR/TaqMan, and assembly
  planning workflows,
- sequence annotation and stateless pasted-DNA inspection,
- gene/locus/genome-context retrieval when supported resources are available
  or can be prepared,
- TFBS/JASPAR, REBASE, ATtRACT-style resource-backed analyses,
- isoform/protein gel and 2D-gel visualizations,
- lab-assistant cloning instruction exports,
- visual review of generated PNG/SVG/report artifacts.

For broad planning questions, prefer a read-only GENtle status/guide/consult
route first. For long-running prepare/sync/download actions, present the
ClawBio suggested action and require user confirmation.

## Known Local Deployment Shape

The concrete deployment reported by the ClawBio/OpenClaw host used:

```text
GENtle checkout:
  /home/clawbio/GENtle

ClawBio checkout:
  /home/clawbio/ClawBio

ClawBio skill symlink:
  /home/clawbio/ClawBio/skills/gentle-cloning
    -> /home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/

GENtle launcher:
  /home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh
```

Treat those paths as a known host profile, not a portable default. On another
machine, resolve the ClawBio checkout, GENtle checkout, and symlink first.

## Preferred Invocation

Run GENtle through ClawBio's runner:

```bash
cd /home/clawbio/ClawBio
GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh \
  python clawbio.py run gentle-cloning \
  --input skills/gentle-cloning/examples/request_runtime_version.json \
  --output /tmp/gentle_runtime_version
```

For other tasks, keep the same `GENTLE_CLI_CMD` and replace `--input` with a
specific `skills/gentle-cloning/examples/*.json` request or an explicitly
constructed `gentle.clawbio_skill_request.v1` file.

Prefer this route over calling `gentle_cli` directly because the ClawBio runner
normalizes:

- `report.md`,
- `result.json`,
- preferred PNG/SVG artifacts,
- `commands.sh`,
- `environment.yml`,
- `checksums.sha256`,
- `suggested_actions[]` follow-ups.

## Result Inspection Order

After a run:

1. Read `result.json` first for `schema`, `status`, `summary_lines`,
   `expected_artifacts`, `preferred_artifacts`, `suggested_actions[]`, and
   structured errors.
2. Read `report.md` for the human-readable explanation.
3. Prefer PNG artifacts for chat display when present; keep SVG artifacts for
   auditability and re-rendering.
4. Inspect `commands.sh`, `environment.yml`, and `checksums.sha256` when
   reproducing or debugging a run.
5. Preserve `suggested_actions[]` as the executable follow-up surface. Do not
   turn prose menu text into commands when structured actions are present.

## Readiness Discipline

Do not hard-code readiness as permanent truth. If a user asks what is available,
run a status route through the skill.

Useful request examples include:

- `request_runtime_version.json`
- `request_skill_info.json`
- `request_services_status.json`
- `request_resources_status.json`
- `request_genomes_status_grch38.json`
- `request_helpers_status_puc19.json`
- `request_services_telegram_guide.json`

The known host smoke-test snapshot reported:

- Available: core GENtle CLI, JASPAR builtin snapshot, REBASE builtin snapshot,
  ATtRACT runtime snapshot.
- Not yet prepared: `Human GRCh38 Ensembl 116`, `Plasmid pUC19 (online)`.
- Missing optional tools: RNAfold/ViennaRNA and rnapkin.

Treat that list as a snapshot. If the user asks now, run `services status`,
`resources status`, `genomes status`, or `helpers status` and quote the fresh
payload.

For missing references or helper tools, probe before proposing installation or
preparation. The agent should first ask GENtle/ClawBio what is already present
and then, only if the user wants setup help, report:

- the exact prepare/install command that would be run,
- the target directory or cache location,
- expected download/build size when known,
- whether the action needs network, elevated permissions, or a package manager,
- platform-specific alternatives for the detected host,
- whether the action is optional for the user's immediate request.

Do not trigger GRCh38, pUC19, RNAfold/ViennaRNA, rnapkin, or other large
downloads/builds from this handoff note. Long-running prepare/sync/download
actions must remain explicit suggested actions with confirmation.

Useful first probes include:

```bash
cd /home/clawbio/ClawBio
GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh \
  python clawbio.py run gentle-cloning \
  --input skills/gentle-cloning/examples/request_services_status.json \
  --output /tmp/gentle_services_status

GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh \
  python clawbio.py run gentle-cloning \
  --input skills/gentle-cloning/examples/request_resources_status.json \
  --output /tmp/gentle_resources_status

GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh \
  python clawbio.py run gentle-cloning \
  --input skills/gentle-cloning/examples/request_genomes_status_grch38.json \
  --output /tmp/gentle_grch38_status

GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh \
  python clawbio.py run gentle-cloning \
  --input skills/gentle-cloning/examples/request_helpers_status_puc19.json \
  --output /tmp/gentle_puc19_status
```

## Minimal Smoke Test

A wrapper installation is healthy when this direct run succeeds:

```bash
cd /home/clawbio/ClawBio
GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh \
  python clawbio.py run gentle-cloning \
  --input skills/gentle-cloning/examples/request_runtime_version.json \
  --output /tmp/gentle_runtime_version
```

Expected:

- exit code `0`,
- `/tmp/gentle_runtime_version/result.json`,
- `/tmp/gentle_runtime_version/report.md`,
- `result.json.status` is `ok`,
- the report mentions the GENtle ClawBio skill wrapper invocation.

For a graphics smoke test, use a demo request that produces PNG/SVG artifacts,
then inspect both `preferred_artifacts` and the actual output directory.

## Copyable Agent Prompt

Use this prompt when registering a thin native wrapper for Codex, Claude, or
OpenClaw:

```text
You have access to GENtle through an already-working ClawBio skill named
gentle-cloning. Do not invent a separate GENtle interface. Use ClawBio's runner
and the existing GENtle skill request examples.

Known local host profile:
- GENtle checkout: /home/clawbio/GENtle
- ClawBio checkout: /home/clawbio/ClawBio
- ClawBio skill symlink:
  /home/clawbio/ClawBio/skills/gentle-cloning ->
  /home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/
- Required launcher env:
  GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh

When a user asks for GENtle, cloning design, plasmid/vector work,
restriction/Gibson/assembly planning, sequence annotation, primer/qPCR design,
isoform/protein-gel figures, TFBS/resource-backed analyses, or lab-assistant
handoffs, route through:

  cd /home/clawbio/ClawBio
  GENTLE_CLI_CMD=/home/clawbio/GENtle/integrations/clawbio/skills/gentle-cloning/gentle_local_checkout_cli.sh \
    python clawbio.py run gentle-cloning --input <request.json> --output <output-dir>

Choose request JSONs from:
  /home/clawbio/ClawBio/skills/gentle-cloning/examples/

After execution, inspect result.json first, then report.md, then preferred
PNG/SVG artifacts and reproducibility files. Preserve structured
suggested_actions[] as the follow-up surface and require confirmation for
prepare/sync/download actions.

Readiness is dynamic. For availability questions, run services/resources/
genomes/helpers status requests instead of answering from memory. Probe first:
do not prepare GRCh38, pUC19, RNAfold/ViennaRNA, rnapkin, or other large
resources/tools until the user has seen the expected command, target location,
download/build size when known, platform-specific install route, and permission
requirements.
```
