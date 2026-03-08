---
name: gentle-cloning
description: >-
  Deterministic cloning/genome-analysis execution wrapper for GENtle CLI with
  reproducibility bundle output.
version: 0.1.0
author: GENtle project
license: MIT
tags: [cloning, genome-analysis, reproducibility, automation]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env:
        - GENTLE_CLI_CMD
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ulfschubert/gentle_rs
    os: [macos, linux]
    install: []
    trigger_keywords:
      - gentle
      - cloning workflow
      - prepare genome
      - primer design
      - blast sequence
      - genome anchor
---

# 🧬 GENtle Cloning

You are **GENtle Cloning**, a specialised ClawBio skill that routes structured
requests into deterministic `gentle_cli` commands and writes reproducible
artifacts.

## Why This Exists

Direct LLM prompting is useful for planning, but weak for reproducibility.

- **Without it**: command arguments/state paths are often implicit and hard to
  replay.
- **With it**: each run is executed through fixed `gentle_cli` commands and
  emits `report.md`, `result.json`, and a reproducibility bundle.
- **Why ClawBio/OpenClaw**: the skill can be orchestrated with other skills
  while preserving strict command-level auditability.

## Core Capabilities

1. Run deterministic GENtle commands from structured request JSON.
2. Support command modes for `capabilities`, `state-summary`, `shell`, `op`,
   `workflow`, and raw direct argument lists.
3. Emit reproducibility artifacts (commands + environment + checksums).

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Skill request JSON (`gentle.clawbio_skill_request.v1`) | `.json` | `schema`, `mode`; plus mode-specific fields | `examples/request_capabilities.json` |

## Workflow

When the user asks for a GENtle operation:

1. **Validate**: parse request JSON and required mode-specific fields.
2. **Resolve executable**: use `--gentle-cli`, then `GENTLE_CLI_CMD`, then
   `gentle_cli` on `PATH`, then repository-local `cargo run` fallback.
3. **Execute**: run one deterministic `gentle_cli` command with timeout.
4. **Generate**: write `report.md`, `result.json`, and reproducibility bundle.

## CLI Reference

```bash
# Standard usage
python skills/gentle-cloning/gentle_cloning.py \
  --input skills/gentle-cloning/examples/request_capabilities.json \
  --output /tmp/gentle_clawbio_run

# Demo mode
python skills/gentle-cloning/gentle_cloning.py \
  --demo --output /tmp/gentle_clawbio_demo

# Optional explicit command resolver
python skills/gentle-cloning/gentle_cloning.py \
  --input request.json --output /tmp/out \
  --gentle-cli "gentle_cli"
```

## Demo

```bash
python skills/gentle-cloning/gentle_cloning.py --demo --output /tmp/gentle_clawbio_demo
```

Expected output: one command execution report for `gentle_cli capabilities` (or
an explicit degraded-demo message if no runnable GENtle CLI resolver exists on
that machine).

## Algorithm / Methodology

1. Parse request JSON into one execution mode.
2. Canonicalize request mode into one CLI argument vector.
3. Execute with fixed timeout and no hidden retries.
4. Persist structured run metadata plus a human-readable report.

Key parameters:

- `timeout_secs`: execution timeout for subprocess call (default `180`).
- `state_path`: optional explicit GENtle state file path.

## Example Queries

- "Run GENtle capabilities and give me the machine-readable output."
- "Execute this workflow JSON through GENtle and capture reproducibility files."
- "Apply this GENtle operation JSON to a state file and summarize the result."

## Output Structure

```
output_directory/
├── report.md
├── result.json
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

Required:

- Python 3.10+
- A resolvable GENtle CLI command:
  - direct binary on PATH (`gentle_cli`), or
  - `GENTLE_CLI_CMD`, or
  - repository-local `cargo run --bin gentle_cli --` fallback.

Optional:

- Rust toolchain (`cargo`) when using repository fallback mode.

## Safety

- Local-first execution; no mandatory data upload.
- Wrapper performs no hidden mutating steps beyond the explicit command.
- Every run records exact invocation in reproducibility artifacts.
- This skill is for in-silico cloning/genome-analysis tooling and is not a
  clinical diagnostic workflow.

## Integration with Bio Orchestrator

Trigger conditions:

- Requests containing "GENtle", "cloning workflow", "prepare genome",
  "primer design", "BLAST in GENtle", or "genome anchor".

Chaining partners:

- `bio-orchestrator`: route high-level user asks into this deterministic tool.
- `profile-report`/other skills: consume `result.json` output where needed.

## Citations

- [GENtle repository](https://github.com/ulfschubert/gentle_rs)
- [ClawBio](https://github.com/ClawBio/ClawBio)
