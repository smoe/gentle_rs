# gentle-cloning skill scaffold

This folder is a ClawBio/OpenClaw-ready skill scaffold for GENtle.

## Files

- `SKILL.md`: skill metadata + routing/instructions
- `gentle_cloning.py`: wrapper executable
- `catalog_entry.json`: ready-to-paste object for ClawBio `skills/catalog.json`
- `examples/*.json`: request payload examples
- `tests/test_gentle_cloning.py`: minimal wrapper tests

## Local smoke test

```bash
python gentle_cloning.py --demo --output /tmp/gentle_clawbio_demo
```

## Request schema

- `schema`: `gentle.clawbio_skill_request.v1`
- `mode`: `capabilities|state-summary|shell|op|workflow|raw`
- optional: `state_path`, `timeout_secs`

Mode-specific fields:

- `shell`: `shell_line`
- `op`: `operation`
- `workflow`: `workflow` or `workflow_path`
- `raw`: `raw_args[]`
