# ClawBio integration scaffold

This directory contains a copy-ready skill scaffold for exposing GENtle through
the ClawBio/OpenClaw skill system.

Current scaffold:

- `skills/gentle-cloning/SKILL.md`
- `skills/gentle-cloning/gentle_cloning.py`
- `skills/gentle-cloning/examples/*.json`
- `skills/gentle-cloning/tests/*`
- `skills/gentle-cloning/catalog_entry.json` (ready-to-paste
  `skills/catalog.json` entry)

Intended use:

1. Copy `integrations/clawbio/skills/gentle-cloning/` into a ClawBio checkout
   under `skills/gentle-cloning/`.
2. Install or otherwise resolve `gentle_cli`.
3. Run:
   - `python clawbio.py run gentle-cloning --demo`
   - `python clawbio.py run gentle-cloning --input <request.json> --output <dir>`

Catalog registration:

1. Copy `integrations/clawbio/skills/gentle-cloning/catalog_entry.json`
   into your ClawBio checkout.
2. Add the object under `skills[]` in `skills/catalog.json` (or regenerate
   catalog via ClawBio's `scripts/generate_catalog.py` flow).

The wrapper emits the standard ClawBio-style bundle:

- `report.md`
- `result.json`
- `reproducibility/commands.sh`
- `reproducibility/environment.yml`
- `reproducibility/checksums.sha256`
