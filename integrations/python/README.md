# Python interface scaffold

This folder contains a thin Python interface to GENtle that wraps `gentle_cli`
without duplicating biology logic.

Path:

- `integrations/python/gentle_py/`

Design intent:

- deterministic command execution through `gentle_cli`
- structured JSON handling for protocol-first routes (`capabilities`,
  `state-summary`, `op`, `workflow`)
- stable error object for scripting automation

Convenience wrappers:

- `render_dotplot_svg(seq_id, dotplot_id, output_svg, ...)`
  - submits engine operation `RenderDotplotSvg` through deterministic `op(...)`
  - optional arguments:
    - `flex_track_id`
    - `display_density_threshold`
    - `display_intensity_gain`

Quick start:

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

Install as local package:

```bash
python3 -m pip install -e integrations/python
```

Test:

```bash
python3 -m unittest -q integrations/python/tests/test_client.py
```
