# Python CLI adapter

This folder contains `gentle_py`, a thin Python adapter over `gentle_cli`.
It is a process-level scripting/notebook convenience layer, not a native Python
biology engine and not a natural-language or ClawBio/Telegram parser.

All biology logic remains in GENtle's Rust engine and shared CLI routes. The
Python package builds deterministic `gentle_cli` invocations, parses JSON when
the selected route is expected to return JSON, and raises structured Python
errors when the subprocess fails.

Path:

- `integrations/python/gentle_py/`

Requirements:

- Python >= 3.10
- a resolvable `gentle_cli`

CLI resolution order:

1. `GentleClient(cli_cmd=...)`
2. `GENTLE_CLI_CMD`
3. `gentle_cli` on `PATH`
4. repository fallback: `cargo run --quiet --bin gentle_cli --`

Core API surface:

- `GentleClient.capabilities()`
- `GentleClient.state_summary()`
- `GentleClient.op(operation)`
- `GentleClient.workflow(workflow|workflow_path=...)`
- `GentleClient.shell(line, expect_json=False)`
- `GentleClient.run(args, expect_json=False)`

Error handling:

- failing subprocesses raise `GentleCliError`
- the exception carries `code`, `command`, `exit_code`, `stdout`, and `stderr`

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

Install as an editable local package:

```bash
python3 -m pip install -e integrations/python
```

Test:

```bash
python3 -m unittest -q integrations/python/tests/test_client.py
```
