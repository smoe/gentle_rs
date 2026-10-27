# Python Wrapper

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

## Python module wrapper (`gentle_py`)

GENtle also ships a thin Python wrapper over `gentle_cli`:

- module path: `integrations/python/gentle_py/`
- quick docs: `integrations/python/README.md`

Quick usage:

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

Deterministic API surface:

- `capabilities()`
- `state_summary()`
- `op(operation)`
- `workflow(workflow|workflow_path=...)`
- `shell(line, expect_json=False)`
- `render_dotplot_svg(seq_id, dotplot_id, output_svg, ...)`

CLI resolution order:

1. constructor `cli_cmd`
2. `GENTLE_CLI_CMD`
3. `gentle_cli` on `PATH`
4. repository fallback `cargo run --quiet --bin gentle_cli --`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Lua Shell](lua.md) | [Next: MCP Server](mcp.md)
<!-- manual-page-nav:end -->
