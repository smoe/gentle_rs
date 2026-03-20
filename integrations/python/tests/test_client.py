from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
import sys

_THIS_FILE = Path(__file__).resolve()
sys.path.insert(0, str(_THIS_FILE.parents[1]))

from gentle_py import GentleCliError, GentleClient


_FAKE_CLI = """#!/usr/bin/env bash
set -euo pipefail
LOG_PATH="$1"
shift
printf '%s\\n' "$@" >> "$LOG_PATH"

if [[ "${1:-}" == "--state" ]]; then
  shift
  shift
fi

cmd="${1:-}"
case "$cmd" in
  capabilities)
    echo '{"protocol_version":"v1","deterministic_operation_log":true}'
    ;;
  state-summary)
    echo '{"sequence_count":2}'
    ;;
  op)
    echo '{"state_changed":true}'
    ;;
  workflow)
    echo '{"state_changed":true,"results":[]}'
    ;;
  shell)
    shift
    if [[ "${1:-}" == "state-summary" ]]; then
      echo '{"shell":"ok"}'
    else
      echo 'shell text output'
    fi
    ;;
  fail)
    echo 'AGENT_SCHEMA_VALIDATION: synthetic failure' >&2
    exit 7
    ;;
  *)
    echo "UNKNOWN_COMMAND: $cmd" >&2
    exit 9
    ;;
esac
"""


class GentleClientTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        root = Path(self.tempdir.name)
        self.log_path = root / "args.log"
        self.script = root / "fake_gentle_cli.sh"
        self.script.write_text(_FAKE_CLI, encoding="utf-8")
        self.script.chmod(0o755)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _client(self) -> GentleClient:
        return GentleClient(
            cli_cmd=[str(self.script), str(self.log_path)],
            state_path=".gentle_state.json",
            default_timeout_secs=10,
        )

    def test_capabilities_and_state_summary(self) -> None:
        client = self._client()
        caps = client.capabilities()
        self.assertEqual(caps["protocol_version"], "v1")
        summary = client.state_summary()
        self.assertEqual(summary["sequence_count"], 2)

    def test_op_and_workflow(self) -> None:
        client = self._client()
        op_result = client.op({"SetDisplayVisibility": {"target": "Features", "visible": False}})
        self.assertTrue(op_result["state_changed"])
        wf_result = client.workflow({"steps": []})
        self.assertIn("results", wf_result)

    def test_shell_json_and_plain(self) -> None:
        client = self._client()
        json_run = client.shell("state-summary", expect_json=True)
        self.assertEqual(json_run.json["shell"], "ok")
        text_run = client.shell("help")
        self.assertIn("shell text output", text_run.stdout)

    def test_error_mapping_contains_code(self) -> None:
        client = self._client()
        with self.assertRaises(GentleCliError) as cm:
            client.run(["fail"])
        self.assertEqual(cm.exception.code, "AGENT_SCHEMA_VALIDATION")
        self.assertEqual(cm.exception.exit_code, 7)

    def test_state_flag_is_forwarded(self) -> None:
        client = self._client()
        client.capabilities()
        log = self.log_path.read_text(encoding="utf-8")
        self.assertIn("--state", log)
        self.assertIn(".gentle_state.json", log)

    def test_render_dotplot_svg_wrapper_uses_op_payload(self) -> None:
        client = self._client()
        out = client.render_dotplot_svg(
            "seq_a",
            "dotplot_primary",
            "/tmp/dotplot.svg",
            flex_track_id="flex_primary",
            display_density_threshold=0.2,
            display_intensity_gain=1.8,
        )
        self.assertTrue(out["state_changed"])
        log = self.log_path.read_text(encoding="utf-8")
        self.assertIn("\"RenderDotplotSvg\"", log)
        self.assertIn("\"seq_id\":\"seq_a\"", log)
        self.assertIn("\"dotplot_id\":\"dotplot_primary\"", log)
        self.assertIn("\"path\":\"/tmp/dotplot.svg\"", log)
        self.assertIn("\"flex_track_id\":\"flex_primary\"", log)
        self.assertIn("\"display_density_threshold\":0.2", log)
        self.assertIn("\"display_intensity_gain\":1.8", log)

    def test_render_dotplot_svg_wrapper_validates_inputs(self) -> None:
        client = self._client()
        with self.assertRaises(GentleCliError):
            client.render_dotplot_svg("", "dp", "/tmp/out.svg")
        with self.assertRaises(GentleCliError):
            client.render_dotplot_svg("seq", "dp", "/tmp/out.svg", display_density_threshold=float("nan"))


if __name__ == "__main__":
    unittest.main()
