import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BRIDGE = REPO_ROOT / "scripts" / "codex-agent-bridge"


def load_bridge_module():
    loader = importlib.machinery.SourceFileLoader("codex_agent_bridge", str(BRIDGE))
    spec = importlib.util.spec_from_loader("codex_agent_bridge", loader)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def request_payload():
    return {
        "schema": "gentle.agent_request.v1",
        "system_id": "codex_local_stdio",
        "prompt": "What can GENtle do?",
        "sent_at_unix_ms": 0,
        "state_summary": {"sequence_count": 0},
    }


def test_extracts_json_from_fenced_codex_output():
    bridge = load_bridge_module()
    parsed = bridge.parse_codex_response(
        textwrap.dedent(
            """\
            Here is the response:

            ```json
            {
              "schema": "gentle.agent_response.v1",
              "assistant_message": "ready",
              "suggested_commands": [
                {"command": "state-summary"}
              ]
            }
            ```

            trailing text
            """
        )
    )

    assert parsed["assistant_message"] == "ready"
    assert parsed["suggested_commands"][0]["execution"] == "ask"


def test_rejects_unknown_agent_response_fields():
    bridge = load_bridge_module()
    try:
        bridge.parse_codex_response(
            json.dumps(
                {
                    "schema": "gentle.agent_response.v1",
                    "assistant_message": "ready",
                    "commands": [],
                }
            )
        )
    except ValueError as exc:
        assert "unsupported fields" in str(exc)
    else:
        raise AssertionError("unknown response field should be rejected")


def test_bridge_keeps_stdout_json_only_when_codex_is_chatty(tmp_path):
    fake_codex = tmp_path / "codex"
    args_path = tmp_path / "args.json"
    prompt_path = tmp_path / "prompt.txt"
    fake_codex.write_text(
        textwrap.dedent(
            f"""\
            #!{sys.executable}
            import json
            import pathlib
            import sys

            args = sys.argv[1:]
            pathlib.Path({str(args_path)!r}).write_text(json.dumps(args), encoding="utf-8")
            pathlib.Path({str(prompt_path)!r}).write_text(sys.stdin.read(), encoding="utf-8")
            output_path = pathlib.Path(args[args.index("--output-last-message") + 1])
            output_path.write_text('''Preamble
            ```json
            {{"schema":"gentle.agent_response.v1","assistant_message":"Codex ready","questions":[],"suggested_commands":[{{"command":"state-summary","execution":"auto"}}]}}
            ```
            trailing''', encoding="utf-8")
            print("codex progress on stdout")
            print("codex progress on stderr", file=sys.stderr)
            """
        ),
        encoding="utf-8",
    )
    fake_codex.chmod(0o755)

    env = os.environ.copy()
    env["CODEX_BIN"] = str(fake_codex)
    env["GENTLE_AGENT_READ_TIMEOUT_SECS"] = "30"
    completed = subprocess.run(
        [sys.executable, str(BRIDGE)],
        input=json.dumps(request_payload()),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stderr == ""
    payload = json.loads(completed.stdout)
    assert payload["assistant_message"] == "Codex ready"
    assert payload["suggested_commands"][0]["execution"] == "auto"

    args = json.loads(args_path.read_text(encoding="utf-8"))
    assert "--sandbox" in args
    assert "read-only" in args
    assert "--skip-git-repo-check" in args
    assert "--ephemeral" in args
    assert "--ignore-user-config" in args
    assert "--output-schema" in args
    assert "GENtle agent request" in prompt_path.read_text(encoding="utf-8")


def test_bridge_reports_login_failure_from_codex_stderr(tmp_path):
    fake_codex = tmp_path / "codex"
    fake_codex.write_text(
        textwrap.dedent(
            f"""\
            #!{sys.executable}
            import sys
            print("please login with your Codex account", file=sys.stderr)
            sys.exit(2)
            """
        ),
        encoding="utf-8",
    )
    fake_codex.chmod(0o755)

    env = os.environ.copy()
    env["CODEX_BIN"] = str(fake_codex)
    completed = subprocess.run(
        [sys.executable, str(BRIDGE)],
        input=json.dumps(request_payload()),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        check=False,
    )

    assert completed.returncode != 0
    assert completed.stdout == ""
    error = json.loads(completed.stderr)
    assert error["schema"] == "gentle.codex_agent_bridge_error.v1"
    assert error["code"] == "codex_not_logged_in"
