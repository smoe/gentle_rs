import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
BRIDGE = REPO_ROOT / "scripts" / "pi-agent-bridge"


def load_bridge_module():
    loader = importlib.machinery.SourceFileLoader("pi_agent_bridge", str(BRIDGE))
    spec = importlib.util.spec_from_loader("pi_agent_bridge", loader)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def request_payload():
    return {
        "schema": "gentle.agent_request.v1",
        "system_id": "pi_local_stdio",
        "prompt": "What can GENtle do?",
        "sent_at_unix_ms": 0,
        "state_summary": {"sequence_count": 0},
    }


class PiAgentBridgeTests(unittest.TestCase):
    def test_pi_command_is_ephemeral_and_tool_free(self):
        bridge = load_bridge_module()
        with mock.patch.dict(
            os.environ, {"GENTLE_AGENT_MODEL": "openai-codex/gpt-5.4"}
        ):
            command = bridge.pi_command("pi")

        for flag in (
            "--print",
            "--no-session",
            "--no-tools",
            "--no-extensions",
            "--no-skills",
            "--no-prompt-templates",
            "--no-context-files",
            "--no-approve",
        ):
            self.assertIn(flag, command)
        model_flag = command.index("--model")
        self.assertEqual(command[model_flag + 1], "openai-codex/gpt-5.4")

    def test_parse_pi_response_accepts_one_json_fence(self):
        bridge = load_bridge_module()
        parsed = bridge.parse_pi_response(
            textwrap.dedent(
                """\
                ```json
                {
                  "schema": "gentle.agent_response.v1",
                  "assistant_message": "ready",
                  "questions": [],
                  "suggested_commands": []
                }
                ```
                """
            )
        )

        self.assertEqual(parsed["assistant_message"], "ready")

    def test_parse_pi_response_rejects_prose_around_fence(self):
        bridge = load_bridge_module()
        with self.assertRaises((json.JSONDecodeError, ValueError)):
            bridge.parse_pi_response(
                """Here is the answer:
```json
{"schema":"gentle.agent_response.v1","assistant_message":"ready"}
```"""
            )

    def test_bridge_invokes_fake_pi_without_tools_and_preserves_context(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_pi = temp_path / "pi"
            args_path = temp_path / "args.json"
            prompt_path = temp_path / "prompt.txt"
            cwd_path = temp_path / "cwd.txt"
            fake_pi.write_text(
                textwrap.dedent(
                    f"""\
                    #!{sys.executable}
                    import json
                    import os
                    import pathlib
                    import sys

                    pathlib.Path({str(args_path)!r}).write_text(json.dumps(sys.argv[1:]), encoding="utf-8")
                    pathlib.Path({str(prompt_path)!r}).write_text(sys.stdin.read(), encoding="utf-8")
                    pathlib.Path({str(cwd_path)!r}).write_text(os.getcwd(), encoding="utf-8")
                    print('''```json
                    {{"schema":"gentle.agent_response.v1","assistant_message":"Pi ready","questions":[],"suggested_commands":[{{"command":"state-summary","execution":"ask"}}]}}
                    ```''')
                    """
                ),
                encoding="utf-8",
            )
            fake_pi.chmod(0o755)
            request = request_payload()
            request["x_conversation"] = {
                "schema": "gentle.agent_conversation.v1",
                "turns": [{"user_message": "Use homo_sapiens."}],
            }
            request["x_local_documents"] = {
                "schema": "gentle.agent_local_documents.v1",
                "documents": [
                    {
                        "requested_path": "/docs/roadmap.md",
                        "content": "Run the GUI smoke checklist.",
                    }
                ],
            }

            env = os.environ.copy()
            env["PI_BIN"] = str(fake_pi)
            completed = subprocess.run(
                [sys.executable, str(BRIDGE)],
                input=json.dumps(request),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(completed.stderr, "")
            self.assertEqual(json.loads(completed.stdout)["assistant_message"], "Pi ready")
            args = json.loads(args_path.read_text(encoding="utf-8"))
            self.assertIn("--no-tools", args)
            self.assertIn("--no-session", args)
            prompt = prompt_path.read_text(encoding="utf-8")
            self.assertIn("x_conversation", prompt)
            self.assertIn("Use homo_sapiens.", prompt)
            self.assertIn("x_local_documents", prompt)
            self.assertIn("Run the GUI smoke checklist.", prompt)
            self.assertIn("do not ask the user to paste", prompt)
            self.assertTrue(
                Path(cwd_path.read_text(encoding="utf-8")).name.startswith(
                    "gentle-pi-agent-"
                )
            )

    def test_bridge_reports_missing_pi_login(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fake_pi = Path(temp_dir) / "pi"
            fake_pi.write_text(
                textwrap.dedent(
                    f"""\
                    #!{sys.executable}
                    import sys
                    print("No models available. Use /login to log into a provider.", file=sys.stderr)
                    sys.exit(1)
                    """
                ),
                encoding="utf-8",
            )
            fake_pi.chmod(0o755)
            env = os.environ.copy()
            env["PI_BIN"] = str(fake_pi)

            completed = subprocess.run(
                [sys.executable, str(BRIDGE)],
                input=json.dumps(request_payload()),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
                check=False,
            )

            self.assertNotEqual(completed.returncode, 0)
            self.assertEqual(completed.stdout, "")
            error = json.loads(completed.stderr)
            self.assertEqual(error["schema"], "gentle.pi_agent_bridge_error.v1")
            self.assertEqual(error["code"], "pi_not_logged_in")


if __name__ == "__main__":
    unittest.main()
