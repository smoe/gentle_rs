import importlib.machinery
import importlib.util
import json
import os
import shutil
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
    def test_render_prompt_includes_host_command_contract_and_local_first_rule(self):
        bridge = load_bridge_module()
        request = request_payload()
        request["x_local_references"] = {
            "schema": "gentle.agent_local_reference_context.v1",
            "references": [
                {
                    "genome_id": "Human GRCh38 Ensembl 116",
                    "gene_extraction_ready": True,
                }
            ],
        }
        request["x_helper_catalog"] = {
            "schema": "gentle.agent_helper_catalog_context.v1",
            "cards": [
                {
                    "helper_id": "Reporter Promega Luciferase AY738222 (online)",
                    "catalog_number": "E6651",
                }
            ],
        }
        request["x_gui_context"] = {
            "schema": "gentle.agent_gui_context.v1",
            "host_available": True,
            "recent_projects": [
                {
                    "item_id": "recent-deadbeef",
                    "open_command": "ui open recent-project recent-deadbeef",
                }
            ],
            "tutorial_projects": [],
            "configuration_sections": [],
        }
        host_contract = (
            "Use /fetch ensembl TP73 --species homo_sapiens --assembly GRCh38 "
            "--flank-bp 10000 --id tp73_grch38."
        )
        with mock.patch.dict(
            os.environ,
            {bridge.HOST_SYSTEM_PROMPT_ENV: host_contract},
        ):
            prompt = bridge.render_pi_prompt(request)

        self.assertIn(host_contract, prompt)
        self.assertIn("x_local_references", prompt)
        self.assertIn("gene_extraction_ready local reference", prompt)
        self.assertIn("local reference ids", prompt)
        self.assertIn("x_helper_catalog", prompt)
        self.assertIn("E6651", prompt)
        self.assertIn("before guessing a reporter/vector identity", prompt)
        self.assertIn("x_gui_context", prompt)
        self.assertIn("authoritative GUI-host catalogs", prompt)
        self.assertIn("ui open recent-project recent-deadbeef", prompt)

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

    def test_pi_command_enables_only_co_shipped_public_web_tools(self):
        bridge = load_bridge_module()
        command = bridge.pi_command("pi", allow_web_research=True)

        self.assertIn("--no-builtin-tools", command)
        self.assertNotIn("--no-tools", command)
        self.assertIn("--no-extensions", command)
        extension_flag = command.index("--extension")
        self.assertEqual(
            Path(command[extension_flag + 1]),
            bridge.WEB_EXTENSION,
        )
        tools_flag = command.index("--tools")
        self.assertEqual(command[tools_flag + 1], bridge.WEB_TOOL_NAMES)
        for forbidden in ("read", "bash", "edit", "write", "grep", "find", "ls"):
            self.assertNotIn(forbidden, command[tools_flag + 1].split(","))

    def test_image_attachment_is_forwarded_as_pi_file_argument(self):
        bridge = load_bridge_module()
        with tempfile.TemporaryDirectory() as temp_dir:
            image_path = Path(temp_dir) / "problem.png"
            image_path.write_bytes(b"synthetic png fixture")
            request = request_payload()
            request["x_attachments"] = [
                {"kind": "image", "path": str(image_path)}
            ]

            paths = bridge.attachment_paths(request)
            command = bridge.pi_command("pi", paths)

            self.assertIn(f"@{image_path}", command)
            prompt = bridge.render_pi_prompt(request)
            self.assertNotIn(str(image_path), prompt)
            self.assertIn("<local attachment>", prompt)

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

    def test_parse_pi_response_accepts_bounded_screenshot_request(self):
        bridge = load_bridge_module()
        parsed = bridge.parse_pi_response(
            json.dumps(
                {
                    "schema": "gentle.agent_response.v1",
                    "assistant_message": "",
                    "questions": [],
                    "suggested_commands": [],
                    "screenshot_request": {
                        "id": "inspect-tp73-map-1",
                        "reason": "Inspect the visible TP73 feature lanes.",
                    },
                }
            )
        )

        self.assertEqual(parsed["screenshot_request"]["id"], "inspect-tp73-map-1")
        self.assertIn("genuinely needed", bridge.render_pi_prompt(request_payload()))

    def test_parse_pi_response_rejects_screenshot_target_fields(self):
        bridge = load_bridge_module()
        response = {
            "schema": "gentle.agent_response.v1",
            "assistant_message": "",
            "questions": [],
            "suggested_commands": [],
            "screenshot_request": {
                "id": "inspect-map",
                "reason": "Inspect the map.",
                "path": "/tmp/problem.png",
            },
        }

        with self.assertRaisesRegex(ValueError, "unsupported fields"):
            bridge.validate_agent_response(response)

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

    def test_bridge_web_opt_in_passes_restricted_tools_and_returns_audit(self):
        bridge = load_bridge_module()
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_pi = temp_path / "pi"
            args_path = temp_path / "args.json"
            fake_pi.write_text(
                textwrap.dedent(
                    f"""\
                    #!{sys.executable}
                    import json
                    import os
                    import pathlib
                    import sys

                    pathlib.Path({str(args_path)!r}).write_text(json.dumps(sys.argv[1:]), encoding="utf-8")
                    log_path = pathlib.Path(os.environ["GENTLE_AGENT_WEB_LOG_PATH"])
                    log_path.write_text(json.dumps({{
                        "kind": "search",
                        "query": "Promega pGL4.10 E6651",
                        "retrieved_at_unix_ms": 123,
                        "results": [{{
                            "title": "pGL4.10 Vector Protocol",
                            "url": "https://www.promega.com/pgl410"
                        }}]
                    }}) + "\\n" + json.dumps({{
                        "kind": "page",
                        "requested_url": "https://www.promega.com/pgl410",
                        "final_url": "https://www.promega.com/pgl410",
                        "title": "pGL4.10 Vector Protocol",
                        "retrieved_at_unix_ms": 124,
                        "content_sha256": "a" * 64,
                        "included_char_count": 321,
                        "truncated": False
                    }}) + "\\n" + json.dumps({{
                        "kind": "warning",
                        "message": "One optional source could not be read"
                    }}) + "\\n", encoding="utf-8")
                    print(json.dumps({{
                        "schema": "gentle.agent_response.v1",
                        "assistant_message": "Promega identifies pGL4.10 as E6651.",
                        "questions": [],
                        "suggested_commands": []
                    }}))
                    """
                ),
                encoding="utf-8",
            )
            fake_pi.chmod(0o755)
            request = request_payload()
            request["x_web_access"] = {
                "schema": "gentle.agent_web_access.v1",
                "enabled": True,
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
            response = json.loads(completed.stdout)
            research = response["x_web_research"]
            self.assertEqual(research["schema"], "gentle.agent_web_research.v1")
            self.assertEqual(research["searches"][0]["query"], "Promega pGL4.10 E6651")
            self.assertEqual(research["pages"][0]["content_sha256"], "a" * 64)
            self.assertEqual(research["warnings"], ["One optional source could not be read"])
            args = json.loads(args_path.read_text(encoding="utf-8"))
            self.assertIn("--no-builtin-tools", args)
            self.assertNotIn("--no-tools", args)
            self.assertEqual(args[args.index("--tools") + 1], bridge.WEB_TOOL_NAMES)

    def test_disabled_web_request_discards_model_claimed_web_audit(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fake_pi = Path(temp_dir) / "pi"
            fake_pi.write_text(
                textwrap.dedent(
                    f"""\
                    #!{sys.executable}
                    import json
                    print(json.dumps({{
                        "schema": "gentle.agent_response.v1",
                        "assistant_message": "No live research performed.",
                        "questions": [],
                        "suggested_commands": [],
                        "x_web_research": {{"schema": "invented"}}
                    }}))
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

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertNotIn("x_web_research", json.loads(completed.stdout))

    @unittest.skipUnless(shutil.which("node"), "Node.js is required for Pi web core tests")
    def test_web_core_blocks_private_targets_and_parses_public_search_results(self):
        core = REPO_ROOT / "scripts" / "pi-extensions" / "gentle-web-research-core.mjs"
        script = textwrap.dedent(
            f"""\
            import {{ isBlockedAddress, normalizePublicUrl, parseBingRssResults, searchResultsMatchQuery }} from {core.as_uri()!r};
            if (!isBlockedAddress('127.0.0.1') || !isBlockedAddress('10.1.2.3') || !isBlockedAddress('192.0.2.1') || !isBlockedAddress('198.51.100.1') || !isBlockedAddress('203.0.113.1') || !isBlockedAddress('::ffff:7f00:1') || !isBlockedAddress('64:ff9b::7f00:1') || isBlockedAddress('93.184.216.34') || isBlockedAddress('2606:2800:220:1:248:1893:25c8:1946')) process.exit(2);
            try {{ normalizePublicUrl('http://localhost:11434/api'); process.exit(3); }} catch {{}}
            const xml = `<rss><channel><item><title><![CDATA[Primary paper]]></title><link>https://example.org/paper?a=1&amp;b=2</link><description>Evidence &amp; methods</description></item></channel></rss>`;
            const results = parseBingRssResults(xml, 3);
            if (results.length !== 1 || results[0].title !== 'Primary paper' || results[0].url !== 'https://example.org/paper?a=1&b=2' || results[0].snippet !== 'Evidence & methods') process.exit(4);
            if (!searchResultsMatchQuery('primary evidence', results)) process.exit(5);
            if (searchResultsMatchQuery('Promega pGL4.10 E6651', results)) process.exit(6);
            """
        )
        completed = subprocess.run(
            [shutil.which("node"), "--input-type=module", "--eval", script],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

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
