#!/usr/bin/env python3
"""Offline headless-container regressions with synthetic executable/receipt data.

Temporary shell stubs exercise the real entrypoint without Docker, GUI, scientific
inputs or network access. Run: python3 -m unittest scripts.test_container -v
Actual Linux linking, dependencies and image execution remain container CI gates.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
import shlex
import subprocess
from tempfile import TemporaryDirectory
import textwrap
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
BINARIES = ["gentle_cli", "gentle_mcp", "gentle_examples_docs"]


class HeadlessEntrypointTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bin_dir = self.root / "bin"
        self.bin_dir.mkdir()
        for binary in BINARIES:
            path = self.bin_dir / binary
            path.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "${0##*/}" "$@"\n')
            path.chmod(0o755)
        source = (ROOT / "docker/entrypoint.sh").read_text()
        binding = 'readonly GENTLE_BIN_DIR="/opt/gentle/bin"'
        self.assertEqual(source.count(binding), 1)
        self.script = self.root / "entrypoint.sh"
        # Relocate only the hard-coded installation path, keeping production
        # dispatch intact and avoiding a production binary-directory override.
        self.script.write_text(source.replace(
            binding, f"readonly GENTLE_BIN_DIR={shlex.quote(str(self.bin_dir))}",
        ))

    def invoke(self, *args: str, **extra_env: str) -> subprocess.CompletedProcess:
        env = {key: value for key, value in os.environ.items()
               if not key.startswith(("APPTAINER_", "SINGULARITY_", "GENTLE_"))}
        return subprocess.run(
            ["bash", str(self.script), *args], cwd=self.root,
            env={**env, **extra_env}, text=True, capture_output=True, timeout=10,
        )

    def test_no_arguments_defaults_to_cli_help_in_every_runtime(self) -> None:
        for env in ({}, {"APPTAINER_NAME": "synthetic.sif"}, {"SINGULARITY_NAME": "synthetic.sif"}):
            with self.subTest(env=env):
                result = self.invoke(**env)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.splitlines(), ["gentle_cli", "--help"])

    def test_supported_modes_preserve_arguments_and_clean_stdout(self) -> None:
        for mode, binary in zip(("cli", "mcp", "examples-docs"), BINARIES):
            with self.subTest(mode=mode):
                result = self.invoke(mode, "--state", "/work/project with spaces.json",
                                     "$(touch unexpected)")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.splitlines(), [
                    binary, "--state", "/work/project with spaces.json", "$(touch unexpected)",
                ])
                self.assertEqual(result.stderr, "")
                self.assertFalse((self.root / "unexpected").exists())

    def test_removed_modes_refuse_even_with_legacy_gui_environment(self) -> None:
        for mode in ("gui", "gui-web", "gui-x11", "js", "lua", "gentle", "gentle_js", "gentle_lua"):
            with self.subTest(mode=mode):
                result = self.invoke(mode, GENTLE_CONTAINER_FLAVOR="gui", APPTAINER_NAME="old.sif")
                self.assertEqual(result.returncode, 64)
                self.assertEqual(result.stdout, "")
                self.assertIn("does not include GUI, JavaScript or Lua", result.stderr)

    def test_apptainer_and_singularity_keep_cli_subcommand_dispatch(self) -> None:
        for key in ("APPTAINER_NAME", "APPTAINER_CONTAINER", "SINGULARITY_NAME", "SINGULARITY_CONTAINER"):
            with self.subTest(key=key):
                result = self.invoke("capabilities", "--json", **{key: "synthetic.sif"})
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.splitlines(), ["gentle_cli", "capabilities", "--json"])

    def test_explicit_external_command_remains_available_in_docker(self) -> None:
        result = self.invoke("/usr/bin/printf", "%s", "external command argument")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "external command argument")

    def test_help_does_not_advertise_removed_modes(self) -> None:
        result = self.invoke("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        for mode in ("cli", "mcp", "examples-docs"):
            self.assertIn(f"gentle-image {mode} ", result.stdout)
        for mode in ("gui", "gui-web", "js", "lua"):
            self.assertNotIn(f"gentle-image {mode} ", result.stdout)


class ContainerContractTests(unittest.TestCase):
    def test_builder_copies_embedded_resource_directories_before_compilation(self) -> None:
        build_script = (ROOT / "build.rs").read_text()
        required = re.search(r"let required_files = \[(.*?)\];", build_script, re.S)
        self.assertIsNotNone(required, "Locate the build script's required resource list")
        resources = re.findall(r'"([^"]+)"', required.group(1))
        self.assertTrue(resources)
        docker = (ROOT / "Dockerfile").read_text()
        self.assertIn("RUN cargo build ", docker)
        before_build = docker.split("RUN cargo build ", 1)[0]
        for directory in sorted({Path(path).parts[0] for path in resources}):
            with self.subTest(directory=directory):
                self.assertIn(f"COPY {directory} ./{directory}\n", before_build)
        for resource in resources:
            with self.subTest(resource=resource):
                self.assertTrue((ROOT / resource).is_file())

    def test_build_and_payload_are_headless_without_scripting(self) -> None:
        docker = (ROOT / "Dockerfile").read_text().replace("\\\n", "")
        command = next(line for line in docker.splitlines() if line.startswith("RUN cargo build "))
        tokens = shlex.split(command)
        self.assertIn("--locked", tokens)
        self.assertIn("--no-default-features", tokens)
        self.assertIn("-j1", tokens)
        self.assertNotIn("--features", tokens)
        self.assertNotIn("--all-features", tokens)
        self.assertNotIn("--bins", tokens)
        self.assertEqual([tokens[i + 1] for i, token in enumerate(tokens) if token == "--bin"], BINARIES)
        self.assertNotIn("runtime-gui", docker)
        self.assertNotIn("gentle-dist-gui", docker)
        for package in ("libgtk-3-dev", "libx11-dev", "novnc", "xvfb", "openbox"):
            self.assertNotIn(package, docker)
        self.assertIn("cp -a assets /opt/gentle-dist-cli/assets", docker)
        self.assertIn("fonts-dejavu-core", docker)
        self.assertIn("integrations/python", docker)
        for binary in BINARIES:
            self.assertIn(f"/opt/gentle-dist-cli/bin/{binary}", docker)

    def test_dependency_guard_recognizes_desktop_and_scripting_packages(self) -> None:
        docker = (ROOT / "Dockerfile").read_text()
        self.assertIn("cargo tree --locked --no-default-features --edges normal,build", docker)
        pattern = re.search(r"grep -E '([^']+)' /tmp/gentle-dependencies.txt", docker).group(1)
        for package in ("eframe", "egui", "gentle-gui", "winit", "arboard", "rfd",
                        "deno_core", "deno_error", "v8", "mlua", "mlua-sys", "lua-src", "luajit-src"):
            with self.subTest(package=package):
                self.assertRegex(f"{package} v1.2.3", pattern)
        for package in ("gentle-engine", "gentle-render", "resvg", "svg", "serde"):
            self.assertNotRegex(f"{package} v1.2.3", pattern)

        command = next(line for line in docker.replace("\\\n", "").splitlines()
                       if line.startswith("RUN cargo tree "))
        guard = command.split("&&", 1)[1]
        with TemporaryDirectory() as folder:
            dependencies = Path(folder) / "dependencies.txt"
            guard = guard.replace("/tmp/gentle-dependencies.txt", shlex.quote(str(dependencies)))
            for extra in ("", "eframe v1.2.3\n", "v8 v1.2.3\n", "mlua v1.2.3\n"):
                with self.subTest(dependency=extra):
                    dependencies.write_text("gentle-engine v1.2.3\n" + extra)
                    result = subprocess.run(["sh", "-c", guard], capture_output=True,
                                            text=True, timeout=10)
                    self.assertEqual(result.returncode, 1 if extra else 0, result.stderr)
                    if extra:
                        self.assertIn("dependency leaked", result.stderr)

    def test_workflow_builds_and_publishes_only_the_headless_target(self) -> None:
        workflow = (ROOT / ".github/workflows/container.yml").read_text()
        self.assertEqual(workflow.count("target: runtime-cli"), 2)
        for old in ("runtime-gui", "check-gui", "push-gui", "GUI_IMAGE", "script-interfaces"):
            self.assertNotIn(old, workflow)
        self.assertIn("--network none", workflow)
        self.assertIn('ldd "/opt/gentle/bin/$bin"', workflow)
        self.assertIn("for bin in gentle gentle_js gentle_lua", workflow)
        self.assertIn("test -s /opt/gentle/assets/genomes.json", workflow)
        self.assertIn('python3 -c "import gentle_py"', workflow)
        for tag in ("cli", "${{ needs.candidate.outputs.tag }}-cli", "${{ needs.candidate.outputs.tag }}", "latest"):
            self.assertIn(f"type=raw,value={tag}\n", workflow)
        self.assertIn("python3 -m unittest scripts.test_container -v",
                      (ROOT / ".github/workflows/ci.yml").read_text())

    def test_receipt_records_exact_headless_features_binaries_and_identity(self) -> None:
        workflow = (ROOT / ".github/workflows/container.yml").read_text()
        script = textwrap.dedent(re.search(
            r"          python3 - <<'PY'\n(.*?)          PY\n", workflow, re.S,
        ).group(1))
        env = {
            "CLI_IMAGE": "sha256:" + "a" * 64, "CLI_DIGEST": "sha256:" + "b" * 64,
            "RELEASE_TAG": "v0.1.0-internal.10", "EXPECTED_REVISION": "c" * 40,
            "EXPECTED_LOCK_SHA256": "d" * 64, "WORKFLOW_REVISION": "e" * 40,
            "CANDIDATE_MODE": "validate_only",
        }
        records = []
        with patch.dict(os.environ, env, clear=True), \
             patch("subprocess.check_output", return_value="synthetic-buildx\n"), \
             patch.object(Path, "write_text", side_effect=lambda text: records.append(json.loads(text))):
            exec(compile(script, "container-receipt", "exec"), {})
        record, = records
        self.assertFalse(record["default_features"])
        self.assertEqual(record["features"], [])
        self.assertEqual(record["binaries"], BINARIES)
        self.assertFalse(record["pushed"])
        self.assertEqual(record["revision"], env["EXPECTED_REVISION"])
        self.assertEqual(record["cargo_lock_sha256"], env["EXPECTED_LOCK_SHA256"])
        self.assertEqual(record["images"], [{
            "target": "runtime-cli", "image_id": env["CLI_IMAGE"],
            "digest": env["CLI_DIGEST"], "entrypoint_smoke": "passed",
        }])


if __name__ == "__main__":
    unittest.main()
