#!/usr/bin/env python3
"""Offline headless-container regressions with synthetic executable/receipt data.

Temporary shell stubs exercise the real entrypoint without Docker, GUI, scientific
inputs or network access. Run: python3 -m unittest scripts.test_container -v
Actual Linux linking, dependencies and image execution remain container CI gates.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
from tempfile import TemporaryDirectory
import textwrap
import tomllib
import unittest
from unittest.mock import patch

from scripts import release_candidate as policy


ROOT = Path(__file__).resolve().parents[1]
BINARIES = ["gentle_cli", "gentle_mcp", "gentle_examples_docs"]


def synthetic_container_outputs(image_id: str, revision: str, tag: str) -> tuple[dict, dict]:
    records, outputs = {}, {}
    for index, binary in enumerate(BINARIES, start=1):
        path = f"/opt/gentle/bin/{binary}"
        version = f"GENtle {tag[1:]}\nBuild synthetic\nSource revision {tag[1:]}+git.{revision}\n"
        digest = f"{index:064x}"
        records[binary] = {"version": version.strip(), "sha256": digest}
        outputs[("docker", "run", "--rm", "--network", "none", "--entrypoint",
                 path, image_id, "--version")] = version
        outputs[("docker", "run", "--rm", "--network", "none", "--entrypoint",
                 "/usr/bin/sha256sum", image_id, path)] = f"{digest}  {path}\n"
    return records, outputs


class ContainerBinaryIdentityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.image_id = "sha256:" + "a" * 64
        self.revision = "c" * 40
        self.tag = "v0.1.0-internal.12"
        self.records, self.outputs = synthetic_container_outputs(
            self.image_id, self.revision, self.tag,
        )

    def verify(self) -> dict:
        return policy.container_binary_identities(self.image_id, self.revision, self.tag)

    def test_all_three_actual_commands_use_the_loaded_image_without_network(self) -> None:
        with patch("subprocess.check_output", side_effect=lambda command, **_: self.outputs[tuple(command)]) as run:
            self.assertEqual(self.verify(), self.records)
        self.assertEqual([call.args[0] for call in run.call_args_list],
                         [list(command) for command in self.outputs])
        for call in run.call_args_list:
            self.assertEqual(call.kwargs, {"text": True, "timeout": 60})

    def test_missing_stale_wrong_version_or_duplicate_source_revision_refuses(self) -> None:
        expected = f"Source revision {self.tag[1:]}+git.{self.revision}"
        for binary in BINARIES:
            command = next(command for command in self.outputs
                           if command[-1] == "--version" and command[-3].endswith(f"/{binary}"))
            for source in ("", "Source revision unknown", f"Source revision {self.tag[1:]}",
                           f"Source revision {self.tag[1:]}+git.{'d' * 40}",
                           f"Source revision 0.1.0-internal.11+git.{self.revision}",
                           f"{expected}\n{expected}"):
                with self.subTest(binary=binary, source=source):
                    outputs = {**self.outputs, command: f"GENtle synthetic\n{source}\n"}
                    with patch("subprocess.check_output", side_effect=lambda args, **_: outputs[tuple(args)]):
                        with self.assertRaisesRegex(ValueError, binary):
                            self.verify()

    def test_malformed_digest_or_wrong_binary_path_refuses(self) -> None:
        command = list(self.outputs)[-1]
        path = command[-1]
        for digest in ("", f"{'a' * 63}  {path}", f"{'A' * 64}  {path}",
                       f"{'a' * 64}  /opt/gentle/bin/gentle_cli",
                       f"{'a' * 64}  {path}\n{'a' * 64}  {path}"):
            with self.subTest(digest=digest):
                outputs = {**self.outputs, command: digest}
                with patch("subprocess.check_output", side_effect=lambda args, **_: outputs[tuple(args)]):
                    with self.assertRaisesRegex(ValueError, "gentle_examples_docs"):
                        self.verify()

    def test_failed_binary_or_hash_command_does_not_return_partial_identity(self) -> None:
        for failed in list(self.outputs)[-2:]:
            def execute(command, **_):
                if tuple(command) == failed:
                    raise subprocess.CalledProcessError(1, command)
                return self.outputs[tuple(command)]
            with self.subTest(command=failed):
                with patch("subprocess.check_output", side_effect=execute):
                    with self.assertRaises(subprocess.CalledProcessError):
                        self.verify()

    def test_timeout_is_a_failure_not_an_identity(self) -> None:
        with patch("subprocess.check_output", side_effect=subprocess.TimeoutExpired("docker", 60)):
            with self.assertRaises(subprocess.TimeoutExpired):
                self.verify()

    def test_invalid_image_revision_or_label_refuses_before_docker(self) -> None:
        inputs = [(image, self.revision, self.tag) for image in (
            "build-check-cli", "sha256:short", "sha256:" + "A" * 64,
            "sha256:" + "a" * 64 + " --privileged",
        )]
        inputs += [(self.image_id, revision, self.tag) for revision in ("main", "c" * 7, "C" * 40)]
        inputs += [(self.image_id, self.revision, tag) for tag in ("", "0.1.0-internal.12", "v1\n")]
        for args in inputs:
            with self.subTest(args=args), patch("subprocess.check_output") as run:
                with self.assertRaises(ValueError):
                    policy.container_binary_identities(*args)
                run.assert_not_called()


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
    def test_candidate_revision_is_passed_only_to_the_gentle_build(self) -> None:
        docker = (ROOT / "Dockerfile").read_text()
        binding = 'ARG GENTLE_GIT_COMMIT=""\nENV GENTLE_GIT_COMMIT=${GENTLE_GIT_COMMIT}\n'
        self.assertEqual(docker.count(binding), 1)
        self.assertLess(docker.index("cargo install --locked --debug"), docker.index(binding))
        self.assertLess(docker.index("RUN cargo tree "), docker.index(binding))
        self.assertLess(docker.index(binding), docker.index("RUN cargo build "))
        workflow = (ROOT / ".github/workflows/container.yml").read_text()
        build_args = re.findall(r"          build-args: \|\n((?:            .*\n)+)", workflow)
        self.assertEqual(len(build_args), 2)
        for args in build_args:
            self.assertIn("GENTLE_GIT_COMMIT=${{ needs.candidate.outputs.revision }}\n", args)

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
        self.assertEqual(tokens[tokens.index("--profile") + 1], "package-opt1")
        self.assertNotIn("--release", tokens)
        self.assertNotIn("GENTLE_CARGO_PROFILE", docker)
        self.assertIn("CARGO_INCREMENTAL=0", docker)
        self.assertIn("CARGO_PROFILE_DEV_DEBUG=0", docker)
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
            self.assertIn(f'install -Dm755 "target/package-opt1/{binary}"', docker)

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

    def test_rnapkin_has_build_and_runtime_fonts_and_installs_before_gentle(self) -> None:
        docker = (ROOT / "Dockerfile").read_text().replace("\\\n", "")
        builder, runtime = docker.split("FROM debian:${DEBIAN_SUITE}-slim AS runtime-cli", 1)
        for stage, packages in (
            (builder, ("pkg-config", "curl", "libfontconfig1-dev", "libfreetype6-dev",
                       "fonts-dejavu-core")),
            (runtime, ("libfontconfig1", "libfreetype6", "fonts-dejavu-core")),
        ):
            install = next(line for line in stage.splitlines() if "apt-get install" in line)
            for package in packages:
                with self.subTest(package=package):
                    self.assertIn(package, shlex.split(install))
        command = re.search(r"cargo install [^&\n]+", builder).group(0).strip()
        self.assertEqual(shlex.split(command), [
            "cargo", "install", "--locked", "--debug", "--path", "/opt/rnapkin-src",
            "--root", "/opt/rnapkin", "-j1",
        ])
        self.assertIn("https://static.crates.io/crates/rnapkin/rnapkin-0.3.9.crate", builder)
        self.assertIn(
            'echo "4495690197e1cced9b16234d6a66b40ebf190e9613b9c1c6aea837adeed00f17  '
            '/tmp/rnapkin.crate" | sha256sum -c -', builder,
        )
        self.assertLess(builder.index("sha256sum -c -"), builder.index("tar -xzf"))
        self.assertIn("COPY docker/rnapkin/Cargo.lock /tmp/rnapkin.Cargo.lock", builder)
        self.assertLess(builder.index("cp /tmp/rnapkin.Cargo.lock /opt/rnapkin-src/Cargo.lock"),
                        builder.index(command))
        self.assertLess(builder.index(command), builder.index("COPY Cargo.toml"))
        self.assertLess(builder.index(command), builder.index("RUN cargo build "))
        self.assertIn("COPY --from=build /opt/rnapkin/bin/rnapkin /usr/local/bin/rnapkin", runtime)

    def test_rnapkin_lock_keeps_the_fixed_bitmap_backend_and_existing_application(self) -> None:
        lock = tomllib.loads((ROOT / "docker/rnapkin/Cargo.lock").read_text())
        packages = {package["name"]: package for package in lock["package"]}
        for name, version in (("rnapkin", "0.3.9"), ("plotters", "0.3.4"),
                              ("plotters-bitmap", "0.3.3"), ("plotters-backend", "0.3.7"),
                              ("gif", "0.12.0")):
            with self.subTest(package=name):
                self.assertEqual(packages[name]["version"], version)
        self.assertEqual(packages["plotters-bitmap"]["checksum"],
                         "0cebbe1f70205299abc69e8b295035bb52a6a70ee35474ad10011f0a4efb8543")
        self.assertEqual(packages["plotters-bitmap"]["source"],
                         "registry+https://github.com/rust-lang/crates.io-index")

    def test_rnapkin_smoke_precedes_gentle_build_without_disabling_checks(self) -> None:
        docker = (ROOT / "Dockerfile").read_text()
        before_gentle = docker.split("COPY Cargo.toml", 1)[0]
        self.assertIn('printf "%s\\n" "GGGAAACCC" "(((...)))"', before_gentle)
        for extension in ("svg", "png"):
            self.assertIn(
                f'timeout 30 /opt/rnapkin/bin/rnapkin --height 128 '
                f'-o "$smoke_dir/hairpin.{extension}" "$smoke_dir/hairpin.dbn"',
                before_gentle,
            )
            self.assertIn(f'test -s "$smoke_dir/hairpin.{extension}"', before_gentle)
        self.assertNotIn("CARGO_PROFILE_DEV_DEBUG_ASSERTIONS", docker)
        self.assertNotIn("CARGO_PROFILE_DEV_OVERFLOW_CHECKS", docker)
        self.assertNotIn("debug-assertions=off", docker)

    def test_container_smoke_renders_rna_with_fonts_without_network(self) -> None:
        workflow = (ROOT / ".github/workflows/container.yml").read_text()
        smoke = workflow.split("      - name: Smoke headless container without network access\n", 1)[1]
        smoke = smoke.split("      - name: Retain container build identity\n", 1)[0]
        self.assertIn('--network none --entrypoint /bin/sh "$CLI_IMAGE" -ec', smoke)
        self.assertIn("ldd /usr/local/bin/rnapkin", smoke)
        self.assertIn("rnapkin --version", smoke)
        self.assertIn('printf "%s\\n" "GGGAAACCC" "(((...)))" > hairpin.dbn', smoke)
        for extension in ("svg", "png"):
            with self.subTest(extension=extension):
                self.assertIn(f"timeout 30 rnapkin --height 128 -o hairpin.{extension} hairpin.dbn", smoke)
                self.assertIn(f"test -s hairpin.{extension}", smoke)

    def test_workflow_builds_and_publishes_only_the_headless_target(self) -> None:
        workflow = (ROOT / ".github/workflows/container.yml").read_text()
        verify = workflow.split("- name: Verify immutable candidate identity\n", 1)[1].split("\n      - name:", 1)[0]
        self.assertIn('test "$NATIVE_PROFILE:$NATIVE_TARGET_SUBDIR" = "package-opt1:package-opt1"', verify)
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
        helper_lock = b"# Synthetic helper lockfile\nversion = 4\n"
        identities, outputs = synthetic_container_outputs(
            env["CLI_IMAGE"], env["EXPECTED_REVISION"], env["RELEASE_TAG"],
        )
        outputs[("docker", "buildx", "version")] = "synthetic-buildx\n"
        with patch.dict(os.environ, env, clear=True), \
             patch("subprocess.check_output", side_effect=lambda args, **_: outputs[tuple(args)]) as run, \
             patch.object(Path, "read_bytes", return_value=helper_lock) as read_bytes, \
             patch.object(Path, "write_text", side_effect=lambda text: records.append(json.loads(text))):
            exec(compile(script, "container-receipt", "exec"), {})
        read_bytes.assert_called_once()
        record, = records
        self.assertEqual(record["profile"], "package-opt1")
        self.assertEqual(record["opt_level"], 1)
        self.assertEqual(record["lto"], "off")
        self.assertEqual(record["codegen_units"], 256)
        self.assertIs(record["debug_assertions"], True)
        self.assertIs(record["overflow_checks"], True)
        self.assertEqual(record["panic"], "unwind")
        self.assertIs(record["incremental"], False)
        self.assertEqual(record["debug"], 0)
        self.assertFalse(record["default_features"])
        self.assertEqual(record["features"], [])
        self.assertEqual(record["binaries"], BINARIES)
        self.assertFalse(record["pushed"])
        self.assertEqual(record["build_args"], {
            "DEBIAN_SUITE": "forky", "GENTLE_GIT_COMMIT": env["EXPECTED_REVISION"],
        })
        self.assertEqual(record["revision"], env["EXPECTED_REVISION"])
        self.assertEqual(record["cargo_lock_sha256"], env["EXPECTED_LOCK_SHA256"])
        self.assertEqual(record["helpers"], {
            "rnapkin": {"version": "0.3.9", "profile": "dev",
                        "cargo_lock_path": "docker/rnapkin/Cargo.lock",
                        "cargo_lock_sha256": hashlib.sha256(helper_lock).hexdigest()},
        })
        self.assertEqual(record["images"], [{
            "target": "runtime-cli", "image_id": env["CLI_IMAGE"],
            "digest": env["CLI_DIGEST"], "entrypoint_smoke": "passed",
            "binary_identities": identities,
        }])
        self.assertEqual([call.args[0] for call in run.call_args_list],
                         [list(command) for command in outputs])

    def test_identity_refusal_prevents_writing_a_container_receipt(self) -> None:
        workflow = (ROOT / ".github/workflows/container.yml").read_text()
        script = textwrap.dedent(re.search(
            r"          python3 - <<'PY'\n(.*?)          PY\n", workflow, re.S,
        ).group(1))
        with patch.dict(os.environ, {"CLI_IMAGE": "synthetic", "EXPECTED_REVISION": "c" * 40,
                                     "RELEASE_TAG": "v0.1.0-internal.12"}, clear=True), \
             patch.object(policy, "container_binary_identities", side_effect=ValueError("unbound binary")) as verify, \
             patch.object(Path, "write_text") as write:
            with self.assertRaisesRegex(ValueError, "unbound binary"):
                exec(compile(script, "container-receipt", "exec"), {})
            verify.assert_called_once_with("synthetic", "c" * 40, "v0.1.0-internal.12")
            write.assert_not_called()


if __name__ == "__main__":
    unittest.main()
