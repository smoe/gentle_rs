"""Offline package-contract tests using disposable, hand-crafted Git repositories.

All resources, app shells and executable stand-ins are synthetic and recreated
by setUp; they test packaging, never native GENtle or scientific acceptance.
Run with: python3 -m unittest scripts.test_package_desktop -v
"""

from __future__ import annotations

import os
from pathlib import Path
import plistlib
import shlex
import shutil
import subprocess
import tempfile
import textwrap
import tomllib
import unittest
from unittest.mock import patch

from scripts import package_desktop as package


class DesktopPackageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="gentle package ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.repo = self.root / "checkout"
        self.repo.mkdir()
        self.git("init", "-q")
        self.git("config", "user.name", "Synthetic packaging test")
        self.git("config", "user.email", "package-test@example.invalid")
        for relative, content in {
            "assets/genomes.json": "{}",
            "icons/icon.png": "synthetic icon",
            "docs/tutorial/manifest.json": '{"chapters": []}',
            "docs/linux_tarball.md": "Synthetic quick start",
            "test_files/tiny.fa": ">synthetic\nACGT\n",
            "integrations/python/gentle_py/__init__.py": "# synthetic adapter\n",
            "data/resources/affymetrix/platform_registry.json": "{}",
            "README.md": "Synthetic distribution",
            "Credits.rtf": "Synthetic credits",
            "unrelated.txt": "Tracked but outside the distribution inventory",
        }.items():
            target = self.repo / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content, encoding="utf-8")
        self.git("add", ".")
        self.git("commit", "-qm", "synthetic distribution")
        self.revision = self.git("rev-parse", "HEAD").strip()
        (self.repo / "assets/private-cache.txt").write_text("must not be packaged")
        self.binaries = self.root / "release"
        self.binaries.mkdir()
        for name in package.BINARIES:
            for suffix in ("", ".exe"):
                path = self.binaries / (name + suffix)
                path.write_text("synthetic executable\n")
                path.chmod(0o755)
        self.app = self.root / "cargo bundle/GENtle.app"
        (self.app / "Contents/Resources").mkdir(parents=True)
        (self.app / "Contents/Info.plist").write_bytes(plistlib.dumps({
            "CFBundleExecutable": "gentle", "CFBundleIconFile": "GENtle.icns",
        }))
        (self.app / "Contents/Resources/GENtle.icns").write_bytes(b"synthetic generated icon")
        (self.app / "Contents/Resources/private-cache.txt").write_text("not tracked")

    def git(self, *args: str) -> str:
        return subprocess.check_output(["git", *args], cwd=self.repo, text=True)

    def stage(self, platform: str) -> Path:
        root = self.root / ("GENtle.app" if platform == "macos" else platform)
        real_run = subprocess.run

        def run(args, **kwargs):
            if args[0] == "ldd":
                return subprocess.CompletedProcess(args, 0, "libc => /system/libc\n", "")
            return real_run(args, **kwargs)

        with patch.object(package.subprocess, "run", side_effect=run):
            package.stage(self.repo, self.binaries, root, platform, "v0.1.0-test", self.app)
        return root

    def smoke(self, root: Path, platform: str) -> None:
        package.smoke(root, platform, self.revision, self.repo)

    def test_all_layouts_include_entrypoints_and_only_tracked_resources(self) -> None:
        for platform in ("macos", "windows", "linux"):
            with self.subTest(platform=platform):
                root = self.stage(platform)
                binary_root, resource_root = package.layout(root, platform)
                suffix = ".exe" if platform == "windows" else ""
                for name in package.BINARIES:
                    self.assertTrue((binary_root / (name + suffix)).is_file())
                self.assertTrue((resource_root / "docs/tutorial/manifest.json").is_file())
                self.assertTrue((resource_root / "test_files/tiny.fa").is_file())
                self.assertFalse((resource_root / "assets/private-cache.txt").exists())
                self.assertFalse((resource_root / "private-cache.txt").exists())
                self.assertFalse((resource_root / "unrelated.txt").exists())
                if platform == "macos":
                    self.assertEqual((resource_root / "GENtle.icns").read_bytes(),
                                     b"synthetic generated icon")
                with patch.object(package.subprocess, "run", return_value=
                                  subprocess.CompletedProcess([], 0, "ok", "")) as run:
                    self.smoke(root, platform)
                self.assertEqual(run.call_count, 6)
                calls = [call.args[0] for call in run.call_args_list]
                self.assertIn([str(binary_root / ("gentle_cli" + suffix)), "capabilities"], calls)
                self.assertIn([str(binary_root / ("gentle_examples_docs" + suffix)),
                               "tutorial-manifest-check"], calls)
                for call in run.call_args_list:
                    self.assertEqual(call.kwargs["cwd"], resource_root)
                    self.assertFalse(resource_root.is_relative_to(self.repo))
                    self.assertEqual(call.kwargs["timeout"], 120)

    def test_staging_excludes_scripting_and_reproduction_binaries(self) -> None:
        extras = ("gentle_js", "gentle_lua", "gentle_dna_viewer_repro", "gentle_egui_window_repro")
        bundle_binaries = self.app / "Contents/MacOS"
        bundle_binaries.mkdir()
        for name in (*package.BINARIES, *extras):
            (bundle_binaries / name).write_text("stale cargo-bundle executable\n")
        for name in extras:
            for suffix in ("", ".exe"):
                (self.binaries / (name + suffix)).write_text("unrequested executable\n")
        for platform in ("macos", "windows", "linux"):
            with self.subTest(platform=platform):
                root = self.stage(platform)
                binary_root, _ = package.layout(root, platform)
                suffix = ".exe" if platform == "windows" else ""
                self.assertEqual({path.name for path in binary_root.iterdir()},
                                 {name + suffix for name in package.BINARIES})
                for name in package.BINARIES:
                    self.assertEqual((binary_root / (name + suffix)).read_bytes(),
                                     (self.binaries / (name + suffix)).read_bytes())

    def test_smoke_rejects_scripting_even_when_recorded_in_checksums(self) -> None:
        for platform in ("macos", "windows", "linux"):
            root = self.stage(platform)
            binary_root, _ = package.layout(root, platform)
            suffix = ".exe" if platform == "windows" else ""
            for name in ("gentle_js", "gentle_lua"):
                with self.subTest(platform=platform, binary=name):
                    extra = binary_root / (name + suffix)
                    extra.write_text("unexpected scripting executable\n")
                    (root / "SHA256SUMS").write_text("".join(
                        f"{package.digest(path)}  {path.relative_to(root).as_posix()}\n"
                        for path in package.files(root)
                    ), encoding="utf-8")
                    with patch.object(package.subprocess, "run") as run:
                        with self.assertRaisesRegex(ValueError, "Unexpected packaged scripting binary"):
                            self.smoke(root, platform)
                        run.assert_not_called()
                    extra.unlink()

    def test_zip_and_tar_survive_extraction_with_complete_inventory(self) -> None:
        for platform, archive_format in (("windows", "zip"), ("linux", "gztar")):
            root = self.stage(platform)
            archive = shutil.make_archive(str(self.root / f"archive-{platform}"),
                                          archive_format, root.parent, root.name)
            extracted = self.root / f"extracted-{platform}"
            shutil.unpack_archive(archive, extracted)
            with patch.object(package.subprocess, "run", return_value=
                              subprocess.CompletedProcess([], 0, "ok", "")):
                self.smoke(extracted / root.name, platform)

    @unittest.skipIf(os.name == "nt", "POSIX executable stand-ins; native Windows smoke runs in CI")
    def test_smoke_executes_relocated_entrypoints_without_original_checkout(self) -> None:
        for name in package.BINARIES:
            (self.binaries / name).write_text(
                "#!/usr/bin/env python3\n"
                "import json, sys\nfrom pathlib import Path\n"
                "assert Path('assets/genomes.json').is_file()\n"
                "assert json.loads(Path('docs/tutorial/manifest.json').read_text()) == {'chapters': []}\n"
                "print(sys.argv)\n"
            )
        root = self.stage("macos")
        relocated = self.root / "mounted/GENtle.app"
        shutil.copytree(root, relocated)
        shutil.rmtree(self.repo)
        shutil.rmtree(self.binaries)
        self.smoke(relocated, "macos")

    def test_staging_rejects_missing_binary_before_creating_package(self) -> None:
        (self.binaries / "gentle_mcp.exe").unlink()
        with self.assertRaisesRegex(ValueError, "Missing or empty release binary"):
            self.stage("windows")
        self.assertFalse((self.root / "windows").exists())

    def test_staging_never_overwrites_existing_destination(self) -> None:
        self.stage("macos")
        with self.assertRaisesRegex(ValueError, "destination must be new"):
            self.stage("macos")

    def test_macos_declared_icon_must_exist_inside_resources(self) -> None:
        for icon in ("missing.icns", "../Info.plist"):
            with self.subTest(icon=icon):
                (self.app / "Contents/Info.plist").write_bytes(
                    plistlib.dumps({"CFBundleIconFile": icon}))
                with self.assertRaisesRegex(ValueError, "bundle icon|resource filename"):
                    self.stage("macos")
                self.assertFalse((self.root / "GENtle.app").exists())

    @unittest.skipIf(os.name == "nt", "Symlink creation needs privileges on Windows")
    def test_staging_rejects_tracked_resource_escaping_checkout(self) -> None:
        outside = self.root / "private.txt"
        outside.write_text("private")
        (self.repo / "assets/escape").symlink_to(outside)
        self.git("add", "assets/escape")
        with self.assertRaisesRegex(ValueError, "escaping tracked resource"):
            self.stage("windows")

    def test_missing_manifest_fails_before_executing_any_binary(self) -> None:
        self.git("rm", "-q", "docs/tutorial/manifest.json")
        root = self.stage("windows")
        with patch.object(package.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "Missing packaged resource"):
                self.smoke(root, "windows")
            run.assert_not_called()

    def test_tampered_or_extra_files_and_wrong_revision_fail_closed(self) -> None:
        root = self.stage("windows")
        with self.assertRaisesRegex(ValueError, "revision does not match"):
            package.smoke(root, "windows", "0" * 40, self.repo)
        extra = root / "unrecorded.txt"
        extra.write_text("unexpected")
        with self.assertRaisesRegex(ValueError, "file inventory differs"):
            self.smoke(root, "windows")
        extra.unlink()
        (root / "assets/genomes.json").write_text("tampered")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.smoke(root, "windows")

    def test_checkout_cannot_masquerade_as_extracted_package(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside the checkout"):
            self.smoke(self.repo / "dist", "windows")

    def test_failing_entrypoint_is_not_swallowed(self) -> None:
        root = self.stage("windows")
        with patch.object(package.subprocess, "run", return_value=
                          subprocess.CompletedProcess([], 9, "", "missing DLL")):
            with self.assertRaisesRegex(ValueError, "gentle --version failed.*9"):
                self.smoke(root, "windows")

    def test_stale_manifest_failure_is_not_swallowed(self) -> None:
        root = self.stage("windows")

        def run(args, **kwargs):
            code = 1 if args[1] == "tutorial-manifest-check" else 0
            return subprocess.CompletedProcess(args, code, "", "stale manifest" if code else "")

        with patch.object(package.subprocess, "run", side_effect=run):
            with self.assertRaisesRegex(ValueError, "tutorial-manifest-check failed"):
                self.smoke(root, "windows")


class WorkflowWiringTests(unittest.TestCase):
    def test_native_release_disables_all_lto_without_changing_custom_profiles(self) -> None:
        repo = Path(__file__).resolve().parents[1]
        manifest = tomllib.loads((repo / "Cargo.toml").read_text())
        profiles = manifest["profile"]
        self.assertEqual(profiles["release"], {"lto": "off"},
                         "Disable all native LTO; retain other Cargo defaults")
        self.assertEqual(profiles["release-fast"], {
            "inherits": "release", "lto": "thin", "codegen-units": 16,
            "opt-level": 2, "panic": "abort", "strip": True,
        })
        self.assertEqual(profiles["bench-audit"], {
            "inherits": "release", "lto": "thin", "codegen-units": 16,
            "panic": "unwind", "strip": True,
        })

    def test_native_release_builds_only_the_five_packaged_binaries(self) -> None:
        self.assertEqual(package.BINARIES, (
            "gentle", "gentle_cli", "gentle_mcp",
            "gentle_examples_docs", "gentle_publication_report",
        ))
        repo = Path(__file__).resolve().parents[1]
        workflow = (repo / ".github/workflows/release.yml").read_text()
        build = workflow.split("- name: Build locked release binaries\n", 1)[1].split("\n      - name:", 1)[0]
        tokens = shlex.split(build.split("build=(", 1)[1].split(")", 1)[0])
        self.assertEqual(tokens[:2], ["cargo", "build"])
        self.assertIn("--locked", tokens)
        self.assertIn("--release", tokens)
        self.assertIn("-j1", tokens)
        self.assertEqual([tokens[i + 1] for i, token in enumerate(tokens) if token == "--bin"],
                         list(package.BINARIES))
        for flag in ("--bins", "--all-targets", "--workspace", "--features", "--all-features", "--no-default-features"):
            self.assertNotIn(flag, tokens)
        self.assertIn("run: cargo bundle --release --bin gentle --format", workflow)
        self.assertNotIn("script-interfaces", workflow)
        self.assertNotIn("CARGO_PROFILE_RELEASE_", workflow)
        self.assertNotIn("RUSTFLAGS:", workflow)
        self.assertIn('"features": []', workflow)
        self.assertIn('"default_features": True, "binaries": list(BINARIES)', workflow)
        self.assertIn("from scripts.package_desktop import BINARIES", workflow)
        for name in package.BINARIES:
            self.assertIn(f'/release/{name}${{suffix}}"', workflow)
        for name in ("gentle_js", "gentle_lua"):
            self.assertNotIn(f'/release/{name}${{suffix}}"', workflow)

        manifest = tomllib.loads((repo / "Cargo.toml").read_text())
        self.assertEqual(set(manifest["features"]["default"]), {"desktop-gui", "screenshot-capture"})
        for name, feature in (("gentle_js", "js-interface"), ("gentle_lua", "lua-interface")):
            binary = next(binary for binary in manifest["bin"] if binary["name"] == name)
            self.assertEqual(binary["required-features"], [feature])
            self.assertIn(feature, manifest["features"]["script-interfaces"])

    def test_release_build_retains_diagnostics_even_on_failure(self) -> None:
        repo = Path(__file__).resolve().parents[1]
        workflow = (repo / ".github/workflows/release.yml").read_text()
        build = workflow.split("- name: Build locked release binaries\n", 1)[1].split("\n      - name:", 1)[0]
        self.assertIn("shell: bash", build)
        self.assertIn("set -euo pipefail", build)
        self.assertIn('build=(/usr/bin/time -v "${build[@]}")', build)
        self.assertIn('build=(/usr/bin/time -l "${build[@]}")', build)
        self.assertIn('"${build[@]}" 2>&1 | tee -a "$log"', build)
        upload = workflow.split("- name: Retain release build diagnostics\n", 1)[1].split("\n      - name:", 1)[0]
        self.assertIn("if: ${{ always() }}", upload)
        self.assertIn("uses: actions/upload-artifact@v6", upload)
        self.assertIn("${{ runner.temp }}/gentle-release-build.log", upload)
        self.assertIn("${{ matrix.platform }}", upload)
        self.assertIn("${{ github.run_attempt }}", upload)

    @unittest.skipUnless(os.name == "posix" and shutil.which("bash"),
                         "Synthetic executable stand-ins require a POSIX host and Bash")
    def test_release_build_logging_preserves_cargo_exit_status(self) -> None:
        repo = Path(__file__).resolve().parents[1]
        workflow = (repo / ".github/workflows/release.yml").read_text()
        build = workflow.split("- name: Build locked release binaries\n", 1)[1].split("\n      - name:", 1)[0]
        script = textwrap.dedent(build.split("run: |\n", 1)[1])
        with tempfile.TemporaryDirectory(prefix="synthetic release log ") as directory:
            root = Path(directory)
            # No real compiler runs: these stand-ins reproduce the process boundary.
            (root / "rustc").write_text("#!/bin/sh\nprintf 'synthetic rustc\\n'\n")
            (root / "cargo").write_text(
                "#!/bin/sh\n"
                "if [ \"$1\" = --version ]; then printf 'synthetic cargo\\n'; exit 0; fi\n"
                "printf 'synthetic build stdout\\n'\n"
                "printf 'synthetic compiler stderr\\n' >&2\n"
                "exit \"$SYNTHETIC_CARGO_EXIT\"\n"
            )
            for name in ("cargo", "rustc"):
                (root / name).chmod(0o755)
            for code in (0, 101, 143):
                with self.subTest(exit_code=code):
                    result = subprocess.run(
                        [shutil.which("bash"), "-c", script], cwd=root,
                        env={**os.environ, "PATH": str(root) + os.pathsep + os.environ["PATH"],
                             "RUNNER_TEMP": str(root), "RUNNER_OS": "Windows",
                             "EXPECTED_REVISION": "synthetic-revision",
                             "SYNTHETIC_CARGO_EXIT": str(code)},
                        capture_output=True, text=True, timeout=30,
                    )
                    self.assertEqual(result.returncode, code, result.stdout + result.stderr)
                    log = (root / "gentle-release-build.log").read_text()
                    self.assertIn("synthetic build stdout", log)
                    self.assertIn("synthetic compiler stderr", log)
                    self.assertIn("revision=synthetic-revision platform=Windows", log)
                    self.assertEqual(result.stdout, log)

    def test_release_stages_and_checks_each_platform_after_extraction(self) -> None:
        repo = Path(__file__).resolve().parents[1]
        workflow = (repo / ".github/workflows/release.yml").read_text()
        for platform in ("macos", "windows", "linux"):
            self.assertIn(f"scripts/package_desktop.py stage --platform {platform}", workflow)
            self.assertIn(f"scripts/package_desktop.py smoke --platform {platform}", workflow)
        self.assertIn("trap 'hdiutil detach", workflow)
        self.assertIn('if ($LASTEXITCODE -ne 0)', workflow)
        self.assertIn("scripts.test_package_desktop", (repo / ".github/workflows/ci.yml").read_text())


if __name__ == "__main__":
    unittest.main()
