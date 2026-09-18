"""Offline package-contract tests using disposable, hand-crafted Git repositories.

All resources, app shells and executable stand-ins are synthetic and recreated
by setUp; they test packaging, never native GENtle or scientific acceptance.
Run with: python3 -m unittest scripts.test_package_desktop -v
"""

from __future__ import annotations

import os
from pathlib import Path
import plistlib
import shutil
import subprocess
import tempfile
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
                self.assertEqual(run.call_count, 8)
                calls = [call.args[0] for call in run.call_args_list]
                self.assertIn([str(binary_root / ("gentle_cli" + suffix)), "capabilities"], calls)
                self.assertIn([str(binary_root / ("gentle_examples_docs" + suffix)),
                               "tutorial-manifest-check"], calls)
                for call in run.call_args_list:
                    self.assertEqual(call.kwargs["cwd"], resource_root)
                    self.assertFalse(resource_root.is_relative_to(self.repo))
                    self.assertEqual(call.kwargs["timeout"], 120)

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
