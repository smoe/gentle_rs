#!/usr/bin/env python3
"""Offline release-policy regressions; all Git/installer fixtures are synthetic.

Fixtures are recreated in temporary directories by these tests. No published
tag, private input, network service, Cargo build or real installer is used.
Run with: python3 -m unittest scripts.test_release_candidate -v
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import textwrap
import tomllib
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts import release_candidate as policy
from scripts import check_tutorial_checkouts as checkouts


class NativeBuildSettingsTests(unittest.TestCase):
    def test_only_internal_version_labels_select_unoptimized_installers(self) -> None:
        for tag in ("v0.1.0-internal.11", "v0.1.0-internal.12", "v2.3.4-internal.1+build.7"):
            with self.subTest(tag=tag):
                self.assertEqual(policy.native_build_settings(tag), {
                    "native_profile": "dev", "native_target_subdir": "debug",
                })
        for tag in ("v0.1.0", "v2.3.4+internal.11", "v0.1.0-rc.1", "v0.1.0-notinternal.11"):
            with self.subTest(tag=tag):
                self.assertEqual(policy.native_build_settings(tag), {
                    "native_profile": "release", "native_target_subdir": "release",
                })
        for tag in ("main", "v0.1.0-internal.11\n", "../v0.1.0-internal.11"):
            with self.subTest(tag=tag), self.assertRaises(ValueError):
                policy.native_build_settings(tag)


class ReleaseCandidateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.git("init", "-q")
        self.git("config", "user.email", "synthetic@example.invalid")
        self.git("config", "user.name", "Synthetic release test")
        (self.root / "Cargo.toml").write_text('[workspace.package]\nversion = "0.1.0-internal.10"\n')
        (self.root / "Cargo.lock").write_text("# Synthetic lockfile\nversion = 4\n")
        self.git("add", ".")
        self.git("commit", "-qm", "synthetic candidate")
        self.sha = self.git("rev-parse", "HEAD")
        self.tag = "v0.1.0-internal.10"
        self.git("remote", "add", "origin", str(self.root))
        self.env = {
            "CANDIDATE_EVENT": "workflow_dispatch", "CANDIDATE_ACTION": "",
            "CANDIDATE_SHA": self.sha, "PUBLISH_REQUESTED": "false",
            "RELEASE_TAG": self.tag, "WORKFLOW_REVISION": self.sha,
            "CANDIDATE_EVENT_SHA": self.sha,
        }

    def git(self, *args: str) -> str:
        return subprocess.check_output(["git", "-C", str(self.root), *args], text=True).strip()

    def test_validation_needs_no_tag_and_does_not_fetch_or_create_one(self) -> None:
        with patch.object(policy, "tag_revision", side_effect=AssertionError("must not fetch a tag")):
            record = policy.prepare(self.root, self.env)
        self.assertFalse(record["publish"])
        self.assertEqual(record["mode"], "validate_only")
        self.assertEqual(record["revision"], self.sha)
        self.assertEqual(self.git("tag", "--list"), "")

    def test_manual_runs_reject_empty_short_branch_and_injected_sha(self) -> None:
        for value in ("", "main", self.sha[:8], "refs/tags/" + self.tag, "$(touch bad)"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                policy.prepare(self.root, {**self.env, "CANDIDATE_SHA": value})

    def test_checkout_must_match_explicit_sha(self) -> None:
        with self.assertRaisesRegex(ValueError, "expected candidate SHA"):
            policy.prepare(self.root, {**self.env, "CANDIDATE_SHA": "0" * 40})

    def test_package_label_and_tag_input_validation(self) -> None:
        for tag in ("v0.1.0-internal.11", "main", "../v0.1", "v0.1\noutput=true", "v0.1;bad"):
            with self.subTest(tag=tag), self.assertRaises(ValueError):
                policy.prepare(self.root, {**self.env, "RELEASE_TAG": tag})

    def test_validation_ignores_older_tag_but_publication_rejects_it(self) -> None:
        self.git("tag", self.tag)
        (self.root / "change").write_text("next synthetic revision\n")
        self.git("add", ".")
        self.git("commit", "-qm", "synthetic later candidate")
        later = self.git("rev-parse", "HEAD")
        env = {**self.env, "CANDIDATE_SHA": later}
        self.assertEqual(policy.prepare(self.root, env)["revision"], later)
        with self.assertRaisesRegex(ValueError, "approved candidate SHA"):
            policy.prepare(self.root, {**env, "PUBLISH_REQUESTED": "true"})
        self.assertEqual(self.git("rev-parse", self.tag), self.sha)

    def test_explicit_publication_accepts_matching_annotated_tag(self) -> None:
        self.git("tag", "-am", "synthetic tag", self.tag)
        record = policy.prepare(self.root, {**self.env, "PUBLISH_REQUESTED": "true"})
        self.assertTrue(record["publish"])
        self.assertEqual(record["mode"], "publish")

    def test_internal_profile_does_not_depend_on_publication_mode(self) -> None:
        self.git("tag", self.tag)
        for event, action, publish in (("workflow_dispatch", "", "false"),
                                       ("workflow_dispatch", "", "true"),
                                       ("release", "published", "false")):
            with self.subTest(event=event, publish=publish):
                record = policy.prepare(self.root, {
                    **self.env, "CANDIDATE_EVENT": event,
                    "CANDIDATE_ACTION": action, "PUBLISH_REQUESTED": publish,
                })
                self.assertEqual(record["native_profile"], "dev")
                self.assertEqual(record["native_target_subdir"], "debug")

    def test_publication_requires_an_existing_tag(self) -> None:
        with self.assertRaises(subprocess.CalledProcessError):
            policy.prepare(self.root, {**self.env, "PUBLISH_REQUESTED": "true"})

    def test_only_explicit_publication_events_are_allowed(self) -> None:
        self.assertFalse(policy.publication_allowed("workflow_dispatch", "", "false"))
        self.assertTrue(policy.publication_allowed("workflow_dispatch", "", "true"))
        self.assertFalse(policy.publication_allowed("push", "", "false"))
        self.assertTrue(policy.publication_allowed("release", "published", "false"))
        for event, action, requested in (
            ("push", "", "true"), ("release", "created", "false"),
            ("pull_request", "", "false"), ("workflow_dispatch", "", "yes"),
        ):
            with self.subTest(event=event, action=action), self.assertRaises(ValueError):
                policy.publication_allowed(event, action, requested)

    def test_tag_push_is_validation_only_even_on_a_tag_ref(self) -> None:
        record = policy.prepare(self.root, {
            **self.env, "CANDIDATE_EVENT": "push", "CANDIDATE_REF": f"refs/tags/{self.tag}",
        })
        self.assertFalse(record["publish"])
        with self.assertRaises(ValueError):
            policy.prepare(self.root, {**self.env, "CANDIDATE_EVENT": "push", "CANDIDATE_REF": "refs/heads/main"})

    def test_release_event_is_bound_to_its_original_sha(self) -> None:
        with self.assertRaisesRegex(ValueError, "expected candidate SHA"):
            policy.prepare(self.root, {
                **self.env, "CANDIDATE_EVENT": "release", "CANDIDATE_ACTION": "published",
                "CANDIDATE_EVENT_SHA": "0" * 40,
            })

    def test_dirty_metadata_and_mismatched_lock_fail(self) -> None:
        record = policy.prepare(self.root, self.env)
        with self.assertRaisesRegex(ValueError, "Cargo.lock"):
            policy.validate_checkout(self.root, self.tag, self.sha, "0" * 64)
        (self.root / "Cargo.lock").write_text("# changed\nversion = 4\n")
        with self.assertRaisesRegex(ValueError, "uncommitted"):
            policy.prepare(self.root, self.env)
        with self.assertRaisesRegex(ValueError, "Cargo.lock"):
            policy.validate_checkout(self.root, self.tag, self.sha, record["cargo_lock_sha256"])

    def checkout_policy_fixture(self, attributes: bytes) -> dict:
        self.git("config", "core.autocrlf", "false")
        (self.root / "Cargo.lock").write_bytes(b"# Synthetic lockfile\nversion = 4\n")
        (self.root / ".gitattributes").write_bytes(attributes)
        self.git("add", "Cargo.lock", ".gitattributes")
        self.git("commit", "-qm", "synthetic checkout policy")
        self.sha = self.git("rev-parse", "HEAD")
        return policy.prepare(self.root, {**self.env, "CANDIDATE_SHA": self.sha})

    def test_repo_checkout_policy_preserves_candidate_lockfile_on_both_modes(self) -> None:
        attributes = (Path(__file__).resolve().parents[1] / ".gitattributes").read_bytes()
        candidate = self.checkout_policy_fixture(attributes)
        for mode in checkouts.MODES:
            with self.subTest(mode=mode[0]):
                checkout = self.root / f"checkout-{mode[0]}"
                checkouts.prepare_checkout(self.root, checkout, mode, revision=self.sha)
                record = policy.validate_checkout(
                    checkout, self.tag, self.sha, candidate["cargo_lock_sha256"]
                )
                self.assertEqual(record["cargo_lock_sha256"], candidate["cargo_lock_sha256"])
                with (checkout / "Cargo.lock").open("ab") as stream:
                    stream.write(b"# changed after checkout\n")
                with self.assertRaisesRegex(ValueError, "Cargo.lock"):
                    policy.validate_checkout(checkout, self.tag, self.sha, candidate["cargo_lock_sha256"])

    def test_unprotected_crlf_lockfile_is_not_silently_normalized(self) -> None:
        candidate = self.checkout_policy_fixture(b"# Deliberately missing LF rule\n")
        checkout = self.root / "checkout-crlf"
        checkouts.prepare_checkout(self.root, checkout, checkouts.MODES[1], revision=self.sha)
        self.assertIn(b"\r\n", (checkout / "Cargo.lock").read_bytes())
        with self.assertRaisesRegex(ValueError, "Cargo.lock no longer matches"):
            policy.validate_checkout(checkout, self.tag, self.sha, candidate["cargo_lock_sha256"])

    def test_workflow_cli_writes_receipt_and_typed_job_outputs(self) -> None:
        script = Path(policy.__file__).resolve()
        output = self.root / "receipt.json"
        job_outputs = self.root / "job-outputs.txt"
        env = {**os.environ, **self.env, "GITHUB_OUTPUT": str(job_outputs)}
        subprocess.run(
            [sys.executable, str(script), "prepare", "--root", str(self.root), "--output", str(output)],
            env=env, check=True, capture_output=True, text=True,
        )
        record = json.loads(output.read_text())
        values = dict(line.split("=", 1) for line in job_outputs.read_text().splitlines())
        self.assertEqual(values["publish"], "false")
        self.assertEqual(values["revision"], self.sha)
        self.assertEqual(values["mode"], "validate_only")
        self.assertEqual(values["native_profile"], "dev")
        self.assertEqual(values["native_target_subdir"], "debug")
        self.assertEqual(record["native_profile"], values["native_profile"])
        verify_env = {
            **env, "EXPECTED_REVISION": self.sha, "EXPECTED_LOCK_SHA256": record["cargo_lock_sha256"],
            "CANDIDATE_MODE": "validate_only",
        }
        subprocess.run(
            [sys.executable, str(script), "verify", "--root", str(self.root)],
            env=verify_env, check=True, capture_output=True, text=True,
        )
        failed = subprocess.run(
            [sys.executable, str(script), "verify", "--root", str(self.root)],
            env={**verify_env, "EXPECTED_LOCK_SHA256": ""}, capture_output=True, text=True,
        )
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("lockfile SHA-256 is required", failed.stderr)

    def test_failed_prepare_does_not_write_receipt_or_job_outputs(self) -> None:
        output = self.root / "receipt.json"
        job_outputs = self.root / "job-outputs.txt"
        result = subprocess.run(
            [sys.executable, str(Path(policy.__file__).resolve()), "prepare", "--root", str(self.root), "--output", str(output)],
            env={**os.environ, **self.env, "CANDIDATE_SHA": "main", "GITHUB_OUTPUT": str(job_outputs)},
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(output.exists())
        self.assertFalse(job_outputs.exists())

    def installers(self) -> tuple[Path, dict, list[Path]]:
        candidate = policy.prepare(self.root, self.env)
        folder = self.root / "artifacts"
        folder.mkdir()
        receipts = []
        for platform, extension in (("linux", "tar.gz"), ("macos", "dmg"), ("windows", "zip")):
            name = f"gentle-{self.tag}-{platform}-x64"
            suffix = "-dev" if candidate["native_profile"] == "dev" else ""
            artifact = folder / f"{name}{suffix}.{extension}"
            artifact.write_bytes(b"synthetic installer, not executable")
            receipt = {
                **{k: candidate[k] for k in ("tag", "revision", "cargo_lock_sha256", "workflow_revision", "mode")},
                "schema": "gentle.release_build.v1", "platform": platform, "arch": "x64",
                "features": [], "default_features": True, "profile": candidate["native_profile"],
                "incremental": False, "debug": 0,
                "binaries": ["gentle", "gentle_cli", "gentle_mcp",
                             "gentle_examples_docs", "gentle_publication_report"],
                "rustc": "synthetic rustc", "cargo": "synthetic cargo",
                "artifact_sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
            }
            path = folder / f"{name}.build.json"
            path.write_text(json.dumps(receipt))
            receipts.append(path)
        return folder, candidate, receipts

    def test_three_platform_collection_preserves_validation_only_state(self) -> None:
        folder, candidate, _ = self.installers()
        result = policy.collect_installers(folder, candidate)
        self.assertEqual(result["revision"], self.sha)
        self.assertEqual(result["mode"], "validate_only")
        self.assertEqual(result["profile"], "dev")
        self.assertEqual(len(result["artifacts"]), 3)
        self.assertTrue(all(row["bytes"] > 0 for row in result["artifacts"]))

    def test_final_release_collection_requires_optimized_receipts(self) -> None:
        (self.root / "Cargo.toml").write_text('[workspace.package]\nversion = "0.1.0"\n')
        self.git("add", "Cargo.toml")
        self.git("commit", "-qm", "synthetic final release")
        self.sha = self.git("rev-parse", "HEAD")
        self.tag = "v0.1.0"
        self.env.update(CANDIDATE_SHA=self.sha, CANDIDATE_EVENT_SHA=self.sha,
                        RELEASE_TAG=self.tag, WORKFLOW_REVISION=self.sha)
        folder, candidate, paths = self.installers()
        self.assertEqual(candidate["native_profile"], "release")
        result = policy.collect_installers(folder, candidate)
        self.assertEqual(result["profile"], "release")
        self.assertTrue(all("-dev." not in item["name"] for item in result["artifacts"]))
        for path in paths:
            receipt = json.loads(path.read_text())
            receipt["profile"] = "dev"
            path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "profile"):
            policy.collect_installers(folder, candidate)

    def test_internal_collection_rejects_mixed_and_consistently_wrong_profiles(self) -> None:
        folder, candidate, paths = self.installers()
        for path in paths:
            receipt = json.loads(path.read_text())
            receipt["profile"] = "release"
            path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, "profile"):
                policy.collect_installers(folder, candidate)
        for key, value in (("native_profile", "release"), ("native_target_subdir", "release")):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, key):
                policy.collect_installers(folder, {**candidate, key: value})

    def test_internal_archive_name_exposes_unoptimized_profile(self) -> None:
        folder, candidate, _ = self.installers()
        path = next(folder.glob("*.zip"))
        self.assertTrue(path.name.endswith("-dev.zip"))
        path.rename(path.with_name(path.name.replace("-dev.zip", ".zip")))
        with self.assertRaisesRegex(ValueError, "named"):
            policy.collect_installers(folder, candidate)

    def test_standalone_collector_uses_the_desktop_binary_contract(self) -> None:
        folder, candidate, _ = self.installers()
        output = self.root / "release-attributes.json"
        subprocess.run(
            [sys.executable, str(Path(policy.__file__).resolve()), "collect",
             "--root", str(self.root), "--artifacts", str(folder), "--output", str(output)],
            env={**os.environ, **self.env, "EXPECTED_REVISION": self.sha,
                 "EXPECTED_LOCK_SHA256": candidate["cargo_lock_sha256"],
                 "CANDIDATE_MODE": "validate_only"},
            check=True, capture_output=True, text=True,
        )
        self.assertEqual(json.loads(output.read_text())["revision"], self.sha)

    def test_consistently_wrong_revision_still_fails_against_candidate(self) -> None:
        folder, candidate, paths = self.installers()
        for path in paths:
            receipt = json.loads(path.read_text())
            receipt["revision"] = "0" * 40
            path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "revision"):
            policy.collect_installers(folder, candidate)

    def test_incompatible_or_missing_receipts_fail_closed(self) -> None:
        folder, candidate, paths = self.installers()
        original = json.loads(paths[0].read_text())
        for key, value in (
            ("schema", "unknown"), ("tag", "v0.1.0-internal.11"), ("revision", "0" * 40),
            ("cargo_lock_sha256", "0" * 64), ("workflow_revision", "0" * 40),
            ("mode", "publish"), ("profile", "release-fast"), ("features", ["script-interfaces"]),
            ("incremental", True), ("incremental", None), ("debug", 1), ("debug", None),
            ("default_features", False), ("default_features", None),
            ("binaries", None), ("binaries", original["binaries"][:-1]),
            ("binaries", [*original["binaries"], "gentle_js", "gentle_lua"]),
            ("rustc", ""), ("cargo", ""),
        ):
            changed = copy.deepcopy(original)
            changed[key] = value
            paths[0].write_text(json.dumps(changed))
            with self.subTest(field=key), self.assertRaises(ValueError):
                policy.collect_installers(folder, candidate)
        for key in ("features", "default_features", "binaries", "incremental", "debug"):
            changed = copy.deepcopy(original)
            del changed[key]
            paths[0].write_text(json.dumps(changed))
            with self.subTest(missing=key), self.assertRaises(ValueError):
                policy.collect_installers(folder, candidate)
        paths[0].unlink()
        with self.assertRaisesRegex(ValueError, "exactly Linux"):
            policy.collect_installers(folder, candidate)

    def test_modified_or_empty_installer_fails(self) -> None:
        folder, candidate, _ = self.installers()
        path = next(folder.glob("*.zip"))
        path.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "digest"):
            policy.collect_installers(folder, candidate)
        path.write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "non-empty"):
            policy.collect_installers(folder, candidate)

    def test_duplicate_platform_or_misnamed_archive_fails(self) -> None:
        folder, candidate, paths = self.installers()
        original = json.loads(paths[0].read_text())
        paths[0].write_text(json.dumps({**original, "platform": "macos"}))
        with self.assertRaisesRegex(ValueError, "exactly Linux"):
            policy.collect_installers(folder, candidate)
        paths[0].write_text(json.dumps(original))
        next(folder.glob("*.zip")).rename(folder / "unrelated.zip")
        with self.assertRaisesRegex(ValueError, "non-empty candidate artifact"):
            policy.collect_installers(folder, candidate)


class WorkflowWiringTests(unittest.TestCase):
    def test_replayed_tutorial_reports_use_current_package_version(self) -> None:
        root = Path(__file__).resolve().parents[1]
        version = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["package"]["version"]
        checked = 0

        def check_versions(value, path):
            nonlocal checked
            if isinstance(value, dict):
                for key, item in value.items():
                    if key in ("gentle_version", "selection_audit_generator_revision"):
                        checked += 1
                        self.assertEqual(item, version,
                                         f"{path}: regenerate tutorial artifacts after a version bump; "
                                         "do not relabel provenance without replaying the workflow")
                    check_versions(item, path)
            elif isinstance(value, list):
                for item in value:
                    check_versions(item, path)

        for path in sorted((root / "docs/tutorial/generated/artifacts").rglob("*.report.json")):
            with self.subTest(path=path.relative_to(root)):
                check_versions(json.loads(path.read_bytes()), path.relative_to(root))
        self.assertGreater(checked, 0, "expected version-bound retained tutorial reports")

    def test_windows_is_unconditional_and_unix_jobs_follow_selection(self) -> None:
        root = Path(__file__).resolve().parents[1]
        text = (root / ".github/workflows/ci.yml").read_text()
        windows = text.split("\n  windows:\n", 1)[1].split("    steps:\n", 1)[0]
        self.assertIn("runs-on: windows-latest", windows)
        self.assertNotIn("needs:", windows)
        self.assertNotIn("if:", windows)
        for platform in ("macos", "linux"):
            job = text.split(f"\n  {platform}:\n", 1)[1].split("    steps:\n", 1)[0]
            self.assertIn("needs: select-platform", job)
            self.assertIn(f"if: needs.select-platform.outputs.platform == '{platform}'", job)

    def test_windows_checks_adapter_processes_before_long_test_suites(self) -> None:
        root = Path(__file__).resolve().parents[1]
        text = (root / ".github/workflows/ci.yml").read_text()
        windows = text.split("\n  windows:\n", 1)[1].split("\n  ci-summary:\n", 1)[0]
        step = "      - name: CLI and MCP process-boundary regressions\n"
        self.assertIn(step, windows)
        body = windows.split(step, 1)[1].split("\n      - name:", 1)[0]
        self.assertIn(
            "run: cargo test -q --locked --test adapter_error_contract "
            "--test reporter_construct_handoff_cli", body)
        self.assertNotIn("continue-on-error:", body)
        self.assertNotIn("if:", body)
        for later in ("Examples and tutorial schema/drift checks",
                      "Workflow example runtime tests", "Full test suite"):
            self.assertLess(windows.index(step), windows.index(f"- name: {later}\n"))

    def test_native_jobs_check_tutorial_discovery_before_long_runtime_suites(self) -> None:
        root = Path(__file__).resolve().parents[1]
        text = (root / ".github/workflows/ci.yml").read_text()
        step = "      - name: Agent tutorial discovery regressions\n"
        for platform in ("macos", "linux", "windows"):
            with self.subTest(platform=platform):
                following = text.split(f"\n  {platform}:\n", 1)[1]
                job = re.split(r"\n  [a-z][a-z0-9-]*:\n", following, maxsplit=1)[0]
                self.assertIn(step, job)
                body = job.split(step, 1)[1].split("\n      - name:", 1)[0]
                self.assertIn(
                    "run: cargo test -q --locked --lib -j1 app::tests::agent_gui_context_ -- --test-threads=1", body)
                self.assertNotIn("continue-on-error:", body)
                self.assertNotIn("if:", body)
                for later in ("Examples and tutorial schema/drift checks",
                              "Workflow example runtime tests", "Full test suite"):
                    self.assertLess(job.index(step), job.index(f"- name: {later}\n"))

    def test_sampled_platform_is_unix_and_manual_selection_is_preserved(self) -> None:
        root = Path(__file__).resolve().parents[1]
        text = (root / ".github/workflows/ci.yml").read_text()
        selection = text.split("\n  select-platform:\n", 1)[1].split("\n  headless-linux:\n", 1)[0]
        script = textwrap.dedent(selection.split("        run: |\n", 1)[1])
        cases = [(requested, sha, "macos" if sha % 2 == 0 else "linux")
                 for requested in ("", "sampled") for sha in (0, 1, 2, 3, 0xffffffff)]
        cases.extend((platform, 0, platform) for platform in ("macos", "linux", "windows"))
        with TemporaryDirectory() as directory:
            output = Path(directory) / "output"
            for requested, sha, expected in cases:
                with self.subTest(requested=requested, sha=sha):
                    output.write_text("")
                    completed = subprocess.run(
                        ["bash", "-c", script],
                        env={**os.environ, "REQUESTED_PLATFORM": requested,
                             "GITHUB_SHA": f"{sha:08x}" + "0" * 32, "GITHUB_OUTPUT": str(output)},
                        capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                    self.assertEqual(output.read_text().strip(), f"platform={expected}")

    def test_workflows_use_reviewed_node24_compatible_actions(self) -> None:
        # Upstream action.yml runtimes were checked on 2026-09-16. The two
        # composite actions use shell steps or Node 24 sub-actions already.
        reviewed_actions = {
            "actions/checkout@v5",
            "actions/setup-python@v6",
            "actions/cache@v5",
            "actions/upload-artifact@v6",
            "actions/download-artifact@v7",
            "actions/attest-build-provenance@v3",
            "docker/setup-buildx-action@v4",
            "docker/build-push-action@v7",
            "docker/login-action@v4",
            "docker/metadata-action@v6",
            "softprops/action-gh-release@v3",
            "dtolnay/rust-toolchain@stable",
        }
        root = Path(__file__).resolve().parents[1]
        for workflow in sorted((root / ".github/workflows").iterdir()):
            if workflow.suffix not in (".yml", ".yaml"):
                continue
            text = workflow.read_text()
            with self.subTest(workflow=workflow.name):
                self.assertNotRegex(
                    text,
                    r"(?m)^\s*(?:ACTIONS_ALLOW_USE_UNSECURE_NODE_VERSION|"
                    r"FORCE_JAVASCRIPT_ACTIONS_TO_NODE24):",
                )
                for reference in re.findall(r"(?m)^\s*(?:-\s+)?uses:\s+(\S+)", text):
                    reference = reference.strip("\"'")
                    if not reference.startswith("./"):
                        self.assertIn(reference, reviewed_actions,
                                      "Review upstream runs.using before changing action versions")

    def test_publication_jobs_are_explicit_and_default_permissions_are_read_only(self) -> None:
        root = Path(__file__).resolve().parents[1]
        for name, job in (("release.yml", "publish-release"), ("container.yml", "container-publish")):
            with self.subTest(workflow=name):
                text = (root / ".github/workflows" / name).read_text()
                before_jobs = text.split("\njobs:")[0]
                self.assertIn("permissions:\n  contents: read", before_jobs)
                # Input details are checked separately to catch a true default.
                self.assertNotIn("default: true", before_jobs)
                self.assertIn("        default: false", before_jobs)
                self.assertIn("      candidate_sha:", before_jobs)
                publication = text.split(f"\n  {job}:\n")[1]
                self.assertIn("needs.candidate.outputs.publish == 'true'", publication)
                self.assertIn("inputs.publish == true", publication)
                self.assertIn("github.event.action == 'published'", publication)
                self.assertIn("ref: ${{ needs.candidate.outputs.revision }}", publication)
                self.assertIn("scripts/release_candidate.py", publication)
                self.assertIn("uses: ./.github/workflows/release-candidate.yml", text)

    def test_candidate_checks_run_on_every_push_and_pr(self) -> None:
        root = Path(__file__).resolve().parents[1]
        text = (root / ".github/workflows/ci.yml").read_text()
        self.assertIn("python3 -m unittest scripts.test_release_candidate -v", text)

    def test_ci_summary_requires_release_policy_and_build_success(self) -> None:
        root = Path(__file__).resolve().parents[1]
        text = (root / ".github/workflows/ci.yml").read_text()
        summary = text.split("\n  ci-summary:\n", 1)[1]
        dependencies = summary.split("    needs:\n", 1)[1].split("    if:", 1)[0]
        self.assertIn("      - release-policy\n", dependencies)
        self.assertIn("      - windows\n", dependencies)
        self.assertIn("    if: always()", summary)
        self.assertIn("RELEASE_POLICY_RESULT: ${{ needs.release-policy.result }}", summary)
        script = textwrap.dedent(summary.split("        run: |\n", 1)[1])

        for platform in ("macos", "linux", "windows"):
            selected = f"{platform.upper()}_RESULT"
            results = {
                "PLATFORM": platform,
                "RELEASE_POLICY_RESULT": "success",
                "HEADLESS_RESULT": "success",
                "MACOS_RESULT": "skipped",
                "LINUX_RESULT": "skipped",
                "WINDOWS_RESULT": "success",
                selected: "success",
            }
            cases = [(results, True)]
            for required in ("RELEASE_POLICY_RESULT", "HEADLESS_RESULT", "WINDOWS_RESULT", selected):
                for outcome in ("failure", "cancelled", "skipped", ""):
                    cases.append(({**results, required: outcome}, False))
            for case, should_pass in cases:
                with self.subTest(results=case):
                    completed = subprocess.run(
                        ["bash", "-c", script], env={**os.environ, **case},
                        capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(
                        completed.returncode == 0, should_pass,
                        completed.stdout + completed.stderr,
                    )


if __name__ == "__main__":
    unittest.main()
