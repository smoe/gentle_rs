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
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts import release_candidate as policy


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
            artifact = folder / f"{name}.{extension}"
            artifact.write_bytes(b"synthetic installer, not executable")
            receipt = {
                **{k: candidate[k] for k in ("tag", "revision", "cargo_lock_sha256", "workflow_revision", "mode")},
                "schema": "gentle.release_build.v1", "platform": platform, "arch": "x64",
                "features": ["script-interfaces"], "profile": "release",
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
        self.assertEqual(len(result["artifacts"]), 3)
        self.assertTrue(all(row["bytes"] > 0 for row in result["artifacts"]))

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
            ("mode", "publish"), ("profile", "release-fast"), ("features", []),
            ("rustc", ""), ("cargo", ""),
        ):
            changed = copy.deepcopy(original)
            changed[key] = value
            paths[0].write_text(json.dumps(changed))
            with self.subTest(field=key), self.assertRaises(ValueError):
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


if __name__ == "__main__":
    unittest.main()
