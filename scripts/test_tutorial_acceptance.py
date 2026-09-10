"""Offline gate regressions using hand-crafted temporary receipts and fake binaries.

All sequences/projects/images here are synthetic placeholders constructed below;
they exercise provenance validation, never claim a live GUI or biological pass.
Run with python3 -m unittest scripts.test_tutorial_acceptance.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts import tutorial_acceptance as gate
from scripts import tutorial_acceptance_evidence as evidence
from scripts import publish_tutorial_gui_screenshots as publisher


REVISION = "a" * 40
REPO = Path(__file__).resolve().parents[1]


def write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")


def synthetic_run(root: Path) -> dict:
    step = {
        "id": "action", "prose_step": 1, "target": "control", "window": "window.main",
        "interaction": {"kind": "click"}, "timeout_class": "interactive",
        "scientific_effect": True, "verifiers": [{"kind": "state"}],
        "evidence": {"snapshot": "required", "screenshot": "required"},
    }
    contract = {"schema": "gentle.tutorial_gui_acceptance.v1", "profile": "smoke",
                "network": "offline", "steps": [step],
                "starter": {"example_id": "start"}, "oracle": {"example_id": "finish"}}
    candidate = {
        "schema": evidence.CANDIDATE_SCHEMA, "source_revision": REVISION,
        "cargo_lock_sha256": "c" * 64, "tutorial_manifest_sha256": "d" * 64,
        "binaries": {name: {"sha256": "b" * 64} for name in (
            "gentle", "gentle_cli", "gentle_examples_docs")},
        "chapters": [{"id": "chapter", "gui_acceptance": contract}],
        "workflow_bindings": {"chapter": {phase: {"sha256": "e" * 64}
                                           for phase in ("starter", "oracle")}},
    }
    environment = {"schema": "gentle.tutorial_acceptance_environment.v1",
                   "network_enforcement": "linux_network_namespace",
                   "network_namespace": {"current": "net:[2]", "parent": "net:[1]"},
                   "source_revision": REVISION, "cargo_lock_sha256": "c" * 64,
                   "git_status": "", "binaries": candidate["binaries"]}
    write(root / "environment.json", environment)
    folder = root / "chapter"
    for phase in ("starter", "oracle"):
        write(folder / f"{phase}.project.gentle.json", {"synthetic": phase})
    checkpoints = folder / "checkpoints"
    snapshot = {"generation": 2, "items": []}
    write(checkpoints / "action.snapshot.json", snapshot)
    (checkpoints / "action.raw.png").write_bytes(b"synthetic image bytes")
    for role in ("orientation", "context"):
        (checkpoints / f"action.{role}.svg").write_text('<svg><image href="action.raw.png"/></svg>')
    record = {
        "schema": publisher.EVIDENCE_SCHEMA, "source_revision": REVISION,
        "gentle_binary_sha256": "b" * 64, "tutorial_manifest_sha256": "d" * 64,
        "acceptance_contract_sha256": evidence.content_digest(contract),
        "chapter_id": "chapter", "step_id": "action", "prose_step": 1,
        "step_sha256": evidence.content_digest(step),
        "requested_target": {"semantic_id": "control", "window_id": "window.main"},
        "capture": {"captured_at_unix_ms": 1, "raw": {"sha256": evidence.digest(checkpoints / "action.raw.png")}},
        "semantic_snapshot": {"canonical_sha256": evidence.content_digest(snapshot),
                              "retained_file_sha256": evidence.digest(checkpoints / "action.snapshot.json")},
        "derived_views": [{"role": role, "path": str(checkpoints / f"action.{suffix}.svg"),
                           "sha256": evidence.digest(checkpoints / f"action.{suffix}.svg")}
                          for role, suffix in (("orientation", "orientation"), ("interaction_context", "context"))],
    }
    write(checkpoints / "action.screenshot.json", record)
    actual = {
        "id": "action", "status": "pass", "step_sha256": evidence.content_digest(step),
        "requested_action": {"kind": "click"}, "semantic_target": "control", "window": "window.main",
        "x11_event": {"synthetic": True}, "resolved_target": {"semantic_id": "control"},
        "before_generation": 1, "after_generation": 2,
        "before_fact": {"truth": "unsatisfied"}, "after_fact": {"truth": "satisfied"},
        "verifiers": [{"kind": "state", "status": "pass"}],
        "evidence": {"snapshot": {"sha256": evidence.digest(checkpoints / "action.snapshot.json")},
                     "screenshot": {"record_sha256": evidence.digest(checkpoints / "action.screenshot.json")}},
    }
    ledger = {
        "schema": "gentle.tutorial_gui_acceptance_ledger.v1", "chapter_id": "chapter",
        "status": "pass", "environment": environment, "manifest_sha256": "d" * 64,
        "acceptance_contract_sha256": evidence.content_digest(contract),
        "preflight": {"status": "pass", "starter_completion_truth": "unsatisfied", "oracle_completion_truth": "satisfied"},
        "completion": {"truth": "satisfied"}, "steps": [actual],
        "final_project_sha256": evidence.digest(folder / "starter.project.gentle.json"),
        "starter": {"example_id": "start", "workflow_source_sha256": "e" * 64},
        "oracle": {"example_id": "finish", "workflow_source_sha256": "e" * 64,
                   "project_sha256": evidence.digest(folder / "oracle.project.gentle.json")},
    }
    write(folder / "acceptance-ledger.json", ledger)
    write(root / "acceptance-report.json", {
        "schema": "gentle.tutorial_gui_acceptance_run.v1", "status": "pass", "chapter_count": 1,
        "chapters": [{"chapter_id": "chapter", "status": "pass",
                      "ledger_sha256": evidence.digest(folder / "acceptance-ledger.json")}],
    })
    return candidate


def change_ledger(root: Path, change) -> None:
    path = root / "chapter/acceptance-ledger.json"
    ledger = evidence.load(path)
    change(ledger)
    write(path, ledger)
    report = evidence.load(root / "acceptance-report.json")
    report["chapters"][0]["ledger_sha256"] = evidence.digest(path)
    write(root / "acceptance-report.json", report)


class TutorialAcceptanceTests(unittest.TestCase):
    def test_inventory_is_deterministic_and_does_not_invent_passes(self):
        first = evidence.coverage_inventory(REPO)
        self.assertEqual(first, evidence.coverage_inventory(REPO))
        self.assertTrue(all(r["execution_status"] == "not_run" for r in first["rows"]))
        self.assertTrue(any(r["coverage_kind"] == "manual_uncovered" for r in first["rows"]))
        self.assertTrue(any(r["coverage_kind"] == "gui_view_only" for r in first["rows"]))
        self.assertEqual(first["gui_contracts"], sum(r["gui_profile"] is not None for r in first["rows"]))

    def test_profiles_are_explicit_not_implicitly_cumulative(self):
        manifest = evidence.load(REPO / "docs/tutorial/manifest.json")
        smoke = gate.select_chapters(manifest, ["smoke"])
        offline = gate.select_chapters(manifest, ["offline-core"])
        self.assertEqual(gate.select_chapters(manifest, ["smoke", "offline-core", "smoke"]), smoke + offline)
        self.assertFalse(set(c["id"] for c in smoke) & set(c["id"] for c in offline))
        with self.assertRaises(gate.gui.AcceptanceFailure):
            gate.select_chapters(manifest, ["unknown"])
        online = {"chapters": [{"id": "online", "gui_acceptance": {"profile": "full", "network": "online"}}]}
        with self.assertRaises(ValueError):
            gate.select_chapters(online, ["full"])

    def test_clean_checkout_requires_exact_revision_and_no_untracked_inputs(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("Cargo.lock", "docs/tutorial/catalog.json", "docs/tutorial/manifest.json"):
                write(root / name, {})
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            subprocess.run(["git", "-C", str(root), "add", "."], check=True)
            subprocess.run(["git", "-C", str(root), "-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                            "-c", "commit.gpgsign=false", "commit", "-qm", "synthetic"], check=True)
            revision = evidence.git(root, "rev-parse", "HEAD")
            self.assertEqual(evidence.checkout_identity(root, revision)["source_revision"], revision)
            for bad in ("main", revision[:8], REVISION):
                with self.assertRaises(ValueError):
                    evidence.checkout_identity(root, bad)
            write(root / "untracked", {})
            with self.assertRaises(ValueError):
                evidence.checkout_identity(root, revision)

    def test_binary_probe_checks_reported_revision_not_just_path(self):
        with TemporaryDirectory() as tmp:
            binary = Path(tmp) / "binary"
            binary.write_bytes(b"fake executable identity")
            for revision, exit_code, valid in ((REVISION, 0, True), ("b" * 40, 0, False), (REVISION, 1, False)):
                result = subprocess.CompletedProcess([], exit_code,
                    f"GENtle test\nSource revision test+git.{revision}\n".encode(), b"")
                with patch.object(evidence.subprocess, "run", return_value=result):
                    if valid:
                        self.assertEqual(evidence.binary_identity(binary, REVISION)["sha256"], evidence.digest(binary))
                    else:
                        with self.assertRaises(ValueError):
                            evidence.binary_identity(binary, REVISION)

    def test_valid_receipts_derive_selection_without_reading_current_head(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = synthetic_run(root)
            selection = evidence.verify_gui_run(root, candidate)
            self.assertEqual(selection["chapters"], [{"chapter_id": "chapter", "checkpoints": ["action"]}])

    def test_empty_missing_duplicate_or_wrong_chapters_never_pass(self):
        for rows in ([], [{"chapter_id": "other"}], [{"chapter_id": "chapter"}] * 2):
            with self.subTest(rows=rows), TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidate = synthetic_run(root)
                report = evidence.load(root / "acceptance-report.json")
                report["chapters"] = rows
                write(root / "acceptance-report.json", report)
                with self.assertRaises(ValueError):
                    evidence.verify_gui_run(root, candidate)

    def test_stale_revision_manifest_contract_and_binary_each_fail(self):
        for field in ("source_revision", "tutorial_manifest_sha256", "contract", "binary"):
            with self.subTest(field=field), TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidate = synthetic_run(root)
                if field == "contract":
                    candidate["chapters"][0]["gui_acceptance"]["steps"][0]["prose_step"] = 2
                elif field == "binary":
                    candidate["binaries"]["gentle_cli"]["sha256"] = "f" * 64
                else:
                    candidate[field] = "f" * (40 if field == "source_revision" else 64)
                with self.assertRaises(ValueError):
                    evidence.verify_gui_run(root, candidate)

    def test_validly_rehashed_but_incomplete_ledgers_fail(self):
        changes = [
            lambda l: l.update(steps=[]),
            lambda l: l["steps"][0].update(verifiers=[]),
            lambda l: l["steps"][0].update(verifiers=[{"kind": "visible_claim", "status": "pass"}]),
            lambda l: l["steps"][0].update(semantic_target="wrong-control"),
            lambda l: l["steps"][0].update(x11_event=None),
            lambda l: l["steps"][0].update(after_generation=1),
            lambda l: l["steps"][0].update(after_fact={"truth": "unknown"}),
            lambda l: l["preflight"].update(starter_completion_truth="satisfied"),
            lambda l: l["oracle"].update(workflow_source_sha256="f" * 64),
            lambda l: l.update(completion={"truth": "unknown"}),
        ]
        for index, change in enumerate(changes):
            with self.subTest(index=index), TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidate = synthetic_run(root)
                change_ledger(root, change)
                with self.assertRaises(ValueError):
                    evidence.verify_gui_run(root, candidate)

    def test_altered_project_image_snapshot_or_sidecar_fails(self):
        for name in ("starter.project.gentle.json", "oracle.project.gentle.json", "acceptance-ledger.json",
                     "checkpoints/action.raw.png", "checkpoints/action.snapshot.json",
                     "checkpoints/action.screenshot.json", "checkpoints/action.context.svg"):
            with self.subTest(name=name), TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidate = synthetic_run(root)
                with (root / "chapter" / name).open("ab") as stream:
                    stream.write(b" ")
                with self.assertRaises(ValueError):
                    evidence.verify_gui_run(root, candidate)

    def test_missing_required_screenshot_fails(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = synthetic_run(root)
            (root / "chapter/checkpoints/action.raw.png").unlink()
            with self.assertRaises(FileNotFoundError):
                evidence.verify_gui_run(root, candidate)

    def test_capture_with_updated_hashes_but_wrong_revision_still_fails(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = synthetic_run(root)
            record_path = root / "chapter/checkpoints/action.screenshot.json"
            record = evidence.load(record_path)
            record["source_revision"] = "f" * 40
            write(record_path, record)
            change_ledger(root, lambda ledger: ledger["steps"][0]["evidence"]["screenshot"].update(
                record_sha256=evidence.digest(record_path)))
            with self.assertRaisesRegex(ValueError, "candidate/contract/step"):
                evidence.verify_gui_run(root, candidate)

    def test_missing_offline_enforcement_is_not_a_candidate_pass(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = synthetic_run(root)
            environment = evidence.load(root / "environment.json")
            environment["network_namespace"]["current"] = "net:[1]"
            write(root / "environment.json", environment)
            change_ledger(root, lambda ledger: ledger.update(environment=environment))
            with self.assertRaisesRegex(ValueError, "offline namespace"):
                evidence.verify_gui_run(root, candidate)

    def test_publication_is_explicit_and_requires_bound_green_run(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = synthetic_run(root / "run")
            write(root / "candidate.json", candidate)
            selection = evidence.verify_gui_run(root / "run", candidate)
            selection["capture_date"] = "2026-09-09"
            write(root / "selection.json", selection)
            publisher.publish(root / "run", root / "selection.json", root / "published",
                              root / "candidate.json", REVISION)
            publisher.check(root / "published")  # Archive check has no moving-HEAD assumption.
            with self.assertRaises(SystemExit):
                publisher.publish(root / "run", root / "selection.json", root / "published",
                                  root / "candidate.json", REVISION)
            with self.assertRaises(SystemExit):
                publisher.publish(root / "run", root / "selection.json", root / "other",
                                  root / "candidate.json", "b" * 40)
            (root / "published").rename(root / "relocated")
            publisher.check(root / "relocated")

    def test_output_receipts_keep_nonzero_exit_and_exact_bytes(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = gate.execute([sys.executable, "-c", "import os; os.write(1,b'raw\\xff'); os.write(2,b'err'); raise SystemExit(7)"],
                                  root, dict(os.environ), root / "logs", timeout=10)
            self.assertEqual(result["status"], "fail")
            self.assertEqual(result["exit_code"], 7)
            self.assertEqual((root / "logs/stdout").read_bytes(), b"raw\xff")
            self.assertEqual(result["stdout"]["sha256"], evidence.digest(root / "logs/stdout"))

    def test_timeout_has_receipt_and_does_not_leave_active_process(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = gate.execute([sys.executable, "-c", "import time; time.sleep(60)"],
                                  root, dict(os.environ), root / "logs", timeout=0.05)
            self.assertEqual(result["status"], "timeout")
            self.assertLess(result["exit_code"], 0)
            self.assertTrue((root / "logs/receipt.json").is_file())

    def test_unspawnable_command_retains_a_missing_dependency_receipt(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = gate.execute([str(root / "absent")], root, dict(os.environ), root / "logs", timeout=10)
            self.assertEqual(result["status"], "missing_dependency")
            self.assertIn("message", result)
            self.assertTrue((root / "logs/receipt.json").is_file())

    def test_environment_removes_api_keys_and_inherited_gentle_configuration(self):
        with TemporaryDirectory() as tmp, patch.dict(os.environ, {
            "GENTLE_TEST_ONLINE": "1", "MISTRAL_API_KEY": "synthetic-secret", "GENTLE_RESOURCE_ROOT": "/private"
        }):
            env = gate.isolated_environment(Path(tmp))
            self.assertNotIn("MISTRAL_API_KEY", env)
            self.assertNotIn("GENTLE_RESOURCE_ROOT", env)
            self.assertNotIn("GENTLE_TEST_ONLINE", env)
            self.assertEqual(env["GENTLE_SKIP_REMOTE_TESTS"], "1")
            self.assertTrue(env["HOME"].startswith(tmp))

    def test_preflight_failure_is_retained_and_cannot_overwrite_previous_run(self):
        with TemporaryDirectory() as tmp:
            output = Path(tmp) / "run"
            args = argparse.Namespace(repo_root=REPO, evidence_dir=output, candidate=REVISION,
                profile=["smoke"], supersedes=None, gentle=Path("missing"), gentle_cli=Path("missing"),
                examples_docs=Path("missing"), parent_network_namespace=None)
            with patch.object(gate, "bind_candidate", side_effect=FileNotFoundError("synthetic missing binary")):
                result = gate.run(args)
            self.assertEqual(result["status"], "incomplete")
            self.assertEqual(result["failure_class"], "missing_dependency")
            self.assertTrue(all(c["status"] == "not_run" for c in result["checks"]))
            before = (output / "tutorial-acceptance-report.json").read_bytes()
            with self.assertRaises(ValueError):
                gate.run(args)
            self.assertEqual(before, (output / "tutorial-acceptance-report.json").read_bytes())

    def test_missing_x11_tool_leaves_a_structured_runner_report(self):
        with TemporaryDirectory() as tmp, patch.dict(os.environ, {"DISPLAY": ":synthetic"}), \
                patch.object(gate.gui.platform, "system", return_value="Linux"), \
                contextlib.redirect_stderr(io.StringIO()):
            output = Path(tmp) / "run"
            code = gate.gui.main(["--repo-root", str(REPO), "--evidence-dir", str(output),
                                  "--xdotool", str(Path(tmp) / "absent")])
            self.assertEqual(code, 1)
            self.assertEqual(evidence.load(output / "acceptance-report.json")["failure_class"], "missing_dependency")

    def test_coordinator_composes_checks_selection_and_supersedes_without_old_passes(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = synthetic_run(root / "synthetic-source")
            for name, binary in candidate["binaries"].items():
                binary["path"] = str(root / name)
            prior = root / "prior-report.json"
            write(prior, {"schema": gate.RUN_SCHEMA, "status": "fail", "candidate_revision": "b" * 40})
            output = root / "new-run"
            args = argparse.Namespace(repo_root=REPO, evidence_dir=output, candidate=REVISION,
                profile=["smoke"], supersedes=prior, gentle=Path("gentle"), gentle_cli=Path("gentle_cli"),
                examples_docs=Path("helper"), parent_network_namespace="net:[1]")
            calls = []
            readlink = os.readlink

            def execute(argv, _root, _env, _folder, _timeout=None, **kwargs):
                calls.append(argv)
                if "--evidence-dir" in argv:
                    synthetic_run(Path(argv[argv.index("--evidence-dir") + 1]))
                return {"status": "pass", "exit_code": 0, "argv": argv}

            with patch.object(gate, "bind_candidate", return_value=candidate), \
                    patch.object(gate, "select_chapters", return_value=candidate["chapters"]), \
                    patch.object(gate, "recheck_candidate"), \
                    patch.object(gate, "execute", side_effect=execute), \
                    patch.object(gate.platform, "system", return_value="Linux"), \
                    patch.object(gate.os, "readlink", side_effect=lambda path, **kw: (
                        "net:[2]" if str(path) == "/proc/self/ns/net" else readlink(path, **kw)
                    )), \
                    patch.dict(os.environ, {"DISPLAY": ":synthetic"}):
                result = gate.run(args)
            self.assertEqual(result["status"], "pass", result)
            self.assertEqual(result["supersedes"]["sha256"], evidence.digest(prior))
            self.assertEqual(len(calls), 4)
            self.assertEqual([c[1] for c in calls[:3]], list(gate.CHECKS))
            self.assertEqual(result["publication"], "not_requested")
            self.assertEqual(result["scientific_study_acceptance"], "not_run")
            self.assertTrue(all(c["status"] == "pass" for c in result["checks"]))
            self.assertEqual(evidence.load(output / "selection.json")["chapters"],
                             [{"chapter_id": "chapter", "checkpoints": ["action"]}])

    def test_end_check_rejects_replaced_file_at_the_same_path(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "binary"
            path.write_bytes(b"first")
            candidate = {"source_revision": REVISION, "inputs": [],
                         "binaries": {"gentle": evidence.file_record(path)}}
            with patch.object(evidence, "checkout_identity", return_value={}):
                gate.recheck_candidate(Path(tmp), candidate)
                path.write_bytes(b"replacement")
                with self.assertRaises(ValueError):
                    gate.recheck_candidate(Path(tmp), candidate)


if __name__ == "__main__":
    unittest.main()
