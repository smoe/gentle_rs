"""Offline identity and evidence checks for the external tutorial acceptance gate.

These checks compare retained evidence; they neither drive the GUI nor infer
biological results. Historical screenshot integrity is a separate concern.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess


CANDIDATE_SCHEMA = "gentle.tutorial_acceptance_candidate.v1"
COVERAGE_SCHEMA = "gentle.tutorial_acceptance_coverage.v1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"Expected a JSON object: {path}")
    return value


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def content_digest(value: object) -> str:
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")).hexdigest()


def file_record(path: Path) -> dict:
    return {"path": str(path), "sha256": digest(path)}


def child_path(root: Path, name: str) -> Path:
    path = (root / name).resolve()
    require(path.is_relative_to(root.resolve()), f"Evidence path escapes its root: {name}")
    return path


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True, timeout=15
    ).strip()


def checkout_identity(root: Path, revision: str) -> dict:
    require(bool(re.fullmatch(r"[0-9a-f]{40}", revision)), "Use a full candidate commit SHA")
    require(git(root, "rev-parse", "HEAD") == revision, "Checkout is not the selected candidate")
    require(not git(root, "status", "--porcelain"),
            "Candidate checkout must be clean; retain evidence outside the checkout")
    return {
        "source_revision": revision,
        "cargo_lock_sha256": digest(root / "Cargo.lock"),
        "tutorial_manifest_sha256": digest(root / "docs/tutorial/manifest.json"),
        "tutorial_catalog_sha256": digest(root / "docs/tutorial/catalog.json"),
    }


def binary_identity(path: Path, revision: str) -> dict:
    before = digest(path)
    result = subprocess.run([str(path), "--version"], capture_output=True, timeout=20)
    text = result.stdout.decode("utf-8", errors="replace")
    matches = re.findall(r"^Source revision \S+\+git\.([0-9a-f]{40})$", text, re.MULTILINE)
    require(result.returncode == 0 and matches == [revision],
            f"Binary does not report the full candidate revision: {path}")
    require(before == digest(path), f"Binary changed during identity probe: {path}")
    return {"path": str(path), "sha256": before, "version_output": text}


def coverage_inventory(root: Path) -> dict:
    catalog_path = root / "docs/tutorial/catalog.json"
    manifest_path = root / "docs/tutorial/manifest.json"
    catalog, manifest = load(catalog_path), load(manifest_path)
    require(catalog.get("schema") == "gentle.tutorial_catalog.v2", "Unsupported catalog")
    require(manifest.get("schema") == "gentle.tutorial_manifest.v2", "Unsupported manifest")
    chapters = {row["id"]: row for row in manifest["chapters"]}
    entries = {row["id"]: row for row in catalog["entries"]}
    require(len(chapters) == len(manifest["chapters"]), "Duplicate chapter IDs")
    require(len(entries) == len(catalog["entries"]), "Duplicate catalog IDs")
    workflows = {}
    for path in sorted((root / "docs/examples/workflows").glob("*.json")):
        workflow = load(path)
        require(workflow["id"] not in workflows, "Duplicate workflow IDs")
        workflows[workflow["id"]] = (path, workflow)
    rows = []
    for identifier in sorted(entries.keys() | chapters.keys()):
        entry, chapter = entries.get(identifier, {}), chapters.get(identifier, {})
        contract = chapter.get("gui_acceptance")
        workflow_id = chapter.get("example_id")
        workflow = None
        if workflow_id:
            require(workflow_id in workflows, f"Missing workflow: {workflow_id}")
            path, definition = workflows[workflow_id]
            workflow = {
                "id": workflow_id, "source": file_record(path),
                "test_mode": definition.get("test_mode"),
                "required_files": [
                    {"path": name, "available": (root / name).is_file()}
                    for name in definition.get("required_files", [])
                ],
            }
        if contract:
            mode = "gui_view_only" if contract.get("view_only", False) else "gui_scientific"
            reason = "Typed GUI contract; execution still required on the selected candidate"
        elif chapter:
            mode, reason = "workflow_only", "No typed GUI acceptance contract"
        elif entry.get("type") in ("operational_reference", "executable_collection"):
            mode, reason = "reference", "Reference/collection entry, not a GUI execution contract"
        else:
            mode, reason = "manual_uncovered", "Manual walkthrough without a typed GUI contract"
        rows.append({
            "id": identifier, "title": entry.get("title", chapter.get("title")),
            "catalog_type": entry.get("type"), "path": entry.get("path"),
            "chapter_present": bool(chapter), "coverage_kind": mode, "reason": reason,
            "workflow": workflow, "gui_profile": contract.get("profile") if contract else None,
            "acceptance_contract_sha256": content_digest(contract) if contract else None,
            "execution_status": "not_run", "human_review_status": entry.get("review_status"),
        })
    return {
        "schema": COVERAGE_SCHEMA, "catalog": file_record(catalog_path),
        "manifest": file_record(manifest_path), "catalog_entries": len(entries),
        "manifest_chapters": len(chapters),
        "gui_contracts": sum(bool(c.get("gui_acceptance")) for c in chapters.values()),
        "rows": rows,
        "non_claim": "Coverage declarations and available files are not test passes.",
    }


def verify_gui_run(root: Path, candidate: dict) -> dict:
    """Reject mixed, incomplete or altered runs before offering publication choices."""
    require(candidate.get("schema") == CANDIDATE_SCHEMA, "Unsupported candidate binding")
    require(bool(re.fullmatch(r"[0-9a-f]{40}", candidate.get("source_revision", ""))),
            "Candidate needs a full revision")
    require(set(candidate.get("binaries", {})) == {"gentle", "gentle_cli", "gentle_examples_docs"},
            "Candidate must bind all three binaries")
    report = load(root / "acceptance-report.json")
    require(report.get("schema") == "gentle.tutorial_gui_acceptance_run.v1"
            and report.get("status") == "pass", "GUI run is not a passing acceptance report")
    expected = candidate["chapters"]
    rows = report.get("chapters", [])
    require(bool(expected) and [r["chapter_id"] for r in rows] == [c["id"] for c in expected],
            "GUI run does not contain the exact ordered required chapter set")
    require(report.get("chapter_count") == len(expected), "Incorrect chapter count")
    environment = load(root / "environment.json")
    require(environment.get("schema") == "gentle.tutorial_acceptance_environment.v1",
            "Unsupported GUI environment schema")
    require(environment.get("source_revision") == candidate["source_revision"],
            "GUI environment revision mismatch")
    require(environment.get("cargo_lock_sha256") == candidate["cargo_lock_sha256"],
            "GUI environment lockfile mismatch")
    require(not environment.get("git_status"), "GUI evidence describes a dirty checkout")
    namespace = environment.get("network_namespace", {})
    require(environment.get("network_enforcement") == "linux_network_namespace"
            and all(re.fullmatch(r"net:\[\d+\]", namespace.get(key) or "") for key in ("current", "parent"))
            and namespace["current"] != namespace["parent"], "GUI run has no verified offline namespace")
    for name, binary in candidate["binaries"].items():
        require(environment.get("binaries", {}).get(name, {}).get("sha256") == binary["sha256"],
                f"Mixed {name} binary in GUI evidence")
    selected = []
    for chapter, row in zip(expected, rows):
        identifier, contract = chapter["id"], chapter["gui_acceptance"]
        ledger_path = child_path(root, f"{identifier}/acceptance-ledger.json")
        require(digest(ledger_path) == row["ledger_sha256"], "GUI ledger hash mismatch")
        ledger = load(ledger_path)
        require(ledger.get("schema") == "gentle.tutorial_gui_acceptance_ledger.v1"
                and ledger.get("chapter_id") == identifier, "GUI ledger identity mismatch")
        require(row.get("status") == ledger.get("status") == "pass", "Chapter did not pass")
        require(ledger.get("environment") == environment, "Mixed chapter environments")
        require(ledger.get("manifest_sha256") == candidate["tutorial_manifest_sha256"],
                "Tutorial manifest mismatch")
        require(ledger.get("acceptance_contract_sha256") == content_digest(contract),
                "GUI contract mismatch")
        require(ledger.get("preflight", {}).get("status") == "pass"
                and ledger.get("completion", {}).get("truth") == "satisfied",
                "Missing starter/oracle or completion verification")
        require(ledger["preflight"].get("starter_completion_truth") == (
            "satisfied" if contract.get("view_only", False) else "unsatisfied"
        ) and ledger["preflight"].get("oracle_completion_truth") == "satisfied",
                "Starter was not distinct from the completed scientific oracle")
        for phase in ("starter", "oracle"):
            prepared = ledger.get(phase, {})
            require(prepared.get("example_id") == contract[phase]["example_id"]
                    and prepared.get("workflow_source_sha256") ==
                    candidate["workflow_bindings"][identifier][phase]["sha256"],
                    "Starter/oracle workflow identity mismatch")
        require(digest(child_path(root, f"{identifier}/starter.project.gentle.json"))
                == ledger.get("final_project_sha256"), "Saved final project hash mismatch")
        require(digest(child_path(root, f"{identifier}/oracle.project.gentle.json"))
                == ledger["oracle"].get("project_sha256"), "Saved oracle hash mismatch")
        steps = ledger.get("steps", [])
        require([s["id"] for s in steps] == [s["id"] for s in contract["steps"]],
                "Missing, reordered or extra GUI steps")
        checkpoints = []
        for step, actual in zip(contract["steps"], steps):
            require(actual.get("status") == "pass"
                    and actual.get("step_sha256") == content_digest(step), "Step binding mismatch")
            require(actual.get("requested_action") == step["interaction"]
                    and actual.get("semantic_target") == step["target"]
                    and actual.get("window") == step["window"]
                    and actual.get("resolved_target", {}).get("semantic_id") == step["target"],
                    "Recorded GUI action/target differs from the contract")
            verifiers = actual.get("verifiers", [])
            require([v.get("kind") for v in verifiers] == [v["kind"] for v in step.get("verifiers", [])]
                    and all(v.get("status") == "pass" for v in verifiers),
                    "Missing or failed step verifiers")
            if step.get("scientific_effect"):
                require(actual.get("after_fact", {}).get("truth") == "satisfied",
                        "Missing scientific after-fact")
                allowed = {"unsatisfied", "satisfied"} if step.get("allow_preexisting") else {"unsatisfied"}
                require(actual.get("before_fact", {}).get("truth") in allowed,
                        "Missing scientific before-fact")
            require(actual.get("x11_event") is not None
                    and actual.get("resolved_target") is not None
                    and actual.get("after_generation", -1) > actual.get("before_generation", -1),
                    "Missing GUI interaction/generation receipt")
            policy, evidence = step.get("evidence", {}), actual.get("evidence", {})
            def checkpoint_path(suffix: str) -> Path:
                return child_path(root, f"{identifier}/checkpoints/{step['id']}.{suffix}")

            snapshot_path = checkpoint_path("snapshot.json")
            if policy.get("snapshot") == "required" or "snapshot" in evidence:
                require(digest(snapshot_path) == evidence.get("snapshot", {}).get("sha256"),
                        "Snapshot hash mismatch")
            if policy.get("screenshot") == "required" or "screenshot" in evidence:
                screenshot = evidence.get("screenshot", {})
                record_path = checkpoint_path("screenshot.json")
                require(digest(record_path) == screenshot.get("record_sha256"),
                        "Screenshot sidecar hash mismatch")
                record = load(record_path)
                bindings = {
                    "schema": "gentle.tutorial_gui_screenshot_evidence.v1",
                    "source_revision": candidate["source_revision"],
                    "gentle_binary_sha256": candidate["binaries"]["gentle"]["sha256"],
                    "tutorial_manifest_sha256": candidate["tutorial_manifest_sha256"],
                    "acceptance_contract_sha256": content_digest(contract),
                    "chapter_id": identifier, "step_id": step["id"],
                    "prose_step": step["prose_step"],
                    "step_sha256": content_digest(step),
                }
                require(all(record.get(k) == v for k, v in bindings.items()),
                        "Screenshot does not belong to this candidate/contract/step")
                require(record.get("requested_target", {}).get("semantic_id") == step["target"]
                        and record.get("requested_target", {}).get("window_id") == step["window"],
                        "Screenshot requested target mismatch")
                snapshot = load(snapshot_path)
                require(record["semantic_snapshot"]["retained_file_sha256"] == digest(snapshot_path)
                        and record["semantic_snapshot"]["canonical_sha256"] == content_digest(snapshot),
                        "Screenshot semantic snapshot mismatch")
                require(digest(checkpoint_path("raw.png")) == record["capture"]["raw"]["sha256"],
                        "Screenshot raw image mismatch")
                views = record.get("derived_views", [])
                require({v["role"] for v in views} == {"orientation", "interaction_context"}
                        and len(views) == 2, "Missing screenshot teaching views")
                for view in views:
                    suffix = "orientation" if view["role"] == "orientation" else "context"
                    require(digest(checkpoint_path(f"{suffix}.svg")) == view["sha256"],
                            "Screenshot teaching-view hash mismatch")
                checkpoints.append(step["id"])
        selected.append({"chapter_id": identifier, "checkpoints": checkpoints})
    return {
        "schema": "gentle.tutorial_gui_screenshot_selection.v1",
        "source_revision": candidate["source_revision"],
        "chapters": selected,
    }
