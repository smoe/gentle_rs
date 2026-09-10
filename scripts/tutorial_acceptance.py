#!/usr/bin/env python3
"""Coordinate candidate-bound offline tutorial acceptance for Glen or another auditor.

No fetching, building, publication or repair is implicit. Run inside the same
isolated Linux/X11 session as tutorial_gui_acceptance.py. Inventory works offline
on any host and does not execute GENtle or claim acceptance.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

if __package__:
    from . import tutorial_gui_acceptance as gui
    from . import tutorial_acceptance_evidence as evidence
else:
    import tutorial_gui_acceptance as gui
    import tutorial_acceptance_evidence as evidence


RUN_SCHEMA = "gentle.tutorial_acceptance.v1"
CHECKS = ("tutorial-catalog-check", "tutorial-manifest-check", "tutorial-check")


def select_chapters(manifest: dict, profiles: list[str]) -> list[dict]:
    selected = []
    for profile in dict.fromkeys(profiles):
        selected.extend(gui.selected_chapters(manifest, [], profile))
    evidence.require(all(c["gui_acceptance"].get("network") == "offline" for c in selected),
                     "This gate is offline-only; online/private studies require separate authorization")
    return selected


def isolated_environment(directory: Path) -> dict:
    _, cleared = gui.redacted_environment(dict(os.environ))
    env = {key: value for key, value in os.environ.items() if key not in cleared}
    for key in ("HOME", "XDG_CONFIG_HOME", "XDG_CACHE_HOME", "XDG_DATA_HOME", "TMPDIR"):
        path = directory / key.lower()
        path.mkdir(parents=True)
        env[key] = str(path)
    env.update(LANG="C.UTF-8", LC_ALL="C.UTF-8", TZ="UTC",
               GENTLE_SKIP_REMOTE_TESTS="1", CARGO_NET_OFFLINE="true")
    return env


def execute(argv: list[str], root: Path, env: dict, folder: Path, timeout: int) -> dict:
    """Keep exact output bytes and bounded process-group cleanup in the receipt."""
    folder.mkdir(parents=True)
    stdout, stderr = folder / "stdout", folder / "stderr"
    receipt = {"argv": argv, "timeout_seconds": timeout, "status": "not_run"}
    started = time.monotonic()
    process = None
    with stdout.open("wb") as out, stderr.open("wb") as err:
        try:
            process = subprocess.Popen(argv, cwd=root, env=env, stdout=out, stderr=err,
                                       start_new_session=True)
            receipt["exit_code"] = process.wait(timeout=timeout)
            receipt["status"] = "pass" if process.returncode == 0 else "fail"
        except OSError as error:
            receipt.update(status="missing_dependency" if isinstance(error, (FileNotFoundError, PermissionError))
                           else "harness_gap", message=str(error))
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            receipt.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "timeout")
        finally:
            if process is not None and process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                except ProcessLookupError:
                    process.wait()
                receipt["exit_code"] = process.returncode
    receipt.update(elapsed_ms=round((time.monotonic() - started) * 1000),
                   stdout=evidence.file_record(stdout), stderr=evidence.file_record(stderr))
    gui.atomic_write_json(folder / "receipt.json", receipt)
    return receipt


def bind_candidate(root: Path, revision: str, binaries: dict, chapters: list[dict]) -> dict:
    identity = evidence.checkout_identity(root, revision)
    harness_sources = []
    for module in (sys.modules[__name__], evidence, gui):
        running = Path(module.__file__).resolve()
        frozen = root / "scripts" / running.name
        evidence.require(evidence.digest(running) == evidence.digest(frozen),
                         "Running harness is not from the frozen candidate")
        harness_sources.append(evidence.file_record(frozen))
    inventory = evidence.coverage_inventory(root)
    workflows = {
        evidence.load(path)["id"]: path
        for path in (root / "docs/examples/workflows").glob("*.json")
    }
    inputs = set()
    workflow_bindings = {}
    for chapter in chapters:
        workflow_bindings[chapter["id"]] = {}
        for phase in ("starter", "oracle"):
            path = workflows[chapter["gui_acceptance"][phase]["example_id"]]
            inputs.add(path)
            workflow_bindings[chapter["id"]][phase] = evidence.file_record(path)
            inputs.update(root / name for name in evidence.load(path).get("required_files", []))
    return {
        "schema": evidence.CANDIDATE_SCHEMA, **identity,
        "binaries": {name: evidence.binary_identity(path, revision) for name, path in binaries.items()},
        "chapters": [{"id": c["id"], "gui_acceptance": c["gui_acceptance"]} for c in chapters],
        "inputs": [evidence.file_record(path) for path in sorted(inputs)],
        "workflow_bindings": workflow_bindings,
        "harness_sources": harness_sources,
        "coverage_sha256": evidence.content_digest(inventory),
        "scope": "Linux/X11 offline tutorial acceptance; not whole-release or biological acceptance",
    }


def recheck_candidate(root: Path, candidate: dict) -> None:
    current = evidence.checkout_identity(root, candidate["source_revision"])
    evidence.require(all(candidate[key] == value for key, value in current.items()),
                     "Candidate source inputs changed during acceptance")
    for row in [*candidate["inputs"], *candidate["binaries"].values(), *candidate.get("harness_sources", [])]:
        evidence.require(evidence.digest(Path(row["path"])) == row["sha256"],
                         f"Bound input changed during acceptance: {row['path']}")


def run(args: argparse.Namespace) -> dict:
    root, output = args.repo_root.resolve(), args.evidence_dir.resolve()
    evidence.require(not output.exists(), "Use a new evidence directory; prior runs are immutable")
    evidence.require(not output.is_relative_to(root), "Keep acceptance evidence outside the checkout")
    output.mkdir(parents=True)
    report = {
        "schema": RUN_SCHEMA, "status": "incomplete", "candidate_revision": args.candidate,
        "scope": "Linux/X11 offline tutorial acceptance", "profiles": args.profile or ["smoke"],
        "checks": [{"id": name, "status": "not_run"} for name in (*CHECKS, "gui", "evidence-binding")],
        "scientific_study_acceptance": "not_run", "publication": "not_requested",
    }
    report_path = output / "tutorial-acceptance-report.json"
    gui.atomic_write_json(report_path, report)
    try:
        if args.supersedes:
            previous = evidence.load(args.supersedes)
            evidence.require(previous.get("schema") == RUN_SCHEMA, "Not a prior tutorial acceptance report")
            report["supersedes"] = evidence.file_record(args.supersedes.resolve())
        inventory = evidence.coverage_inventory(root)
        gui.atomic_write_json(output / "coverage.json", inventory)
        report["coverage"] = evidence.file_record(output / "coverage.json")
        manifest = evidence.load(root / "docs/tutorial/manifest.json")
        chapters = select_chapters(manifest, report["profiles"])
        report["required_chapters"] = [c["id"] for c in chapters]
        report["chapter_results"] = [{"chapter_id": c["id"], "status": "not_run"} for c in chapters]
        candidate = bind_candidate(root, args.candidate, {
            name: (root / getattr(args, flag)).resolve()
            for name, flag in (("gentle", "gentle"), ("gentle_cli", "gentle_cli"),
                               ("gentle_examples_docs", "examples_docs"))
        }, chapters)
        gui.atomic_write_json(output / "candidate.json", candidate)
        report["candidate"] = evidence.file_record(output / "candidate.json")
        if platform.system() != "Linux" or not os.environ.get("DISPLAY"):
            raise gui.AcceptanceFailure("missing_dependency", "Live acceptance needs Linux/X11 with DISPLAY")
        gui.validate_parent_network_namespace(os.readlink("/proc/self/ns/net"), args.parent_network_namespace)
        env = isolated_environment(output / "check-profile")
        redacted, _ = gui.redacted_environment(dict(os.environ))
        report["cleared_inherited_variables"] = redacted
        helper = candidate["binaries"]["gentle_examples_docs"]["path"]
        for check in report["checks"][:3]:
            recheck_candidate(root, candidate)
            check.update(execute([helper, check["id"]], root, env,
                                 output / "checks" / check["id"], timeout=1200))
            gui.atomic_write_json(report_path, report)
            if check["status"] == "interrupted":
                raise KeyboardInterrupt
        if any(c["status"] != "pass" for c in report["checks"][:3]):
            raise gui.AcceptanceFailure("product_failure", "Tutorial source/artifact checks did not all pass")
        argv = [sys.executable, str(root / "scripts/tutorial_gui_acceptance.py"),
                "--repo-root", str(root), "--evidence-dir", str(output / "gui"),
                "--network-enforcement", "linux_network_namespace",
                "--parent-network-namespace", args.parent_network_namespace]
        for name, flag in (("gentle", "--gentle"), ("gentle_cli", "--gentle-cli"),
                           ("gentle_examples_docs", "--examples-docs")):
            argv.extend([flag, candidate["binaries"][name]["path"]])
        for chapter in chapters:
            argv.extend(["--chapter", chapter["id"]])
        # Each chapter owns per-step budgets; the outer bound only detects a stuck harness.
        timeout = 1200 * len(chapters) + sum(
            6 * gui.TIMEOUT_DEFAULTS[s["timeout_class"]]
            for c in chapters for s in c["gui_acceptance"]["steps"]
        )
        report["checks"][3].update(execute(argv, root, env, output / "checks/gui", int(timeout)))
        gui_report = output / "gui/acceptance-report.json"
        if gui_report.is_file():
            report["gui_report"] = evidence.file_record(gui_report)
            by_id = {row["chapter_id"]: row for row in evidence.load(gui_report).get("chapters", [])}
            report["chapter_results"] = [by_id.get(row["chapter_id"], row) for row in report["chapter_results"]]
            if report["checks"][3]["status"] != "pass":
                classes = {row.get("failure_class") for row in by_id.values() if row.get("status") != "pass"}
                report["checks"][3]["failure_class"] = evidence.load(gui_report).get("failure_class") or (
                    next(iter(classes)) if len(classes) == 1 else None
                )
        if report["checks"][3]["status"] != "pass":
            failure = report["checks"][3].get("failure_class") or "gui_acceptance_failed"
            if report["checks"][3]["status"] == "interrupted":
                failure = "interrupted"
            raise gui.AcceptanceFailure(failure, "See retained GUI/checkpoint receipts")
        report["checks"][4]["status"] = "running"
        recheck_candidate(root, candidate)
        selection = evidence.verify_gui_run(output / "gui", candidate)
        selection["capture_date"] = datetime.datetime.now(datetime.timezone.utc).date().isoformat()
        gui.atomic_write_json(output / "selection.json", selection)
        report["screenshot_selection"] = evidence.file_record(output / "selection.json")
        report["checks"][4]["status"] = "pass"
        report["status"] = "pass"
    except (ValueError, KeyError, TypeError, OSError, subprocess.SubprocessError,
            gui.AcceptanceFailure, KeyboardInterrupt) as error:
        failure = getattr(error, "failure_class", "evidence_mismatch")
        if isinstance(error, FileNotFoundError):
            failure = "missing_dependency"
        elif isinstance(error, KeyboardInterrupt):
            failure = "interrupted"
        report.update(failure_class=failure, message=str(error))
        for check in report["checks"]:
            if check["status"] == "running":
                check.update(status="fail", failure_class=failure, message=str(error))
        report["status"] = "incomplete" if failure in (
            "missing_dependency", "harness_gap", "interrupted"
        ) else "fail"
    finally:
        gui.atomic_write_json(report_path, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("inventory", "run", "verify"))
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, help="Inventory JSON; stdout when omitted")
    parser.add_argument("--candidate", help="Full frozen commit SHA")
    parser.add_argument("--candidate-binding", type=Path, help="candidate.json for offline evidence verification")
    parser.add_argument("--profile", action="append", choices=("smoke", "offline-core"))
    parser.add_argument("--evidence-dir", type=Path)
    parser.add_argument("--parent-network-namespace")
    parser.add_argument("--gentle", type=Path, default=Path("target/debug/gentle"))
    parser.add_argument("--gentle-cli", type=Path, default=Path("target/debug/gentle_cli"))
    parser.add_argument("--examples-docs", type=Path, default=Path("target/debug/gentle_examples_docs"))
    parser.add_argument("--supersedes", type=Path, help="Prior report to reference, never overwrite or reuse as a pass")
    args = parser.parse_args()
    try:
        if args.command == "inventory":
            result = evidence.coverage_inventory(args.repo_root.resolve())
            if args.output:
                gui.atomic_write_json(args.output, result)
        elif args.command == "verify":
            evidence.require(args.candidate_binding is not None and args.evidence_dir is not None
                             and args.candidate is not None,
                             "verify requires --candidate, --candidate-binding and --evidence-dir (GUI run)")
            candidate = evidence.load(args.candidate_binding)
            evidence.require(candidate["source_revision"] == args.candidate, "Candidate revision mismatch")
            selection = evidence.verify_gui_run(args.evidence_dir, candidate)
            result = {"status": "pass", "scope": "retained GUI evidence binding",
                      "candidate": evidence.file_record(args.candidate_binding), "selection": selection}
        else:
            evidence.require(args.candidate is not None and args.evidence_dir is not None,
                             "run requires --candidate and --evidence-dir")
            result = run(args)
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0 if result.get("status", "pass") == "pass" else 1
    except (ValueError, KeyError, TypeError, OSError, gui.AcceptanceFailure) as error:
        print(json.dumps({"status": "fail", "message": str(error)}), file=sys.stderr)
        return 1


if __name__ == "__main__":
    def interrupt(_signal: int, _frame: object) -> None:
        raise KeyboardInterrupt

    for sig in (signal.SIGTERM, getattr(signal, "SIGHUP", None)):
        if sig is not None:
            signal.signal(sig, interrupt)
    raise SystemExit(main())
