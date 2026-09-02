#!/usr/bin/env python3
"""Drive a typed GENtle tutorial GUI contract through ordinary X11 events.

The runner deliberately has no synthetic GUI input API. Rust validates the
closed semantic-control contract; this process reads only the opt-in semantic
snapshot, translates rectangles into xdotool events, saves through Ctrl+S, and
checks scientific outcomes against the persisted project with gentle_cli.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import signal
import struct
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable, Iterable, Optional


LEDGER_SCHEMA = "gentle.tutorial_gui_acceptance_run.v1"
SNAPSHOT_SCHEMA = "gentle.gui_semantic_snapshot.v2"
CONTRACT_SCHEMA = "gentle.tutorial_gui_acceptance.v1"
TIMEOUT_SECONDS = {
    "instant": 3.0,
    "interactive": 30.0,
    "compute": 180.0,
    "io": 90.0,
}


class AuditFailure(RuntimeError):
    def __init__(self, result: str, message: str):
        super().__init__(message)
        self.result = result


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AuditFailure("contract_error", f"Could not read JSON '{path}': {error}") from error


def subject_scope(subject: dict[str, str]) -> Optional[str]:
    if not subject:
        return None
    if set(subject) != {"sequence"} or not subject["sequence"]:
        raise AuditFailure(
            "contract_error",
            "The X11 runner currently supports only a non-empty 'sequence' subject; "
            f"received keys {sorted(subject)}",
        )
    identity = bytearray(b"gentle.gui.subject_scope.v1\0")
    part = subject["sequence"].encode("utf-8")
    identity.extend(struct.pack(">Q", len(part)))
    identity.extend(part)
    return "subject-" + sha256_bytes(bytes(identity))[:32]


def json_path(value: Any, path: str) -> Any:
    current = value
    for segment in path.split("."):
        if isinstance(current, list):
            try:
                current = current[int(segment)]
            except (ValueError, IndexError) as error:
                raise KeyError(path) from error
        elif isinstance(current, dict) and segment in current:
            current = current[segment]
        else:
            raise KeyError(path)
    return current


def is_non_empty(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, (str, list, dict, tuple, set)):
        return len(value) > 0
    return True


def compare_numbers(left: Any, op: str, right: Any) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return False
    if not isinstance(left, (int, float)) or not isinstance(right, (int, float)):
        return False
    comparisons = {
        "eq": lambda: left == right,
        "=": lambda: left == right,
        "==": lambda: left == right,
        "ne": lambda: left != right,
        "!=": lambda: left != right,
        "gt": lambda: left > right,
        ">": lambda: left > right,
        "gte": lambda: left >= right,
        ">=": lambda: left >= right,
        "lt": lambda: left < right,
        "<": lambda: left < right,
        "lte": lambda: left <= right,
        "<=": lambda: left <= right,
    }
    return op in comparisons and comparisons[op]()


def assert_report(report: Any, verifier: dict[str, Any]) -> list[str]:
    checks: list[str] = []
    if report.get("schema") != verifier["schema"]:
        raise AuditFailure(
            "product_failure",
            f"Report '{verifier['report_id']}' schema is {report.get('schema')!r}, "
            f"expected {verifier['schema']!r}",
        )
    for path in verifier.get("required_fields", []):
        try:
            json_path(report, path)
        except KeyError as error:
            raise AuditFailure(
                "product_failure", f"Report is missing required field '{path}'"
            ) from error
        checks.append(f"required:{path}")
    for assertion in verifier.get("assertions", []):
        kind = assertion["kind"]
        if kind == "value":
            path = assertion["path"]
            try:
                actual = json_path(report, path)
            except KeyError as error:
                raise AuditFailure(
                    "product_failure", f"Report assertion field '{path}' is missing"
                ) from error
            if "equals" in assertion and actual != assertion["equals"]:
                raise AuditFailure(
                    "product_failure",
                    f"Report field '{path}' is {actual!r}, expected {assertion['equals']!r}",
                )
            if assertion.get("non_empty") and not is_non_empty(actual):
                raise AuditFailure("product_failure", f"Report field '{path}' is empty")
            if "compare" in assertion:
                comparison = assertion["compare"]
                if not compare_numbers(actual, comparison["op"], comparison["value"]):
                    raise AuditFailure(
                        "product_failure",
                        f"Report field '{path}' value {actual!r} failed "
                        f"{comparison['op']} {comparison['value']!r}",
                    )
            checks.append(f"value:{path}")
        elif kind == "relation":
            left = json_path(report, assertion["left_path"])
            right = json_path(report, assertion["right_path"])
            if not compare_numbers(left, assertion["op"], right):
                raise AuditFailure(
                    "product_failure",
                    f"Report relation {assertion['left_path']}={left!r} "
                    f"{assertion['op']} {assertion['right_path']}={right!r} failed",
                )
            checks.append(
                f"relation:{assertion['left_path']}:{assertion['op']}:{assertion['right_path']}"
            )
        else:
            raise AuditFailure("contract_error", f"Unsupported report assertion kind '{kind}'")
    return checks


class TutorialGuiRunner:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.repo_root = args.repo_root.resolve()
        self.output_dir = args.output_dir.resolve()
        self.ledger_path = self.output_dir / "acceptance-ledger.json"
        self.logs_dir = self.output_dir / "logs"
        self.checkpoints_dir = self.output_dir / "checkpoints"
        self.runtime_dir = self.output_dir / "runtime"
        self.projects_dir = self.output_dir / "projects"
        self.queries_dir = self.output_dir / "queries"
        self.snapshot_path = self.runtime_dir / "semantic-snapshot.json"
        self.command_index = 0
        self.gui_process: Optional[subprocess.Popen[str]] = None
        self.screenshot_tool: Optional[str] = None
        self.evidence_initialized = False
        self.ledger: dict[str, Any] = {
            "schema": LEDGER_SCHEMA,
            "tutorial_id": args.tutorial_id,
            "result": "running",
            "acceptance_executed": not args.validate_only,
            "started_at": utc_now(),
            "finished_at": None,
            "source": {},
            "environment": {},
            "projects": {},
            "commands": [],
            "checkpoints": [],
            "failure": None,
        }

    def fail(self, result: str, message: str) -> None:
        raise AuditFailure(result, message)

    def write_ledger(self) -> None:
        if not self.evidence_initialized:
            return
        payload = copy.deepcopy(self.ledger)
        payload.pop("payload_sha256", None)
        self.ledger["payload_sha256"] = sha256_bytes(canonical_json_bytes(payload))
        atomic_write_json(self.ledger_path, self.ledger)

    def relative(self, path: Path) -> str:
        try:
            return path.relative_to(self.output_dir).as_posix()
        except ValueError:
            return path.name

    def command(
        self,
        label: str,
        argv: list[str],
        *,
        cwd: Optional[Path] = None,
        env: Optional[dict[str, str]] = None,
        timeout: float = 300.0,
        failure_result: str = "product_failure",
    ) -> subprocess.CompletedProcess[str]:
        self.command_index += 1
        safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
        stdout_path = self.logs_dir / f"{self.command_index:02d}-{safe_label}.stdout.txt"
        stderr_path = self.logs_dir / f"{self.command_index:02d}-{safe_label}.stderr.txt"
        started = time.monotonic()
        try:
            completed = subprocess.run(
                argv,
                cwd=str(cwd or self.repo_root),
                env=env,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            self.fail(failure_result, f"Could not run {label}: {error}")
        stdout_path.write_text(completed.stdout, encoding="utf-8")
        stderr_path.write_text(completed.stderr, encoding="utf-8")
        record = {
            "label": label,
            "program": Path(argv[0]).name,
            "arguments": [Path(argv[0]).name, *argv[1:]],
            "cwd_role": "repo" if (cwd or self.repo_root) == self.repo_root else "isolated_run",
            "exit_code": completed.returncode,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "stdout": self.relative(stdout_path),
            "stderr": self.relative(stderr_path),
        }
        self.ledger["commands"].append(record)
        self.write_ledger()
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic"
            self.fail(failure_result, f"{label} failed ({completed.returncode}): {detail}")
        return completed

    def command_json(self, label: str, argv: list[str], **kwargs: Any) -> Any:
        completed = self.command(label, argv, **kwargs)
        try:
            return json.loads(completed.stdout)
        except json.JSONDecodeError as error:
            self.fail("harness_gap", f"{label} did not emit one JSON value: {error}")

    def tool_version(self, argv: list[str]) -> str:
        try:
            completed = subprocess.run(
                argv,
                cwd=str(self.repo_root),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=15.0,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            return f"unavailable: {error}"
        return completed.stdout.strip().splitlines()[0] if completed.stdout.strip() else "unknown"

    def prepare_directories(self) -> None:
        if self.output_dir.exists() and any(self.output_dir.iterdir()):
            self.fail(
                "harness_gap",
                f"Refusing to overwrite non-empty audit directory '{self.output_dir}'",
            )
        for path in (
            self.logs_dir,
            self.checkpoints_dir,
            self.runtime_dir,
            self.projects_dir,
            self.queries_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)
        self.evidence_initialized = True

    def require_file(self, path: Path, label: str) -> None:
        if not path.is_file():
            self.fail("harness_gap", f"{label} does not exist: '{path}'")

    def establish_provenance(self) -> None:
        manifest = self.args.manifest.resolve()
        source_dir = self.args.tutorial_sources.resolve()
        examples_dir = self.args.examples.resolve()
        for path, label in (
            (manifest, "tutorial manifest"),
            (self.args.gentle_cli.resolve(), "gentle_cli"),
            (self.args.examples_bin.resolve(), "gentle_examples_docs"),
        ):
            self.require_file(path, label)
        if not self.args.validate_only:
            self.require_file(self.args.gentle_bin.resolve(), "GENtle GUI binary")
        revision = self.command(
            "git-revision", ["git", "rev-parse", "HEAD"], failure_result="harness_gap"
        ).stdout.strip()
        self.ledger["source"] = {
            "revision": revision,
            "manifest_sha256": sha256_file(manifest),
            "tutorial_sources_role": source_dir.name,
            "examples_role": examples_dir.name,
            "binaries": {
                "gentle_cli": sha256_file(self.args.gentle_cli.resolve()),
                "gentle_examples_docs": sha256_file(self.args.examples_bin.resolve()),
                "gentle": None
                if self.args.validate_only
                else sha256_file(self.args.gentle_bin.resolve()),
            },
        }
        self.ledger["environment"] = {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "gentle_cli": self.tool_version([str(self.args.gentle_cli.resolve()), "--version"]),
            "gentle": None
            if self.args.validate_only
            else self.tool_version([str(self.args.gentle_bin.resolve()), "--version"]),
            "display_present": bool(os.environ.get("DISPLAY")),
            "network_policy": "contract-recorded; OS-level isolation remains auditor-owned",
            "profile_isolation": "temporary HOME and XDG directories; inherited GENTLE_* and secret variables removed",
        }
        self.write_ledger()

    def validate_manifest(self) -> None:
        self.command(
            "tutorial-manifest-check",
            [
                str(self.args.examples_bin.resolve()),
                "tutorial-manifest-check",
                "--catalog-meta",
                str(self.args.catalog_meta.resolve()),
                "--tutorial-sources",
                str(self.args.tutorial_sources.resolve()),
                "--manifest",
                str(self.args.manifest.resolve()),
            ],
            failure_result="contract_error",
        )

    def load_contract(self) -> tuple[dict[str, Any], dict[str, Any]]:
        manifest = load_json(self.args.manifest.resolve())
        chapters = [
            chapter
            for chapter in manifest.get("chapters", [])
            if chapter.get("id") == self.args.tutorial_id
        ]
        if len(chapters) != 1:
            self.fail(
                "contract_error",
                f"Expected one tutorial chapter '{self.args.tutorial_id}', found {len(chapters)}",
            )
        chapter = chapters[0]
        contract = chapter.get("gui_acceptance")
        if not isinstance(contract, dict):
            self.fail("contract_error", "Tutorial chapter has no GUI acceptance contract")
        if contract.get("schema") != CONTRACT_SCHEMA:
            self.fail(
                "contract_error",
                f"Unsupported GUI acceptance schema {contract.get('schema')!r}",
            )
        self.ledger["source"]["contract_sha256"] = sha256_bytes(
            canonical_json_bytes(contract)
        )
        self.ledger["contract"] = {
            "schema": contract["schema"],
            "profile": contract["profile"],
            "network": contract["network"],
            "starter_example_id": contract["starter"]["example_id"],
            "oracle_example_id": contract["oracle"]["example_id"],
            "step_count": len(contract["steps"]),
        }
        self.write_ledger()
        return chapter, contract

    def materialize_project(self, role: str, example_id: str) -> Path:
        project = self.projects_dir / f"{role}.project.gentle.json"
        run_dir = self.output_dir / role
        result = self.command_json(
            f"materialize-{role}",
            [
                str(self.args.examples_bin.resolve()),
                "example-project",
                example_id,
                str(project),
                "--run-dir",
                str(run_dir),
                "--source",
                str(self.args.examples.resolve()),
                "--repo-root",
                str(self.repo_root),
            ],
            failure_result="product_failure",
            timeout=600.0,
        )
        materializer = {
            key: result.get(key)
            for key in ("status", "mode", "example_id", "sequence_count")
        }
        self.ledger["projects"][role] = {
            "example_id": example_id,
            "project": self.relative(project),
            "sha256": sha256_file(project),
            "materializer": materializer,
        }
        self.write_ledger()
        return project

    def shell_json(self, project: Path, label: str, tokens: Iterable[str]) -> Any:
        token_list = list(tokens)
        for token in token_list:
            if not token or any(ord(character) < 32 or ord(character) == 127 for character in token):
                self.fail("contract_error", f"Unsafe token in typed verifier: {token!r}")
        command_text = shlex.join(token_list)
        return self.command_json(
            label,
            [str(self.args.gentle_cli.resolve()), "--project", str(project), "shell", command_text],
            timeout=180.0,
        )

    def fact_evaluation(self, project: Path, expression: dict[str, Any], label: str) -> Any:
        digest = sha256_bytes(canonical_json_bytes(expression))[:16]
        query = self.queries_dir / f"fact-{digest}.json"
        if not query.exists():
            atomic_write_json(query, expression)
        return self.shell_json(project, label, ["facts", "eval", "@" + str(query)])

    def state_summary(self, project: Path, label: str) -> Any:
        return self.command_json(
            label,
            [str(self.args.gentle_cli.resolve()), "--project", str(project), "state-summary"],
            timeout=120.0,
        )

    def report(self, project: Path, report_id: str, label: str) -> Any:
        payload = self.shell_json(project, label, ["primers", "show-report", report_id])
        report = payload.get("report") if isinstance(payload, dict) else None
        if not isinstance(report, dict):
            self.fail("product_failure", f"Report '{report_id}' did not yield a report object")
        return report

    def verify_nonvisual(
        self,
        project: Path,
        verifier: dict[str, Any],
        label: str,
        artifact_root: Path,
    ) -> dict[str, Any]:
        kind = verifier["kind"]
        if kind == "facts":
            result = self.fact_evaluation(project, verifier["expression"], label)
            if result.get("truth") != "satisfied":
                self.fail("product_failure", f"Typed fact verifier is {result.get('truth')!r}")
            return {"kind": kind, "truth": result.get("truth")}
        if kind == "expected_effects":
            tokens = ["introspect", "verify-effects", verifier["capability_id"]]
            for key, value in sorted(verifier.get("args", {}).items()):
                tokens.extend(["--arg", f"{key}={value}"])
            result = self.shell_json(project, label, tokens)
            if result.get("verified") is not True:
                self.fail(
                    "product_failure",
                    f"Expected effects for '{verifier['capability_id']}' are not verified",
                )
            return {"kind": kind, "status": result.get("status")}
        if kind == "report":
            if verifier["schema"] != "gentle.primer_design_report.v1":
                self.fail(
                    "harness_gap",
                    "No typed CLI report adapter is registered for "
                    f"schema {verifier['schema']!r}",
                )
            report = self.report(project, verifier["report_id"], label)
            checks = assert_report(report, verifier)
            return {
                "kind": kind,
                "report_id": verifier["report_id"],
                "schema": report.get("schema"),
                "checks": checks,
            }
        if kind == "artifact":
            relative_path = Path(verifier["path"])
            if relative_path.is_absolute() or ".." in relative_path.parts:
                self.fail("contract_error", f"Unsafe artifact path '{relative_path}'")
            artifact = artifact_root / relative_path
            if not artifact.is_file():
                self.fail("product_failure", f"Typed artifact does not exist: '{relative_path}'")
            digest = sha256_file(artifact)
            if verifier.get("sha256") and digest.lower() != verifier["sha256"].lower():
                self.fail("product_failure", f"Artifact '{relative_path}' digest differs")
            structured = None
            if verifier.get("schema") or verifier.get("required_attributes"):
                structured = load_json(artifact)
                if not isinstance(structured, dict):
                    self.fail(
                        "product_failure",
                        f"Typed artifact '{relative_path}' is not a JSON object",
                    )
            if verifier.get("schema") and structured.get("schema") != verifier["schema"]:
                self.fail("product_failure", f"Artifact '{relative_path}' schema differs")
            for path in verifier.get("required_attributes", []):
                try:
                    json_path(structured, path)
                except KeyError as error:
                    self.fail("product_failure", f"Artifact '{relative_path}' lacks '{path}'")
            return {"kind": kind, "path": relative_path.as_posix(), "sha256": digest}
        if kind == "state":
            summary = self.state_summary(project, label)
            actual = {row["id"] for row in summary.get("sequences", [])}
            missing = sorted(set(verifier["seq_ids"]) - actual)
            if missing:
                self.fail("product_failure", f"Project state is missing sequence ids {missing}")
            return {"kind": kind, "sequence_ids": sorted(actual)}
        if kind == "visible_claim":
            self.fail("contract_error", "Visible claims require a live semantic snapshot")
        self.fail("contract_error", f"Unsupported verifier kind '{kind}'")

    def preflight_projects(
        self, contract: dict[str, Any], starter: Path, oracle: Path
    ) -> None:
        starter_completion = self.fact_evaluation(
            starter, contract["completion_condition"], "starter-completion"
        )
        if starter_completion.get("truth") != "unsatisfied":
            self.fail(
                "contract_error",
                "Starter completion condition must be unsatisfied, not "
                f"{starter_completion.get('truth')!r}",
            )
        oracle_completion = self.fact_evaluation(
            oracle, contract["completion_condition"], "oracle-completion"
        )
        if oracle_completion.get("truth") != "satisfied":
            self.fail(
                "contract_error",
                f"Oracle completion condition is {oracle_completion.get('truth')!r}",
            )
        starter_ids = {
            row["id"]
            for row in self.state_summary(starter, "starter-state").get("sequences", [])
        }
        oracle_ids = {
            row["id"]
            for row in self.state_summary(oracle, "oracle-state").get("sequences", [])
        }
        for source, target in contract["oracle"].get("seq_id_map", {}).items():
            if source not in starter_ids or target not in oracle_ids:
                self.fail(
                    "contract_error",
                    f"Oracle seq_id_map entry {source!r}->{target!r} is not present",
                )
        oracle_checks: list[dict[str, Any]] = []
        for step in contract["steps"]:
            if not step.get("scientific_effect"):
                continue
            after = self.fact_evaluation(
                oracle, step["after"], f"oracle-after-{step['id']}"
            )
            if after.get("truth") != "satisfied":
                self.fail(
                    "contract_error",
                    f"Oracle postcondition for '{step['id']}' is {after.get('truth')!r}",
                )
            for index, verifier in enumerate(step.get("verifiers", []), 1):
                if verifier["kind"] == "visible_claim":
                    continue
                oracle_checks.append(
                    self.verify_nonvisual(
                        oracle,
                        verifier,
                        f"oracle-{step['id']}-{index}",
                        self.output_dir / "oracle",
                    )
                )
        self.ledger["oracle_preflight"] = {
            "starter_completion": "unsatisfied",
            "oracle_completion": "satisfied",
            "typed_checks": oracle_checks,
        }
        self.write_ledger()

    def isolated_gui_environment(self) -> dict[str, str]:
        env = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith("GENTLE_")
            and not key.endswith("_API_KEY")
            and not key.endswith("_TOKEN")
        }
        profile = self.runtime_dir / "profile"
        env.update(
            {
                "HOME": str(profile / "home"),
                "XDG_CONFIG_HOME": str(profile / "config"),
                "XDG_CACHE_HOME": str(profile / "cache"),
                "XDG_DATA_HOME": str(profile / "data"),
                "GENTLE_GUI_TEST_SNAPSHOT": str(self.snapshot_path),
            }
        )
        for path in (env["HOME"], env["XDG_CONFIG_HOME"], env["XDG_CACHE_HOME"], env["XDG_DATA_HOME"]):
            Path(path).mkdir(parents=True, exist_ok=True)
        return env

    def discover_screenshot_tool(self) -> Optional[str]:
        for candidate in ("scrot", "gnome-screenshot", "import"):
            if shutil.which(candidate):
                return candidate
        return None

    def require_x11(self, contract: dict[str, Any]) -> None:
        if sys.platform != "linux":
            self.fail("harness_gap", "The tutorial GUI executor currently requires Linux/X11")
        if not os.environ.get("DISPLAY"):
            self.fail("harness_gap", "DISPLAY is not set; run under xvfb-run or an X11 session")
        if not shutil.which("xdotool"):
            self.fail("harness_gap", "xdotool is required for ordinary X11 input delivery")
        needs_screenshot = any(
            step["evidence"]["screenshot"] == "required" for step in contract["steps"]
        )
        self.screenshot_tool = self.discover_screenshot_tool()
        if needs_screenshot and self.screenshot_tool is None:
            self.fail(
                "harness_gap",
                "The contract requires screenshots, but scrot, gnome-screenshot, and ImageMagick import are unavailable",
            )
        self.ledger["environment"].update(
            {
                "display": os.environ["DISPLAY"],
                "xdotool": self.tool_version(["xdotool", "version"]),
                "screenshot_tool": self.screenshot_tool,
            }
        )
        self.write_ledger()

    def start_gui(self, project: Path) -> None:
        if self.snapshot_path.exists():
            self.snapshot_path.unlink()
        stdout_path = self.logs_dir / "gui.stdout.txt"
        stderr_path = self.logs_dir / "gui.stderr.txt"
        stdout = stdout_path.open("w", encoding="utf-8")
        stderr = stderr_path.open("w", encoding="utf-8")
        try:
            self.gui_process = subprocess.Popen(
                [str(self.args.gentle_bin.resolve()), "--project", str(project)],
                cwd=str(self.runtime_dir),
                env=self.isolated_gui_environment(),
                text=True,
                stdout=stdout,
                stderr=stderr,
            )
        except OSError as error:
            stdout.close()
            stderr.close()
            self.fail("harness_gap", f"Could not launch GENtle GUI: {error}")
        stdout.close()
        stderr.close()
        self.ledger["gui"] = {
            "program": self.args.gentle_bin.name,
            "stdout": self.relative(stdout_path),
            "stderr": self.relative(stderr_path),
        }
        self.write_ledger()

    def stop_gui(self) -> None:
        if self.gui_process is None or self.gui_process.poll() is not None:
            return
        self.gui_process.terminate()
        try:
            self.gui_process.wait(timeout=10.0)
        except subprocess.TimeoutExpired:
            self.gui_process.kill()
            self.gui_process.wait(timeout=5.0)

    def read_snapshot(self) -> Optional[dict[str, Any]]:
        if self.gui_process is not None and self.gui_process.poll() is not None:
            self.fail(
                "product_failure",
                f"GENtle exited before acceptance completed (exit {self.gui_process.returncode})",
            )
        if not self.snapshot_path.is_file():
            return None
        try:
            snapshot = json.loads(self.snapshot_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if snapshot.get("schema") != SNAPSHOT_SCHEMA:
            self.fail(
                "product_failure",
                f"GUI emitted unsupported semantic snapshot schema {snapshot.get('schema')!r}",
            )
        return snapshot

    def wait_snapshot(
        self,
        description: str,
        timeout: float,
        predicate: Callable[[dict[str, Any]], bool],
    ) -> dict[str, Any]:
        deadline = time.monotonic() + timeout * self.args.timeout_scale
        last_snapshot: Optional[dict[str, Any]] = None
        while time.monotonic() < deadline:
            snapshot = self.read_snapshot()
            if snapshot is not None:
                last_snapshot = snapshot
                if snapshot.get("settled") and predicate(snapshot):
                    return snapshot
            time.sleep(0.1)
        generation = None if last_snapshot is None else last_snapshot.get("generation")
        self.fail(
            "product_failure",
            f"Timed out waiting for {description}; last semantic generation={generation}",
        )

    @staticmethod
    def matching_items(
        snapshot: dict[str, Any],
        semantic_id: str,
        scope: Optional[str] = None,
        window: Optional[str] = None,
    ) -> list[dict[str, Any]]:
        return [
            item
            for item in snapshot.get("items", [])
            if item.get("semantic_id") == semantic_id
            and (scope is None or item.get("subject_scope") == scope)
            and (window is None or item.get("window_id") == window)
        ]

    def unique_item(
        self,
        snapshot: dict[str, Any],
        semantic_id: str,
        scope: Optional[str] = None,
        window: Optional[str] = None,
    ) -> dict[str, Any]:
        items = self.matching_items(snapshot, semantic_id, scope, window)
        if len(items) != 1:
            self.fail(
                "product_failure",
                f"Expected one semantic item '{semantic_id}' (scope={scope}, window={window}), found {len(items)}",
            )
        return items[0]

    @staticmethod
    def item_ready(item: dict[str, Any], enabled: bool = True) -> bool:
        state = item.get("state", {})
        return state.get("visible") is True and (not enabled or state.get("enabled") is True)

    def wait_item(
        self,
        semantic_id: str,
        scope: Optional[str],
        window: Optional[str],
        timeout: float,
        *,
        enabled: bool = True,
        min_generation: int = 0,
        allow_scroll: bool = False,
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        deadline = time.monotonic() + timeout * self.args.timeout_scale
        last_snapshot: Optional[dict[str, Any]] = None
        scroll_attempts = 0
        while time.monotonic() < deadline:
            snapshot = self.read_snapshot()
            if snapshot is not None:
                last_snapshot = snapshot
                generation = int(snapshot.get("generation", 0))
                items = self.matching_items(snapshot, semantic_id, scope, window)
                if (
                    snapshot.get("settled")
                    and generation >= min_generation
                    and len(items) == 1
                ):
                    if self.item_ready(items[0], enabled):
                        return snapshot, items[0]
                    if (
                        allow_scroll
                        and items[0].get("state", {}).get("visible") is False
                        and scroll_attempts < 16
                        and self.scroll_owning_window(
                            snapshot,
                            items[0],
                            f"reveal-{semantic_id}-{scroll_attempts + 1}",
                        )
                    ):
                        scroll_attempts += 1
            time.sleep(0.15)
        generation = None if last_snapshot is None else last_snapshot.get("generation")
        self.fail(
            "product_failure",
            f"Timed out waiting for {semantic_id}; last semantic generation={generation}; "
            f"scroll attempts={scroll_attempts}",
        )

    def xdotool(self, label: str, args: list[str]) -> None:
        self.command(
            label,
            ["xdotool", *args],
            cwd=self.runtime_dir,
            timeout=20.0,
            failure_result="harness_gap",
        )

    def point_for_item(self, item: dict[str, Any]) -> tuple[int, int]:
        rect = item["rect_logical_points"]
        scale = float(item["pixels_per_point"])
        x = round((float(rect["min_x"]) + float(rect["max_x"])) * 0.5 * scale)
        y = round((float(rect["min_y"]) + float(rect["max_y"])) * 0.5 * scale)
        return x, y

    def scroll_owning_window(
        self, snapshot: dict[str, Any], item: dict[str, Any], label: str
    ) -> bool:
        window_id = item.get("window_id")
        if window_id != "window.pcr_design":
            return False
        windows = self.matching_items(snapshot, window_id, None, window_id)
        if len(windows) != 1 or not self.item_ready(windows[0], enabled=False):
            return False
        rect = windows[0]["rect_logical_points"]
        scale = float(windows[0]["pixels_per_point"])
        min_x = float(rect["min_x"])
        max_x = float(rect["max_x"])
        min_y = float(rect["min_y"])
        max_y = float(rect["max_y"])
        # The first supported contract uses the right-hand PCR constraint
        # scroll area. Ordinary wheel input keeps clipping authoritative while
        # allowing the harness to reveal a catalogued control before using it.
        x = round((min_x + (max_x - min_x) * 0.75) * scale)
        y = round((min_y + (max_y - min_y) * 0.60) * scale)
        self.xdotool(
            label,
            [
                "mousemove",
                "--sync",
                str(x),
                str(y),
                "click",
                "--repeat",
                "5",
                "--delay",
                "80",
                "5",
            ],
        )
        return True

    def apply_interaction(self, step: dict[str, Any], item: dict[str, Any]) -> None:
        interaction = step["interaction"]
        kind = interaction["kind"]
        x, y = self.point_for_item(item)
        self.xdotool(f"move-{step['id']}", ["mousemove", "--sync", str(x), str(y)])
        if kind == "click" or kind == "select_tab":
            self.xdotool(f"click-{step['id']}", ["click", "1"])
        elif kind == "right_click":
            self.xdotool(f"right-click-{step['id']}", ["click", "3"])
        elif kind == "double_click":
            self.xdotool(
                f"double-click-{step['id']}", ["click", "--repeat", "2", "--delay", "100", "1"]
            )
        elif kind == "replace_text":
            self.xdotool(f"focus-{step['id']}", ["click", "1"])
            self.xdotool(f"select-{step['id']}", ["key", "--clearmodifiers", "ctrl+a"])
            self.xdotool(
                f"type-{step['id']}",
                ["type", "--clearmodifiers", "--delay", "1", interaction["text"]],
            )
        elif kind == "set_checkbox":
            if item.get("state", {}).get("selected") != interaction["selected"]:
                self.xdotool(f"checkbox-{step['id']}", ["click", "1"])
        elif kind == "press_key":
            key = interaction["key"]
            if not re.fullmatch(r"[A-Za-z0-9_+-]+", key):
                self.fail("contract_error", f"Unsupported X11 key token {key!r}")
            self.xdotool(f"key-{step['id']}", ["key", "--clearmodifiers", key])
        else:
            self.fail("contract_error", f"Unsupported interaction kind '{kind}'")

    def visible_claim_matches(
        self, snapshot: dict[str, Any], verifier: dict[str, Any], scope: Optional[str]
    ) -> bool:
        items = self.visible_claim_items(snapshot, verifier, scope)
        if len(items) != 1:
            return False
        item = items[0]
        state = item.get("state", {})
        for key in ("visible", "enabled", "selected"):
            if key in verifier and state.get(key) != verifier[key]:
                return False
        if "outcome_role" in verifier and item.get("outcome_role") != verifier["outcome_role"]:
            return False
        return True

    def visible_claim_items(
        self, snapshot: dict[str, Any], verifier: dict[str, Any], scope: Optional[str]
    ) -> list[dict[str, Any]]:
        all_items = self.matching_items(snapshot, verifier["semantic_id"])
        items = [item for item in all_items if item.get("subject_scope") == scope]
        if not items and scope is not None:
            items = [item for item in all_items if item.get("subject_scope") is None]
        return items

    def first_hidden_visible_claim(
        self,
        snapshot: dict[str, Any],
        claims: list[dict[str, Any]],
        scope: Optional[str],
    ) -> Optional[dict[str, Any]]:
        for claim in claims:
            if claim.get("visible") is not True:
                continue
            items = self.visible_claim_items(snapshot, claim, scope)
            if len(items) == 1 and items[0].get("state", {}).get("visible") is False:
                return items[0]
        return None

    def wait_visible_claims(
        self,
        step: dict[str, Any],
        scope: Optional[str],
        min_generation: int,
        timeout: float,
    ) -> dict[str, Any]:
        claims = [v for v in step.get("verifiers", []) if v["kind"] == "visible_claim"]
        deadline = time.monotonic() + timeout * self.args.timeout_scale
        last_snapshot: Optional[dict[str, Any]] = None
        scroll_attempts = 0
        while time.monotonic() < deadline:
            snapshot = self.read_snapshot()
            if snapshot is not None:
                last_snapshot = snapshot
                generation = int(snapshot.get("generation", 0))
                if snapshot.get("settled") and generation >= min_generation:
                    if all(
                        self.visible_claim_matches(snapshot, claim, scope)
                        for claim in claims
                    ):
                        return snapshot
                    if scroll_attempts < 16:
                        hidden_item = self.first_hidden_visible_claim(
                            snapshot, claims, scope
                        )
                        if hidden_item is not None and self.scroll_owning_window(
                            snapshot,
                            hidden_item,
                            f"reveal-{step['id']}-claim-{scroll_attempts + 1}",
                        ):
                            scroll_attempts += 1
            time.sleep(0.15)
        generation = None if last_snapshot is None else last_snapshot.get("generation")
        self.fail(
            "product_failure",
            f"Timed out waiting for step '{step['id']}' visible claims; "
            f"last semantic generation={generation}; scroll attempts={scroll_attempts}",
        )

    def wait_save_state(self, role: str, timeout: float, min_generation: int = 0) -> dict[str, Any]:
        def predicate(snapshot: dict[str, Any]) -> bool:
            if snapshot.get("generation", 0) < min_generation:
                return False
            items = self.matching_items(snapshot, "main.project.save_state", None, "window.main")
            return len(items) == 1 and items[0].get("outcome_role") == role

        return self.wait_snapshot(f"project save state '{role}'", timeout, predicate)

    def capture_evidence(
        self, index: int, step: dict[str, Any], snapshot: dict[str, Any]
    ) -> dict[str, Any]:
        stem = f"{index:02d}-{step['id']}"
        snapshot_output = self.checkpoints_dir / f"{stem}.snapshot.json"
        atomic_write_json(snapshot_output, snapshot)
        evidence: dict[str, Any] = {
            "snapshot": self.relative(snapshot_output),
            "snapshot_sha256": sha256_file(snapshot_output),
            "screenshot": None,
            "screenshot_sha256": None,
        }
        screenshot_policy = step["evidence"]["screenshot"]
        should_capture = screenshot_policy == "required" or (
            screenshot_policy == "optional" and self.args.capture_optional_screenshots
        )
        if should_capture:
            if self.screenshot_tool is None:
                if screenshot_policy == "required":
                    self.fail("harness_gap", f"Step '{step['id']}' requires a screenshot")
            else:
                screenshot = self.checkpoints_dir / f"{stem}.png"
                if self.screenshot_tool == "scrot":
                    argv = ["scrot", str(screenshot)]
                elif self.screenshot_tool == "gnome-screenshot":
                    argv = ["gnome-screenshot", "-f", str(screenshot)]
                else:
                    argv = ["import", "-window", "root", str(screenshot)]
                self.command(
                    f"screenshot-{step['id']}",
                    argv,
                    cwd=self.runtime_dir,
                    timeout=30.0,
                    failure_result="harness_gap",
                )
                self.require_file(screenshot, "captured screenshot")
                evidence["screenshot"] = self.relative(screenshot)
                evidence["screenshot_sha256"] = sha256_file(screenshot)
        return evidence

    def bootstrap_sequence_window(self, contract: dict[str, Any]) -> None:
        first_subject = contract["steps"][0].get("subject", {})
        scope = subject_scope(first_subject)
        snapshot, item = self.wait_item(
            "main.project.sequence.open", scope, "window.main", 90.0
        )
        x, y = self.point_for_item(item)
        self.xdotool("bootstrap-sequence-move", ["mousemove", "--sync", str(x), str(y)])
        self.xdotool(
            "bootstrap-sequence-open", ["click", "--repeat", "2", "--delay", "100", "1"]
        )
        self.wait_item(
            contract["steps"][0]["target"],
            scope,
            contract["steps"][0]["window"],
            90.0,
            min_generation=int(snapshot.get("generation", 0)) + 1,
        )

    def execute_steps(self, contract: dict[str, Any], project: Path) -> None:
        self.bootstrap_sequence_window(contract)
        for index, step in enumerate(contract["steps"], 1):
            timeout = TIMEOUT_SECONDS.get(step["timeout_class"])
            if timeout is None:
                self.fail(
                    "contract_error", f"Unknown timeout class {step['timeout_class']!r}"
                )
            scope = subject_scope(step.get("subject", {}))
            if step.get("scientific_effect"):
                before = self.fact_evaluation(
                    project, step["before"], f"before-{step['id']}"
                )
                accepted_truths = {"unsatisfied"}
                if step.get("allow_preexisting"):
                    accepted_truths.add("satisfied")
                if before.get("truth") not in accepted_truths:
                    self.fail(
                        "product_failure",
                        f"Before fact for '{step['id']}' is {before.get('truth')!r}",
                    )
            snapshot, item = self.wait_item(
                step["target"], scope, step["window"], timeout, allow_scroll=True
            )
            before_generation = int(snapshot.get("generation", 0))
            self.apply_interaction(step, item)
            snapshot = self.wait_visible_claims(
                step, scope, before_generation + 1, timeout
            )
            if step.get("persists_project_state"):
                snapshot = self.wait_save_state(
                    "unsaved", timeout, min_generation=before_generation + 1
                )
                save_generation = int(snapshot.get("generation", 0))
                self.xdotool(f"save-{step['id']}", ["key", "--clearmodifiers", "ctrl+s"])
                snapshot = self.wait_save_state(
                    "saved", timeout, min_generation=save_generation + 1
                )
            typed_checks: list[dict[str, Any]] = []
            if step.get("scientific_effect"):
                after = self.fact_evaluation(project, step["after"], f"after-{step['id']}")
                if after.get("truth") != "satisfied":
                    self.fail(
                        "product_failure",
                        f"After fact for '{step['id']}' is {after.get('truth')!r}",
                    )
            for verifier_index, verifier in enumerate(step.get("verifiers", []), 1):
                if verifier["kind"] == "visible_claim":
                    typed_checks.append(
                        {
                            "kind": "visible_claim",
                            "semantic_id": verifier["semantic_id"],
                            "matched": True,
                        }
                    )
                else:
                    typed_checks.append(
                        self.verify_nonvisual(
                            project,
                            verifier,
                            f"verify-{step['id']}-{verifier_index}",
                            self.runtime_dir,
                        )
                    )
            evidence = self.capture_evidence(index, step, snapshot)
            checkpoint = {
                "step_id": step["id"],
                "prose_step": step["prose_step"],
                "interaction": step["interaction"]["kind"],
                "semantic_target": step["target"],
                "subject_scope": scope,
                "generation": snapshot.get("generation"),
                "project_sha256": sha256_file(project),
                "typed_checks": typed_checks,
                "evidence": evidence,
                "result": "passed",
            }
            self.ledger["checkpoints"].append(checkpoint)
            self.write_ledger()
        completion = self.fact_evaluation(
            project, contract["completion_condition"], "runtime-completion"
        )
        if completion.get("truth") != "satisfied":
            self.fail(
                "product_failure",
                f"Final completion condition is {completion.get('truth')!r}",
            )
        self.ledger["runtime_completion"] = completion

    def run(self) -> None:
        self.prepare_directories()
        self.write_ledger()
        self.establish_provenance()
        self.validate_manifest()
        _chapter, contract = self.load_contract()
        starter = self.materialize_project("starter", contract["starter"]["example_id"])
        oracle = self.materialize_project("oracle", contract["oracle"]["example_id"])
        self.preflight_projects(contract, starter, oracle)
        if self.args.validate_only:
            self.ledger["result"] = "validated_only"
            return
        self.require_x11(contract)
        self.start_gui(starter)
        self.execute_steps(contract, starter)
        self.ledger["result"] = "passed"


def self_test() -> None:
    expected = "subject-f33b2c9685525d19be7ecd5efd1f518f"
    actual = subject_scope({"sequence": "tp73_locus"})
    assert actual == expected, (actual, expected)
    report = {
        "schema": "demo.v1",
        "count": 2,
        "items": [{"start": 3}, {"start": 8}],
        "enabled": True,
    }
    checks = assert_report(
        report,
        {
            "report_id": "demo",
            "schema": "demo.v1",
            "required_fields": ["items.0.start"],
            "assertions": [
                {"kind": "value", "path": "enabled", "equals": True},
                {
                    "kind": "value",
                    "path": "count",
                    "compare": {"op": "gt", "value": 0},
                },
                {
                    "kind": "relation",
                    "left_path": "items.0.start",
                    "op": "lt",
                    "right_path": "items.1.start",
                },
            ],
        },
    )
    assert len(checks) == 4
    print("tutorial_gui_acceptance self-test: ok")


def parser() -> argparse.ArgumentParser:
    repo_default = Path(__file__).resolve().parents[1]
    result = argparse.ArgumentParser(
        description="Execute one typed GENtle tutorial GUI acceptance contract on Linux/X11."
    )
    result.add_argument("--tutorial-id", default="simple_pcr_selection_gui")
    result.add_argument("--repo-root", type=Path, default=repo_default)
    result.add_argument("--manifest", type=Path)
    result.add_argument("--tutorial-sources", type=Path)
    result.add_argument("--catalog-meta", type=Path)
    result.add_argument("--examples", type=Path)
    result.add_argument("--gentle-bin", type=Path)
    result.add_argument("--gentle-cli", type=Path)
    result.add_argument("--examples-bin", type=Path)
    result.add_argument("--output-dir", type=Path)
    result.add_argument("--timeout-scale", type=float, default=1.0)
    result.add_argument("--capture-optional-screenshots", action="store_true")
    result.add_argument("--validate-only", action="store_true")
    result.add_argument("--self-test", action="store_true")
    return result


def resolve_defaults(args: argparse.Namespace) -> None:
    repo = args.repo_root.resolve()
    args.repo_root = repo
    args.manifest = (args.manifest or repo / "docs/tutorial/manifest.json").resolve()
    args.tutorial_sources = (
        args.tutorial_sources or repo / "docs/tutorial/sources"
    ).resolve()
    args.catalog_meta = (
        args.catalog_meta or repo / "docs/tutorial/sources/catalog_meta.json"
    ).resolve()
    args.examples = (args.examples or repo / "docs/examples/workflows").resolve()
    args.gentle_bin = (args.gentle_bin or repo / "target/debug/gentle").resolve()
    args.gentle_cli = (args.gentle_cli or repo / "target/debug/gentle_cli").resolve()
    args.examples_bin = (
        args.examples_bin or repo / "target/debug/gentle_examples_docs"
    ).resolve()
    if args.output_dir is None:
        raise AuditFailure("harness_gap", "--output-dir is required for retained audit evidence")
    args.output_dir = args.output_dir.resolve()
    if args.timeout_scale <= 0:
        raise AuditFailure("harness_gap", "--timeout-scale must be greater than zero")


def main() -> int:
    args = parser().parse_args()
    if args.self_test:
        self_test()
        return 0
    runner: Optional[TutorialGuiRunner] = None
    try:
        resolve_defaults(args)
        runner = TutorialGuiRunner(args)
        runner.run()
        return 0
    except AuditFailure as error:
        if runner is not None:
            runner.ledger["result"] = error.result
            runner.ledger["failure"] = {"kind": error.result, "message": str(error)}
        print(f"{error.result}: {error}", file=sys.stderr)
        return {"product_failure": 1, "contract_error": 2, "harness_gap": 3}.get(
            error.result, 1
        )
    except KeyboardInterrupt:
        if runner is not None:
            runner.ledger["result"] = "interrupted"
            runner.ledger["failure"] = {"kind": "interrupted", "message": "Interrupted"}
        print("interrupted", file=sys.stderr)
        return 130
    finally:
        if runner is not None:
            runner.stop_gui()
            runner.ledger["finished_at"] = utc_now()
            runner.write_ledger()


def interrupt_on_signal(_signum: int, _frame: Any) -> None:
    raise KeyboardInterrupt


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, interrupt_on_signal)
    raise SystemExit(main())
