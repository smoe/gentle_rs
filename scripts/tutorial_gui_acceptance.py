#!/usr/bin/env python3
"""Run typed GENtle tutorial GUI acceptance contracts through ordinary X11 input.

The tutorial manifest supplies intent and typed postconditions. This runner owns
input delivery, independent verification, and the final pass/fail verdict. It
never executes prose or command strings from a tutorial.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable


SNAPSHOT_SCHEMA = "gentle.gui_semantic_snapshot.v2"
CONTRACT_SCHEMA = "gentle.tutorial_gui_acceptance.v1"
LEDGER_SCHEMA = "gentle.tutorial_gui_acceptance_ledger.v1"
RUN_SCHEMA = "gentle.tutorial_gui_acceptance_run.v1"
ENVIRONMENT_SCHEMA = "gentle.tutorial_acceptance_environment.v1"
PREPARATION_SCHEMA = "gentle.tutorial_gui_project_preparation.v1"

TIMEOUT_DEFAULTS = {
    "instant": 5.0,
    "interactive": 25.0,
    "io": 75.0,
    "compute": 300.0,
}

SAFE_KEY_NAMES = {
    "enter": "Return",
    "return": "Return",
    "escape": "Escape",
    "tab": "Tab",
    "space": "space",
    "up": "Up",
    "down": "Down",
    "left": "Left",
    "right": "Right",
}


class AcceptanceFailure(RuntimeError):
    def __init__(self, failure_class: str, message: str):
        super().__init__(message)
        self.failure_class = failure_class


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AcceptanceFailure(
            "tutorial_ambiguity", f"Could not read JSON '{path}': {error}"
        ) from error


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


def numeric_compare(actual: float, op: str, expected: float) -> bool:
    if op in {"eq", "=", "=="}:
        return actual == expected
    if op in {"ne", "!="}:
        return actual != expected
    if op in {"gt", ">"}:
        return actual > expected
    if op in {"gte", ">="}:
        return actual >= expected
    if op in {"lt", "<"}:
        return actual < expected
    if op in {"lte", "<="}:
        return actual <= expected
    raise AcceptanceFailure("tutorial_ambiguity", f"Unsupported comparison '{op}'")


def assert_report_contract(report: dict[str, Any], verifier: dict[str, Any]) -> None:
    if report.get("schema") != verifier["schema"]:
        raise AcceptanceFailure(
            "product_failure",
            f"Report schema is {report.get('schema')!r}, expected {verifier['schema']!r}",
        )
    for path in verifier.get("required_fields", []):
        try:
            json_path(report, path)
        except KeyError as error:
            raise AcceptanceFailure(
                "product_failure", f"Required report field '{path}' is absent"
            ) from error
    for assertion in verifier.get("assertions", []):
        kind = assertion.get("kind")
        if kind == "value":
            path = assertion["path"]
            try:
                actual = json_path(report, path)
            except KeyError as error:
                raise AcceptanceFailure(
                    "product_failure", f"Report assertion field '{path}' is absent"
                ) from error
            if "equals" in assertion and actual != assertion["equals"]:
                raise AcceptanceFailure(
                    "product_failure",
                    f"Report field '{path}' is {actual!r}, expected {assertion['equals']!r}",
                )
            comparison = assertion.get("compare")
            if comparison is not None:
                if not isinstance(actual, (int, float)) or isinstance(actual, bool):
                    raise AcceptanceFailure(
                        "product_failure", f"Report field '{path}' is not numeric"
                    )
                if not numeric_compare(
                    float(actual), comparison["op"], float(comparison["value"])
                ):
                    raise AcceptanceFailure(
                        "product_failure",
                        f"Report comparison failed: {path}={actual!r} "
                        f"{comparison['op']} {comparison['value']!r}",
                    )
            if assertion.get("non_empty") and not value_is_non_empty(actual):
                raise AcceptanceFailure(
                    "product_failure", f"Report field '{path}' is empty"
                )
        elif kind == "relation":
            try:
                left = json_path(report, assertion["left_path"])
                right = json_path(report, assertion["right_path"])
            except KeyError as error:
                raise AcceptanceFailure(
                    "product_failure", f"Report relation field '{error.args[0]}' is absent"
                ) from error
            if (
                not isinstance(left, (int, float))
                or isinstance(left, bool)
                or not isinstance(right, (int, float))
                or isinstance(right, bool)
            ):
                raise AcceptanceFailure(
                    "product_failure", "Report relation operands must be numeric"
                )
            if not numeric_compare(float(left), assertion["op"], float(right)):
                raise AcceptanceFailure(
                    "product_failure",
                    f"Report relation failed: {assertion['left_path']}={left!r} "
                    f"{assertion['op']} {assertion['right_path']}={right!r}",
                )
        else:
            raise AcceptanceFailure(
                "tutorial_ambiguity", f"Unsupported report assertion kind '{kind}'"
            )


def value_is_non_empty(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, (str, list, dict)):
        return len(value) > 0
    return True


def sanitized_label(value: str) -> str:
    label = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in value)
    return label[:80] or "command"


def fixed_shell_command(argv: list[str]) -> str:
    """Quote runner-owned shared-shell arguments without accepting command text."""
    return shlex.join(argv)


def selected_chapters(
    manifest: dict[str, Any], chapter_ids: list[str], profile: str | None
) -> list[dict[str, Any]]:
    chapters = manifest.get("chapters", [])
    if chapter_ids:
        by_id = {chapter.get("id"): chapter for chapter in chapters}
        missing = [chapter_id for chapter_id in chapter_ids if chapter_id not in by_id]
        if missing:
            raise AcceptanceFailure(
                "tutorial_ambiguity",
                f"Unknown tutorial chapter(s): {', '.join(missing)}",
            )
        selected = [by_id[chapter_id] for chapter_id in chapter_ids]
    else:
        selected = [
            chapter
            for chapter in chapters
            if chapter.get("gui_acceptance", {}).get("profile") == profile
        ]
    missing_contract = [
        chapter.get("id", "<unknown>")
        for chapter in selected
        if not isinstance(chapter.get("gui_acceptance"), dict)
    ]
    if missing_contract:
        raise AcceptanceFailure(
            "tutorial_ambiguity",
            f"Selected chapter(s) have no GUI acceptance contract: {', '.join(missing_contract)}",
        )
    if not selected:
        raise AcceptanceFailure(
            "tutorial_ambiguity", f"No GUI acceptance chapters matched profile '{profile}'"
        )
    return selected


def redacted_environment(environment: dict[str, str]) -> tuple[list[dict[str, str]], set[str]]:
    sensitive_names = {
        name
        for name in environment
        if name.startswith("GENTLE_") or name.endswith("_API_KEY")
    }
    rows = [
        {
            "name": name,
            "value_sha256": sha256_bytes(environment[name].encode("utf-8")),
        }
        for name in sorted(sensitive_names)
    ]
    return rows, sensitive_names


@dataclass
class CommandResult:
    payload: Any
    receipt: dict[str, Any]


class CommandRecorder:
    def __init__(self, output_dir: Path, cwd: Path):
        self.output_dir = output_dir
        self.cwd = cwd
        self.receipts: list[dict[str, Any]] = []
        self.counter = 0

    def run_json(
        self,
        argv: list[str],
        label: str,
        *,
        timeout: float = 180.0,
        env: dict[str, str] | None = None,
        failure_class: str = "product_failure",
    ) -> CommandResult:
        self.counter += 1
        stem = f"{self.counter:03d}-{sanitized_label(label)}"
        started = time.monotonic()
        try:
            completed = subprocess.run(
                argv,
                cwd=self.cwd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise AcceptanceFailure(
                failure_class, f"Command '{label}' could not complete: {error}"
            ) from error
        elapsed_ms = round((time.monotonic() - started) * 1000)
        stdout_path = self.output_dir / "commands" / f"{stem}.stdout"
        stderr_path = self.output_dir / "commands" / f"{stem}.stderr"
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        stdout_path.write_bytes(completed.stdout)
        stderr_path.write_bytes(completed.stderr)
        receipt = {
            "label": label,
            "argv": argv,
            "cwd": str(self.cwd),
            "exit_code": completed.returncode,
            "elapsed_ms": elapsed_ms,
            "stdout_path": str(stdout_path),
            "stdout_sha256": sha256_bytes(completed.stdout),
            "stderr_path": str(stderr_path),
            "stderr_sha256": sha256_bytes(completed.stderr),
        }
        self.receipts.append(receipt)
        if completed.returncode != 0:
            diagnostic = completed.stderr.decode("utf-8", errors="replace").strip()
            raise AcceptanceFailure(
                failure_class,
                f"Command '{label}' exited {completed.returncode}: {diagnostic[-1000:]}",
            )
        try:
            payload = json.loads(completed.stdout)
        except json.JSONDecodeError as error:
            raise AcceptanceFailure(
                failure_class, f"Command '{label}' did not emit one JSON document: {error}"
            ) from error
        return CommandResult(payload=payload, receipt=receipt)


class TutorialAcceptanceRun:
    def __init__(
        self,
        args: argparse.Namespace,
        repo_root: Path,
        chapter: dict[str, Any],
        environment_record: dict[str, Any],
    ):
        self.args = args
        self.repo_root = repo_root
        self.chapter = chapter
        self.acceptance = chapter["gui_acceptance"]
        self.chapter_dir = args.evidence_dir / chapter["id"]
        if self.chapter_dir.exists() and any(self.chapter_dir.iterdir()):
            raise AcceptanceFailure(
                "harness_gap", f"Evidence directory is not empty: {self.chapter_dir}"
            )
        self.chapter_dir.mkdir(parents=True, exist_ok=True)
        self.recorder = CommandRecorder(self.chapter_dir, repo_root)
        self.snapshot_path = self.chapter_dir / "live-semantic-snapshot.json"
        self.gui_stdout_path = self.chapter_dir / "gui.stdout"
        self.gui_stderr_path = self.chapter_dir / "gui.stderr"
        self.gui_process: subprocess.Popen[bytes] | None = None
        self.gui_stdout = None
        self.gui_stderr = None
        self.isolation_paths = self.prepare_isolation_paths()
        self.process_environment = self.isolated_process_environment()
        self.starter_preparation: dict[str, Any] = {}
        self.oracle_preparation: dict[str, Any] = {}
        self.sequence_scopes: dict[str, str] = {}
        self.steps: list[dict[str, Any]] = []
        self.ledger: dict[str, Any] = {
            "schema": LEDGER_SCHEMA,
            "chapter_id": chapter["id"],
            "chapter_title": chapter.get("title", ""),
            "profile": self.acceptance.get("profile"),
            "network_policy": self.acceptance.get("network"),
            "network_enforcement": args.network_enforcement,
            "manifest_path": str(args.manifest),
            "manifest_sha256": sha256_file(args.manifest),
            "acceptance_contract_sha256": sha256_bytes(
                canonical_json_bytes(self.acceptance)
            ),
            "environment": environment_record,
            "isolation_paths": {
                name: str(path) for name, path in self.isolation_paths.items()
            },
            "steps": self.steps,
            "status": "running",
            "failure_class": None,
            "message": "",
        }

    @property
    def ledger_path(self) -> Path:
        return self.chapter_dir / "acceptance-ledger.json"

    def write_ledger(self) -> None:
        self.ledger["command_receipts"] = self.recorder.receipts
        atomic_write_json(self.ledger_path, self.ledger)

    def prepare_project(self, phase: str) -> dict[str, Any]:
        project_path = self.chapter_dir / f"{phase}.project.gentle.json"
        run_dir = self.chapter_dir / f"{phase}-workflow"
        result = self.recorder.run_json(
            [
                str(self.args.examples_docs),
                "tutorial-gui-project",
                "--chapter",
                self.chapter["id"],
                "--phase",
                phase,
                "--project-output",
                str(project_path),
                "--run-dir",
                str(run_dir),
                "--source",
                str(self.args.workflow_source),
                "--manifest",
                str(self.args.manifest),
                "--repo-root",
                str(self.repo_root),
            ],
            f"prepare-{phase}",
            env=self.process_environment,
            failure_class="tutorial_ambiguity",
        ).payload
        if result.get("schema") != PREPARATION_SCHEMA:
            raise AcceptanceFailure(
                "tutorial_ambiguity",
                f"Unexpected {phase} preparation schema {result.get('schema')!r}",
            )
        if result.get("project_sha256") != sha256_file(project_path):
            raise AcceptanceFailure(
                "harness_gap", f"{phase} preparation project hash does not match its file"
            )
        return result

    def fact_eval(
        self, project_path: Path, expression: dict[str, Any], label: str
    ) -> dict[str, Any]:
        expression_path = self.chapter_dir / "verifiers" / f"{sanitized_label(label)}.fact.json"
        atomic_write_json(expression_path, expression)
        output = self.recorder.run_json(
            [
                str(self.args.gentle_cli),
                "--project",
                str(project_path),
                "shell",
                fixed_shell_command(["facts", "eval", f"@{expression_path}"]),
            ],
            label,
            env=self.process_environment,
        ).payload
        if output.get("schema") != "gentle.fact_evaluation.v1":
            raise AcceptanceFailure(
                "product_failure", f"Fact verifier '{label}' returned an unexpected schema"
            )
        return output

    def show_report(
        self, project_path: Path, verifier: dict[str, Any], label: str
    ) -> dict[str, Any]:
        schema = verifier["schema"]
        if schema != "gentle.primer_design_report.v1":
            raise AcceptanceFailure(
                "harness_gap",
                f"No fixed tutorial verifier route is registered for report schema '{schema}'",
            )
        output = self.recorder.run_json(
            [
                str(self.args.gentle_cli),
                "--project",
                str(project_path),
                "primers",
                "show-report",
                verifier["report_id"],
            ],
            label,
            env=self.process_environment,
        ).payload
        report = output.get("report")
        if not isinstance(report, dict):
            raise AcceptanceFailure("product_failure", f"Report '{label}' is absent")
        assert_report_contract(report, verifier)
        return {
            "status": "pass",
            "schema": report.get("schema"),
            "report_id": report.get("report_id"),
            "content_sha256": sha256_bytes(canonical_json_bytes(report)),
        }

    def state_verify(
        self, project_path: Path, verifier: dict[str, Any], phase: str, label: str
    ) -> dict[str, Any]:
        output = self.recorder.run_json(
            [str(self.args.gentle_cli), "--project", str(project_path), "state-summary"],
            label,
            env=self.process_environment,
        ).payload
        present = {row.get("id") for row in output.get("sequences", [])}
        mapping = self.binding_map(
            self.starter_preparation if phase == "starter" else self.oracle_preparation
        )
        required = [mapping.get(seq_id, seq_id) for seq_id in verifier.get("seq_ids", [])]
        missing = [seq_id for seq_id in required if seq_id not in present]
        if missing:
            raise AcceptanceFailure(
                "product_failure", f"State verifier is missing sequence(s): {', '.join(missing)}"
            )
        return {"status": "pass", "sequence_ids": required}

    def expected_effects_verify(
        self, project_path: Path, verifier: dict[str, Any], label: str
    ) -> dict[str, Any]:
        argv = [
            str(self.args.gentle_cli),
            "--project",
            str(project_path),
            "introspect",
            "verify-effects",
            verifier["capability_id"],
        ]
        for name, value in sorted(verifier.get("args", {}).items()):
            argv.extend(["--arg", f"{name}={value}"])
        output = self.recorder.run_json(
            argv, label, env=self.process_environment
        ).payload
        if output.get("schema") != "gentle.introspection.v1" or not output.get("verified"):
            raise AcceptanceFailure(
                "product_failure",
                f"Expected effects were not verified for {verifier['capability_id']}: "
                f"{output.get('status')}",
            )
        return {
            "status": "pass",
            "capability_id": output.get("canonical_capability_id"),
            "verification_status": output.get("status"),
        }

    def artifact_verify(
        self, verifier: dict[str, Any], artifact_root: Path
    ) -> dict[str, Any]:
        relative = Path(verifier["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise AcceptanceFailure(
                "tutorial_ambiguity", f"Artifact path is not confined: {relative}"
            )
        path = artifact_root / relative
        if not path.is_file():
            raise AcceptanceFailure("product_failure", f"Artifact is absent: {path}")
        actual_sha = sha256_file(path)
        if verifier.get("sha256") and actual_sha != verifier["sha256"]:
            raise AcceptanceFailure(
                "product_failure", f"Artifact hash differs for '{relative}'"
            )
        content: Any = None
        if verifier.get("schema") or verifier.get("required_attributes"):
            try:
                content = json.loads(path.read_text(encoding="utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError):
                content = path.read_text(encoding="utf-8", errors="replace")
        if verifier.get("schema"):
            if not isinstance(content, dict) or content.get("schema") != verifier["schema"]:
                raise AcceptanceFailure(
                    "product_failure", f"Artifact schema differs for '{relative}'"
                )
        for attribute in verifier.get("required_attributes", []):
            if isinstance(content, dict):
                try:
                    json_path(content, attribute)
                except KeyError as error:
                    raise AcceptanceFailure(
                        "product_failure",
                        f"Artifact '{relative}' lacks attribute '{attribute}'",
                    ) from error
            elif attribute not in str(content):
                raise AcceptanceFailure(
                    "product_failure",
                    f"Artifact '{relative}' lacks marker '{attribute}'",
                )
        return {"status": "pass", "path": str(path), "sha256": actual_sha}

    @staticmethod
    def binding_map(preparation: dict[str, Any]) -> dict[str, str]:
        return {
            row["source_sequence_id"]: row["resolved_sequence_id"]
            for row in preparation.get("sequence_bindings", [])
        }

    @staticmethod
    def scope_map(preparation: dict[str, Any]) -> dict[str, str]:
        return {
            row["source_sequence_id"]: row["subject_scope"]
            for row in preparation.get("sequence_bindings", [])
        }

    def preflight_contract(self) -> None:
        self.starter_preparation = self.prepare_project("starter")
        self.oracle_preparation = self.prepare_project("oracle")
        self.sequence_scopes = self.scope_map(self.starter_preparation)
        starter_path = Path(self.starter_preparation["project_path"])
        oracle_path = Path(self.oracle_preparation["project_path"])
        starter_completion = self.fact_eval(
            starter_path, self.acceptance["completion_condition"], "starter-completion"
        )
        if starter_completion.get("truth") != "unsatisfied":
            failure_class = (
                "starter_precompleted"
                if starter_completion.get("truth") == "satisfied"
                else "tutorial_ambiguity"
            )
            raise AcceptanceFailure(
                failure_class,
                f"Starter completion fact is {starter_completion.get('truth')!r}, expected 'unsatisfied'",
            )
        oracle_completion = self.fact_eval(
            oracle_path, self.acceptance["completion_condition"], "oracle-completion"
        )
        if oracle_completion.get("truth") != "satisfied":
            raise AcceptanceFailure(
                "tutorial_ambiguity",
                f"Oracle completion fact is {oracle_completion.get('truth')!r}, expected 'satisfied'",
            )
        for step in self.acceptance["steps"]:
            if not step.get("scientific_effect"):
                continue
            before = self.fact_eval(
                starter_path, step["before"], f"oracle-check-{step['id']}-before"
            )
            expected_before = {"unsatisfied"}
            if step.get("allow_preexisting"):
                expected_before.add("satisfied")
            if before.get("truth") not in expected_before:
                raise AcceptanceFailure(
                    "tutorial_ambiguity",
                    f"Starter fact for step '{step['id']}' is {before.get('truth')!r}",
                )
            after = self.fact_eval(
                oracle_path, step["after"], f"oracle-check-{step['id']}-after"
            )
            if after.get("truth") != "satisfied":
                raise AcceptanceFailure(
                    "tutorial_ambiguity",
                    f"Oracle fact for step '{step['id']}' is not satisfied",
                )
            for index, verifier in enumerate(step.get("verifiers", []), start=1):
                if verifier["kind"] == "report":
                    self.show_report(
                        oracle_path,
                        verifier,
                        f"oracle-check-{step['id']}-report-{index}",
                    )
                elif verifier["kind"] == "facts":
                    result = self.fact_eval(
                        oracle_path,
                        verifier["expression"],
                        f"oracle-check-{step['id']}-facts-{index}",
                    )
                    if result.get("truth") != "satisfied":
                        raise AcceptanceFailure(
                            "tutorial_ambiguity", "Oracle typed fact verifier is not satisfied"
                        )
        self.ledger["starter"] = self.starter_preparation
        self.ledger["oracle"] = self.oracle_preparation
        self.ledger["preflight"] = {
            "status": "pass",
            "starter_completion_truth": starter_completion.get("truth"),
            "oracle_completion_truth": oracle_completion.get("truth"),
        }
        self.write_ledger()

    def prepare_isolation_paths(self) -> dict[str, Path]:
        profile_root = self.chapter_dir / "profile"
        paths = {
            "home": profile_root / "home",
            "xdg_config": profile_root / "xdg-config",
            "xdg_cache": profile_root / "xdg-cache",
            "xdg_data": profile_root / "xdg-data",
            "tmpdir": profile_root / "tmp",
        }
        for path in paths.values():
            path.mkdir(parents=True, exist_ok=True)
        return paths

    def isolated_process_environment(self) -> dict[str, str]:
        inherited = dict(os.environ)
        _, sensitive_names = redacted_environment(inherited)
        environment = {
            key: value
            for key, value in inherited.items()
            if key not in sensitive_names
        }
        environment.update(
            {
                "HOME": str(self.isolation_paths["home"]),
                "XDG_CONFIG_HOME": str(self.isolation_paths["xdg_config"]),
                "XDG_CACHE_HOME": str(self.isolation_paths["xdg_cache"]),
                "XDG_DATA_HOME": str(self.isolation_paths["xdg_data"]),
                "TMPDIR": str(self.isolation_paths["tmpdir"]),
                "LANG": "C.UTF-8",
                "LC_ALL": "C.UTF-8",
                "TZ": "UTC",
            }
        )
        return environment

    def isolated_gui_environment(self) -> dict[str, str]:
        environment = dict(self.process_environment)
        environment["GENTLE_GUI_TEST_SNAPSHOT"] = str(self.snapshot_path)
        return environment

    def launch_gui(self) -> None:
        self.snapshot_path.unlink(missing_ok=True)
        self.gui_stdout = self.gui_stdout_path.open("wb")
        self.gui_stderr = self.gui_stderr_path.open("wb")
        self.gui_process = subprocess.Popen(
            [
                str(self.args.gentle),
                "--project",
                self.starter_preparation["project_path"],
            ],
            cwd=self.repo_root,
            env=self.isolated_gui_environment(),
            stdout=self.gui_stdout,
            stderr=self.gui_stderr,
            start_new_session=True,
        )
        self.ledger["gui_process"] = {
            "argv": [
                str(self.args.gentle),
                "--project",
                self.starter_preparation["project_path"],
            ],
            "pid": self.gui_process.pid,
            "stdout_path": str(self.gui_stdout_path),
            "stderr_path": str(self.gui_stderr_path),
        }

    def stop_gui(self) -> None:
        if self.gui_process is not None and self.gui_process.poll() is None:
            try:
                os.killpg(self.gui_process.pid, signal.SIGTERM)
                self.gui_process.wait(timeout=10)
            except (ProcessLookupError, subprocess.TimeoutExpired):
                try:
                    os.killpg(self.gui_process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                self.gui_process.wait(timeout=5)
        if self.gui_stdout is not None:
            self.gui_stdout.close()
        if self.gui_stderr is not None:
            self.gui_stderr.close()
        process_info = self.ledger.get("gui_process")
        if isinstance(process_info, dict):
            process_info["exit_code"] = (
                None if self.gui_process is None else self.gui_process.returncode
            )
            if self.gui_stdout_path.is_file():
                process_info["stdout_sha256"] = sha256_file(self.gui_stdout_path)
            if self.gui_stderr_path.is_file():
                process_info["stderr_sha256"] = sha256_file(self.gui_stderr_path)

    def read_snapshot(self) -> dict[str, Any] | None:
        try:
            snapshot = json.loads(self.snapshot_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if snapshot.get("schema") != SNAPSHOT_SCHEMA:
            raise AcceptanceFailure(
                "harness_gap", f"Unexpected GUI snapshot schema {snapshot.get('schema')!r}"
            )
        return snapshot

    def wait_snapshot(
        self,
        predicate: Callable[[dict[str, Any]], bool],
        timeout: float,
        description: str,
        *,
        after_generation: int | None = None,
    ) -> dict[str, Any]:
        deadline = time.monotonic() + timeout
        observed_ids: set[str] = set()
        last_snapshot: dict[str, Any] | None = None
        while time.monotonic() < deadline:
            if self.gui_process is not None and self.gui_process.poll() is not None:
                raise AcceptanceFailure(
                    "product_failure",
                    f"GENtle exited {self.gui_process.returncode} while waiting for {description}",
                )
            snapshot = self.read_snapshot()
            if snapshot is not None:
                last_snapshot = snapshot
                observed_ids.update(
                    item.get("semantic_id", "") for item in snapshot.get("items", [])
                )
                generation_ok = after_generation is None or snapshot.get("generation", 0) > after_generation
                if snapshot.get("settled") and generation_ok and predicate(snapshot):
                    return snapshot
            time.sleep(0.1)
        failure_class = "harness_gap" if not observed_ids else "product_failure"
        generation = None if last_snapshot is None else last_snapshot.get("generation")
        raise AcceptanceFailure(
            failure_class,
            f"Timed out waiting for {description}; last generation={generation}, "
            f"observed semantic ids={sorted(observed_ids)}",
        )

    def item_for(
        self,
        snapshot: dict[str, Any],
        semantic_id: str,
        *,
        window_id: str | None = None,
        subject_scope: str | None = None,
        allow_unscoped_fallback: bool = False,
    ) -> dict[str, Any] | None:
        candidates = [
            item
            for item in snapshot.get("items", [])
            if item.get("semantic_id") == semantic_id
            and (window_id is None or item.get("window_id") == window_id)
        ]
        if subject_scope is not None:
            scoped = [
                item for item in candidates if item.get("subject_scope") == subject_scope
            ]
            if scoped:
                candidates = scoped
            elif allow_unscoped_fallback:
                candidates = [item for item in candidates if item.get("subject_scope") is None]
            else:
                candidates = []
        if len(candidates) > 1:
            raise AcceptanceFailure(
                "harness_gap",
                f"Semantic target '{semantic_id}' is ambiguous ({len(candidates)} matches)",
            )
        return candidates[0] if candidates else None

    def scope_for_step(self, step: dict[str, Any]) -> str | None:
        sequence_id = step.get("subject", {}).get("sequence")
        if sequence_id is None:
            return None
        if sequence_id not in self.sequence_scopes:
            raise AcceptanceFailure(
                "tutorial_ambiguity",
                f"Step '{step['id']}' sequence '{sequence_id}' has no prepared subject scope",
            )
        return self.sequence_scopes[sequence_id]

    def target_ready(
        self, snapshot: dict[str, Any], step: dict[str, Any], scope: str | None
    ) -> bool:
        item = self.item_for(
            snapshot,
            step["target"],
            window_id=step["window"],
            subject_scope=scope,
        )
        return bool(item and item.get("state", {}).get("visible") and item.get("state", {}).get("enabled"))

    def visible_claim_holds(
        self,
        snapshot: dict[str, Any],
        verifier: dict[str, Any],
        scope: str | None,
    ) -> bool:
        semantic_id = verifier["semantic_id"]
        item = self.item_for(
            snapshot,
            semantic_id,
            subject_scope=scope,
            allow_unscoped_fallback=semantic_id.startswith("window."),
        )
        if item is None:
            return False
        state = item.get("state", {})
        for field in ("visible", "enabled", "selected"):
            if field in verifier and state.get(field) != verifier[field]:
                return False
        if "outcome_role" in verifier and item.get("outcome_role") != verifier["outcome_role"]:
            return False
        return True

    def emit_x11(self, step: dict[str, Any], item: dict[str, Any]) -> dict[str, Any]:
        rectangle = item["rect_logical_points"]
        scale = float(item["pixels_per_point"])
        x = round((float(rectangle["min_x"]) + float(rectangle["max_x"])) / 2 * scale)
        y = round((float(rectangle["min_y"]) + float(rectangle["max_y"])) / 2 * scale)
        interaction = step["interaction"]
        kind = interaction["kind"]
        commands: list[list[str]] = [
            [str(self.args.xdotool), "mousemove", "--sync", str(x), str(y)]
        ]
        if kind in {"click", "select_tab"}:
            commands.append([str(self.args.xdotool), "click", "1"])
        elif kind == "right_click":
            commands.append([str(self.args.xdotool), "click", "3"])
        elif kind == "double_click":
            commands.append(
                [
                    str(self.args.xdotool),
                    "click",
                    "--repeat",
                    "2",
                    "--delay",
                    "100",
                    "1",
                ]
            )
        elif kind == "replace_text":
            commands.extend(
                [
                    [str(self.args.xdotool), "click", "1"],
                    [str(self.args.xdotool), "key", "--clearmodifiers", "ctrl+a"],
                    [
                        str(self.args.xdotool),
                        "type",
                        "--clearmodifiers",
                        "--delay",
                        "1",
                        "--",
                        interaction["text"],
                    ],
                ]
            )
        elif kind == "set_checkbox":
            if bool(item.get("state", {}).get("selected")) != bool(interaction["selected"]):
                commands.append([str(self.args.xdotool), "click", "1"])
        elif kind == "press_key":
            key = SAFE_KEY_NAMES.get(str(interaction["key"]).lower())
            if key is None:
                raise AcceptanceFailure(
                    "tutorial_ambiguity", f"Unsupported typed key '{interaction['key']}'"
                )
            commands = [[str(self.args.xdotool), "key", "--clearmodifiers", key]]
        else:
            raise AcceptanceFailure(
                "tutorial_ambiguity", f"Unsupported interaction kind '{kind}'"
            )
        for command in commands:
            completed = subprocess.run(
                command,
                cwd=self.repo_root,
                env=self.process_environment,
                timeout=10,
                check=False,
            )
            if completed.returncode != 0:
                raise AcceptanceFailure(
                    "harness_gap",
                    f"X11 event command exited {completed.returncode}: {command}",
                )
        return {"kind": kind, "screen_x": x, "screen_y": y, "argv": commands}

    def wait_after_interaction(
        self,
        step: dict[str, Any],
        scope: str | None,
        prior_generation: int,
        timeout: float,
    ) -> dict[str, Any]:
        visible_verifiers = [
            verifier
            for verifier in step.get("verifiers", [])
            if verifier["kind"] == "visible_claim"
        ]

        def ready(snapshot: dict[str, Any]) -> bool:
            if not all(
                self.visible_claim_holds(snapshot, verifier, scope)
                for verifier in visible_verifiers
            ):
                return False
            if step.get("persists_project_state"):
                save_state = self.item_for(snapshot, "main.project.save_state")
                if not save_state or save_state.get("outcome_role") != "unsaved":
                    return False
            if step.get("scientific_effect"):
                target = self.item_for(
                    snapshot,
                    step["target"],
                    window_id=step["window"],
                    subject_scope=scope,
                )
                return bool(target and target.get("state", {}).get("enabled"))
            return True

        return self.wait_snapshot(
            ready,
            timeout,
            f"postconditions for step '{step['id']}'",
            after_generation=prior_generation,
        )

    def save_project(self, prior_generation: int, timeout: float) -> dict[str, Any]:
        completed = subprocess.run(
            [str(self.args.xdotool), "key", "--clearmodifiers", "ctrl+s"],
            cwd=self.repo_root,
            env=self.process_environment,
            timeout=10,
            check=False,
        )
        if completed.returncode != 0:
            raise AcceptanceFailure("harness_gap", "Could not emit Ctrl+S through X11")

        def saved(snapshot: dict[str, Any]) -> bool:
            item = self.item_for(snapshot, "main.project.save_state")
            return bool(item and item.get("outcome_role") == "saved")

        snapshot = self.wait_snapshot(
            saved,
            timeout,
            "known-path project save",
            after_generation=prior_generation,
        )
        return {
            "event": "ctrl+s",
            "generation": snapshot["generation"],
            "project_sha256": sha256_file(Path(self.starter_preparation["project_path"])),
        }

    def verify_step(
        self,
        step: dict[str, Any],
        snapshot: dict[str, Any],
        scope: str | None,
    ) -> list[dict[str, Any]]:
        project_path = Path(self.starter_preparation["project_path"])
        results = []
        for index, verifier in enumerate(step.get("verifiers", []), start=1):
            kind = verifier["kind"]
            label = f"runtime-{step['id']}-{kind}-{index}"
            if kind == "visible_claim":
                if not self.visible_claim_holds(snapshot, verifier, scope):
                    raise AcceptanceFailure(
                        "product_failure", f"Visible claim failed in step '{step['id']}'"
                    )
                semantic_id = verifier["semantic_id"]
                item = self.item_for(
                    snapshot,
                    semantic_id,
                    subject_scope=scope,
                    allow_unscoped_fallback=semantic_id.startswith("window."),
                )
                results.append(
                    {
                        "kind": kind,
                        "status": "pass",
                        "semantic_id": verifier["semantic_id"],
                        "observed": item,
                    }
                )
            elif kind == "facts":
                output = self.fact_eval(project_path, verifier["expression"], label)
                if output.get("truth") != "satisfied":
                    raise AcceptanceFailure(
                        "product_failure", f"Fact verifier failed in step '{step['id']}'"
                    )
                results.append({"kind": kind, "status": "pass", "truth": "satisfied"})
            elif kind == "expected_effects":
                results.append(
                    {"kind": kind, **self.expected_effects_verify(project_path, verifier, label)}
                )
            elif kind == "report":
                results.append({"kind": kind, **self.show_report(project_path, verifier, label)})
            elif kind == "state":
                results.append(
                    {
                        "kind": kind,
                        **self.state_verify(project_path, verifier, "starter", label),
                    }
                )
            elif kind == "artifact":
                results.append(
                    {
                        "kind": kind,
                        **self.artifact_verify(verifier, self.chapter_dir),
                    }
                )
            else:
                raise AcceptanceFailure(
                    "tutorial_ambiguity", f"Unsupported verifier kind '{kind}'"
                )
        return results

    def retain_step_evidence(
        self, step: dict[str, Any], snapshot: dict[str, Any]
    ) -> dict[str, Any]:
        evidence: dict[str, Any] = {}
        policy = step["evidence"]
        snapshot_requirement = policy["snapshot"]
        if snapshot_requirement != "omitted":
            snapshot_path = self.chapter_dir / "checkpoints" / f"{step['id']}.snapshot.json"
            atomic_write_json(snapshot_path, snapshot)
            evidence["snapshot"] = {
                "requirement": snapshot_requirement,
                "path": str(snapshot_path),
                "sha256": sha256_file(snapshot_path),
            }
        screenshot_requirement = policy["screenshot"]
        if screenshot_requirement != "omitted" and self.args.scrot is not None:
            screenshot_path = self.chapter_dir / "checkpoints" / f"{step['id']}.png"
            screenshot_path.parent.mkdir(parents=True, exist_ok=True)
            completed = subprocess.run(
                [str(self.args.scrot), str(screenshot_path)],
                cwd=self.repo_root,
                env=self.process_environment,
                timeout=20,
                check=False,
            )
            if completed.returncode != 0 or not screenshot_path.is_file():
                if screenshot_requirement == "required":
                    raise AcceptanceFailure(
                        "harness_gap", f"Required screenshot failed for step '{step['id']}'"
                    )
            else:
                evidence["screenshot"] = {
                    "requirement": screenshot_requirement,
                    "path": str(screenshot_path),
                    "sha256": sha256_file(screenshot_path),
                }
        elif screenshot_requirement == "required":
            raise AcceptanceFailure(
                "missing_dependency", f"Step '{step['id']}' requires scrot"
            )
        return evidence

    def execute_steps(self) -> None:
        starter_path = Path(self.starter_preparation["project_path"])
        for step in self.acceptance["steps"]:
            started = time.monotonic()
            timeout = self.args.timeouts[step["timeout_class"]]
            scope = self.scope_for_step(step)
            step_record: dict[str, Any] = {
                "id": step["id"],
                "step_sha256": sha256_bytes(canonical_json_bytes(step)),
                "prose_step": step["prose_step"],
                "requested_action": step["interaction"],
                "semantic_target": step["target"],
                "window": step["window"],
                "subject_scope": scope,
                "timeout_class": step["timeout_class"],
                "status": "running",
            }
            self.steps.append(step_record)
            self.write_ledger()
            before_snapshot = self.wait_snapshot(
                lambda snapshot: self.target_ready(snapshot, step, scope),
                timeout,
                f"enabled target for step '{step['id']}'",
            )
            target = self.item_for(
                before_snapshot,
                step["target"],
                window_id=step["window"],
                subject_scope=scope,
            )
            if target is None:
                raise AcceptanceFailure(
                    "harness_gap", f"Target vanished for step '{step['id']}'"
                )
            step_record["before_generation"] = before_snapshot["generation"]
            step_record["resolved_target"] = target
            if step.get("scientific_effect"):
                before = self.fact_eval(
                    starter_path, step["before"], f"runtime-{step['id']}-before"
                )
                expected = {"unsatisfied"}
                if step.get("allow_preexisting"):
                    expected.add("satisfied")
                if before.get("truth") not in expected:
                    raise AcceptanceFailure(
                        "product_failure",
                        f"Before fact for step '{step['id']}' is {before.get('truth')!r}",
                    )
                step_record["before_fact"] = before
            step_record["x11_event"] = self.emit_x11(step, target)
            after_snapshot = self.wait_after_interaction(
                step, scope, before_snapshot["generation"], timeout
            )
            step_record["after_generation"] = after_snapshot["generation"]
            if step.get("persists_project_state"):
                step_record["save"] = self.save_project(
                    after_snapshot["generation"], self.args.timeouts["io"]
                )
                after_snapshot = self.wait_snapshot(
                    lambda snapshot: snapshot.get("generation", 0)
                    >= step_record["save"]["generation"],
                    self.args.timeouts["instant"],
                    "post-save semantic snapshot",
                )
            if step.get("scientific_effect"):
                after = self.fact_eval(
                    starter_path, step["after"], f"runtime-{step['id']}-after"
                )
                if after.get("truth") != "satisfied":
                    raise AcceptanceFailure(
                        "product_failure",
                        f"After fact for step '{step['id']}' is {after.get('truth')!r}",
                    )
                step_record["after_fact"] = after
            step_record["verifiers"] = self.verify_step(step, after_snapshot, scope)
            verifier_kinds = {
                verifier["kind"] for verifier in step.get("verifiers", [])
            }
            step_record["visual_verdict"] = (
                "pass" if "visible_claim" in verifier_kinds else "not_requested"
            )
            step_record["scientific_verdict"] = (
                "pass"
                if any(kind != "visible_claim" for kind in verifier_kinds)
                else "not_requested"
            )
            step_record["evidence"] = self.retain_step_evidence(step, after_snapshot)
            step_record["elapsed_ms"] = round((time.monotonic() - started) * 1000)
            step_record["status"] = "pass"
            self.write_ledger()

    def run(self) -> dict[str, Any]:
        try:
            self.preflight_contract()
            self.launch_gui()
            self.execute_steps()
            completion = self.fact_eval(
                Path(self.starter_preparation["project_path"]),
                self.acceptance["completion_condition"],
                "runtime-completion",
            )
            if completion.get("truth") != "satisfied":
                raise AcceptanceFailure(
                    "product_failure", "Final tutorial completion fact is not satisfied"
                )
            self.ledger["completion"] = completion
            self.ledger["final_project_sha256"] = sha256_file(
                Path(self.starter_preparation["project_path"])
            )
            self.ledger["status"] = "pass"
            self.ledger["message"] = "GUI actions and typed scientific verification passed"
        except AcceptanceFailure as error:
            self.ledger["status"] = "fail"
            self.ledger["failure_class"] = error.failure_class
            self.ledger["message"] = str(error)
            if self.steps and self.steps[-1].get("status") == "running":
                self.steps[-1]["status"] = "fail"
                self.steps[-1]["failure_class"] = error.failure_class
                self.steps[-1]["message"] = str(error)
        except Exception as error:  # retain unexpected harness failures as evidence
            self.ledger["status"] = "fail"
            self.ledger["failure_class"] = "harness_gap"
            self.ledger["message"] = f"Unexpected runner failure: {type(error).__name__}: {error}"
        finally:
            self.stop_gui()
            self.write_ledger()
        return self.ledger


def tool_version(argv: list[str]) -> dict[str, Any]:
    try:
        completed = subprocess.run(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"argv": argv, "available": False, "diagnostic": str(error)}
    output = completed.stdout.decode("utf-8", errors="replace").strip()
    return {
        "argv": argv,
        "available": completed.returncode == 0,
        "exit_code": completed.returncode,
        "output": output[:1000],
    }


def detect_window_manager(xprop: Path) -> dict[str, Any]:
    root = tool_version([str(xprop), "-root", "_NET_SUPPORTING_WM_CHECK"])
    output = str(root.get("output", ""))
    match = re.search(r"0x[0-9a-fA-F]+", output)
    if not root.get("available") or match is None:
        return {
            "available": False,
            "diagnostic": "No EWMH window manager was detected on the X11 root window",
            "root_probe": root,
        }
    window_id = match.group(0)
    name = tool_version([str(xprop), "-id", window_id, "_NET_WM_NAME"])
    if not name.get("available"):
        return {
            "available": False,
            "diagnostic": f"Could not read the window-manager name from {window_id}",
            "root_probe": root,
            "name_probe": name,
        }
    return {
        "available": True,
        "window_id": window_id,
        "name": name.get("output", ""),
        "root_probe": root,
        "name_probe": name,
    }


def build_environment_record(
    args: argparse.Namespace, repo_root: Path
) -> dict[str, Any]:
    inherited, cleared_names = redacted_environment(dict(os.environ))
    git_revision = tool_version(["git", "-C", str(repo_root), "rev-parse", "HEAD"])
    git_status = tool_version(["git", "-C", str(repo_root), "status", "--short"])
    tools = {
        "xdotool": tool_version([str(args.xdotool), "version"]),
        "xdpyinfo": tool_version([str(args.xdpyinfo)]),
        "xprop": tool_version([str(args.xprop), "-version"]),
        "window_manager": detect_window_manager(args.xprop),
        "scrot": (
            tool_version([str(args.scrot), "--version"])
            if args.scrot is not None
            else {"available": False}
        ),
        "rustc": tool_version(["rustc", "--version"]),
        "python": {"available": True, "version": sys.version},
        "primer3": tool_version(
            [str(args.gentle_cli), "primers", "preflight", "--backend", "primer3"]
        ),
        "blastn": optional_tool_version("blastn", ["-version"]),
        "blastdbcmd": optional_tool_version("blastdbcmd", ["-version"]),
        "pandoc": optional_tool_version("pandoc", ["--version"]),
        "ghostscript": optional_tool_version("gs", ["--version"]),
    }
    return {
        "schema": ENVIRONMENT_SCHEMA,
        "source_revision": git_revision.get("output"),
        "git_status": git_status.get("output", ""),
        "os": platform.platform(),
        "kernel": platform.release(),
        "machine": platform.machine(),
        "display": os.environ.get("DISPLAY"),
        "display_geometry": tools["xdpyinfo"],
        "locale": {"LANG": "C.UTF-8", "LC_ALL": "C.UTF-8", "TZ": "UTC"},
        "profile_isolation": [
            "HOME",
            "XDG_CONFIG_HOME",
            "XDG_CACHE_HOME",
            "XDG_DATA_HOME",
            "TMPDIR",
        ],
        "cleared_inherited_variables": inherited,
        "cleared_variable_names": sorted(cleared_names),
        "network_enforcement": args.network_enforcement,
        "binaries": {
            "gentle": {
                "path": str(args.gentle),
                "sha256": sha256_file(args.gentle),
            },
            "gentle_cli": {
                "path": str(args.gentle_cli),
                "sha256": sha256_file(args.gentle_cli),
            },
            "gentle_examples_docs": {
                "path": str(args.examples_docs),
                "sha256": sha256_file(args.examples_docs),
            },
        },
        "cargo_lock_sha256": sha256_file(repo_root / "Cargo.lock"),
        "tools": tools,
    }


def resolve_tool(value: str | None, name: str, required: bool) -> Path | None:
    resolved = value or shutil.which(name)
    if resolved is None:
        if required:
            raise AcceptanceFailure("missing_dependency", f"Required tool '{name}' is absent")
        return None
    path = Path(resolved).expanduser().resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise AcceptanceFailure("missing_dependency", f"Tool is not executable: {path}")
    return path


def optional_tool_version(name: str, version_args: list[str]) -> dict[str, Any]:
    resolved = shutil.which(name)
    if resolved is None:
        return {"available": False, "name": name}
    return tool_version([resolved, *version_args])


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, default=Path("docs/tutorial/manifest.json"))
    parser.add_argument(
        "--workflow-source", type=Path, default=Path("docs/examples/workflows")
    )
    parser.add_argument("--chapter", action="append", default=[])
    parser.add_argument("--profile", default="smoke")
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--gentle", type=Path, default=Path("target/debug/gentle"))
    parser.add_argument(
        "--gentle-cli", type=Path, default=Path("target/debug/gentle_cli")
    )
    parser.add_argument(
        "--examples-docs", type=Path, default=Path("target/debug/gentle_examples_docs")
    )
    parser.add_argument("--xdotool")
    parser.add_argument("--xdpyinfo")
    parser.add_argument("--xprop")
    parser.add_argument("--scrot")
    parser.add_argument(
        "--network-enforcement",
        choices=["not_enforced", "linux_network_namespace", "container_none", "external_policy"],
        default="not_enforced",
        help="How the caller enforces an offline tutorial's no-network boundary.",
    )
    for timeout_class, default in TIMEOUT_DEFAULTS.items():
        parser.add_argument(
            f"--timeout-{timeout_class}", type=float, default=default, metavar="SECONDS"
        )
    parsed = parser.parse_args(argv)
    parsed.repo_root = parsed.repo_root.resolve()
    parsed.manifest = (parsed.repo_root / parsed.manifest).resolve() if not parsed.manifest.is_absolute() else parsed.manifest.resolve()
    parsed.workflow_source = (parsed.repo_root / parsed.workflow_source).resolve() if not parsed.workflow_source.is_absolute() else parsed.workflow_source.resolve()
    parsed.evidence_dir = parsed.evidence_dir.resolve()
    for field in ("gentle", "gentle_cli", "examples_docs"):
        path = getattr(parsed, field)
        path = (parsed.repo_root / path).resolve() if not path.is_absolute() else path.resolve()
        setattr(parsed, field, path)
    parsed.xdotool = resolve_tool(parsed.xdotool, "xdotool", required=True)
    parsed.xdpyinfo = resolve_tool(parsed.xdpyinfo, "xdpyinfo", required=True)
    parsed.xprop = resolve_tool(parsed.xprop, "xprop", required=True)
    parsed.scrot = resolve_tool(parsed.scrot, "scrot", required=False)
    parsed.timeouts = {
        timeout_class: getattr(parsed, f"timeout_{timeout_class}")
        for timeout_class in TIMEOUT_DEFAULTS
    }
    return parsed


def validate_runtime(args: argparse.Namespace, manifest: dict[str, Any]) -> list[dict[str, Any]]:
    if platform.system() != "Linux":
        raise AcceptanceFailure(
            "missing_dependency", "Tutorial GUI acceptance currently supports Linux/X11 only"
        )
    if not os.environ.get("DISPLAY"):
        raise AcceptanceFailure(
            "missing_dependency", "DISPLAY is unset; run this under Xvfb/X11"
        )
    window_manager = detect_window_manager(args.xprop)
    if not window_manager.get("available"):
        raise AcceptanceFailure(
            "harness_gap",
            "No named EWMH window manager is active; start one inside the Xvfb session",
        )
    for path in (args.gentle, args.gentle_cli, args.examples_docs):
        if not path.is_file() or not os.access(path, os.X_OK):
            raise AcceptanceFailure("missing_dependency", f"Binary is not executable: {path}")
    chapters = selected_chapters(manifest, args.chapter, args.profile)
    if any(chapter["gui_acceptance"].get("network") == "offline" for chapter in chapters):
        if args.network_enforcement == "not_enforced":
            raise AcceptanceFailure(
                "harness_gap",
                "Offline acceptance requires an explicit network-enforcement method",
            )
        if args.network_enforcement == "linux_network_namespace":
            try:
                own_namespace = os.readlink("/proc/self/ns/net")
                init_namespace = os.readlink("/proc/1/ns/net")
            except OSError as error:
                raise AcceptanceFailure(
                    "harness_gap", f"Could not verify Linux network namespace: {error}"
                ) from error
            if own_namespace == init_namespace:
                raise AcceptanceFailure(
                    "harness_gap",
                    "linux_network_namespace was declared but the runner shares PID 1's network namespace",
                )
    if any(
        step.get("evidence", {}).get("screenshot") == "required"
        for chapter in chapters
        for step in chapter["gui_acceptance"].get("steps", [])
    ) and args.scrot is None:
        raise AcceptanceFailure(
            "missing_dependency", "Selected acceptance contract requires scrot"
        )
    return chapters


def main(argv: Iterable[str] | None = None) -> int:
    args: argparse.Namespace | None = None
    try:
        args = parse_args(argv)
        manifest = load_json(args.manifest)
        if manifest.get("schema") != "gentle.tutorial_manifest.v2":
            raise AcceptanceFailure(
                "tutorial_ambiguity", f"Unsupported tutorial manifest schema {manifest.get('schema')!r}"
            )
        chapters = validate_runtime(args, manifest)
        args.evidence_dir.mkdir(parents=True, exist_ok=True)
        environment_record = build_environment_record(args, args.repo_root)
        atomic_write_json(args.evidence_dir / "environment.json", environment_record)
        ledgers = []
        for chapter in chapters:
            run = TutorialAcceptanceRun(args, args.repo_root, chapter, environment_record)
            ledgers.append(run.run())
        status = "pass" if all(ledger["status"] == "pass" for ledger in ledgers) else "fail"
        report = {
            "schema": RUN_SCHEMA,
            "status": status,
            "profile": args.profile,
            "chapter_count": len(ledgers),
            "chapters": [
                {
                    "chapter_id": ledger["chapter_id"],
                    "status": ledger["status"],
                    "failure_class": ledger.get("failure_class"),
                    "message": ledger.get("message"),
                    "ledger_path": str(args.evidence_dir / ledger["chapter_id"] / "acceptance-ledger.json"),
                    "ledger_sha256": sha256_file(
                        args.evidence_dir / ledger["chapter_id"] / "acceptance-ledger.json"
                    ),
                }
                for ledger in ledgers
            ],
        }
        atomic_write_json(args.evidence_dir / "acceptance-report.json", report)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if status == "pass" else 1
    except AcceptanceFailure as error:
        report = {
            "schema": RUN_SCHEMA,
            "status": "fail",
            "failure_class": error.failure_class,
            "message": str(error),
        }
        evidence_dir = None if args is None else args.evidence_dir
        if evidence_dir is not None:
            evidence_dir.mkdir(parents=True, exist_ok=True)
            atomic_write_json(evidence_dir / "acceptance-report.json", report)
        print(json.dumps(report, indent=2, sort_keys=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
