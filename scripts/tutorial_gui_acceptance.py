#!/usr/bin/env python3
"""Run typed GENtle tutorial GUI acceptance contracts through ordinary X11 input.

The tutorial manifest supplies intent and typed postconditions. This runner owns
input delivery, independent verification, and the final pass/fail verdict. It
never executes prose or command strings from a tutorial.
"""

from __future__ import annotations

import argparse
import hashlib
import html
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
SCREENSHOT_EVIDENCE_SCHEMA = "gentle.tutorial_gui_screenshot_evidence.v1"
SEMANTIC_COORDINATE_SPACE = (
    "egui logical points translated by the viewport outer rectangle; "
    "native window decorations are excluded"
)

TIMEOUT_DEFAULTS = {
    "instant": 5.0,
    "interactive": 25.0,
    "io": 75.0,
    "compute": 600.0,
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


def should_flush_pending_save_before_step(
    pending_project_save: bool, step: dict[str, Any]
) -> bool:
    """Keep an old dirty flag from masquerading as scientific completion."""
    return pending_project_save and bool(step.get("scientific_effect"))


class AcceptanceFailure(RuntimeError):
    def __init__(self, failure_class: str, message: str):
        super().__init__(message)
        self.failure_class = failure_class


def expected_starter_completion_truth(acceptance_contract: dict[str, Any]) -> str:
    """Return the typed starter invariant for mutating versus view-only tutorials."""
    return "satisfied" if acceptance_contract.get("view_only", False) else "unsatisfied"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def png_dimensions(path: Path) -> tuple[int, int]:
    """Read PNG canvas dimensions without adding an image-library dependency."""
    header = path.read_bytes()[:24]
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise AcceptanceFailure("harness_gap", f"Screenshot is not a readable PNG: {path}")
    width = int.from_bytes(header[16:20], "big")
    height = int.from_bytes(header[20:24], "big")
    if width <= 0 or height <= 0:
        raise AcceptanceFailure("harness_gap", f"Screenshot has invalid dimensions: {path}")
    return width, height


def semantic_pixel_rect(item: dict[str, Any]) -> dict[str, int]:
    rectangle = item.get("rect_logical_points", {})
    scale = float(item.get("pixels_per_point", 1.0))
    min_x = round(float(rectangle["min_x"]) * scale)
    min_y = round(float(rectangle["min_y"]) * scale)
    max_x = round(float(rectangle["max_x"]) * scale)
    max_y = round(float(rectangle["max_y"]) * scale)
    return {
        "min_x": min(min_x, max_x),
        "min_y": min(min_y, max_y),
        "max_x": max(min_x, max_x),
        "max_y": max(min_y, max_y),
    }


def padded_crop_rect(
    rectangle: dict[str, int], canvas_width: int, canvas_height: int
) -> dict[str, int]:
    target_width = max(1, rectangle["max_x"] - rectangle["min_x"])
    target_height = max(1, rectangle["max_y"] - rectangle["min_y"])
    pad_x = max(180, target_width)
    pad_y = max(120, target_height * 2)
    min_x = max(0, rectangle["min_x"] - pad_x)
    min_y = max(0, rectangle["min_y"] - pad_y)
    max_x = min(canvas_width, rectangle["max_x"] + pad_x)
    max_y = min(canvas_height, rectangle["max_y"] + pad_y)
    return {
        "min_x": min_x,
        "min_y": min_y,
        "max_x": max(min_x + 1, max_x),
        "max_y": max(min_y + 1, max_y),
    }


def write_screenshot_view_svg(
    output: Path,
    raw_png: Path,
    canvas_width: int,
    canvas_height: int,
    crop: dict[str, int],
    focus: dict[str, int],
    label: str,
) -> None:
    """Create a lossless teaching view over one immutable raw X11 capture."""
    crop_width = crop["max_x"] - crop["min_x"]
    crop_height = crop["max_y"] - crop["min_y"]
    relative_raw = os.path.relpath(raw_png, output.parent)
    focus_x = focus["min_x"] - crop["min_x"]
    focus_y = focus["min_y"] - crop["min_y"]
    focus_width = max(1, focus["max_x"] - focus["min_x"])
    focus_height = max(1, focus["max_y"] - focus["min_y"])
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{crop_width}" height="{crop_height}" viewBox="0 0 {crop_width} {crop_height}">
  <title>{html.escape(label)}</title>
  <image href="{html.escape(relative_raw)}" x="{-crop['min_x']}" y="{-crop['min_y']}" width="{canvas_width}" height="{canvas_height}"/>
  <rect x="{focus_x}" y="{focus_y}" width="{focus_width}" height="{focus_height}" fill="none" stroke="#d62728" stroke-width="3"/>
  <circle cx="{focus_x + 12}" cy="{focus_y + 12}" r="11" fill="#d62728"/>
  <text x="{focus_x + 12}" y="{focus_y + 17}" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="bold" fill="white">1</text>
</svg>
'''
    output.write_text(svg, encoding="utf-8")


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def sequence_content_identity(project: dict[str, Any], seq_id: str) -> str:
    """Compare persisted biological content, not volatile runtime caches or labels."""
    try:
        dna = project["sequences"][seq_id]
        sequence = dna["seq"]
        content = {
            "bases": sequence["seq"],
            "topology": sequence["topology"],
            "features": sequence["features"],
            "molecule_type": sequence.get("molecule_type"),
            "overhang": dna["overhang"],
        }
    except (KeyError, TypeError) as error:
        raise AcceptanceFailure("product_failure", f"Missing sequence content for {seq_id}") from error
    return sha256_bytes(canonical_json_bytes(content))


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


@dataclass(frozen=True)
class NativeWindowGeometry:
    """Root-screen geometry for one exact X11 client window."""

    window_id: int
    root_x: int
    root_y: int
    width: int
    height: int


def parse_xwininfo_geometry(output: str, window_id: int) -> NativeWindowGeometry:
    """Parse the stable, locale-independent numeric fields used from xwininfo."""

    fields: dict[str, int] = {}
    patterns = {
        "root_x": r"^\s*Absolute upper-left X:\s*(-?\d+)\s*$",
        "root_y": r"^\s*Absolute upper-left Y:\s*(-?\d+)\s*$",
        "width": r"^\s*Width:\s*(\d+)\s*$",
        "height": r"^\s*Height:\s*(\d+)\s*$",
    }
    for name, pattern in patterns.items():
        match = re.search(pattern, output, flags=re.MULTILINE)
        if match is None:
            raise AcceptanceFailure(
                "harness_gap", f"xwininfo omitted {name} for X11 window {window_id}"
            )
        fields[name] = int(match.group(1))
    if fields["width"] <= 0 or fields["height"] <= 0:
        raise AcceptanceFailure(
            "harness_gap", f"X11 window {window_id} has invalid client geometry"
        )
    return NativeWindowGeometry(window_id=window_id, **fields)


def native_screen_rect(
    item: dict[str, Any],
    semantic_window: dict[str, Any],
    native: NativeWindowGeometry,
) -> dict[str, int]:
    """Map one egui semantic rectangle into X11 root-screen physical pixels."""

    item_scale = float(item.get("pixels_per_point", 1.0))
    window_scale = float(semantic_window.get("pixels_per_point", 1.0))
    if abs(item_scale - window_scale) > 1e-6:
        raise AcceptanceFailure(
            "harness_gap", "Target and semantic window use different pixel scales"
        )
    target = semantic_pixel_rect(item)
    viewport = semantic_pixel_rect(semantic_window)
    local = {
        "min_x": target["min_x"] - viewport["min_x"],
        "min_y": target["min_y"] - viewport["min_y"],
        "max_x": target["max_x"] - viewport["min_x"],
        "max_y": target["max_y"] - viewport["min_y"],
    }
    tolerance = 2
    if (
        local["min_x"] < -tolerance
        or local["min_y"] < -tolerance
        or local["max_x"] > native.width + tolerance
        or local["max_y"] > native.height + tolerance
    ):
        raise AcceptanceFailure(
            "harness_gap",
            f"Semantic target lies outside bound X11 client window {native.window_id}",
        )
    return {
        "min_x": native.root_x + local["min_x"],
        "min_y": native.root_y + local["min_y"],
        "max_x": native.root_x + local["max_x"],
        "max_y": native.root_y + local["max_y"],
    }


def validate_parent_network_namespace(own: str, parent: str | None) -> None:
    """Prove that the runner entered a namespace distinct from its caller."""

    namespace_pattern = re.compile(r"^net:\[\d+\]$")
    if parent is None or not namespace_pattern.fullmatch(parent):
        raise AcceptanceFailure(
            "harness_gap",
            "linux_network_namespace requires the parent namespace identity "
            "captured before unshare",
        )
    if not namespace_pattern.fullmatch(own):
        raise AcceptanceFailure(
            "harness_gap", f"Unexpected current network namespace identity: {own!r}"
        )
    if own == parent:
        raise AcceptanceFailure(
            "harness_gap",
            "linux_network_namespace was declared but the runner shares its "
            "caller's network namespace",
        )


def parse_ewmh_client_ids(output: str) -> list[int]:
    """Parse the root window's ordered EWMH client inventory."""

    if "not found" in output.lower() or "no such atom" in output.lower():
        return []
    marker = "window id #"
    if marker not in output:
        raise AcceptanceFailure(
            "harness_gap", "Window manager did not expose a valid _NET_CLIENT_LIST"
        )
    values = output.split(marker, 1)[1]
    ids: list[int] = []
    for value in values.split(","):
        value = value.strip()
        if value:
            try:
                ids.append(int(value, 0))
            except ValueError as error:
                raise AcceptanceFailure(
                    "harness_gap", f"Invalid X11 client id in _NET_CLIENT_LIST: {value!r}"
                ) from error
    return ids


def parse_ewmh_process_id(output: str, window_id: int) -> int | None:
    """Parse one EWMH client PID without relying on xdotool's process lookup."""

    if "not found" in output.lower():
        return None
    match = re.search(r"_NET_WM_PID\(CARDINAL\)\s*=\s*(\d+)\s*$", output)
    if match is None:
        raise AcceptanceFailure(
            "harness_gap", f"X11 client {window_id} has malformed _NET_WM_PID"
        )
    return int(match.group(1))


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
        self.native_window_bindings: dict[str, int] = {}
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

    def run_x11(
        self, argv: list[str], *, timeout: float = 10.0
    ) -> subprocess.CompletedProcess[bytes]:
        try:
            completed = subprocess.run(
                argv,
                cwd=self.repo_root,
                env=self.process_environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise AcceptanceFailure(
                "harness_gap", f"X11 command could not complete: {argv}: {error}"
            ) from error
        if completed.returncode != 0:
            diagnostic = completed.stderr.decode("utf-8", errors="replace").strip()
            raise AcceptanceFailure(
                "harness_gap",
                f"X11 command exited {completed.returncode}: {argv}: {diagnostic[-500:]}",
            )
        return completed

    def visible_native_window_ids(self) -> list[int]:
        if self.gui_process is None:
            raise AcceptanceFailure("harness_gap", "GENtle process is not active")
        deadline = time.monotonic() + 3.0
        while time.monotonic() < deadline:
            root = self.run_x11(
                [str(self.args.xprop), "-root", "_NET_CLIENT_LIST"]
            )
            candidates = parse_ewmh_client_ids(
                root.stdout.decode("utf-8", errors="replace")
            )
            matches: list[int] = []
            for window_id in candidates:
                owner = self.run_x11(
                    [str(self.args.xprop), "-id", str(window_id), "_NET_WM_PID"]
                )
                if (
                    parse_ewmh_process_id(
                        owner.stdout.decode("utf-8", errors="replace"), window_id
                    )
                    == self.gui_process.pid
                ):
                    geometry = self.query_native_geometry(window_id)
                    matches.append(geometry.window_id)
            if matches:
                return matches
            time.sleep(0.1)
        return []

    @staticmethod
    def semantic_window_ids(snapshot: dict[str, Any]) -> list[str]:
        return sorted(
            {
                str(item["semantic_id"])
                for item in snapshot.get("items", [])
                if item.get("widget_kind") == "window"
                and item.get("semantic_id") == item.get("window_id")
                and item.get("state", {}).get("visible")
            }
        )

    def refresh_native_window_bindings(self, snapshot: dict[str, Any]) -> None:
        native_ids = self.visible_native_window_ids()
        semantic_ids = self.semantic_window_ids(snapshot)
        visible_native = set(native_ids)
        visible_semantic = set(semantic_ids)
        self.native_window_bindings = {
            semantic_id: native_id
            for semantic_id, native_id in self.native_window_bindings.items()
            if semantic_id in visible_semantic and native_id in visible_native
        }
        viewport_by_semantic = {
            str(item["semantic_id"]): str(item.get("egui_viewport_id", ""))
            for item in snapshot.get("items", [])
            if item.get("widget_kind") == "window"
            and item.get("semantic_id") == item.get("window_id")
        }
        for semantic_id in semantic_ids:
            if semantic_id in self.native_window_bindings:
                continue
            viewport_id = viewport_by_semantic.get(semantic_id)
            inherited = {
                native_id
                for bound_semantic, native_id in self.native_window_bindings.items()
                if viewport_id
                and viewport_by_semantic.get(bound_semantic) == viewport_id
            }
            if len(inherited) == 1:
                self.native_window_bindings[semantic_id] = inherited.pop()
            elif len(inherited) > 1:
                raise AcceptanceFailure(
                    "harness_gap",
                    f"Egui viewport {viewport_id!r} maps to multiple X11 clients",
                )
        bound_native = set(self.native_window_bindings.values())
        unbound_native = [window_id for window_id in native_ids if window_id not in bound_native]
        unbound_semantic = [
            semantic_id
            for semantic_id in semantic_ids
            if semantic_id not in self.native_window_bindings
        ]
        if not unbound_semantic:
            return
        if len(unbound_native) != len(unbound_semantic):
            raise AcceptanceFailure(
                "harness_gap",
                "Could not bind semantic viewports to exact GENtle X11 clients: "
                f"semantic={unbound_semantic}, native={unbound_native}",
            )
        if len(unbound_semantic) > 1:
            raise AcceptanceFailure(
                "harness_gap",
                "Multiple semantic and native GENtle windows appeared simultaneously; "
                "their identities are ambiguous",
            )
        self.native_window_bindings[unbound_semantic[0]] = unbound_native[0]

    def query_native_geometry(self, window_id: int) -> NativeWindowGeometry:
        completed = self.run_x11(
            [str(self.args.xwininfo), "-id", str(window_id), "-stats"]
        )
        return parse_xwininfo_geometry(
            completed.stdout.decode("utf-8", errors="replace"), window_id
        )

    def active_native_window_id(self) -> int:
        return int(
            self.run_x11([str(self.args.xdotool), "getactivewindow"])
            .stdout.decode("ascii")
            .strip()
        )

    def focused_native_window_id(self) -> int:
        return int(
            self.run_x11([str(self.args.xdotool), "getwindowfocus"])
            .stdout.decode("ascii")
            .strip()
        )

    def ensure_native_window_active(self, window_id: int) -> list[list[str]]:
        raise_window = [str(self.args.xdotool), "windowraise", str(window_id)]
        self.run_x11(raise_window)
        commands: list[list[str]] = [raise_window]
        if self.active_native_window_id() != window_id:
            activate = [
                str(self.args.xdotool),
                "windowactivate",
                "--sync",
                str(window_id),
            ]
            self.run_x11(activate)
            commands.append(activate)
        if self.focused_native_window_id() != window_id:
            focus = [
                str(self.args.xdotool),
                "windowfocus",
                "--sync",
                str(window_id),
            ]
            self.run_x11(focus)
            commands.append(focus)
        deadline = time.monotonic() + 2.0
        while time.monotonic() < deadline:
            if (
                self.active_native_window_id() == window_id
                and self.focused_native_window_id() == window_id
            ):
                time.sleep(0.1)
                return commands
            time.sleep(0.05)
        raise AcceptanceFailure(
            "harness_gap",
            f"X11 client {window_id} did not become both active and input-focused",
        )

    def stable_native_geometry(
        self, window_id: int, timeout: float = 3.0
    ) -> NativeWindowGeometry:
        deadline = time.monotonic() + timeout
        previous: NativeWindowGeometry | None = None
        while time.monotonic() < deadline:
            current = self.query_native_geometry(window_id)
            if current == previous:
                return current
            previous = current
            time.sleep(0.1)
        raise AcceptanceFailure(
            "harness_gap", f"X11 client geometry did not settle for window {window_id}"
        )

    def native_target_binding(
        self, snapshot: dict[str, Any], item: dict[str, Any]
    ) -> dict[str, Any]:
        if snapshot.get("coordinate_space") != SEMANTIC_COORDINATE_SPACE:
            raise AcceptanceFailure(
                "harness_gap",
                "Semantic snapshot does not declare the client-relative "
                "coordinate contract required for native X11 input",
            )
        self.refresh_native_window_bindings(snapshot)
        semantic_window_id = str(item.get("window_id", ""))
        semantic_window = self.item_for(
            snapshot, semantic_window_id, window_id=semantic_window_id
        )
        if semantic_window is None:
            viewport_id = item.get("egui_viewport_id")
            candidates = [
                candidate
                for candidate in snapshot.get("items", [])
                if candidate.get("widget_kind") == "window"
                and candidate.get("semantic_id") == candidate.get("window_id")
                and viewport_id
                and candidate.get("egui_viewport_id") == viewport_id
            ]
            if len(candidates) != 1:
                raise AcceptanceFailure(
                    "harness_gap",
                    f"Semantic surface '{semantic_window_id}' has no unique owning viewport",
                )
            semantic_window = candidates[0]
        native_owner_id = str(semantic_window["semantic_id"])
        native_window_id = self.native_window_bindings.get(native_owner_id)
        if native_window_id is None:
            raise AcceptanceFailure(
                "harness_gap",
                f"Semantic viewport '{native_owner_id}' has no exact X11 client binding",
            )
        native = self.stable_native_geometry(native_window_id)
        screen_rect = native_screen_rect(item, semantic_window, native)
        return {
            "semantic_window_id": semantic_window_id,
            "native_owner_semantic_window_id": native_owner_id,
            "egui_viewport_id": semantic_window.get("egui_viewport_id"),
            "x11_client_window_id": native_window_id,
            "native_client_geometry": {
                "root_x": native.root_x,
                "root_y": native.root_y,
                "width": native.width,
                "height": native.height,
            },
            "semantic_window_rect_logical_points": semantic_window.get(
                "rect_logical_points"
            ),
            "coordinate_transform": (
                "root_screen_px = native_client_origin_px + "
                "(target_logical - semantic_window_logical_origin) * pixels_per_point"
            ),
            "screen_rect_physical_pixels": screen_rect,
        }

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
        identities = {}
        if verifier.get("compare_with_oracle"):
            current = load_json(project_path)
            oracle = load_json(Path(self.oracle_preparation["project_path"]))
            oracle_mapping = self.binding_map(self.oracle_preparation)
            for source, actual in zip(verifier["seq_ids"], required):
                expected = oracle_mapping.get(source, source)
                identity = sequence_content_identity(current, actual)
                if identity != sequence_content_identity(oracle, expected):
                    raise AcceptanceFailure("product_failure", f"Sequence content differs from oracle: {actual}")
                identities[actual] = identity
        return {"status": "pass", "sequence_ids": required, "sequence_content_sha256": identities}

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
        expected_starter_truth = expected_starter_completion_truth(self.acceptance)
        if starter_completion.get("truth") != expected_starter_truth:
            failure_class = (
                "starter_precompleted"
                if not self.acceptance.get("view_only", False)
                and starter_completion.get("truth") == "satisfied"
                else "tutorial_ambiguity"
            )
            raise AcceptanceFailure(
                failure_class,
                f"Starter completion fact is {starter_completion.get('truth')!r}, "
                f"expected {expected_starter_truth!r}",
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
                elif verifier["kind"] == "state":
                    self.state_verify(oracle_path, verifier, "oracle", f"oracle-check-{step['id']}-state-{index}")
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

    def emit_x11(
        self,
        step: dict[str, Any],
        item: dict[str, Any],
        snapshot: dict[str, Any],
        timeout: float,
    ) -> dict[str, Any]:
        binding = self.native_target_binding(snapshot, item)
        native_window_id = binding["x11_client_window_id"]
        commands: list[list[str]] = []
        focus_commands = self.ensure_native_window_active(native_window_id)
        if focus_commands:
            commands.extend(focus_commands)
            binding = self.native_target_binding(snapshot, item)
        rectangle = binding["screen_rect_physical_pixels"]
        x = round((rectangle["min_x"] + rectangle["max_x"]) / 2)
        y = round((rectangle["min_y"] + rectangle["max_y"]) / 2)
        interaction = step["interaction"]
        kind = interaction["kind"]
        location = self.run_x11(
            [str(self.args.xdotool), "getmouselocation", "--shell"]
        ).stdout.decode("ascii", errors="strict")
        pointer = {
            key: int(value)
            for key, value in re.findall(r"^(X|Y)=(-?\d+)$", location, flags=re.MULTILINE)
        }
        if pointer.get("X") != x or pointer.get("Y") != y:
            move = [str(self.args.xdotool), "mousemove", str(x), str(y)]
            self.run_x11(move)
            commands.append(move)
            moved = self.run_x11(
                [str(self.args.xdotool), "getmouselocation", "--shell"]
            ).stdout.decode("ascii", errors="strict")
            moved_pointer = {
                key: int(value)
                for key, value in re.findall(
                    r"^(X|Y)=(-?\d+)$", moved, flags=re.MULTILINE
                )
            }
            if moved_pointer.get("X") != x or moved_pointer.get("Y") != y:
                raise AcceptanceFailure(
                    "harness_gap",
                    f"X11 pointer did not reach verified target ({x}, {y})",
                )
        hover_generation = None
        if kind != "press_key" and not (
            kind == "set_checkbox"
            and bool(item.get("state", {}).get("selected"))
            == bool(interaction["selected"])
        ):
            semantic_id = item["semantic_id"]
            window_id = item["window_id"]
            subject_scope = item.get("subject_scope")

            def target_hovered(candidate_snapshot: dict[str, Any]) -> bool:
                candidate = self.item_for(
                    candidate_snapshot,
                    semantic_id,
                    window_id=window_id,
                    subject_scope=subject_scope,
                )
                return bool(candidate and candidate.get("state", {}).get("hovered"))

            hovered_snapshot = self.wait_snapshot(
                target_hovered,
                min(timeout, self.args.timeouts["instant"]),
                f"semantic hover confirmation for step '{step['id']}'",
            )
            hovered_item = self.item_for(
                hovered_snapshot,
                semantic_id,
                window_id=window_id,
                subject_scope=subject_scope,
            )
            if hovered_item is None:
                raise AcceptanceFailure(
                    "harness_gap",
                    f"Target vanished while confirming hover for step '{step['id']}'",
                )
            binding = self.native_target_binding(hovered_snapshot, hovered_item)
            rectangle = binding["screen_rect_physical_pixels"]
            if not (
                rectangle["min_x"] <= x <= rectangle["max_x"]
                and rectangle["min_y"] <= y <= rectangle["max_y"]
            ):
                raise AcceptanceFailure(
                    "harness_gap",
                    f"Target geometry moved after pointer placement for step '{step['id']}'",
                )
            hover_generation = hovered_snapshot["generation"]
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
                    # Leaving the field commits the replacement and gives the
                    # child viewport a deterministic semantic publication.
                    [str(self.args.xdotool), "key", "--clearmodifiers", "Tab"],
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
            commands.append([str(self.args.xdotool), "key", "--clearmodifiers", key])
        elif kind == "scroll":
            direction = interaction.get("direction")
            button = {"up": "4", "down": "5"}.get(direction)
            clicks = interaction.get("clicks")
            if button is None or not isinstance(clicks, int) or not 1 <= clicks <= 20:
                raise AcceptanceFailure(
                    "tutorial_ambiguity",
                    f"Unsupported bounded scroll request {interaction!r}",
                )
            commands.append(
                [
                    str(self.args.xdotool),
                    "click",
                    "--repeat",
                    str(clicks),
                    "--delay",
                    "50",
                    button,
                ]
            )
        else:
            raise AcceptanceFailure(
                "tutorial_ambiguity", f"Unsupported interaction kind '{kind}'"
            )
        for command in commands:
            if command[1] not in {
                "windowactivate",
                "windowfocus",
                "windowraise",
                "mousemove",
            }:
                self.run_x11(command)
        return {
            "kind": kind,
            "screen_x": x,
            "screen_y": y,
            "native_window_binding": binding,
            "hover_confirmed_generation": hover_generation,
            "argv": commands,
        }

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
            return True

        return self.wait_snapshot(
            ready,
            timeout,
            f"postconditions for step '{step['id']}'",
            after_generation=prior_generation,
        )

    def save_project(self, snapshot: dict[str, Any], timeout: float) -> dict[str, Any]:
        save_item = self.item_for(snapshot, "main.project.save_state")
        if save_item is None:
            raise AcceptanceFailure(
                "harness_gap", "Project save status is absent from the semantic snapshot"
            )
        binding = self.native_target_binding(snapshot, save_item)
        native_window_id = binding["x11_client_window_id"]
        commands: list[list[str]] = []
        focus_commands = self.ensure_native_window_active(native_window_id)
        if focus_commands:
            commands.extend(focus_commands)
        save = [str(self.args.xdotool), "key", "--clearmodifiers", "ctrl+s"]
        self.run_x11(save)
        commands.append(save)

        def saved(snapshot: dict[str, Any]) -> bool:
            item = self.item_for(snapshot, "main.project.save_state")
            return bool(item and item.get("outcome_role") == "saved")

        snapshot = self.wait_snapshot(
            saved,
            timeout,
            "known-path project save",
            after_generation=snapshot["generation"],
        )
        return {
            "event": "ctrl+s",
            "native_window_binding": binding,
            "argv": commands,
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
            checkpoint_dir = self.chapter_dir / "checkpoints"
            screenshot_path = checkpoint_dir / f"{step['id']}.raw.png"
            screenshot_path.parent.mkdir(parents=True, exist_ok=True)
            # The semantic snapshot is written during egui's frame. Give the
            # X11 surface one bounded interval to present that verified frame
            # before capturing its pixels.
            time.sleep(0.15)
            completed = subprocess.run(
                [str(self.args.scrot), str(screenshot_path)],
                cwd=self.repo_root,
                env=self.process_environment,
                timeout=20,
                check=False,
            )
            # scrot briefly grabs the X server. Let the compositor and egui
            # observe the released grab before the next ordinary input event.
            time.sleep(0.5)
            if completed.returncode != 0 or not screenshot_path.is_file():
                if screenshot_requirement == "required":
                    raise AcceptanceFailure(
                        "harness_gap", f"Required screenshot failed for step '{step['id']}'"
                    )
            else:
                canvas_width, canvas_height = png_dimensions(screenshot_path)
                focus_item = self.screenshot_focus_item(step, snapshot)
                if focus_item is None:
                    raise AcceptanceFailure(
                        "harness_gap",
                        f"Screenshot step '{step['id']}' has no visible semantic focus item",
                    )
                focus_binding = self.native_target_binding(snapshot, focus_item)
                focus_rect = focus_binding["screen_rect_physical_pixels"]
                focus_rect = {
                    "min_x": max(0, min(canvas_width - 1, focus_rect["min_x"])),
                    "min_y": max(0, min(canvas_height - 1, focus_rect["min_y"])),
                    "max_x": max(1, min(canvas_width, focus_rect["max_x"])),
                    "max_y": max(1, min(canvas_height, focus_rect["max_y"])),
                }
                context_rect = padded_crop_rect(focus_rect, canvas_width, canvas_height)
                orientation_rect = {
                    "min_x": 0,
                    "min_y": 0,
                    "max_x": canvas_width,
                    "max_y": canvas_height,
                }
                orientation_path = checkpoint_dir / f"{step['id']}.orientation.svg"
                context_path = checkpoint_dir / f"{step['id']}.context.svg"
                write_screenshot_view_svg(
                    orientation_path,
                    screenshot_path,
                    canvas_width,
                    canvas_height,
                    orientation_rect,
                    focus_rect,
                    f"{self.chapter['id']} / {step['id']} — orientation",
                )
                write_screenshot_view_svg(
                    context_path,
                    screenshot_path,
                    canvas_width,
                    canvas_height,
                    context_rect,
                    focus_rect,
                    f"{self.chapter['id']} / {step['id']} — interaction context",
                )
                snapshot_sha256 = sha256_bytes(canonical_json_bytes(snapshot))
                screenshot_record = {
                    "schema": SCREENSHOT_EVIDENCE_SCHEMA,
                    "requirement": screenshot_requirement,
                    "chapter_id": self.chapter["id"],
                    "step_id": step["id"],
                    "prose_step": step["prose_step"],
                    "step_sha256": sha256_bytes(canonical_json_bytes(step)),
                    "acceptance_contract_sha256": self.ledger[
                        "acceptance_contract_sha256"
                    ],
                    "tutorial_manifest_sha256": self.ledger["manifest_sha256"],
                    "source_revision": self.ledger["environment"].get("source_revision"),
                    "gentle_binary_sha256": self.ledger["environment"]["binaries"][
                        "gentle"
                    ]["sha256"],
                    "project_state": {
                        "path": self.starter_preparation["project_path"],
                        "sha256": sha256_file(
                            Path(self.starter_preparation["project_path"])
                        ),
                    },
                    "capture": {
                        "backend": "scrot",
                        "backend_version": self.ledger["environment"]["tools"].get(
                            "scrot"
                        ),
                        "scope": "x11_root",
                        "captured_at_unix_ms": int(time.time() * 1000),
                        "display": self.ledger["environment"].get("display"),
                        "gui_process_pid": self.ledger.get("gui_process", {}).get("pid"),
                        "raw": {
                            "path": str(screenshot_path),
                            "sha256": sha256_file(screenshot_path),
                            "width_px": canvas_width,
                            "height_px": canvas_height,
                        },
                    },
                    "semantic_snapshot": {
                        "schema": snapshot.get("schema"),
                        "generation": snapshot.get("generation"),
                        "settled": snapshot.get("settled"),
                        "canonical_sha256": snapshot_sha256,
                        "retained_path": evidence.get("snapshot", {}).get("path"),
                        "retained_file_sha256": evidence.get("snapshot", {}).get("sha256"),
                    },
                    "requested_target": {
                        "semantic_id": step["target"],
                        "window_id": step["window"],
                        "subject_scope": self.scope_for_step(step),
                    },
                    "visual_focus": {
                        "semantic_id": focus_item.get("semantic_id"),
                        "window_id": focus_item.get("window_id"),
                        "subject_scope": focus_item.get("subject_scope"),
                        "rect_logical_points": focus_item.get("rect_logical_points"),
                        "pixels_per_point": focus_item.get("pixels_per_point"),
                        "rect_physical_pixels": focus_rect,
                        "native_window_binding": focus_binding,
                    },
                    "derived_views": [
                        {
                            "role": "orientation",
                            "path": str(orientation_path),
                            "sha256": sha256_file(orientation_path),
                            "crop_physical_pixels": orientation_rect,
                            "annotation": "Numbered outline derived from the semantic focus rectangle",
                        },
                        {
                            "role": "interaction_context",
                            "path": str(context_path),
                            "sha256": sha256_file(context_path),
                            "crop_physical_pixels": context_rect,
                            "annotation": "Padded crop around the same semantic focus rectangle",
                        },
                    ],
                }
                record_path = checkpoint_dir / f"{step['id']}.screenshot.json"
                atomic_write_json(record_path, screenshot_record)
                evidence["screenshot"] = {
                    "schema": SCREENSHOT_EVIDENCE_SCHEMA,
                    "requirement": screenshot_requirement,
                    "record_path": str(record_path),
                    "record_sha256": sha256_file(record_path),
                    "raw_path": str(screenshot_path),
                    "raw_sha256": screenshot_record["capture"]["raw"]["sha256"],
                    "derived_views": screenshot_record["derived_views"],
                }
        elif screenshot_requirement == "required":
            raise AcceptanceFailure(
                "missing_dependency", f"Step '{step['id']}' requires scrot"
            )
        return evidence

    def screenshot_focus_item(
        self, step: dict[str, Any], snapshot: dict[str, Any]
    ) -> dict[str, Any] | None:
        """Prefer the visible post-action result, then the requested control."""
        scope = self.scope_for_step(step)
        visible_ids = [
            verifier["semantic_id"]
            for verifier in step.get("verifiers", [])
            if verifier.get("kind") == "visible_claim"
        ]
        # Prefer a concrete control over a whole-window identity.
        visible_ids.sort(key=lambda semantic_id: semantic_id.startswith("window."))
        for semantic_id in visible_ids:
            item = self.item_for(
                snapshot,
                semantic_id,
                subject_scope=scope,
                allow_unscoped_fallback=semantic_id.startswith("window."),
            )
            if item is not None and item.get("state", {}).get("visible"):
                return item
        return self.item_for(
            snapshot,
            step["target"],
            window_id=step["window"],
            subject_scope=scope,
        )

    def execute_steps(self) -> None:
        starter_path = Path(self.starter_preparation["project_path"])
        pending_project_save = False
        steps = self.acceptance["steps"]
        for step_index, step in enumerate(steps):
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
            if should_flush_pending_save_before_step(pending_project_save, step):
                step_record["pre_scientific_save"] = self.save_project(
                    before_snapshot, self.args.timeouts["io"]
                )
                pending_project_save = False
                before_snapshot = self.wait_snapshot(
                    lambda snapshot: self.target_ready(snapshot, step, scope),
                    timeout,
                    f"reacquired target after pre-scientific save for step '{step['id']}'",
                    after_generation=step_record["pre_scientific_save"]["generation"],
                )
                target = self.item_for(
                    before_snapshot,
                    step["target"],
                    window_id=step["window"],
                    subject_scope=scope,
                )
                if target is None:
                    raise AcceptanceFailure(
                        "harness_gap",
                        f"Target vanished after pre-scientific save for step '{step['id']}'",
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
            step_record["x11_event"] = self.emit_x11(
                step, target, before_snapshot, timeout
            )
            after_snapshot = self.wait_after_interaction(
                step, scope, before_snapshot["generation"], timeout
            )
            step_record["after_generation"] = after_snapshot["generation"]
            pending_project_save = pending_project_save or bool(
                step.get("persists_project_state")
            )
            should_save = pending_project_save and (
                step.get("scientific_effect") or step_index == len(steps) - 1
            )
            if should_save:
                step_record["save"] = self.save_project(
                    after_snapshot, self.args.timeouts["io"]
                )
                pending_project_save = False
                after_snapshot = self.wait_snapshot(
                    lambda snapshot: snapshot.get("generation", 0)
                    >= step_record["save"]["generation"],
                    self.args.timeouts["instant"],
                    "post-save semantic snapshot",
                )
            elif pending_project_save:
                step_record["save"] = {
                    "event": "deferred_until_scientific_checkpoint",
                    "generation": after_snapshot["generation"],
                }
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
            if self.args.scrot is not None and self.gui_process is not None:
                failure_path = self.chapter_dir / "failure-diagnostic.raw.png"
                completed = subprocess.run(
                    [str(self.args.scrot), str(failure_path)],
                    cwd=self.repo_root,
                    env=self.process_environment,
                    timeout=20,
                    check=False,
                )
                if completed.returncode == 0 and failure_path.is_file():
                    self.ledger["failure_diagnostic"] = {
                        "role": "diagnostic_only",
                        "path": str(failure_path),
                        "sha256": sha256_file(failure_path),
                    }
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
        "xwininfo": tool_version([str(args.xwininfo), "-version"]),
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
        "network_namespace": {
            "current": (
                os.readlink("/proc/self/ns/net")
                if args.network_enforcement == "linux_network_namespace"
                else None
            ),
            "parent": args.parent_network_namespace,
        },
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
    parser.add_argument("--xwininfo")
    parser.add_argument("--scrot")
    parser.add_argument(
        "--parent-network-namespace",
        help=(
            "Parent network namespace identity captured with "
            "readlink /proc/self/ns/net before entering unshare"
        ),
    )
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
    parsed.xwininfo = resolve_tool(parsed.xwininfo, "xwininfo", required=True)
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
            except OSError as error:
                raise AcceptanceFailure(
                    "harness_gap", f"Could not verify Linux network namespace: {error}"
                ) from error
            validate_parent_network_namespace(
                own_namespace, args.parent_network_namespace
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
