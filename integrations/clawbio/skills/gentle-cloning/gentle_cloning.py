#!/usr/bin/env python3
"""ClawBio/OpenClaw wrapper skill for deterministic GENtle CLI execution."""

from __future__ import annotations

import argparse
import base64
import dataclasses
import datetime as dt
import hashlib
import html
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
from typing import Any
import xml.etree.ElementTree as ET

REQUEST_SCHEMA = "gentle.clawbio_skill_request.v1"
RESULT_SCHEMA = "gentle.clawbio_skill_result.v1"
SVG_DIMENSION_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)")


class SkillError(RuntimeError):
    """Base skill error with deterministic message formatting."""


@dataclasses.dataclass
class Request:
    mode: str
    timeout_secs: int = 180
    state_path: str | None = None
    raw_args: list[str] | None = None
    shell_line: str | None = None
    operation: Any = None
    workflow: Any = None
    workflow_path: str | None = None
    expected_artifacts: list[str] | None = None
    ensure_reference_prepared: Any = None


@dataclasses.dataclass
class EnsureReferencePrepared:
    genome_id: str
    catalog_path: str | None = None
    cache_dir: str | None = None
    status_timeout_secs: int = 300
    prepare_timeout_secs: int = 7200


@dataclasses.dataclass
class CliResolution:
    argv_prefix: list[str]
    cwd: str | None
    label: str


def _format_command_text(command: list[str] | None) -> str:
    if not command:
        return "(none)"
    return " ".join(shlex.quote(v) for v in command)


def _now_utc_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as e:
        raise SkillError(f"request JSON file '{path}' does not exist") from e
    except json.JSONDecodeError as e:
        raise SkillError(f"invalid JSON in '{path}': {e}") from e


def _repo_root_candidates(script_path: Path) -> list[Path]:
    candidates: list[Path] = [Path.cwd()]
    candidates.extend(script_path.resolve().parents)
    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def _find_repo_root(script_path: Path) -> Path | None:
    for candidate in _repo_root_candidates(script_path):
        if (candidate / "Cargo.toml").exists() and (
            candidate / "src" / "bin" / "gentle_cli.rs"
        ).exists():
            return candidate
    return None


def _find_repo_root_from_path(path: Path) -> Path | None:
    candidate = path.resolve()
    if candidate.is_file():
        candidate = candidate.parent
    while True:
        if (candidate / "Cargo.toml").exists() and (
            candidate / "src" / "bin" / "gentle_cli.rs"
        ).exists():
            return candidate
        if candidate.parent == candidate:
            return None
        candidate = candidate.parent


def _request_path_candidates(raw_path: str, script_path: Path) -> list[Path]:
    path = Path(raw_path)
    if path.is_absolute():
        return [path]

    candidates: list[Path] = [Path.cwd() / path]
    for base in script_path.resolve().parents:
        candidates.append(base / path)
    configured_repo_root = os.environ.get("GENTLE_REPO_ROOT", "").strip()
    if configured_repo_root:
        candidates.append(Path(configured_repo_root) / path)
    repo_root = _find_repo_root(script_path)
    if repo_root is not None:
        candidates.append(repo_root / path)

    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key not in seen:
            seen.add(key)
            unique.append(candidate)
    return unique


def _resolve_request_path(raw_path: str, script_path: Path) -> str:
    trimmed = raw_path.strip()
    if not trimmed:
        return raw_path
    for candidate in _request_path_candidates(trimmed, script_path):
        if candidate.exists():
            return str(candidate.resolve())
    return raw_path


def _resolve_existing_request_file(raw_path: str, script_path: Path) -> Path:
    trimmed = raw_path.strip()
    candidates = _request_path_candidates(trimmed, script_path)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    tried = ", ".join(f"`{candidate}`" for candidate in candidates)
    raise SkillError(
        f"request JSON file '{raw_path}' was not found. "
        f"Current cwd: `{Path.cwd()}`. Tried: {tried}"
    )


def _resolve_cli(explicit: str | None, script_path: Path) -> CliResolution:
    if explicit:
        return CliResolution(
            argv_prefix=shlex.split(explicit),
            cwd=None,
            label=f"explicit --gentle-cli: {explicit}",
        )
    env_cmd = os.environ.get("GENTLE_CLI_CMD", "").strip()
    if env_cmd:
        return CliResolution(
            argv_prefix=shlex.split(env_cmd),
            cwd=None,
            label=f"GENTLE_CLI_CMD: {env_cmd}",
        )
    path_hit = shutil.which("gentle_cli")
    if path_hit:
        return CliResolution(
            argv_prefix=[path_hit],
            cwd=None,
            label=f"PATH gentle_cli: {path_hit}",
        )
    repo_root = _find_repo_root(script_path)
    if repo_root is not None and shutil.which("cargo"):
        return CliResolution(
            argv_prefix=["cargo", "run", "--quiet", "--bin", "gentle_cli", "--"],
            cwd=str(repo_root),
            label=f"cargo run fallback in repo: {repo_root}",
        )
    raise SkillError(
        "Could not resolve gentle_cli executable. "
        "Recommended: set GENTLE_CLI_CMD to the included "
        "'./gentle_local_checkout_cli.sh' launcher (typically with "
        "GENTLE_REPO_ROOT=/absolute/path/to/GENtle), or to a Docker/OCI or "
        "Apptainer/Singularity command such as "
        "'docker run --rm -i -v \"$PWD\":/work -w /work "
        "ghcr.io/smoe/gentle_rs:cli' or "
        "'./gentle_apptainer_cli.sh /absolute/path/to/gentle.sif'. "
        "Alternatives: use --gentle-cli, install gentle_cli on PATH, or run "
        "from a local repo with cargo available."
    )


def _coerce_request(payload: dict[str, Any]) -> Request:
    schema = payload.get("schema")
    if schema != REQUEST_SCHEMA:
        raise SkillError(
            f"unsupported request schema '{schema}', expected '{REQUEST_SCHEMA}'"
        )
    mode = str(payload.get("mode", "")).strip()
    if not mode:
        raise SkillError("request missing required field 'mode'")
    timeout_secs = int(payload.get("timeout_secs", 180))
    if timeout_secs <= 0:
        raise SkillError("timeout_secs must be > 0")
    request = Request(
        mode=mode,
        timeout_secs=timeout_secs,
        state_path=payload.get("state_path"),
        raw_args=payload.get("raw_args"),
        shell_line=payload.get("shell_line"),
        operation=payload.get("operation"),
        workflow=payload.get("workflow"),
        workflow_path=payload.get("workflow_path"),
        expected_artifacts=payload.get("expected_artifacts"),
        ensure_reference_prepared=payload.get("ensure_reference_prepared"),
    )
    if request.expected_artifacts is not None:
        if not isinstance(request.expected_artifacts, list):
            raise SkillError("expected_artifacts must be a string array when present")
        if not all(isinstance(v, str) and v.strip() for v in request.expected_artifacts):
            raise SkillError("expected_artifacts must contain non-empty strings")
    if request.mode == "raw":
        if not isinstance(request.raw_args, list) or not request.raw_args:
            raise SkillError("mode=raw requires non-empty string array 'raw_args'")
        if not all(isinstance(v, str) and v for v in request.raw_args):
            raise SkillError("mode=raw 'raw_args' must contain non-empty strings")
    elif request.mode == "shell":
        if not isinstance(request.shell_line, str) or not request.shell_line.strip():
            raise SkillError("mode=shell requires non-empty string field 'shell_line'")
    elif request.mode == "op":
        if request.operation is None:
            raise SkillError("mode=op requires field 'operation'")
    elif request.mode == "workflow":
        if request.workflow is None and not request.workflow_path:
            raise SkillError("mode=workflow requires 'workflow' or 'workflow_path'")
    elif request.mode not in ("capabilities", "state-summary"):
        raise SkillError(
            "unsupported mode. Use one of: capabilities, state-summary, shell, "
            "op, workflow, raw"
        )
    if request.ensure_reference_prepared is not None:
        ensure = request.ensure_reference_prepared
        if not isinstance(ensure, dict):
            raise SkillError(
                "ensure_reference_prepared must be an object when present"
            )
        genome_id = str(ensure.get("genome_id", "")).strip()
        if not genome_id:
            raise SkillError(
                "ensure_reference_prepared requires non-empty string field 'genome_id'"
            )
        status_timeout_secs = int(ensure.get("status_timeout_secs", 300))
        prepare_timeout_secs = int(ensure.get("prepare_timeout_secs", 7200))
        if status_timeout_secs <= 0:
            raise SkillError(
                "ensure_reference_prepared.status_timeout_secs must be > 0"
            )
        if prepare_timeout_secs <= 0:
            raise SkillError(
                "ensure_reference_prepared.prepare_timeout_secs must be > 0"
            )
        for optional_path_field in ("catalog_path", "cache_dir"):
            value = ensure.get(optional_path_field)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise SkillError(
                    f"ensure_reference_prepared.{optional_path_field} must be a non-empty string when present"
                )
        request.ensure_reference_prepared = EnsureReferencePrepared(
            genome_id=genome_id,
            catalog_path=ensure.get("catalog_path"),
            cache_dir=ensure.get("cache_dir"),
            status_timeout_secs=status_timeout_secs,
            prepare_timeout_secs=prepare_timeout_secs,
        )
    return request


def _json_arg(value: Any) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=True, separators=(",", ":"))


def _preview_text(text: str, *, max_lines: int = 6, max_chars: int = 600) -> str | None:
    trimmed = text.strip()
    if not trimmed:
        return None
    lines = trimmed.splitlines()
    preview = "\n".join(lines[:max_lines])
    if len(preview) > max_chars:
        preview = preview[: max_chars - 3].rstrip() + "..."
    elif len(lines) > max_lines or len(trimmed) > len(preview):
        preview = preview.rstrip() + "..."
    return preview


def _infer_command_compatibility_note(
    command: list[str] | None,
    stderr_text: str,
) -> str | None:
    if not command:
        return None

    stderr = stderr_text.strip()
    if not stderr:
        return None

    requested = ""
    if len(command) >= 3 and command[1] == "shell":
        requested = command[2]
    elif len(command) >= 3:
        requested = " ".join(command[1:3])
    elif len(command) >= 2:
        requested = command[1]

    requested = requested.strip()
    if not requested:
        return None

    command_unknown = (
        "Unknown shell command" in stderr
        or "unrecognized subcommand" in stderr
        or "unexpected argument" in stderr
        or "Found argument" in stderr
    )
    if not command_unknown:
        return None

    if requested.startswith(("services ", "resources ", "helpers status", "genomes status")):
        return (
            "This GENtle CLI looks older than the current ClawBio skill scaffold and "
            "may not support the requested status route yet. Update the installed "
            "`gentle_cli` on PATH, or set `GENTLE_CLI_CMD` to "
            "`gentle_local_checkout_cli.sh` for an updated GENtle checkout."
        )
    return (
        "This GENtle CLI may be older than the current ClawBio skill scaffold and "
        "may not support the requested command yet. Update the installed "
        "`gentle_cli` on PATH, or point `GENTLE_CLI_CMD` at the updated GENtle "
        "checkout launcher."
    )


def _build_failure_summary(
    *,
    stage: str,
    step: dict[str, Any] | None,
    execution_cwd: Path | None,
    note: str | None = None,
) -> dict[str, Any]:
    command = step.get("command") if isinstance(step, dict) else None
    effective_note = note
    if effective_note is None and isinstance(step, dict):
        effective_note = _infer_command_compatibility_note(
            command,
            step.get("stderr", ""),
        )
    summary = {
        "stage": stage,
        "note": effective_note,
        "command": command,
        "command_text": _format_command_text(command),
        "execution_cwd": (str(execution_cwd) if execution_cwd is not None else None),
        "exit_code": (step.get("exit_code") if isinstance(step, dict) else None),
        "stderr_preview": (
            _preview_text(step.get("stderr", "")) if isinstance(step, dict) else None
        ),
        "stdout_preview": (
            _preview_text(step.get("stdout", "")) if isinstance(step, dict) else None
        ),
    }
    return summary


def _build_failure_message(
    *,
    headline: str,
    failure_summary: dict[str, Any] | None,
) -> str:
    if not failure_summary:
        return headline

    parts = [headline]
    note = failure_summary.get("note")
    if isinstance(note, str) and note.strip():
        parts.append(note.strip())
    command_text = failure_summary.get("command_text")
    if isinstance(command_text, str) and command_text and command_text != "(none)":
        parts.append(f"Command: `{command_text}`.")
    execution_cwd = failure_summary.get("execution_cwd")
    if isinstance(execution_cwd, str) and execution_cwd:
        parts.append(f"Cwd: `{execution_cwd}`.")
    exit_code = failure_summary.get("exit_code")
    if exit_code is not None:
        parts.append(f"Exit code: `{exit_code}`.")
    stderr_preview = failure_summary.get("stderr_preview")
    stdout_preview = failure_summary.get("stdout_preview")
    if isinstance(stderr_preview, str) and stderr_preview:
        parts.append(f"Stderr preview: `{stderr_preview}`.")
    elif isinstance(stdout_preview, str) and stdout_preview:
        parts.append(f"Stdout preview: `{stdout_preview}`.")
    parts.append("Full stdout/stderr are recorded in `report.md` and `result.json`.")
    return " ".join(parts)


def _resolve_execution_cwd(
    request: Request,
    resolution: CliResolution | None,
    script_path: Path,
) -> Path:
    if request.mode == "workflow" and request.workflow_path:
        resolved_workflow_path = Path(_resolve_request_path(request.workflow_path, script_path))
        if resolved_workflow_path.exists():
            repo_root = _find_repo_root_from_path(resolved_workflow_path)
            if repo_root is not None:
                return repo_root
            return resolved_workflow_path.parent
    if resolution and resolution.cwd:
        return Path(resolution.cwd)
    return Path.cwd()


def _resolve_expected_artifact(raw_path: str, execution_cwd: Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return execution_cwd / path


def _safe_artifact_destination(raw_path: str) -> Path:
    normalized = raw_path.replace("\\", "/")
    parts: list[str] = []
    for part in Path(normalized).parts:
        if part in ("", ".", "..", "/", "\\"):
            continue
        if part.endswith(":"):
            continue
        parts.append(part)
    if not parts:
        return Path("artifact")
    return Path(*parts)


def _copy_collected_artifacts(
    request: Request,
    output_dir: Path,
    execution_cwd: Path,
) -> list[dict[str, str]]:
    collected: list[dict[str, str]] = []
    if not request.expected_artifacts:
        return collected

    generated_dir = output_dir / "generated"
    for raw_path in request.expected_artifacts:
        source_path = _resolve_expected_artifact(raw_path, execution_cwd)
        if not source_path.exists() or not source_path.is_file():
            continue

        destination = generated_dir / _safe_artifact_destination(raw_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination)
        collected.append(
            {
                "declared_path": raw_path,
                "source_path": str(source_path.resolve()),
                "copied_path": str(destination.resolve()),
            }
        )
    return collected


def _parse_svg_dimension(value: str | None) -> float | None:
    if value is None:
        return None
    match = SVG_DIMENSION_RE.match(value)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def _storyboard_caption_for_artifact(raw_path: str, preferred: dict[str, Any] | None) -> str:
    if preferred is not None and isinstance(preferred.get("caption"), str):
        caption = preferred["caption"].strip()
        if caption:
            return caption
    stem = Path(raw_path).stem.replace("_", " ").replace(".", " ")
    return " ".join(part for part in stem.split() if part)


def _storyboard_panel_title(raw_path: str, preferred: dict[str, Any] | None) -> str:
    lowered = raw_path.lower()
    if "context" in lowered:
        return "Genomic context"
    if "reference" in lowered and ("reporter" in lowered or "construct" in lowered):
        return "Reference allele reporter"
    if "alternate" in lowered and ("reporter" in lowered or "construct" in lowered):
        return "Alternate allele reporter"
    return _storyboard_caption_for_artifact(raw_path, preferred)


def _storyboard_panel_priority(raw_path: str) -> tuple[int, str]:
    lowered = raw_path.lower()
    if "context" in lowered:
        return (0, lowered)
    if "reference" in lowered:
        return (1, lowered)
    if "alternate" in lowered:
        return (2, lowered)
    return (3, lowered)


def _storyboard_panel_record(
    copied_artifact: dict[str, str],
    preferred: dict[str, Any] | None,
) -> dict[str, Any] | None:
    copied_path = Path(copied_artifact["copied_path"])
    if copied_path.suffix.lower() != ".svg":
        return None
    try:
        text = copied_path.read_text(encoding="utf-8")
    except OSError:
        return None
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return None

    width = None
    height = None
    view_box = root.get("viewBox")
    if view_box:
        parts = view_box.replace(",", " ").split()
        if len(parts) == 4:
            try:
                width = float(parts[2])
                height = float(parts[3])
            except ValueError:
                width = None
                height = None
    if width is None or height is None or width <= 0.0 or height <= 0.0:
        width = _parse_svg_dimension(root.get("width")) or 800.0
        height = _parse_svg_dimension(root.get("height")) or 600.0
    if width <= 0.0 or height <= 0.0:
        return None

    return {
        "declared_path": copied_artifact["declared_path"],
        "copied_path": copied_artifact["copied_path"],
        "artifact_id": (
            preferred.get("artifact_id")
            if preferred is not None and isinstance(preferred.get("artifact_id"), str)
            else Path(copied_artifact["declared_path"]).stem
        ),
        "caption": _storyboard_caption_for_artifact(copied_artifact["declared_path"], preferred),
        "panel_title": _storyboard_panel_title(copied_artifact["declared_path"], preferred),
        "priority": _storyboard_panel_priority(copied_artifact["declared_path"]),
        "width": width,
        "height": height,
        "data_uri": (
            "data:image/svg+xml;base64,"
            + base64.b64encode(text.encode("utf-8")).decode("ascii")
        ),
    }


def _storyboard_frame_rects(
    panel_count: int,
    variant_layout: bool,
) -> list[tuple[float, float, float, float]]:
    if variant_layout and panel_count >= 3:
        return [
            (64.0, 172.0, 900.0, 860.0),
            (1000.0, 172.0, 536.0, 404.0),
            (1000.0, 628.0, 536.0, 404.0),
        ]
    columns = 2 if panel_count > 1 else 1
    rows = max(1, math.ceil(panel_count / columns))
    margin = 64.0
    gutter = 36.0
    top = 172.0
    canvas_width = 1600.0
    canvas_height = 1100.0
    frame_width = (canvas_width - margin * 2 - gutter * (columns - 1)) / columns
    frame_height = (canvas_height - top - 72.0 - gutter * (rows - 1)) / rows
    rects: list[tuple[float, float, float, float]] = []
    for idx in range(panel_count):
        row = idx // columns
        col = idx % columns
        rects.append(
            (
                margin + col * (frame_width + gutter),
                top + row * (frame_height + gutter),
                frame_width,
                frame_height,
            )
        )
    return rects


def _write_storyboard_svg(
    title: str,
    subtitle: str,
    panels: list[dict[str, Any]],
    output_path: Path,
) -> None:
    variant_layout = (
        len(panels) == 3
        and "context" in panels[0]["declared_path"].lower()
        and any("reference" in panel["declared_path"].lower() for panel in panels)
        and any("alternate" in panel["declared_path"].lower() for panel in panels)
    )
    frame_rects = _storyboard_frame_rects(len(panels), variant_layout)

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1100" viewBox="0 0 1600 1100" role="img" aria-labelledby="title subtitle">',
        "  <defs>",
        '    <linearGradient id="bg" x1="0%" y1="0%" x2="100%" y2="100%">',
        '      <stop offset="0%" stop-color="#f8fafc" />',
        '      <stop offset="100%" stop-color="#e2e8f0" />',
        "    </linearGradient>",
        '    <filter id="card-shadow" x="-10%" y="-10%" width="130%" height="140%">',
        '      <feDropShadow dx="0" dy="10" stdDeviation="14" flood-color="#0f172a" flood-opacity="0.12" />',
        "    </filter>",
        "  </defs>",
        '  <rect width="1600" height="1100" fill="url(#bg)" />',
        '  <text id="title" x="64" y="78" font-family="Inter,Segoe UI,sans-serif" font-size="34" font-weight="700" fill="#0f172a">'
        + html.escape(title)
        + "</text>",
        '  <text id="subtitle" x="64" y="114" font-family="Inter,Segoe UI,sans-serif" font-size="16" fill="#334155">'
        + html.escape(subtitle)
        + "</text>",
        '  <text x="64" y="144" font-family="Inter,Segoe UI,sans-serif" font-size="13" fill="#475569">Deterministic GENtle figures bundled for a ClawBio/OpenClaw variant follow-up reply.</text>',
    ]

    for panel, (x, y, width, height) in zip(panels, frame_rects):
        header_height = 76.0
        image_pad = 20.0
        image_box_x = x + image_pad
        image_box_y = y + header_height
        image_box_width = width - image_pad * 2
        image_box_height = height - header_height - image_pad
        scale = min(
            image_box_width / panel["width"],
            image_box_height / panel["height"],
        )
        render_width = panel["width"] * scale
        render_height = panel["height"] * scale
        render_x = image_box_x + (image_box_width - render_width) / 2.0
        render_y = image_box_y + (image_box_height - render_height) / 2.0
        lines.extend(
            [
                f'  <g filter="url(#card-shadow)">',
                f'    <rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{height:.1f}" rx="28" fill="#ffffff" />',
                f'  </g>',
                f'  <text x="{x + 24.0:.1f}" y="{y + 32.0:.1f}" font-family="Inter,Segoe UI,sans-serif" font-size="20" font-weight="700" fill="#0f172a">{html.escape(panel["panel_title"])}</text>',
                f'  <text x="{x + 24.0:.1f}" y="{y + 56.0:.1f}" font-family="Inter,Segoe UI,sans-serif" font-size="13" fill="#64748b">{html.escape(panel["caption"])}</text>',
                f'  <image href="{panel["data_uri"]}" x="{render_x:.1f}" y="{render_y:.1f}" width="{render_width:.1f}" height="{render_height:.1f}" preserveAspectRatio="xMidYMid meet" />',
            ]
        )

    lines.extend(
        [
            '  <text x="64" y="1070" font-family="Inter,Segoe UI,sans-serif" font-size="12" fill="#64748b">Typical answer shape: locus context, assayable promoter fragment logic, and engineered allele-specific reporter constructs for synthetic-biology follow-up.</text>',
            "</svg>",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _augment_artifacts_with_storyboard(
    request: Request,
    output_dir: Path,
    collected_artifacts: list[dict[str, str]],
    preferred_artifacts: list[dict[str, Any]] | None,
) -> tuple[list[dict[str, str]], list[dict[str, Any]] | None]:
    preferred_by_path = {
        artifact["path"]: artifact
        for artifact in (preferred_artifacts or [])
        if isinstance(artifact, dict) and isinstance(artifact.get("path"), str)
    }
    panels = []
    for artifact in collected_artifacts:
        preferred = preferred_by_path.get(artifact["declared_path"])
        panel = _storyboard_panel_record(artifact, preferred)
        if panel is not None:
            panels.append(panel)
    panels.sort(key=lambda panel: panel["priority"])
    if len(panels) < 2:
        return collected_artifacts, preferred_artifacts

    workflow_hint = " ".join(
        value
        for value in (
            request.workflow_path or "",
            request.shell_line or "",
        )
        if value
    ).lower()
    variant_story = (
        any("context" in panel["declared_path"].lower() for panel in panels)
        and any("reference" in panel["declared_path"].lower() for panel in panels)
        and any("alternate" in panel["declared_path"].lower() for panel in panels)
    ) or ("variant" in workflow_hint or "luciferase" in workflow_hint)
    if variant_story:
        title = "Variant-to-Synthetic-Biology assay storyboard"
        subtitle = (
            "GENtle bridges genomic variant context to engineered allele-specific "
            "reporter constructs for reproducible functional follow-up."
        )
    else:
        title = "GENtle graphical storyboard"
        subtitle = "One shareable figure assembled from the deterministic graphics in this run."

    storyboard_path = output_dir / "generated" / "clawbio_storyboard.svg"
    storyboard_path.parent.mkdir(parents=True, exist_ok=True)
    _write_storyboard_svg(title, subtitle, panels[:4], storyboard_path)

    storyboard_artifact = {
        "declared_path": "generated/clawbio_storyboard.svg",
        "source_path": str(storyboard_path.resolve()),
        "copied_path": str(storyboard_path.resolve()),
    }
    updated_collected = [*collected_artifacts, storyboard_artifact]

    storyboard_entry = {
        "artifact_id": "clawbio_storyboard_svg",
        "path": "generated/clawbio_storyboard.svg",
        "caption": title,
        "recommended_use": "best_first_figure",
        "presentation_rank": 0,
        "is_best_first_artifact": True,
    }
    merged_preferred: list[dict[str, Any]] = [storyboard_entry]
    seen_paths = {storyboard_entry["path"]}

    source_entries = preferred_artifacts or []
    if not source_entries:
        source_entries = [
            {
                "artifact_id": panel["artifact_id"],
                "path": panel["declared_path"],
                "caption": panel["caption"],
                "recommended_use": "supporting_figure",
                "presentation_rank": idx + 1,
                "is_best_first_artifact": False,
            }
            for idx, panel in enumerate(panels)
        ]

    for idx, artifact in enumerate(source_entries, start=1):
        if not isinstance(artifact, dict):
            continue
        path = artifact.get("path")
        if not isinstance(path, str) or path in seen_paths:
            continue
        updated = dict(artifact)
        updated["is_best_first_artifact"] = False
        if not isinstance(updated.get("presentation_rank"), int):
            updated["presentation_rank"] = idx
        elif updated["presentation_rank"] <= 0:
            updated["presentation_rank"] = idx
        seen_paths.add(path)
        merged_preferred.append(updated)

    return updated_collected, merged_preferred


def _build_cli_args(request: Request, script_path: Path) -> list[str]:
    args: list[str] = []
    if request.state_path:
        args.extend(["--state", request.state_path])
    if request.mode == "capabilities":
        args.append("capabilities")
    elif request.mode == "state-summary":
        args.append("state-summary")
    elif request.mode == "shell":
        args.extend(["shell", request.shell_line.strip()])
    elif request.mode == "op":
        args.extend(["op", _json_arg(request.operation)])
    elif request.mode == "workflow":
        if request.workflow_path:
            resolved_workflow_path = _resolve_request_path(request.workflow_path, script_path)
            args.extend(["workflow", f"@{resolved_workflow_path}"])
        else:
            args.extend(["workflow", _json_arg(request.workflow)])
    elif request.mode == "raw":
        args.extend(request.raw_args or [])
    else:
        raise SkillError(f"Unsupported mode '{request.mode}'")
    return args


def _build_reference_status_args(reference: EnsureReferencePrepared) -> list[str]:
    args = ["genomes", "status", reference.genome_id]
    if reference.catalog_path:
        args.extend(["--catalog", reference.catalog_path])
    if reference.cache_dir:
        args.extend(["--cache-dir", reference.cache_dir])
    return args


def _build_reference_prepare_args(reference: EnsureReferencePrepared) -> list[str]:
    args = [
        "genomes",
        "prepare",
        reference.genome_id,
        "--timeout-secs",
        str(reference.prepare_timeout_secs),
    ]
    if reference.catalog_path:
        args.extend(["--catalog", reference.catalog_path])
    if reference.cache_dir:
        args.extend(["--cache-dir", reference.cache_dir])
    return args


def _run_cli_command(
    resolution: CliResolution,
    cli_args: list[str],
    execution_cwd: Path,
    timeout_secs: int,
) -> tuple[subprocess.CompletedProcess[str], dict[str, Any]]:
    command = resolution.argv_prefix + cli_args
    started_utc = _now_utc_iso()
    run_result = subprocess.run(
        command,
        cwd=execution_cwd,
        capture_output=True,
        text=True,
        timeout=timeout_secs,
        check=False,
    )
    ended_utc = _now_utc_iso()
    step = {
        "command": command,
        "started_utc": started_utc,
        "ended_utc": ended_utc,
        "exit_code": run_result.returncode,
        "stdout": run_result.stdout,
        "stderr": run_result.stderr,
        "status": ("ok" if run_result.returncode == 0 else "command_failed"),
    }
    return run_result, step


def _reference_preflight_record(reference: EnsureReferencePrepared) -> dict[str, Any]:
    return {
        "genome_id": reference.genome_id,
        "catalog_path": reference.catalog_path,
        "cache_dir": reference.cache_dir,
        "status_timeout_secs": reference.status_timeout_secs,
        "prepare_timeout_secs": reference.prepare_timeout_secs,
        "prepared_before": None,
        "prepared_after": None,
        "prepare_attempted": False,
        "status_before": None,
        "status_after": None,
        "steps": [],
        "last_failure": None,
        "status": "pending",
    }


def _parse_json_stdout(
    stdout: str,
    context: str,
    *,
    failure_summary: dict[str, Any] | None = None,
) -> Any:
    try:
        return json.loads(stdout)
    except json.JSONDecodeError as e:
        detail = _build_failure_message(
            headline=f"{context} did not return valid JSON: {e}",
            failure_summary=failure_summary,
        )
        raise SkillError(detail) from e


def _ensure_reference_prepared(
    reference: EnsureReferencePrepared,
    resolution: CliResolution,
    execution_cwd: Path,
    record: dict[str, Any],
) -> None:
    status_args = _build_reference_status_args(reference)
    status_result, status_step = _run_cli_command(
        resolution,
        status_args,
        execution_cwd,
        reference.status_timeout_secs,
    )
    record["steps"].append(status_step)
    if status_result.returncode != 0:
        record["status"] = "failed"
        failure_summary = _build_failure_summary(
            stage="reference_status_preflight",
            step=status_step,
            execution_cwd=execution_cwd,
        )
        record["last_failure"] = failure_summary
        raise SkillError(
            _build_failure_message(
                headline=(
                    "reference-status preflight failed while checking whether the requested "
                    "reference is already prepared."
                ),
                failure_summary=failure_summary,
            )
        )
    status_payload = _parse_json_stdout(
        status_result.stdout,
        f"genomes status {reference.genome_id}",
        failure_summary=_build_failure_summary(
            stage="reference_status_preflight",
            step=status_step,
            execution_cwd=execution_cwd,
        ),
    )
    record["status_before"] = status_payload
    record["prepared_before"] = bool(status_payload.get("prepared"))
    if record["prepared_before"]:
        record["prepared_after"] = True
        record["status_after"] = status_payload
        record["status"] = "already_prepared"
        return

    prepare_args = _build_reference_prepare_args(reference)
    prepare_result, prepare_step = _run_cli_command(
        resolution,
        prepare_args,
        execution_cwd,
        reference.prepare_timeout_secs,
    )
    record["steps"].append(prepare_step)
    record["prepare_attempted"] = True
    if prepare_result.returncode != 0:
        record["status"] = "failed"
        failure_summary = _build_failure_summary(
            stage="reference_prepare_preflight",
            step=prepare_step,
            execution_cwd=execution_cwd,
        )
        record["last_failure"] = failure_summary
        raise SkillError(
            _build_failure_message(
                headline=(
                    "reference-prepare preflight failed while trying to prepare the requested "
                    "reference locally."
                ),
                failure_summary=failure_summary,
            )
        )

    status_after_result, status_after_step = _run_cli_command(
        resolution,
        status_args,
        execution_cwd,
        reference.status_timeout_secs,
    )
    record["steps"].append(status_after_step)
    if status_after_result.returncode != 0:
        record["status"] = "failed"
        failure_summary = _build_failure_summary(
            stage="post_prepare_reference_status_preflight",
            step=status_after_step,
            execution_cwd=execution_cwd,
        )
        record["last_failure"] = failure_summary
        raise SkillError(
            _build_failure_message(
                headline=(
                    "post-prepare reference-status preflight failed after the prepare step "
                    "completed."
                ),
                failure_summary=failure_summary,
            )
        )
    status_after_payload = _parse_json_stdout(
        status_after_result.stdout,
        f"post-prepare genomes status {reference.genome_id}",
        failure_summary=_build_failure_summary(
            stage="post_prepare_reference_status_preflight",
            step=status_after_step,
            execution_cwd=execution_cwd,
        ),
    )
    record["status_after"] = status_after_payload
    record["prepared_after"] = bool(status_after_payload.get("prepared"))
    if not record["prepared_after"]:
        record["status"] = "failed"
        raise SkillError(
            "reference-prepare completed but the requested reference still is not reported as prepared"
        )
    record["status"] = "prepared_during_run"


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_repro_environment(path: Path) -> None:
    lines = [
        "name: gentle-clawbio-skill",
        "channels:",
        "  - defaults",
        "dependencies:",
        f"  - python={sys.version_info.major}.{sys.version_info.minor}",
        "variables:",
        f'  PYTHON_EXECUTABLE: "{sys.executable}"',
        f'  PLATFORM: "{platform.platform()}"',
        f'  PYTHON_VERSION: "{platform.python_version()}"',
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parse_stdout_json(stdout: str) -> Any | None:
    trimmed = stdout.strip()
    if not trimmed:
        return None
    try:
        return json.loads(trimmed)
    except json.JSONDecodeError:
        return None


def _summary_lines_from_sequence_context_view(candidate: Any) -> list[str] | None:
    if not isinstance(candidate, dict):
        return None
    if candidate.get("schema") != "gentle.sequence_context_view.v1":
        return None
    raw_lines = candidate.get("summary_lines")
    if not isinstance(raw_lines, list):
        return None
    lines = [line.strip() for line in raw_lines if isinstance(line, str) and line.strip()]
    return lines or None


def _extract_sequence_context_bundle(stdout_json: Any) -> dict[str, Any] | None:
    if not isinstance(stdout_json, dict):
        return None
    if stdout_json.get("schema") == "gentle.sequence_context_bundle.v1":
        return stdout_json
    candidate = stdout_json.get("sequence_context_bundle")
    if isinstance(candidate, dict) and candidate.get("schema") == "gentle.sequence_context_bundle.v1":
        return candidate
    return None


def _extract_preferred_artifacts(stdout_json: Any) -> list[dict[str, Any]] | None:
    bundle = _extract_sequence_context_bundle(stdout_json)
    if not isinstance(bundle, dict):
        return None
    raw_artifacts = bundle.get("artifacts")
    if not isinstance(raw_artifacts, list):
        return None
    artifacts = [
        artifact
        for artifact in raw_artifacts
        if isinstance(artifact, dict) and isinstance(artifact.get("path"), str)
    ]
    if not artifacts:
        return None
    artifacts.sort(
        key=lambda artifact: (
            int(artifact.get("presentation_rank", 10**9))
            if isinstance(artifact.get("presentation_rank"), int)
            else 10**9,
            str(artifact.get("artifact_id", "")),
        )
    )
    return artifacts


def _extract_chat_summary_lines(stdout_json: Any) -> list[str] | None:
    if not isinstance(stdout_json, dict):
        return None
    bundle = _extract_sequence_context_bundle(stdout_json)
    for candidate in (
        stdout_json,
        stdout_json.get("sequence_context_view"),
        bundle,
        bundle.get("sequence_context_view") if isinstance(bundle, dict) else None,
    ):
        lines = _summary_lines_from_sequence_context_view(candidate)
        if lines:
            return lines
    return None


def _write_report(
    path: Path,
    request: Request | None,
    resolution: CliResolution | None,
    execution_cwd: Path,
    command: list[str] | None,
    run_result: subprocess.CompletedProcess[str] | None,
    stdout_json: Any | None,
    chat_summary_lines: list[str] | None,
    collected_artifacts: list[dict[str, str]],
    reference_preflight: dict[str, Any] | None,
    started_utc: str,
    ended_utc: str,
    status: str,
    error_message: str | None,
    failure_summary: dict[str, Any] | None,
    preferred_artifacts: list[dict[str, Any]] | None,
) -> None:
    command_text = _format_command_text(command)
    stdout = run_result.stdout if run_result else ""
    stderr = run_result.stderr if run_result else ""
    exit_code = run_result.returncode if run_result else None
    lines = [
        "# GENtle ClawBio Skill Report",
        "",
        f"- Started (UTC): `{started_utc}`",
        f"- Ended (UTC): `{ended_utc}`",
        f"- Status: `{status}`",
        f"- Mode: `{request.mode if request is not None else 'unknown'}`",
    ]
    if request is not None and request.state_path:
        lines.append(f"- State path: `{request.state_path}`")
    if resolution is not None:
        lines.append(f"- Resolver: `{resolution.label}`")
    lines.append(f"- Execution cwd: `{execution_cwd}`")
    if exit_code is not None:
        lines.append(f"- Exit code: `{exit_code}`")
    if isinstance(stdout_json, dict) and isinstance(stdout_json.get("schema"), str):
        lines.append(f"- Parsed stdout JSON schema: `{stdout_json['schema']}`")
    if error_message:
        lines.append(f"- Error: `{error_message}`")
    if failure_summary:
        lines.append(
            f"- Failure stage: `{failure_summary.get('stage', 'unknown')}`"
        )
    if collected_artifacts:
        lines.append("- Collected artifacts:")
        for artifact in collected_artifacts:
            lines.append(
                f"  - `{artifact['declared_path']}` -> `{artifact['copied_path']}`"
            )
    if preferred_artifacts:
        best_first = next(
            (
                artifact
                for artifact in preferred_artifacts
                if artifact.get("is_best_first_artifact") is True
            ),
            preferred_artifacts[0],
        )
        lines.append(
            f"- Best first artifact: `{best_first.get('path', '(unknown)')}`"
        )
    if reference_preflight:
        lines.append(
            f"- Reference preflight: `{reference_preflight.get('status', 'unknown')}`"
        )
    if chat_summary_lines:
        lines.extend(["", "## Chat Summary", ""])
        lines.extend(f"- {line}" for line in chat_summary_lines)
    if preferred_artifacts:
        lines.extend(["", "## Preferred Artifacts", ""])
        for artifact in preferred_artifacts:
            marker = " (best first)" if artifact.get("is_best_first_artifact") is True else ""
            lines.append(
                f"- `{artifact.get('artifact_id', 'artifact')}`{marker}: "
                f"`{artifact.get('path', '(unknown)')}`"
            )
            if artifact.get("caption"):
                lines.append(f"  Caption: `{artifact['caption']}`")
            if artifact.get("recommended_use"):
                lines.append(f"  Recommended use: `{artifact['recommended_use']}`")
    if failure_summary:
        lines.extend(
            [
                "",
                "## Failure Summary",
                "",
                f"- Stage: `{failure_summary.get('stage', 'unknown')}`",
                f"- Command: `{failure_summary.get('command_text', '(none)')}`",
            ]
        )
        if failure_summary.get("execution_cwd"):
            lines.append(
                f"- Cwd: `{failure_summary['execution_cwd']}`"
            )
        if failure_summary.get("exit_code") is not None:
            lines.append(
                f"- Exit code: `{failure_summary['exit_code']}`"
            )
        if failure_summary.get("note"):
            lines.append(f"- Note: `{failure_summary['note']}`")
        if failure_summary.get("stderr_preview"):
            lines.extend(
                [
                    "",
                    "### Stderr Preview",
                    "",
                    "```text",
                    str(failure_summary["stderr_preview"]),
                    "```",
                ]
            )
        elif failure_summary.get("stdout_preview"):
            lines.extend(
                [
                    "",
                    "### Stdout Preview",
                    "",
                    "```text",
                    str(failure_summary["stdout_preview"]),
                    "```",
                ]
            )
    lines.extend(
        [
            "",
        ]
    )
    if reference_preflight:
        lines.extend(
            [
                "## Reference Preflight",
                "",
                f"- Genome: `{reference_preflight['genome_id']}`",
                f"- Prepared before run: `{reference_preflight['prepared_before']}`",
                f"- Prepare attempted: `{reference_preflight['prepare_attempted']}`",
                f"- Prepared after run: `{reference_preflight['prepared_after']}`",
                "",
            ]
        )
        status_before = reference_preflight.get("status_before")
        if isinstance(status_before, dict):
            lines.extend(
                [
                    "### Status Before",
                    "",
                    "```json",
                    json.dumps(status_before, indent=2, ensure_ascii=True),
                    "```",
                    "",
                ]
            )
        status_after = reference_preflight.get("status_after")
        if (
            isinstance(status_after, dict)
            and status_after != status_before
        ):
            lines.extend(
                [
                    "### Status After",
                    "",
                    "```json",
                    json.dumps(status_after, indent=2, ensure_ascii=True),
                    "```",
                    "",
                ]
            )
        lines.extend(["### Preflight Commands", ""])
        for idx, step in enumerate(reference_preflight.get("steps", []), start=1):
            step_command = " ".join(
                shlex.quote(v) for v in step.get("command", [])
            ) or "(none)"
            lines.extend(
                [
                    f"{idx}. `{step.get('status', 'unknown')}`",
                    f"   Command: `{step_command}`",
                    f"   Exit code: `{step.get('exit_code')}`",
                ]
            )
        lines.extend([""])
    lines.extend(
        [
            "## Command",
            "",
            "```bash",
            command_text,
            "```",
            "",
            "## Stdout",
            "",
            "```text",
            stdout.rstrip(),
            "```",
            "",
            "## Stderr",
            "",
            "```text",
            stderr.rstrip(),
            "```",
        ]
    )
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _default_demo_request() -> Request:
    return Request(mode="capabilities", timeout_secs=180)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "GENtle ClawBio skill wrapper. Executes deterministic gentle_cli "
            "commands from structured request JSON, supports sequence-grounded "
            "follow-up to biological observations, and writes reproducibility artifacts."
        )
    )
    parser.add_argument("--input", help="Path to request JSON")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run built-in demo request")
    parser.add_argument(
        "--gentle-cli",
        help="Explicit command used to invoke gentle_cli. Recommended runtime "
        "is Docker/OCI or Apptainer/Singularity via GENTLE_CLI_CMD; examples "
        "include 'gentle_cli', './gentle_apptainer_cli.sh /path/to/gentle.sif', "
        "or 'cargo run --bin gentle_cli --'.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    started = _now_utc_iso()
    request: Request | None = None
    run_result: subprocess.CompletedProcess[str] | None = None
    command: list[str] | None = None
    resolution: CliResolution | None = None
    execution_cwd = Path.cwd()
    status = "failed"
    error_message: str | None = None
    failure_summary: dict[str, Any] | None = None
    collected_artifacts: list[dict[str, str]] = []
    reference_preflight: dict[str, Any] | None = None
    stdout_json: Any | None = None
    chat_summary_lines: list[str] | None = None
    preferred_artifacts: list[dict[str, Any]] | None = None

    try:
        if args.demo:
            request = _default_demo_request()
        else:
            if not args.input:
                raise SkillError("--input is required unless --demo is used")
            payload = _read_json(_resolve_existing_request_file(args.input, Path(__file__)))
            request = _coerce_request(payload)

        try:
            resolution = _resolve_cli(args.gentle_cli, Path(__file__))
        except SkillError as e:
            if args.demo:
                status = "degraded_demo"
                error_message = str(e)
                request = _default_demo_request()
                raise
            raise

        execution_cwd = _resolve_execution_cwd(request, resolution, Path(__file__))
        if request.ensure_reference_prepared is not None:
            reference_preflight = _reference_preflight_record(
                request.ensure_reference_prepared
            )
            _ensure_reference_prepared(
                request.ensure_reference_prepared,
                resolution,
                execution_cwd,
                reference_preflight,
            )

        cli_args = _build_cli_args(request, Path(__file__))
        run_result, main_step = _run_cli_command(
            resolution,
            cli_args,
            execution_cwd,
            request.timeout_secs,
        )
        command = main_step["command"]
        status = "ok" if run_result.returncode == 0 else "command_failed"
        if run_result.returncode != 0:
            failure_summary = _build_failure_summary(
                stage="main_command",
                step=main_step,
                execution_cwd=execution_cwd,
            )
            error_message = _build_failure_message(
                headline="gentle_cli exited with a non-zero status.",
                failure_summary=failure_summary,
            )
        stdout_json = _parse_stdout_json(run_result.stdout)
        chat_summary_lines = _extract_chat_summary_lines(stdout_json)
        preferred_artifacts = _extract_preferred_artifacts(stdout_json)
        collected_artifacts = _copy_collected_artifacts(
            request, output_dir, execution_cwd
        )
        collected_artifacts, preferred_artifacts = _augment_artifacts_with_storyboard(
            request,
            output_dir,
            collected_artifacts,
            preferred_artifacts,
        )
    except subprocess.TimeoutExpired as e:
        request = request if request is not None else _default_demo_request()
        error_message = f"command timed out after {e.timeout} seconds"
        status = "timeout"
    except SkillError as e:
        if (
            failure_summary is None
            and reference_preflight is not None
            and isinstance(reference_preflight.get("last_failure"), dict)
        ):
            failure_summary = reference_preflight["last_failure"]
        error_message = str(e)
        if status != "degraded_demo":
            status = "failed"
    except Exception as e:  # pragma: no cover - defensive boundary
        request = request if request is not None else _default_demo_request()
        error_message = f"unexpected error: {type(e).__name__}: {e}"
        status = "failed"

    ended = _now_utc_iso()

    report_path = output_dir / "report.md"
    result_path = output_dir / "result.json"
    commands_path = repro_dir / "commands.sh"
    env_path = repro_dir / "environment.yml"
    checksums_path = repro_dir / "checksums.sha256"

    _write_report(
        path=report_path,
        request=request,
        resolution=resolution,
        execution_cwd=execution_cwd,
        command=command,
        run_result=run_result,
        stdout_json=stdout_json,
        chat_summary_lines=chat_summary_lines,
        collected_artifacts=collected_artifacts,
        reference_preflight=reference_preflight,
        started_utc=started,
        ended_utc=ended,
        status=status,
        error_message=error_message,
        failure_summary=failure_summary,
        preferred_artifacts=preferred_artifacts,
    )

    command_lines = [
        " ".join(shlex.quote(v) for v in step.get("command", []))
        for step in (reference_preflight or {}).get("steps", [])
        if step.get("command")
    ]
    if command:
        command_lines.append(" ".join(shlex.quote(v) for v in command))
    if not command_lines:
        command_lines = ["# no command executed"]
    commands_text = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"# generated_utc: {ended}",
            f"cd {shlex.quote(str(execution_cwd))}",
            *command_lines,
            "",
        ]
    )
    commands_path.write_text(commands_text, encoding="utf-8")
    os.chmod(commands_path, 0o755)

    _write_repro_environment(env_path)

    result_payload = {
        "schema": RESULT_SCHEMA,
        "status": status,
        "request": (dataclasses.asdict(request) if request is not None else None),
        "started_utc": started,
        "ended_utc": ended,
        "resolver": (dataclasses.asdict(resolution) if resolution else None),
        "command": command,
        "exit_code": (run_result.returncode if run_result else None),
        "stdout": (run_result.stdout if run_result else ""),
        "stdout_json": stdout_json,
        "stderr": (run_result.stderr if run_result else ""),
        "chat_summary_lines": chat_summary_lines,
        "preferred_artifacts": preferred_artifacts,
        "error": error_message,
        "failure_summary": failure_summary,
        "preflight": {
            "reference_preparation": reference_preflight,
        },
        "artifacts": {
            "report_md": str(report_path),
            "result_json": str(result_path),
            "repro_commands": str(commands_path),
            "repro_environment": str(env_path),
            "repro_checksums": str(checksums_path),
            "collected": collected_artifacts,
        },
    }
    result_path.write_text(
        json.dumps(result_payload, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )

    checksums: list[tuple[str, str]] = []
    for artifact in [report_path, result_path, commands_path, env_path]:
        checksums.append((artifact.name, _sha256_file(artifact)))
    checksums_lines = [f"{digest}  {name}" for name, digest in checksums]
    checksums_path.write_text("\n".join(checksums_lines) + "\n", encoding="utf-8")

    print(json.dumps(result_payload, indent=2, ensure_ascii=True))

    return 0 if status in ("ok", "degraded_demo") else 1


if __name__ == "__main__":
    sys.exit(main())
