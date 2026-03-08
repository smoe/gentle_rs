#!/usr/bin/env python3
"""ClawBio/OpenClaw wrapper skill for deterministic GENtle CLI execution."""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
from typing import Any

REQUEST_SCHEMA = "gentle.clawbio_skill_request.v1"
RESULT_SCHEMA = "gentle.clawbio_skill_result.v1"


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


@dataclasses.dataclass
class CliResolution:
    argv_prefix: list[str]
    cwd: str | None
    label: str


def _now_utc_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
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
        "Set --gentle-cli, or GENTLE_CLI_CMD, or install gentle_cli on PATH."
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
    )
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
    return request


def _json_arg(value: Any) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=True, separators=(",", ":"))


def _build_cli_args(request: Request) -> list[str]:
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
            args.extend(["workflow", f"@{request.workflow_path}"])
        else:
            args.extend(["workflow", _json_arg(request.workflow)])
    elif request.mode == "raw":
        args.extend(request.raw_args or [])
    else:
        raise SkillError(f"Unsupported mode '{request.mode}'")
    return args


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


def _write_report(
    path: Path,
    request: Request,
    resolution: CliResolution | None,
    command: list[str] | None,
    run_result: subprocess.CompletedProcess[str] | None,
    started_utc: str,
    ended_utc: str,
    status: str,
    error_message: str | None,
) -> None:
    command_text = " ".join(shlex.quote(v) for v in command) if command else "(none)"
    stdout = run_result.stdout if run_result else ""
    stderr = run_result.stderr if run_result else ""
    exit_code = run_result.returncode if run_result else None
    lines = [
        "# GENtle ClawBio Skill Report",
        "",
        f"- Started (UTC): `{started_utc}`",
        f"- Ended (UTC): `{ended_utc}`",
        f"- Status: `{status}`",
        f"- Mode: `{request.mode}`",
    ]
    if request.state_path:
        lines.append(f"- State path: `{request.state_path}`")
    if resolution is not None:
        lines.append(f"- Resolver: `{resolution.label}`")
        lines.append(f"- Resolver cwd: `{resolution.cwd or os.getcwd()}`")
    if exit_code is not None:
        lines.append(f"- Exit code: `{exit_code}`")
    if error_message:
        lines.append(f"- Error: `{error_message}`")
    lines.extend(
        [
            "",
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
            "commands from structured request JSON and writes reproducibility artifacts."
        )
    )
    parser.add_argument("--input", help="Path to request JSON")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run built-in demo request")
    parser.add_argument(
        "--gentle-cli",
        help="Explicit command used to invoke gentle_cli (e.g. 'gentle_cli' or "
        "'cargo run --bin gentle_cli --')",
    )
    args = parser.parse_args()

    output_dir = Path(args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    started = _now_utc_iso()
    request: Request
    run_result: subprocess.CompletedProcess[str] | None = None
    command: list[str] | None = None
    resolution: CliResolution | None = None
    status = "failed"
    error_message: str | None = None

    try:
        if args.demo:
            request = _default_demo_request()
        else:
            if not args.input:
                raise SkillError("--input is required unless --demo is used")
            payload = _read_json(Path(args.input))
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

        cli_args = _build_cli_args(request)
        command = resolution.argv_prefix + cli_args
        run_result = subprocess.run(
            command,
            cwd=resolution.cwd,
            capture_output=True,
            text=True,
            timeout=request.timeout_secs,
            check=False,
        )
        status = "ok" if run_result.returncode == 0 else "command_failed"
        if run_result.returncode != 0:
            error_message = (
                f"gentle_cli exited with {run_result.returncode}; inspect stderr in report.md"
            )
    except subprocess.TimeoutExpired as e:
        request = request if "request" in locals() else _default_demo_request()
        error_message = f"command timed out after {e.timeout} seconds"
        status = "timeout"
    except SkillError as e:
        request = request if "request" in locals() else _default_demo_request()
        error_message = str(e)
        if status != "degraded_demo":
            status = "failed"
    except Exception as e:  # pragma: no cover - defensive boundary
        request = request if "request" in locals() else _default_demo_request()
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
        command=command,
        run_result=run_result,
        started_utc=started,
        ended_utc=ended,
        status=status,
        error_message=error_message,
    )

    command_line = " ".join(shlex.quote(v) for v in command) if command else "# no command executed"
    commands_text = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"# generated_utc: {ended}",
            (f"cd {shlex.quote(resolution.cwd)}" if resolution and resolution.cwd else ""),
            command_line,
            "",
        ]
    )
    commands_path.write_text(commands_text, encoding="utf-8")
    os.chmod(commands_path, 0o755)

    _write_repro_environment(env_path)

    result_payload = {
        "schema": RESULT_SCHEMA,
        "status": status,
        "request": dataclasses.asdict(request),
        "started_utc": started,
        "ended_utc": ended,
        "resolver": (dataclasses.asdict(resolution) if resolution else None),
        "command": command,
        "exit_code": (run_result.returncode if run_result else None),
        "stdout": (run_result.stdout if run_result else ""),
        "stderr": (run_result.stderr if run_result else ""),
        "error": error_message,
        "artifacts": {
            "report_md": str(report_path),
            "result_json": str(result_path),
            "repro_commands": str(commands_path),
            "repro_environment": str(env_path),
            "repro_checksums": str(checksums_path),
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
