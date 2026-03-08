"""Deterministic Python adapter over ``gentle_cli`` command routes.

This module intentionally keeps biology logic in GENtle's Rust engine and
provides a process-level adapter for Python scripting/notebook automation.
"""

from __future__ import annotations

import dataclasses
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
from typing import Any, Sequence

_ERROR_CODE_RE = re.compile(r"\b([A-Z][A-Z0-9_]{2,})\b")


@dataclasses.dataclass(frozen=True)
class GentleRunResult:
    """Structured subprocess result for one ``gentle_cli`` invocation."""

    command: tuple[str, ...]
    cwd: str | None
    exit_code: int
    stdout: str
    stderr: str
    json: Any | None = None


class GentleCliError(RuntimeError):
    """Raised when a ``gentle_cli`` command fails or returns invalid JSON."""

    def __init__(
        self,
        message: str,
        *,
        code: str | None = None,
        command: Sequence[str] | None = None,
        exit_code: int | None = None,
        stdout: str = "",
        stderr: str = "",
    ) -> None:
        super().__init__(message)
        self.code = code
        self.command = tuple(command or [])
        self.exit_code = exit_code
        self.stdout = stdout
        self.stderr = stderr


def _find_repo_root(seed: Path) -> Path | None:
    for candidate in [Path.cwd(), *seed.resolve().parents]:
        if (candidate / "Cargo.toml").exists() and (
            candidate / "src" / "bin" / "gentle_cli.rs"
        ).exists():
            return candidate
    return None


def _normalize_cli_cmd(cli_cmd: str | Sequence[str]) -> list[str]:
    if isinstance(cli_cmd, str):
        cmd = shlex.split(cli_cmd)
    else:
        cmd = [str(x) for x in cli_cmd]
    if not cmd:
        raise GentleCliError("Empty CLI command")
    return cmd


def _extract_error_code(stderr: str, stdout: str) -> str | None:
    for text in (stderr, stdout):
        if not text:
            continue
        # Prefer explicit prefixed codes, e.g. AGENT_SCHEMA_VALIDATION: ...
        prefix = text.strip().split(":", 1)[0].strip()
        if prefix and re.fullmatch(r"[A-Z][A-Z0-9_]{2,}", prefix):
            return prefix
        match = _ERROR_CODE_RE.search(text)
        if match:
            return match.group(1)
    return None


class GentleClient:
    """Thin deterministic Python interface to GENtle CLI routes.

    Resolution order for CLI executable:
    1. ``cli_cmd`` constructor argument
    2. ``GENTLE_CLI_CMD`` environment variable
    3. ``gentle_cli`` on ``PATH``
    4. repository fallback: ``cargo run --quiet --bin gentle_cli --``
    """

    def __init__(
        self,
        *,
        state_path: str | None = None,
        cli_cmd: str | Sequence[str] | None = None,
        cwd: str | None = None,
        default_timeout_secs: int = 180,
    ) -> None:
        self.state_path = state_path
        self.default_timeout_secs = int(default_timeout_secs)
        if self.default_timeout_secs <= 0:
            raise GentleCliError("default_timeout_secs must be > 0")

        self._cwd = cwd
        self._cli_prefix, self._resolved_cwd = self._resolve_cli(cli_cmd, cwd)

    @property
    def cli_prefix(self) -> tuple[str, ...]:
        return tuple(self._cli_prefix)

    @property
    def cwd(self) -> str | None:
        return self._resolved_cwd

    def capabilities(self, *, timeout_secs: int | None = None) -> dict[str, Any]:
        result = self.run(["capabilities"], expect_json=True, timeout_secs=timeout_secs)
        return result.json if isinstance(result.json, dict) else {}

    def state_summary(self, *, timeout_secs: int | None = None) -> dict[str, Any]:
        result = self.run(["state-summary"], expect_json=True, timeout_secs=timeout_secs)
        return result.json if isinstance(result.json, dict) else {}

    def op(self, operation: dict[str, Any] | str, *, timeout_secs: int | None = None) -> dict[str, Any]:
        payload = operation if isinstance(operation, str) else json.dumps(operation, separators=(",", ":"), ensure_ascii=True)
        result = self.run(["op", payload], expect_json=True, timeout_secs=timeout_secs)
        return result.json if isinstance(result.json, dict) else {}

    def workflow(
        self,
        workflow: dict[str, Any] | str | None = None,
        *,
        workflow_path: str | None = None,
        timeout_secs: int | None = None,
    ) -> dict[str, Any]:
        if workflow is None and not workflow_path:
            raise GentleCliError("workflow(...) requires 'workflow' or 'workflow_path'")
        if workflow_path:
            arg = f"@{workflow_path}"
        elif isinstance(workflow, str):
            arg = workflow
        else:
            arg = json.dumps(workflow, separators=(",", ":"), ensure_ascii=True)
        result = self.run(["workflow", arg], expect_json=True, timeout_secs=timeout_secs)
        return result.json if isinstance(result.json, dict) else {}

    def shell(
        self,
        line: str,
        *,
        expect_json: bool = False,
        timeout_secs: int | None = None,
    ) -> GentleRunResult:
        if not isinstance(line, str) or not line.strip():
            raise GentleCliError("shell(...) requires non-empty command line text")
        return self.run(["shell", line.strip()], expect_json=expect_json, timeout_secs=timeout_secs)

    def run(
        self,
        args: Sequence[str],
        *,
        include_state: bool = True,
        expect_json: bool = False,
        timeout_secs: int | None = None,
    ) -> GentleRunResult:
        cmd = list(self._cli_prefix)
        if include_state and self.state_path:
            cmd.extend(["--state", self.state_path])
        cmd.extend(str(a) for a in args)
        timeout = self.default_timeout_secs if timeout_secs is None else int(timeout_secs)
        if timeout <= 0:
            raise GentleCliError("timeout_secs must be > 0")

        try:
            proc = subprocess.run(
                cmd,
                cwd=self._resolved_cwd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as e:
            raise GentleCliError(
                f"gentle_cli timed out after {e.timeout} seconds",
                code="GENTLE_TIMEOUT",
                command=cmd,
            ) from e
        except OSError as e:
            raise GentleCliError(
                f"failed to spawn gentle_cli command: {e}",
                code="GENTLE_SPAWN_FAILED",
                command=cmd,
            ) from e

        if proc.returncode != 0:
            code = _extract_error_code(proc.stderr, proc.stdout)
            message = (
                f"gentle_cli command failed with exit code {proc.returncode}"
                + (f" ({code})" if code else "")
            )
            raise GentleCliError(
                message,
                code=code,
                command=cmd,
                exit_code=proc.returncode,
                stdout=proc.stdout,
                stderr=proc.stderr,
            )

        parsed: Any | None = None
        if expect_json:
            raw = proc.stdout.strip()
            if not raw:
                raise GentleCliError(
                    "gentle_cli returned empty stdout when JSON was expected",
                    code="GENTLE_PARSE_ERROR",
                    command=cmd,
                    exit_code=proc.returncode,
                    stdout=proc.stdout,
                    stderr=proc.stderr,
                )
            try:
                parsed = json.loads(raw)
            except json.JSONDecodeError as e:
                raise GentleCliError(
                    f"gentle_cli returned invalid JSON: {e}",
                    code="GENTLE_PARSE_ERROR",
                    command=cmd,
                    exit_code=proc.returncode,
                    stdout=proc.stdout,
                    stderr=proc.stderr,
                ) from e

        return GentleRunResult(
            command=tuple(cmd),
            cwd=self._resolved_cwd,
            exit_code=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            json=parsed,
        )

    def _resolve_cli(
        self,
        cli_cmd: str | Sequence[str] | None,
        cwd: str | None,
    ) -> tuple[list[str], str | None]:
        if cli_cmd is not None:
            return _normalize_cli_cmd(cli_cmd), cwd

        env_cmd = os.environ.get("GENTLE_CLI_CMD", "").strip()
        if env_cmd:
            return _normalize_cli_cmd(env_cmd), cwd

        path_hit = shutil.which("gentle_cli")
        if path_hit:
            return [path_hit], cwd

        repo_root = _find_repo_root(Path(__file__))
        if repo_root is not None and shutil.which("cargo"):
            return ["cargo", "run", "--quiet", "--bin", "gentle_cli", "--"], str(repo_root)

        raise GentleCliError(
            "Could not resolve gentle_cli executable; set cli_cmd or GENTLE_CLI_CMD."
        )
