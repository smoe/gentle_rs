#!/usr/bin/env python3
"""Cargo test runner that removes the stack override at the final process boundary.

Cargo can restore RUST_MIN_STACK from .cargo/config.toml after a shell unset.
This wrapper preserves Cargo's runtime environment except that one variable.
Use only for the focused default-stack regressions, not the full suite.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


def cargo_runner_config() -> str:
    command = [sys.executable, str(Path(__file__).resolve())]
    return 'target."cfg(all())".runner = ' + json.dumps(command)


def run_test(command: list[str]) -> int:
    if not command:
        raise ValueError("Cargo must supply a test executable")
    environment = os.environ.copy()
    inherited_stack = environment.pop("RUST_MIN_STACK", None)
    print(json.dumps({
        "schema": "gentle.default_stack_test_runner.v1",
        "command": command,
        "inherited_rust_min_stack": inherited_stack,
        "child_rust_min_stack": environment.get("RUST_MIN_STACK"),
    }), flush=True)
    result = subprocess.run(command, env=environment, check=False)
    return result.returncode if result.returncode >= 0 else 128 - result.returncode


def main(arguments: list[str]) -> int:
    if arguments == ["--cargo-config"]:
        print(cargo_runner_config())
        return 0
    return run_test(arguments)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
