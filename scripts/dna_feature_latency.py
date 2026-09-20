#!/usr/bin/env python3
"""Build once, then smoke/audit the exact DNA latency binary without Cargo.

Receipts bind source, toolchain/profile, lockfile and executable bytes. These are
CPU observations, never a native-GUI or release performance verdict. Uses only
the Python standard library; the build is explicitly offline.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
PREFIX = "GENTLE_DNA_LATENCY_DIAGNOSTICS="


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path: Path, value: dict) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write("\n")


def command_text(argv: list[str], root: Path) -> str:
    return subprocess.check_output(argv, cwd=root, timeout=30).decode("utf-8").strip()


def source_identity(root: Path) -> dict:
    # Hash dirty input without copying any source/private prose into the receipt.
    delta = subprocess.check_output(["git", "diff", "HEAD", "--binary"], cwd=root, timeout=30)
    untracked = subprocess.check_output(
        ["git", "ls-files", "--others", "--exclude-standard", "-z", "--",
         "src", "crates", "packages", "vendor", "assets", "scripts", "build.rs", ".cargo"],
        cwd=root, timeout=30,
    )
    status = subprocess.check_output(
        ["git", "status", "--porcelain", "-z"], cwd=root, timeout=30,
    )
    extra = {os.fsdecode(name): sha256(root / os.fsdecode(name))
             for name in untracked.split(b"\0") if name}
    return {
        "revision": command_text(["git", "rev-parse", "HEAD"], root),
        "tracked_diff_sha256": hashlib.sha256(delta).hexdigest(),
        "untracked_build_inputs": extra,
        "dirty": bool(status),
        "cargo_lock_sha256": sha256(root / "Cargo.lock"),
    }


def binary_identity(binary: Path) -> dict:
    result = subprocess.run([str(binary), "--gentle-identity"], capture_output=True, check=True, timeout=30)
    identity = json.loads(result.stdout)
    if identity.get("schema") != "gentle.dna_feature_latency_binary.v1":
        raise ValueError("Not a DNA feature latency benchmark binary")
    return identity


def prepare(root: Path, output: Path, profile: str, timeout: int) -> dict:
    output.mkdir(parents=True, exist_ok=False)
    source = source_identity(root)
    argv = ["cargo", "test", "--locked", "--offline", "--profile", profile,
            "-p", "gentle-benchmarks", "--bench", "dna_feature_latency",
            "--no-run", "--message-format=json"]
    started = time.monotonic()
    with (output / "build.stdout.jsonl").open("wb") as stdout, (output / "build.stderr.log").open("wb") as stderr:
        result = subprocess.run(argv, cwd=root, stdout=stdout, stderr=stderr, timeout=timeout)
    if result.returncode:
        raise RuntimeError(f"Build failed ({result.returncode}); retained logs in {output}")
    binaries = []
    for line in (output / "build.stdout.jsonl").read_text(encoding="utf-8").splitlines():
        event = json.loads(line)
        if (event.get("reason") == "compiler-artifact"
                and event.get("target", {}).get("name") == "dna_feature_latency"
                and event.get("executable")):
            binaries.append(Path(event["executable"]).resolve())
    if len(binaries) != 1 or source != source_identity(root):
        raise ValueError("Ambiguous executable or source changed during build; no receipt issued")
    binary = binaries[0]
    identity = binary_identity(binary)
    if identity.get("revision") != source["revision"]:
        raise ValueError("Compiled revision differs from source revision")
    record = {
        "schema": "gentle.dna_feature_latency_build.v1", "source": source,
        "binary": str(binary), "binary_sha256": sha256(binary),
        "profile": profile, "rustc": command_text(["rustc", "-Vv"], root),
        "cargo": command_text(["cargo", "-V"], root), "build_command": argv,
        "build_wall_seconds": time.monotonic() - started,
        "build_environment": {key: os.environ.get(key) for key in
                              ("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "CARGO_BUILD_TARGET", "RUST_MIN_STACK")},
    }
    write_json(output / "build.json", record)
    return record


def load_build(receipt: Path, mode: str) -> tuple[dict, Path]:
    build = json.loads(receipt.read_bytes())
    if build.get("schema") != "gentle.dna_feature_latency_build.v1":
        raise ValueError("Unsupported build receipt")
    if mode == "audit" and build.get("profile") != "bench-audit":
        raise ValueError("Timed audit requires the recorded bench-audit profile")
    if mode == "audit" and build.get("source", {}).get("dirty") is not False:
        raise ValueError("Timed audit requires a clean source receipt")
    binary = Path(build["binary"])
    if sha256(binary) != build["binary_sha256"]:
        raise ValueError("Benchmark binary changed since preparation")
    if binary_identity(binary).get("revision") != build["source"]["revision"]:
        raise ValueError("Binary identity does not match the build receipt")
    return build, binary


def collect_observations(log: Path, revision: str) -> list[dict]:
    rows = []
    fixtures = {}
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(PREFIX):
            row = json.loads(line[len(PREFIX):])
            if row.get("schema") != "gentle.dna_feature_latency_observation.v1" or row.get("revision") != revision:
                raise ValueError("Unbound benchmark observation")
            digest = row.get("fixture_sha256", "")
            if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                raise ValueError("Missing fixture content digest")
            fixture = (row["length_bp"], row["feature_count"])
            if fixtures.setdefault(fixture, digest) != digest:
                raise ValueError("Fixture content changed between interactions")
            rows.append(row)
    keys = {(row["length_bp"], row["feature_count"], row["interaction"]) for row in rows}
    expected = {(length, count, interaction) for length in (20_000, 250_000, 2_000_000)
                for count in (100, 1_000, 10_000)
                for interaction in ("steady", "pan_1bp", "zoom", "toggle_mrna", "select", "hover",
                                    "resize_compact", "resize_desktop", "resize_fullhd")}
    if len(keys) != len(rows) or keys != expected:
        raise ValueError("Missing or duplicated density/interaction observations")
    return rows


def run(receipt: Path, output: Path, mode: str, timeout: int) -> dict:
    build, binary = load_build(receipt, mode)
    output.mkdir(parents=True, exist_ok=False)
    output = output.resolve()
    env = os.environ.copy()
    for key in list(env):
        if key.startswith("GENTLE_") or key.startswith("XDG_"):
            del env[key]
    isolated = {}
    for key, directory in (("HOME", "home"), ("USERPROFILE", "home"),
                           ("XDG_CONFIG_HOME", "config"), ("XDG_CACHE_HOME", "cache"),
                           ("XDG_DATA_HOME", "data"), ("TMPDIR", "tmp"),
                           ("TMP", "tmp"), ("TEMP", "tmp")):
        path = output / directory
        path.mkdir(exist_ok=True)
        env[key] = str(path)
        isolated[key] = str(path)
    # Always isolate Criterion output; do not inherit unrelated audit destinations.
    env["CRITERION_HOME"] = str(output / "criterion")
    env["GENTLE_DNA_CACHE_DIAGNOSTICS"] = "0"
    env["GENTLE_GUI_PROFILE"] = "0"
    env["RUST_MIN_STACK"] = "16777216"
    argv = [str(binary), "--test"] if mode == "smoke" else [str(binary), "--bench", "--quick", "--noplot"]
    log = output / "run.log"
    started = time.monotonic()
    exit_code = None
    failure = None
    try:
        with log.open("wb") as stream:
            exit_code = subprocess.run(argv, cwd=output, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        failure = "timeout"
    except OSError as error:
        failure = f"launch failed: {error}"
    record = {
        "schema": "gentle.dna_feature_latency_run.v1", "build": build,
        "build_receipt_sha256": sha256(receipt), "mode": mode,
        "runner_sha256": sha256(Path(__file__)),
        "host": {"system": platform.system(), "release": platform.release(),
                 "version": platform.version(), "machine": platform.machine()},
        "isolated_environment": isolated,
        "runtime_environment": {key: env.get(key) for key in
                                ("LANG", "LC_ALL", "LC_CTYPE", "RAYON_NUM_THREADS", "RUST_MIN_STACK")},
        "command": argv, "exit_code": exit_code, "failure": failure,
        "process_wall_seconds": time.monotonic() - started,
        "log_sha256": sha256(log), "native_gui_measured": False,
        "observations": [], "statistics": [],
    }
    try:
        if exit_code != 0:
            raise ValueError(f"Benchmark did not complete: {failure or exit_code}")
        record["observations"] = collect_observations(log, build["source"]["revision"])
        for path in sorted((output / "criterion").rglob("new/estimates.json")):
            record["statistics"].append({
                "case": str(path.parent.parent.relative_to(output / "criterion")),
                "estimates": json.loads(path.read_bytes()), "sha256": sha256(path),
            })
        if mode == "audit" and len(record["statistics"]) != 117:
            raise ValueError("Timed audit did not produce all 117 Criterion estimates")
    except (ValueError, KeyError) as error:
        record["failure"] = str(error)
    write_json(output / "run.json", record)
    columns = ["length_bp", "feature_count", "interaction", "tree_builds", "layer_builds", "layer_gc_bases", "layouts"]
    with (output / "work.tsv").open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("\t".join(columns) + "\n")
        for row in record["observations"]:
            before, after = row["before"], row["after"]
            values = [row[key] for key in columns[:3]]
            values.extend(after[key] - before[key] for key in columns[3:6])
            values.append(after["linear"]["layouts"] - before["linear"]["layouts"])
            stream.write("\t".join(map(str, values)) + "\n")
    if record["failure"]:
        raise ValueError(f"{record['failure']}; receipt retained in {output}")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("prepare")
    build.add_argument("--profile", choices=("dev", "bench-audit"), default="bench-audit")
    build.add_argument("--output", type=Path, required=True)
    build.add_argument("--timeout", type=int, default=7200)
    execute = commands.add_parser("run")
    execute.add_argument("--receipt", type=Path, required=True)
    execute.add_argument("--output", type=Path, required=True)
    execute.add_argument("--mode", choices=("smoke", "audit"), default="smoke")
    execute.add_argument("--repeats", type=int, choices=(1, 2), default=2)
    execute.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    if args.command == "prepare":
        prepare(ROOT, args.output, args.profile, args.timeout)
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        for index in range(args.repeats):
            run(args.receipt.resolve(), args.output / f"repeat-{index + 1}", args.mode, args.timeout)


if __name__ == "__main__":
    main()
