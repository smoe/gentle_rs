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
LEGACY_WORKLOAD = "density_ladder_v1"
BOUNDARY_WORKLOAD = "density_boundary_v1"
INTERACTIVE_BOUNDARY = (250_000, 5_000)
INTERACTIONS = ("steady", "pan_1bp", "zoom", "toggle_mrna", "select", "hover",
                "resize_compact", "resize_desktop", "resize_fullhd")
CASES = ("constructor", "hydration", "first_frame_tree_deferred", "first_frame_tree_loaded", *INTERACTIONS)


def workload_fixtures(workload: str) -> list[tuple[int, int]]:
    if workload not in (LEGACY_WORKLOAD, BOUNDARY_WORKLOAD):
        raise ValueError(f"Unsupported DNA latency workload: {workload}")
    fixtures = [(length, count) for length in (20_000, 250_000, 2_000_000)
                for count in (100, 1_000, 10_000)]
    return fixtures + ([INTERACTIVE_BOUNDARY] if workload == BOUNDARY_WORKLOAD else [])


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
    workload_fixtures(identity.get("workload", LEGACY_WORKLOAD))
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
        "workload": identity.get("workload", LEGACY_WORKLOAD),
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
    identity = binary_identity(binary)
    if identity.get("revision") != build["source"]["revision"]:
        raise ValueError("Binary identity does not match the build receipt")
    workload = build.get("workload", LEGACY_WORKLOAD)
    workload_fixtures(workload)
    if identity.get("workload", LEGACY_WORKLOAD) != workload:
        raise ValueError("Binary workload does not match the build receipt")
    return build, binary


def validate_layers(row: dict) -> None:
    boundary = (row["length_bp"], row["feature_count"]) == INTERACTIVE_BOUNDARY
    if row.get("derived_layers") is not boundary:
        raise ValueError("Wrong derived-layer configuration")
    for side in ("layers_before", "layers_after"):
        layers = row.get(side)
        if not isinstance(layers, dict) or layers.get("annotation_features") != row["feature_count"]:
            raise ValueError("Missing or inconsistent layer inventory")
        viewport = layers.get("viewport")
        if (not isinstance(viewport, list) or len(viewport) != 2
                or any(type(value) is not int for value in viewport)
                or not 0 <= viewport[0] < viewport[1] <= row["length_bp"]):
            raise ValueError("Invalid layer inventory viewport")
        if type(layers.get("gc_bin_size_bp")) is not int or layers["gc_bin_size_bp"] <= 0:
            raise ValueError("Invalid GC-bin configuration")
        for name in ("restriction_groups", "gc_bins", "orfs", "methylation_sites"):
            layer = layers.get(name)
            if (not isinstance(layer, dict) or type(layer.get("enabled")) is not bool
                    or any(type(layer.get(key)) is not int for key in ("total", "viewport_eligible"))
                    or not 0 <= layer["viewport_eligible"] <= layer["total"]):
                raise ValueError(f"Invalid derived-layer inventory: {name}")
            if boundary and (not layer["enabled"] or layer["viewport_eligible"] == 0):
                raise ValueError(f"Boundary derived layer unavailable: {name}")


def collect_observations(log: Path, revision: str, workload: str = LEGACY_WORKLOAD) -> list[dict]:
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
            if workload == BOUNDARY_WORKLOAD:
                if row.get("workload") != workload:
                    raise ValueError("Observation workload does not match the build receipt")
                validate_layers(row)
            rows.append(row)
    keys = {(row["length_bp"], row["feature_count"], row["interaction"]) for row in rows}
    expected = {(length, count, interaction) for length, count in workload_fixtures(workload)
                for interaction in INTERACTIONS}
    if len(keys) != len(rows) or keys != expected:
        raise ValueError("Missing or duplicated density/interaction observations")
    return rows


def collect_statistics(directory: Path, rows: list[dict]) -> list[dict]:
    expected = {
        f"dna_feature_latency/{row['length_bp']}bp_{row['feature_count']}features_{row['fixture_sha256'][:12]}/{case}"
        for row in rows for case in CASES
    }
    statistics = []
    for path in sorted(directory.rglob("new/estimates.json")):
        identity_path = path.with_name("benchmark.json")
        if not identity_path.is_file():
            raise ValueError("Missing Criterion case identity")
        identity = json.loads(identity_path.read_bytes())
        if not isinstance(identity.get("full_id"), str):
            raise ValueError("Invalid Criterion case identity")
        statistics.append({
            "case": identity.get("full_id"), "benchmark_sha256": sha256(identity_path),
            "estimates": json.loads(path.read_bytes()), "sha256": sha256(path),
        })
    actual = {row["case"] for row in statistics}
    if len(statistics) != len(actual) or actual != expected:
        raise ValueError(f"Timed audit did not produce all {len(expected)} Criterion estimates with matching case identities")
    return statistics


def run(receipt: Path, output: Path, mode: str, timeout: int) -> dict:
    build, binary = load_build(receipt, mode)
    workload = build.get("workload", LEGACY_WORKLOAD)
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
        "workload": workload, "interactive_boundary_exercised": False,
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
        record["observations"] = collect_observations(log, build["source"]["revision"], workload)
        if mode == "audit":
            record["statistics"] = collect_statistics(output / "criterion", record["observations"])
        record["interactive_boundary_exercised"] = workload == BOUNDARY_WORKLOAD
    except (ValueError, KeyError) as error:
        record["failure"] = str(error)
    write_json(output / "run.json", record)
    columns = ["length_bp", "feature_count", "interaction", "tree_builds", "layer_builds", "layer_gc_bases", "layouts"]
    layer_names = ("restriction_groups", "gc_bins", "orfs", "methylation_sites")
    with (output / "work.tsv").open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("\t".join([*columns, *(f"{name}_viewport_eligible" for name in layer_names)]) + "\n")
        for row in record["observations"]:
            before, after = row["before"], row["after"]
            values = [row[key] for key in columns[:3]]
            values.extend(after[key] - before[key] for key in columns[3:6])
            values.append(after["linear"]["layouts"] - before["linear"]["layouts"])
            layers = row.get("layers_after", {})
            values.extend(layers.get(name, {}).get("viewport_eligible", "") for name in layer_names)
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
