#!/usr/bin/env python3
"""Run a release build while streaming sampled child-process resource data.

The monitor is intentionally dependency-free.  Its samples go both to the
Actions stream and to the retained build log, so an abrupt runner shutdown can
still leave useful evidence even when the artifact-upload step never starts.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import threading
import time
from typing import Iterable


def parse_process_table(text: str, root_pid: int) -> dict[str, int]:
    """Return sampled RSS totals for root_pid's process tree."""
    rows: list[tuple[int, int, int, str]] = []
    for line in text.splitlines():
        fields = line.strip().split(maxsplit=3)
        if len(fields) != 4:
            continue
        try:
            pid, parent, rss = map(int, fields[:3])
        except ValueError:
            continue
        rows.append((pid, parent, rss, fields[3]))

    descendants = {root_pid}
    for _ in range(len(rows) + 1):
        added = {pid for pid, parent, _rss, _name in rows if parent in descendants}
        if added <= descendants:
            break
        descendants.update(added)

    tree_rss = 0
    rustc_rss = 0
    rustc_processes = 0
    largest_rustc = 0
    for pid, _parent, rss, name in rows:
        if pid not in descendants:
            continue
        tree_rss += rss
        if Path(name).name == "rustc":
            rustc_rss += rss
            rustc_processes += 1
            largest_rustc = max(largest_rustc, rss)
    return {
        "tree_rss_kib": tree_rss,
        "rustc_rss_kib": rustc_rss,
        "rustc_processes": rustc_processes,
        "largest_rustc_rss_kib": largest_rustc,
    }


def process_sample(root_pid: int) -> dict[str, int]:
    completed = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss=,comm="],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return {}
    return parse_process_table(completed.stdout, root_pid)


def _read_linux_value(path: Path) -> str | None:
    try:
        return (
            path.read_text(encoding="utf-8")
            .strip()
            .replace(" ", ":")
            .replace("\n", ",")
        )
    except OSError:
        return None


def linux_memory_sample(root: Path = Path("/")) -> dict[str, str | int]:
    values: dict[str, str | int] = {}
    meminfo = root / "proc/meminfo"
    try:
        for line in meminfo.read_text(encoding="utf-8").splitlines():
            key, value = line.split(":", 1)
            if key in {"MemAvailable", "SwapFree", "SwapTotal"}:
                values[f"linux_{key.lower()}_kib"] = int(value.split()[0])
    except (OSError, ValueError, IndexError):
        pass
    try:
        pressure = (root / "proc/pressure/memory").read_text(encoding="utf-8")
    except OSError:
        pressure = None
    if pressure is not None:
        match = re.search(r"some avg10=([0-9.]+)", pressure)
        if match:
            values["linux_pressure_some_avg10"] = match.group(1)
        match = re.search(r"full avg10=([0-9.]+)", pressure)
        if match:
            values["linux_pressure_full_avg10"] = match.group(1)
    cgroup = root / "sys/fs/cgroup"
    for filename in ("memory.current", "memory.max", "memory.events"):
        value = _read_linux_value(cgroup / filename)
        if value is not None:
            values[f"cgroup_{filename.replace('.', '_')}"] = value
    return values


def macos_memory_sample() -> dict[str, str | int]:
    values: dict[str, str | int] = {}
    vm = subprocess.run(["vm_stat"], check=False, capture_output=True, text=True)
    if vm.returncode == 0:
        page_match = re.search(r"page size of (\d+) bytes", vm.stdout)
        page_size = int(page_match.group(1)) if page_match else 4096
        pages = 0
        for label in ("Pages free", "Pages inactive", "Pages speculative"):
            match = re.search(rf"^{label}:\s+(\d+)\.", vm.stdout, re.MULTILINE)
            if match:
                pages += int(match.group(1))
        values["macos_available_estimate_kib"] = pages * page_size // 1024
    swap = subprocess.run(
        ["sysctl", "-n", "vm.swapusage"], check=False, capture_output=True, text=True
    )
    if swap.returncode == 0:
        free = re.search(r"free = ([0-9.]+)([MG])", swap.stdout)
        if free:
            amount = float(free.group(1)) * (1024 if free.group(2) == "G" else 1)
            values["macos_swap_free_mib"] = int(amount)
    return values


def system_memory_sample() -> dict[str, str | int]:
    if sys.platform == "linux":
        return linux_memory_sample()
    if sys.platform == "darwin":
        return macos_memory_sample()
    return {}


def format_sample(kind: str, values: dict[str, object]) -> str:
    timestamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    fields = " ".join(f"{key}={value}" for key, value in sorted(values.items()))
    return f"release_resource_{kind} timestamp={timestamp} {fields}\n"


def run(command: list[str], log_path: Path, sample_seconds: float) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("ab", buffering=0) as log:
        lock = threading.RLock()

        def emit(data: bytes) -> None:
            with lock:
                sys.stdout.buffer.write(data)
                sys.stdout.buffer.flush()
                log.write(data)

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )

        received_signal: list[int] = []

        def forward(signum: int, _frame: object) -> None:
            received_signal.append(signum)
            emit(format_sample("signal", {"signal": signum}).encode())
            try:
                os.killpg(process.pid, signum)
            except ProcessLookupError:
                pass

        previous = {
            signum: signal.signal(signum, forward)
            for signum in (signal.SIGINT, signal.SIGTERM)
        }

        def stream_output() -> None:
            assert process.stdout is not None
            for chunk in iter(process.stdout.readline, b""):
                emit(chunk)

        reader = threading.Thread(target=stream_output, name="release-build-output")
        reader.start()
        peak = {
            "peak_tree_rss_kib": 0,
            "peak_rustc_rss_kib": 0,
            "peak_largest_rustc_rss_kib": 0,
        }
        try:
            while process.poll() is None:
                sample: dict[str, object] = {"root_pid": process.pid}
                sample.update(process_sample(process.pid))
                sample.update(system_memory_sample())
                peak["peak_tree_rss_kib"] = max(
                    peak["peak_tree_rss_kib"], int(sample.get("tree_rss_kib", 0))
                )
                peak["peak_rustc_rss_kib"] = max(
                    peak["peak_rustc_rss_kib"], int(sample.get("rustc_rss_kib", 0))
                )
                peak["peak_largest_rustc_rss_kib"] = max(
                    peak["peak_largest_rustc_rss_kib"],
                    int(sample.get("largest_rustc_rss_kib", 0)),
                )
                emit(format_sample("sample", sample).encode())
                time.sleep(sample_seconds)
            returncode = process.wait()
        finally:
            reader.join(timeout=30)
            for signum, handler in previous.items():
                signal.signal(signum, handler)

        peak["command_returncode"] = returncode
        peak["received_signal"] = received_signal[-1] if received_signal else "none"
        emit(format_sample("peak", peak).encode())
        if received_signal:
            return 128 + received_signal[-1]
        return 128 - returncode if returncode < 0 else returncode


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--sample-seconds", type=float, default=2.0)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    if args.sample_seconds <= 0:
        parser.error("--sample-seconds must be positive")
    return args


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    return run(args.command, args.log, args.sample_seconds)


if __name__ == "__main__":
    raise SystemExit(main())
