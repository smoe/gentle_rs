#!/usr/bin/env python3
"""Summarize one saved GUI startup trace offline, without running or timing GENtle.

Preserves the raw input and reports inclusive span durations and same-subject
boundaries, never a sum of nested/parallel work or a native performance verdict.
Uses only the Python standard library. Output must be a new directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

MAX_BYTES = 1024 * 1024
MAX_EVENTS = 512
U64_MAX = 2**64 - 1
# Keep this vocabulary in step with src/gui_profiler/startup_trace.rs (tested).
CHECKPOINTS = frozenset((
    "main_entered", "root_first_frame", "root_workspace_frame",
    "dna_worker_scheduled", "dna_worker_result", "dna_load_failed",
    "dna_native_content_frame", "dna_embedded_content_frame",
))
SPANS = frozenset((
    "native_run", "app_initialize", "app_defaults", "help_preparation",
    "configuration_load", "credential_refresh", "initial_state_preparation",
    "project_load", "project_read_decode", "project_install", "dna_open_dispatch",
    "dna_construct", "dna_placeholder_construct", "dna_engine_read_lock",
    "dna_sequence_clone", "dna_hydrate",
    "help_manuals", "help_shell_reference", "help_tutorial_discovery",
    "help_tutorial_selected_load", "help_open", "help_tutorial_open",
    "help_tutorial_menu_discovery", "help_tutorial_switch",
))
HELP_IMAGE_COUNTERS = (
    "svg_references", "cache_hits", "preparation_failures", "rasterization_attempts",
    "rasterization_completed", "rasterization_failures", "rasterization_us",
)
METADATA_STRINGS = ("source_revision", "git_commit", "os", "architecture", "clock", "scope")
METADATA_FLAGS = ("debug_assertions", "gui_test_support", "gui_profiler")
LIMITATIONS = [
    "CPU-side monotonic wall time since Rust main; excludes OS process loading.",
    "Inclusive spans may overlap or nest; do not sum them or treat them as CPU usage.",
    "native_run is the application loop's lifetime, not a startup duration.",
    "Content markers are not fully ready, subject-verified or compositor-presented content.",
    "No input responsiveness, native GUI acceptance or performance verdict is established.",
    "Revision and build flags are producer-reported, not independently verified binary identity.",
    "Retain binary/input hashes, effective profile, host/toolchain and cold/warm conditions separately.",
    "Process-local subjects cannot identify a biological sequence or link separate runs.",
    "Help image counters cover the whole session, including reloads; not unique files or startup-only work.",
    "Rasterization time includes file I/O and font loading; it overlaps help spans, not an additional cost.",
    "First help/tutorial spans measure handler work; menu discovery excludes other menu work and presentation.",
]


def unsigned(value: object, label: str) -> int:
    if type(value) is not int or not 0 <= value <= U64_MAX:
        raise ValueError(f"{label} must be an unsigned 64-bit integer")
    return value


def unique_object(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate JSON field")
        result[key] = value
    return result


def help_image_work(trace: dict) -> dict:
    """Optional v1 extension: absence is unavailable, loss is not measured zero."""
    if "help_image_work" not in trace and "dropped_help_image_observations" not in trace:
        return {"status": "unavailable", "reason": "producer_did_not_record", "recorded": None}
    work = trace.get("help_image_work")
    if (not isinstance(work, dict) or set(work) != {*HELP_IMAGE_COUNTERS, "saturated"}
            or type(work["saturated"]) is not bool):
        raise ValueError("Invalid help image work fields")
    for key in HELP_IMAGE_COUNTERS:
        unsigned(work[key], key)
    dropped = unsigned(trace.get("dropped_help_image_observations"), "dropped_help_image_observations")
    reason = "dropped_observations" if dropped else "saturated" if work["saturated"] else None
    if reason is None:
        if (work["rasterization_failures"] > work["rasterization_completed"]
                or work["rasterization_completed"] > work["rasterization_attempts"]
                or work["cache_hits"] + work["preparation_failures"] + work["rasterization_attempts"] > work["svg_references"]
                or (work["rasterization_completed"] == 0 and work["rasterization_us"] != 0)):
            raise ValueError("Inconsistent help image work counts")
        if (work["rasterization_completed"] != work["rasterization_attempts"]
                or work["cache_hits"] + work["preparation_failures"] + work["rasterization_attempts"] != work["svg_references"]):
            reason = "unfinished_work"
    return {
        "status": "incomplete" if reason else "no_recorded_loss", "reason": reason,
        "recorded": work, "dropped_observations": dropped,
        "rasterization_total_us": work["rasterization_us"] if reason is None else None,
    }


def summarize(raw: bytes) -> dict:
    if len(raw) > MAX_BYTES:
        raise ValueError("Startup trace exceeds the 1 MiB input limit")
    try:
        trace = json.loads(raw, object_pairs_hook=unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError) as error:
        raise ValueError("Invalid startup trace JSON") from error
    if not isinstance(trace, dict) or trace.get("schema") != "gentle.gui_startup_trace.v1":
        raise ValueError("Unsupported startup trace schema")
    for key in METADATA_STRINGS:
        value = trace.get(key)
        if (not isinstance(value, str) or not value or len(value) > 512
                or any(ord(char) < 32 for char in value)):
            raise ValueError(f"Invalid trace metadata: {key}")
    for key in METADATA_FLAGS:
        if type(trace.get(key)) is not bool:
            raise ValueError(f"Invalid trace metadata: {key}")
    limit = unsigned(trace.get("event_limit"), "event_limit")
    dropped = unsigned(trace.get("dropped_events"), "dropped_events")
    events = trace.get("events")
    if not 0 < limit <= MAX_EVENTS or not isinstance(events, list) or len(events) > limit:
        raise ValueError("Invalid event limit or event count")
    subjects = {}
    spans = {}
    previous_time = 0
    for event in events:
        if not isinstance(event, dict):
            raise ValueError("Invalid event object")
        elapsed = unsigned(event.get("elapsed_us"), "elapsed_us")
        subject = unsigned(event.get("subject"), "subject")
        phase, kind = event.get("phase"), event.get("kind")
        if not isinstance(phase, str) or not isinstance(kind, str):
            raise ValueError("Invalid event phase or kind")
        if elapsed < previous_time:
            raise ValueError("Event timestamps are not monotonic")
        previous_time = elapsed
        row = subjects.setdefault(subject, {"subject": subject, "checkpoints": {}, "spans": []})
        if kind == "checkpoint":
            if phase not in CHECKPOINTS or "span" in event:
                raise ValueError("Invalid checkpoint phase or span")
            if phase in row["checkpoints"]:
                raise ValueError("Duplicate checkpoint for one subject")
            row["checkpoints"][phase] = elapsed
        elif kind in ("begin", "completed", "failed", "interrupted"):
            if phase not in SPANS:
                raise ValueError("Invalid span phase")
            span_id = unsigned(event.get("span"), "span")
            if span_id == 0:
                raise ValueError("Span identity must be positive")
            entry = spans.setdefault(span_id, {"span": span_id, "subject": subject, "phase": phase})
            if (entry["subject"], entry["phase"]) != (subject, phase):
                raise ValueError("Span identity has conflicting subjects or phases")
            side = "begin" if kind == "begin" else "end"
            if side in entry:
                raise ValueError("Duplicate span endpoint")
            if side == "begin" and "end" in entry:
                raise ValueError("Span end recorded before its begin")
            entry[side] = elapsed
            if side == "end":
                entry["outcome"] = kind
        else:
            raise ValueError("Unknown event kind")

    warnings = []
    images = help_image_work(trace)
    if images["status"] != "no_recorded_loss":
        warnings.append(f"Help image work {images['status']}: {images['reason']}; not a measured complete total.")
    if dropped:
        warnings.append(f"{dropped} events dropped; all derived durations are unavailable (loss cannot be localized).")
    incomplete = bool(dropped) or images["status"] == "incomplete"
    for entry in spans.values():
        begin, end = entry.get("begin"), entry.get("end")
        if begin is not None and end is not None and end < begin:
            raise ValueError("Span ends before it begins")
        reason = None
        if begin is None or end is None:
            reason = "missing_begin" if begin is None else "missing_end"
            warnings.append(f"Span {entry['span']} has {reason}; no duration inferred.")
            incomplete = True
        elif dropped:
            reason = "dropped_events"
        outcome = entry.get("outcome", "unclosed")
        if outcome == "interrupted":
            incomplete = True
        subjects[entry["subject"]]["spans"].append({
            "span": entry["span"], "phase": entry["phase"], "outcome": outcome,
            "begin_us": begin, "end_us": end,
            "duration_us": end - begin if reason is None else None,
            "unavailable_reason": reason,
        })

    main = subjects.get(0, {}).get("checkpoints", {}).get("main_entered")
    if main is None:
        warnings.append("main_entered checkpoint absent; trace completeness is not established.")
        incomplete = True
    if subjects.get(0, {}).get("checkpoints", {}).get("root_workspace_frame") is None:
        warnings.append("No root workspace CPU marker recorded; a splash or dispatch is not usable content.")
    for row in subjects.values():
        row["spans"].sort(key=lambda span: span["span"])
        row["boundaries"] = boundaries(row, dropped)
        if "dna_load_failed" in row["checkpoints"]:
            warnings.append(f"Subject {row['subject']} recorded dna_load_failed; do not infer successful loading.")
        if any(span["phase"] == "dna_construct" for span in row["spans"]) and not any(
                phase in row["checkpoints"] for phase in ("dna_native_content_frame", "dna_embedded_content_frame")):
            warnings.append(f"Subject {row['subject']} has no DNA content CPU marker; absence is not zero latency.")
    failures = sum(entry.get("outcome") == "failed" for entry in spans.values())
    if failures:
        warnings.append(f"{failures} failed spans retained; failed work is not successful startup.")
    return {
        "schema": "gentle.gui_startup_summary.v1",
        "input": {"schema": trace["schema"], "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)},
        "producer": {key: trace[key] for key in (*METADATA_STRINGS, *METADATA_FLAGS)},
        "event_count": len(events), "dropped_events": dropped,
        "trace_completeness": "incomplete" if incomplete else "no_recorded_loss",
        "failed_spans": failures, "native_gui_measured": False,
        "help_image_work": images,
        "performance_verdict": "not_assessed", "limitations": LIMITATIONS,
        "warnings": warnings, "subjects": [subjects[key] for key in sorted(subjects)],
    }


def boundaries(row: dict, dropped: int) -> list[dict]:
    """Only pair explicit unique endpoints within one subject; never cross windows."""
    def endpoint(phase: str, side: str) -> tuple[int | None, str | None]:
        if side == "checkpoint":
            value = row["checkpoints"].get(phase)
            return value, None if value is not None else "missing_endpoint"
        candidates = [span for span in row["spans"] if span["phase"] == phase]
        if len(candidates) != 1:
            return None, "missing_endpoint" if not candidates else "ambiguous_endpoint"
        span = candidates[0]
        if span["outcome"] != "completed" or span["begin_us"] is None or span["end_us"] is None:
            return None, "unsuccessful_or_incomplete_span"
        return span[f"{side}_us"], None

    specs = [
        ("native_setup_before_app_constructor", "native_run", "begin", "app_initialize", "begin"),
        ("worker_schedule_to_lock_attempt", "dna_worker_scheduled", "checkpoint", "dna_engine_read_lock", "begin"),
        ("worker_result_to_hydration", "dna_worker_result", "checkpoint", "dna_hydrate", "begin"),
        ("hydration_to_native_content_cpu", "dna_hydrate", "end", "dna_native_content_frame", "checkpoint"),
        ("hydration_to_embedded_content_cpu", "dna_hydrate", "end", "dna_embedded_content_frame", "checkpoint"),
        ("construct_to_native_content_cpu", "dna_construct", "begin", "dna_native_content_frame", "checkpoint"),
        ("construct_to_embedded_content_cpu", "dna_construct", "begin", "dna_embedded_content_frame", "checkpoint"),
    ]
    result = []
    phases = set(row["checkpoints"]) | {span["phase"] for span in row["spans"]}
    for name, start_phase, start_side, end_phase, end_side in specs:
        if start_phase not in phases and end_phase not in phases:
            continue
        start, start_reason = endpoint(start_phase, start_side)
        end, end_reason = endpoint(end_phase, end_side)
        reason = "dropped_events" if dropped else start_reason or end_reason
        if reason is None and end < start:
            reason = "endpoint_order"
        result.append({"name": name, "start_us": start, "end_us": end,
                       "duration_us": end - start if reason is None else None,
                       "unavailable_reason": reason})
    return result


def markdown(report: dict) -> str:
    def safe(value: object) -> str:
        return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("|", "&#124;").replace("`", "&#96;")

    def ms(value: int | None) -> str:
        return "unavailable" if value is None else f"{value / 1000:.3f}"

    lines = ["# GUI Startup CPU Attribution", "", "Performance verdict: **not assessed**.", "",
             f"Input SHA-256: `{report['input']['sha256']}`",
             f"Producer revision: {safe(report['producer']['source_revision'])}",
             f"Recorded events: {report['event_count']}; dropped: {report['dropped_events']}; "
             f"completeness: {report['trace_completeness']}.", ""]
    lines.extend(f"- {item}" for item in report["limitations"])
    lines.extend(["", "## Warnings", ""])
    lines.extend(f"- {item}" for item in report["warnings"] or ["No recorded loss; this is not a successful-startup or responsiveness verdict."])
    images = report["help_image_work"]
    lines.extend(["", "## Help Image Work (Whole Session, Not Additive)", "",
                  f"Status: {images['status']}; reason: {images['reason'] or 'none'}."])
    if images["recorded"] is not None:
        lines.extend([f"Dropped observations: {images['dropped_observations']}; "
                      f"saturated: {images['recorded']['saturated']}.", "",
                      "| Recorded counter (may be partial) | Value |", "| --- | ---: |"])
        lines.extend(f"| {key} | {images['recorded'][key]} |" for key in HELP_IMAGE_COUNTERS)
        lines.extend(["", f"Complete rasterization total (ms): {ms(images['rasterization_total_us'])}."])
    for row in report["subjects"]:
        lines.extend(["", f"## Subject {row['subject']} (Process-Local Ordinal)", "",
                      "### CPU Markers", "", "| Marker | Since Rust main (ms) |", "| --- | ---: |"])
        lines.extend(f"| {phase} | {ms(value)} |" for phase, value in row["checkpoints"].items())
        lines.extend(["", "### Inclusive Spans (Do Not Sum)", "",
                      "| Phase / span | Outcome | Begin (ms) | End (ms) | Duration (ms) | Unavailable reason |",
                      "| --- | --- | ---: | ---: | ---: | --- |"])
        for span in row["spans"]:
            lines.append(f"| {span['phase']} / {span['span']} | {span['outcome']} | {ms(span['begin_us'])} | "
                         f"{ms(span['end_us'])} | {ms(span['duration_us'])} | {span['unavailable_reason'] or ''} |")
        if row["boundaries"]:
            lines.extend(["", "### Same-Subject Boundaries (Not Additive)", "",
                          "| Boundary | Elapsed (ms) | Unavailable reason |", "| --- | ---: | --- |"])
            lines.extend(f"| {item['name']} | {ms(item['duration_us'])} | {item['unavailable_reason'] or ''} |"
                         for item in row["boundaries"])
    return "\n".join(lines) + "\n"


def write_summary(source: Path, output: Path) -> dict:
    with source.open("rb") as stream:
        raw = stream.read(MAX_BYTES + 1)
    report = summarize(raw)
    report["summarizer_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    # Reserve a new directory only after validation. Never replace prior evidence.
    output.mkdir(parents=True, exist_ok=False)
    with (output / "trace.json").open("xb") as stream:
        stream.write(raw)
    for name, content in (("summary.json", json.dumps(report, indent=2, sort_keys=True) + "\n"),
                          ("summary.md", markdown(report))):
        with (output / name).open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(content)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = write_summary(args.trace, args.output_dir)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Cannot summarize startup trace: {error}\n")
    print(f"Retained {report['event_count']} events ({report['trace_completeness']}); "
          f"performance not assessed. Output: {args.output_dir}")


if __name__ == "__main__":
    main()
