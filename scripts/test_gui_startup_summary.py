"""Hand-crafted CPU traces; no GUI launch, private data or timing assertions.

Recreate with python3 -m unittest scripts.test_gui_startup_summary -v. The
synthetic events exercise the real CLI/file boundary, raw hashes and span math;
they do not establish native performance or replace Glen's recorded traces.
"""

import copy
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

from scripts import gui_startup_summary as summary

ROOT = Path(__file__).resolve().parents[1]


def event(time, phase, kind="checkpoint", subject=0, span=None):
    result = {"elapsed_us": time, "subject": subject, "phase": phase, "kind": kind}
    if span is not None:
        result["span"] = span
    return result


def fixture():
    return {
        "schema": "gentle.gui_startup_trace.v1", "source_revision": "a" * 40,
        "git_commit": "a" * 40, "os": "synthetic", "architecture": "synthetic",
        "clock": "monotonic microseconds since Rust main entry; excludes process loader",
        "scope": "CPU checkpoints only; not compositor presentation or native usability acceptance",
        "debug_assertions": False, "gui_test_support": False, "gui_profiler": False,
        "event_limit": 512, "dropped_events": 0,
        "events": [
            event(0, "main_entered"),
            event(10, "native_run", "begin", span=1),
            event(20, "app_initialize", "begin", span=2),
            event(25, "project_load", "begin", span=3),
            event(30, "project_read_decode", "begin", span=4),
            event(40, "project_read_decode", "completed", span=4),
            event(45, "project_install", "begin", span=5),
            event(50, "project_install", "completed", span=5),
            event(60, "project_load", "completed", span=3),
            event(70, "app_initialize", "completed", span=2),
            event(80, "root_first_frame"),
            event(90, "root_workspace_frame"),
            event(100, "dna_open_dispatch", "begin", span=6),
            event(105, "dna_construct", "begin", subject=7, span=8),
            event(110, "dna_worker_scheduled", subject=7),
            event(115, "dna_placeholder_construct", "begin", subject=7, span=9),
            event(120, "dna_engine_read_lock", "begin", subject=7, span=10),
            event(125, "dna_placeholder_construct", "completed", subject=7, span=9),
            event(130, "dna_construct", "completed", subject=7, span=8),
            event(135, "dna_open_dispatch", "completed", span=6),
            event(140, "dna_engine_read_lock", "completed", subject=7, span=10),
            event(145, "dna_sequence_clone", "begin", subject=7, span=11),
            event(150, "dna_sequence_clone", "completed", subject=7, span=11),
            event(155, "dna_worker_result", subject=7),
            event(170, "dna_hydrate", "begin", subject=7, span=12),
            event(190, "dna_hydrate", "completed", subject=7, span=12),
            event(210, "dna_native_content_frame", subject=7),
            event(10_000, "native_run", "completed", span=1),
        ],
    }


def encode(trace):
    return (json.dumps(trace, indent=2) + "\n").encode()


def help_fixture():
    trace = fixture()
    trace["help_image_work"] = dict.fromkeys(summary.HELP_IMAGE_COUNTERS, 0)
    trace["help_image_work"].update(svg_references=7, cache_hits=2, preparation_failures=1,
                                    rasterization_attempts=4, rasterization_completed=4,
                                    rasterization_failures=1, rasterization_us=35, saturated=False)
    trace["dropped_help_image_observations"] = 0
    trace["events"].extend([
        event(21, "help_preparation", "begin", span=30),
        event(65, "help_preparation", "completed", span=30),
    ])
    for i, phase in enumerate(("help_manuals", "help_shell_reference", "help_tutorial_discovery",
                               "help_tutorial_selected_load", "help_open", "help_tutorial_open",
                               "help_tutorial_menu_discovery", "help_tutorial_switch")):
        start = 22 + i * 10 if i < 4 else 300 + i * 10
        trace["events"].extend([event(start, phase, "begin", span=31 + i),
                                event(start + 5, phase, "completed", span=31 + i)])
    trace["events"].sort(key=lambda row: row["elapsed_us"])
    return trace


class GuiStartupSummaryTests(unittest.TestCase):
    def analyze(self, trace=None):
        return summary.summarize(encode(fixture() if trace is None else trace))

    def test_help_extension_is_optional_not_zero_in_old_traces(self):
        self.assertEqual(self.analyze()["help_image_work"], {
            "status": "unavailable", "reason": "producer_did_not_record", "recorded": None,
        })
        trace = help_fixture()
        trace["help_image_work"] = dict.fromkeys(summary.HELP_IMAGE_COUNTERS, 0)
        trace["help_image_work"]["saturated"] = False
        images = self.analyze(trace)["help_image_work"]
        self.assertEqual(images["status"], "no_recorded_loss")
        self.assertEqual(images["rasterization_total_us"], 0)

    def test_help_subphases_keep_nested_intervals_and_first_use_separate(self):
        report = self.analyze(help_fixture())
        spans = {row["phase"]: row for row in report["subjects"][0]["spans"]}
        parent = spans["help_preparation"]
        self.assertEqual(parent["duration_us"], 44)
        for phase in ("help_manuals", "help_shell_reference", "help_tutorial_discovery", "help_tutorial_selected_load"):
            child = spans[phase]
            self.assertLess(parent["begin_us"], child["begin_us"])
            self.assertLess(child["end_us"], parent["end_us"])
            self.assertEqual(child["duration_us"], 5)
        self.assertGreater(spans["help_open"]["begin_us"], parent["end_us"])
        images = report["help_image_work"]
        self.assertEqual(images["rasterization_total_us"], 35)
        self.assertEqual(images["recorded"]["rasterization_failures"], 1)
        self.assertIn("Whole Session, Not Additive", summary.markdown(report))

    def test_help_losses_and_unfinished_work_are_not_complete_totals(self):
        for change, reason in (({"rasterization_completed": 3}, "unfinished_work"),
                               ({"saturated": True}, "saturated"),
                               ({"svg_references": 0}, "dropped_observations")):
            trace = help_fixture()
            trace["help_image_work"].update(change)
            if reason == "dropped_observations":
                trace["dropped_help_image_observations"] = 1
            with self.subTest(reason=reason):
                report = self.analyze(trace)
                images = report["help_image_work"]
                self.assertEqual(images["status"], "incomplete")
                self.assertEqual(images["reason"], reason)
                self.assertIsNone(images["rasterization_total_us"])
                self.assertIsNotNone(report["subjects"][0]["spans"][0]["duration_us"])
        trace = help_fixture()
        trace["dropped_events"] = 1
        self.assertEqual(self.analyze(trace)["help_image_work"]["rasterization_total_us"], 35)

    def test_help_extension_rejects_bad_types_and_unexplained_inconsistency(self):
        for key, value in (("svg_references", -1), ("cache_hits", True), ("rasterization_us", 0.5),
                           ("saturated", 1), ("rasterization_completed", 5),
                           ("rasterization_failures", 5), ("cache_hits", 8), ("path", "private")):
            trace = help_fixture()
            trace["help_image_work"][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.analyze(trace)
        for field in ("help_image_work", "dropped_help_image_observations"):
            trace = help_fixture()
            del trace[field]
            with self.subTest(missing=field), self.assertRaises(ValueError):
                self.analyze(trace)

    def test_span_pairing_and_cpu_only_limits(self):
        report = self.analyze()
        self.assertEqual(report["trace_completeness"], "no_recorded_loss")
        self.assertFalse(report["native_gui_measured"])
        self.assertEqual(report["performance_verdict"], "not_assessed")
        root, window = report["subjects"]
        self.assertEqual(root["checkpoints"]["main_entered"], 0)
        spans = {span["phase"]: span for span in root["spans"]}
        self.assertEqual(spans["project_load"]["duration_us"], 35)
        self.assertEqual(spans["project_read_decode"]["duration_us"], 10)
        self.assertEqual(spans["project_install"]["duration_us"], 5)
        self.assertEqual(spans["native_run"]["duration_us"], 9990)
        self.assertEqual(root["boundaries"][0]["duration_us"], 10)
        intervals = {row["name"]: row for row in window["boundaries"]}
        self.assertEqual(intervals["worker_schedule_to_lock_attempt"]["duration_us"], 10)
        self.assertEqual(intervals["worker_result_to_hydration"]["duration_us"], 15)
        self.assertEqual(intervals["hydration_to_native_content_cpu"]["duration_us"], 20)
        self.assertEqual(intervals["construct_to_native_content_cpu"]["duration_us"], 105)
        self.assertIsNone(intervals["construct_to_embedded_content_cpu"]["duration_us"])
        text = summary.markdown(report)
        self.assertIn("Do Not Sum", text)
        self.assertIn("not a startup duration", text)
        self.assertIn("producer-reported", text)
        self.assertNotIn("total_duration", json.dumps(report))

    def test_concurrent_subjects_are_not_joined_by_phase(self):
        trace = fixture()
        trace["events"].extend([
            event(160, "dna_construct", "begin", subject=20, span=21),
            event(160, "dna_construct", "completed", subject=20, span=21),
            event(220, "dna_embedded_content_frame", subject=20),
        ])
        trace["events"].sort(key=lambda row: row["elapsed_us"])
        report = self.analyze(trace)
        window = report["subjects"][2]
        self.assertEqual(window["spans"][0]["duration_us"], 0)  # A real observed zero.
        intervals = {row["name"]: row for row in window["boundaries"]}
        self.assertEqual(intervals["construct_to_embedded_content_cpu"]["duration_us"], 60)
        self.assertIsNone(intervals["hydration_to_embedded_content_cpu"]["duration_us"])
        self.assertNotIn("dna_native_content_frame", window["checkpoints"])

    def test_any_dropped_event_suppresses_all_derived_durations(self):
        trace = fixture()
        trace["dropped_events"] = 1
        report = self.analyze(trace)
        self.assertEqual(report["trace_completeness"], "incomplete")
        for subject in report["subjects"]:
            self.assertTrue(all(span["duration_us"] is None for span in subject["spans"]))
            self.assertTrue(all(row["duration_us"] is None for row in subject["boundaries"]))
        self.assertEqual(report["subjects"][1]["checkpoints"]["dna_native_content_frame"], 210)

    def test_orphan_endpoints_are_explicit_not_zero(self):
        for kind, reason in (("begin", "missing_begin"), ("completed", "missing_end")):
            trace = fixture()
            trace["events"] = [row for row in trace["events"]
                               if not (row.get("span") == 12 and row["kind"] == kind)]
            with self.subTest(kind=kind):
                report = self.analyze(trace)
                span = report["subjects"][1]["spans"][-1]
                self.assertIsNone(span["duration_us"])
                self.assertEqual(span["unavailable_reason"], reason)
                self.assertEqual(report["trace_completeness"], "incomplete")

    def test_failed_and_interrupted_work_is_not_completed_startup(self):
        for kind in ("failed", "interrupted"):
            trace = fixture()
            next(row for row in trace["events"] if row.get("span") == 12 and row["kind"] == "completed")["kind"] = kind
            report = self.analyze(trace)
            span = report["subjects"][1]["spans"][-1]
            self.assertEqual(span["outcome"], kind)
            self.assertEqual(span["duration_us"], 20)  # Time spent failing is still observed.
            interval = next(row for row in report["subjects"][1]["boundaries"]
                            if row["name"] == "hydration_to_native_content_cpu")
            self.assertIsNone(interval["duration_us"])
            self.assertEqual(interval["unavailable_reason"], "unsuccessful_or_incomplete_span")
            self.assertEqual(report["performance_verdict"], "not_assessed")

    def test_missing_content_and_empty_trace_remain_unavailable(self):
        trace = fixture()
        trace["events"] = [row for row in trace["events"] if row["phase"] not in
                           ("dna_native_content_frame", "root_workspace_frame")]
        trace["events"].insert(-1, event(215, "dna_load_failed", subject=7))
        report = self.analyze(trace)
        self.assertIn("no DNA content", " ".join(report["warnings"]))
        self.assertIn("dna_load_failed", " ".join(report["warnings"]))
        self.assertIn("No root workspace", " ".join(report["warnings"]))
        trace["events"] = []
        report = self.analyze(trace)
        self.assertEqual(report["trace_completeness"], "incomplete")
        self.assertEqual(report["subjects"], [])

    def test_repeated_spans_do_not_invent_a_boundary_pair(self):
        trace = fixture()
        trace["events"][-1:-1] = [
            event(230, "dna_hydrate", "begin", subject=7, span=25),
            event(240, "dna_hydrate", "completed", subject=7, span=25),
        ]
        report = self.analyze(trace)
        row = next(row for row in report["subjects"][1]["boundaries"]
                   if row["name"] == "worker_result_to_hydration")
        self.assertIsNone(row["duration_us"])
        self.assertEqual(row["unavailable_reason"], "ambiguous_endpoint")

    def test_invalid_trace_structure_fails_closed(self):
        invalid = []
        for key, value in (("schema", "future"), ("events", {}), ("dropped_events", -1),
                           ("event_limit", 513), ("debug_assertions", "false"), ("os", "\n"),
                           ("event_limit", 1), ("dropped_events", True)):
            trace = fixture()
            trace[key] = value
            invalid.append(trace)
        for key, value in (("elapsed_us", -1), ("subject", True), ("elapsed_us", 0.5),
                           ("subject", 2**64), ("kind", "invented"), ("phase", "unknown"),
                           ("phase", []), ("span", 1)):
            trace = fixture()
            trace["events"][0][key] = value
            invalid.append(trace)
        for key, value in (("span", 0), ("span", True), ("span", None), ("phase", "root_first_frame")):
            trace = fixture()
            trace["events"][1][key] = value
            invalid.append(trace)
        trace = fixture()
        trace["events"][2]["elapsed_us"] = 0
        invalid.append(trace)
        trace = fixture()
        trace["events"].insert(1, copy.deepcopy(trace["events"][0]))
        invalid.append(trace)
        trace = fixture()
        trace["events"].insert(2, copy.deepcopy(trace["events"][1]))
        invalid.append(trace)
        trace = fixture()
        trace["events"][-1]["subject"] = 7
        invalid.append(trace)
        trace = fixture()
        trace["events"][-1]["phase"] = "project_load"
        invalid.append(trace)
        trace = fixture()
        trace["events"][1]["kind"] = "completed"
        trace["events"][-1]["kind"] = "begin"
        invalid.append(trace)
        for index, trace in enumerate(invalid):
            with self.subTest(index=index), self.assertRaises(ValueError):
                self.analyze(trace)

    def test_invalid_json_and_oversized_input(self):
        for raw in (b"[]", b"{", b"\xff", b"[" * 2000, b" " * (summary.MAX_BYTES + 1),
                    encode(fixture()).replace(b'"dropped_events": 0', b'"dropped_events": 0, "dropped_events": 1')):
            with self.subTest(raw_size=len(raw)), self.assertRaises(ValueError):
                summary.summarize(raw)

    def test_inverted_boundary_is_unavailable_not_negative_duration(self):
        trace = fixture()
        next(row for row in trace["events"] if row["phase"] == "dna_worker_result")["elapsed_us"] = 195
        trace["events"].sort(key=lambda row: row["elapsed_us"])
        report = self.analyze(trace)
        row = next(row for row in report["subjects"][1]["boundaries"]
                   if row["name"] == "worker_result_to_hydration")
        self.assertIsNone(row["duration_us"])
        self.assertEqual(row["unavailable_reason"], "endpoint_order")

    def test_unknown_metadata_stays_in_raw_input_only(self):
        trace = fixture()
        trace["private_note"] = "not part of the producer contract"
        trace["source_revision"] = "<script>|revision`"
        report = self.analyze(trace)
        self.assertNotIn("private_note", report["producer"])
        self.assertNotIn("<script>", summary.markdown(report))
        self.assertIn("&lt;script&gt;&#124;revision&#96;", summary.markdown(report))

    def test_phase_vocabulary_matches_rust_producer(self):
        source = (ROOT / "src/gui_profiler/startup_trace.rs").read_text(encoding="utf-8")
        enum = source.split("pub enum Phase {", 1)[1].split("}", 1)[0]
        phases = {re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()
                  for name in re.findall(r"\b([A-Z][A-Za-z]+),", enum)}
        self.assertEqual(phases, summary.CHECKPOINTS | summary.SPANS)
        self.assertIn(f"const MAX_EVENTS: usize = {summary.MAX_EVENTS};", source)

    def test_cli_retains_raw_lf_crlf_hashes_and_never_overwrites(self):
        with tempfile.TemporaryDirectory(prefix="startup space ") as directory:
            root = Path(directory)
            hashes = []
            for name, raw in (("LF", encode(help_fixture())), ("CRLF", encode(help_fixture()).replace(b"\n", b"\r\n"))):
                source = root / f"{name} input.json"
                output = root / f"{name} report"
                source.write_bytes(raw)
                command = [sys.executable, str(ROOT / "scripts/gui_startup_summary.py"),
                           str(source), "--output-dir", str(output)]
                result = subprocess.run(command, capture_output=True, timeout=10)
                self.assertEqual(result.returncode, 0, result.stderr)
                record = json.loads((output / "summary.json").read_bytes())
                hashes.append(record["input"]["sha256"])
                self.assertEqual(hashes[-1], hashlib.sha256(raw).hexdigest())
                self.assertEqual(record["summarizer_sha256"], hashlib.sha256(Path(summary.__file__).read_bytes()).hexdigest())
                self.assertEqual((output / "trace.json").read_bytes(), raw)
                self.assertEqual(source.read_bytes(), raw)
                saved = {path.name: path.read_bytes() for path in output.iterdir()}
                self.assertNotIn(b"\r", saved["summary.md"])
                self.assertNotIn(b"\r", saved["summary.json"])
                rerun = subprocess.run(command, capture_output=True, timeout=10)
                self.assertEqual(rerun.returncode, 2)
                self.assertEqual(saved, {path.name: path.read_bytes() for path in output.iterdir()})
            self.assertNotEqual(*hashes)

    def test_cli_invalid_input_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "malformed.json"
            source.write_bytes(b"{}")
            output = root / "not-created"
            result = subprocess.run([sys.executable, str(ROOT / "scripts/gui_startup_summary.py"),
                                     str(source), "--output-dir", str(output)], capture_output=True, timeout=10)
            self.assertEqual(result.returncode, 2)
            self.assertFalse(output.exists())
            self.assertEqual(source.read_bytes(), b"{}")


if __name__ == "__main__":
    unittest.main()
