"""Fast source/tutorial checks only; not Rust execution or scientific acceptance."""

import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class PrimerTraceContractTests(unittest.TestCase):
    def test_tutorial_preparation_is_synthetic_and_digest_bound(self):
        guide = (ROOT / "docs/tutorial/10-04_targeted_nanopore_fastq_panel.md").read_text()
        addendum = guide.split("## Primer Trace Addendum: Synthetic Only", 1)[1]
        code = addendum.split("<<'PY'\n", 1)[1].split("\nPY\n", 1)[0]
        with tempfile.TemporaryDirectory(prefix="primer trace's synthetic ") as scratch:
            result = subprocess.run([sys.executable, "-c", code, scratch], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            root = Path(scratch)
            panel = json.loads((root / "panel.json").read_text())
            digest = panel["content_sha256"]
            panel["content_sha256"] = ""
            encoded = json.dumps(panel, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
            self.assertEqual(digest, "sha256:" + hashlib.sha256(encoded).hexdigest())
            primer = panel["primers"][0]
            self.assertNotEqual(len(primer["sequence"]), 44)
            self.assertEqual(primer["core_orientation"], "sense_5prime")
            self.assertEqual(primer["expected_targets"], [])
            binding = json.loads((root / "binding.json").read_text())
            request = json.loads((root / "trace.json").read_text())
            self.assertEqual(binding, request["primer_panel"])
            self.assertEqual(binding["expected_sha256"], digest)
            fastq = (root / "reads.fastq").read_bytes()
            self.assertEqual(request["files"][0]["expected_sha256"], "sha256:" + hashlib.sha256(fastq).hexdigest())
            lines = fastq.decode().splitlines()
            self.assertEqual(len(lines), 8)
            forward, reverse = lines[1], lines[5]
            rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            full = primer["sequence"]
            self.assertEqual(forward.index(full), 6)
            self.assertEqual(len(reverse) - reverse.index(rc(full)) - len(full), 10)
            self.assertEqual(forward.split(full)[1], rc(reverse[:reverse.index(rc(full))]))
            self.assertTrue(set(lines[3] + lines[7]) == {"I"})
            self.assertFalse(request["export_payload_fasta"])
            self.assertNotIn("reference_panel", request)

    def test_all_routes_share_engine_contract_and_confirmation(self):
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())
        commands = {c["path"]: c for c in glossary["commands"]}
        expected = {
            "rna-reads primer-panel-validate": ["ValidateRnaPrimerPanel"],
            "rna-reads primer-trace": ["TracePrimerReads", "CancelPrimerTrace"],
            "rna-reads primer-trace-show": ["InspectPrimerTrace"],
        }
        engine = (ROOT / "src/engine.rs").read_text()
        registry = (ROOT / "crates/gentle-protocol/src/lib.rs").read_text()
        handlers = (ROOT / "src/engine/ops/operation_handlers.rs").read_text()
        shell = (ROOT / "src/engine_shell.rs").read_text()
        cli = (ROOT / "docs/cli.md").read_text()
        for route, operations in expected.items():
            self.assertEqual(commands[route]["engine_operations"], operations)
            self.assertNotIn("gui-menu", commands[route]["interfaces"])
            self.assertIn(route, cli)
            self.assertIn(f'"{route}"', registry)
            mcp = (ROOT / "src/mcp_server.rs").read_text()
            command_paths = mcp.split('"op" => &[', 1)[1].split("]", 1)[0]
            self.assertIn(f'"{route}"', command_paths)
            for op in operations:
                self.assertRegex(engine, r"\b" + op + r"\s*\{")
                self.assertIn(f'"{op}"', registry)
                self.assertIn("Operation::" + op, handlers)
                self.assertIn(f'primer_trace_descriptor("{op}")', shell)
        raw_mcp = (ROOT / "src/mcp_server.rs").read_text().split("fn op_tool_result(", 1)[1].split("fn workflow_tool_result(", 1)[0]
        self.assertIn('require_confirm_true(&args, "op")', raw_mcp)
        for operation in (op for ops in expected.values() for op in ops):
            self.assertIn("Operation::" + operation, raw_mcp)
        tests = (ROOT / "src/mcp_server.rs").read_text()
        self.assertIn("fn primer_trace_mcp_and_agent_cannot_bypass_artifact_confirmation", tests)
        self.assertIn("fn primer_trace_shell_contract_keeps_observation_read_only_and_writes_confirmed", (ROOT / "src/engine_shell/tests.rs").read_text())

    def test_quality_extension_does_not_filter_existing_nominations(self):
        source = (ROOT / "src/engine/analysis/rna_panel/input.rs").read_text()
        record = source.split("struct FastqRecord", 1)[1].split("#[derive", 1)[0]
        self.assertIn("pub quality: Vec<u8>", record)
        panel = (ROOT / "src/engine/analysis/rna_panel.rs").read_text()
        self.assertNotIn("record.quality", panel)
        legacy = (ROOT / "src/engine/analysis/rna_reads.rs").read_text()
        self.assertNotIn("TracePrimerReads", legacy)
        types = (ROOT / "src/engine/protocol/primer_trace.rs").read_text()
        row = types.split("pub struct RnaPrimerTraceRow", 1)[1].split("#[derive", 1)[0]
        self.assertNotIn("sequence:", row, "read-level inspection must not retain read sequence")

    def test_managed_trace_and_feedback_regressions_are_wired(self):
        host = (ROOT / "src/command_execution.rs").read_text()
        self.assertIn("mod primer_trace_tests;", host)
        tests = (ROOT / "src/command_execution/primer_trace_tests.rs").read_text()
        for name in (
            "primer_trace_managed_partial_reports_keep_outcomes_without_project_lock",
            "primer_trace_managed_compact_result_survives_late_cancel_and_project_replacement",
        ):
            self.assertIn("fn " + name, tests)
        feedback = (ROOT / "src/agent_feedback.rs").read_text()
        self.assertIn("fn primer_trace_feedback_keeps_full_and_compact_publication_outcomes", feedback)
        interpreter = feedback.split("pub fn from_shell_output(", 1)[1].split("pub fn after_sequence_import(", 1)[0]
        for marker in (
            "gentle.command_primer_trace_result.v1",
            "gentle.rna_primer_trace_result_binding.v1",
            "/result/rna_primer_trace",
            "gentle.rna_primer_trace_report.v1",
        ):
            self.assertIn(marker, interpreter)

    def test_source_modules_and_regressions_are_explicit(self):
        for parent, child in [
            ("src/engine/analysis/rna_panel.rs", "rna_panel/primer_trace.rs"),
            ("src/engine/analysis/rna_panel/primer_trace.rs", "primer_trace/matching.rs"),
            ("src/engine/analysis/rna_panel/primer_trace.rs", "primer_trace/tests.rs"),
            ("src/engine/protocol.rs", "protocol/primer_trace.rs"),
        ]:
            text = (ROOT / parent).read_text()
            self.assertIn(f'#[path = "{child}"]', text)
            module = ROOT / parent
            self.assertTrue((module.parent / child).read_text().startswith("//!"))
        types = (ROOT / "src/engine/protocol/primer_trace.rs").read_text()
        for text in (types, types.replace("\n", "\r\n")):
            self.assertIn('rename = "sense_5prime"', text)
            self.assertIn('rename = "antisense_3prime"', text)
            self.assertIn("pub core_and_probe_contiguous: bool", text)
        tests = (ROOT / "src/engine/analysis/rna_panel/primer_trace/tests.rs").read_text()
        for name in (
            "four_orientation_geometries_keep_the_boundary_adjacent_probe",
            "substitutions_insertions_deletions_and_core_limit",
            "handle_no_anchor_short_and_quality_are_not_target_calls",
            "multiple_anchors_do_not_become_unchecked_molecules",
            "reference_ambiguity_is_retained_and_bad_flanks_stay_unknown",
            "policy_binding_budgets_and_cancel_are_explicit",
            "publication_inspection_tamper_and_digest_failure",
            "shared_fastq_keeps_quality_and_raw_lf_crlf_hashes",
        ):
            self.assertIn("fn primer_trace_" + name, tests)

    def test_parity_projection_counts_cover_new_rows_without_new_gaps(self):
        text = (ROOT / "docs/gui_cli_mcp_parity.md").read_text()
        rows = []
        for line in text.splitlines():
            cells = [c.strip() for c in re.split(r"(?<!\\)\|", line.strip("|"))]
            if len(cells) == 9 and cells[1] in {"glossary-command", "engine-operation", "mcp-tool"}:
                rows.append(cells)
        for i, adapter in enumerate(["GUI", "gentle_cli", "MCP", "JS", "Lua", "ClawBio"]):
            match = re.search(r"^\| " + adapter + r" \| (\d+) \| (\d+) \| (\d+) \|$", text, re.M)
            self.assertIsNotNone(match)
            actual = tuple(sum(row[i + 2] == kind for row in rows) for kind in ["prominent", "shell-only", "gap"])
            self.assertEqual(tuple(map(int, match.groups())), actual, adapter)


if __name__ == "__main__":
    unittest.main()
