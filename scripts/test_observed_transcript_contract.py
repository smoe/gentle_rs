"""Fast source/documentation guards, not execution of the Rust structural analyser."""

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "src/engine/analysis/rna_panel/primer_payload"


class ObservedTranscriptContractTests(unittest.TestCase):
    def test_strict_new_wire_records_do_not_extend_phase_two_schema(self):
        source = (ROOT / "src/engine/protocol/observed_transcripts.rs").read_text()
        self.assertEqual(source.count("pub struct "), source.count("#[serde(deny_unknown_fields)]"))
        self.assertNotIn("ObservedTranscript", (ROOT / "src/engine/protocol/primer_payload.rs").read_text())
        result = (ROOT / "src/engine/protocol.rs").read_text()
        for field in ("observed_transcript_structures", "observed_transcript_structures_inspection"):
            self.assertRegex(result, r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s+pub ' + field + r':')
        for path in ("src/engine/ops/operation_handlers.rs", "src/engine/analysis/rna_reads.rs",
                     "src/app.rs", "src/app/tests.rs", "src/main_area_dna/tests.rs"):
            lines = [line.strip() for line in (ROOT / path).read_text().splitlines()]
            for index, line in enumerate(lines):
                if line == "primer_payload_inspection: None,":
                    self.assertEqual(lines[index + 1:index + 3], [
                        "observed_transcript_structures: None,",
                        "observed_transcript_structures_inspection: None,",
                    ], path)

    def test_shared_routes_have_explicit_effect_and_lifecycle_metadata(self):
        commands = {c["path"]: c for c in json.loads((ROOT / "docs/glossary.json").read_text())["commands"]}
        surfaces = ["src/engine.rs", "src/engine_shell.rs", "src/command_execution/rna_artifacts.rs",
                    "src/engine/operation_lifecycle.rs", "crates/gentle-protocol/src/lib.rs"]
        routes = {
            "rna-reads transcript-structures": "SummarizeObservedTranscriptStructures",
            "rna-reads transcript-structures-show": "InspectObservedTranscriptStructures",
            "rna-reads transcript-structures-cancel": "CancelObservedTranscriptStructures",
        }
        for path, operation in routes.items():
            self.assertEqual(commands[path]["engine_operations"], [operation])
            self.assertIn(path, (ROOT / "docs/cli/observed-transcripts.md").read_text())
            self.assertIn(path, (ROOT / "src/mcp_server.rs").read_text())
            for surface in surfaces:
                self.assertIn(operation, (ROOT / surface).read_text(), surface)

    def test_intake_is_saved_evidence_not_an_aligner_or_gene_specific_replay(self):
        source = (BASE / "observed.rs").read_text()
        for call in ("binding::queries", "reference::catalog", "VerifiedReader::open", "kernel::summarize_with_control"):
            self.assertIn(call, source)
        for forbidden in ("Command::", "reqwest", "verify_fasta(", "paf::row(", "TP53", "TP73", "GRCh38", "barcode13"):
            self.assertNotIn(forbidden, source)
        self.assertIn("target_length_unavailable", source)
        self.assertIn("unassessed_records", source)

    def test_saved_result_has_separate_outcomes_and_bounded_inspection(self):
        source = (BASE / "observed.rs").read_text()
        for word in ("unavailable", "control_error", "counts.is_none()", "reader.finish()?", "128_000_000", "4096", "BoundedOutput"):
            self.assertIn(word, source)
        self.assertIn("source/cohort", (ROOT / "docs/observed_transcript_structures.md").read_text())

    def test_source_bound_runtime_regressions_have_an_explicit_execution_scope(self):
        tests = (BASE / "observed_tests.rs").read_text()
        for case in ("keeps_phase2_bytes", "rehashed_duplicate_missing_or_foreign", "incomplete_batch_counts",
                     "late_cancel_and_control_fault", "shared_parser", "raw_lifecycle", "crlf"):
            self.assertIn(case, tests)
        docs = (ROOT / "docs/observed_transcript_structures.md").read_text()
        self.assertIn("## Verification Scope", docs)
        self.assertIn("--lib observed_transcripts -- --test-threads=1", docs)
        self.assertIn("Glen/CI", docs)

    def test_new_source_modules_have_module_summaries(self):
        for path in [BASE / "observed.rs", BASE / "observed_tests.rs", ROOT / "src/engine/protocol/observed_transcripts.rs"]:
            self.assertTrue(path.read_text().startswith("//!"), str(path))

    def test_annotation_routes_name_their_effect_lifecycle_and_shared_parser(self):
        commands = {c["path"]: c for c in json.loads((ROOT / "docs/glossary.json").read_text())["commands"]}
        routes = {
            "rna-reads transcript-annotation-prepare": "PrepareObservedTranscriptAnnotation",
            "rna-reads transcript-annotation-cancel": "CancelObservedTranscriptAnnotation",
        }
        for route, operation in routes.items():
            self.assertEqual(commands[route]["engine_operations"], [operation])
            self.assertIn(route, (ROOT / "docs/cli/observed-transcripts.md").read_text())
            for path in ("src/engine.rs", "src/engine/ops/operation_handlers.rs",
                         "src/engine/operation_lifecycle.rs", "src/engine_shell/command_parsers.rs",
                         "src/command_execution/rna_artifacts.rs", "crates/gentle-protocol/src/lib.rs"):
                self.assertIn(operation, (ROOT / path).read_text(), path)
        protocol = (ROOT / "src/engine/protocol.rs").read_text()
        self.assertRegex(protocol, r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s+pub observed_transcript_annotation:')
        self.assertNotIn("observed_transcript_annotation", (ROOT / "src/engine/protocol/primer_payload.rs").read_text())
        for path in (BASE / "annotation_prepare.rs", BASE / "annotation_prepare/adapter_tests.rs",
                     ROOT / "src/engine/protocol/observed_annotation.rs"):
            self.assertTrue(path.read_text().startswith("//!"), str(path))


if __name__ == "__main__":
    unittest.main()
