"""Fast source/doc guards and an explicit opt-in synthetic CLI/MCP process smoke."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "src/engine/analysis/rna_panel/primer_payload"


class ObservedTranscriptContractTests(unittest.TestCase):
    def test_junction_followup_contract_keeps_independent_geometry_and_unknown_legacy_causes(self):
        wire = (ROOT / "src/engine/protocol/observed_junctions.rs").read_text()
        guide = (ROOT / "docs/observed_splice_junction_assessment.md").read_text()
        for field in ("search_not_declared_complete_records", "query_views_unassessed_records", "gene_unassessed_records"):
            self.assertRegex(wire, r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s+pub ' + field + r': Option<usize>')
            self.assertIn(field, guide)
        for concept in ("outside_gene_locus", "both ends", "any", "not summed", "remain unknown"):
            self.assertIn(concept, guide)
        tests = (BASE / "junction_assessment/tests.rs").read_text()
        for case in ("same_and_cross_contig_off_gene", "exact_shared_end_is_independent", "unassessed_causes_are_overlapping", "legacy_denominator_missing_causes"):
            self.assertIn(case, tests)

    def test_aligner_policy_and_vm_checks_are_explicit_without_reopening_phase_two(self):
        decisions = (ROOT / "docs/decisions.md").read_text()
        for word in ("B-i", "possible-only", "definite_or_possible", "not an aligner-specific"):
            self.assertIn(word, decisions)
        producer = (ROOT / "docs/primer_payload_evidence.md").read_text()
        for word in ("EXACT_VERSION", "ACTUAL_N", "ACTUAL_MIN_CHAIN_SCORE", "template, not a recorded run"):
            self.assertIn(word, producer)
        for language in ("js", "lua"):
            source = (ROOT / f"src/{language}_interface.rs").read_text()
            for route in ("AssessObservedSpliceJunctions", "InspectObservedSpliceJunctions"):
                self.assertIn(route, source)
            self.assertIn(f"junction_assessment_{language}_vm_assess_and_inspect_match_engine_without_project_writes", source)

    @unittest.skipUnless(os.environ.get("GENTLE_CONTRACT_CLI") and os.environ.get("GENTLE_CONTRACT_MCP"),
                         "set both GENTLE_CONTRACT_CLI and GENTLE_CONTRACT_MCP to same-source binaries")
    def test_actual_sam_cli_and_framed_mcp_saved_inspection(self):
        """Hand-crafted temporary SAM/normalized annotation, not provider or prepared evidence."""
        cli = str(Path(os.environ["GENTLE_CONTRACT_CLI"]).resolve())
        mcp = str(Path(os.environ["GENTLE_CONTRACT_MCP"]).resolve())
        with tempfile.TemporaryDirectory(prefix="gentle-sam-public-") as directory:
            root = Path(directory)

            def put(name, value):
                path = root / name
                raw = value.encode() if isinstance(value, str) else json.dumps(value, sort_keys=True).encode()
                path.write_bytes(raw)
                return {"path": str(path), "sha256": "sha256:" + hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}

            declared = "sha256:" + "1" * 64
            sam = put("query's input.sam", "@HD\tVN:1.6\n@SQ\tSN:ctg_a\tLN:1000\nq1\t0\tctg_a\t101\t60\t10M10N20M\t*\t0\t0\t*\t*\n")
            dictionary = put("dictionary.json", {
                "schema": "gentle.reference_contig_dictionary.v1", "reference_id": "synthetic_reference",
                "assembly": "synthetic", "reference_fasta_sha256": declared,
                "contigs": [{"name": "ctg_a", "length": 1000}],
            })
            manifest = put("query units.json", {
                "schema": "gentle.observed_sam_query_units.v1", "sam_sha256": sam["sha256"],
                "units": [{"unit_id": "synthetic_unit", "unclipped_axis_length": 30,
                           "rows": [{"ordinal": 1, "qname": "q1", "read_group": None}], "ends": [[], []],
                           "masks": {"coordinate_convention": "zero_based_half_open_unclipped_query",
                                     "unclipped_axis_length": 30, "intervals": [], "source_sha256": sam["sha256"],
                                     "derivation": "Hand-crafted fixture with no excluded bases"},
                           "search": {"search_complete": True, "views_complete": True, "limitations": [],
                                      "producer": "hand-crafted fixture", "version": "1", "arguments": [],
                                      "selection_scope": "only supplied synthetic alternatives"}}],
            })
            annotation = put("annotation.json", {
                "schema": "gentle.primer_payload_annotation.v1", "reference_id": "synthetic_reference",
                "assembly": "synthetic", "release": "synthetic_snapshot", "reference_fasta_sha256": declared,
                "coordinate_convention": "zero_based_half_open", "transcripts": [{
                    "transcript_id": "synthetic_transcript", "gene_ids": ["synthetic_gene"],
                    "target_name": "ctg_a", "target_length": 1000, "genomic_contig": "ctg_a",
                    "genomic_strand": "+", "genomic_exons": [[100, 110], [120, 140]],
                }, {
                    "transcript_id": "synthetic_reverse_transcript", "gene_ids": ["synthetic_other_gene"],
                    "target_name": "ctg_a", "target_length": 1000, "genomic_contig": "ctg_a",
                    "genomic_strand": "-", "genomic_exons": [[500, 510], [520, 540]],
                }],
            })
            selector = {"source_id": "source_a", "gene_id": "synthetic_gene"}
            request = {
                "schema": "gentle.observed_sam_structures_request.v1", "sam": sam, "dictionary": dictionary,
                "manifest": manifest, "conversion": {"kind": "direct_sam", "producer": "synthetic process smoke",
                    "version": "1", "arguments": [], "selection_scope": "one hand-crafted unit"},
                "reference_id": "synthetic_reference", "reference": {"assembly": "synthetic", "fasta_sha256": declared},
                "reference_gene": selector, "comparison_genes": [], "annotation_sources": [{
                    "kind": "normalized", "selector": selector, "provider": "synthetic", "annotation": annotation}],
                "version_policy": "strict16", "unpadded_reference": True, "required_flank_bases": 2,
                "text_budgets": {"max_input_bytes": 1000000, "max_line_bytes": 65536, "max_lines": 100,
                    "max_rows": 100, "max_dictionary_entries": 100, "max_header_ids": 10, "max_fields_per_line": 64,
                    "max_operations": 10000, "max_sa_entries": 100, "max_retained_bytes": 4194304, "max_work": 100000},
                "kernel_budgets": {"max_queries": 100, "max_records": 100, "max_segments_per_query": 100,
                    "max_operations": 10000, "max_skips_per_chain": 100, "max_chain_families": 100,
                    "max_alternatives_per_record": 10, "max_annotation_rows": 100,
                    "max_model_comparisons_per_family": 100, "max_total_work": 1000000},
                "budgets": {"max_input_bytes": 1000000, "max_annotation_bytes": 65536, "max_output_bytes": 2097152,
                    "max_output_row_bytes": 65536, "max_units": 100, "max_retained_bytes": 8388608,
                    "max_projection_work": 1000000}, "output_dir": str(root / "gene's evidence"), "cancel_path": None,
            }
            request_file = put("request's SAM.json", request)
            state = root / "must not be created.json"

            def run_cli(*args):
                process = subprocess.run([cli, "--state", str(state), *args], cwd=root,
                                         capture_output=True, text=True, timeout=60)
                self.assertEqual(process.returncode, 0, process.stderr)
                return json.loads(process.stdout)

            def run_mcp(operation):
                call = {"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {
                    "name": "op", "arguments": {"confirm": True, "operation": operation}}}
                body = json.dumps(call).encode()
                framed = f"Content-Length: {len(body)}\r\n\r\n".encode() + body
                process = subprocess.run([mcp, "--state", str(state)], cwd=root, input=framed,
                                         capture_output=True, timeout=60)
                self.assertEqual(process.returncode, 0, process.stderr.decode(errors="replace"))
                header, body = process.stdout.split(b"\r\n\r\n", 1)
                self.assertEqual(int(header.split(b":", 1)[1]), len(body))
                result = json.loads(body)["result"]
                self.assertFalse(result.get("isError"), result)
                return result["structuredContent"]

            def saved_data(report):
                return json.loads((Path(report["resolved_output_dir"]) / "artifacts/analysis.json").read_bytes())

            result = run_cli("rna-reads", "transcript-structures-from-sam", "@" + request_file["path"])
            report = result["result"]["observed_sam_structures"]
            self.assertEqual(report["status"], "completed")
            self.assertEqual(report["counts"]["query_units"], 1)
            self.assertEqual(report["coordinate_authority"], "declared")
            self.assertFalse(state.exists())
            binding = {"schema": "gentle.observed_sam_structures_result_binding.v1",
                "request_sha256": report["request_sha256"], "report_content_sha256": report["content_sha256"],
                "report_path": str(Path(report["resolved_output_dir"]) / "report.json"), "outcome": "completed"}
            first = saved_data(report)["analysis"]
            other_request = json.loads(json.dumps(request))
            other_request["reference_gene"]["gene_id"] = "synthetic_other_gene"
            other_request["annotation_sources"][0]["selector"]["gene_id"] = "synthetic_other_gene"
            other_request["output_dir"] = str(root / "reverse gene's evidence")
            other_report = run_cli("rna-reads", "transcript-structures-from-sam", json.dumps(other_request))["result"]["observed_sam_structures"]
            other = saved_data(other_report)["analysis"]
            self.assertEqual(other["gene_resolutions"][0]["comparison"]["gene_id"], "synthetic_other_gene")
            self.assertEqual(other["gene_resolutions"][0]["comparison"]["resolution"], "resolved")
            self.assertNotEqual(first["context_sha256"], other["context_sha256"])
            self.assertEqual(first["chain_families"], other["chain_families"])
            self.assertNotEqual(first["observations"][0]["comparisons"], other["observations"][0]["comparisons"])

            second_catalog = json.loads(Path(annotation["path"]).read_bytes())
            second_catalog["release"] = "synthetic_second_snapshot"
            second_annotation = put("second annotation.json", second_catalog)
            multi_request = json.loads(json.dumps(other_request))
            multi_request["reference_gene"]["source_id"] = "source_b"
            multi_request["comparison_genes"] = [selector]
            multi_request["annotation_sources"] = [
                {"kind": "normalized", "selector": multi_request["reference_gene"], "provider": "synthetic",
                 "annotation": second_annotation}, request["annotation_sources"][0]]
            multi_request["output_dir"] = str(root / "multi-source evidence")
            multi = run_mcp({"SummarizeObservedTranscriptStructuresFromSam": {"request": multi_request}})
            multi_binding = multi["observed_sam_result"]
            self.assertEqual(multi_binding["outcome"], "completed")
            self.assertEqual([r["source_id"] for r in multi["gene_resolutions"]], ["source_a", "source_b"])
            self.assertEqual(multi["gene_resolutions"][1]["comparison"]["gene_id"], "synthetic_other_gene")
            for file in [sam, dictionary, manifest, annotation, second_annotation]:
                Path(file["path"]).unlink()
            inspected = run_cli("rna-reads", "transcript-structures-from-sam-show", binding["report_path"],
                "--table", "junction_catalog.tsv", "--rows", "1", "--verify-artifacts", "--binding", json.dumps(binding))
            rows = inspected["result"]["observed_sam_inspection"]["rows"]
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["counts"]["definite"]["qualified_only"], 1)
            inspected = run_mcp({"InspectObservedTranscriptStructuresFromSam": {
                "path": multi_binding["report_path"], "expected_result": multi_binding,
                "table": "junction_catalog.tsv", "rows": 1, "verify_artifacts": True}})
            multi_rows = inspected["result"]["observed_sam_inspection"]["rows"]
            self.assertEqual(len(multi_rows), 1)
            self.assertEqual(multi_rows[0]["counts"]["definite"]["qualified_only"], 1)
            self.assertFalse(state.exists())

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
                    self.assertEqual(lines[index + 1:index + 5], [
                        "observed_transcript_structures: None,",
                        "observed_sam_structures: None,",
                        "observed_sam_inspection: None,",
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
                         "src/command_execution/rna_artifacts.rs", "crates/gentle-protocol/src/lib.rs",
                         "src/mcp_server.rs"):
                self.assertIn(operation, (ROOT / path).read_text(), path)
        protocol = (ROOT / "src/engine/protocol.rs").read_text()
        self.assertRegex(protocol, r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s+pub observed_transcript_annotation:')
        self.assertNotIn("observed_transcript_annotation", (ROOT / "src/engine/protocol/primer_payload.rs").read_text())
        for path in (BASE / "annotation_prepare.rs", BASE / "annotation_prepare/adapter_tests.rs",
                     ROOT / "src/engine/protocol/observed_annotation.rs"):
            self.assertTrue(path.read_text().startswith("//!"), str(path))

    def test_sam_decoder_has_no_io_and_public_ingress_is_separate(self):
        root = ROOT / "src/observed_transcript_structures"
        for name in ("sam_text.rs", "sam_text_tests.rs"):
            self.assertTrue((root / name).read_text().startswith("//!"), name)
        module = (root / "sam_text.rs").read_text()
        for forbidden in ("std::fs", "fs::", "Command::", "reqwest"):
            self.assertNotIn(forbidden, module)
        self.assertIn("sam::validate_decoded_row", module)
        self.assertIn("mod sam_text;", (ROOT / "src/observed_transcript_structures.rs").read_text())
        commands = {c["path"]: c for c in json.loads((ROOT / "docs/glossary.json").read_bytes())["commands"]}
        for suffix, op in (("", "Summarize"), ("-show", "Inspect"), ("-cancel", "Cancel")):
            name = op + "ObservedTranscriptStructuresFromSam"
            route = "rna-reads transcript-structures-from-sam" + suffix
            self.assertEqual(commands[route]["engine_operations"], [name])
            for path in ("src/engine.rs", "src/engine/operation_lifecycle.rs", "src/mcp_server.rs",
                         "src/command_execution/rna_artifacts.rs", "src/engine_shell/command_parsers.rs"):
                self.assertIn(name, (ROOT / path).read_text(), path)
            self.assertIn(route, (ROOT / "docs/cli/observed-transcripts.md").read_text())
        self.assertIn("verified_json", (root / "sam_intake.rs").read_text())


if __name__ == "__main__":
    unittest.main()
