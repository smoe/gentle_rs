"""Source/metadata checks for 10.04; these do not execute GENtle or certify biology."""

import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT / "docs/tutorial/10-04_targeted_nanopore_fastq_panel.md"
SOURCE = ROOT / "docs/tutorial/sources/10-04_targeted_nanopore_fastq_panel.json"
HELPER = ROOT / "scripts/prepare_rna_panel_tutorial.py"


class RnaPanelTutorialTests(unittest.TestCase):
    def test_table_verifier_uses_the_declared_streaming_hash_backend(self):
        dependencies = tomllib.loads((ROOT / "Cargo.toml").read_text())["dependencies"]
        self.assertIn("ring", dependencies)
        source = (ROOT / "src/engine/analysis/rna_panel/verification.rs").read_text()
        for text in (source, source.replace("\n", "\r\n")):
            self.assertIn("use ring::digest::{Context, SHA256};", text)
            self.assertIn("Context::new(&SHA256)", text)
            self.assertIn("hasher.finish().as_ref()", text)
            self.assertNotIn("use sha2::", text)
            self.assertIn("file.take(artifact.bytes + 1)", text)
            self.assertIn("panel_verify_artifact_hash_or_size", text)

    def test_agent_confirmation_uses_parsed_commands_before_execution(self):
        source = (ROOT / "src/agent_execution.rs").read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            guard = text.split("pub fn execute_agent_plan_candidate_with_feedback(", 1)[1].split(
                "status = AgentExecutionStatus::Failed;", 1)[0]
            self.assertIn("parsed.is_state_mutating()", guard)
            self.assertIn("selected.requires_confirmation", guard)
            self.assertIn("explicit --confirm", guard)
            self.assertRegex(guard, r"AgentPlanCandidateKind::Op => Some\(ShellCommand::Op")
            self.assertRegex(guard, r"AgentPlanCandidateKind::Workflow => Some\(ShellCommand::Workflow")
        self.assertIn("fn rna_artifact_writes_require_confirm_even_when_candidate_flags_are_false()", source)
        self.assertIn("fn read_only_shell_and_rna_inspection_do_not_require_confirm()", source)

    def test_demo_saves_reconciled_state_before_the_mcp_unchanged_state_check(self):
        source = (ROOT / "src/engine/analysis/rna_panel/tools.rs").read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            self.assertRegex(text, r"GentleEngine::from_state\(state\)\s*\.state\(\)\s*\.save_to_path\(")
        mcp = (ROOT / "src/mcp_server.rs").read_text(encoding="utf-8")
        self.assertIn("the synthetic starter must already contain reconciled lineage and containers", mcp)

    def test_path_loaded_panel_children_resolve_to_their_tracked_sources(self):
        source_path = ROOT / "src/engine/analysis/rna_panel.rs"
        source = source_path.read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            for child in ("input", "control", "verification", "tools", "tests", "summary", "summary_tests", "project"):
                match = re.search(r'#\[path = "([^"]+)"\]\s*mod ' + child + r';', text)
                self.assertIsNotNone(match, f"path-loaded panel needs an explicit {child} path")
                self.assertEqual(match.group(1), f"rna_panel/{child}.rs")
                self.assertTrue((source_path.parent / match.group(1)).is_file())

    def test_native_preparation_reuses_gene_sets_and_shared_routes(self):
        guide = GUIDE.read_text(encoding="utf-8")
        for phrase in ("Gene sets...", "Resolve new set > Explicit members",
                       "Bind loaded gene-set members", "Keep panel in project graph",
                       "parent_panel_sha256", "expected_source_record_sha256",
                       "register_in_project=true", "20", "ui open rna-read-panel"):
            self.assertIn(phrase, guide)
        source = (ROOT / "src/app/rna_panel_ui.rs").read_text(encoding="utf-8")
        for phrase in ("ShellCommand::RnaReadPanelSources", "Operation::BuildRnaReadPanel",
                       "execute_on_engine_snapshot", "try_open_rna_panel_graph_node",
                       "rna_panel_gene_set_registration_graph_save_reload_and_undo_agree",
                       "rna_panel_gene_set_binding_rejects_stale_unresolved_and_relabelled_sources"):
            self.assertIn(phrase, source)
        project = (ROOT / "src/engine/analysis/rna_panel/project.rs").read_text(encoding="utf-8")
        for phrase in ("get_gene_set_resolution_artifact", "resolution_sha256",
                       "panel_gene_set_ambiguous_or_missing_source", "panel_gene_set_changed"):
            self.assertIn(phrase, project)
        bridge = (ROOT / "src/agent_bridge.rs").read_text(encoding="utf-8")
        self.assertIn("--list-sources [GENE_SET_ID]", bridge)
        self.assertIn('Builds remain execution="ask"', bridge)
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())["commands"]
        entry, = [row for row in glossary if row["path"] == "rna-reads panel-build"]
        self.assertIn("gui-menu", entry["interfaces"])

    def test_source_catalog_discovery_and_unique_placement(self):
        source = json.loads(SOURCE.read_text())
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_text())
        entry, = [e for e in catalog["entries"] if e["id"] == source["id"]]
        self.assertEqual(entry["title"], source["title"])
        self.assertEqual(entry["notes"], source["catalog"]["notes"])
        self.assertEqual(entry["produces"], source["catalog"]["produces"])
        self.assertEqual(entry["decimal_id"], "10.04")
        self.assertEqual(entry["review_status"], "unreviewed")
        for path in SOURCE.parent.glob("*.json"):
            other = json.loads(path.read_text())
            if other.get("id") != source["id"] and "catalog" in other:
                self.assertNotEqual(other["catalog"]["order"], source["catalog"]["order"])
        self.assertIn("10-04_targeted_nanopore_fastq_panel.md",
                      (ROOT / "docs/tutorial/README.md").read_text())

    def test_literal_bash_is_syntactically_valid_in_lf_and_crlf(self):
        payload = GUIDE.read_bytes()
        for raw in (payload.replace(b"\r\n", b"\n"), payload.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")):
            text = raw.decode().replace("\r\n", "\n")
            blocks = re.findall(r"```bash\n(.*?)```", text, re.S)
            self.assertGreaterEqual(len(blocks), 5)
            for block in blocks:
                result = subprocess.run(["bash", "-n"], input=block, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_request_helper_binds_explicit_samples_without_opening_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "request with spaces.json"
            # Intentionally nonexistent references/reads: this is metadata preparation only.
            args = [sys.executable, str(HELPER), "request", "--panel", str(root / "panel.json"),
                    "--md5-manifest", str(root / "md5_checksums.txt"), "--output-dir", str(root / "run"),
                    "--output", str(output), "--declared-complete", "--write-locus-fastas",
                    "--input", "26R266-1", "barcode13", "custom_only", str(root / "renamed folder/a.fastq.gz"),
                    "--input", "26R266-2", "barcode14", "custom_only", str(root / "b.fastq.gz")]
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            request = json.loads(output.read_text())
            self.assertEqual([(f["sample_id"], f["barcode"]) for f in request["files"]],
                             [("26R266-1", "barcode13"), ("26R266-2", "barcode14")])
            self.assertTrue(all(f["declared_complete"] for f in request["files"]))
            self.assertTrue(request["policy"]["write_locus_fastas"])
            self.assertNotIn(b"\r", output.read_bytes())
            before = output.read_bytes()
            self.assertNotEqual(subprocess.run(args, capture_output=True).returncode, 0)
            self.assertEqual(output.read_bytes(), before)
            args.remove("--declared-complete")
            self.assertNotEqual(subprocess.run(args, capture_output=True).returncode, 0)
            self.assertFalse((root / "run").exists())

    def test_all_tutorial_stages_have_registered_inner_agent_routes(self):
        operations = ["BuildRnaReadPanel", "ScreenRnaReadPanel", "PrepareRnaReadPanelDemo",
                      "InspectRnaReadPanel", "CancelRnaReadPanel", "ExportRnaReadPanelSummary"]
        routes = ["rna-reads panel-" + suffix for suffix in ("build", "screen", "demo", "show", "cancel", "export-summary")]
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())["commands"]
        protocol = (ROOT / "crates/gentle-protocol/src/lib.rs").read_text()
        mcp = (ROOT / "src/mcp_server.rs").read_text()
        bridge = (ROOT / "src/agent_bridge.rs").read_text()
        matrix = (ROOT / "docs/gui_cli_mcp_parity.md").read_text()
        for operation, route in zip(operations, routes):
            entry, = [row for row in glossary if row["path"] == route]
            self.assertEqual(entry["engine_operations"], [operation])
            self.assertIn('"' + operation + '"', protocol)
            self.assertIn('"' + route + '"', mcp)
            self.assertIn(route.split(" ", 1)[1], bridge)
            self.assertIn("| " + route + " | glossary-command |", matrix)
        guide = GUIDE.read_text()
        for phrase in ("Follow This Through The Inner Agent", "panel-show", "panel-cancel",
                       "26R266-1", "barcode13", "26R266-2", "barcode14",
                       "not a validated assay", "unresolved", "artifact_content_verified=false"):
            self.assertIn(phrase, guide)
        self.assertIn("rna_panel_tutorial_is_executable_through_shared_agent_commands",
                      (ROOT / "src/engine/analysis/rna_panel/tests.rs").read_text())
        shell = (ROOT / "src/engine_shell.rs").read_text()
        self.assertIn("request_template_is_example_not_resolved_input", shell)
        for schema in ("gentle.rna_read_panel_cancel_request.v1", "gentle.rna_read_panel_inspection_command.v1"):
            self.assertIn(schema, shell)

    def test_parity_table_totals_match_its_rows(self):
        text = (ROOT / "docs/gui_cli_mcp_parity.md").read_text()
        lines = text.split("## Open Gaps", 1)[0].splitlines()
        rows = [line.split("|")[1:-1] for line in lines
                if re.match(r"\| .* \| (glossary-command|engine-operation|mcp-tool) \|", line)]
        for column, adapter in enumerate(("GUI", "gentle_cli", "MCP", "JS", "Lua", "ClawBio"), 2):
            expected = [sum(row[column].strip() == state for row in rows)
                        for state in ("prominent", "shell-only", "gap")]
            match = re.search(r"^\| " + re.escape(adapter) + r" \| (\d+) \| (\d+) \| (\d+) \|$", text, re.M)
            self.assertEqual([int(v) for v in match.groups()], expected)

    def test_summary_lesson_preserves_ambiguity_and_has_an_agent_executable_request(self):
        for guide in (GUIDE.read_bytes(), GUIDE.read_bytes().replace(b"\n", b"\r\n")):
            text = guide.decode().replace("\r\n", "\n")
            requests = [json.loads(block) for block in re.findall(r"```json\n(.*?)```", text, re.S)]
            summary, = [r for r in requests if r.get("schema") == "gentle.rna_read_panel_summary_request.v1"]
            self.assertEqual(set(summary), {"schema", "report_paths", "output_dir"})
            self.assertEqual(len(summary["report_paths"]), 2)
            for phrase in ("single_candidate_reads", "co_nominated_reads", "member_ids_json",
                           "not original molecule copy numbers", "completeness: partial", "no numerical rows",
                           "panel-export-summary @$RUN/artifacts/summary.json", "panel-show $RUN/summary/report.json"):
                self.assertIn(phrase, text)
        implementation = (ROOT / "src/engine/analysis/rna_panel/tools.rs").read_text()
        self.assertIn('join("summary.json")', implementation)
        self.assertIn("synthetic_library={number}", implementation)
        tests = (ROOT / "src/engine/analysis/rna_panel/summary_tests.rs").read_text()
        self.assertIn("without_resolving_ambiguous_isoforms_or_reopening_inputs", tests)
        self.assertIn("fs::remove_file", tests)

    def test_repaired_lesson_explains_limits_and_all_new_shared_controls(self):
        text = GUIDE.read_text()
        for phrase in ("demo_1", "supported_junction_offsets: [40]", "whole-read",
                       "legacy_class_unresolved", "class_member_counts_complete=false",
                       "sample_coverage", "unassessed", "inputs_not_completed", "MD5SUMS",
                       "merged delivery FASTQ", "constituent chunks", "decoder read-ahead",
                       "file_checksum_verified", "pre-start", "record index",
                       "--verify-artifacts", "256 MB", "co_nomination_sets.tsv"):
            self.assertIn(phrase, text)
        bridge = (ROOT / "src/agent_bridge.rs").read_text()
        for phrase in ("--verify-artifacts", "legacy_class_unresolved", "co_nominated_reads",
                       "sample_coverage", "decoder read-ahead"):
            self.assertIn(phrase, bridge)
        # These are source-contract assertions, not execution of the Rust paths.
        control = (ROOT / "src/engine/analysis/rna_panel/control.rs").read_text()
        self.assertIn("panel_cancel_marker_mismatch", control)
        source = (ROOT / "src/engine/analysis/rna_panel.rs").read_text()
        polling = source.split("let mut gate = control::PollGate::default();", 1)[1]
        self.assertLess(polling.index("if !gate.due"), polling.index("control::marker_requested(request.cancel_path"))

    def test_control_errors_are_taught_as_failed_analysis_with_preserved_units(self):
        for raw in (GUIDE.read_bytes(), GUIDE.read_bytes().replace(b"\n", b"\r\n")):
            text = re.sub(r"\s+", " ", raw.decode())
            for phrase in ("status: control_error", "panel_control_marker_error", "marker_path",
                           "file_boundary", "streaming", "finalization", "source_runs[].status",
                           "not your authorized cancellation", "unfinished unit is discarded",
                           "fresh cancel path", "no automatic restart", "before any unit completes",
                           "only for post-start marker errors", "process exit", "competitors from the index",
                           "16 MB across", "all files in one screen request"):
                self.assertIn(phrase, text)
        protocol = (ROOT / "src/engine/protocol/rna_panel.rs").read_text()
        self.assertEqual(len(re.findall(r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s*pub control_error: Option<RnaReadPanelControlError>', protocol)), 2)
        bridge = (ROOT / "src/agent_bridge.rs").read_text()
        self.assertIn("FASTQ control-error rule", bridge)
        self.assertIn("Never silently overwrite the marker/report or auto-retry", bridge)
        feedback = (ROOT / "src/agent_feedback.rs").read_text()
        self.assertIn('"control_error" => Self::Failed', feedback)
        tests = (ROOT / "src/engine/analysis/rna_panel/tests.rs").read_text()
        for name in ("rna_panel_mid_run_marker_errors_preserve_verified_units_with_explicit_outcomes",
                     "rna_panel_control_repair_keeps_checksum_and_resource_failures_hard",
                     "rna_panel_control_error_binding_cannot_be_relabelled_as_success_or_legacy"):
            self.assertIn("fn " + name, tests)


if __name__ == "__main__":
    unittest.main()
