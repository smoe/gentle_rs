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
    def test_historical_source_and_budget_review_repairs_have_scoped_contracts(self):
        panel = (ROOT / "src/engine/analysis/rna_panel.rs").read_text(encoding="utf-8")
        legacy = (ROOT / "src/engine/analysis/rna_reads.rs").read_text(encoding="utf-8")
        types = (ROOT / "src/engine/protocol/rna_panel.rs").read_text(encoding="utf-8")
        tests = (ROOT / "src/engine/analysis/rna_panel/tests.rs").read_text(encoding="utf-8")
        for text in (panel, panel.replace("\n", "\r\n")):
            self.assertIn("same_sequence", text)
            self.assertTrue(re.search(r"a\s*\.source_sequence_sha256", text), "DNA comparison missing")
            self.assertIn("removes unselected competitors", text)
            self.assertIn("MAX_RUN_EQUIVALENCE_KEY_BYTES", text)
            self.assertIn("used.checked_add(additional)", re.sub(r"\s+", " ", text))
            screen = text.split("let mut equivalence_bytes = 0usize;", 1)[1]
            self.assertIn("for (ordinal, source) in request.files.iter().enumerate()", screen)
            self.assertIn("admit_run_key_bytes(equivalence_bytes, equivalence.len())", re.sub(r"\s+", " ", screen))
        self.assertRegex(types, r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s*pub source_sequence_sha256: Option<String>')
        for text in (legacy, legacy.replace("\n", "\r\n")):
            collect = text.split("pub(super) fn collect_sparse_target_transcript_lanes", 1)[1].split(
                "pub(super) fn build_isoform_support_accumulators", 1)[0]
            self.assertLess(collect.index("if restrict_to_target_strand && is_reverse"),
                            collect.index("matched_requested.insert"))
            self.assertLess(collect.index("matched_requested.insert"), collect.index("already_present_feature_ids.contains"))
        for name in ("rna_panel_cross_source_identical_dna_without_anchors_is_flagged_before_streaming",
                     "rna_panel_legacy_locus_sequence_hash_is_absent_not_fabricated_and_stays_digest_bound",
                     "rna_panel_shared_isoforms_can_exceed_the_default_anchor_budget_without_publishing",
                     "rna_panel_whole_gzip_member_loss_is_caught_by_the_expected_transfer_digest"):
            self.assertIn("fn " + name, tests)
        for raw in (GUIDE.read_bytes(), GUIDE.read_bytes().replace(b"\n", b"\r\n")):
            text = re.sub(r"\s+", " ", raw.decode())
            for phrase in ("source_sequence_sha256", "Peak anchors on successful reads still need profiling",
                           "Zero-length reads", "alternative exon paths", "not a drop-in replacement"):
                self.assertIn(phrase, text)

    def test_glossary_rna_flags_keep_typed_binding_and_required_draft_mode(self):
        source = (ROOT / "src/engine_shell/tests.rs").read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            self.assertIn('"BINDING_JSON_OR_@FILE" => {', text)
            self.assertIn('"schema": "gentle.rna_read_panel_result_binding.v1"', text)
            self.assertIn('path == "rna-reads panel-screen" && flag == "--save"', text)
            self.assertIn('tokens.insert(2, "--request".into())', text)
            self.assertIn("fn glossary_rna_flag_smokes_retain_typed_binding_and_draft_save_mode()", text)
            self.assertIn("assert!(parse_shell_line(&without_request).is_err())", text)

    def test_rna_parity_fixtures_preserve_exact_path_and_request_hash_bindings(self):
        mcp = (ROOT / "src/mcp_server.rs").read_text(encoding="utf-8")
        summaries = (ROOT / "src/engine/analysis/rna_panel/summary_tests.rs").read_text(encoding="utf-8")
        for text in (mcp, mcp.replace("\n", "\r\n")):
            self.assertIn('"path":binding.report_path, "expected_result":binding', text)
            self.assertIn('contains("panel_result_content_binding")', text)
            self.assertIn("let mut expected_mcp_summary = expected_summary.clone()", text)
            self.assertIn("canonical_hash(serde_json::to_value(&expected_mcp_summary).unwrap())", text)
        for text in (summaries, summaries.replace("\n", "\r\n")):
            self.assertIn("fn expected_summary_for_request(", text)
            self.assertIn("Some(sha256_prefixed_bytes(&encode(request).unwrap()))", text)
            self.assertIn("expected.content_sha256 = sha256_prefixed_bytes(&encode(&expected).unwrap())", text)
            self.assertIn("rna_read_panel_summary_result_binding(&reversed, &repeated)", text)
        binding = (ROOT / "src/engine/analysis/rna_panel/result_binding.rs").read_text(encoding="utf-8")
        self.assertIn("binding.report_path != path", binding)

    def test_native_panel_module_paths_resolve_to_tracked_children(self):
        modules = (
            ("src/app.rs", "rna_panel_ui", "app/rna_panel_ui.rs"),
            ("src/app/rna_panel_ui.rs", "processing", "rna_panel_ui/processing.rs"),
            ("src/app/rna_panel_ui/processing.rs", "evidence", "processing/evidence.rs"),
        )
        for parent, child, relative in modules:
            source_path = ROOT / parent
            source = source_path.read_text(encoding="utf-8")
            for text in (source, source.replace("\n", "\r\n")):
                with self.subTest(parent=parent, child=child, crlf="\r\n" in text):
                    match = re.search(r'#\[path = "([^"]+)"\]\s*mod ' + child + r';', text)
                    self.assertIsNotNone(match, f"{parent} needs an explicit path for {child}")
                    self.assertEqual(match.group(1), relative)
                    self.assertTrue((source_path.parent / match.group(1)).is_file())

    def test_command_admission_baseline_is_not_shadowed_by_a_journal_count(self):
        source = (ROOT / "src/command_execution.rs").read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            worker = text.split("let mut detached = if let Some", 1)[1].split(
                "let mut bytes = bounded_output", 1)[0]
            self.assertIn("let journal_before_execution = detached.engine().journal_len();", worker)
            self.assertIn("saturating_sub(journal_before_execution)", worker)
            self.assertIn("project_commit && (baseline.is_none()", worker)
            self.assertIn("Some(detached.engine().instance_id()) != owner", worker)
            self.assertNotRegex(worker, r"let\s+baseline\s*=")

    def test_draft_file_loading_is_lazy_metadata_observation_in_every_checkout(self):
        sources = {path: (ROOT / path).read_text(encoding="utf-8") for path in (
            "src/engine_shell/command_parsers.rs", "src/command_execution/rna_artifacts.rs",
            "src/engine/analysis/rna_panel/tools.rs", "src/app/rna_panel_ui/processing.rs",
            "src/engine/analysis/rna_panel/result_binding_tests.rs")}
        for crlf in (False, True):
            with self.subTest(crlf=crlf):
                parser, policy, tools, native, tests = (
                    text.replace("\n", "\r\n") if crlf else text for text in sources.values())
                branch = parser.split('"panel-screen" if tokens.get(2)', 1)[1].split(
                    '"panel-build" if', 1)[0]
                self.assertLess(branch.index("ShellCommand::RnaReadPanelLoadDraft"),
                                branch.index("parse_required_json_payload"))
                self.assertRegex(policy, r"RnaReadPanelLoadDraft\s*\{\s*\.\.\s*\}\s*=> Policy::Observe\s*\{\s*expensive: false")
                loader = tools.split("pub fn load_rna_read_panel_screen_draft", 1)[1].split(
                    "pub fn rna_read_panel_screen_draft", 1)[0]
                self.assertIn("file.take(DRAFT_JSON_BYTES + 1)", loader)
                self.assertIn("panel_draft_json_budget", loader)
                self.assertIn('text.starts_with("#!")', loader)
                self.assertIn("self.rna_read_panel_screen_draft(&request, None)", loader)
                native_load = native.split("fn start_rna_draft_load", 1)[1].split(
                    "pub(super) fn render_rna_rehearsal_controls", 1)[0]
                self.assertIn("self.submit_rna_processing(", native_load)
                self.assertIn("ShellCommand::RnaReadPanelLoadDraft", native_load)
                self.assertNotIn(".submit(", native_load)
                self.assertIn("rna_screen_draft_load_uses_observation_with_busy_execution_and_locked_project", native)
                self.assertIn("rna_panel_screen_draft_file_is_lazy_bounded_and_does_not_open_referenced_inputs", tests)

    def test_native_draft_retry_clears_only_a_successful_operations_old_error(self):
        source = (ROOT / "src/app/rna_panel_ui/processing.rs").read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            draft = text.split("Task::Draft { load } => match", 1)[1].split("Task::Inspect(", 1)[0]
            self.assertIn("Ok(draft) => {", draft)
            self.assertIn("state.error = None", draft)
            self.assertLess(draft.index("state.error = None"), draft.index("if load {"))
            self.assertIn("Err(error) => state.error = Some(error)", draft)
            self.assertIn("fn rna_screen_draft_load_retry_preserves_review_on_failure_and_clears_stale_error()", text)
            self.assertIn("state.reviewed_request().unwrap()", text)
            self.assertIn('"successful retry must clear the old failure"', text)

    def test_model_stop_command_cancel_and_project_change_have_distinct_contracts(self):
        source = (ROOT / "src/app/routine_and_agent_assistant_ui.rs").read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            guard = text.split("pub(super) fn cancel_agent_commands_for_project_change", 1)[1].split(
                "pub(super) fn finish_agent_shell_run", 1)[0]
            self.assertIn("receipt.requires_project_commit", guard)
            self.assertIn("self.agent_command_service.cancel(task.job_id)", guard)
            self.assertIn("self.request_agent_command_cancel(task.job_id)", text)
            for name in ("agent_rna_command_cancel_retains_partial_report_with_or_without_marker",
                         "agent_rna_screen_survives_model_stop_and_project_replacement"):
                self.assertIn("fn " + name, text)
            replacement = text.split("fn agent_rna_screen_survives_model_stop_and_project_replacement()", 1)[1].split(
                "fn run_and_prompt_submit_exact_command", 1)[0]
            self.assertIn("for marker in [true, false]", replacement)
            self.assertIn("held_rna_screen(&mut app, td.path(), marker)", replacement)
        for raw in (GUIDE.read_bytes(), GUIDE.read_bytes().replace(b"\n", b"\r\n")):
            text = re.sub(r"\s+", " ", raw.decode())
            for phrase in ("Which Stop Or Cancel Do I Use?", "independent observation lane",
                           "128 MB", "16 MiB", "same-process", "without a marker path",
                           "must not silently fall back", "after publication cannot relabel it"):
                self.assertIn(phrase, text)

    def test_full_native_rehearsal_and_prominence_have_actual_shared_handlers(self):
        source = (ROOT / "src/app/rna_panel_ui/processing.rs").read_text()
        evidence = (ROOT / "src/app/rna_panel_ui/processing/evidence.rs").read_text()
        native = (ROOT / "src/app/window_registry.rs").read_text()
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())["commands"]
        operations = {"demo": "PrepareRnaReadPanelDemo", "screen": "ScreenRnaReadPanel",
                      "show": "InspectRnaReadPanel", "cancel": "CancelRnaReadPanel",
                      "export-summary": "ExportRnaReadPanelSummary"}
        for suffix, operation in operations.items():
            route = "rna-reads panel-" + suffix
            entry, = [row for row in glossary if row["path"] == route]
            self.assertIn("gui-menu", entry["interfaces"])
            self.assertIn('glossary_path: "' + route + '"', native)
            self.assertIn("Operation::" + operation, source + evidence)
        self.assertIn("rna_native_screen_monitor_bound_reopen_and_summary_work_without_the_live_project_lock", evidence)
        guide = re.sub(r"\s+", " ", GUIDE.read_text())
        for phrase in ("GUI-Only Rehearsal", "Typed Checkpoints And Capture Targets",
                       "DEMO_1", "DEMO_2", "harness_gap", "no capture", "verified_tables_only"):
            self.assertIn(phrase, guide)

    def test_native_inspector_and_export_use_shared_records_and_preserve_ambiguity(self):
        source = (ROOT / "src/app/rna_panel_ui/processing/evidence.rs").read_text()
        for phrase in ("Operation::InspectRnaReadPanel", "Operation::ExportRnaReadPanelSummary",
                       "expected_result:", "expected_content", "co_nomination_sets", "member_ids",
                       "legacy_class_unresolved", "artifact_content_verified", "sample_coverage",
                       "control_error", "export_request", "show_rows", "RnaReadPanelSummary"):
            self.assertIn(phrase, source)
        for forbidden in ("File::open", "visit_fastq", "GzDecoder"):
            self.assertNotIn(forbidden, source)
        guide = re.sub(r"\s+", " ", GUIDE.read_text())
        for phrase in ("Open bound result", "Inspect nomination prefix", "Add screen to summary",
                       "Export reviewed summary (writes files)", "Verify summary tables"):
            self.assertIn(phrase, guide)

    def test_native_screen_uses_shared_service_and_explicit_review_without_fastq_io(self):
        source = (ROOT / "src/app/rna_panel_ui/processing.rs").read_text()
        for phrase in ("agent_command_service.submit_parsed", "Operation::ScreenRnaReadPanel",
                       "Operation::CancelRnaReadPanel", "Operation::InspectRnaReadPanel",
                       "expected_panel_sha256", "reviewed_request", "declared_complete: false",
                       "decoder read-ahead", "rna_screen_window_close_does_not_cancel"):
            self.assertIn(phrase, source)
        for forbidden in ("visit_fastq", "GzDecoder", "read_to_end", "File::open"):
            self.assertNotIn(forbidden, source)
        parser = (ROOT / "src/engine_shell/command_parsers.rs").read_text()
        self.assertIn('"--request"', parser)
        self.assertIn("ShellCommand::RnaReadPanelScreenDraft", parser)
        guide = GUIDE.read_text()
        for phrase in ("Native Setup And Monitoring", "Run reviewed screen (writes files)",
                       "File copy finished", "Copy report binding", "--save NEW_FILE"):
            self.assertIn(phrase, guide)

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
            self.assertIn("let mut demo_engine = GentleEngine::from_state(state)", text)
            self.assertRegex(text, r"demo_engine\s*\.state\(\)\s*\.save_to_path\(")
            self.assertIn("upsert_gene_set_resolution_artifact", text)
        mcp = (ROOT / "src/mcp_server.rs").read_text(encoding="utf-8")
        self.assertIn("the synthetic starter must already contain reconciled lineage and containers", mcp)

    def test_path_loaded_panel_children_resolve_to_their_tracked_sources(self):
        source_path = ROOT / "src/engine/analysis/rna_panel.rs"
        source = source_path.read_text(encoding="utf-8")
        for text in (source, source.replace("\n", "\r\n")):
            for child in ("input", "control", "verification", "tools", "tests", "summary", "summary_tests", "project", "result_binding", "result_binding_tests"):
                match = re.search(r'#\[path = "([^"]+)"\]\s*mod ' + child + r';', text)
                self.assertIsNotNone(match, f"path-loaded panel needs an explicit {child} path")
                self.assertEqual(match.group(1), f"rna_panel/{child}.rs")
                self.assertTrue((source_path.parent / match.group(1)).is_file())

    def test_artifact_only_policy_and_request_bindings_are_not_inferred_from_flags(self):
        policy = (ROOT / "src/command_execution/rna_artifacts.rs").read_text(encoding="utf-8")
        for operation in ("ScreenRnaReadPanel", "ExportRnaReadPanelSummary", "InspectRnaReadPanel", "CancelRnaReadPanel"):
            self.assertIn("Operation::" + operation, policy)
        self.assertIn("_ => Policy::Project", policy)
        self.assertNotIn("state_changed", policy)
        types = (ROOT / "src/engine/protocol/rna_panel.rs").read_text(encoding="utf-8")
        for field in ("expected_panel_sha256", "submitted_request_sha256"):
            self.assertRegex(types, r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s*pub ' + field)
        tests = (ROOT / "src/command_execution/rna_tests.rs").read_text(encoding="utf-8")
        for name in ("rna_panel_service_late_cancel_and_live_edit_keep_published_completion",
                     "rna_panel_service_retention_does_not_stale_another_detached_job",
                     "rna_panel_control_and_metadata_do_not_compete_with_execution_slots_or_project_lock",
                     "rna_panel_registered_build_keeps_project_commit_and_mutation_feedback"):
            self.assertIn(name, tests)

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
        for field in ("interfaces", "starting_state", "estimated_minutes"):
            self.assertEqual(entry[field], source["catalog"][field], field)
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
