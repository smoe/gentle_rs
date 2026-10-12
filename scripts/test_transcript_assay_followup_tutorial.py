"""Offline companion for the hand-written transcript-assay follow-up tutorial.

Reuses committed synthetic PATZ1 fixtures and canonical 04.06 reports; origins
and request recreation are in docs/tutorial/inputs/README.md. Source checks do
not execute biology. Set GENTLE_TUTORIAL_BIN_DIR for real CLI/shared-shell
replay, using only temporary state and the internal backend. Never builds,
launches BLAST, executes an approved study, contacts an agent, or submits orders.
"""

import copy
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = ROOT / "docs/tutorial"
STEM = "04-07_transcript_assay_followup_gui_cli"
INPUTS = TUTORIAL / "inputs"
RT_REQUEST = INPUTS / "transcript_assay_followup_rt_pool.json"
STUDY_REQUEST = INPUTS / "transcript_assay_followup_study.json"
FIXTURE = ROOT / "test_files/fixtures/transcript_assay_panel/patz1"
WORKFLOW = ROOT / "docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json"
ARTIFACTS = TUTORIAL / "generated/artifacts/patz1_transcript_assay_panels_cli/artifacts"
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(data):
    return "sha256:" + hashlib.sha256(data).hexdigest()


def comparable_rt_report(report):
    """Remove only top-level execution identity, retaining every scientific field."""
    return {key: value for key, value in report.items()
            if key not in {"generated_at_unix_ms", "op_id", "run_id"}}


class PanelAssertions:
    def assert_panel(self, panel):
        self.assertEqual(panel["schema"], "gentle.transcript_assay_panel.v2")
        self.assertEqual(panel["completion_status"], "complete")
        self.assertEqual(panel["strand"], "-")
        self.assertEqual(panel["transcript_count"], 3)
        self.assertEqual(panel["equivalence_group_count"], 3)
        self.assertGreater(panel["selected_assay_count"], 0)
        self.assertEqual(panel["uncovered_equivalence_group_ids"], [])
        self.assertEqual(panel["provenance"]["primer_backend"], "internal")
        self.assertEqual({row["genomic_confirmation_status"]
                          for row in panel["specificity_followups"]}, {"not_run"})
        self.assertTrue(panel["order_ready_primers"])
        self.assertTrue(all(row["role"] in {"forward", "reverse"}
                            for row in panel["order_ready_primers"]))

    def assert_juc_missing_threshold(self, panel):
        self.assertEqual(panel["selected_assay_count"], 3)
        junction, = panel["junction_evaluations"]
        self.assertEqual(junction["priority"], "required")
        self.assertEqual(junction["status"], "selected_spanning_assay")
        self.assertEqual(junction["resolved_transcript_ids"], ["PATZ1-202"])
        self.assertEqual(junction["local_positions_0based"], [40])
        juc_rows = []
        for assay in panel["selected_assays"]:
            summary = assay["primer_pair_summary"]
            self.assertEqual(summary["selection_audit_status"], "computed")
            self.assertEqual(summary["selection_audit_method"],
                             "binary_detection_leave_one_out_v1")
            self.assertEqual(summary["selection_operation_sha256"], panel["operation_sha256"])
            self.assertFalse(summary["retained_because_of_differential_junction_evidence"])
            for evidence in summary["selection_evidence"]:
                if evidence["evidence_kind"] == "juc":
                    juc_rows.append(evidence)
                    self.assertEqual(summary["selection_evidence_observation_count"], 1)
                    self.assertEqual(summary["selection_evidence_projection_count"], 1)
        row, = juc_rows
        self.assertEqual(row["evidence_id"], "JUC2200054820.hg.1")
        self.assertEqual(row["measured_value"], -1.2)
        self.assertEqual(row["contrast"], "synthetic_case_vs_control")
        self.assertEqual(row["differential_eligibility"], "not_assessed")
        self.assertEqual(row["disposition"], "incomplete_missing_threshold")


class FollowupSourceTests(PanelAssertions, unittest.TestCase):
    def test_source_is_discoverable_without_claiming_generated_or_human_review(self):
        source = read_json(TUTORIAL / "sources" / f"{STEM}.json")
        self.assertEqual(source["schema"], "gentle.tutorial_source.v4")
        self.assertEqual(source["id"], "transcript_assay_followup_gui_cli")
        self.assertNotIn("generated_chapter", source)
        catalog = source["catalog"]
        self.assertEqual(catalog["source"], "hand_written_markdown")
        self.assertEqual(catalog["status"], "manual/hybrid")
        self.assertEqual((catalog["group"], catalog["group_position"]), ("04", 7))
        for path in (TUTORIAL / "sources").glob("*.json"):
            other = read_json(path)
            if other.get("id") == source["id"] or "catalog" not in other:
                continue
            self.assertNotEqual(other["catalog"]["order"], catalog["order"], str(path))
            self.assertNotEqual((other["catalog"].get("group"),
                                 other["catalog"].get("group_position")), ("04", 7))
        # Guide retrieval indexes title, notes (as summary) and audiences, not Markdown.
        indexed = (source["title"] + " " + catalog["notes"]).lower()
        for term in ["terminal-exon", "rt primer", "differential", "junction",
                     "specificity", "gene-isoform", "study", "offline"]:
            self.assertIn(term, indexed)
        guide = ROOT / catalog["path"]
        for target in re.findall(r"\]\(([^)]+)\)", guide.read_text(encoding="utf-8")):
            if "://" not in target:
                self.assertTrue((guide.parent / target.split("#")[0]).exists(), target)

    def test_requests_bind_reused_fixtures_without_fabricating_threshold_provenance(self):
        study = read_json(STUDY_REQUEST)
        self.assertEqual(sha256((ROOT / study["isoform_evidence_path"]).read_bytes()),
                         study["expected_isoform_evidence_sha256"])
        evidence_ref, = study["junction_evidence"]
        path = ROOT / evidence_ref["path"]
        self.assertEqual(sha256(path.read_bytes()), evidence_ref["expected_sha256"])
        evidence = read_json(path)
        for field in ["min_abs_logfc", "threshold_source", "policy_sha256"]:
            self.assertIsNone(evidence.get(field))
        self.assertEqual(study["policy"]["min_abs_regional_effect"], 0.5)
        self.assertEqual(study["profile_override"]["selected_profile"],
                         "targeted_junction_validation")
        self.assertTrue(study["profile_override"]["reason"])
        ledger = read_json(ROOT / study["isoform_evidence_path"])["data"]
        self.assertEqual(ledger["schema"], "gentle.gene_isoform_evidence.v2")
        self.assertEqual(ledger["splicing"]["target_feature_id"], 2)
        readme = (INPUTS / "README.md").read_text(encoding="utf-8")
        for path in [RT_REQUEST, STUDY_REQUEST]:
            self.assertIn(path.name, readme)

    def test_rt_request_matches_documented_minus_strand_geometry(self):
        request = read_json(RT_REQUEST)
        self.assertNotIn("genomic_specificity", request)
        self.assertEqual(request["fixed_adapter_5prime"], "ACTTGCCTGTCGCTCTATCTTC")
        self.assertEqual(request["variable_length_bp"], 22)
        self.assertEqual(request["terminal_exon_search_window_bp"], 40)
        self.assertEqual(request["max_candidates_per_target"], 5)
        self.assertEqual([row["transcript_id"] for row in request["targets"]],
                         ["PATZ1-201", "PATZ1-202"])
        self.assertTrue(all(row["source_feature_id"] == 1 for row in request["targets"]))
        genbank = (FIXTURE / "patz1_assay_minus_strand.gb").read_text(encoding="utf-8")
        bases = re.sub(r"[^acgt]", "", genbank.split("ORIGIN\n")[1].split("//")[0])
        self.assertEqual(len(bases), 240)
        for location in ["complement(join(21..60,101..140,181..220))",
                         "complement(join(21..60,181..220))",
                         "complement(join(21..60,121..140,181..220))"]:
            self.assertIn(location, genbank)
        workflow = read_json(WORKFLOW)["workflow"]["ops"]
        self.assertEqual(workflow[0]["LoadFile"]["as_id"], request["targets"][0]["seq_id"])

    def test_existing_reports_support_the_teaching_observations(self):
        for name in ["patz1_endpoint_end_matrix", "patz1_sybr_juc_panel",
                     "patz1_routine_common_region_screen"]:
            with self.subTest(report=name):
                self.assert_panel(read_json(ARTIFACTS / f"{name}.report.json"))
        self.assert_juc_missing_threshold(read_json(ARTIFACTS / "patz1_sybr_juc_panel.report.json"))

    def test_parity_comparison_keeps_geometry_and_policy_changes_visible(self):
        report = {"op_id": "a", "run_id": "a", "generated_at_unix_ms": 1,
                  "targets": [{"source_start_0based": 20}], "tm_policy": "descriptive"}
        other = dict(report, op_id="b", run_id="b", generated_at_unix_ms=2)
        self.assertEqual(comparable_rt_report(report), comparable_rt_report(other))
        changed = copy.deepcopy(other)
        changed["targets"][0]["source_start_0based"] = 21
        self.assertNotEqual(comparable_rt_report(report), comparable_rt_report(changed))
        changed = dict(other, tm_policy="ranked")
        self.assertNotEqual(comparable_rt_report(report), comparable_rt_report(changed))


@unittest.skipUnless(BIN_DIR, "Set GENTLE_TUTORIAL_BIN_DIR for real offline CLI/Shell replay")
class FollowupReplayTests(PanelAssertions, unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gentle-transcript-followup-")
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)
        self.state = self.out / "followup.gentle.json"
        suffix = ".exe" if os.name == "nt" else ""
        self.cli = Path(BIN_DIR).resolve() / f"gentle_cli{suffix}"
        self.invoke("op", json.dumps(read_json(WORKFLOW)["workflow"]["ops"][0]))
        self.invoke("op", json.dumps({"SetParameter": {
            "name": "primer_design_backend", "value": "internal"}}))

    def invoke(self, *args, expect_success=True):
        result = subprocess.run([str(self.cli), "--state", str(self.state), *map(str, args)],
                                cwd=ROOT, capture_output=True, text=True, timeout=180)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            return json.loads(result.stdout)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def shell(self, *args, **kwargs):
        return self.invoke("shell", shlex.join(list(map(str, args))), **kwargs)

    def test_rt_pool_direct_shell_export_and_invalid_request(self):
        direct = self.invoke("primers", "design-terminal-exon-rt-pool", f"@{RT_REQUEST}")["report"]
        shared = self.shell("primers", "design-terminal-exon-rt-pool", f"@{RT_REQUEST}")["report"]
        self.assertEqual(comparable_rt_report(direct), comparable_rt_report(shared))
        self.assertEqual(direct["schema"], "gentle.terminal_exon_rt_primer_pool.v1")
        self.assertIsNone(direct.get("genomic_specificity"))
        self.assertIn("excluded from ranking", direct["tm_policy"])
        self.assertEqual(len(direct["selected_pool_interactions"]), 1)
        self.assertEqual(len(direct["targets"]), 2)
        for target, transcript_id, start in zip(direct["targets"], ["PATZ1-201", "PATZ1-202"], [80, 40]):
            self.assertEqual(target["resolved_transcript_id"], transcript_id)
            self.assertEqual(target["strand"], "-")
            self.assertEqual((target["terminal_exon_source_start_0based"],
                              target["terminal_exon_source_end_0based_exclusive"]), (20, 60))
            self.assertEqual((target["terminal_exon_transcript_start_0based"],
                              target["terminal_exon_transcript_end_0based_exclusive"]), (start, start + 40))
            self.assertEqual(target["evaluated_candidate_count"], 19)
            self.assertEqual(target["ambiguous_candidate_count"], 0)
            self.assertEqual(len(target["candidates"]), 5)
            selected, = [row for row in target["candidates"] if row["selected"]]
            self.assertEqual(len(selected["variable_primer_5_to_3"]), 22)
            self.assertEqual(selected["variable_primer_5_to_3"],
                             selected["target_segment_5_to_3"].translate(str.maketrans("ACGT", "TGCA"))[::-1])
            self.assertEqual(selected["full_oligo_5_to_3"],
                             direct["fixed_adapter_5prime"] + selected["variable_primer_5_to_3"])
        export = self.out / "rt_pool.json"
        self.shell("primers", "export-report", direct["report_id"], export)
        self.assertEqual(comparable_rt_report(direct), comparable_rt_report(read_json(export)))
        before = self.state.read_bytes()
        invalid = dict(read_json(RT_REQUEST), invent_specificity_pass=True)
        refused = self.invoke("primers", "design-terminal-exon-rt-pool", json.dumps(invalid), expect_success=False)
        self.assertIn("unknown field", refused.stdout + refused.stderr)
        self.assertEqual(self.state.read_bytes(), before)

    def test_existing_panel_operations_and_persisted_shell_export(self):
        operations = read_json(WORKFLOW)["workflow"]["ops"]
        for index in [2, 3, 4]:
            operation = copy.deepcopy(operations[index])
            request = operation["DesignTranscriptAssayPanel"]
            report_id = request["report_id"]
            output = self.out / f"{report_id}.json"
            request["path"] = str(output)
            operation_path = self.out / f"{report_id}.operation.json"
            operation_path.write_text(json.dumps(operation), encoding="utf-8")
            self.invoke("primers", "design-transcript-assay-panel", f"@{operation_path}", "--backend", "internal")
            report = read_json(output)
            self.assert_panel(report)
            if index == 3:
                self.assert_juc_missing_threshold(report)
            export = self.out / f"{report_id}.export.json"
            before = self.state.read_bytes()
            self.shell("primers", "show-transcript-assay-panel", report_id)
            self.shell("primers", "export-transcript-assay-panel", report_id, export)
            self.assertEqual(read_json(export), report)
            self.assertEqual(self.state.read_bytes(), before)

    def test_study_normalization_planning_and_changed_workflow_refusal(self):
        before = self.state.read_bytes()
        normalized_path = self.out / "study.normalized.json"
        normalized = self.invoke("primers", "plan-gene-isoform-study", f"@{STUDY_REQUEST}",
                                 "--normalize-only", "--normalized-request", normalized_path)
        self.assertEqual(normalized, read_json(normalized_path))
        self.assertEqual(normalized["schema"], "gentle.gene_isoform_assay_study_plan_request.v1")
        self.assertTrue(normalized["policy_sha256"].startswith("sha256:"))
        self.assertEqual(normalized["fallback_submission"]["mode"], "never")
        plan_path, workflow_path = self.out / "study.plan.json", self.out / "study.workflow.json"
        plan = self.shell("primers", "plan-gene-isoform-study", f"@{normalized_path}",
                          "--path", plan_path, "--workflow", workflow_path)["report"]
        self.assertEqual(self.state.read_bytes(), before)
        self.assertEqual(plan, read_json(plan_path))
        self.assertEqual(plan["schema"], "gentle.gene_isoform_assay_study_plan.v1")
        self.assertEqual(plan["selected_profile"], "targeted_junction_validation")
        self.assertEqual(plan["approved_workflow_sha256"], sha256(workflow_path.read_bytes()))
        planned = [row["operation"] for row in plan["planned_operations"]]
        self.assertEqual(read_json(workflow_path)["ops"], planned)
        self.assertEqual([row["DesignTranscriptAssayPanel"]["objective"] for row in planned],
                         ["pan_transcript", "minimal_discrimination_panel"])
        self.assertTrue(all(row["DesignTranscriptAssayPanel"]["coverage_policy"] == "require_all"
                            for row in planned))
        self.assertEqual(planned[1]["DesignTranscriptAssayPanel"]["junction_evidence_priority"], "preferred")
        self.assertTrue(any(row["rule_id"] == "junction_interpretation_incomplete_provenance"
                            and row["triggered"] for row in plan["decision_factors"]))
        tampered = self.out / "study.changed.workflow.json"
        tampered.write_bytes(workflow_path.read_bytes() + b"\n")
        refused = self.invoke("primers", "execute-gene-isoform-study-workflow", plan_path, tampered,
                              expect_success=False)
        self.assertIn("Approved workflow digest mismatch", refused.stdout + refused.stderr)
        self.assertEqual(self.state.read_bytes(), before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
