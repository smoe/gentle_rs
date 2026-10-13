"""Offline source checks for tutorial 04.09; no binaries or biological acceptance.

The Markdown JSON examples are hand-crafted learning templates, not gene data.
Tests verify catalog projection, discoverability and teaching contract only.
Rust engine/MCP execution and real-reference acceptance stay with Glen/CI.
"""

import json
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
STEM = "04-09_reference_bound_primer_specificity"


class PrimerSpecificityTutorialTests(unittest.TestCase):
    def setUp(self):
        self.source = json.loads((ROOT / f"docs/tutorial/sources/{STEM}.json").read_text())
        self.page = (ROOT / self.source["catalog"]["path"]).read_text()

    def test_result_literals_initialize_the_new_summary_fields(self):
        # Source-only preflight; the CI compiler remains the definitive check.
        count = 0
        for path in (ROOT / "src").rglob("*.rs"):
            text = path.read_text()
            initializers = re.findall(
                r"primer_specificity_multi_handoff:\s*None,\s*"
                r"primer_specificity_multi_summary:\s*None,\s*"
                r"primer_specificity_multi_summaries:\s*vec!\[\],",
                text,
            )
            handoffs = re.findall(r"primer_specificity_multi_handoff:\s*None,", text)
            self.assertEqual(len(initializers), len(handoffs), str(path))
            count += len(initializers)
        self.assertGreaterEqual(count, 10)

    def test_catalog_and_review_match_source_without_claiming_replay(self):
        sources = [json.loads(p.read_text()) for p in (ROOT / "docs/tutorial/sources").glob("*.json")]
        orders = [s["catalog"]["order"] for s in sources if "catalog" in s]
        self.assertEqual(len(orders), len(set(orders)))
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_text())
        row, = [r for r in catalog["entries"] if r["id"] == self.source["id"]]
        for key, value in self.source["catalog"].items():
            if key != "order":
                self.assertEqual(row[key], value, key)
        self.assertEqual(row["decimal_id"], "04.09")
        self.assertEqual(row["review_status"], "codex_reviewed")
        reviews = json.loads((ROOT / "docs/tutorial/review_manifest.json").read_text())
        review, = [r for r in reviews["entries"] if r["tutorial_id"] == self.source["id"]]
        self.assertEqual(row["codex_reviewed_at"], review["codex_reviewed_at"])
        self.assertIsNone(review["human_reviewed_at"])
        manifest = json.loads((ROOT / "docs/tutorial/manifest.json").read_text())
        self.assertNotIn(self.source["id"], [r["id"] for r in manifest["chapters"]])

    def test_learning_request_and_mapping_have_explicit_different_spaces(self):
        request, mapping = [json.loads(s) for s in re.findall(r"```json\n(.*?)\n```", self.page, re.S)]
        self.assertEqual(request["schema"], "gentle.primer_pair_multi_reference_request.v1")
        self.assertEqual(request["pair"]["kind"], "saved_pair")
        self.assertEqual(request["pair"]["pair_rank"], 1)
        self.assertNotIn("specificity_target_genome_id", request["policy"])
        self.assertEqual([r["required"] for r in request["references"]], [True, True, False])
        self.assertEqual({r["expected_index_kind"] for r in request["references"]}, {"genomic_dna", "transcriptome_cdna"})
        self.assertEqual(mapping["model"], "transcript_set")
        self.assertEqual(mapping["expected_products"][0]["target_space"], "transcriptome_cdna")
        self.assertIn("SYNTHETIC", mapping["expected_products"][0]["subject_id"])

    def test_commands_and_motivations_remain_discoverable(self):
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())
        for suffix in ["handoff", "import", "show", "list"]:
            path = f"primers specificity-multi-{suffix}"
            self.assertIn(path, self.page)
            row, = [r for r in glossary["commands"] if r["path"] == path]
            self.assertIn("mcp", row["interfaces"])
            self.assertEqual(len(row["engine_operations"]), 1)
        for topic in ["Why This Matters", "assembly", "caller-provided", "annealing", "required", "optional", "cancelled", "not_requested", "not_required", "selection receipt", "Legacy", "readiness", "authenticated", "confirm: true"]:
            self.assertIn(topic, self.page)
        for target in re.findall(r"\]\(([^)]+)\)", self.page):
            if "://" not in target:
                self.assertTrue((ROOT / "docs/tutorial" / target.split("#")[0]).exists(), target)

    def test_panel_matrix_teaches_scope_and_read_only_bindings(self):
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())
        row, = [r for r in glossary["commands"] if r["path"] == "primers transcript-assay-specificity-matrix"]
        self.assertEqual(row["engine_operations"], ["InspectTranscriptAssayPanelSpecificityMatrix"])
        self.assertIn("mcp", row["interfaces"])
        for ending in ("\n", "\r\n"):
            page = self.page.replace("\n", ending)
            for term in [row["path"], "transcript_assay_specificity_matrix", "is_superseded",
                         "matches_current_panel_design", "retained_binding_verified",
                         "unassessed_genomic_assay_ids", "unlinked_pair_summaries",
                         "specificity and", "isoform discrimination", "does not combine"]:
                self.assertIn(term, page)
        implementation = (ROOT / "src/engine/analysis/primer_specificity_matrix.rs").read_text()
        for forbidden in ["record_primer_specificity_reference_selection(",
                          "write_primer_design_store(", "inspect_blast_database(",
                          "selected_primer_specificity_report("]:
            self.assertNotIn(forbidden, implementation)
        for path in (ROOT / "src").rglob("*.rs"):
            text = path.read_text()
            for match in re.finditer(r"primer_specificity_multi_summaries:\s*vec!\[\],", text):
                self.assertRegex(text[match.end():], r"^\s*transcript_assay_specificity_matrix:\s*None,", str(path))

    def test_slice_b_is_deferred_and_slice_a_uses_the_printed_command(self):
        plan = (ROOT / "docs/primer_design_specialist_plan.md").read_text()
        for term in ["## Conditional Slice B: Bound Panel Fan-Out", "not implementation authorization",
                     "Motivation And Entry Conditions", "atomic cross-reference import",
                     "128 MiB total across the entire parent", "unknown sizes are not zero",
                     "record_primer_specificity_reference_selection", "Readiness attachment needs a later policy decision"]:
            self.assertIn(term, plan)
        path = "src/engine/analysis/primer_specificity_matrix_tests.rs"
        regression = (ROOT / path).read_text()
        self.assertIn('include_str!("../../../docs/tutorial/04-09_reference_bound_primer_specificity.md")', regression)
        self.assertIn('"gentle_cli --state study.gentle.json "', regression)
        self.assertIn('template.replace(', regression)
        self.assertIn('execute_shell_command(&mut engine, &command)', regression)
        self.assertIn("primers list-transcript-assay-panels", self.page)
        self.assertIn("Hands-on inspection exercise", self.page)


if __name__ == "__main__":
    unittest.main()
