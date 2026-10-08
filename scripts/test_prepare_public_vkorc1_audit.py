"""Hand-crafted inline tests for the public audit helper's narrow assertions."""

import unittest
from pathlib import Path

from scripts.prepare_public_vkorc1_audit import (
    assert_reviewed_fragment, pair_content, report_content, retained_refsnp_url, PREFIX,
)


class PublicAuditAssertionsTests(unittest.TestCase):
    def project(self, reference="aaCaa", alternate="aaTaa", fragment="aacaa"):
        return {"sequences": {PREFIX + role: {"seq": {"seq": list(value.encode())}}
                              for role, value in (("fragment", fragment), ("reference", reference), ("alternate", alternate))}}

    def test_retained_file_override_preserves_raw_spaces_quotes_and_backslashes(self):
        path = Path('/tmp/public response/ref"snp\\record.json')
        self.assertEqual(retained_refsnp_url(path), 'file:///tmp/public response/ref"snp\\record.json')

    def test_engine_recommendation_must_match_reviewed_tutorial_geometry(self):
        report = {"recommended_candidate_id": "chosen", "candidates": [
            {"candidate_id": "chosen", "recommended": True, "start_0based": 10, "end_0based_exclusive": 20}]}
        assert_reviewed_fragment(report, {"ExtractRegion": {"from": 10, "to": 20}})
        with self.assertRaises(AssertionError):
            assert_reviewed_fragment(report, {"ExtractRegion": {"from": 10, "to": 21}})

    def test_only_known_execution_identity_is_excluded_from_report_comparison(self):
        report = {"generated_at_unix_ms": 1, "op_id": "op-1", "run_id": "run",
                  "genomic_alt": "A,G,T", "promoter_windows": [2, 1], "unknown": "retain"}
        self.assertEqual(report_content(report), {"genomic_alt": "A,G,T", "promoter_windows": [2, 1], "unknown": "retain"})

    def test_exact_pair_preserves_raw_nonvariant_bytes_without_scientific_claims(self):
        result = pair_content(self.project())
        self.assertEqual(result["differences"], [{"position_0based": 2, "reference": "C", "alternate": "T"}])
        self.assertEqual(result["raw_sequences"]["fragment"], "aacaa")
        self.assertNotIn("synthetic", result, "Pair shape alone does not establish public provenance")
        self.assertFalse(result["native_screenshot"])
        self.assertFalse(result["human_scientific_approval"])

    def test_wrong_allele_length_second_difference_and_changed_nonvariant_case_fail(self):
        for project in (self.project(alternate="aaGaa"), self.project(alternate="aaTaaa"),
                        self.project(alternate="aaTta"), self.project(reference="AaCaa", alternate="AaTaa")):
            with self.subTest(project=project), self.assertRaises(AssertionError):
                pair_content(project)


if __name__ == "__main__":
    unittest.main()
