"""Inline hand-crafted JSON tests for the synthetic CI base-comparison boundary."""

import unittest

from scripts.ci_vkorc1_gui_audit import pair_evidence


class PairEvidenceTests(unittest.TestCase):
    def project(self, source="aaaaaacaaaaaaaaaaaaa", reference="aaaaaaCaaaaaaaaaaaaa",
                alternate="aaaaaaTaaaaaaaaaaaaa"):
        return {"sequences": {"vkorc1_rs9923231_promoter_" + role: {"seq": {"seq": list(bases.encode())}}
                              for role, bases in (("fragment", source), ("reference", reference), ("alternate", alternate))}}

    def test_base_proof_is_explicitly_synthetic_and_not_a_screenshot(self):
        result = pair_evidence(self.project())
        self.assertEqual(result["differences"], [{"position_0based": 6, "reference": "C", "alternate": "T"}])
        self.assertTrue(result["synthetic"])
        self.assertFalse(result["native_screenshot"])
        self.assertFalse(result["online_vkorc1_accepted"])
        self.assertFalse(result["human_scientific_approval"])
        self.assertEqual(result["raw_sequences"]["fragment"], "aaaaaacaaaaaaaaaaaaa")
        self.assertEqual(result["raw_sequences"]["reference"], "aaaaaaCaaaaaaaaaaaaa")
        self.assertEqual(result["sequences"]["reference"], "AAAAAACAAAAAAAAAAAAA")

    def test_changed_source_missing_output_wrong_allele_and_length_fail_closed(self):
        for project in (self.project(source="A" * 20), self.project(alternate="A" * 20),
                        self.project(reference="AaaaaaCaaaaaaaaaaaaa"),
                        self.project(alternate="aaaaaaTaaaaaaaaaaaaaa"), {"sequences": {}}):
            with self.subTest(project=project), self.assertRaises((AssertionError, KeyError)):
                pair_evidence(project)


if __name__ == "__main__":
    unittest.main()
