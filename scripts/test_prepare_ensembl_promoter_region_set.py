import json
from pathlib import Path
import tempfile
import unittest

from scripts import prepare_ensembl_promoter_region_set as prep


class EnsemblPromoterRegionSetPreparationTests(unittest.TestCase):
    @staticmethod
    def _report(*, verified=True, truncated=False):
        source = {
            "source_id": "ensembl_regulation_2026_08_grch38",
            "species_scientific_name": "Homo sapiens",
            "taxon_id": 9606,
            "assembly_name": "GRCh38",
            "assembly_accession": "GCA_000001405.29",
        }
        return {
            "schema": prep.REPORT_SCHEMA,
            "seq_id": "gene_locus",
            "gene_symbol": "GENE1",
            "isoform_evidence": {"chromosome": "17"},
            "ensembl_regulation": {
                "availability": "available",
                "evidence_statement": "Provider annotation only.",
                "non_claims": ["No biosample activity is implied."],
                "source_binding": {
                    "source": source,
                    "content_identity_verified": verified,
                    "truncated": truncated,
                    "index_sha256": "sha256:index",
                    "intervals_sha256": "sha256:intervals",
                },
                "rows": [
                    {
                        "source_id": source["source_id"],
                        "feature_id": "ENSR000001",
                        "feature_type": "promoter",
                        "core_genomic_start_1based": 101,
                        "core_genomic_end_1based": 180,
                    },
                    {
                        "source_id": source["source_id"],
                        "feature_id": "ENSR000002",
                        "feature_type": "enhancer",
                        "core_genomic_start_1based": 250,
                        "core_genomic_end_1based": 310,
                    },
                ],
            },
        }

    def _write_report(self, root: Path, **kwargs) -> Path:
        path = root / "report.json"
        path.write_text(json.dumps(self._report(**kwargs)), encoding="utf-8")
        return path

    def test_only_requested_feature_type_is_captured_with_exact_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report_path = self._write_report(root)
            requests = prep.capture_requests(report_path, "promoters", {"promoter"})
            self.assertEqual(len(requests), 1)
            request = requests[0]
            self.assertEqual(request["region_id"], "GENE1_ENSR000001_promoter")
            self.assertEqual(request["purpose"], "promoter_region")
            source = request["source"]
            self.assertEqual(source["source_kind"], "gene_locus_ensembl_regulatory_feature")
            self.assertEqual(source["row"]["feature_id"], "ENSR000001")
            self.assertTrue(source["source_binding"]["content_identity_verified"])
            self.assertEqual(source["reference"]["contig_name"], "17")
            self.assertIn("sha256:", request["notes"][1])

    def test_unverified_source_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            report_path = self._write_report(Path(directory), verified=False)
            with self.assertRaisesRegex(RuntimeError, "content identity is not verified"):
                prep.capture_requests(report_path, "promoters", {"promoter"})

    def test_truncated_overlap_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            report_path = self._write_report(Path(directory), truncated=True)
            with self.assertRaisesRegex(RuntimeError, "rows are truncated"):
                prep.capture_requests(report_path, "promoters", {"promoter"})


if __name__ == "__main__":
    unittest.main()
