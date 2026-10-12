import gzip
import json
from pathlib import Path
import tempfile
import unittest

from scripts import prepare_regulatory_region_indexes as prep


class RegulatoryRegionIndexPreparationTests(unittest.TestCase):
    @staticmethod
    def _region_set(region_id="ensembl_feature", assembly="GRCh38", taxon_id=9606):
        return {
            "schema": prep.GENOMIC_REGION_SET_SCHEMA,
            "set_id": "source_regions",
            "content_sha256": "sha256:declared-by-gentle",
            "regions": [{
                "schema": prep.GENOMIC_REGION_SCHEMA,
                "region_id": region_id,
                "label": "Typed promoter candidate",
                "interval": {
                    "reference": {
                        "species_scientific_name": "Homo sapiens",
                        "taxon_id": taxon_id,
                        "assembly_name": assembly,
                        "assembly_accession": "GCA_000001405.29",
                        "contig_name": "17",
                    },
                    "start_0based": 100,
                    "end_0based_exclusive": 180,
                    "strand": "minus",
                    "coordinate_convention": "zero_based_half_open",
                },
                "purpose": "promoter_region",
                "selection_method": "ensembl_regulatory_feature",
                "identity_sha256": "sha256:identity",
                "content_sha256": "sha256:content",
            }],
        }

    def test_canonical_kmers_are_strand_neutral(self):
        forward = prep.canonical_kmers("AACCGGTTACGA", 5)
        reverse = prep.canonical_kmers(prep.reverse_complement("AACCGGTTACGA"), 5)
        self.assertEqual(forward, reverse)

    def test_ambiguous_words_are_excluded(self):
        words = prep.canonical_kmers("AAAANCCCC", 4)
        self.assertIn("AAAA", words)
        self.assertIn("CCCC", words)
        self.assertTrue(all("N" not in word for word in words))

    def test_pairwise_resources_are_deterministic(self):
        records = [("left", "AACCGGTTAACC"), ("right", "AACCGGAAAACC")]
        with tempfile.TemporaryDirectory() as first_dir, tempfile.TemporaryDirectory() as second_dir:
            first = Path(first_dir)
            second = Path(second_dir)
            first_hashes = prep.write_kmer_resources(records, [3, 5], first)
            second_hashes = prep.write_kmer_resources(records, [3, 5], second)
            self.assertEqual(first_hashes, second_hashes)
            self.assertEqual((first / "kmer_pairwise.tsv").read_bytes(),
                             (second / "kmer_pairwise.tsv").read_bytes())
            with gzip.open(first / "kmer_signatures.json.gz", "rb") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["k_values"], [3, 5])
            self.assertIn("not establish", payload["interpretation"])

    def test_single_fasta_rejects_multiple_records(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "two.fa"
            path.write_text(">one\nACGT\n>two\nTGCA\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "one FASTA record"):
                prep.parse_single_fasta(path)

    def test_digest_is_sequence_content_bound(self):
        self.assertEqual(prep.sha256_bytes(b"ACGT"), prep.sha256_bytes(b"ACGT"))
        self.assertNotEqual(prep.sha256_bytes(b"ACGT"), prep.sha256_bytes(b"ACGA"))

    def test_extraction_copy_excludes_runtime_timestamp(self):
        state = {"metadata": {"provenance": {"genome_extractions": [{
            "seq_id": "roi", "chromosome": "7", "recorded_at_unix_ms": 123,
        }]}}}
        extraction = prep.latest_genome_extraction(state, "roi")
        self.assertEqual(extraction["chromosome"], "7")
        self.assertNotIn("recorded_at_unix_ms", extraction)
        self.assertIn("recorded_at_unix_ms", state["metadata"]["provenance"]["genome_extractions"][0])

    def test_equivalent_regions_are_indexed_once_but_all_mapped(self):
        indexed, classes, class_by_region = prep.collapse_equivalent_records([
            ("tx_a", "AACCGG"), ("tx_b", "AACCGG"), ("tx_c", "TTGGCC"),
        ])
        self.assertEqual([record_id for record_id, _ in indexed], ["tx_a", "tx_c"])
        self.assertEqual(classes[0]["member_region_ids"], ["tx_a", "tx_b"])
        self.assertEqual(class_by_region["tx_a"], class_by_region["tx_b"])
        self.assertNotEqual(class_by_region["tx_a"], class_by_region["tx_c"])

    def test_duplicate_region_ids_fail_closed(self):
        with self.assertRaisesRegex(RuntimeError, "region_id values must be unique"):
            prep.collapse_equivalent_records([("dup", "AAAA"), ("dup", "CCCC")])

    def test_canonical_region_set_is_filtered_and_coordinates_are_converted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "ensembl.json"
            source.write_text(json.dumps(self._region_set()), encoding="utf-8")
            manifest_path = root / "comparison.json"
            manifest = {
                "region_sets": [{
                    "path": "ensembl.json",
                    "comparison_class": "ensembl_regulation",
                    "id_prefix": "ens_",
                    "purposes": ["promoter_region"],
                    "selection_methods": ["ensembl_regulatory_feature"],
                }]
            }
            tasks = prep.canonical_region_tasks(manifest, manifest_path, {
                "expected_reference": {"taxon_id": 9606, "assembly_names": ["GRCh38"]}
            })
            self.assertEqual(len(tasks), 1)
            self.assertEqual(tasks[0]["region_id"], "ens_ensembl_feature")
            self.assertEqual(tasks[0]["start_1based"], 101)
            self.assertEqual(tasks[0]["end_1based"], 180)
            self.assertEqual(tasks[0]["source_region"]["interval"]["strand"], "minus")
            self.assertEqual(tasks[0]["source_region_set"]["file_sha256"], prep.sha256_file(source))

    def test_canonical_region_set_rejects_wrong_assembly(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "wrong.json"
            source.write_text(json.dumps(self._region_set(assembly="GRCh37")), encoding="utf-8")
            manifest = {"region_sets": [{
                "path": str(source), "comparison_class": "self_defined",
            }]}
            with self.assertRaisesRegex(RuntimeError, "reference assembly is not allowed"):
                prep.canonical_region_tasks(manifest, root / "manifest.json", {
                    "expected_reference": {"taxon_id": 9606, "assembly_names": ["GRCh38"]}
                })

    def test_canonical_region_sets_require_machine_checkable_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "regions.json"
            source.write_text(json.dumps(self._region_set()), encoding="utf-8")
            manifest = {"region_sets": [{
                "path": str(source), "comparison_class": "self_defined",
            }]}
            with self.assertRaisesRegex(RuntimeError, "expected_reference is required"):
                prep.canonical_region_tasks(manifest, root / "manifest.json", {})


if __name__ == "__main__":
    unittest.main()
