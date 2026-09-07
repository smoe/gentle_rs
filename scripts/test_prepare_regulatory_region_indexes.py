import gzip
import json
from pathlib import Path
import tempfile
import unittest

from scripts import prepare_regulatory_region_indexes as prep


class RegulatoryRegionIndexPreparationTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
