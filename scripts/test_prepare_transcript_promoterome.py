from pathlib import Path
import tempfile
import unittest

from scripts import prepare_transcript_promoterome as prep


class TranscriptPromoteromePreparationTests(unittest.TestCase):
    def test_plus_and_minus_intervals_match_gentle_geometry(self):
        plus = {"strand": "+", "transcript_start_1based": 1000, "transcript_end_1based": 1500}
        minus = {"strand": "-", "transcript_start_1based": 1000, "transcript_end_1based": 1500}
        self.assertEqual(prep.promoter_interval(plus, 100, 20, 5000), (899, 1020, 1000, False))
        self.assertEqual(prep.promoter_interval(minus, 100, 20, 5000), (1479, 1600, 1500, False))

    def test_boundary_clipping_is_explicit(self):
        plus = {"strand": "+", "transcript_start_1based": 50, "transcript_end_1based": 100}
        self.assertEqual(prep.promoter_interval(plus, 100, 20, 1000), (0, 70, 50, True))

    def test_shared_tss_window_is_indexed_once_but_maps_every_transcript(self):
        transcripts = [
            {"chromosome": "1", "strand": "+", "transcript_start_1based": 1000,
             "transcript_end_1based": 1500, "gene_id": "G1", "gene_name": "ONE", "transcript_id": "T1"},
            {"chromosome": "1", "strand": "+", "transcript_start_1based": 1000,
             "transcript_end_1based": 1700, "gene_id": "G1", "gene_name": "ONE", "transcript_id": "T2"},
            {"chromosome": "1", "strand": "+", "transcript_start_1based": 2000,
             "transcript_end_1based": 2200, "gene_id": "G2", "gene_name": "TWO", "transcript_id": "T3"},
        ]
        windows, mappings, excluded = prep.build_windows(transcripts, {"1": 5000}, 100, 20)
        self.assertEqual(len(windows), 2)
        self.assertEqual(len(mappings), 3)
        shared = next(row for row in windows if row["tss_1based"] == 1000)
        self.assertEqual(shared["gene_count"], 1)
        self.assertEqual(shared["transcript_count"], 2)
        self.assertEqual(excluded, {"invalid_record": 0, "contig_absent": 0})

    def test_bedtools_headers_are_normalized_and_validated(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "out.fa"
            path.write_text(">promoter_a(+)\naacg\n>promoter_b(-)\nttgc\n", encoding="utf-8")
            prep.normalize_bedtools_fasta(path, {"promoter_a", "promoter_b"})
            self.assertEqual(path.read_text(encoding="utf-8"), ">promoter_a\nAACG\n>promoter_b\nTTGC\n")


if __name__ == "__main__":
    unittest.main()
