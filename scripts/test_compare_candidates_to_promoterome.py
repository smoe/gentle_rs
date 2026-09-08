import unittest

from scripts import compare_candidates_to_promoterome as compare


class CandidatePromoteromeComparisonTests(unittest.TestCase):
    def test_merge_intervals_collapses_overlap_and_touching_spans(self):
        self.assertEqual(compare.merge_intervals([(5, 10), (1, 4), (4, 7), (20, 22)]),
                         [(1, 10), (20, 22)])

    def test_query_interval_supports_canonical_and_extracted_regions(self):
        canonical = {"input_kind": "canonical_region", "source_region": {"interval": {
            "reference": {"contig_name": "7"}, "start_0based": 100, "end_0based_exclusive": 180,
        }}}
        extracted = {"input_kind": "transcript_promoter", "genome_extraction": {
            "chromosome": "11", "start_1based": 501, "end_1based": 700,
        }}
        self.assertEqual(compare.query_interval(canonical), ("7", 100, 180))
        self.assertEqual(compare.query_interval(extracted), ("11", 500, 700))

    def test_equivalent_index_query_is_expanded_to_every_declared_region(self):
        rows = [{"qseqid": "roi_a", "sseqid": "target"}]
        expanded = compare.expand_equivalent_query_rows(rows, [{
            "representative_region_id": "roi_a",
            "member_region_ids": ["roi_a", "roi_b"],
        }], {"roi_a"})
        self.assertEqual([row["qseqid"] for row in expanded], ["roi_a", "roi_b"])
        self.assertTrue(all(row["indexed_representative_query_id"] == "roi_a" for row in expanded))

    def test_contig_alias_normalization_handles_chr_prefix_only(self):
        self.assertEqual(compare.normalized_contig("chr7"), "7")
        self.assertEqual(compare.normalized_contig("7"), "7")

    def test_coverage_counts_windows_genes_and_transcripts_separately(self):
        segments = compare.coverage_segments(
            {"p1": [(0, 10)], "p2": [(5, 15)]},
            {"p1": {"g1"}, "p2": {"g1", "g2"}},
            {"p1": {"t1", "t2"}, "p2": {"t3"}},
        )
        self.assertEqual(segments, [
            {"query_start_0based": 0, "query_end_0based_exclusive": 5,
             "distinct_promoter_windows": 1, "distinct_genes": 1, "distinct_transcripts": 2},
            {"query_start_0based": 5, "query_end_0based_exclusive": 10,
             "distinct_promoter_windows": 2, "distinct_genes": 2, "distinct_transcripts": 3},
            {"query_start_0based": 10, "query_end_0based_exclusive": 15,
             "distinct_promoter_windows": 1, "distinct_genes": 2, "distinct_transcripts": 1},
        ])


if __name__ == "__main__":
    unittest.main()
