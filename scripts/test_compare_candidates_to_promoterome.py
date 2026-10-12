import unittest
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts import compare_candidates_to_promoterome as compare


class CandidatePromoteromeComparisonTests(unittest.TestCase):
    """Hand-crafted tiny HSPs exercise accounting only, not biological similarity."""

    def test_caps_are_detected_before_filters(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "hits.tsv"
            # Both targets and two HSPs for t1 reach the limits, but every row
            # fails the scientific length/identity filter.
            path.write_text("\n".join([
                "q\tt1\t20\t5\t0\t0\t1\t5\t1\t5\t1\t10",
                "q\tt1\t20\t5\t0\t0\t6\t10\t6\t10\t1\t10",
                "q\tt2\t20\t5\t0\t0\t1\t5\t1\t5\t1\t10",
            ]) + "\n")
            audit = {}
            self.assertEqual(compare.parse_blast(path, 40, 80, 1e-5, limits=(2, 2), audit=audit), [])
            self.assertEqual(audit["q"], {"raw_target_count": 2, "max_hsps_observed": 2,
                                          "target_cap_reached": True, "hsp_cap_reached": True,
                                          "counts_are_lower_bounds": True})

    def test_shared_gene_windows_are_excluded_from_counts_tiers_and_strips(self):
        queries = {"q": {"sequence_length_bp": 100, "genome_extraction": {"gene_id": "g1"}}}
        windows = {key: {"chromosome": "2", "start_0based": "1000", "end_0based_exclusive": "1100",
                         "tss_1based": "1001", "strand": "+"} for key in ("own", "shared", "other")}
        mappings = [{"promoter_id": key, "gene_id": gene, "gene_name": "same_display_name", "transcript_id": tx}
                    for key, gene, tx in [("own", "g1", "t1"), ("shared", "g1", "t2"),
                                          ("shared", "g2", "t3"), ("other", "g2", "t4")]]
        rows = [{"qseqid": "q", "sseqid": key, "qstart": 1, "qend": 50, "pident": 90, "bitscore": 100}
                for key in windows]
        with TemporaryDirectory() as directory:
            target = Path(directory) / "matches.tsv"
            summary = compare.summarize_task(rows, queries, windows, mappings, target)[0]
            with target.open() as handle:
                matches = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(summary["other_promoters"], {"distinct_promoter_windows": 1, "distinct_genes": 1,
                                                      "distinct_transcripts": 1})
        self.assertEqual(summary["other_promoters_by_min_query_coverage"]["0.50"]["distinct_genes"], 1)
        self.assertEqual(summary["recurrent_query_segments"][0]["distinct_promoter_windows"], 1)
        self.assertEqual({row["promoter_id"] for row in matches if row["same_gene"] == "False"}, {"other"})
        self.assertTrue(summary["other_gene_exclusion_verified"])

    def test_unassigned_query_gene_is_not_marked_verified(self):
        with TemporaryDirectory() as directory:
            summary = compare.summarize_task([], {"q": {"sequence_length_bp": 100}}, {}, [],
                                             Path(directory) / "matches.tsv")[0]
        self.assertFalse(summary["other_gene_exclusion_verified"])
        self.assertEqual(summary["other_promoters"]["distinct_genes"], 0)

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
