"""Synthetic inline fixture tests for the BAM/FASTQ-to-shared-contract adapter.

Every sequence, quality, SAM/BAM row, GTF and metadata file is invented below in
private temporary directories. Cell-line/gene labels only exercise fixture
selection; they do not supply or predict biological results. Recreate with
`python -m unittest scripts.test_prepare_cdna_fixture_comparison`; optional
pysam/CLI process coverage needs those explicitly supplied external tools.
"""

import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from scripts import prepare_cdna_fixture_comparison as prep

try:
    import pysam
except ImportError:
    pysam = None


class PreparationHelpers(unittest.TestCase):
    def test_inventory_checks_raw_lf_crlf_and_exact_population(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            path = root / "input's.txt"
            for raw in (b"x\ny\n", b"x\r\ny\r\n"):
                path.write_bytes(raw)
                (root / "SHA256SUMS").write_text(hashlib.sha256(raw).hexdigest() + "  ./input's.txt\n", encoding="ascii")
                bound = prep.verify_inventory(root)[0]
                self.assertEqual(bound["sha256"], prep.sha_bytes(raw))
                (root / "input's.txt").write_bytes(raw + b" ")
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    prep.verify_inventory(root)
                (root / "input's.txt").unlink()

    def test_inventory_rejects_unsafe_and_unlisted_paths(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            (root / "a").write_bytes(b"a")
            digest = hashlib.sha256(b"a").hexdigest()
            for path in ("../a", "/a", "x//a", "x/./a", "x\\a"):
                (root / "SHA256SUMS").write_text(digest + "  " + path + "\n")
                with self.assertRaises(ValueError):
                    prep.verify_inventory(root)
            (root / "SHA256SUMS").write_text(digest + "  a\n")
            (root / "unlisted").write_text("extra")
            with self.assertRaisesRegex(ValueError, "exactly"):
                prep.verify_inventory(root)

    def test_region_headers_and_raw_coordinate_offsets_are_not_contig_aliases(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch) / "regions.fa"
            path.write_bytes(b">17:101-104\r\nAcgT\r\n>1:201-204\r\nNNNN\r\n")
            self.assertEqual(prep.region_sequences(path), {"17:101-104": "ACGT", "1:201-204": "NNNN"})
            path.write_text(">17:101-104\nACGT\n>17:101-104\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                prep.region_sequences(path)

    def test_fastq_quality_identity_and_duplicate_refuse(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch) / "a.fastq.gz"
            for raw in (b"@r\nAC\n+\nI\n", b"@r\nAC\n+\nII\n@r\nAC\n+\nII\n", b"@r\nAC\n+\nI \n"):
                path.write_bytes(gzip.compress(raw))
                with self.assertRaises(ValueError):
                    prep.read_fastq(path)
            path.write_bytes(gzip.compress(b"@r\r\nAC\r\n+\r\nI&\r\n"))
            self.assertEqual(prep.read_fastq(path), {"r": (1, "AC", "I&")})

    def test_original_axis_and_quality_equality_are_mandatory(self):
        row = SimpleNamespace(query_length=2, infer_read_length=lambda: 2,
                              get_forward_sequence=lambda: "AC", get_forward_qualities=lambda: [40, 5])
        self.assertEqual(prep.matched_original(row, (1, "AC", "I&")), ("AC", "I&"))
        for record, fastq in ((row, (1, "AC", "II")), (row, (1, "GT", "I&")),
                              (SimpleNamespace(**{**row.__dict__, "infer_read_length": lambda: 3}), (1, "AC", "I&")),
                              (SimpleNamespace(**{**row.__dict__, "get_forward_qualities": lambda: None}), (1, "AC", "I&"))):
            with self.assertRaises(ValueError):
                prep.matched_original(record, fastq)

    def test_external_command_uses_argv_and_preserves_caller_environment(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            response = SimpleNamespace(returncode=0, stdout='{"result":{}}', stderr="")
            with patch.dict(os.environ, {"RUST_MIN_STACK": "99999999"}), patch.object(prep.subprocess, "run", return_value=response) as run:
                self.assertEqual(prep.invoke(root / "tool's α", root, "rna-reads", "@a b"), {"result": {}})
            args, kwargs = run.call_args
            self.assertEqual(args[0][-2:], ["rna-reads", "@a b"])
            self.assertNotIn("env", kwargs)
            self.assertNotIn("shell", kwargs)


@unittest.skipUnless(pysam, "requires explicitly supplied external pysam")
class BamEvidenceTests(unittest.TestCase):
    def test_partition_units_preserve_original_ordinals_and_reset_only_sam_row_numbers(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            sam = root / "source.sam"
            sam.write_text("@HD\tVN:1.6\n@SQ\tSN:17\tLN:1000\n" + "".join(
                f"q{i}\t0\t17\t101\t60\t2M\t*\t0\t0\tAC\tII\n" for i in range(5)))
            originals = [{"ordinal": i + 1, "qname": f"q{i}", "unit_id": f"bam_record_{i + 1:08d}",
                          "sequence": "AC", "qualities_phred33": "II"} for i in range(5)]
            parts = prep.partition_sam(sam, originals, root, 2)
            self.assertEqual([len(rows) for _, _, rows in parts], [2, 2, 1])
            self.assertEqual([r["original_bam_ordinal"] for _, _, rows in parts for r in rows], [1, 2, 3, 4, 5])
            self.assertEqual([[r["ordinal"] for r in rows] for _, _, rows in parts], [[1, 2], [1, 2], [1]])
            self.assertEqual([r["unit_id"] for _, _, rows in parts for r in rows], [r["unit_id"] for r in originals])
            for directory, path, rows in parts:
                proof = json.loads((directory / "partition-proof.json").read_bytes())
                self.assertEqual(proof["sam"], prep.bound(path))
                self.assertEqual(len(proof["records"]), len(rows))
            with self.assertRaisesRegex(ValueError, "budget"):
                prep.partition_sam(sam, originals * 13, root, 1)

    def test_library_reverse_orientation_preserves_original_quality_order_and_sa(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            bam, fastq, sam = root / "a.bam", root / "a.fastq.gz", root / "a.sam"
            with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "17", "LN": 1000}]}) as out:
                row = pysam.AlignedSegment(out.header)
                row.query_name, row.flag, row.reference_id, row.reference_start = "reverse", 16, 0, 100
                row.query_sequence, row.query_qualities, row.mapping_quality, row.cigarstring = "GT", [5, 40], 60, "2M"
                row.set_tag("SA", "17,201,+,2M,60,0;")
                out.write(row)
            fastq.write_bytes(gzip.compress(b"@reverse\nAC\n+\nI&\n"))
            contigs, rows, proof = prep.convert_bam(bam, fastq, sam)
            self.assertEqual(contigs, [{"name": "17", "length": 1000}])
            self.assertEqual((rows[0]["sequence"], rows[0]["qualities_phred33"]), ("AC", "I&"))
            self.assertEqual(proof[0]["ordinal"], 1)
            self.assertIn("SA:Z:17,201,+,2M,60,0;", sam.read_text())

    def test_hard_clips_and_duplicate_primary_names_are_not_reconstructed_or_counted(self):
        with tempfile.TemporaryDirectory() as scratch:
            root = Path(scratch)
            for index, (cigar, copies) in enumerate((("1H2M", 1), ("2M", 2))):
                bam, fastq = root / f"{index}.bam", root / f"{index}.fq.gz"
                with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "17", "LN": 1000}]}) as out:
                    for _ in range(copies):
                        row = pysam.AlignedSegment(out.header)
                        row.query_name, row.reference_id, row.reference_start = "q", 0, 100
                        row.query_sequence, row.query_qualities, row.mapping_quality, row.cigarstring = "AC", [40, 40], 60, cigar
                        out.write(row)
                fastq.write_bytes(gzip.compress(b"@q\nAC\n+\nII\n"))
                with self.assertRaises(ValueError):
                    prep.convert_bam(bam, fastq, root / f"{index}.sam")


@unittest.skipUnless(pysam and os.environ.get("GENTLE_CDNA_CLI"), "requires external pysam and source-matched GENTLE_CDNA_CLI")
class SharedEngineReplay(unittest.TestCase):
    def test_fixture_preparation_uses_gtf_engine_and_preserves_incomplete_search_scope(self):
        with tempfile.TemporaryDirectory(prefix="cDNA prepared's ") as scratch:
            root = Path(scratch).resolve()
            fixture, output = root / "fixture's α", root / "output's α"
            data = fixture / "data"
            for directory in ("reads", "annotation", "reference"):
                (data / directory).mkdir(parents=True, exist_ok=True)
            (data / "samples.tsv").write_text("barcode\tcell_line\nbarcode15\tSK-MEL-28\nbarcode16\tWM793\n")
            (data / "loci.tsv").write_text("gene\tgene_id\tchromosome\tstart_1based\tend_1based\tstrand\tassembly\n"
                "TP53\tENSG00000141510.21\t17\t101\t140\t-\tGRCh38\n")
            (data / "source_bindings.tsv").write_text("role\tsha256_or_bound_digest\nreference_fasta\t" + "1" * 64 + "\n")
            (data / "reference/ensembl116.tp53_tp73.genomic_regions.fa").write_text(">17:101-140\n" + "A" * 40 + "\n")
            attrs = 'gene_id "ENSG00000141510"; gene_version "21";'
            tx = attrs + ' transcript_id "invented_tx"; transcript_version "1";'
            gtf = f'17\tinvented\tgene\t101\t140\t.\t-\t.\t{attrs}\n'
            for feature, start, end in (("transcript", 101, 140), ("exon", 101, 110), ("exon", 121, 140)):
                gtf += f'17\tinvented\t{feature}\t{start}\t{end}\t.\t-\t.\t{tx}\n'
            (data / "annotation/ensembl116.tp53_tp73.gtf").write_text(gtf)
            for barcode in ("barcode15", "barcode16"):
                prefix = data / "reads" / f"{barcode}.TP53.primary"
                records = []
                with pysam.AlignmentFile(str(prefix) + ".bam", "wb", header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "17", "LN": 1000}]}) as out:
                    for index in range(3):
                        sequence = "AAAC" + "A" * 11 + "G" + "A" * 14 if barcode == "barcode15" else "A" * 30
                        row = pysam.AlignedSegment(out.header)
                        row.query_name, row.reference_id, row.reference_start = f"invented_{barcode}_{index}", 0, 100
                        row.query_sequence, row.query_qualities, row.mapping_quality, row.cigarstring = sequence, [40] * 30, 60, "10M10N20M"
                        out.write(row)
                        records.append(f"@{row.query_name}\n{sequence}\n+\n{'I' * 30}\n")
                Path(str(prefix) + ".fastq.gz").write_bytes(gzip.compress("".join(records).encode()))
            sums = "".join(hashlib.sha256(path.read_bytes()).hexdigest() + "  ./" + path.relative_to(data).as_posix() + "\n" for path in sorted(data.rglob("*")) if path.is_file())
            (data / "SHA256SUMS").write_text(sums)
            receipt = prep.prepare(fixture, Path(os.environ["GENTLE_CDNA_CLI"]), output, "a" * 40, ["TP53"], True,
                                   records_per_partition=1)
            self.assertEqual([row["units"] for row in receipt["genes"][0]["counts"]], [3, 3])
            self.assertEqual([row["partitions"] for row in receipt["genes"][0]["counts"]], [3, 3])
            self.assertEqual([row["units_without_complete_search_declaration"] for row in receipt["genes"][0]["counts"]], [3, 3])
            self.assertEqual([row["units_without_complete_view_declaration"] for row in receipt["genes"][0]["counts"]], [3, 3])
            self.assertEqual([row["units_with_unassessed_supplied_views"] for row in receipt["genes"][0]["counts"]], [0, 0])
            result = json.loads((output / "TP53/comparison.json").read_bytes())
            self.assertEqual([s["position_0based"] for s in result["comparison"]["sites"]], [103, 125])
            annotation = json.loads(Path(result["request"]["normalized_annotation"]["path"]).read_bytes())
            self.assertEqual(annotation["transcripts"][0]["genomic_exons"], [[100, 110], [120, 140]])
            self.assertEqual(annotation["transcripts"][0]["genomic_strand"], "-")
            with self.assertRaises(FileExistsError):
                prep.prepare(fixture, Path(os.environ["GENTLE_CDNA_CLI"]), output, "a" * 40, ["TP53"], True)


if __name__ == "__main__":
    unittest.main()
