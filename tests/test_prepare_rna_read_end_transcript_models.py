"""GTF-to-read-end-model preparation regressions."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("prepare_rna_read_end_transcript_models", ROOT / "scripts" / "prepare_rna_read_end_transcript_models.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
ANNOTATION = "sha256:" + "a" * 64


class PrepareTests(unittest.TestCase):
    def test_gtf_coordinates_and_provenance_are_frozen(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf, fasta, output = root / "models.gtf", root / "ref.fa", root / "models.json"
            gtf.write_text(
                'chr1\tsrc\texon\t101\t110\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
                'chr1\tsrc\tCDS\t106\t110\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
                'chr1\tsrc\texon\t121\t140\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
                'chr1\tsrc\tCDS\t121\t130\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n',
                encoding="utf-8",
            )
            fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
            MODULE.run(gtf, fasta, "synthetic-v1", "annotation1", ANNOTATION, output)
            document = json.loads(output.read_text(encoding="utf-8"))
            transcript = document["transcripts"][0]
            self.assertEqual(transcript["genomic_exons_0based_half_open"], [[100, 110], [120, 140]])
            self.assertEqual(transcript["genomic_cds_0based_half_open"], [[105, 110], [120, 130]])
            self.assertEqual(document["annotation_sha256"], ANNOTATION)
            self.assertEqual(document["model_source"]["name"], "models.gtf")

    def test_inconsistent_transcript_identity_fails(self):
        text = (
            'chr1\ts\texon\t1\t10\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
            'chr1\ts\tCDS\t2\t9\t.\t+\t0\tgene_id "g2"; transcript_id "t1";\n'
        )
        with self.assertRaisesRegex(ValueError, "inconsistent transcript identity"):
            MODULE.parse_gtf(text)


if __name__ == "__main__":
    unittest.main()
