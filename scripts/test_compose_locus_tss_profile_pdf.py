#!/usr/bin/env python3

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compose_locus_tss_profile_pdf as target


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class CompositeLocusTssPdfTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.locus_svg = self.root / "GENE_locus.svg"
        self.locus_svg.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" '
            'data-gentle-schema="gentle.gene_locus_evidence_display.v1" '
            'data-gentle-panel-id="GENE_panel" width="1400" height="100">'
            '<rect data-gentle-tss-stretch-band="GENE_tss_stretch_1" '
            'data-gentle-genomic-start="900" data-gentle-genomic-end="1100"/>'
            '</svg>'
        )
        self.locus_receipt = self.root / "locus.receipt.json"
        self.locus_receipt.write_text(json.dumps({
            "schema": target.LOCUS_RECEIPT_SCHEMA,
            "gene": "GENE",
            "outputs": {self.locus_svg.name: f"sha256:{digest(self.locus_svg)}"},
        }))
        self.selected_svg = self.root / "selected.svg"
        self.selected_svg.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="30" '
            'data-gentle-plot-left="255" data-gentle-plot-right="1050"/>'
        )
        self.other_svg = self.root / "other.svg"
        self.other_svg.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="30" '
            'data-gentle-plot-left="255" data-gentle-plot-right="1050"/>'
        )
        self.sequences = {"selected": "A" * 701, "other": "C" * 701}
        self.report = self.root / "report.json"
        self.report.write_text(json.dumps({
            "schema": target.TSS_REPORT_SCHEMA,
            "score_policy": {"score_kind": "llr_background_tail_log10"},
            "reference": {"assembly": "GRCh38"},
            "verification": "bundle_consistency_verified; prepared_reference_not_assessed",
            "non_claims": "predictions are not measured binding",
            "windows": [self.window("selected", 1000, True), self.window("other", 1050, False)],
        }))
        self.index = self.root / "index.json"
        self.index.write_text(json.dumps({
            "schema": target.TSS_INDEX_SCHEMA,
            "report_sha256": digest(self.report),
            "genes": [{
                "gene_symbol": "GENE",
                "pages": [
                    {"promoter_ids": ["selected"], "files": [self.selected_svg.name]},
                    {"promoter_ids": ["other"], "files": [self.other_svg.name]},
                ],
            }],
        }))
        self.fasta = self.root / "GENE_TSS_minus500_plus200.fasta"
        self.fasta.write_text("".join(
            f">GENE|gene_id=ENSG_TEST|promoter_id={promoter_id}|assembly=GRCh38|"
            f"chromosome=1|strand=+|tss_1based={tss}|genomic_1based={tss - 500}-{tss + 200}|"
            "window=minus500_plus200|orientation=transcript_5prime_to_3prime|"
            f"transcripts=ENST_TEST|sequence_sha256={hashlib.sha256(self.sequences[promoter_id].encode()).hexdigest()}\n"
            + self.sequences[promoter_id] + "\n"
            for promoter_id, tss in [("selected", 1000), ("other", 1050)]
        ))
        self.checksums = self.root / "SHA256SUMS"
        self.checksums.write_text(f"{digest(self.fasta)}  {self.fasta.name}\n")
        self.manifest = self.root / "manifest.json"
        self.manifest.write_text(json.dumps({
            "schema": target.TARGET_FASTA_SCHEMA,
            "sha256sums_sha256": digest(self.checksums),
            "files": [{
                "filename": self.fasta.name,
                "gene_id": "ENSG_TEST",
                "gene_symbol": "GENE",
                "sha256": digest(self.fasta),
                "record_count": 2,
                "records": [
                    self.manifest_record("selected", 1000),
                    self.manifest_record("other", 1050),
                ],
            }],
        }))
        self.tss_receipt = self.root / "tss.receipt.json"
        self.write_tss_receipt()
        self.cli = self.root / "gentle_cli"
        self.cli.write_bytes(b"synthetic executable")

    def tearDown(self):
        self.temp.cleanup()

    def window(self, promoter_id, tss, selected):
        return {
            "selected": selected,
            "selection_evidence": ({
                "label": "Selected in integrated report",
                "factor": "TP73",
                "criterion": "synthetic criterion",
            } if selected else None),
            "record": {
                "promoter_id": promoter_id,
                "gene_id": "ENSG_TEST",
                "gene_symbol": "GENE",
                "geometry": {
                    "chromosome": "1", "strand": "+", "tss_1based": tss,
                    "start_1based": tss - 500, "end_1based": tss + 200,
                },
                "sequence_sha256": hashlib.sha256(self.sequences[promoter_id].encode()).hexdigest(),
                "transcripts": ["ENST_TEST"],
            },
        }

    def manifest_record(self, promoter_id, tss):
        return {
            "promoter_id": promoter_id,
            "gene_id": "ENSG_TEST",
            "chromosome": "1",
            "strand": "+",
            "tss_1based": tss,
            "genomic_start_1based": tss - 500,
            "genomic_end_1based": tss + 200,
            "sequence_length_bp": 701,
            "sequence_sha256": hashlib.sha256(self.sequences[promoter_id].encode()).hexdigest(),
            "transcript_ids": ["ENST_TEST"],
        }

    def write_tss_receipt(self):
        self.tss_receipt.write_text(json.dumps({
            "schema": target.TSS_RECEIPT_SCHEMA,
            "inputs": [
                {"role": "bundle_manifest", "name": self.manifest.name,
                 "sha256": digest(self.manifest)},
                {"role": "bundle_checksums", "name": self.checksums.name,
                 "sha256": digest(self.checksums)},
                {"role": "fasta", "name": self.fasta.name,
                 "sha256": digest(self.fasta)},
            ],
            "outputs": {
                self.report.name: digest(self.report),
                self.index.name: digest(self.index),
                self.selected_svg.name: digest(self.selected_svg),
                self.other_svg.name: digest(self.other_svg),
            },
        }))

    def args(self):
        return argparse.Namespace(
            gene="GENE", locus_svg=self.locus_svg, locus_receipt=self.locus_receipt,
            tss_report=self.report, tss_index=self.index, tss_receipt=self.tss_receipt,
            tss_bundle_manifest=self.manifest,
            gentle_cli=self.cli, producer_revision="test-revision",
            output_pdf=self.root / "combined.pdf",
            output_fasta=self.root / "selected.fasta",
            output_receipt=self.root / "combined.receipt.json",
        )

    def test_compose_binds_context_and_only_selected_tss(self):
        def fake_run(command, **kwargs):
            Path(command[2]).write_bytes(b"%PDF-1.4\nsynthetic")
            return subprocess.CompletedProcess(
                command, 0,
                stdout=json.dumps({"page_count": 2, "pages": [{}, {}]}), stderr="",
            )
        with mock.patch.object(target.subprocess, "run", side_effect=fake_run):
            receipt = target.compose(self.args())
        self.assertEqual(receipt["page_order"], ["locus_context", "selected_tss_tfbs"])
        self.assertEqual(receipt["selected_tss_count"], 1)
        self.assertEqual(receipt["selected_tss_bindings"][0]["promoter_id"], "selected")
        self.assertEqual(receipt["selected_tss_bindings"][0]["covering_tss_bands"],
                         ["GENE_tss_stretch_1"])
        self.assertTrue(self.args().output_pdf.is_file())
        self.assertTrue(self.args().output_fasta.is_file())
        self.assertTrue(self.args().output_receipt.is_file())
        fasta = self.args().output_fasta.read_text()
        self.assertIn("promoter_id=selected", fasta)
        self.assertNotIn("promoter_id=other", fasta)
        self.assertEqual(receipt["horizontal_alignment"], {
            "page_width": 1400.0,
            "plot_left": 255.0,
            "plot_right": 1050.0,
            "policy": "context and detail pages share the canonical locus plot frame",
        })

    def test_rejects_selected_tss_outside_bound_locus_bands(self):
        value = json.loads(self.report.read_text())
        value["windows"][0]["record"]["geometry"]["tss_1based"] = 2000
        self.report.write_text(json.dumps(value))
        self.index.write_text(json.dumps({**json.loads(self.index.read_text()),
                                          "report_sha256": digest(self.report)}))
        self.write_tss_receipt()
        with self.assertRaisesRegex(ValueError, "outside every locus band"):
            target.compose(self.args())

    def test_rejects_tampered_selected_page_before_rendering(self):
        self.selected_svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="1401" height="30"/>')
        with self.assertRaisesRegex(ValueError, "output hash mismatch"):
            target.compose(self.args())

    def test_rejects_tampered_source_fasta_before_rendering(self):
        self.fasta.write_text(self.fasta.read_text().replace("A", "T", 1))
        with self.assertRaisesRegex(ValueError, "FASTA hash does not match"):
            target.compose(self.args())


if __name__ == "__main__":
    unittest.main()
