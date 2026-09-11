#!/usr/bin/env python3
"""Compositor regressions using only hand-crafted temporary JSON/SVG/FASTA.

The synthetic MA9991.1 identity is not a JASPAR measurement. Recreate the
fixtures by running this module. Rendering is mocked to exercise joins and
publication; Rust SVG/PDF tests cover real rasterization and font identities.
"""

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
        self.locus_report = self.root / "locus.report.json"
        self.locus_report.write_text(json.dumps({
            "schema": target.LOCUS_SVG_SCHEMA, "gene_symbol": "GENE",
            "panel_id": "GENE_panel", "gene_strand": "+",
            "isoform_evidence": {"assembly": "GRCh38", "chromosome": "chr1"},
            "regulatory_score_tracks": [{"track_id": "synthetic_matrix",
                "provider_kind": "jaspar_pwm", "provider_version": "MA9991.1",
                "source_ids": ["MA9991.1"]}],
        }))
        self.write_locus_receipt()
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
            "panel_resolution": {"matrices": [{"specification": {"source_id": "MA9991.1"}}]},
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
            "assembly_id": "GRCh38", "upstream_bp": 500, "downstream_bp": 200,
            "sequence_orientation": "transcript_5prime_to_3prime",
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
                **{name: digest(self.root / name)
                   for page in json.loads(self.index.read_text())["genes"][0]["pages"]
                   for name in page["files"]},
            },
        }))

    def write_locus_receipt(self):
        self.locus_receipt.write_text(json.dumps({
            "schema": target.LOCUS_RECEIPT_SCHEMA, "gene": "GENE",
            "inputs": {"locus_report": digest(self.locus_report)},
            "outputs": {self.locus_svg.name: digest(self.locus_svg)},
        }))

    def rebind(self):
        """Refresh file hashes only, so semantic inconsistencies reach validators."""
        self.checksums.write_text(f"{digest(self.fasta)}  {self.fasta.name}\n")
        manifest = json.loads(self.manifest.read_text())
        manifest["files"][0]["sha256"] = digest(self.fasta)
        manifest["sha256sums_sha256"] = digest(self.checksums)
        self.manifest.write_text(json.dumps(manifest))
        index = json.loads(self.index.read_text())
        index["report_sha256"] = digest(self.report)
        self.index.write_text(json.dumps(index))
        self.write_tss_receipt()
        self.write_locus_receipt()

    @staticmethod
    def fake_run(command, **kwargs):
        Path(command[2]).write_bytes(b"%PDF-1.4\nsynthetic")
        return subprocess.CompletedProcess(command, 0, stdout=json.dumps({
            "page_count": len(command) - 3,
            "pages": [{"font_identity_status": target.FONT_IDENTITY_STATUS,
                       "font_identities": []} for _ in command[3:]],
        }), stderr="")

    def assert_no_outputs(self):
        for path in (self.args().output_pdf, self.args().output_fasta, self.args().output_receipt):
            self.assertFalse(path.exists(), path)
            self.assertFalse(path.with_name(path.name + ".partial").exists(), path)

    def args(self):
        return argparse.Namespace(
            gene="GENE", locus_svg=self.locus_svg, locus_receipt=self.locus_receipt,
            locus_report=self.locus_report,
            tss_report=self.report, tss_index=self.index, tss_receipt=self.tss_receipt,
            tss_bundle_manifest=self.manifest,
            gentle_cli=self.cli, producer_revision="test-revision",
            output_pdf=self.root / "combined.pdf",
            output_fasta=self.root / "selected.fasta",
            output_receipt=self.root / "combined.receipt.json",
        )

    def test_compose_binds_context_and_only_selected_tss(self):
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
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
        self.assertEqual(receipt["locus_matrix_bindings"], [{
            "source_id": "MA9991.1", "locus_track_id": "synthetic_matrix",
        }])

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

    def test_rejects_changed_bases_even_with_rebound_file_hashes(self):
        self.fasta.write_text(self.fasta.read_text().replace("\n" + "A" * 701 + "\n",
                                                           "\n" + "T" * 701 + "\n"))
        self.rebind()
        with self.assertRaisesRegex(ValueError, "sequence digest does not match its bases"):
            target.compose(self.args())
        self.assert_no_outputs()

    def test_header_cannot_override_sequence_payload(self):
        self.fasta.write_text(self.fasta.read_text().replace("|gene_id=", "|sequence=AAAA|gene_id=", 1))
        self.rebind()
        with self.assertRaisesRegex(ValueError, "reserved field"):
            target.compose(self.args())

    def test_rejects_disagreeing_fasta_header_metadata(self):
        original = self.fasta.read_text()
        for old, new, diagnostic in [
            ("genomic_1based=500-1200", "genomic_1based=400-1100", "genomic span"),
            ("orientation=transcript_5prime_to_3prime", "orientation=genomic", "orientation"),
            ("window=minus500_plus200", "window=minus600_plus100", "window"),
            ("transcripts=ENST_TEST", "transcripts=ENST_OTHER", "transcript membership"),
            ("assembly=GRCh38", "assembly=GRCh37", "assembly"),
        ]:
            with self.subTest(field=old):
                self.fasta.write_text(original.replace(old, new, 1))
                self.rebind()
                with self.assertRaisesRegex(ValueError, diagnostic):
                    target.compose(self.args())
                self.assert_no_outputs()

    def test_rejects_disagreeing_report_or_manifest_span(self):
        for path, mutate in [
            (self.report, lambda value: value["windows"][0]["record"]["geometry"].update(start_1based=400)),
            (self.manifest, lambda value: value["files"][0]["records"][0].update(genomic_end_1based=1100)),
        ]:
            original = path.read_text()
            with self.subTest(path=path.name):
                value = json.loads(original)
                mutate(value)
                path.write_text(json.dumps(value))
                self.rebind()
                with self.assertRaisesRegex(ValueError, "genomic span"):
                    target.compose(self.args())
            path.write_text(original)
            self.rebind()

    def test_rejects_duplicate_manifest_promoters(self):
        value = json.loads(self.manifest.read_text())
        value["files"][0]["records"].append(value["files"][0]["records"][0])
        self.manifest.write_text(json.dumps(value))
        self.rebind()
        with self.assertRaisesRegex(ValueError, "unique and equal"):
            target.compose(self.args())

    def test_accepts_minus_strand_without_reversing_transcript_oriented_dna(self):
        fasta = self.fasta.read_text().replace("strand=+", "strand=-")
        report = json.loads(self.report.read_text())
        manifest = json.loads(self.manifest.read_text())
        for window, declared in zip(report["windows"], manifest["files"][0]["records"]):
            geometry = window["record"]["geometry"]
            tss = geometry["tss_1based"]
            geometry.update(strand="-", start_1based=tss - 200, end_1based=tss + 500)
            declared.update(strand="-", genomic_start_1based=tss - 200, genomic_end_1based=tss + 500)
            fasta = fasta.replace(f"genomic_1based={tss - 500}-{tss + 200}",
                                  f"genomic_1based={tss - 200}-{tss + 500}")
        locus = json.loads(self.locus_report.read_text())
        locus["gene_strand"] = "-"
        self.locus_report.write_text(json.dumps(locus))
        self.fasta.write_text(fasta)
        self.report.write_text(json.dumps(report))
        self.manifest.write_text(json.dumps(manifest))
        self.rebind()
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            receipt = target.compose(self.args())
        self.assertEqual(receipt["selected_tss_bindings"][0]["start_1based"], 800)
        self.assertEqual(target.parse_fasta(self.args().output_fasta)["selected"]["sequence"], "A" * 701)

    def test_rejects_locus_matrix_alias_collapse(self):
        value = json.loads(self.locus_report.read_text())
        value["regulatory_score_tracks"][0].update(
            source_ids=["MA9993.7"], provider_version="MA9993.7")
        self.locus_report.write_text(json.dumps(value))
        self.rebind()
        with self.assertRaisesRegex(ValueError, "exactly one resolved MA9991.1"):
            target.compose(self.args())
        self.assert_no_outputs()

    def test_rejects_unbound_locus_report(self):
        self.locus_report.write_text(self.locus_report.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "locus report hash mismatch"):
            target.compose(self.args())

    def test_rejects_cross_chromosome_or_strand_locus_join(self):
        original = self.locus_report.read_text()
        for field, changed in [("chromosome", "chr2"), ("gene_strand", "-")]:
            with self.subTest(field=field):
                value = json.loads(original)
                if field == "chromosome":
                    value["isoform_evidence"][field] = changed
                else:
                    value[field] = changed
                self.locus_report.write_text(json.dumps(value))
                self.rebind()
                with self.assertRaisesRegex(ValueError, "chromosome or strand mismatch"):
                    target.compose(self.args())

    def test_keeps_continuation_pages_but_exports_one_fasta_record(self):
        continued = self.root / "continued.svg"
        continued.write_text(self.selected_svg.read_text().replace('height="30"', 'height="31"'))
        index = json.loads(self.index.read_text())
        index["genes"][0]["pages"].insert(1, {"promoter_ids": ["selected"], "files": [continued.name]})
        self.index.write_text(json.dumps(index))
        self.rebind()
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            receipt = target.compose(self.args())
        self.assertEqual(receipt["page_count"], 3)
        self.assertEqual(receipt["selected_tss_count"], 1)
        self.assertEqual([Path(item["path"]).name for item in receipt["inputs"]["pages"]],
                         [self.locus_svg.name, self.selected_svg.name, continued.name])
        self.assertEqual(self.args().output_fasta.read_text().count(">"), 1)

    def test_rejects_duplicate_svg_not_continuation(self):
        index = json.loads(self.index.read_text())
        index["genes"][0]["pages"].append(index["genes"][0]["pages"][0])
        self.index.write_text(json.dumps(index))
        self.rebind()
        with self.assertRaisesRegex(ValueError, "repeats the same SVG page"):
            target.compose(self.args())

    def test_creates_receipt_parent_before_publication(self):
        args = self.args()
        args.output_receipt = self.root / "new_receipts" / "receipt.json"
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            target.compose(args)
        self.assertTrue(args.output_receipt.is_file())

    def test_stale_partial_receipt_does_not_publish_or_remove_foreign_file(self):
        partial = self.args().output_receipt.with_name(self.args().output_receipt.name + ".partial")
        partial.write_text("existing incomplete run")
        with mock.patch.object(target.subprocess, "run") as renderer:
            with self.assertRaisesRegex(ValueError, "stale partial"):
                target.compose(self.args())
            renderer.assert_not_called()
        self.assertFalse(self.args().output_pdf.exists())
        self.assertFalse(self.args().output_fasta.exists())
        self.assertEqual(partial.read_text(), "existing incomplete run")

    def test_receipt_write_failure_rolls_back_and_allows_retry(self):
        original_write = Path.write_text
        def fail_receipt(path, *args, **kwargs):
            if path.name == "combined.receipt.json.partial":
                raise OSError("synthetic receipt write failure")
            return original_write(path, *args, **kwargs)
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            with mock.patch.object(Path, "write_text", new=fail_receipt):
                with self.assertRaisesRegex(OSError, "receipt write failure"):
                    target.compose(self.args())
            self.assert_no_outputs()
            target.compose(self.args())

    def test_each_promotion_failure_rolls_back_and_allows_retry(self):
        original_link = target.os.link
        for failure in (1, 2, 3):
            calls = 0
            def fail_link(source, destination):
                nonlocal calls
                calls += 1
                if calls == failure:
                    raise OSError("synthetic promotion failure")
                return original_link(source, destination)
            with self.subTest(failure=failure), \
                 mock.patch.object(target.subprocess, "run", side_effect=self.fake_run), \
                 mock.patch.object(target.os, "link", side_effect=fail_link):
                with self.assertRaisesRegex(OSError, "promotion failure"):
                    target.compose(self.args())
                self.assert_no_outputs()
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            target.compose(self.args())

    def test_racing_output_is_not_overwritten_or_deleted(self):
        original_link = target.os.link
        def race(source, destination):
            destination.write_bytes(b"another producer's output")
            return original_link(source, destination)
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run), \
             mock.patch.object(target.os, "link", side_effect=race):
            with self.assertRaises(FileExistsError):
                target.compose(self.args())
        self.assertEqual(self.args().output_pdf.read_bytes(), b"another producer's output")
        self.assertFalse(self.args().output_fasta.exists())
        self.assertFalse(self.args().output_receipt.exists())

    def test_rejects_renderer_without_used_font_audit(self):
        def missing_audit(command, **kwargs):
            result = self.fake_run(command, **kwargs)
            value = json.loads(result.stdout)
            value["pages"][0] = {"font_face_count": 205}
            result.stdout = json.dumps(value)
            return result
        with mock.patch.object(target.subprocess, "run", side_effect=missing_audit):
            with self.assertRaisesRegex(ValueError, "lacks used-font audit"):
                target.compose(self.args())
        self.assert_no_outputs()

    def test_detail_context_must_bind_the_same_locus_as_the_overview(self):
        report = json.loads(self.report.read_text())
        report["windows"][0]["detail_context"] = {
            "schema": "gentle.tss_detail_context.v1",
            "locus_report_sha256": "f" * 64,
        }
        self.report.write_text(json.dumps(report))
        self.rebind()
        with self.assertRaisesRegex(ValueError, "different locus reports"):
            target.compose(self.args())
        self.assert_no_outputs()
        report["windows"][0]["detail_context"]["locus_report_sha256"] = digest(self.locus_report)
        self.report.write_text(json.dumps(report))
        self.rebind()
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            target.compose(self.args())


if __name__ == "__main__":
    unittest.main()
