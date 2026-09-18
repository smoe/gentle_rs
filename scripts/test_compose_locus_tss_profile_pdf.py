#!/usr/bin/env python3
"""Compositor regressions using hand-crafted temporary JSON/SVG/FASTA/GenBank/EMBL.

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
                "selection_window": {"upstream_bp": 2000, "downstream_bp": 200,
                                     "length_bp": 2201, "sequence_sha256": "a" * 64},
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
                **{name: digest(self.root / name)
                   for kind in ("genbank", "embl")
                   for name in json.loads(self.index.read_text())["genes"][0].get(kind, {}).values()},
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

    def test_ordered_svg_bundle_preserves_hover_and_exact_pdf_input_bytes(self):
        import zipfile
        args = self.args()
        args.output_svg_directory = self.root / "svg-pages"
        args.output_svg_zip = self.root / "svg-pages.zip"
        observed = []

        def render(command, **kwargs):
            observed.extend(Path(path).read_bytes() for path in command[3:])
            return self.fake_run(command, **kwargs)

        with mock.patch.object(target.subprocess, "run", side_effect=render):
            receipt = target.compose(args)
        bundle = receipt["output"]["svg_bundle"]
        self.assertEqual(len(observed), receipt["page_count"])
        for index, page in enumerate(bundle["pages"]):
            path = args.output_svg_directory / page["file"]
            self.assertEqual(path.read_bytes(), observed[index])
            self.assertEqual(digest(path), page["sha256"])
            self.assertEqual(page["role"], receipt["page_order"][index])
        for name, expected in bundle["files"].items():
            self.assertEqual(digest(args.output_svg_directory / name), expected)
        self.assertEqual(digest(args.output_svg_zip), bundle["zip"]["sha256"])
        with zipfile.ZipFile(args.output_svg_zip) as archive:
            self.assertEqual(archive.namelist(), sorted(bundle["files"]))
            for info in archive.infolist():
                self.assertEqual(info.date_time, (1980, 1, 1, 0, 0, 0))
                self.assertEqual(archive.read(info.filename), (args.output_svg_directory / info.filename).read_bytes())
        html = (args.output_svg_directory / "index.html").read_text()
        self.assertIn('type="image/svg+xml"', html)
        self.assertIn('selected_tss.fasta', html)

    def test_svg_bundle_failure_rolls_back_without_a_success_receipt(self):
        args = self.args()
        args.output_svg_directory = self.root / "svg-pages"
        args.output_svg_zip = self.root / "svg-pages.zip"
        with mock.patch.object(target.subprocess, "run", side_effect=RuntimeError("renderer failed")):
            with self.assertRaises(RuntimeError):
                target.compose(args)
        self.assert_no_outputs(args)
        self.assertFalse(args.output_svg_directory.exists())
        self.assertFalse(args.output_svg_zip.exists())

    def test_svg_bundle_rejects_page_tamper_and_keeps_title_bytes(self):
        args = self.args()
        # Explicit synthetic hover text is an annotation, not an instruction.
        self.locus_svg.write_text(self.locus_svg.read_text().replace('</svg>', '<title>source E.1 and R.2</title></svg>'))
        self.rebind()
        args.output_svg_directory = self.root / "svg-pages"
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            receipt = target.compose(args)
        first = receipt["output"]["svg_bundle"]["pages"][0]["file"]
        self.assertIn(b'<title>source E.1 and R.2</title>', (args.output_svg_directory / first).read_bytes())
        self.locus_svg.write_bytes(self.locus_svg.read_bytes() + b'\n<!-- changed -->')
        with self.assertRaises(ValueError):
            target.svg_sequence_bundle(receipt, [self.locus_svg], "A", self.root / "new", None)

    def source_coherent_fixture(self):
        # Synthetic binding token only: the Rust tests verify geometry; the
        # compositor must compare canonical content, not reimplement grouping.
        presentation = {"schema": "gentle.transcript_structure_presentation.v1", "content_sha256": "a" * 64}
        locus = json.loads(self.locus_report.read_text())
        locus["transcript_presentation"] = presentation
        self.locus_report.write_text(json.dumps(locus))
        report = json.loads(self.report.read_text())
        for window in report["windows"]:
            if window.get("selected"):
                window["detail_context"] = {"schema": "gentle.tss_detail_context.v1",
                    "locus_report_sha256": digest(self.locus_report), "transcript_presentation": presentation}
        self.report.write_text(json.dumps(report))
        for path in self.root.glob("*.svg"):
            root = target.ET.fromstring(path.read_bytes())
            target.ET.SubElement(root, "g", {"data-role": "source-coherent-transcripts", "data-content-sha256": "a" * 64})
            path.write_bytes(target.ET.tostring(root))
        self.rebind()

    def test_composite_requires_the_same_transcript_presentation_and_svg_binding(self):
        self.source_coherent_fixture()
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            receipt = target.compose(self.args())
        self.assertEqual(receipt["inputs"]["transcript_presentation_sha256"], "a" * 64)

    def test_rebound_different_detail_presentation_is_rejected_before_pdf(self):
        self.source_coherent_fixture()
        report = json.loads(self.report.read_text())
        next(w for w in report["windows"] if w.get("selected"))["detail_context"]["transcript_presentation"]["content_sha256"] = "b" * 64
        self.report.write_text(json.dumps(report))
        self.rebind()
        with mock.patch.object(target.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "transcript presentations differ"):
                target.compose(self.args())
            run.assert_not_called()
        self.assert_no_outputs()

    def test_rebound_svg_with_missing_shared_structure_layer_is_rejected(self):
        self.source_coherent_fixture()
        self.locus_svg.write_text(self.locus_svg.read_text().replace('data-role="source-coherent-transcripts"', 'data-role="omitted"'))
        self.rebind()
        with mock.patch.object(target.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "SVG page does not bind"):
                target.compose(self.args())
            run.assert_not_called()
        self.assert_no_outputs()

    def assert_no_outputs(self, args=None):
        args = args or self.args()
        for path in (args.output_pdf, args.output_fasta, args.output_receipt,
                     getattr(args, "output_genbank", None), getattr(args, "output_embl", None)):
            if path is None:
                continue
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
        self.assertEqual(receipt["selected_tss_bindings"][0]["historical_selection_window"]["length_bp"], 2201)
        self.assertEqual(receipt["selected_tss_bindings"][0]["selection_window_join_policy"],
                         "provenance_only_not_display_geometry")
        self.assertEqual(receipt["selected_tss_bindings"][0]["covering_tss_bands"],
                         ["GENE_tss_stretch_1"])
        self.assertTrue(self.args().output_pdf.is_file())
        self.assertTrue(self.args().output_fasta.is_file())
        self.assertTrue(self.args().output_receipt.is_file())
        fasta = self.args().output_fasta.read_text()
        self.assertIn("promoter_id=selected", fasta)
        self.assertNotIn("promoter_id=other", fasta)
        self.assertEqual(receipt["producer"]["pdf_representation"], "raster")
        self.assertEqual(receipt["horizontal_alignment"], {
            "page_width": 1400.0,
            "plot_left": 255.0,
            "plot_right": 1050.0,
            "policy": "context and detail pages share the canonical locus plot frame",
        })

    def test_vector_pdf_representation_uses_typed_renderer_and_binds_contract(self):
        args = self.args()
        args.pdf_representation = "vector"

        def render(command, **kwargs):
            self.assertEqual(command[1], "svg-vector-pdf-set")
            Path(command[2]).write_bytes(b"%PDF-1.7\nsynthetic vector")
            return subprocess.CompletedProcess(command, 0, stdout=json.dumps({
                "page_count": len(command) - 3,
                "pdf_representation": "static multipage vector PDF with embedded selectable text",
                "embedded_text": True,
                "svg_interactivity_preserved": False,
                "svg_uri_links_preserved": False,
                "pages": [{"font_identity_status": target.FONT_IDENTITY_STATUS,
                           "font_identities": []} for _ in command[3:]],
            }), stderr="")

        with mock.patch.object(target.subprocess, "run", side_effect=render):
            receipt = target.compose(args)
        self.assertEqual(receipt["producer"]["pdf_representation"], "vector")
        self.assertTrue(args.output_pdf.read_bytes().startswith(b"%PDF-1.7"))
        self.assertEqual(receipt["locus_matrix_bindings"], [{
            "source_id": "MA9991.1", "locus_track_id": "synthetic_matrix",
        }])

    def genbank_fixture(self):
        index = json.loads(self.index.read_text())
        files = {}
        for promoter, bases in self.sequences.items():
            file = self.root / f"{promoter}.gb"
            file.write_text(f"LOCUS       {promoter} 701 bp DNA linear\nFEATURES             Location/Qualifiers\n     misc_feature    501\n                     /label=\"Synthetic TSS\"\nORIGIN\n        1 {bases.lower()}\n//\n")
            files[promoter] = file.name
        index["genes"][0]["genbank"] = files
        self.index.write_text(json.dumps(index))
        self.write_tss_receipt()

    def embl_fixture(self):
        index = json.loads(self.index.read_text())
        files = {}
        for promoter, bases in self.sequences.items():
            file = self.root / f"{promoter}.embl"
            sequence_lines = []
            for offset in range(0, len(bases), 60):
                groups = " ".join(bases[i:i + 10].lower()
                                  for i in range(offset, min(offset + 60, len(bases)), 10))
                sequence_lines.append(f"     {groups:<65} {min(offset + 60, len(bases))}\n")
            file.write_text(
                f"ID   {promoter}; SV 1; linear; DNA; STD; SYN; 701 BP.\n"
                "XX\nFH   Key             Location/Qualifiers\nFH\n"
                "FT   misc_feature    501\nFT                   /label=\"Synthetic TSS\"\nXX\n"
                f"SQ   Sequence 701 BP; {bases.count('A')} A; {bases.count('C')} C; "
                f"{bases.count('G')} G; {bases.count('T')} T; {bases.count('N')} other;\n"
                + "".join(sequence_lines) + "//\n"
            )
            files[promoter] = file.name
        index["genes"][0]["embl"] = files
        self.index.write_text(json.dumps(index))
        self.write_tss_receipt()

    def test_genbank_selected_records_are_copied_and_receipt_bound(self):
        self.genbank_fixture()
        args = self.args()
        args.output_genbank = self.root / "annotated.gb"
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            receipt = target.compose(args)
        self.assertEqual(args.output_genbank.read_bytes(), (self.root / "selected.gb").read_bytes())
        self.assertEqual(receipt["output"]["selected_tss_genbank_sha256"], digest(args.output_genbank))
        self.assertEqual(receipt["inputs"]["selected_tss_genbank"][0]["promoter_id"], "selected")

    def test_genbank_missing_tampered_or_wrong_bases_fail_before_publication(self):
        args = self.args()
        args.output_genbank = self.root / "annotated.gb"
        with self.assertRaisesRegex(ValueError, "lacks indexed GenBank"):
            target.compose(args)
        self.genbank_fixture()
        path = self.root / "selected.gb"
        path.write_text(path.read_text().replace("a" * 701, "t" * 701))
        with self.assertRaisesRegex(ValueError, "output hash mismatch"):
            target.compose(args)
        receipt = json.loads(self.tss_receipt.read_text())
        receipt["outputs"][path.name] = digest(path)
        self.tss_receipt.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "GenBank bases differ"):
            target.compose(args)
        self.assert_no_outputs()
        self.assertFalse(args.output_genbank.exists())

    def test_genbank_publication_failure_rolls_back_whole_new_bundle(self):
        self.genbank_fixture()
        args = self.args()
        args.output_genbank = self.root / "annotated.gb"
        link = target.os.link
        def fail_receipt(src, dst):
            if dst == args.output_receipt.resolve():
                raise OSError("synthetic receipt failure")
            return link(src, dst)
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run), mock.patch.object(target.os, "link", side_effect=fail_receipt):
            with self.assertRaisesRegex(OSError, "synthetic receipt failure"):
                target.compose(args)
        self.assert_no_outputs()
        self.assertFalse(args.output_genbank.exists())

    def test_embl_selected_record_is_copied_and_receipt_bound_without_genbank(self):
        self.embl_fixture()
        path = self.root / "selected.embl"
        # Preserve source bytes, including CRLF and annotations; ignore unselected files.
        path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
        self.write_tss_receipt()
        (self.root / "other.embl").unlink()
        args = self.args()
        args.output_embl = self.root / "annotated.embl"
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            receipt = target.compose(args)
        self.assertEqual(args.output_embl.read_bytes(), path.read_bytes())
        self.assertEqual(receipt["inputs"]["selected_tss_embl"], [{
            "promoter_id": "selected", "path": str(path.resolve()), "sha256": digest(path),
        }])
        self.assertEqual(receipt["output"]["selected_tss_embl"], str(args.output_embl.resolve()))
        self.assertEqual(receipt["output"]["selected_tss_embl_sha256"], digest(args.output_embl))
        self.assertEqual(receipt["output"]["selected_tss_embl_bytes"], args.output_embl.stat().st_size)
        self.assertNotIn("selected_tss_genbank", receipt["inputs"])
        self.assertNotIn("selected_tss_genbank", receipt["output"])

    def test_embl_missing_tampered_or_wrong_bases_fail_before_publication(self):
        args = self.args()
        args.output_embl = self.root / "annotated.embl"
        with mock.patch.object(target.subprocess, "run") as renderer:
            with self.assertRaisesRegex(ValueError, "lacks indexed EMBL"):
                target.compose(args)
            self.embl_fixture()
            path = self.root / "selected.embl"
            original = path.read_bytes()
            path.unlink()
            with self.assertRaisesRegex(ValueError, "EMBL must be a direct regular file"):
                target.compose(args)
            path.write_bytes(original)
            receipt = json.loads(self.tss_receipt.read_text())
            del receipt["outputs"][path.name]
            self.tss_receipt.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, "receipt does not bind output selected.embl"):
                target.compose(args)
            self.write_tss_receipt()
            path.write_bytes(original.replace(b"a", b"t"))
            with self.assertRaisesRegex(ValueError, "output hash mismatch"):
                target.compose(args)
            self.write_tss_receipt()
            with self.assertRaisesRegex(ValueError, "EMBL bases differ"):
                target.compose(args)
            renderer.assert_not_called()
        self.assert_no_outputs(args)

    def test_embl_requires_exactly_one_complete_record_and_sequence_below_sq(self):
        self.embl_fixture()
        path = self.root / "selected.embl"
        original = path.read_text()
        sq = original.index("SQ   ")
        sequence_start = original.index("\n", sq) + 1
        invalid = {
            "missing_id": original.replace("ID   ", "XX   ", 1),
            "duplicate_id": original.replace("XX\n", "ID   duplicate\n", 1),
            "missing_sq": original.replace("SQ   ", "XX   "),
            "duplicate_sq": original.replace("SQ   ", "SQ   duplicate\nSQ   "),
            "two_records": original + original,
            "missing_terminator": original.removesuffix("//\n"),
            "early_terminator": original.replace("XX\n", "//\n", 1),
            "trailing_record_text": original + "XX\n",
            "no_bases": original[:sequence_start] + "//\n",
            "invalid_base": original[:sequence_start] + "     Z 1\n//\n",
        }
        args = self.args()
        args.output_embl = self.root / "annotated.embl"
        for label, value in invalid.items():
            with self.subTest(label=label), mock.patch.object(target.subprocess, "run") as renderer:
                path.write_text(value)
                self.write_tss_receipt()
                with self.assertRaisesRegex(ValueError, "complete engine-exported EMBL record|EMBL bases differ"):
                    target.compose(args)
                renderer.assert_not_called()
                self.assert_no_outputs(args)

    def test_embl_hash_covers_the_exact_copied_bytes_not_an_earlier_read(self):
        self.embl_fixture()
        args = self.args()
        args.output_embl = self.root / "annotated.embl"
        read = Path.read_bytes
        def changed_annotations(path):
            raw = read(path)
            if path == (self.root / "selected.embl").resolve():
                return raw.replace(b"Synthetic TSS", b"Changed TSS label")
            return raw
        with mock.patch.object(target.subprocess, "run") as renderer, \
             mock.patch.object(Path, "read_bytes", new=changed_annotations):
            with self.assertRaisesRegex(ValueError, "output hash mismatch"):
                target.compose(args)
            renderer.assert_not_called()
        self.assert_no_outputs(args)

    def test_embl_source_paths_must_be_direct_regular_files(self):
        self.embl_fixture()
        index = json.loads(self.index.read_text())
        args = self.args()
        args.output_embl = self.root / "annotated.embl"
        for name in ("../selected.embl", str(self.root / "selected.embl"), "selected.gb"):
            with self.subTest(name=name):
                index["genes"][0]["embl"]["selected"] = name
                self.index.write_text(json.dumps(index))
                receipt = json.loads(self.tss_receipt.read_text())
                receipt["outputs"][self.index.name] = digest(self.index)
                self.tss_receipt.write_text(json.dumps(receipt))
                with self.assertRaisesRegex(ValueError, "lacks indexed EMBL"):
                    target.compose(args)
                self.assert_no_outputs(args)
        self.embl_fixture()
        path = self.root / "selected.embl"
        path.unlink()
        path.symlink_to(self.root / "other.embl")
        with self.assertRaisesRegex(ValueError, "EMBL must be a direct regular file"):
            target.compose(args)
        self.assert_no_outputs(args)

    def test_both_annotated_formats_follow_selected_report_order_and_publish_receipt_last(self):
        report = json.loads(self.report.read_text())
        report["windows"] = [self.window("other", 1050, True), self.window("selected", 1000, True)]
        self.report.write_text(json.dumps(report))
        self.rebind()
        self.genbank_fixture()
        self.embl_fixture()
        args = self.args()
        args.output_genbank = self.root / "annotated.gb"
        args.output_embl = self.root / "annotated.embl"
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run), \
             mock.patch.object(target.os, "link", wraps=target.os.link) as publish:
            receipt = target.compose(args)
        self.assertEqual([call.args[1] for call in publish.call_args_list], [path.resolve() for path in (
            args.output_pdf, args.output_fasta, args.output_genbank, args.output_embl, args.output_receipt,
        )])
        for kind, suffix in (("genbank", ".gb"), ("embl", ".embl")):
            output = getattr(args, f"output_{kind}")
            self.assertEqual(output.read_bytes(), b"".join(
                (self.root / (promoter + suffix)).read_bytes() for promoter in ("other", "selected")))
            self.assertEqual([row["promoter_id"] for row in receipt["inputs"][f"selected_tss_{kind}"]],
                             ["other", "selected"])
            self.assertEqual(receipt["output"][f"selected_tss_{kind}_sha256"], digest(output))
        self.assertEqual(json.loads(args.output_receipt.read_text()), receipt)

    def test_embl_and_both_formats_roll_back_at_each_publication_step(self):
        self.genbank_fixture()
        self.embl_fixture()
        link = target.os.link
        for both in (False, True):
            args = self.args()
            args.output_embl = self.root / "annotated.embl"
            if both:
                args.output_genbank = self.root / "annotated.gb"
            for failure in range(1, 6 if both else 5):
                calls = 0
                def fail_link(src, dst):
                    nonlocal calls
                    calls += 1
                    if calls == failure:
                        raise OSError("synthetic annotated promotion failure")
                    return link(src, dst)
                with self.subTest(both=both, failure=failure), \
                     mock.patch.object(target.subprocess, "run", side_effect=self.fake_run), \
                     mock.patch.object(target.os, "link", side_effect=fail_link):
                    with self.assertRaisesRegex(OSError, "annotated promotion failure"):
                        target.compose(args)
                    self.assert_no_outputs(args)
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run):
            target.compose(args)

    def test_embl_write_failure_rolls_back_all_staged_outputs(self):
        self.genbank_fixture()
        self.embl_fixture()
        args = self.args()
        args.output_genbank = self.root / "annotated.gb"
        args.output_embl = self.root / "annotated.embl"
        write = Path.write_bytes
        def fail_write(path, data):
            if path.name == "annotated.embl.partial":
                raise OSError("synthetic EMBL write failure")
            return write(path, data)
        with mock.patch.object(target.subprocess, "run") as renderer, \
             mock.patch.object(Path, "write_bytes", new=fail_write):
            with self.assertRaisesRegex(OSError, "EMBL write failure"):
                target.compose(args)
            renderer.assert_not_called()
        self.assert_no_outputs(args)

    def test_embl_output_collisions_and_stale_partials_are_not_overwritten(self):
        self.genbank_fixture()
        self.embl_fixture()
        args = self.args()
        args.output_genbank = self.root / "annotated.gb"
        for path in (args.output_pdf, args.output_fasta, args.output_receipt, args.output_genbank):
            with self.subTest(path=path), mock.patch.object(target.subprocess, "run") as renderer:
                args.output_embl = path
                with self.assertRaisesRegex(ValueError, "output paths must be distinct"):
                    target.compose(args)
                renderer.assert_not_called()
                self.assert_no_outputs(args)
        args.output_embl = self.root / "annotated.embl"
        for path in (args.output_embl, args.output_embl.with_suffix(".embl.partial")):
            with self.subTest(path=path), mock.patch.object(target.subprocess, "run") as renderer:
                path.write_text("existing unrelated output")
                with self.assertRaisesRegex(ValueError, "output or stale partial"):
                    target.compose(args)
                renderer.assert_not_called()
                self.assertEqual(path.read_text(), "existing unrelated output")
                path.unlink()
                self.assert_no_outputs(args)

    def test_cli_accepts_either_or_both_annotated_outputs(self):
        base = ["compose_locus_tss_profile_pdf.py"]
        for name, value in vars(self.args()).items():
            base.extend(["--" + name.replace("_", "-"), str(value)])
        for kinds in (("embl",), ("genbank",), ("genbank", "embl")):
            with self.subTest(kinds=kinds):
                argv = base + [arg for kind in kinds for arg in (f"--output-{kind}", f"output.{kind}")]
                with mock.patch.object(sys, "argv", argv):
                    parsed = target.parse_args()
                for kind in ("genbank", "embl"):
                    self.assertEqual(getattr(parsed, f"output_{kind}"),
                                     Path(f"output.{kind}") if kind in kinds else None)

    def test_racing_embl_output_is_preserved_while_other_outputs_roll_back(self):
        self.genbank_fixture()
        self.embl_fixture()
        args = self.args()
        args.output_genbank = self.root / "annotated.gb"
        args.output_embl = self.root / "annotated.embl"
        link = target.os.link
        def race(src, dst):
            if dst == args.output_embl.resolve():
                dst.write_bytes(b"another producer's EMBL output")
            return link(src, dst)
        with mock.patch.object(target.subprocess, "run", side_effect=self.fake_run), \
             mock.patch.object(target.os, "link", side_effect=race):
            with self.assertRaises(FileExistsError):
                target.compose(args)
        self.assertEqual(args.output_embl.read_bytes(), b"another producer's EMBL output")
        args.output_embl.unlink()
        self.assert_no_outputs(args)

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
