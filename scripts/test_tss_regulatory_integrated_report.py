#!/usr/bin/env python3
"""Focused deterministic tests for the TSS-local integrated report helpers."""

from __future__ import annotations

import importlib.util
from contextlib import redirect_stdout
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import tss_regulatory_report_binding as binding
from test_tp73_cutrun_promoter_comparison import reference_fixture, write_json, write_tsv


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PREPARE = load_module(
    "prepare_tss_regulatory_similarity_candidates",
    "prepare_tss_regulatory_similarity_candidates.py",
)
RENDER = load_module(
    "render_integrated_tss_regulatory_report",
    "render_integrated_tss_regulatory_report.py",
)
APPEND = load_module(
    "append_tss_similarity_to_locus_report",
    "append_tss_similarity_to_locus_report.py",
)
BUNDLE = ROOT.parent / "docs" / "examples" / "regulatory_region_comparison" / \
    "tp73_tss_local_integrated"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_fixture(root, strand="+"):
    """Hand-crafted chr1/TOY inputs; recreated in tempdirs, never scientific evidence.

    Used by preparation and rendering rejection tests. The two-promoter reference
    is shared with the existing comparison tests; no downloads or BLAST are used.
    """
    reference, windows, mappings = reference_fixture(root / "reference")
    catalog = root / "reference/genomes.json"
    write_json(catalog, {reference["genome_id"]: {
        "ensembl_template": {"file_stem": "Homo_sapiens.GRCh38"}}})
    reference.update(catalog_path=str(catalog), catalog_sha256=f"sha256:{sha256(catalog)}")
    write_json(root / "reference/receipt.json", reference)
    window = next(row for row in windows if row["strand"] == strand)
    mapping = next(row for row in mappings if row["promoter_id"] == window["promoter_id"])
    extraction = dict(genome_id=reference["genome_id"], chromosome="1", strand=strand,
                      start_1based=window["start_0based"] + 1,
                      end_1based=window["end_0based_exclusive"], tss_1based=window["tss_1based"],
                      promoter_upstream_bp=2000, promoter_downstream_bp=200,
                      transcript_ids=[mapping["transcript_id"]], gene_id="G1", gene_name="TOY")
    selected = {
        "schema": PREPARE.SCHEMA,
        "source_bindings": {"promoterome_receipt_sha256": sha256(root / "reference/receipt.json")},
        "regions": [{"region_id": "selected", "promoterome_id": window["promoter_id"],
                     "gene_query": "TOY", "transcript_ids": extraction["transcript_ids"][:],
                     "genome_extraction": extraction, "sequence_length_bp": 2201,
                     "sequence_orientation": "biological_5prime_to_3prime",
                     "sequence_sha256": hashlib.sha256(("A" * 2201).encode()).hexdigest()}],
    }
    tss = window["tss_1based"]
    report = {
        "schema": binding.LOCUS_SCHEMA, "seq_id": "toy", "gene_symbol": "TOY",
        "panel_id": "toy_panel", "gene_strand": strand,
        "axis_left_genomic_1based": 1 if strand == "+" else 10000,
        "axis_right_genomic_1based": 10000 if strand == "+" else 1,
        "sequence_binding": {"genome_anchor": {"genome_id": "GRCh38",
                             "chromosome": "1", "start_1based": 1, "end_1based": 10000}},
        "isoform_evidence": {"chromosome": "1"},
        "ensembl_regulation": {
            "source_binding": {"content_identity_verified": True, "truncated": False},
            "rows": [{"feature_id": "F1", "feature_type": "promoter",
                      "core_genomic_start_1based": tss - 100, "core_genomic_end_1based": tss + 100,
                      "assembly_name": "GRCh38", "assembly_accession": "synthetic",
                      "canonical_feature_url": "https://example.invalid/synthetic/F1",
                      "source_id": "synthetic", "annotation_release": "synthetic"}],
        },
    }
    write_json(root / "selected.json", selected)
    write_json(root / "report.json", report)
    svg = ('<svg data-gentle-schema="' + binding.LOCUS_SCHEMA + '" '
           'data-gentle-panel-id="toy_panel" height="200" viewBox="0 0 1400 200" width="1400">'
           '<rect fill="#ffffff" height="200" width="1400" x="0" y="0"/>'
           '<text x="34" y="20">existing evidence</text>'
           '<text x="34" y="80">Transcript models and annotation-derived metrics</text>'
           '<line data-gentle-transcript="T1" x1="255" x2="1050" y1="110" y2="110"/>'
           '<text data-gentle-overlay-non-claims="true" x="34" y="170">'
           'Reporter interpretation boundaries</text>'
           '<text x="34" y="190">Evidence provenance</text></svg>')
    (root / "base.svg").write_text(svg)
    return selected, report, svg


def run_preparation(root, assembly_id="GRCh38"):
    argv = ["prepare", "--selected-tss", str(root / "selected.json"),
            "--locus-report", str(root / "report.json"), "--locus-svg", f"TOY={root / 'base.svg'}",
            "--promoterome", str(root / "reference"), "--source-revision", "synthetic-revision",
            "--assembly-id", assembly_id,
            "--output", str(root / "output")]
    with patch.object(sys, "argv", argv), patch.object(
            PREPARE.subprocess, "check_output", return_value="synthetic-revision\n"), redirect_stdout(io.StringIO()):
        PREPARE.main()
    return json.loads((root / "output/candidate_regions.json").read_text())


def plot_fixture(strand="+", match_counts=(2,)):
    """Synthetic display-only rows with known query offsets and no biological claim."""
    candidates = {
        "stretches": [{"gene": "TOY", "stretch_id": "TOY_tss_stretch_1",
                       "start_1based": 100, "end_1based": 800,
                       "tss_windows": [{"strand": strand, "chromosome": "1", "tss_1based": 600}]}],
        "regions": [],
    }
    matches, hsps, summaries = {}, {}, {}
    for index, count in enumerate(match_counts):
        query = f"query_{index}"
        length = 20 if count < 0 else 200
        candidates["regions"].append({
            "region_id": query, "sequence_length_bp": length,
            "source_region": {"stretch_id": "TOY_tss_stretch_1", "region_id": query,
                              "feature_type": "promoter",
                              "interval": {"start_0based": 299, "end_0based_exclusive": 299 + length}},
        })
        summaries[query] = {"other_promoters": {"distinct_genes": max(0, count)},
                            "other_promoters_by_min_query_coverage": {
                                key: {"distinct_genes": 0} for key in ("0.25", "0.50")},
                            "counts_are_lower_bounds": False}
        matches[query] = []
        for rank in range(count):
            target = f"target_{rank}"
            matches[query].append({"promoter_id": target, "gene_ids": f"g{rank}",
                                   "gene_names": f"gene{rank}", "aligned_query_fraction": "0.1",
                                   "best_bitscore": "80", "query_intervals_0based_half_open": "0-20"})
            hsps[query, target] = [{"qstart": 1, "qend": 20, "sstart": 11, "send": 30,
                                   "bitscore": 80, "pident": 90}]
    thresholds = {"minimum_bp": 40, "minimum_identity": 80.0, "maximum_evalue": 1e-5}
    return candidates, matches, hsps, summaries, thresholds


class SourceBindingTests(unittest.TestCase):
    def test_preparation_binds_sources_and_extracts_both_strands(self):
        for strand, expected_base in [("+", "A"), ("-", "T")]:
            with self.subTest(strand=strand), TemporaryDirectory() as tmp:
                root = Path(tmp)
                source_fixture(root, strand)
                candidates = run_preparation(root)
                report, _ = binding.load_bound_locus_report(candidates, root / "report.json", "TOY")
                binding.load_bound_locus_svg(candidates, root / "base.svg", report)
                self.assertEqual(len(candidates["regions"]), 1)
                sequence = "".join((root / "output/candidate_regions.fa").read_text().splitlines()[1:])
                self.assertEqual(sequence, expected_base * 201)

    def test_selected_identity_mismatches_fail_before_output(self):
        for defect in ("chromosome", "strand", "tss", "span", "transcript", "empty_transcripts",
                       "duplicate_transcripts", "gene", "digest", "orientation"):
            with self.subTest(defect=defect), TemporaryDirectory() as tmp:
                root = Path(tmp)
                selected, _, _ = source_fixture(root)
                row = selected["regions"][0]
                extraction = row["genome_extraction"]
                if defect == "chromosome":
                    extraction["chromosome"] = "2"
                elif defect == "strand":
                    extraction["strand"] = "-"
                elif defect == "tss":
                    extraction["tss_1based"] += 1
                elif defect == "span":
                    extraction["start_1based"] += 1
                elif defect in {"transcript", "empty_transcripts", "duplicate_transcripts"}:
                    ids = {"transcript": ["T2"], "empty_transcripts": [],
                           "duplicate_transcripts": ["T1", "T1"]}[defect]
                    row["transcript_ids"] = extraction["transcript_ids"] = ids
                elif defect == "gene":
                    extraction["gene_id"] = "G2"
                elif defect == "digest":
                    row["sequence_sha256"] = "0" * 64
                else:
                    row["sequence_orientation"] = "assembly_reference_forward"
                write_json(root / "selected.json", selected)
                with self.assertRaises(RuntimeError):
                    run_preparation(root)
                self.assertFalse((root / "output").exists())

    def test_duplicate_reference_inventory_fails_even_when_rehashed(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            selected, _, _ = source_fixture(root)
            path = root / "reference/promoter_transcripts.tsv"
            rows = PREPARE.load_tsv(path)
            rows.append(rows[0])
            write_tsv(path, rows)
            receipt_path = root / "reference/receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["artifacts"][path.name] = f"sha256:{sha256(path)}"
            receipt["included_transcript_count"] += 1
            write_json(receipt_path, receipt)
            selected["source_bindings"]["promoterome_receipt_sha256"] = sha256(receipt_path)
            write_json(root / "selected.json", selected)
            with self.assertRaisesRegex(RuntimeError, "inventory"):
                run_preparation(root)
            self.assertFalse((root / "output").exists())

    def test_locus_reference_mismatches_fail_before_output(self):
        for defect in ("genome", "chromosome", "strand", "range", "missing_binding", "svg_panel"):
            with self.subTest(defect=defect), TemporaryDirectory() as tmp:
                root = Path(tmp)
                _, report, svg = source_fixture(root)
                anchor = report["sequence_binding"]["genome_anchor"]
                if defect == "genome":
                    anchor["genome_id"] = "another assembly"
                elif defect == "chromosome":
                    anchor["chromosome"] = "2"
                elif defect == "strand":
                    report["gene_strand"] = "-"
                elif defect == "range":
                    anchor["end_1based"] = 100
                elif defect == "missing_binding":
                    del report["sequence_binding"]
                else:
                    (root / "base.svg").write_text(svg.replace("toy_panel", "different_panel"))
                write_json(root / "report.json", report)
                with self.assertRaises(RuntimeError):
                    run_preparation(root)
                self.assertFalse((root / "output").exists())

    def test_assembly_identifier_is_exact_and_independent_of_catalog_words(self):
        rejected = ("Human", "Ensembl", "116", "GRCh3", "GRCh37",
                    "Human GRCh38 Ensembl 115")
        for assembly_id in rejected:
            with self.subTest(assembly_id=assembly_id), TemporaryDirectory() as tmp:
                root = Path(tmp)
                source_fixture(root)
                argv = ["prepare", "--selected-tss", str(root / "selected.json"),
                        "--locus-report", str(root / "report.json"),
                        "--locus-svg", f"TOY={root / 'base.svg'}",
                        "--promoterome", str(root / "reference"),
                        "--assembly-id", assembly_id,
                        "--source-revision", "synthetic-revision",
                        "--output", str(root / "output")]
                with patch.object(sys, "argv", argv), patch.object(
                        PREPARE.subprocess, "check_output",
                        return_value="synthetic-revision\n"), redirect_stdout(io.StringIO()):
                    with self.assertRaisesRegex(RuntimeError, "declared assembly"):
                        PREPARE.main()
                self.assertFalse((root / "output").exists())

    def test_reference_assembly_is_bound_to_catalog_content_not_matching_report_labels(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, report, _ = source_fixture(root)
            reference = json.loads((root / "reference/receipt.json").read_text())
            # A report and CLI can agree with one another while using the wrong reference.
            report["sequence_binding"]["genome_anchor"]["genome_id"] = "GRCh37"
            for row in report["ensembl_regulation"]["rows"]:
                row["assembly_name"] = "GRCh37"
            write_json(root / "report.json", report)
            with self.assertRaisesRegex(RuntimeError, "receipt-bound prepared genome"):
                run_preparation(root, "GRCh37")
            self.assertFalse((root / "output").exists())
            relocated = root / "relocated.json"
            relocated.write_bytes((root / "reference/genomes.json").read_bytes())
            reference["catalog_path"] = "/unavailable/original/catalog.json"
            PREPARE.validate_reference_assembly(reference, "GRCh38", relocated)
            relocated.write_text('{}')
            with self.assertRaisesRegex(RuntimeError, "catalog hash mismatch"):
                PREPARE.validate_reference_assembly(reference, "GRCh38", relocated)

    def test_report_and_svg_are_bound_not_just_gene_labels(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, report, svg = source_fixture(root)
            candidates = run_preparation(root)
            with self.assertRaisesRegex(RuntimeError, "another gene"):
                binding.load_bound_locus_report(candidates, root / "report.json", "OTHER")
            changed = {**report, "new_annotation": True}
            write_json(root / "report.json", changed)
            with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
                binding.load_bound_locus_report(candidates, root / "report.json")
            for replacement in (svg + "\n", svg.replace("toy_panel", "other_panel")):
                (root / "base.svg").write_text(replacement)
                with self.assertRaisesRegex(RuntimeError, "not preparation-bound"):
                    binding.load_bound_locus_svg(candidates, root / "base.svg", report)
            (root / "base.svg").write_text(svg)
            del candidates["source_bindings"]["locus_svgs"]
            with self.assertRaisesRegex(RuntimeError, "reprepare"):
                binding.load_bound_locus_svg(candidates, root / "base.svg", report)

    def test_per_region_report_binding_cannot_disagree(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_fixture(root)
            candidates = run_preparation(root)
            candidates["regions"][0]["source_region"]["source_report_sha256"] = "0" * 64
            with self.assertRaisesRegex(RuntimeError, "per-region"):
                binding.load_bound_locus_report(candidates, root / "report.json")

    def test_already_extended_svg_cannot_be_registered_as_original(self):
        with TemporaryDirectory() as tmp:
            _, report, svg = source_fixture(Path(tmp))
            svg = svg.replace("</svg>", '<g data-gentle-panel="tss-local-promoter-similarity"/></svg>')
            with self.assertRaisesRegex(RuntimeError, "already contains"):
                binding.validate_locus_svg(svg, report)

    def test_renderers_validate_before_writing_or_loading_plotting(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_fixture(root)
            candidates = run_preparation(root)
            regions = {row["region_id"]: row for row in candidates["regions"]}
            bound = (candidates, {}, regions, {}, {}, {}, {})
            common = ["--candidates-json", "unused", "--comparison", "unused",
                      "--matches", "unused", "--hits", "unused", "--locus-report", str(root / "report.json")]
            duplicate = ["render", *common, "--locus-report", str(root / "report.json"),
                         "--output", str(root / "compact")]
            with patch.object(sys, "argv", duplicate), patch.object(
                    RENDER, "load_bound_comparison", return_value=bound):
                with self.assertRaisesRegex(RuntimeError, "duplicate locus report"):
                    RENDER.main()
            self.assertFalse((root / "compact").exists())
            (root / "base.svg").write_text("wrong source")
            argv = ["append", *common, "--base-svg", str(root / "base.svg"), "--gene", "TOY",
                    "--output-svg", str(root / "tall/out.svg"), "--output-pdf", str(root / "tall/out.pdf"),
                    "--output-png", str(root / "tall/out.png")]
            with patch.object(sys, "argv", argv), patch.object(
                    APPEND, "load_bound_comparison", return_value=bound), patch.object(APPEND, "render_derivatives") as render:
                with self.assertRaisesRegex(RuntimeError, "not preparation-bound"):
                    APPEND.main()
                render.assert_not_called()
            self.assertFalse((root / "tall").exists())

    def test_output_paths_cannot_overwrite_sources_or_each_other(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.svg"
            source.write_text("bound source")
            for paths in ([source], [root / "new", root / "new"]):
                with self.assertRaisesRegex(ValueError, "distinct and absent"):
                    APPEND.require_new_outputs(paths)
            self.assertEqual(source.read_text(), "bound source")


class TssWindowTests(unittest.TestCase):
    def test_transcript_oriented_windows(self) -> None:
        self.assertEqual(PREPARE.selected_window(1_000, "+", 500, 200), (500, 1_200))
        self.assertEqual(PREPARE.selected_window(1_000, "-", 500, 200), (800, 1_500))
        with self.assertRaises(ValueError):
            PREPARE.selected_window(1_000, ".", 500, 200)

    def test_connected_stretches_merge_only_touching_windows(self) -> None:
        windows = [
            {"start_1based": 900, "end_1based": 1_200, "name": "b"},
            {"start_1based": 500, "end_1based": 900, "name": "a"},
            {"start_1based": 1_202, "end_1based": 1_500, "name": "c"},
        ]
        stretches = PREPARE.connected_stretches(windows)
        self.assertEqual([(row["start_1based"], row["end_1based"])
                          for row in stretches], [(500, 1_200), (1_202, 1_500)])
        self.assertEqual([row["name"] for row in stretches[0]["tss_windows"]], ["a", "b"])

    def test_receipt_bound_window_slice_is_assembly_forward(self) -> None:
        window = {"start_0based": "100", "end_0based_exclusive": "110", "strand": "+"}
        self.assertEqual(
            PREPARE.assembly_forward_slice("AACCGGTTAA", window, 103, 106),
            "CCGG",
        )
        assembly = "AACCGGTTAA"
        transcript_oriented = PREPARE.reverse_complement(assembly)
        window["strand"] = "-"
        self.assertEqual(
            PREPARE.assembly_forward_slice(transcript_oriented, window, 103, 106),
            "CCGG",
        )
        with self.assertRaises(RuntimeError):
            PREPARE.assembly_forward_slice("AACCGGTTAA", window, 99, 106)


class FrequencyTests(unittest.TestCase):
    def test_query_interval_projection_mirrors_without_changing_width(self):
        self.assertEqual(APPEND.query_interval_pixels(0, 20, 100, 255, 795, "+"), (255, 159))
        self.assertEqual(APPEND.query_interval_pixels(0, 20, 100, 255, 795, "-"), (891, 159))
        for strand in ("+", "-"):
            self.assertEqual(APPEND.query_interval_pixels(0, 100, 100, 255, 795, strand), (255, 795))
            self.assertEqual(APPEND.query_interval_pixels(0, 1, 1, 255, 795, strand), (255, 795))
        for start, end, strand in [(-1, 10, "+"), (0, 101, "-"), (2, 2, "+"), (0, 10, ".")]:
            with self.assertRaises(ValueError):
                APPEND.query_interval_pixels(start, end, 100, 255, 795, strand)

    def test_drawn_frequency_and_hsps_follow_negative_strand(self):
        with TemporaryDirectory() as tmp:
            _, _, base = source_fixture(Path(tmp))
            starts = {}
            for strand in ("+", "-"):
                fixture = plot_fixture(strand)
                output, _ = APPEND.append_section(base, "TOY", *fixture)
                tree = ET.fromstring(output)
                blocks = [node for node in tree.iter("rect")
                          if node.get("height") == "12" and node.get("width") != "795"]
                self.assertEqual(len(blocks), 3)  # One frequency segment, two ranked HSPs.
                self.assertEqual(len({node.get("x") for node in blocks}), 1)
                starts[strand] = float(blocks[0].get("x"))
                stretch = fixture[0]["stretches"][0]
                ends = sorted(APPEND.x_for(pos, stretch, strand) for pos in (300, 499))
                expected = ends[0] if strand == "+" else ends[0] + 0.9 * (ends[1] - ends[0])
                self.assertAlmostEqual(starts[strand], expected, places=2)
            self.assertGreater(starts["-"], starts["+"])

    def test_footer_clears_actual_content_for_dense_and_empty_layouts(self):
        with TemporaryDirectory() as tmp:
            _, _, base = source_fixture(Path(tmp))
            for rows in [(), (0,), (-1,), (2,), (0,) * 12, (2,) * 12, (0, 1, 2, -1) * 4]:
                with self.subTest(rows=rows):
                    fixture = plot_fixture(match_counts=rows)
                    output, height = APPEND.append_section(base, "TOY", *fixture)
                    self.assertEqual(APPEND.append_section(base, "TOY", *fixture), (output, height))
                    tree = ET.fromstring(output)
                    panel = next(node for node in tree if node.get("data-gentle-panel"))
                    last_baseline = max(float(node.get("y")) for node in panel.iter("text"))
                    footer = next(node for node in tree if node.get("data-gentle-shifted-footer"))
                    shift = float(footer.get("transform").split()[1].rstrip(")"))
                    footer_top = min(float(node.get("y")) - float(node.get("font-size", 13))
                                     for node in footer.iter("text")) + shift
                    self.assertGreaterEqual(footer_top - last_baseline, 19)
                    self.assertLess(max(float(node.get("y")) for node in footer.iter("text")) + shift, height)
                    self.assertEqual(tree.get("viewBox"), f"0 0 1400 {height}")
                    self.assertEqual(tree.find("rect").get("height"), str(height))

    def test_upper_stretch_references_share_ids_colours_and_bound_axis(self):
        with TemporaryDirectory() as tmp:
            _, report, base = source_fixture(Path(tmp))
            fixture = plot_fixture(match_counts=(0, 2))
            # Two separated stretches keep independent labels and rows.
            second = deepcopy(fixture[0]["stretches"][0])
            second.update(stretch_id="TOY_tss_stretch_2", start_1based=1000, end_1based=1200)
            second["tss_windows"][0]["tss_1based"] = 1100
            fixture[0]["stretches"].append(second)
            lower, old_height = APPEND.append_section(base, "TOY", *fixture)
            for left, right in [(1, 10000), (10000, 1)]:
                report.update(axis_left_genomic_1based=left, axis_right_genomic_1based=right)
                output, height = APPEND.add_stretch_overview(lower, "TOY", fixture[0], report)
                tree = ET.fromstring(output)
                overview = next(node for node in tree if node.get("data-gentle-tss-stretch-overview"))
                shifted = next(node for node in tree if node.get("data-gentle-shifted-locus"))
                backgrounds = shifted[0]
                self.assertEqual(backgrounds.get("data-gentle-tss-stretch-backgrounds"), "true")
                self.assertEqual(backgrounds.get("pointer-events"), "none")
                self.assertEqual(len(backgrounds), len(fixture[0]["stretches"]))
                similarity = next(node for node in shifted if node.get("data-gentle-panel"))
                similarity_y = float(similarity.find("text").get("y"))
                self.assertGreater(height, old_height)
                self.assertLess(output.index('data-gentle-tss-stretch-overview'),
                                output.index("Transcript models and annotation-derived metrics"))
                for index, stretch in enumerate(fixture[0]["stretches"]):
                    name = stretch["stretch_id"]
                    upper_link = next(node for node in overview if node.get("id") == f"overview-{name}")
                    lower_link = next(node for node in shifted.iter("a") if node.get("id") == f"similarity-{name}")
                    self.assertEqual(upper_link.get("href"), f"#similarity-{name}")
                    self.assertEqual(lower_link.get("href"), f"#overview-{name}")
                    bar = upper_link.find("rect")
                    self.assertEqual(bar.get("fill"), lower_link.find("text").get("fill"))
                    band = backgrounds[index]
                    self.assertEqual(band.get("data-gentle-tss-stretch-band"), name)
                    for key in ("x", "width", "fill", "data-gentle-genomic-start", "data-gentle-genomic-end"):
                        self.assertEqual(band.get(key), bar.get(key))
                    self.assertEqual(band.get("fill-opacity"), "0.12")
                    self.assertLess(float(band.get("y")), 110)  # Original transcript lane.
                    self.assertGreater(float(band.get("y")) + float(band.get("height")), 110)
                    self.assertLess(float(band.get("y")) + float(band.get("height")), similarity_y)
                    expected = min(255 + (pos - left) / (right - left) * 795
                                   for pos in (stretch["start_1based"], stretch["end_1based"]))
                    self.assertAlmostEqual(float(bar.get("x")), expected, places=2)
                original_line = ET.fromstring(base).find("line")
                self.assertIn(ET.tostring(original_line), ET.tostring(shifted))
                self.assertEqual((output, height), APPEND.add_stretch_overview(lower, "TOY", fixture[0], report))
            report.update(axis_left_genomic_1based=1, axis_right_genomic_1based=150)
            with self.assertRaisesRegex(ValueError, "outside"):
                APPEND.add_stretch_overview(lower, "TOY", fixture[0], report)

    def test_background_bands_preserve_native_feature_and_occupancy_elements(self):
        with TemporaryDirectory() as tmp:
            _, report, base = source_fixture(Path(tmp))
            # Synthetic glyphs exercise layering/preservation, not experimental signal validity.
            glyphs = ('<rect data-gentle-exon="E1" x="280" y="104" width="20" height="12" fill="#2563eb"/>'
                      '<rect data-gentle-regulatory-feature="F1" x="290" y="122" width="30" height="6" fill="#61d9a8"/>'
                      '<line data-gentle-occupancy-lane="O1" data-gentle-occupancy-state="available" '
                      'x1="255" x2="1050" y1="140" y2="140"/>'
                      '<rect data-gentle-occupancy-interval="I1" x="300" y="134" width="8" height="6" fill="#475569"/>')
            base = base.replace('<text data-gentle-overlay-non-claims=', glyphs + '<text data-gentle-overlay-non-claims=')
            fixture = plot_fixture()
            lower, _ = APPEND.append_section(base, "TOY", *fixture)
            # The real renderer emits a namespaced SVG; support that as well as tiny test SVGs.
            lower = lower.replace('<svg ', '<svg xmlns="http://www.w3.org/2000/svg" ', 1)
            output, _ = APPEND.add_stretch_overview(lower, "TOY", fixture[0], report)
            before = ET.fromstring(lower)
            after = ET.fromstring(output)
            original = list(before)[2:]  # Root background and leading title remain outside the shift.
            shifted = next(node for node in after if node.get("data-gentle-shifted-locus"))
            self.assertEqual([ET.tostring(node) for node in original],
                             [ET.tostring(node) for node in list(shifted)[1:]])
            self.assertIn(glyphs, output)

    def test_background_bands_require_a_valid_similarity_boundary(self):
        with TemporaryDirectory() as tmp:
            _, report, base = source_fixture(Path(tmp))
            fixture = plot_fixture()
            with self.assertRaisesRegex(ValueError, "similarity boundary"):
                APPEND.add_stretch_overview(base, "TOY", fixture[0], report)
            lower, _ = APPEND.append_section(base, "TOY", *fixture)
            lower = lower.replace('x="34.00" y="170.00"', 'x="34.00" y="60.00"', 1)
            with self.assertRaisesRegex(ValueError, "overview height"):
                APPEND.add_stretch_overview(lower, "TOY", fixture[0], report)

    def test_frequency_segments_count_distinct_genes_not_transcripts(self) -> None:
        rows = [
            {
                "promoter_id": "p1",
                "gene_ids": "g1",
                "query_intervals_0based_half_open": "0-10;20-30",
            },
            {
                "promoter_id": "p2",
                "gene_ids": "g1;g2",
                "query_intervals_0based_half_open": "5-25",
            },
        ]
        self.assertEqual(
            RENDER.frequency_segments(rows),
            [(0, 5, 1), (5, 10, 2), (10, 20, 2), (20, 25, 2), (25, 30, 1)],
        )

    def test_interval_merge_is_half_open(self) -> None:
        self.assertEqual(
            RENDER.merge_intervals([(8, 12), (0, 5), (5, 9), (20, 21)]),
            [(0, 12), (20, 21)],
        )

    def test_reverse_chain_decreases_in_query_without_false_break(self) -> None:
        reverse = [
            {"sstart": 10, "send": 20, "qstart": 90, "qend": 80, "bitscore": 20},
            {"sstart": 30, "send": 40, "qstart": 60, "qend": 50, "bitscore": 20},
        ]
        self.assertEqual([broken for _, broken in APPEND.ordered_blocks(reverse)],
                         [False, False])
        reordered = [reverse[0], {**reverse[1], "qstart": 100, "qend": 95}]
        self.assertEqual([broken for _, broken in APPEND.ordered_blocks(reordered)],
                         [False, True])

    def test_similarity_panel_is_inserted_before_interpretation_footer(self) -> None:
        base = (
            '<svg height="100" viewBox="0 0 1400 100" width="1400">\n'
            '<rect fill="#ffffff" height="100" width="1400" x="0" y="0"/>\n'
            '<text x="34" y="20">existing evidence</text>\n'
            '<text data-gentle-overlay-non-claims="true" x="34" y="70">\n'
            'Reporter interpretation boundaries\n</text>\n'
            '<text x="34" y="90">Evidence provenance</text>\n</svg>\n'
        )
        candidates = {
            "stretches": [{
                "stretch_id": "G_tss_stretch_1", "gene": "G",
                "start_1based": 100, "end_1based": 800,
                "tss_windows": [{"strand": "+", "chromosome": "1", "tss_1based": 600}],
            }],
            "regions": [],
        }
        output, height = APPEND.append_section(
            base, "G", candidates, {}, {}, {},
            {"minimum_bp": 40, "minimum_identity": 80.0, "maximum_evalue": 1e-5},
        )
        self.assertGreater(height, 100)
        self.assertLess(output.index("TSS-local Ensembl-feature promoterome similarity"),
                        output.index("Reporter interpretation boundaries"))
        self.assertIn('data-gentle-shifted-footer="true"', output)
        self.assertIn(f'height="{height}"', output)


class RetainedBundleTests(unittest.TestCase):
    def test_manifest_and_derivative_receipts(self) -> None:
        entries = {}
        for line in (BUNDLE / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            entries[name] = digest
        self.assertGreaterEqual(len(entries), 20)
        for name, digest in entries.items():
            self.assertEqual(sha256(BUNDLE / name), digest, name)
        for gene in ("CD44", "TGFB1", "SERPINE1"):
            stem = f"{gene}_luciferase_planning_EnsemblReg_15TF_with_TSS_similarity"
            receipt = json.loads((BUNDLE / f"{stem}.receipt.json").read_text())
            for name, expected in receipt["outputs"].items():
                self.assertEqual(f"sha256:{sha256(BUNDLE / name)}", expected, name)

    def test_retained_comparison_uses_corrected_policy_and_limit_audit(self) -> None:
        comparison = json.loads((BUNDLE / "comparison.json").read_text())
        self.assertEqual(comparison["schema"], RENDER.comparison_tools.SCHEMA)
        self.assertEqual(comparison["counting_policy_id"],
                         RENDER.comparison_tools.COUNTING_POLICY)
        queries = comparison["tasks"]["blastn"]["queries"]
        self.assertEqual(len(queries), 15)
        for row in queries:
            self.assertTrue(row["other_gene_exclusion_verified"])
            self.assertEqual(row["counts_are_lower_bounds"],
                             row["target_cap_reached"] or row["hsp_cap_reached"])

    def test_retained_reports_have_validated_nonempty_cutrun_lanes(self) -> None:
        validation = json.loads((BUNDLE / "cutrun_lane_validation.json").read_text())
        self.assertEqual(validation["schema"],
                         "gentle.tp73_locus_cutrun_lane_validation.v1")
        self.assertEqual(validation["assembly_id"], "GRCh38")
        self.assertEqual(validation["total_lane_count"], 36)
        self.assertGreater(validation["total_rendered_interval_count"], 0)
        self.assertEqual({row["gene"] for row in validation["genes"]},
                         {"CD44", "TGFB1", "SERPINE1"})
        for row in validation["genes"]:
            self.assertEqual(row["lane_count"], 12)
            self.assertGreater(row["rendered_interval_count"], 0)
            self.assertEqual(len(row["lanes"]), 12)
            for lane in row["lanes"]:
                self.assertEqual(lane["source_assembly"], "GRCh38")
                self.assertGreater(lane["source_overlapping_interval_count"], 0)
                self.assertGreater(lane["rendered_interval_count"], 0)
                self.assertRegex(lane["source_sha256"], r"^sha256:[0-9a-f]{64}$")
            stem = f'{row["gene"]}_luciferase_planning_EnsemblReg_15TF_with_TSS_similarity.svg'
            svg = (BUNDLE / stem).read_text()
            self.assertEqual(svg.count('data-gentle-occupancy-state="available"'), 12)
            self.assertNotIn('data-gentle-occupancy-state="no_compatible_interval"', svg)

    def test_tall_reports_place_similarity_before_original_footer(self) -> None:
        for gene in ("CD44", "TGFB1", "SERPINE1"):
            svg = (BUNDLE / f"{gene}_luciferase_planning_EnsemblReg_15TF_with_TSS_similarity.svg").read_text()
            self.assertLess(svg.index("TSS-local Ensembl-feature promoterome similarity"),
                            svg.index("Reporter interpretation boundaries"))
            self.assertIn('data-gentle-shifted-footer="true"', svg)


if __name__ == "__main__":
    unittest.main()
