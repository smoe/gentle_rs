"""Source/fixture checks for bounded manual tutorial additions.

These tests do not execute GENtle, DuckDB, a GUI or a biological analysis.
The arithmetic oracle and drawn PNG are explicitly synthetic teaching inputs.
"""

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import unittest
import zlib


ROOT = Path(__file__).resolve().parents[1]
GEL = ROOT / "docs/tutorial/inputs/gel_sizing"
GENERATOR = ROOT / "docs/tutorial/inputs/make_gel_sizing_fixture.py"
spec = importlib.util.spec_from_file_location("gel_teaching_generator", GENERATOR)
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
motif_spec = importlib.util.spec_from_file_location(
    "regulatory_teaching_generator",
    ROOT / "test_files/fixtures/regulatory_motif_subset/make_fixture.py",
)
motif_generator = importlib.util.module_from_spec(motif_spec)
motif_spec.loader.exec_module(motif_generator)


class TutorialGapLessonTests(unittest.TestCase):
    def test_png_is_original_eight_bit_grayscale_with_six_drawn_bands(self):
        payload = (GEL / "gel_sizing_synthetic.png").read_bytes()
        self.assertEqual(payload, generator.png_bytes())
        self.assertEqual(payload[:8], b"\x89PNG\r\n\x1a\n")
        offset, compressed, kinds = 8, b"", []
        while offset < len(payload):
            length = struct.unpack_from(">I", payload, offset)[0]
            kind = payload[offset + 4:offset + 8]
            data = payload[offset + 8:offset + 8 + length]
            crc = struct.unpack_from(">I", payload, offset + 8 + length)[0]
            self.assertEqual(crc, zlib.crc32(kind + data))
            kinds.append(kind)
            if kind == b"IHDR":
                self.assertEqual(struct.unpack(">IIBBBBB", data), (160, 220, 8, 0, 0, 0, 0))
            elif kind == b"IDAT":
                compressed += data
            offset += length + 12
        self.assertEqual(kinds, [b"IHDR", b"IDAT", b"IEND"])
        rows = zlib.decompress(compressed)
        self.assertEqual(len(rows), 220 * 161)
        centers = [(20, 20), (20, 80), (20, 140), (20, 200), (90, 110), (90, 210)]
        for y in range(220):
            self.assertEqual(rows[y * 161], 0, "unfiltered original pixels")
            for x in range(160):
                band = any(abs(x - bx) <= 8 and abs(y - by) <= 1 for bx, by in centers)
                self.assertEqual(rows[y * 161 + 1 + x], 210 if band else 12)

    def test_fixture_requests_and_raw_hashes_recreate_exactly(self):
        files = generator.fixture_files()
        for name, expected in files.items():
            self.assertEqual((GEL / name).read_bytes(), expected, name)
        for line in files["SHA256SUMS"].decode().splitlines():
            digest, name = line.split("  ", 1)
            self.assertEqual(hashlib.sha256(files[name]).hexdigest(), digest)
        image_hash = hashlib.sha256(files["gel_sizing_synthetic.png"]).hexdigest()
        for name in ("dna-analysis.json", "protein-analysis.json"):
            request = json.loads(files[name])
            self.assertEqual(request["image_sha256"], image_hash)
            self.assertEqual(request["image_id"], "gel_teaching")
            self.assertEqual(request["migration"], "down")
        self.assertTrue(json.loads(files["protein-analysis.json"])["ladder"]["prestained"])

    def test_generator_refuses_clobber_and_check_never_repairs_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "gel fixture with 'quote'"
            command = [sys.executable, str(GENERATOR), str(output)]
            self.assertEqual(subprocess.run(command, capture_output=True).returncode, 0)
            before = {p.name: p.read_bytes() for p in output.iterdir()}
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)
            self.assertEqual({p.name: p.read_bytes() for p in output.iterdir()}, before)
            self.assertEqual(subprocess.run(command + ["--check"], capture_output=True).returncode, 0)
            damaged = output / "dna-analysis.json"
            damaged.write_bytes(b"not the original request")
            self.assertNotEqual(subprocess.run(command + ["--check"], capture_output=True).returncode, 0)
            self.assertEqual(damaged.read_bytes(), b"not the original request")

    def test_independent_log_size_oracle_and_outside_range(self):
        dna = json.loads((GEL / "dna-analysis.json").read_bytes())
        protein = json.loads((GEL / "protein-analysis.json").read_bytes())
        for request, expected in [(dna, math.sqrt(1000 * 500)), (protein, math.sqrt(50 * 25))]:
            refs = request["ladder"]["bands"]
            a, b = refs[1], refs[2]
            y = request["sample_bands"][0]["center"]["y"]
            weight = (y - a["center"]["y"]) / (b["center"]["y"] - a["center"]["y"])
            estimate = 10 ** ((1 - weight) * math.log10(a["size"]) + weight * math.log10(b["size"]))
            self.assertAlmostEqual(estimate, expected, places=9)
            outside = request["sample_bands"][1]["center"]["y"]
            self.assertGreater(outside, max(ref["center"]["y"] for ref in refs))
        self.assertNotAlmostEqual(math.sqrt(1000 * 500), 750)

    def test_microarray_rows_and_coordinate_oracle_are_vendor_subset_not_measurements(self):
        folder = ROOT / "test_files/fixtures/microarray_tracks"
        manifest = json.loads((folder / "clariomd.tp73_vendor_subset.manifest.json").read_bytes())
        self.assertIn("synthetic", " ".join(manifest["warnings"]))
        self.assertEqual(manifest["coordinate_system"], "hg38")
        rows = []
        for contrast in manifest["contrasts"]:
            with (folder / contrast["path"]).open(encoding="utf-8", newline="") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(parsed), contrast["row_count"])
            rows.extend(parsed)
        self.assertEqual(len(rows), 5)
        first = rows[0]
        self.assertEqual(first["feature_id"], "PSR0100145779.hg.1")
        self.assertEqual((int(first["start_1based"]) - 3652516, int(first["end_1based"]) - 3652516 + 1), (4, 49))

    def test_new_manual_sources_do_not_invent_executable_gui_acceptance(self):
        manifest = json.loads((ROOT / "docs/tutorial/manifest.json").read_bytes())
        for source_name in ("07-04_microarray_evidence_gui_cli.json", "10-03_measured_gel_sizing_gui_cli.json",
                            "08-18_repeat_environment_cohorts.json"):
            source = json.loads((ROOT / "docs/tutorial/sources" / source_name).read_bytes())
            self.assertEqual(source["catalog"]["status"], "manual/hybrid")
            self.assertNotIn("generated", source)
            self.assertNotIn("gui_acceptance", source)
            self.assertNotIn(source["id"], [chapter["id"] for chapter in manifest["chapters"]])

    def test_lesson_local_links_and_visible_pending_limits(self):
        for filename in ("07-04_microarray_evidence_gui_cli.md", "10-03_measured_gel_sizing_gui_cli.md",
                         "08-16_tss_regulatory_view_gui.md", "08-18_repeat_environment_cohorts.md"):
            path = ROOT / "docs/tutorial" / filename
            prose = path.read_text(encoding="utf-8")
            for link in re.findall(r"\]\(([^)]+)\)", prose):
                if not link.startswith(("https://", "http://", "#")):
                    self.assertTrue((path.parent / link.split("#")[0]).exists(), link)
            self.assertIn("synthetic", prose)
            self.assertIn("pending", prose.lower())

    def test_repeat_fixture_keeps_identical_dna_at_distinct_locations(self):
        fixture = ROOT / "test_files/fixtures/repeat_cohort/local_project.rmsk.txt"
        rows = [line.split() for line in fixture.read_text().splitlines()]
        self.assertEqual(len(rows), 3)
        selected = [row for row in rows if row[11:13] == ["LINE", "L1"]]
        self.assertEqual(len(selected), 2)
        sequence = "".join(line.strip() for line in
                           (ROOT / "test_files/fixtures/genomes/AB011549.2.fa").read_text().splitlines()
                           if not line.startswith(">"))
        intervals = [(int(row[6]), int(row[7])) for row in selected]
        self.assertEqual(intervals, [(1843, 1849), (3877, 3883)])
        self.assertEqual([sequence[a:b].upper() for a, b in intervals], ["AGCCAG", "AGCCAG"])
        self.assertEqual([row[9] for row in selected], ["+", "-"])
        prose = (ROOT / "docs/tutorial/08-18_repeat_environment_cohorts.md").read_text()
        self.assertIn("RNA-read ranking is not implemented", prose)
        self.assertIn("not evidence of mammalian LINE biology", prose)

    def test_tss_companion_footprints_match_original_geometry_without_faking_receipts(self):
        original = json.loads((ROOT / "docs/tutorial/reproducibility/tss_regulatory_view_gui/synthetic-minus.profile-report.json").read_bytes())
        geometry = original["windows"][0]["record"]["geometry"]
        self.assertEqual(original["panel_resolution"]["panel"]["score_kind"], "llr_bits")
        self.assertEqual((geometry["start_1based"], geometry["end_1based"], geometry["tss_1based"]),
                         (293, 303, 300))
        hits = motif_generator.tss_view_hits()
        self.assertEqual(len(hits), 5)
        low, high = geometry["start_1based"] - 1, geometry["end_1based"]
        local = [(high - (start + 6), high - start, score) for start, _, score, _ in hits
                 if start < high and start + 6 > low]
        self.assertEqual(local, [(3, 9, 1.25), (1, 7, 0.0), (7, 13, -0.5)])
        with tempfile.TemporaryDirectory() as tmp:
            nonexistent = Path(tmp) / "must not exist"
            with self.assertRaises(ValueError):
                motif_generator.build(nonexistent, "not-invoked", "unknown")
            self.assertFalse(nonexistent.exists())
            with self.assertRaises(FileExistsError):
                motif_generator.build(Path(tmp), "not-invoked", "tss-view")
        draft = ROOT / "docs/tutorial/drafts/flna_splicing_interpretation.md"
        self.assertIn("Inputs Pending", draft.read_text())
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_bytes())
        self.assertNotIn(str(draft.relative_to(ROOT)), json.dumps(catalog))

    def test_navigation_does_not_duplicate_review_badges_or_planned_existing_variant_guide(self):
        readme = (ROOT / "docs/tutorial/README.md").read_text(encoding="utf-8")
        self.assertNotRegex(readme, r"review `(?:unreviewed|codex_reviewed|human_reviewed)`")
        self.assertIn("review_manifest.json", readme)
        landscape = (ROOT / "docs/tutorial/landscape_overview.md").read_text(encoding="utf-8")
        self.assertNotIn("Planned: SNP retrieval", landscape)
        self.assertIn("08-04_vkorc1_warfarin_promoter_luciferase_gui.md", landscape)


if __name__ == "__main__":
    unittest.main()
