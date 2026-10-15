"""Offline semantic checks for the existing synthetic rack documentation figures.

No Rust, renderer, network or native GUI runs. Existing fixture provenance and
regeneration commands live in docs/figures/README.md. Temporary hero outputs
exercise the Python presentation transform, not current Rust-export freshness.
"""

from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "docs/figures"


class RackFigureReviewTests(unittest.TestCase):
    def test_scad_openings_retain_bottom_floor_in_all_templates(self):
        cases = (("cell_culture_plate", 6, 1.0),
                 ("gibson_single_insert_storage_rack", 24, 1.2),
                 ("gibson_single_insert_pipetting_rack", 24, 2.0))
        for stem, count, floor in cases:
            text = (FIGURES / f"{stem}.scad").read_text()
            for source in (text, text.replace("\n", "\r\n")):
                with self.subTest(stem=stem, crlf="\r" in source):
                    dimensions = dict(re.findall(r"^(\w+) = ([0-9.]+);", source, re.M))
                    self.assertEqual(float(dimensions["floor_thickness"]), floor)
                    self.assertLess(floor, float(dimensions["rack_height"]))
                    cuts = re.findall(r"translate\(\[([0-9.]+), ([0-9.]+), floor_thickness\]\) cylinder\(h=rack_height - floor_thickness \+ 0.10, d=opening_diameter\);", source)
                    self.assertEqual(len(cuts), count)
                    self.assertNotIn(", -0.05]) cylinder", source)
                    for x, y in cuts:
                        radius = float(dimensions["opening_diameter"]) / 2
                        self.assertGreaterEqual(float(x) - radius, 0)
                        self.assertLessEqual(float(x) + radius, float(dimensions["outer_width"]))
                        self.assertGreaterEqual(float(y) - radius, 0)
                        self.assertLessEqual(float(y) + radius, float(dimensions["outer_depth"]))
        exporter = (ROOT / "src/engine/state/lineage_containers.rs").read_text()
        self.assertIn(", floor_thickness]) cylinder(h=rack_height - floor_thickness + 0.10", exporter)
        self.assertNotIn(", -0.05]) cylinder", exporter)

    def test_hero_transform_preserves_existing_isometric_surface_coordinates(self):
        ns = {"s": "http://www.w3.org/2000/svg"}
        for stem in ("gibson_single_insert_rack_isometric", "cell_culture_plate_isometric"):
            source = FIGURES / f"{stem}.svg"
            original = ET.parse(source).getroot()
            with tempfile.TemporaryDirectory() as temp:
                output = Path(temp) / "hero.svg"
                result = subprocess.run([sys.executable, str(FIGURES / "render_rack_isometric_hero.py"), str(source), str(output)], capture_output=True, text=True, timeout=15)
                self.assertEqual(result.returncode, 0, result.stderr)
                hero = ET.parse(output).getroot()
            self.assertEqual(hero.attrib["data-rack-isometric-template"], original.attrib["data-rack-isometric-template"])
            self.assertEqual([e.attrib["points"] for e in original.findall("s:polygon", ns)][:3], [e.attrib["points"] for e in hero.findall("s:polygon", ns)][:3])
            self.assertIn("data-rack-hero-background", ET.tostring(hero, encoding="unicode"))


if __name__ == "__main__":
    unittest.main()
