"""Offline SCAD-template checks, not execution of GENtle's Rust exporter.

Literal dimensions are hand-transcribed from the drawing sources documented in
data/labware/README.md. The actual body template is read from the Rust renderer;
no duplicate model implementation is used. OpenSCAD is optional, not installed
by this test. Mesh bounds are visualization checks, never physical fit evidence.
"""

import math
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CASES = (
    ("0_5ml", 30.0, 6.7, 7.8, 16.5, False),
    ("1_5ml", 38.9, 9.85, 10.8, 20.0, False),
    ("2_0ml", 40.0, 9.85, 10.8, 32.3, True),
)


def body_template(source):
    match = re.search(r'const TUBE_BODY_SCAD: &str = r#"(.*?)"#;', source, re.S)
    if not match:
        raise ValueError("SCAD body template not found; update the harness explicitly")
    return match.group(1)


def openscad_binary():
    explicit = os.environ.get("GENTLE_OPENSCAD_BIN")
    if explicit:
        # A misspelled explicit path must fail rather than silently skip.
        return explicit
    native = Path("/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD")
    return str(native) if native.is_file() else shutil.which("openscad")


def carrier_template(source, rows, columns, width, depth, height, opening, edge, front,
                     pitch=20, outer_wall=None):
    """Exercise actual SCAD format strings, not the unexecuted Rust derivation."""
    source = source.split("pub(super) fn export_rack_openscad(", 1)[1].split(
        "pub(super) fn export_rack_simulation_json(", 1)[0]
    fragments = re.findall(
        r'(?:let mut scad = format!|scad.push_str\(&format!)\(\s*("(?:\\.|[^"\\])*")', source)
    if len(fragments) != 3:
        raise ValueError("carrier SCAD format strings changed; review the harness")
    header, bore, strip = [json.loads(value).replace("{:.3}", "{:.3f}") for value in fragments]
    end = re.search(r'scad.push_str\(("(?:\\.|[^"\\])*")\)', source)
    if not end:
        raise ValueError("carrier SCAD terminator missing")
    result = header.format("synthetic", "unverified", "tube_rack", "reference_only",
                           width, depth, height, opening, 2.0, front, 6.0, 0.8, 2.0)
    for row in range(rows):
        for column in range(columns):
            result += bore.format(edge + opening / 2 + column * pitch,
                                  front + opening / 2 + row * pitch)
    outer_wall = edge if outer_wall is None else outer_wall
    result += strip.format(outer_wall / 2, 0.8, width - outer_wall, max(front - 1.6, 0))
    return result + json.loads(end.group(1))


class LabwareTemplateTests(unittest.TestCase):
    def test_compact_carrier_scad_retains_floor_and_variable_grid_coordinates(self):
        # Synthetic geometry, not output from Rust or fit to a storage box.
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text()
        for rows, columns in ((9, 9), (8, 8), (10, 10), (5, 12)):
            edge, opening, pitch, front = 4.775, 10.45, 14, 10
            width = 2 * edge + opening + (columns - 1) * pitch
            depth = front + edge + opening + (rows - 1) * pitch
            scad = carrier_template(source, rows, columns, width, depth, 24,
                                    opening, edge, front, pitch)
            self.assertEqual(scad.count(" cylinder("), rows * columns)
            self.assertEqual(scad.count(", floor_thickness]) cylinder(h=rack_height - floor_thickness + 0.10"), rows * columns)
            self.assertIn(f"outer_width = {width:.3f};", scad)
            if (rows, columns) == (9, 9):
                self.assertEqual(width, 132)
                self.assertIn("translate([10.000, 15.225, floor_thickness])", scad)
                self.assertIn("translate([122.000, 127.225, floor_thickness])", scad)

    def test_carrier_scad_fragments_preserve_floor_open_top_and_newline_semantics(self):
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text()
        args = (4, 9, 180.0, 85.075, 44.0, 10.15, 4.925, 10.0)
        scad = carrier_template(source, *args)
        self.assertEqual(scad, carrier_template(source.replace("\n", "\r\n"), *args))
        self.assertEqual(scad.count(", floor_thickness]) cylinder(h=rack_height - floor_thickness + 0.10"), 36)
        self.assertIn("translate([10.000, 15.075, floor_thickness])", scad)
        self.assertIn("translate([170.000, 75.075, floor_thickness])", scad)
        self.assertNotIn(", -0.05]) cylinder", scad)
        pcr = carrier_template(source, 8, 12, 127.76, 85.48, 15.0, 6.6, 11.08, 7.94, 9, 2)
        self.assertEqual(pcr.count(" cylinder("), 96)
        self.assertIn("translate([14.380, 11.240, floor_thickness])", pcr)
        self.assertIn("translate([113.380, 74.240, floor_thickness])", pcr)

    def test_carrier_meshes_have_expected_outer_bounds(self):
        binary = openscad_binary()
        if not binary:
            self.skipTest("OpenSCAD unavailable; source-matched exporter acceptance remains separate")
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text()
        # Hand-derived engineering assumptions, not generated Rust output or fit evidence.
        cases = [("collar", 4, 9, 180.0, 85.075, 44.0, 10.15, 4.925, 10.0),
                 ("floor", 4, 9, 180.0, 85.225, 24.0, 10.45, 4.775, 10.0),
                 ("pcr96", 8, 12, 127.76, 85.48, 15.0, 6.6, 11.08, 7.94, 9, 2)]
        with tempfile.TemporaryDirectory(prefix="gentle carrier models ") as directory:
            for name, *args in cases:
                with self.subTest(carrier=name):
                    scad = Path(directory) / f"{name}.scad"
                    mesh = scad.with_suffix(".off")
                    scad.write_text(carrier_template(source, *args), encoding="utf-8")
                    result = subprocess.run([binary, "-o", str(mesh), str(scad)],
                                            capture_output=True, text=True, timeout=120)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertNotIn("ERROR:", result.stderr)
                    tokens = mesh.read_text(encoding="utf-8").split()
                    self.assertEqual(tokens[0], "OFF")
                    count, faces, _ = map(int, tokens[1:4])
                    self.assertGreater(count, 0)
                    self.assertGreater(faces, 0)
                    points = [tuple(map(float, tokens[start:start + 3]))
                              for start in range(4, 4 + 3 * count, 3)]
                    self.assertTrue(all(math.isfinite(value) for point in points for value in point))
                    for axis, expected in enumerate(args[2:5]):
                        self.assertAlmostEqual(min(point[axis] for point in points), 0, places=4)
                        self.assertAlmostEqual(max(point[axis] for point in points), expected, places=3)

    def test_template_extraction_accepts_both_checkout_newlines(self):
        source = (ROOT / "crates/gentle-render/src/labware.rs").read_text()
        lf = body_template(source)
        crlf = body_template(source.replace("\n", "\r\n"))
        self.assertEqual(lf, crlf.replace("\r\n", "\n"))
        self.assertIn("module gentle_tube_outer_body()", lf)
        self.assertTrue(lf.endswith("gentle_tube_outer_body();\n"))
        with self.assertRaises(ValueError):
            body_template("missing")

    def test_openscad_meshes_have_the_drawing_outer_bounds(self):
        binary = openscad_binary()
        if not binary:
            self.skipTest("OpenSCAD unavailable; Rust exporter and native acceptance are separate")
        source = (ROOT / "crates/gentle-render/src/labware.rs").read_text()
        body = body_template(source)
        with tempfile.TemporaryDirectory(prefix="gentle tube models ") as directory:
            for name, length, diameter, collar, cylindrical, rounded in CASES:
                with self.subTest(profile=name):
                    scad = Path(directory) / f"{name}.scad"
                    mesh = scad.with_suffix(".off")
                    scad.write_text(
                        f"$fn=64; body_length_mm={length}; body_diameter_mm={diameter};\n"
                        f"collar_diameter_mm={collar}; cylindrical_length_mm={cylindrical};\n"
                        f"rounded_bottom={str(rounded).lower()};\n"
                        "visual_collar_thickness_mm=1.0; visual_tip_radius_mm=0.6;\n"
                        + body,
                        encoding="utf-8",
                    )
                    result = subprocess.run(
                        [binary, "-o", str(mesh), str(scad)],
                        capture_output=True, text=True, timeout=120,
                    )
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertNotIn("ERROR:", result.stderr)
                    # OFF permits counts on the header line or the next line.
                    tokens = mesh.read_text(encoding="utf-8").split()
                    self.assertEqual(tokens[0], "OFF")
                    vertex_count, face_count, _ = map(int, tokens[1:4])
                    self.assertGreater(vertex_count, 0)
                    self.assertGreater(face_count, 0)
                    points = [tuple(map(float, tokens[start:start + 3]))
                              for start in range(4, 4 + 3 * vertex_count, 3)]
                    self.assertEqual(len(points), vertex_count)
                    self.assertTrue(all(len(point) == 3 for point in points))
                    self.assertTrue(all(math.isfinite(value) for point in points for value in point))
                    bounds = [(min(p[axis] for p in points), max(p[axis] for p in points))
                              for axis in range(3)]
                    for axis in (0, 1):
                        self.assertAlmostEqual(bounds[axis][0], -collar / 2, places=4)
                        self.assertAlmostEqual(bounds[axis][1], collar / 2, places=4)
                    # Faceted spheres may stop fractionally above the analytic tip.
                    self.assertGreaterEqual(bounds[2][0], -0.001)
                    self.assertLess(bounds[2][0], 0.02)
                    self.assertAlmostEqual(bounds[2][1], length, places=4)


if __name__ == "__main__":
    unittest.main()
