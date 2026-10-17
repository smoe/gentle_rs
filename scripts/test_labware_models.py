"""Offline SCAD-template checks, not execution of GENtle's Rust exporter.

Literal dimensions are hand-transcribed from the drawing sources documented in
data/labware/README.md. The actual body template is read from the Rust renderer;
no duplicate model implementation is used. OpenSCAD is optional, not installed
by this test. Mesh bounds are visualization checks, never physical fit evidence.
"""

import math
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


class LabwareTemplateTests(unittest.TestCase):
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
