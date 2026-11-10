"""Offline labware source checks, not execution of GENtle's Rust paths.

Literal dimensions are hand-transcribed from the drawing sources documented in
data/labware/README.md. The actual body template is read from the Rust renderer;
no duplicate model implementation is used. OpenSCAD is optional, not installed
by this test. Mesh bounds are visualization checks, never physical fit evidence.
"""

import hashlib
import math
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET


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


def strip_template(source, count):
    """Actual strip SCAD with literal BRAND callouts/declared visual assumptions."""
    match = re.search(r'const PCR_STRIP_SCAD: &str = r#"(.*?)"#;', source, re.S)
    if not match:
        raise ValueError("PCR strip template not found; review the harness")
    return (
        f"$fn=64; tube_count={count}; max_tubes=8; center_pitch_mm=9;\n"
        "body_length_mm=20.75; body_diameter_mm=6; collar_diameter_mm=8;\n"
        "rim_diameter_mm=7.05; tip_diameter_mm=2.75;\n"
        "visual_shoulder_fraction=0.55; visual_collar_drop_mm=4.5;\n"
        "visual_collar_thickness_mm=0.9; visual_rim_thickness_mm=0.6;\n"
        "visual_bridge_width_mm=1.0; visual_bridge_height_mm=0.7;\n"
        + match.group(1)
    )


def openscad_binary():
    explicit = os.environ.get("GENTLE_OPENSCAD_BIN")
    if explicit:
        # A misspelled explicit path must fail rather than silently skip.
        return explicit
    native = Path("/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD")
    return str(native) if native.is_file() else shutil.which("openscad")


def carrier_template(source, rows, columns, width, depth, height, opening, edge, front,
                     pitch=20, outer_wall=None, *, safe_lock=False,
                     label_depth=6.0, corner_radius=2.0):
    """Exercise actual SCAD format strings, not the unexecuted Rust derivation."""
    source = source.replace("\r\n", "\n")
    rounded = re.search(r'const SAFE_LOCK_RACK_BODY_SCAD: &str = r#"(.*?)"#;', source, re.S)
    source = source.split("pub(super) fn export_rack_openscad(", 1)[1].split(
        "pub(super) fn export_rack_simulation_json(", 1)[0]
    # Only the original bore generator; the new storage-box branch is separate.
    source = source.split("let blocked = rack", 1)[1]
    fragments = re.findall(
        r'(?:let mut scad = format!|scad.push_str\(&format!)\(\s*("(?:\\.|[^"\\])*")', source)
    if len(fragments) != 3:
        raise ValueError("carrier SCAD format strings changed; review the harness")
    header, bore, strip = [json.loads(value).replace("{:.3}", "{:.3f}") for value in fragments]
    end = re.search(r'scad.push_str\(("(?:\\.|[^"\\])*")\)', source)
    if not end:
        raise ValueError("carrier SCAD terminator missing")
    body = re.search(
        r'scad.push_str\(if safe_lock \{\s*SAFE_LOCK_RACK_BODY_SCAD\s*\} else \{\s*'
        r'("(?:\\.|[^"\\])*")\s*\}\);', source)
    if not body or not rounded:
        raise ValueError("carrier body selection changed; review the harness")
    result = header.format("synthetic", "unverified", "tube_rack", "reference_only",
                           width, depth, height, opening, 2.0, front,
                           label_depth, 0.8, corner_radius)
    result += rounded.group(1) if safe_lock else json.loads(body.group(1))
    for row in range(rows):
        for column in range(columns):
            result += bore.format(edge + opening / 2 + column * pitch,
                                  front + opening / 2 + row * pitch)
    outer_wall = edge if outer_wall is None else outer_wall
    strip_depth = "front_label_strip_depth" if safe_lock else f"{max(front - 1.6, 0):.3f}"
    result += strip.format(outer_wall / 2, 0.8, width - outer_wall, strip_depth)
    return result + json.loads(end.group(1))


def section_contours(path):
    """Read only OpenSCAD's polygonal SVG sections; restore millimetre Y direction."""
    contours = []
    for element in ET.parse(path).getroot().iter("{http://www.w3.org/2000/svg}path"):
        for subpath in re.findall(r"M(.*?)z", element.attrib["d"], re.S):
            numbers = [float(value) for value in
                       re.findall(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", subpath)]
            if len(numbers) < 6 or len(numbers) % 2 or not all(map(math.isfinite, numbers)):
                raise ValueError("invalid OpenSCAD section polygon")
            points = list(zip(numbers[::2], (-value for value in numbers[1::2])))
            contours.append(points)
    if not contours:
        raise ValueError("OpenSCAD section has no polygon contours")
    return contours


def contour_bounds(points):
    return [(min(p[axis] for p in points), max(p[axis] for p in points)) for axis in (0, 1)]


def storage_box_template(source, rows, columns):
    """Actual box SCAD with hand-derived Safe-Lock engineering fixture values, not Rust execution."""
    match = re.search(r'const STORAGE_BOX_SCAD: &str = r#"(.*?)"#;', source, re.S)
    if not match:
        raise ValueError("storage-box SCAD template missing")
    width, depth = columns * 10.15 + (columns - 1) * 1.2, rows * 10.15 + (rows - 1) * 1.2
    fields = dict(rows=rows, columns=columns, cell_width_mm=10.15, divider_thickness_mm=1.2,
                  pitch_mm=11.35, divider_height_mm=22, wall_thickness_mm=2,
                  interior_width_mm=width, interior_depth_mm=depth, interior_height_mm=45,
                  exterior_width_mm=width + 4, exterior_depth_mm=depth + 4, exterior_height_mm=47,
                  lid_inner_width_mm=width + 4.6, lid_inner_depth_mm=depth + 4.6,
                  lid_exterior_width_mm=width + 8.6, lid_exterior_depth_mm=depth + 8.6,
                  lid_exterior_height_mm=10, lid_skirt_height_mm=8)
    return "".join(f"{key} = {value};\n" for key, value in fields.items()) + match.group(1)


class LabwareTemplateTests(unittest.TestCase):
    def test_safe_lock_source_uses_dimensions_without_moving_holes_in_lf_and_crlf(self):
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text(encoding="utf-8")
        selection = source.split("let safe_lock = matches!(", 1)[1].split(");", 1)[0]
        self.assertEqual(set(re.findall(r"RackPhysicalTemplateKind::(\w+)", selection)), {
            "SafeLockTubeRackCollarSeated", "SafeLockTubeRackFloorSeated",
        })
        for rows, columns, pitch in ((4, 6, 20), (1, 1, 20), (3, 8, 14)):
            for opening, edge, height in ((10.45, 4.775, 24), (10.15, 4.925, 44)):
                args = (rows, columns, 2 * edge + opening + (columns - 1) * pitch,
                        10 + edge + opening + (rows - 1) * pitch,
                        height, opening, edge, 10, pitch)
                scad = carrier_template(source, *args, safe_lock=True)
                self.assertEqual(scad, carrier_template(source.replace("\n", "\r\n"),
                                                       *args, safe_lock=True))
                self.assertIn("linear_extrude(height=rack_height)", scad)
                self.assertIn("circle(r=corner_radius)", scad)
                self.assertIn("outer_width - corner_radius", scad)
                self.assertIn("outer_depth - corner_radius", scad)
                self.assertNotIn("cube([outer_width, outer_depth, rack_height]", scad)
                self.assertIn(", front_label_strip_depth, front_label_strip_recess + 0.05]", scad)
                holes = [line for line in scad.splitlines() if "cylinder(" in line]
                self.assertEqual(len(holes), rows * columns)
                for row in range(rows):
                    for column in range(columns):
                        self.assertIn(
                            f"translate([{10 + column * pitch:.3f}, "
                            f"{10 + opening / 2 + row * pitch:.3f}, floor_thickness])",
                            holes[row * columns + column])

    def test_other_carrier_source_bytes_preserve_pre_fix_baselines(self):
        # Synthetic source-fragment controls from 6adee466, recorded before this fix.
        # No committed CAD fixture, Rust execution or manufacturing claim.
        cases = [
            ((4, 6, 62.8, 44.0, 12.0, 9.2, 1.8, 3.0, 10, 1.2),
             "0039700249401982d6f0d10ad0c879d833d2a81d35c7d45b07fface8033e9e36"),
            ((4, 6, 87.2, 60.2, 18.0, 9.2, 4.0, 5.0, 14, 2.4),
             "a6973a768ce96e0383e9142041a682dbd4fb3734f558695f9645ad9c08076117"),
            ((8, 12, 115.1, 78.5, 4.0, 6.9, 4.6, 4.0, 9, 1.0),
             "b7829da79a2ec1fc758d8506e340419ae3fdf48878e37f10b213b341f3b04294"),
            ((8, 12, 127.76, 85.48, 15.0, 6.6, 11.08, 7.94, 9, 2),
             "f904d670a805a9e226722ec3d2c65b670d6f60f95ed3f1e94177116306216b9d"),
        ]
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text(encoding="utf-8")
        for args, expected in cases:
            for text in (source, source.replace("\n", "\r\n")):
                actual = carrier_template(text, *args).encode("utf-8")
                self.assertEqual(hashlib.sha256(actual).hexdigest(), expected)

    def test_safe_lock_mesh_radius_label_depth_floor_and_d1(self):
        binary = openscad_binary()
        if not binary:
            self.skipTest("OpenSCAD unavailable; new Rust exporter execution remains separate")
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text(encoding="utf-8")
        # Both support modes, a minimal grid, compact spacing, and nondefault
        # SCAD parameter controls: fixed literals cannot pass the latter case.
        cases = [
            (4, 6, 20, 10.45, 4.775, 24, 2, 6),
            (4, 6, 20, 10.15, 4.925, 44, 2, 6),
            (1, 1, 20, 10.45, 4.775, 24, 2, 6),
            (3, 8, 14, 10.45, 4.775, 24, 2, 6),
            (4, 6, 20, 10.45, 4.775, 24, 1, 4),
        ]
        with tempfile.TemporaryDirectory(prefix="gentle Safe-Lock geometry ") as directory:
            for index, (rows, columns, pitch, opening, edge, height, radius, label) in enumerate(cases):
                with self.subTest(case=index):
                    width = 2 * edge + opening + (columns - 1) * pitch
                    depth = 10 + edge + opening + (rows - 1) * pitch
                    scad = Path(directory) / f"carrier-{index}.scad"
                    mesh = scad.with_suffix(".off")
                    scad.write_text(carrier_template(
                        source, rows, columns, width, depth, height, opening, edge, 10, pitch,
                        safe_lock=True, label_depth=label, corner_radius=radius), encoding="utf-8")
                    self.assertIn("$fn = 56;", scad.read_text(encoding="utf-8"))
                    result = subprocess.run([binary, "--hardwarnings", "-o", str(mesh), str(scad)],
                                            capture_output=True, text=True, timeout=120)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertNotIn("ERROR:", result.stderr)
                    cuts = {}
                    for name, z in (("floor", 1.999), ("bores", 2.001),
                                    ("label", height - 0.4), ("below_label", height - 0.801)):
                        wrapper = Path(directory) / f"cut-{index}-{name}.scad"
                        svg = wrapper.with_suffix(".svg")
                        wrapper.write_text(
                            f'projection(cut=true) translate([0,0,-{z}]) import("{mesh.name}");\n',
                            encoding="utf-8")
                        result = subprocess.run([binary, "--hardwarnings", "-o", str(svg), str(wrapper)],
                                                capture_output=True, text=True, timeout=120)
                        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                        self.assertNotIn("ERROR:", result.stderr)
                        cuts[name] = section_contours(svg)
                    self.assertEqual({key: len(value) for key, value in cuts.items()}, {
                        "floor": 1, "bores": rows * columns + 1,
                        "label": rows * columns + 2, "below_label": rows * columns + 1,
                    })
                    outline = cuts["floor"][0]
                    for (lo, hi), expected in zip(contour_bounds(outline), (width, depth)):
                        self.assertAlmostEqual(lo, 0, delta=0.001)
                        self.assertAlmostEqual(hi, expected, delta=0.001)
                    for cx, cy in ((radius, radius), (width - radius, radius),
                                   (radius, depth - radius), (width - radius, depth - radius)):
                        arc = [(x, y) for x, y in outline
                               if (x <= radius + 0.001 if cx == radius else x >= width - radius - 0.001)
                               and (y <= radius + 0.001 if cy == radius else y >= depth - radius - 0.001)]
                        self.assertGreater(len(arc), 2)
                        for x, y in arc:
                            # A mesh-plane intersection may fall anywhere on a
                            # 56-gon chord. Bound its apothem, not just vertices;
                            # 0.001 mm separately covers SVG coordinate rounding.
                            distance = math.hypot(x - cx, y - cy)
                            self.assertGreaterEqual(distance, radius * math.cos(math.pi / 56) - 0.001)
                            self.assertLessEqual(distance, radius + 0.001)
                    hole_bounds = [contour_bounds(points) for points in cuts["bores"]
                                   if contour_bounds(points)[0][1] - contour_bounds(points)[0][0] < opening + 0.01]
                    centers = []
                    for bounds in hole_bounds:
                        for lo, hi in bounds:
                            self.assertAlmostEqual(hi - lo, opening, delta=0.001)
                        centers.append(tuple(sum(bound) / 2 for bound in bounds))
                    centers.sort(key=lambda p: (round(p[1], 3), round(p[0], 3)))
                    expected = [(10 + column * pitch, 10 + opening / 2 + row * pitch)
                                for row in range(rows) for column in range(columns)]
                    self.assertEqual(len(centers), len(expected))
                    for actual, target in zip(centers, expected):
                        for value, wanted in zip(actual, target):
                            self.assertAlmostEqual(value, wanted, delta=0.001)
                    pockets = [contour_bounds(points) for points in cuts["label"]
                               if opening + 0.01 < contour_bounds(points)[0][1] - contour_bounds(points)[0][0] < width - 0.01]
                    self.assertEqual(len(pockets), 1)
                    (xlo, xhi), (ylo, yhi) = pockets[0]
                    self.assertAlmostEqual(xlo, float(f"{edge / 2:.3f}"), delta=0.001)
                    self.assertAlmostEqual(xhi - xlo, width - edge, delta=0.001)
                    self.assertAlmostEqual(ylo, 0.8, delta=0.001)
                    self.assertAlmostEqual(yhi - ylo, label, delta=0.001)
                    self.assertLess(yhi, centers[0][1] - opening / 2)

    def test_storage_box_square_cells_and_shared_template_routes_in_lf_and_crlf(self):
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text()
        for rows, columns in ((9, 9), (4, 12), (1, 1)):
            lf = storage_box_template(source, rows, columns)
            crlf = storage_box_template(source.replace("\n", "\r\n"), rows, columns)
            self.assertEqual(lf, crlf.replace("\r\n", "\n"))
            self.assertNotIn("cylinder(", lf)
            self.assertIn("cube([cell_width_mm, cell_width_mm, divider_height_mm + 0.1])", lf)
            self.assertIn("cell_count = rows * columns", lf)
            self.assertIn("Blocked coordinates are still physical square cells", lf)
            for module in ("gentle_storage_box", "gentle_storage_box_insert", "gentle_storage_box_lid"):
                self.assertIn(f"module {module}()", lf)
        for path in ("src/app/rack_workspace_ui.rs", "src/engine_shell.rs", "crates/gentle-protocol/src/lib.rs"):
            self.assertIn("SafeLockStorageBoxDerived", (ROOT / path).read_text())
        types = (ROOT / "crates/gentle-protocol/src/lib.rs").read_text()
        self.assertRegex(types, r'#\[serde\(default, skip_serializing_if = "Option::is_none"\)\]\s*pub storage_box:')

    def test_storage_box_meshes_have_finite_nonempty_geometry_and_expected_box_bounds(self):
        binary = openscad_binary()
        if not binary:
            self.skipTest("OpenSCAD unavailable; Rust exporter/native/fit acceptance stays separate")
        source = (ROOT / "src/engine/state/lineage_containers.rs").read_text()
        with tempfile.TemporaryDirectory(prefix="gentle derived storage boxes ") as directory:
            for rows, columns, width, depth in ((9, 9, 104.95, 104.95), (4, 12, 139, 48.2), (1, 1, 14.15, 14.15)):
                for part in ("all", "box", "lid"):
                    scad = Path(directory) / f"{rows}x{columns}-{part}.scad"
                    mesh = scad.with_suffix(".off")
                    scad.write_text(storage_box_template(source, rows, columns), encoding="utf-8")
                    result = subprocess.run([binary, "-D", f'part="{part}"', "-o", str(mesh), str(scad)],
                                            capture_output=True, text=True, timeout=120)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertNotIn("ERROR:", result.stderr)
                    tokens = mesh.read_text(encoding="utf-8").split()
                    self.assertEqual(tokens[0], "OFF")
                    vertices, faces, _ = map(int, tokens[1:4])
                    self.assertGreater(vertices, 0)
                    self.assertGreater(faces, 0)
                    points = [tuple(map(float, tokens[start:start + 3]))
                              for start in range(4, 4 + 3 * vertices, 3)]
                    self.assertTrue(all(math.isfinite(v) for point in points for v in point))
                    if part != "all":
                        expected = (width, depth, 47) if part == "box" else (width + 4.6, depth + 4.6, 10)
                        for axis, bound in enumerate(expected):
                            self.assertAlmostEqual(min(p[axis] for p in points), 0, places=4)
                            self.assertAlmostEqual(max(p[axis] for p in points), bound, places=3)

    def test_pcr_strip_template_newlines_and_internal_bridges(self):
        source = (ROOT / "crates/gentle-render/src/labware.rs").read_text()
        for count in range(1, 9):
            lf = strip_template(source, count)
            crlf = strip_template(source.replace("\n", "\r\n"), count)
            self.assertEqual(lf, crlf.replace("\r\n", "\n"))
            self.assertIn("if (tube_count > 1)", lf)
            self.assertIn("[0 : tube_count - 2]", lf)
        with self.assertRaises(ValueError):
            strip_template("missing", 8)

    def test_pcr_strip_meshes_have_one_connected_body_and_cut_length_bounds(self):
        binary = openscad_binary()
        if not binary:
            self.skipTest("OpenSCAD unavailable; Rust exporter acceptance remains separate")
        source = (ROOT / "crates/gentle-render/src/labware.rs").read_text()
        with tempfile.TemporaryDirectory(prefix="gentle PCR strips ") as directory:
            for count in range(1, 9):
                with self.subTest(tubes=count):
                    scad = Path(directory) / f"strip-{count}.scad"
                    mesh = scad.with_suffix(".off")
                    scad.write_text(strip_template(source, count), encoding="utf-8")
                    result = subprocess.run([binary, "-o", str(mesh), str(scad)],
                                            capture_output=True, text=True, timeout=120)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertNotIn("ERROR:", result.stderr)
                    mesh_text = mesh.read_text(encoding="utf-8")
                    tokens = mesh_text.split()
                    self.assertEqual(tokens[0], "OFF")
                    vertices, faces, _ = map(int, tokens[1:4])
                    self.assertGreater(vertices, 0)
                    self.assertGreater(faces, 0)
                    points = [tuple(map(float, tokens[start:start + 3]))
                              for start in range(4, 4 + 3 * vertices, 3)]
                    self.assertTrue(all(math.isfinite(v) for point in points for v in point))
                    self.assertAlmostEqual(min(p[0] for p in points), -4, places=4)
                    self.assertAlmostEqual(max(p[0] for p in points), (count - 1) * 9 + 4, places=4)
                    self.assertAlmostEqual(min(p[1] for p in points), -4, places=4)
                    self.assertAlmostEqual(max(p[1] for p in points), 4, places=4)
                    self.assertGreaterEqual(min(p[2] for p in points), -0.001)
                    self.assertLess(min(p[2] for p in points), 0.02)
                    self.assertAlmostEqual(max(p[2] for p in points), 20.75, places=4)
                    # Surface connectivity catches detached tubes or bridges.
                    neighbors = [set() for _ in points]
                    lines = mesh_text.splitlines()
                    header_lines = 1 if len(lines[0].split()) > 1 else 2
                    face_lines = lines[header_lines + vertices:]
                    self.assertEqual(len(face_lines), faces)
                    for line in face_lines:
                        fields = line.split()
                        size = int(fields[0])
                        # OFF face rows may end with RGB/RGBA colors.
                        indices = list(map(int, fields[1:1 + size]))
                        self.assertEqual(len(indices), size)
                        for index in indices[1:]:
                            neighbors[indices[0]].add(index)
                            neighbors[index].add(indices[0])
                    visited, pending = set(), [0]
                    while pending:
                        index = pending.pop()
                        if index not in visited:
                            visited.add(index)
                            pending.extend(neighbors[index] - visited)
                    self.assertEqual(len(visited), vertices, "disconnected strip geometry")
            for count in (0, 9):
                scad = Path(directory) / f"invalid-{count}.scad"
                mesh = scad.with_suffix(".off")
                scad.write_text(strip_template(source, count), encoding="utf-8")
                result = subprocess.run([binary, "-o", str(mesh), str(scad)],
                                        capture_output=True, text=True, timeout=120)
                self.assertIn("ERROR:", result.stderr)
                self.assertFalse(mesh.exists())

    def test_carrier_pitch_glossary_smoke_keeps_the_numeric_example(self):
        source = (ROOT / "src/engine_shell/tests.rs").read_text(encoding="utf-8")
        glossary = json.loads((ROOT / "docs/glossary.json").read_text(encoding="utf-8"))
        command = next(row for row in glossary["commands"]
                       if row["path"] == "racks set-carrier-pitch")
        literal = command["usage"].split()[-1]
        self.assertEqual(literal, "14")
        for text in (source, source.replace("\n", "\r\n")):
            overrides = text.split("fn smoke_command_override(", 1)[1].split(
                "fn smoke_command_line_for_glossary_command(", 1)[0]
            match = re.search(r'"racks set-carrier-pitch"\s*=>\s*Some\("([^"]+)"\)', overrides)
            self.assertIsNotNone(match, "numeric pitch needs a parser-only smoke fixture")
            self.assertEqual(match.group(1), f"racks set-carrier-pitch demo {literal}")
            self.assertIn("fn glossary_carrier_pitch_smoke_retains_numeric_literal_and_clear_flag()", text)

    def test_rack_mutation_admission_matches_dispatch_in_lf_and_crlf_source(self):
        source = (ROOT / "src/engine_shell.rs").read_text(encoding="utf-8")
        cases = (ROOT / "src/engine_shell/tests.rs").read_text(encoding="utf-8").split(
            "fn rack_mutator_command_cases()", 1)[1].split("#[test]", 1)[0]
        expected = set(re.findall(r'"(Racks\w+)"', cases))
        self.assertEqual(len(expected), 11)
        for text in (source, source.replace("\n", "\r\n")):
            with self.subTest(crlf="\r\n" in text):
                classifier = text.split("pub fn is_state_mutating(&self)", 1)[1].split(
                    "/// Render default text help", 1)[0]
                executor = text.split("fn execute_rack_mutation_command(", 1)[1].split(
                    "fn execute_rack_mutation_command_with_expanded_stack(", 1)[0]
                self.assertEqual(set(re.findall(r"ShellCommand::(Racks\w+)", executor)), expected)
                self.assertEqual(set(re.findall(r"Self::(Racks\w+)", classifier)), expected)
                # A deliberate omission must fail the membership comparison.
                incomplete = classifier.replace("Self::RacksSetCarrierPitch", "Self::Missing", 1)
                self.assertNotEqual(set(re.findall(r"Self::(Racks\w+)", incomplete)), expected)

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
