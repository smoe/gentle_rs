"""Source/fixture checks for bounded manual tutorial additions.

These tests do not execute GENtle, DuckDB, a GUI or a biological analysis.
The arithmetic oracle and drawn PNG are explicitly synthetic teaching inputs.
"""

import csv
from fractions import Fraction
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
    def test_custom_only_reference_lesson_has_native_density_and_no_real_data_subtraction(self):
        prose = (ROOT / "docs/tutorial/10-03_measured_gel_sizing_gui_cli.md").read_text(encoding="utf-8")
        for newline in ("\n", "\r\n"):
            text = prose.replace("\n", newline)
            lesson = text.split("<!-- gel-algebra-reference-lesson:start -->")[1].split(
                "<!-- gel-algebra-reference-lesson:end -->")[0]
            requests = [json.loads(line.split(" workbench '", 1)[1][:-1])
                        for line in lesson.splitlines() if line.startswith("gentle_cli ")]
            self.assertEqual(len(requests), 4)
            profile = requests[0]["profile"]
            self.assertEqual(profile["evidence_kind"], "synthetic")
            self.assertEqual(profile["source_sha256"], hashlib.sha256(b"GENtle synthetic gel reference panel v1").hexdigest())
            self.assertEqual(profile["edges"], [90, 125, 175, 250, 350])
            density = [None if weight is None else weight / (b - a)
                       for weight, a, b in zip(profile["weights"], profile["edges"], profile["edges"][1:])]
            self.assertEqual(density, [2, None, 2, 0])
            self.assertNotEqual(density[0], profile["weights"][0] / (125 - 100))
            self.assertEqual(requests[1]["reference_profile_id"], profile["id"])
            for export in requests[2:]:
                self.assertEqual(export["algebra_ids"], ["teaching_residual"])
                self.assertEqual(export["algebra_reference_profile_id"], profile["id"])
                self.assertNotIn("comparison_display", export)
            for forbidden in ("display_factor", "reference_scale", "fit_scale"):
                self.assertTrue(all(forbidden not in request for request in requests))
        self.assertIn("positions, not amplitudes", prose)
        self.assertIn("No real C1/D1 subtraction", prose)
        self.assertIn('use `"comparison_display":"signal_density_per_size"`', prose)
        self.assertIn("unreviewed manifest status is retained", prose)

    def test_signed_algebra_overlay_lesson_has_independent_weighted_arithmetic(self):
        prose = (ROOT / "docs/tutorial/10-03_measured_gel_sizing_gui_cli.md").read_text(encoding="utf-8")
        for newline in ("\n", "\r\n"):
            text = prose.replace("\n", newline)
            requests = []
            for name in ("gel-algebra-lesson", "gel-algebra-weighted-lesson"):
                block = text.split(f"<!-- {name}:start -->", 1)[1].split(f"<!-- {name}:end -->", 1)[0]
                for line in block.splitlines():
                    if line.startswith("gentle_cli "):
                        payload = line.split("gel-image workbench '", 1)[1].removesuffix("'")
                        requests.append(json.loads(payload))
            self.assertEqual(len(requests), 8)
            profiles = {r["profile"]["id"]: r["profile"] for r in requests if r["action"] == "import_profile"}
            expected = {
                "teaching_residual": ([Fraction(-1), None, Fraction(3)],
                                      [Fraction("-0.02"), None, Fraction("0.03")],
                                      [Fraction("-0.06"), Fraction("-0.02"), Fraction("-0.03")]),
                "teaching_half_b": ([Fraction("0.5"), None, Fraction("4.5")],
                                    [Fraction("0.01"), None, Fraction("0.045")],
                                    [Fraction("-0.03"), Fraction("-0.01"), Fraction("-0.015")]),
            }
            for request in (r for r in requests if r["action"] == "algebra"):
                terms = request["terms"]
                edges = profiles[terms[0]["profile_id"]]["edges"]
                widths = [Fraction(b - a) for a, b in zip(edges, edges[1:])]
                contributions = []
                for term in terms:
                    p = profiles[term["profile_id"]]
                    self.assertEqual(p["edges"], edges)
                    contributions.append([None if w is None else Fraction(str(w)) * Fraction(str(term["coefficient"]))
                                          for w in p["weights"]])
                weights = [None if any(w is None for w in values) else sum(values)
                           for values in zip(*contributions)]
                density = [None if w is None else w / width for w, width in zip(weights, widths)]
                b_density = [None if w is None else w / width for w, width in zip(contributions[1], widths)]
                self.assertEqual((weights, density, b_density), expected[request["output_id"]])
            self.assertIn("does **not** recover the missing residual", text)
            self.assertIn("same weighted curves appear in the PDF", text)

    @staticmethod
    def editing_teaching_data(prose):
        block = prose.split("<!-- editing-teaching-data:start -->")[1].split(
            "<!-- editing-teaching-data:end -->"
        )[0]
        return json.loads(block.split("```json", 1)[1].split("```", 1)[0])

    def test_editing_draft_keeps_literal_edit_and_reverse_coordinate_expectations(self):
        prose = (ROOT / "docs/tutorial_drafts/base_prime_editing.md").read_text(encoding="utf-8")
        data = self.editing_teaching_data(prose)
        reference = data["reference"]
        self.assertEqual(reference, "ACGTAC")
        self.assertEqual([edit["expected"] for edit in data["edits"]],
                         ["ATGTAC", "ACGGGTAC", "ACAC"])
        for edit in data["edits"]:
            self.assertEqual(reference[edit["start"]:edit["end"]], edit["deleted"])
            self.assertEqual(reference[:edit["start"]] + edit["inserted"] + reference[edit["end"]:],
                             edit["expected"])
        complement = str.maketrans("ACGT", "TGCA")
        self.assertEqual(reference.translate(complement)[::-1], data["reverse_display"])
        substitution = data["edits"][0]
        self.assertEqual(substitution["expected"].translate(complement)[::-1],
                         data["reverse_edited_display"])
        local_start, local_end = data["reverse_target_interval"]
        self.assertEqual([len(reference) - substitution["end"], len(reference) - substitution["start"]],
                         [local_start, local_end])
        self.assertEqual(data["reverse_display"][local_start:local_end], "G")
        self.assertEqual(data["reverse_edited_display"][local_start:local_end], "A")
        self.assertEqual(data["genomic_start_0based"] + substitution["start"], 1001)
        self.assertEqual(data["genomic_start_0based"] + substitution["end"], 1002)
        for boundary in (0, 3, len(reference)):
            self.assertEqual(len(reference) - (len(reference) - boundary), boundary)

    def test_editing_draft_counts_bystanders_without_claiming_outcome_probabilities(self):
        prose = (ROOT / "docs/tutorial_drafts/base_prime_editing.md").read_text(encoding="utf-8")
        toy = self.editing_teaching_data(prose)["toy_bystanders"]
        self.assertEqual(toy["sequence"], "ACCCGTAA")
        self.assertEqual([toy["sequence"][i] for i in toy["editable_positions"]], ["C", "C", "C"])
        self.assertIn(toy["intended_position"], toy["editable_positions"])
        self.assertEqual([i for i in toy["editable_positions"] if i != toy["intended_position"]],
                         toy["expected_bystanders"])
        self.assertEqual(toy["expected_bystanders"], [1, 3])
        self.assertIn("no validated spacer, PAM or", prose)
        self.assertIn("Possible edits are not measured", prose)

    def test_editing_draft_preserves_denominators_and_lf_crlf_semantics(self):
        prose = (ROOT / "docs/tutorial_drafts/base_prime_editing.md").read_text(encoding="utf-8")
        data = self.editing_teaching_data(prose)
        self.assertEqual(self.editing_teaching_data(prose.replace("\n", "\r\n")), data)
        counts = data["counts"]
        categories = counts["categories"]
        self.assertEqual(set(categories), {"intended_only", "intended_plus_bystander", "other_substitutions",
                                           "unexpected_indels", "unedited", "unclassified"})
        self.assertTrue(all(type(value) is int and value >= 0 for value in categories.values()))
        self.assertEqual(sum(categories.values()), counts["retained"])
        self.assertEqual(counts["input"] - counts["quality_excluded"], counts["retained"])
        self.assertEqual(counts["retained"] - categories["unclassified"], counts["callable"])
        self.assertAlmostEqual(100 * categories["intended_only"] / counts["callable"], 23.5294117647)
        self.assertAlmostEqual(100 * categories["intended_only"] / counts["retained"], 22.2222222222)
        self.assertEqual(100 * categories["intended_only"] / counts["input"], 20)
        for disclosure in ("23.53%", "22.22%", "20%", "unavailable, not 0%", "No untreated control",
                           "total intended-edit presence cannot be calculated"):
            self.assertIn(disclosure, prose)

    def test_editing_draft_stays_outside_active_discovery_and_resolves_local_links(self):
        draft = ROOT / "docs/tutorial_drafts/base_prime_editing.md"
        prose = draft.read_text(encoding="utf-8")
        self.assertFalse(draft.is_relative_to(ROOT / "docs/tutorial"), "Help fallback scans that tree")
        for registry in ("catalog.json", "manifest.json", "review_manifest.json"):
            self.assertNotIn(str(draft.relative_to(ROOT)),
                             (ROOT / "docs/tutorial" / registry).read_text(encoding="utf-8"))
        self.assertFalse((ROOT / "docs/tutorial/sources/base_prime_editing.json").exists())
        for link in re.findall(r"\]\(([^)]+)\)", prose):
            if not link.startswith(("https://", "http://", "#")):
                self.assertTrue((draft.parent / link.split("#")[0]).is_file(), link)
        for disclosure in ("Inactive teaching draft", "teaching data, not a GENtle request",
                           "same engine report", "not a GENtle", "RNA 5' -> spacer | scaffold | RTT | PBS",
                           "are **pending**", "No external tool is installed"):
            self.assertIn(disclosure, prose)
        self.assertNotRegex(prose, r"(?m)^(?:gentle_cli |cargo run |editing (?:intent|base-feasibility|prime-feasibility))")

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

    def test_gel_density_lesson_has_independent_unequal_bin_arithmetic(self):
        prose = (ROOT / "docs/tutorial/10-03_measured_gel_sizing_gui_cli.md").read_text(encoding="utf-8")
        lesson = prose.split("<!-- gel-density-lesson:start -->")[1].split("<!-- gel-density-lesson:end -->")[0]
        requests = [json.loads(line.split(" workbench '", 1)[1][:-1])
                    for line in lesson.splitlines() if line.startswith("gentle_cli ")]
        self.assertEqual(len(requests), 6)
        profiles = [r["profile"] for r in requests if r["action"] == "import_profile"]
        self.assertEqual(len(profiles), 2)
        edges = sorted({x for profile in profiles for x in profile["edges"]})
        self.assertEqual(edges, [100, 125, 200, 400])
        widths = [b - a for a, b in zip(edges, edges[1:])]
        self.assertEqual(widths, [25, 75, 200])
        for profile in profiles:
            self.assertEqual(profile["evidence_kind"], "synthetic")
            self.assertEqual(profile["source_sha256"], hashlib.sha256(b"GENtle synthetic gel density v1").hexdigest())
            rebinned = []
            for a, b in zip(edges, edges[1:]):
                matches = [(i, hi - lo) for i, (lo, hi) in
                           enumerate(zip(profile["edges"], profile["edges"][1:])) if lo <= a and b <= hi]
                (i, width), = matches
                rebinned.append(profile["weights"][i] * (b - a) / width)
            self.assertEqual(rebinned, [25, 75, 200])
            self.assertEqual([w / width for w, width in zip(rebinned, widths)], [1, 1, 1])
            self.assertEqual(sum(rebinned), 300)
        exports = [r for r in requests if r["action"] == "export"]
        self.assertEqual([r.get("comparison_display") for r in exports],
                         ["signal_density_per_size", "bin_integral", None])
        self.assertIn("Equal units do not establish matched exposure", prose)
        self.assertIn("unreviewed manifest status is retained", prose)

    def test_tss_codex_review_is_scoped_to_current_prose(self):
        reviews = json.loads((ROOT / "docs/tutorial/review_manifest.json").read_bytes())
        review, = [row for row in reviews["entries"] if row["tutorial_id"] == "tss_regulatory_view_gui"]
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_bytes())
        projected, = [row for row in catalog["entries"] if row["id"] == review["tutorial_id"]]
        prose = (ROOT / "docs/tutorial/08-16_tss_regulatory_view_gui.md").read_text(encoding="utf-8")
        self.assertGreaterEqual(review["codex_reviewed_at"], "2026-10-12")
        self.assertEqual(projected["codex_reviewed_at"], review["codex_reviewed_at"])
        self.assertEqual(projected["review_status"], "codex_reviewed")
        self.assertIn(f'**{review["codex_reviewed_at"]}**', prose)
        self.assertRegex(prose, r"Reviewed source: `[0-9a-f]{40}`")
        self.assertIn("including the optional DuckDB companion", prose)
        self.assertIn("readability, pedagogy and command/source consistency", prose)
        self.assertIn("not Rust, native GUI, model or human scientific acceptance", prose)
        self.assertIsNone(review["human_reviewed_at"])
        self.assertIsNone(review["human_reviewer"])

    def test_canonical_tutorial_sources_have_review_manifest_entries(self):
        manifest = json.loads((ROOT / "docs/tutorial/review_manifest.json").read_bytes())
        identifiers = [entry["tutorial_id"] for entry in manifest["entries"]]
        self.assertEqual(len(identifiers), len(set(identifiers)), "duplicate review identities")
        known = {source["id"] for path in (ROOT / "docs/tutorial/sources").glob("*.json")
                 if "id" in (source := json.loads(path.read_bytes()))}
        self.assertEqual(known - set(identifiers), set(), "missing review entries change report.json warnings")
        by_id = {entry["tutorial_id"]: entry for entry in manifest["entries"]}
        catalog = {entry["id"]: entry for entry in
                   json.loads((ROOT / "docs/tutorial/catalog.json").read_bytes())["entries"]}
        for identifier in ("microarray_evidence_gui_cli", "measured_gel_sizing_gui_cli",
                           "repeat_environment_cohorts"):
            with self.subTest(tutorial=identifier):
                entry = by_id[identifier]
                self.assertEqual(entry["tutorial_kind"], "guided_walkthrough")
                self.assertEqual(entry["tutorial_status"], "active")
                self.assertNotEqual(catalog[identifier]["review_status"], "missing_review_manifest_entry")
                if entry["codex_reviewed_at"] is None and entry["human_reviewed_at"] is None:
                    self.assertEqual(catalog[identifier]["review_status"], "unreviewed")

    def test_candidate_handoff_links_resolve_to_historical_and_topic_owned_sections(self):
        linked_from = (ROOT / "docs/roadmap.md",
                       ROOT / "docs/release_notes/release_notes_v0.1.0-internal.12.md",
                       ROOT / "docs/clawbio_gentle_integration_onepager.md",
                       ROOT / "docs/release.md")
        historical_revision = "7c52082c24fca85817d70956146623eef6f24607"
        historical_url = (f"https://github.com/smoe/gentle_rs/blob/{historical_revision}/"
                          "docs/internal_12_candidate_preparation.md")
        historical = subprocess.check_output(
            ["git", "show", f"{historical_revision}:docs/internal_12_candidate_preparation.md"],
            cwd=ROOT,
        ).decode("utf-8")

        def headings(text):
            return {re.sub(r"[^a-z0-9 -]", "", heading.lower()).replace(" ", "-")
                    for heading in re.findall(r"(?m)^#{1,6} (.+)$", text)}

        expected = {"corrected-ci-evidence", "ubuntu-runner-decision", "compatibility-review-inventory"}
        seen = set()
        for source in linked_from:
            for link in re.findall(r"\]\(([^)]+)\)", source.read_text(encoding="utf-8")):
                relative, _, fragment = link.partition("#")
                if "internal_12_candidate_preparation.md" in link and link.startswith("https://"):
                    self.assertEqual(relative, historical_url, "historical evidence must stay revision-pinned")
                    text = historical
                elif "internal_12_candidate_preparation.md" in link or fragment in expected:
                    target = (source.parent / relative).resolve() if relative else source
                    self.assertTrue(target.is_file(), link)
                    text = target.read_text(encoding="utf-8")
                else:
                    continue
                if fragment:
                    self.assertIn(fragment, headings(text), link)
                    seen.add(fragment)
        self.assertEqual(seen, expected)
        restored = (ROOT / "docs/internal_12_candidate_preparation.md").read_text(encoding="utf-8")
        self.assertIn("Historical handoff", restored)
        self.assertTrue({"corrected-ci-evidence", "ubuntu-runner-decision", "clawbio-compatibility-decision"}
                        <= headings(restored))

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
        self.assertIn("Activation Pending", draft.read_text())
        self.assertIn("no Ensembl", draft.read_text())
        self.assertIn("release record", draft.read_text())
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
