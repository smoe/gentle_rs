"""Public real-image contracts plus opt-in replay with an already built GENtle.

No runtime calculation of fragment sizes is mocked or reimplemented. Synthetic
reports below test only comparison failure handling. Real execution requires an
explicit binary; CI provides its existing native CLI without an extra build.
"""

import copy
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

from scripts import gel_tapestation_acceptance as audit
from scripts import check_tutorial_checkouts as checkout


class TapeStationFixtureTests(unittest.TestCase):
    def test_cross_ladder_permission_does_not_approve_real_subtraction_parameters(self):
        decisions = (audit.ROOT / "docs/decisions.md").read_text(encoding="utf-8")
        policy = decisions.split("### Cross-Ladder Residual Policy\n", 1)[1].split("\n## ", 1)[0]
        self.assertIn("2026-10-25", policy)
        self.assertIn("residuals across different ladders are permitted", policy)
        self.assertIn("coefficients, background and marker masks remain undecided", policy)
        self.assertIn("does not approve", policy)
        for relative in (
            "docs/gel_workbench.md",
            "docs/gel_lane_profile_comparison_plan.md",
            "docs/tutorial/10-03_measured_gel_sizing_gui_cli.md",
            "docs/figures/tapestation_priming/README.md",
            "docs/glen_next_verification_handoff.md",
            "docs/roadmap.md",
        ):
            with self.subTest(path=relative):
                text = (audit.ROOT / relative).read_text(encoding="utf-8")
                self.assertIn("decisions.md#cross-ladder-residual-policy", text)

    def test_priming_plan_keeps_both_cells_roles_ladders_and_no_subtraction(self):
        plan = audit.priming_requests(audit.load_fixture())
        self.assertEqual([g["cell_line"] for g in plan], ["cell_1", "cell_2"])
        self.assertEqual([[(p["lane"], p["ladder"]) for p in g["profiles"]] for g in plan],
                         [[("C1", "L1"), ("E1", "L1"), ("G1", "L2")],
                          [("D1", "L1"), ("F1", "L2"), ("H1", "L2")]])
        for group in plan:
            for profile in group["profiles"]:
                request = profile["request"]
                self.assertEqual(request["action"], "extract_profile")
                self.assertTrue(request["dark_bands"])
                self.assertEqual(request["background"], 0.0)
                self.assertEqual(request["masked_ranges"], [])
            self.assertEqual(len(group["comparisons"]), 3)
            for request in group["comparisons"]:
                self.assertEqual(request["action"], "compare")
                self.assertEqual(request["range"], [100.0, 5000.0])
                self.assertIsNone(request["short_threshold"])
                self.assertNotIn("coefficient", json.dumps(request))
                self.assertNotIn("algebra", json.dumps(request))

    def test_priming_plan_refuses_ambiguous_calibration_and_changed_pairing_bytes(self):
        manifest = audit.load_fixture()
        changed = copy.deepcopy(manifest)
        changed["cases"].append(copy.deepcopy(changed["cases"][1]))
        with self.assertRaisesRegex(ValueError, "one explicit retained calibration"):
            audit.priming_requests(changed)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("region_summaries.json", "SHA256SUMS"):
                (root / name).write_bytes((audit.FIXTURE / name).read_bytes())
            (root / "region_summaries.json").write_bytes(b"{}\n")
            with self.assertRaisesRegex(ValueError, "pairing source digest mismatch"):
                audit.priming_requests(manifest, root)
            data = json.loads((audit.FIXTURE / "region_summaries.json").read_bytes())
            data["after_priming_groups"][0]["custom_lane"] = "../escape"
            (root / "region_summaries.json").write_text(json.dumps(data), encoding="utf-8")
            (root / "SHA256SUMS").write_text(
                f'{audit.sha256(root / "region_summaries.json")}  region_summaries.json\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Invalid priming lane identity"):
                audit.priming_requests(manifest, root)

    def test_retained_priming_return_figures_and_tutorial_are_source_bound(self):
        root = audit.ROOT / "docs/figures/tapestation_priming"
        for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            self.assertEqual(audit.sha256(audit.fixture_file(root, name)), digest)
        data = json.loads((root / "observations.json").read_bytes())
        self.assertRegex(data["binary_source_sha"], r"^[0-9a-f]{40}$")
        self.assertRegex(data["binary_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(set(data["matching_gel_source_files"]), set(audit.PRIMING_FILES))
        after = audit.load_fixture()["sources"][1]
        self.assertEqual(data["source_image_sha256"], after["image_sha256"])
        self.assertEqual(data["pairing_source_sha256"], audit.sha256(audit.FIXTURE / "region_summaries.json"))
        self.assertEqual(data["range_bp"], [100, 5000])
        prose = (audit.ROOT / "docs/tutorial/10-03_measured_gel_sizing_gui_cli.md").read_text(encoding="utf-8")
        plan = audit.priming_requests(audit.load_fixture())
        for group, expected in zip(data["cells"], plan):
            self.assertEqual(group["profiles"], [{k: v for k, v in p.items() if k != "request"}
                                              | {"profile_id": p["request"]["profile_id"]}
                                              for p in expected["profiles"]])
            for report, request, label in zip(group["descriptive_results"], expected["comparisons"],
                    ("Custom", "Custom", "Manufacturer")):
                self.assertEqual(report["request"], {k: v for k, v in request.items() if k != "action"})
                self.assertIsNone(report["p_value"])
                self.assertTrue(report["excluded_intervals"])
                right = "manufacturer" if request["right_profile"].find("manufacturer") >= 0 else "combined"
                left_lane = report["request"]["left_profile"].rsplit("_", 1)[1]
                right_lane = report["request"]["right_profile"].rsplit("_", 1)[1]
                row = f'| {group["cell_line"][-1]} | {label} {left_lane} / {right} {right_lane} | {report["ks_distance"]:.3f} | {report["wasserstein_distance"]:.0f} |'
                self.assertIn(row, prose)
            results = group["descriptive_results"]
            self.assertLess(results[2]["ks_distance"], results[1]["ks_distance"])
            self.assertLess(results[2]["wasserstein_distance"], results[1]["wasserstein_distance"])
            for role in ("custom", "manufacturer"):
                name = f'{group["cell_line"]}_{role}_vs_combined.png'
                png = (root / name).read_bytes()
                self.assertEqual(png[:8], b"\x89PNG\r\n\x1a\n")
                self.assertEqual((int.from_bytes(png[16:20], "big"), int.from_bytes(png[20:24], "big")), (1100, 1460))
                self.assertIn(f"../figures/tapestation_priming/{name}", prose)
        for wording in ("--compare-priming", "No subtraction, coefficient", "not screenshots",
                        "not four percent of DNA molecules", "not between"):
            self.assertIn(wording, prose)

    def test_optional_algebra_lesson_is_explicitly_synthetic_and_fail_closed(self):
        prose = (audit.ROOT / "docs/tutorial/10-03_measured_gel_sizing_gui_cli.md").read_text(encoding="utf-8")
        lesson = prose.split("<!-- gel-algebra-lesson:start -->")[1].split("<!-- gel-algebra-lesson:end -->")[0]
        requests = []
        prefix = "gentle_cli --state /tmp/gentle-algebra-lesson/state.json gel-image workbench '"
        for line in lesson.splitlines():
            if line.startswith("gentle_cli "):
                self.assertTrue(line.startswith(prefix))
                self.assertTrue(line.endswith("'"))
                requests.append(json.loads(line[len(prefix):-1]))
        self.assertEqual([r["action"] for r in requests],
                         ["import_profile", "import_profile", "algebra", "inspect_algebra", "export", "export"])
        digest = hashlib.sha256(b"GENtle synthetic gel algebra v1").hexdigest()
        for request in requests[:2]:
            self.assertEqual(request["profile"]["evidence_kind"], "synthetic")
            self.assertEqual(request["profile"]["source_sha256"], digest)
        self.assertEqual(requests[0]["profile"]["weights"], [2, None, 6])
        self.assertEqual(requests[1]["profile"]["weights"], [3, 1, 3])
        self.assertEqual([t["coefficient"] for t in requests[2]["terms"]], [1, -1])
        for request in requests[4:]:
            self.assertEqual(request["algebra_ids"], ["teaching_residual"])
            self.assertNotIn("allow_unverified", request)
        for wording in ("not an approved normalization", "no biological coefficients", "C1/E1 on L1 and G1 on L2", "D1 on L1 and F1/H1 on L2"):
            self.assertIn(wording, prose)

    def test_source_identity_uses_git_blobs_and_rejects_changed_code_in_both_checkout_modes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "source"
            root.mkdir()
            checkout.git(root, "init", "--quiet")
            checkout.git(root, "config", "core.autocrlf", "false")
            checkout.git(root, "config", "user.name", "Synthetic source identity test")
            checkout.git(root, "config", "user.email", "fixture@example.invalid")
            for name in audit.PRIMING_FILES:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"//! Synthetic source identity fixture.\npub fn example() {}\n")
            checkout.git(root, "add", "--all")
            checkout.git(root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "synthetic source")
            revision = checkout.git(root, "rev-parse", "HEAD").decode("ascii").strip()
            expected = audit.verify_source_files(revision, root)
            expected_priming = audit.verify_source_files(revision, root, files=audit.PRIMING_FILES)
            for mode in checkout.MODES:
                clone = root.parent / mode[0]
                checkout.prepare_checkout(root, clone, mode)
                self.assertEqual(audit.verify_source_files(revision, clone), expected)
                self.assertEqual(audit.verify_source_files(revision, clone, files=audit.PRIMING_FILES), expected_priming)
                path = clone / audit.ALGORITHM_FILES[0]
                self.assertEqual(b"\r\n" in path.read_bytes(), mode[0] == "crlf")
                path.write_bytes(path.read_bytes() + b"// changed\n")
                with self.assertRaises(subprocess.CalledProcessError):
                    audit.verify_source_files(revision, clone)

    def test_sources_geometry_and_semantic_boundaries_are_explicit(self):
        fixture = audit.load_fixture()
        self.assertIn("approved public", fixture["publication_permission"])
        self.assertIn("annotation-assisted", fixture["marking_method"])
        self.assertEqual(fixture["nominal_sizing_range_bp"], [100, 5000])
        self.assertEqual(len(fixture["sources"]), 2)
        self.assertEqual(len(fixture["cases"]), 3)
        sources = {row["id"]: row for row in fixture["sources"]}
        self.assertEqual(sources["before"]["region_window_bp"], [400, 5000])
        self.assertEqual(sources["after"]["region_window_bp"], [500, 7000])
        self.assertTrue(all(row["warnings"] for row in sources.values()))
        self.assertEqual([row["id"] for row in sources["after"]["lanes"]],
                         ["L1", "C1", "D1", "E1", "L2", "F1", "G1", "H1"])
        for case in fixture["cases"]:
            request = audit.request_for_case(fixture, case)
            with self.subTest(case=case["id"]):
                self.assertEqual(request["image_sha256"], sources[case["source"]]["image_sha256"])
                self.assertEqual(request["ladder"]["lane_id"], case["calibration_lane"])
                self.assertEqual(len(request["ladder"]["bands"]), 9)
                self.assertTrue(all(100 <= band["size"] <= 5000 for band in request["ladder"]["bands"]))
                self.assertNotIn("vendor_peak_bp", json.dumps(request))
                self.assertNotIn("Average Size", json.dumps(request))
                self.assertEqual(len(request["sample_bands"]), len(case["bands"]))
                self.assertTrue(all(row["center"]["x"] >= 86 for row in request["sample_bands"]))
        self.assertEqual(sum(band.get("vendor_peak_bp") is not None
                             for case in fixture["cases"] for band in case["bands"]
                             if band.get("role", "sample_peak") == "sample_peak"), 13)

    def test_mismatch_paths_and_missing_binary_fail_instead_of_skipping(self):
        for name in ("../before.png", "/before.png", "C:\\before.png", ".", "a/b.pdf"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                audit.fixture_file(audit.FIXTURE, name)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for path in audit.FIXTURE.iterdir():
                if path.is_file():
                    (root / path.name).write_bytes(path.read_bytes())
            manifest = audit.load_fixture(root)
            for bad_id in ("../outside", "C:\\outside", "/outside"):
                changed = copy.deepcopy(manifest)
                changed["cases"][0]["id"] = bad_id
                (root / "manifest.json").write_text(json.dumps(changed), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "Invalid or duplicate cases identity"):
                    audit.load_fixture(root)
            (root / "manifest.json").write_bytes((audit.FIXTURE / "manifest.json").read_bytes())
            (root / "before.png").write_bytes(b"corrupted image")
            with self.assertRaisesRegex(ValueError, "image digest mismatch"):
                audit.load_fixture(root)
            result = subprocess.run([sys.executable, str(audit.ROOT / "scripts/gel_tapestation_acceptance.py"),
                                     "--gentle-cli", str(root / "missing"), "--output-dir", str(root / "audit")],
                                    text=True, capture_output=True, timeout=15)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Supply an existing executable", result.stderr)
            self.assertFalse((root / "audit").exists())

    def test_comparison_does_not_turn_missing_or_invalid_evidence_into_a_pass(self):
        # Literal synthetic result values test the comparison, not GENtle sizing.
        case = {"references": [[1500, 10], [100, 50]], "bands": [
            {"id": "known", "lane": "A", "vendor_peak_bp": 1000},
            {"id": "unknown", "lane": "A", "vendor_peak_bp": None},
            {"id": "marker", "lane": "A", "expected_status": "outside_calibrated_range"},
        ]}
        report = {"estimates": [
            {"band_id": "known", "lane_id": "A", "center": {"y": 30}, "status": "interpolated", "estimated_size": 1000, "reference_band_ids": ["a", "b"]},
            {"band_id": "unknown", "lane_id": "A", "center": {"y": 35}, "status": "interpolated", "estimated_size": 900, "reference_band_ids": ["a", "b"]},
            {"band_id": "marker", "lane_id": "A", "center": {"y": 60}, "status": "outside_calibrated_range", "estimated_size": None, "reference_band_ids": []},
        ]}
        for band, y in zip(case["bands"], (30, 35, 60)):
            band["y"] = y
        result = audit.compare_report(report, case, 0.05)
        self.assertTrue(all(row["check_pass"] for row in result))
        self.assertEqual(result[1]["comparison"], "vendor_peak_unavailable")
        self.assertIsNone(result[1]["relative_error"])
        self.assertIsNone(result[1]["agreement_pass"])
        for value in (float("nan"), None, -1, 1100, True):
            changed = copy.deepcopy(report)
            changed["estimates"][0]["estimated_size"] = value
            self.assertFalse(audit.compare_report(changed, case, 0.05)[0]["check_pass"])
        changed = copy.deepcopy(report)
        changed["estimates"][2]["estimated_size"] = 15
        self.assertFalse(audit.compare_report(changed, case, 0.05)[2]["check_pass"])
        changed = copy.deepcopy(report)
        changed["estimates"][1]["estimated_size"] = 10000
        self.assertFalse(audit.compare_report(changed, case, 0.05)[1]["check_pass"])
        changed["estimates"][1]["lane_id"] = "other_ladder"
        with self.assertRaisesRegex(ValueError, "changed the retained lane"):
            audit.compare_report(changed, case, 0.05)
        for values in (report["estimates"][:-1], report["estimates"] + report["estimates"][:1]):
            with self.assertRaisesRegex(ValueError, "omitted, duplicated or invented"):
                audit.compare_report({"estimates": values}, case, 0.05)

    def test_native_ci_uses_existing_cli_without_an_extra_build(self):
        workflow = (audit.ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
        for job, suffix, python in (("macos", "", "python3"), ("linux", "", "python3"), ("windows", ".exe", "python")):
            block = re.split(r"\n  [a-z][a-z0-9_-]*:\n", workflow.split(f"\n  {job}:\n", 1)[1], maxsplit=1)[0]
            self.assertIn(f"GENTLE_GEL_TEST_CLI: ${{{{ github.workspace }}}}/target/debug/gentle_cli{suffix}", block)
            self.assertIn(f"{python} -m unittest scripts.test_gel_tapestation.TapeStationBinaryTests -v", block)

    def test_real_tutorial_pairs_and_vendor_regions_do_not_claim_full_distributions(self):
        summary_path = audit.FIXTURE / "region_summaries.json"
        data = json.loads(summary_path.read_bytes())
        self.assertEqual(data["schema"], "gentle.fixture.tapestation_region_summaries.v1")
        self.assertIn("not GENtle measurements", data["origin"])
        for source in audit.load_fixture()["sources"]:
            self.assertEqual(data["source_pdf_sha256"][source["id"]], source["pdf_sha256"])
        checksums = dict(line.split("  ", 1)[::-1] for line in
                         (audit.FIXTURE / "SHA256SUMS").read_text().splitlines())
        self.assertEqual(audit.sha256(summary_path), checksums[summary_path.name])
        rows = {(row["source"], row["lane"]): row for row in data["rows"]}
        self.assertEqual(len(rows), 10)
        self.assertTrue(all(row["pdf_pages"] and 0 <= row["percent_of_total"] <= 100
                            for row in rows.values()))
        for pair in data["cleanup_pairs"]:
            before, after = (rows[source, pair[key]] for source, key in
                             (("before", "before_lane"), ("after", "after_lane")))
            self.assertEqual((before["cell_line"], before["priming"]),
                             (after["cell_line"], after["priming"]))
            self.assertEqual((before["from_bp"], before["to_bp"]), (400, 5000))
            self.assertEqual((after["from_bp"], after["to_bp"]), (500, 7000))
        self.assertEqual([(row["custom_lane"], row["manufacturer_lane"], row["combined_lane"])
                          for row in data["after_priming_groups"]],
                         [("C1", "E1", "G1"), ("D1", "F1", "H1")])
        for group in data["after_priming_groups"]:
            for key, priming in (("custom_lane", "custom"), ("manufacturer_lane", "manufacturer_crta"),
                                 ("combined_lane", "custom_and_manufacturer")):
                self.assertEqual((rows["after", group[key]]["cell_line"], rows["after", group[key]]["priming"]),
                                 (group["cell_line"], priming))
        for absent in ("full_distributions", "quantitative_cleanup_effect", "custom_enrichment_from_lane_subtraction"):
            self.assertIsNone(data[absent])
        prose = (audit.ROOT / "docs/tutorial/10-03_measured_gel_sizing_gui_cli.md").read_text(encoding="utf-8")
        for cell in ("cell_1", "cell_2"):
            for priming, label in (("custom", "Custom"), ("manufacturer_crta", "Manufacturer"),
                                   ("custom_and_manufacturer", "Combined")):
                cells = [cell[-1], label]
                for source in ("before", "after"):
                    matches = [row for row in rows.values() if
                               (row["source"], row["cell_line"], row["priming"]) == (source, cell, priming)]
                    cells.extend([str(matches[0]["average_size_bp"]), f'{matches[0]["percent_of_total"]:.2f}']
                                 if matches else ["not supplied", "not supplied"])
                self.assertIn("| " + " | ".join(cells) + " |", prose)
        for wording in ("Real TapeStation Companion", "not matched-window effect estimates",
                        "agreement_pass: null", "does **not** compare full measured lane distributions",
                        "explicitly reviewed compatible profile inputs", "encrypted measurement/model"):
            self.assertIn(wording, prose)


@unittest.skipUnless(os.environ.get("GENTLE_GEL_TEST_CLI"), "Set GENTLE_GEL_TEST_CLI to an existing CLI; no build performed")
class TapeStationBinaryTests(unittest.TestCase):
    def test_real_images_replay_through_cli_persistence_shell_and_exports(self):
        with tempfile.TemporaryDirectory(prefix="gentle tapestation test ") as directory:
            receipt, rows, reports = audit.run_audit(os.environ["GENTLE_GEL_TEST_CLI"], Path(directory) / "evidence")
            self.assertTrue(receipt["all_required_checks_pass"], rows)
            self.assertTrue(all(row["agreement_pass"] is None for row in rows
                                if row["comparison"] == "vendor_peak_unavailable"))
            self.assertEqual(receipt["sample_peak_comparisons"], 13)
            self.assertEqual(receipt["cross_ladder_comparisons"], 14)
            self.assertEqual(receipt["vendor_peak_unavailable"], 15)
            self.assertEqual(receipt["outside_calibration_controls"], 3)
            self.assertEqual(set(reports), {"before_A1", "after_L1", "after_L2"})
            self.assertTrue(all(reports[name]["algorithm"] == "piecewise_log10_size_v1" for name in reports))
            for name in reports:
                self.assertIn(f"{name}.svg", receipt["artifacts"])
                self.assertIn(f"{name}.tsv", receipt["artifacts"])


if __name__ == "__main__":
    unittest.main()
