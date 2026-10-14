"""Source checks for the integrative Synthetic Biology learning route.

These checks validate documentation boundaries and retained visual identities;
they do not execute GENtle or claim native/scientific acceptance.
"""

import json
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/tutorial/sources/03-09_synthetic_biology_design_build_verify.json"
PROSE = ROOT / "docs/tutorial/03-09_synthetic_biology_design_build_verify.md"


class SyntheticBiologyTutorialTests(unittest.TestCase):
    def test_reality_boundary_checks_run_in_ci(self):
        workflow = (ROOT / ".github/workflows/ci.yml").read_text()
        safety = workflow.split("\n  release-policy:\n", 1)[1].split("\n  select-platform:\n", 1)[0]
        self.assertIn("python3 -m unittest scripts.test_synthetic_biology_tutorial -v", safety)

    def test_source_is_manual_route_not_fabricated_executable_contract(self):
        source = json.loads(SOURCE.read_bytes())
        self.assertEqual(source["id"], "synthetic_biology_design_build_verify")
        self.assertEqual(source["catalog"]["type"], "operational_reference")
        self.assertEqual(source["catalog"]["status"], "manual/reference")
        self.assertNotIn("generated_chapter", source)
        self.assertNotIn("gui_acceptance", source)

    def test_links_and_declared_graphics_resolve(self):
        prose = " ".join(PROSE.read_text(encoding="utf-8").split())
        for link in re.findall(r"\]\(([^)#]+)", prose):
            self.assertTrue((PROSE.parent / link).resolve().is_file(), link)
        source = json.loads(SOURCE.read_bytes())
        for graphic in source["graphics"]:
            self.assertTrue((ROOT / graphic["path"]).is_file(), graphic["path"])
            if graphic["kind"] == "screenshot":
                self.assertRegex(graphic["capture_date"], r"^\d{4}-\d{2}-\d{2}$")

    def test_prose_keeps_planning_experiment_and_function_separate(self):
        prose = " ".join(PROSE.read_text(encoding="utf-8").split()).lower()
        required = (
            "do **not** describe one continuous experimental construct",
            "does not prove biological function",
            "not evidence that the wet-lab assembly succeeded",
            "planning, physical execution and biological function are not conflated",
            "synthetic/public inputs",
        )
        for phrase in required:
            self.assertIn(phrase, prose)


if __name__ == "__main__":
    unittest.main()
