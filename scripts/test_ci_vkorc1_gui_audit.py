"""Hand-crafted tests for synthetic base proofs and CI failure retention."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts.ci_vkorc1_gui_audit import pair_evidence, run_logged, write_generation_ledger


class FailureRetentionTests(unittest.TestCase):
    def test_restored_generation_ledger_preserves_unicode_lf_and_checksum_values(self):
        ledger = {"caption": "800\u00d7600", "file_checksums": {"historical.json": "unchanged"}}
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "report.json"
            write_generation_ledger(path, ledger)
            payload = path.read_bytes()
            self.assertIn("800\u00d7600".encode("utf-8"), payload)
            self.assertNotIn(b"\\u00d7", payload)
            self.assertNotIn(b"\r", payload)
            self.assertTrue(payload.endswith(b"\n"))
            self.assertEqual(json.loads(payload), ledger)

    def test_native_audit_installs_dynamically_loaded_x11_keyboard_library(self):
        workflow = (Path(__file__).resolve().parents[1] / ".github/workflows/ci.yml").read_text()
        native_job = workflow.split("  vkorc1-gui-audit:", 1)[1].split("  gui-scroll-audit:", 1)[0]
        self.assertIn("libxkbcommon-x11-dev", native_job)

    def test_native_namespace_drops_privilege_and_retains_before_public_preparation(self):
        root = Path(__file__).resolve().parents[1]
        helper = (root / "scripts/ci_vkorc1_gui_audit.py").read_text()
        self.assertIn('sudo --preserve-env=DISPLAY,XAUTHORITY unshare --net --', helper)
        self.assertIn('setpriv --reuid="$(id -u)" --regid="$(id -g)" --init-groups', helper)
        self.assertLess(helper.index('setpriv --reuid='), helper.index('python3 "$2/scripts/tutorial_gui_acceptance.py"'))
        workflow = (root / ".github/workflows/ci.yml").read_text()
        self.assertLess(workflow.index("- name: Preserve candidate-bound evidence, including failures"),
                        workflow.index("- name: Prepare fresh public 08.04 inputs"))

    def test_vkorc1_companion_declares_required_use_case_context(self):
        root = Path(__file__).resolve().parents[1]
        source = json.loads((root / "docs/tutorial/sources/08-04_vkorc1_warfarin_promoter_luciferase_gui.json").read_bytes())
        use_cases = source["generated_chapter"]["use_cases"]
        self.assertTrue(use_cases)
        self.assertTrue(all(isinstance(value, str) and value.strip() for value in use_cases))
        self.assertIn("synthetic", " ".join(use_cases))
        self.assertIn("genomic-forward T", " ".join(use_cases))
        guide = (root / "docs/tutorial/08-04_vkorc1_warfarin_promoter_luciferase_gui.md").read_text()
        self.assertIn("help variant materialize-allele", guide)
        self.assertIn("does not authorize execution", guide)

    def test_failed_generation_records_exit_and_exposes_retained_diagnostic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            receipt = {"commands": []}
            output = io.StringIO()

            def fail(command, **kwargs):
                kwargs["stdout"].write(b"hand-crafted generation failure\n")
                return subprocess.CompletedProcess(command, 1)

            with patch("scripts.ci_vkorc1_gui_audit.subprocess.run", side_effect=fail), redirect_stdout(output):
                with self.assertRaisesRegex(RuntimeError, "tutorial-generate failed"):
                    run_logged(["unused-test-helper"], "tutorial-generate", root, root, receipt)
            self.assertIn("hand-crafted generation failure", output.getvalue())
            self.assertEqual((root / "tutorial-generate.log").read_bytes(), b"hand-crafted generation failure\n")
            self.assertEqual(receipt["commands"][0]["exit_code"], 1)

    def test_rust_filter_log_names_are_portable_without_changing_filters(self):
        workflow = (Path(__file__).resolve().parents[1] / ".github/workflows/ci.yml").read_text()
        self.assertIn("tutorial_gui_semantics::tests; do", workflow)
        self.assertIn('vkorc1-gui-audit/${filter//:/_}.log', workflow)
        self.assertNotIn('vkorc1-gui-audit/$filter.log', workflow)


class PairEvidenceTests(unittest.TestCase):
    def project(self, source="aaaaaacaaaaaaaaaaaaa", reference="aaaaaaCaaaaaaaaaaaaa",
                alternate="aaaaaaTaaaaaaaaaaaaa"):
        return {"sequences": {"vkorc1_rs9923231_promoter_" + role: {"seq": {"seq": list(bases.encode())}}
                              for role, bases in (("fragment", source), ("reference", reference), ("alternate", alternate))}}

    def test_base_proof_is_explicitly_synthetic_and_not_a_screenshot(self):
        result = pair_evidence(self.project())
        self.assertEqual(result["differences"], [{"position_0based": 6, "reference": "C", "alternate": "T"}])
        self.assertTrue(result["synthetic"])
        self.assertFalse(result["native_screenshot"])
        self.assertFalse(result["online_vkorc1_accepted"])
        self.assertFalse(result["human_scientific_approval"])
        self.assertEqual(result["raw_sequences"]["fragment"], "aaaaaacaaaaaaaaaaaaa")
        self.assertEqual(result["raw_sequences"]["reference"], "aaaaaaCaaaaaaaaaaaaa")
        self.assertEqual(result["sequences"]["reference"], "AAAAAACAAAAAAAAAAAAA")

    def test_changed_source_missing_output_wrong_allele_and_length_fail_closed(self):
        for project in (self.project(source="A" * 20), self.project(alternate="A" * 20),
                        self.project(reference="AaaaaaCaaaaaaaaaaaaa"),
                        self.project(alternate="aaaaaaTaaaaaaaaaaaaaa"), {"sequences": {}}):
            with self.subTest(project=project), self.assertRaises((AssertionError, KeyError)):
                pair_evidence(project)


if __name__ == "__main__":
    unittest.main()
