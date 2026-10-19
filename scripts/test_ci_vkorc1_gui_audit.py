"""Hand-crafted tests for synthetic base proofs and CI failure retention."""

from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts.ci_vkorc1_gui_audit import pair_evidence, run_logged, write_generation_ledger


class PreparedOwnerCommandTests(unittest.TestCase):
    def active_section(self):
        path = Path(__file__).resolve().parents[1] / "docs/glen_audit_followups_20261008.md"
        return path.read_text(encoding="utf-8").split("## Prepared VKORC1 Next Actions", 1)[1].split(
            "## Owner-Requested Main Integration", 1)[0]

    def test_active_guidance_records_platform_passes_without_selecting_a_candidate(self):
        section = self.active_section()
        command = re.search(r"```bash\n(.*?)```", section, re.S).group(1)
        self.assertNotRegex(command, r"(?m)^AUDIT_SHA=")
        self.assertIn('${AUDIT_SHA:?', command)
        self.assertIn('${NORMAL_RUN_ID:?', command)
        self.assertIn('^[0-9a-f]{40}$', command)
        self.assertIn("38774813595", section)
        self.assertIn("116559658104", section)
        self.assertIn("116559599111", section)
        self.assertIn("39066178748", section)
        self.assertIn("two `gene_set_arrangement` tests on macOS, one on Windows", section)
        self.assertIn("Neither failed run in the table qualifies", section)
        self.assertIn("not an atomic SHA", section)
        self.assertNotIn("Windows and the other four filters remain open", section)
        self.assertIn("single next-check handoff", section)

    @unittest.skipUnless(os.name != "nt" and shutil.which("bash"),
                         "Executes a Unix Bash documentation contract, not native Windows acceptance")
    def test_prepared_bash_guards_fail_closed_without_network_or_actual_dispatch(self):
        command = re.search(r"```bash\n(.*?)```", self.active_section(), re.S).group(1)
        # All gh calls are a shell function. An empty PATH also prevents a real
        # network executable from being used if the documented command drifts.
        stub = r'''
gh() {
  printf '%s\n' "$*" >> "$CALLS_FILE"
  case "$1:$2" in
    api:repos/*/git/ref/heads/main) printf '%s\n' "$FAKE_REMOTE_SHA" ;;
    api:repos/*/actions/runs/*)
      printf '%s\t%s\t%s\t%s\t%s\n' "$FAKE_CI_PATH" "$FAKE_CI_SHA" "$FAKE_CI_STATUS" "$FAKE_CI_CONCLUSION" "$FAKE_CI_EVENT" ;;
    workflow:run) printf '%s\n' 'simulated, not dispatched' > "$DISPATCH_FILE" ;;
    *) return 97 ;;
  esac
}
'''
        sha = "a" * 40
        cases = (
            ("matching successful source", {}, True),
            ("missing SHA", {"AUDIT_SHA": None}, False),
            ("abbreviated SHA", {"AUDIT_SHA": "6112095d"}, False),
            ("invalid SHA", {"AUDIT_SHA": "x" * 40}, False),
            ("missing normal run", {"NORMAL_RUN_ID": None}, False),
            ("invalid normal run", {"NORMAL_RUN_ID": "42;exit 0"}, False),
            ("published ref moved", {"FAKE_REMOTE_SHA": "b" * 40}, False),
            ("different CI source", {"FAKE_CI_SHA": "b" * 40}, False),
            ("failed normal CI", {"FAKE_CI_CONCLUSION": "failure"}, False),
            ("incomplete normal CI", {"FAKE_CI_STATUS": "in_progress"}, False),
            ("wrong workflow", {"FAKE_CI_PATH": ".github/workflows/release.yml"}, False),
            ("not push CI", {"FAKE_CI_EVENT": "workflow_dispatch"}, False),
        )
        for name, overrides, allowed in cases:
            with self.subTest(case=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                calls, dispatched = root / "calls", root / "simulated-dispatch"
                env = {"PATH": str(root / "empty-path"), "CALLS_FILE": str(calls),
                       "DISPATCH_FILE": str(dispatched), "AUDIT_SHA": sha, "NORMAL_RUN_ID": "42",
                       "FAKE_REMOTE_SHA": sha, "FAKE_CI_SHA": sha,
                       "FAKE_CI_PATH": ".github/workflows/ci.yml", "FAKE_CI_STATUS": "completed",
                       "FAKE_CI_CONCLUSION": "success", "FAKE_CI_EVENT": "push"}
                env.update({key: value for key, value in overrides.items() if value is not None})
                for key, value in overrides.items():
                    if value is None:
                        env.pop(key, None)
                result = subprocess.run([shutil.which("bash"), "--noprofile", "--norc", "-c", stub + command],
                                        env=env, text=True, capture_output=True, timeout=10)
                self.assertEqual(result.returncode == 0, allowed, result.stderr)
                self.assertEqual(dispatched.exists(), allowed)
                invoked = calls.read_text() if calls.exists() else ""
                self.assertEqual("workflow run" in invoked, allowed)
                if allowed:
                    self.assertIn("--ref main -f platform=linux -f audit_followups=true -f audit_public_vkorc1=true", invoked)
                if name in ("missing SHA", "abbreviated SHA", "invalid SHA", "missing normal run", "invalid normal run"):
                    self.assertEqual(invoked, "", "Invalid input must stop before any GitHub API access")


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

    def test_native_artifact_keeps_its_hash_bound_isolated_gui_settings(self):
        workflow = (Path(__file__).resolve().parents[1] / ".github/workflows/ci.yml").read_text()
        upload = workflow.split("- name: Preserve candidate-bound evidence, including failures", 1)[1].split("- name: Prepare fresh public 08.04 inputs", 1)[0]
        self.assertIn("path: ${{ runner.temp }}/vkorc1-gui-audit/", upload)
        self.assertIn("include-hidden-files: true", upload)
        self.assertNotIn("path: ~/", upload)

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
