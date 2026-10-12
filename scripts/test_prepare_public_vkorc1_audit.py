"""Hand-crafted inline tests for the public audit helper's narrow assertions."""

from contextlib import redirect_stdout
import hashlib
import io
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import URLError
from xml.etree import ElementTree

from scripts import prepare_public_vkorc1_audit as public_audit
from scripts.tutorial_acceptance_evidence import checkout_identity
from scripts.prepare_public_vkorc1_audit import (
    assert_reviewed_fragment, pair_content, public_pair_svg, report_content, retained_refsnp_url, PREFIX,
    REFSNP_URL,
)


class PublicAuditIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="public audit identity ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.repo = self.root / "repo with spaces"
        self.repo.mkdir()
        # Disposable committed inputs and a non-executable stand-in are synthetic.
        files = {
            "Cargo.lock": b"synthetic lock\n",
            "docs/tutorial/catalog.json": b"{}\n",
            "docs/tutorial/manifest.json": b"{}\n",
            "assets/genomes.json": b"{}\n",
            public_audit.SOURCE: b"{}\n",
            public_audit.WORKFLOW: json.dumps({"workflow": {"ops": [
                {"PrepareGenome": {}}, {"FetchGenomeSequence": {}},
            ]}}).encode(),
        }
        for name, payload in files.items():
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        for args in (
            ("init", "--quiet"), ("config", "user.name", "Synthetic audit test"),
            ("config", "user.email", "fixture@example.invalid"),
            ("add", "--all"),
            ("-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "synthetic inputs"),
        ):
            subprocess.run(["git", "-C", str(self.repo), *args], check=True, capture_output=True)
        self.revision = subprocess.check_output(
            ["git", "-C", str(self.repo), "rev-parse", "HEAD"], text=True,
        ).strip()
        self.binary = self.root / "synthetic gentle_cli"
        self.binary.write_bytes(b"synthetic non-executable binary fixture\n")
        self.evidence = self.root / "evidence"
        self.version = f"GENtle synthetic\nSource revision synthetic+git.{self.revision}\n"

    def version_probe(self, version=None, mutate=False):
        original_run = subprocess.run

        def run(command, **kwargs):
            if command == [str(self.binary), "--version"]:
                if mutate:
                    self.binary.write_bytes(b"changed during synthetic probe\n")
                output = self.version if version is None else version
                return subprocess.CompletedProcess(
                    command, 0, stdout=output if kwargs.get("text") else output.encode(),
                )
            return original_run(command, **kwargs)

        return patch.object(subprocess, "run", side_effect=run)

    def test_dirty_inputs_refuse_before_binary_probe_network_or_evidence_creation(self):
        for index, name in enumerate((public_audit.SOURCE, "untracked-input.txt")):
            with self.subTest(name=name):
                evidence = self.root / f"dirty-evidence-{index}"
                path = self.repo / name
                before = path.read_bytes() if path.exists() else None
                path.write_bytes(b"uncommitted synthetic input\n")
                try:
                    with self.version_probe(), redirect_stdout(io.StringIO()), \
                            patch.object(public_audit, "binary_identity") as binary, \
                            patch.object(public_audit.urllib.request, "urlopen", side_effect=URLError("synthetic offline refusal")) as network:
                        with self.assertRaisesRegex(ValueError, "clean"):
                            public_audit.audit(self.repo, self.binary, evidence)
                        binary.assert_not_called()
                        network.assert_not_called()
                        self.assertFalse(evidence.exists())
                finally:
                    if before is None:
                        path.unlink()
                    else:
                        path.write_bytes(before)

    def test_missing_stale_or_duplicate_source_refuses_before_fetch(self):
        for index, version in enumerate((
            "GENtle synthetic\n", "Source revision unknown\n",
            f"Source revision synthetic+git.{'0' * 40}\n",
            self.version + self.version,
        )):
            evidence = self.root / f"invalid-evidence-{index}"
            with self.subTest(version=version), self.version_probe(version), \
                    patch.object(public_audit.urllib.request, "urlopen") as network:
                with self.assertRaisesRegex(ValueError, "full candidate revision"):
                    public_audit.audit(self.repo, self.binary, evidence)
                network.assert_not_called()
                self.assertFalse(evidence.exists())

    def test_replaced_binary_during_probe_refuses_before_fetch(self):
        with self.version_probe(mutate=True), \
                patch.object(public_audit.urllib.request, "urlopen") as network:
            with self.assertRaisesRegex(ValueError, "changed during identity probe"):
                public_audit.audit(self.repo, self.binary, self.evidence)
            network.assert_not_called()
            self.assertFalse(self.evidence.exists())

    def test_clean_inputs_are_bound_before_retaining_an_offline_fetch_failure(self):
        with self.version_probe(), redirect_stdout(io.StringIO()), \
                patch.object(public_audit.urllib.request, "urlopen", side_effect=URLError("synthetic offline refusal")) as network:
            self.assertEqual(public_audit.audit(self.repo, self.binary, self.evidence), 1)
        network.assert_called_once()
        receipt = json.loads((self.evidence / "receipt.json").read_bytes())
        self.assertEqual(receipt["status"], "fail")
        self.assertIn("synthetic offline refusal", receipt["message"])
        self.assertEqual(receipt["source_revision"], self.revision)
        self.assertEqual(receipt["binary_version"], self.version)
        self.assertEqual(receipt["binary_sha256"], hashlib.sha256(self.binary.read_bytes()).hexdigest())
        self.assertEqual(receipt["tutorial_source_sha256"], hashlib.sha256((self.repo / public_audit.SOURCE).read_bytes()).hexdigest())
        self.assertFalse(receipt["native_gui_accepted"])
        self.assertFalse(receipt["human_scientific_approval"])

    def test_completion_refuses_changed_source_or_binary(self):
        identity = checkout_identity(self.repo, self.revision)
        binary_hash = public_audit.sha(self.binary)
        public_audit.assert_audit_inputs_unchanged(self.repo, self.binary, identity, binary_hash)
        source = self.repo / public_audit.SOURCE
        before = source.read_bytes()
        source.write_bytes(b"changed after synthetic work\n")
        with self.assertRaisesRegex(ValueError, "clean"):
            public_audit.assert_audit_inputs_unchanged(self.repo, self.binary, identity, binary_hash)
        source.write_bytes(before)
        self.binary.write_bytes(b"changed after synthetic work\n")
        with self.assertRaisesRegex(ValueError, "binary changed during audit"):
            public_audit.assert_audit_inputs_unchanged(self.repo, self.binary, identity, binary_hash)


class PublicAuditAssertionsTests(unittest.TestCase):
    def test_public_endpoint_matches_engine_default_without_obsolete_beta_path(self):
        engine = (Path(__file__).resolve().parents[1] / "src/engine/analysis/feature_expert_ops.rs").read_text()
        self.assertIn(REFSNP_URL.replace("9923231", "{refsnp_id}"), engine)
        self.assertNotIn("/beta/", REFSNP_URL)

    def project(self, reference="aaCaa", alternate="aaTaa", fragment="aacaa"):
        return {"sequences": {PREFIX + role: {"seq": {"seq": list(value.encode())}}
                              for role, value in (("fragment", fragment), ("reference", reference), ("alternate", alternate))}}

    def test_retained_file_override_preserves_raw_spaces_quotes_and_backslashes(self):
        path = Path('/tmp/public response/ref"snp\\record.json')
        self.assertEqual(retained_refsnp_url(path), 'file:///tmp/public response/ref"snp\\record.json')

    def test_engine_recommendation_must_match_reviewed_tutorial_geometry(self):
        report = {"recommended_candidate_id": "chosen", "candidates": [
            {"candidate_id": "chosen", "recommended": True, "start_0based": 10, "end_0based_exclusive": 20}]}
        assert_reviewed_fragment(report, {"ExtractRegion": {"from": 10, "to": 20}})
        with self.assertRaises(AssertionError):
            assert_reviewed_fragment(report, {"ExtractRegion": {"from": 10, "to": 21}})

    def test_only_known_execution_identity_is_excluded_from_report_comparison(self):
        report = {"generated_at_unix_ms": 1, "op_id": "op-1", "run_id": "run",
                  "genomic_alt": "A,G,T", "promoter_windows": [2, 1], "unknown": "retain"}
        self.assertEqual(report_content(report), {"genomic_alt": "A,G,T", "promoter_windows": [2, 1], "unknown": "retain"})

    def test_exact_pair_preserves_raw_nonvariant_bytes_without_scientific_claims(self):
        result = pair_content(self.project())
        self.assertEqual(result["differences"], [{"position_0based": 2, "reference": "C", "alternate": "T"}])
        self.assertEqual(result["raw_sequences"]["fragment"], "aacaa")
        self.assertNotIn("synthetic", result, "Pair shape alone does not establish public provenance")
        self.assertFalse(result["native_screenshot"])
        self.assertFalse(result["human_scientific_approval"])

    def test_wrong_allele_length_second_difference_and_changed_nonvariant_case_fail(self):
        for project in (self.project(alternate="aaGaa"), self.project(alternate="aaTaaa"),
                        self.project(alternate="aaTta"), self.project(reference="AaCaa", alternate="AaTaa")):
            with self.subTest(project=project), self.assertRaises(AssertionError):
                pair_content(project)

    def test_base_view_is_data_bound_and_not_a_native_or_scientific_claim(self):
        proof = pair_content(self.project())
        proof["synthetic"] = False  # Renderer test only; inline DNA remains hand-crafted.
        source, binary = "a" * 40, "b" * 64
        svg = public_pair_svg(proof, source, binary)
        self.assertEqual(ElementTree.fromstring(svg).tag, "{http://www.w3.org/2000/svg}svg")
        self.assertIn(source, svg)
        self.assertIn(binary, svg)
        self.assertIn("differing base 2 (zero-based)", svg)
        self.assertIn("not a native screenshot", svg)
        self.assertIn("No functional, wet-lab, clinical or human scientific approval claim", svg)
        proof["sequence_sha256"]["alternate"] = "0" * 64
        with self.assertRaisesRegex(AssertionError, "base-view proof changed"):
            public_pair_svg(proof, source, binary)

    def test_base_view_refuses_synthetic_or_wrong_revision_bindings(self):
        proof = pair_content(self.project())
        proof["synthetic"] = True
        with self.assertRaises(AssertionError):
            public_pair_svg(proof, "a" * 40, "b" * 64)
        proof["synthetic"] = False
        with self.assertRaises(AssertionError):
            public_pair_svg(proof, "short", "b" * 64)


if __name__ == "__main__":
    unittest.main()
