"""Hand-crafted inline tests for the public audit helper's narrow assertions."""

from contextlib import redirect_stdout
import copy
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


def synthetic_marker_project(alternates=("G", "T")):
    """Hand-crafted saved-project fragment; never loaded as real VKORC1 data."""
    qualifiers = [
        ["label", "rs9923231"], ["db_xref", "dbSNP:rs9923231"],
        ["gentle_generated", "dbsnp_variant_marker"],
        ["dbsnp_reference_check", "match"], ["dbsnp_assembly_check", "match"],
        ["dbsnp_assembly_fallback", "false"], ["dbsnp_variant_class", "snv"],
        ["vcf_ref", "C"], ["vcf_alt", ",".join(alternates)],
        ["assembly_name", "GRCh38.p14"], ["chromosome", "NC_000016.10"],
        ["genomic_position_1based", "31096368"],
    ]
    qualifiers.extend(["dbsnp_spdi", json.dumps({
        "seq_id": "NC_000016.10", "position": 31096367,
        "deleted_sequence": "C", "inserted_sequence": allele,
    })] for allele in ("C", *alternates))
    return {"sequences": {public_audit.CONTEXT: {"seq": {
        "seq": list(b"aaCaa"), "features": [{"kind": "variation", "qualifiers": qualifiers}],
    }}}}


class PublicMarkerEvidenceTests(unittest.TestCase):
    def checked_marker(self, project):
        evidence = public_audit.fetched_marker_evidence(project)
        public_audit.assert_verified_marker(evidence)
        return evidence

    def test_complete_marker_retains_actual_alternates_and_raw_spdi_without_mutation(self):
        for alternates in (("G", "T"), ("T",), ("A", "T")):
            project = synthetic_marker_project(alternates)
            before = copy.deepcopy(project)
            evidence = self.checked_marker(project)
            self.assertEqual(project, before)
            self.assertEqual(evidence["status"], "verified")
            self.assertEqual(evidence["vcf_alt"], list(alternates))
            self.assertEqual(evidence["vcf_alt_raw"], [",".join(alternates)])
            self.assertEqual(evidence["dbsnp_spdi"], [
                value for key, value in project["sequences"][public_audit.CONTEXT]["seq"]["features"][0]["qualifiers"]
                if key == "dbsnp_spdi"])
            self.assertEqual(evidence["chromosome"], "NC_000016.10")
            self.assertEqual(evidence["genomic_position_1based"], 31096368)

    def test_missing_duplicate_and_unsafe_qualifiers_fail_closed(self):
        bad_values = {
            "dbsnp_reference_check": ("mismatch", "unavailable"),
            "dbsnp_assembly_check": ("mismatch", "unavailable"),
            "dbsnp_assembly_fallback": ("true", None),
            "dbsnp_variant_class": ("non_snv",), "vcf_ref": ("G", "CC"),
            "vcf_alt": ("A,G", "T,N", "T,T", "C,T", "", "TT"),
            "assembly_name": ("GRCh37.p13", "notGRCh38", "GRCh380"),
            "chromosome": ("NC_000015.10", "16", "NC_000016.unknown"),
            "genomic_position_1based": ("31096367", "31096369", "invalid"),
        }
        for key, values in bad_values.items():
            for value in (*values, "__missing__", "__duplicate__"):
                with self.subTest(key=key, value=value):
                    project = synthetic_marker_project()
                    qualifiers = project["sequences"][public_audit.CONTEXT]["seq"]["features"][0]["qualifiers"]
                    if value == "__missing__":
                        qualifiers[:] = [row for row in qualifiers if row[0] != key]
                    elif value == "__duplicate__":
                        qualifiers.append(next(row.copy() for row in qualifiers if row[0] == key))
                    else:
                        next(row for row in qualifiers if row[0] == key)[1] = value
                    with self.assertRaisesRegex(ValueError, key):
                        self.checked_marker(project)

    def test_missing_or_ambiguous_marker_and_missing_spdi_are_not_passes(self):
        project = synthetic_marker_project()
        features = project["sequences"][public_audit.CONTEXT]["seq"]["features"]
        features.append(copy.deepcopy(features[0]))
        for invalid in ({"sequences": {}}, project):
            with self.assertRaises(ValueError):
                self.checked_marker(invalid)
        project = synthetic_marker_project()
        qualifiers = project["sequences"][public_audit.CONTEXT]["seq"]["features"][0]["qualifiers"]
        qualifiers[:] = [row for row in qualifiers if row[0] != "dbsnp_spdi"]
        with self.assertRaisesRegex(ValueError, "dbsnp_spdi"):
            self.checked_marker(project)

    def test_spdi_is_retained_but_invalid_or_incompatible_entries_fail(self):
        for value in ("not json", "[]", json.dumps({"seq_id": "NC_000015.10", "position": 31096367}),
                      json.dumps({"seq_id": "NC_000016.10", "position": 31096368})):
            project = synthetic_marker_project()
            qualifiers = project["sequences"][public_audit.CONTEXT]["seq"]["features"][0]["qualifiers"]
            next(row for row in qualifiers if row[0] == "dbsnp_spdi")[1] = value
            evidence = public_audit.fetched_marker_evidence(project)
            self.assertIn(value, evidence["dbsnp_spdi"])
            with self.assertRaisesRegex(ValueError, "dbsnp_spdi"):
                public_audit.assert_verified_marker(evidence)
            self.assertEqual(evidence["status"], "observed_unverified")

    def test_bare_alternate_refusal_is_conditional_and_preserves_persisted_state(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            state, error = root / "state.json", root / "refusal.stderr.txt"
            state.write_bytes(b"hand-crafted unchanged project bytes\n")
            calls = []
            def refuse():
                calls.append("called")
                error.write_text("multiple alternate alleles 'G,T'; choose one explicitly\n")
            single = public_audit.bare_alternate_probe(state, ["T"], refuse, error)
            self.assertEqual(calls, [])
            self.assertFalse(error.exists())
            self.assertEqual(single["status"], "not_applicable_single_alternate")
            self.assertFalse(single["executed"])
            multiple = public_audit.bare_alternate_probe(state, ["G", "T"], refuse, error)
            self.assertEqual(calls, ["called"])
            self.assertEqual(multiple["status"], "refused")
            self.assertEqual(multiple["state_sha256_before"], multiple["state_sha256_after"])
            self.assertEqual(multiple["observed_alternates"], ["G", "T"])

    def test_unrelated_refusal_changed_state_or_wrong_alternates_do_not_pass(self):
        for case in ("changed_state", "wrong_error", "wrong_alternates"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                state, error = root / "state.json", root / "stderr.txt"
                state.write_bytes(b"synthetic original project\n")
                def invoke():
                    if case == "changed_state":
                        state.write_bytes(b"unexpected mutation\n")
                    error.write_text("unrelated error" if case == "wrong_error" else
                                     "multiple alternate alleles 'A,T'; choose one explicitly" if case == "wrong_alternates" else
                                     "multiple alternate alleles 'G,T'; choose one explicitly")
                with self.assertRaises(ValueError):
                    public_audit.bare_alternate_probe(state, ["G", "T"], invoke, error)


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
                {"PrepareGenome": {}}, {"FetchDbSnpRegion": {}},
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

    def version_probe(self, version=None, mutate=False, fetched_project=None):
        original_run = subprocess.run

        def run(command, **kwargs):
            if command == [str(self.binary), "--version"]:
                if mutate:
                    self.binary.write_bytes(b"changed during synthetic probe\n")
                output = self.version if version is None else version
                return subprocess.CompletedProcess(
                    command, 0, stdout=output if kwargs.get("text") else output.encode(),
                )
            if fetched_project is not None and command[0] == str(self.binary):
                if command[3] != "op":
                    raise AssertionError("Unexpected synthetic CLI command before marker validation")
                operation = json.loads(command[4])
                if "FetchDbSnpRegion" in operation:
                    Path(command[2]).write_text(json.dumps(fetched_project), encoding="utf-8")
                elif "PrepareGenome" not in operation:
                    raise AssertionError("Unexpected synthetic operation before marker validation")
                kwargs["stdout"].write(b"{}\n")
                return subprocess.CompletedProcess(command, 0)
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

    def test_saved_marker_is_validated_before_context_or_allele_work(self):
        project = synthetic_marker_project()
        qualifiers = project["sequences"][public_audit.CONTEXT]["seq"]["features"][0]["qualifiers"]
        next(row for row in qualifiers if row[0] == "dbsnp_reference_check")[1] = "mismatch"
        # This stand-in is neither a network response nor an executable GENtle.
        response = io.BytesIO(b'{"refsnp_id":"9923231"}')
        response.url, response.status, response.headers = REFSNP_URL, 200, {"Date": "synthetic"}
        with self.version_probe(fetched_project=project), redirect_stdout(io.StringIO()), \
                patch.object(public_audit.urllib.request, "urlopen", return_value=response):
            self.assertEqual(public_audit.audit(self.repo, self.binary, self.evidence), 1)
        receipt = json.loads((self.evidence / "receipt.json").read_bytes())
        self.assertEqual(receipt["status"], "fail")
        self.assertIn("dbsnp_reference_check", receipt["message"])
        self.assertEqual(receipt["fetched_marker"]["status"], "observed_unverified")
        self.assertEqual(receipt["fetched_marker"]["vcf_alt"], ["G", "T"])
        self.assertEqual(receipt["fetched_marker"]["qualifiers"]["dbsnp_reference_check"], ["mismatch"])
        self.assertEqual([command["label"] for command in receipt["commands"]],
                         ["prepare-complete-reference", "fetch-retained-public-slice"])
        self.assertIn("source.project.gentle.json", receipt["artifacts"])
        self.assertFalse((self.evidence / "direct.project.gentle.json").exists())

    def test_valid_two_alternate_marker_reaches_next_phase_without_three_allele_assumption(self):
        project = synthetic_marker_project()
        response = io.BytesIO(b'{"refsnp_id":"9923231"}')
        response.url, response.status, response.headers = REFSNP_URL, 200, {"Date": "synthetic"}
        with self.version_probe(fetched_project=project), redirect_stdout(io.StringIO()), \
                patch.object(public_audit.urllib.request, "urlopen", return_value=response):
            self.assertEqual(public_audit.audit(self.repo, self.binary, self.evidence), 1)
        receipt = json.loads((self.evidence / "receipt.json").read_bytes())
        self.assertEqual(receipt["fetched_marker"]["status"], "verified")
        self.assertEqual(receipt["fetched_marker"]["vcf_alt"], ["G", "T"])
        self.assertTrue(receipt["fetched_marker"]["dbsnp_spdi"])
        self.assertIn("agent_parity", receipt["message"], "synthetic source intentionally stops the next phase")
        self.assertTrue((self.evidence / "direct.project.gentle.json").exists())
        self.assertTrue((self.evidence / "shared.project.gentle.json").exists())

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
