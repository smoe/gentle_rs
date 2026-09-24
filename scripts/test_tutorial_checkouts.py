"""Synthetic temporary Git fixtures for the cross-host tutorial checkout gate.

Fixtures are created below from literal bytes, committed only inside temporary
repositories, and consumed by the checkout/replay tests. The byte-preservation
regressions also copy the generated capability matrix and existing
provenance-documented TSS and probe fixtures.
No network, build, real checkout changes, or platform-specific tools are used.
"""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

from scripts import check_tutorial_checkouts as checker


class TutorialCheckoutTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="gentle-checkout-test-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "source"
        self.root.mkdir()
        checker.git(self.root, "init", "--quiet")
        checker.git(self.root, "config", "core.autocrlf", "false")
        checker.git(self.root, "config", "user.name", "Synthetic test")
        checker.git(self.root, "config", "user.email", "fixture@example.invalid")
        self.payload = b'{\n  "schema": "synthetic.fixture.v1"\n}\n'
        (self.root / "evidence.json").write_bytes(self.payload)
        (self.root / "unbound.txt").write_bytes(b"first\nsecond\n")
        (self.root / "binary.bin").write_bytes(b"\x00\xff\n\r\n")
        (self.root / ".gitattributes").write_bytes(b"evidence.json text eol=lf\n")
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "fixture")

    def test_both_modes_preserve_bound_bytes_binary_and_history(self):
        revision = checker.git(self.root, "rev-parse", "HEAD").decode().strip()
        (self.root / "untracked-private.txt").write_text("must not be copied")
        before = checker.git(self.root, "status", "--porcelain")
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / mode[0]
                self.assertEqual(checker.prepare_checkout(self.root, target, mode), revision)
                self.assertEqual((target / "evidence.json").read_bytes(), self.payload)
                self.assertEqual((target / "binary.bin").read_bytes(), b"\x00\xff\n\r\n")
                separator = b"\r\n" if mode[0] == "crlf" else b"\n"
                self.assertEqual((target / "unbound.txt").read_bytes(),
                                 separator.join((b"first", b"second", b"")))
                self.assertFalse((target / "untracked-private.txt").exists())
                self.assertEqual(checker.git(target, "log", "-1", "--format=%H", "--", "evidence.json")
                                 .decode().strip(), revision)
        self.assertEqual(checker.git(self.root, "status", "--porcelain"), before)

    def test_missing_lf_rule_changes_digest_and_explicit_overlay_restores_it(self):
        (self.root / ".gitattributes").write_bytes(b"# No protection\n")
        checker.git(self.root, "add", "--", ".gitattributes")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "regression")
        broken = Path(self.tmp.name) / "broken"
        checker.prepare_checkout(self.root, broken, checker.MODES[1])
        self.assertNotEqual(hashlib.sha256((broken / "evidence.json").read_bytes()).digest(),
                            hashlib.sha256(self.payload).digest())
        fixed = Path(self.tmp.name) / "fixed"
        checker.prepare_checkout(self.root, fixed, checker.MODES[1],
                                 b"evidence.json text eol=lf\n")
        self.assertEqual((fixed / "evidence.json").read_bytes(), self.payload)
        self.assertEqual((self.root / ".gitattributes").read_bytes(), b"# No protection\n")

    def test_generated_tutorial_json_stays_byte_exact_in_both_checkout_modes(self):
        # Use the real checkout policy with synthetic generated JSON, so removing
        # a targeted LF rule reproduces Windows' strict drift-check failure.
        (self.root / ".gitattributes").write_bytes(
            (checker.ROOT / ".gitattributes").read_bytes())
        paths = ("docs/tutorial/catalog.json", "docs/tutorial/manifest.json")
        for relative in paths:
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(self.payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "generated JSON")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"generated-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative in paths:
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), self.payload)

    def test_gene_assay_gui_evidence_hashes_survive_both_checkout_modes(self):
        (self.root / ".gitattributes").write_bytes(
            (checker.ROOT / ".gitattributes").read_bytes())
        evidence_dir = Path("docs/screenshots/gene_assay_study_gui")
        source_dir = checker.ROOT / evidence_dir
        evidence = json.loads((source_dir / "evidence.json").read_bytes())
        expected = {}
        for row in evidence["captures"]:
            for field in ("raw_png", "semantic_snapshot", "context_svg"):
                relative = evidence_dir / row[field]
                payload = (checker.ROOT / relative).read_bytes()
                self.assertEqual(hashlib.sha256(payload).hexdigest(), row[f"{field}_sha256"])
                expected[relative] = (payload, row[f"{field}_sha256"])
                destination = self.root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(payload)
        self.assertEqual(len(expected), 18)
        checker.git(self.root, "add", "--all")
        checker.git(
            self.root,
            "-c",
            "commit.gpgsign=false",
            "commit",
            "--quiet",
            "-m",
            "public GUI evidence",
        )
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if b"docs/screenshots/gene_assay_study_gui/" not in line
        )
        broken = Path(self.tmp.name) / "gene-assay-evidence-unprotected"
        checker.prepare_checkout(self.root, broken, checker.MODES[1], unprotected)
        changed = [
            relative
            for relative, (payload, _digest) in expected.items()
            if (broken / relative).read_bytes() != payload
        ]
        self.assertEqual(
            sorted(changed),
            sorted(relative for relative in expected if relative.suffix in (".json", ".svg")),
        )
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"gene-assay-evidence-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    checked_out = (target / relative).read_bytes()
                    self.assertEqual(checked_out, payload)
                    self.assertEqual(hashlib.sha256(checked_out).hexdigest(), digest)

    def test_tss_regulatory_gui_evidence_hashes_survive_both_checkout_modes(self):
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        evidence = json.loads((checker.ROOT /
            "docs/screenshots/tss_regulatory_view_gui/evidence.json").read_bytes())
        records = list(evidence["inputs"].values())
        records.append(evidence["agent_assistant"]["checkpoint"])
        for capture in evidence["captures"]:
            records.extend(capture[field] for field in ("raw_png", "semantic_snapshot"))
        expected = {}
        for record in records:
            relative = Path(record["path"])
            payload = (checker.ROOT / relative).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), record["sha256"])
            expected[relative] = (payload, record["sha256"])
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        self.assertEqual(len(expected), 8)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "retained synthetic TSS GUI evidence")

        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if b"docs/screenshots/tss_regulatory_view_gui/" not in line
        )
        broken = Path(self.tmp.name) / "tss-evidence-unprotected"
        checker.prepare_checkout(self.root, broken, checker.MODES[1], unprotected)
        self.assertEqual(
            sorted(relative for relative, (payload, _) in expected.items()
                   if (broken / relative).read_bytes() != payload),
            sorted(Path(capture["semantic_snapshot"]["path"])
                   for capture in evidence["captures"]),
        )
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"tss-evidence-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, (payload, digest) in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    checked_out = (target / relative).read_bytes()
                    self.assertEqual(checked_out, payload,
                                     f"{relative} needs a scoped .gitattributes LF rule")
                    self.assertEqual(hashlib.sha256(checked_out).hexdigest(), digest)

    def test_generated_parity_matrix_stays_byte_exact_in_both_checkout_modes(self):
        relative = "docs/gui_cli_mcp_parity.md"
        payload = (checker.ROOT / relative).read_bytes()
        self.assertIn(b"\n", payload)
        self.assertFalse(b"\r" in payload,
                         f"{relative} requires a scoped .gitattributes text eol=lf rule")
        attributes = (checker.ROOT / ".gitattributes").read_bytes()
        (self.root / ".gitattributes").write_bytes(attributes)
        destination = self.root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet",
                    "-m", "generated capability matrix")

        # Negative control: the missing rule must reproduce CRLF checkout drift.
        unprotected = b"\n".join(
            line for line in attributes.split(b"\n")
            if relative.encode() not in line
        )
        broken = Path(self.tmp.name) / "parity-unprotected"
        checker.prepare_checkout(self.root, broken, checker.MODES[1], unprotected)
        self.assertTrue(
            (broken / relative).read_bytes() == payload.replace(b"\n", b"\r\n"),
            "The missing-rule control must reproduce CRLF conversion of the parity matrix",
        )
        for mode in checker.MODES:
            with self.subTest(mode=mode[0]):
                target = Path(self.tmp.name) / f"parity-{mode[0]}"
                checker.prepare_checkout(self.root, target, mode)
                self.assertTrue(
                    (target / relative).read_bytes() == payload,
                    f"{relative} changed during {mode[0]} checkout; preserve the "
                    "generator's exact bytes with .gitattributes text eol=lf",
                )

    def test_replay_uses_existing_binary_and_forces_offline(self):
        with patch.dict(os.environ, {"GENTLE_TEST_ONLINE": "1"}), \
                patch.object(checker.subprocess, "run") as run:
            checker.check_checkout(Path("existing-binary"), self.root, "crlf", 12)
        self.assertEqual([call.args[0] for call in run.call_args_list],
                         [["existing-binary", "parity-matrix-check"],
                          ["existing-binary", "--check"], ["existing-binary", "tutorial-check"]])
        for call in run.call_args_list:
            self.assertNotIn("GENTLE_TEST_ONLINE", call.kwargs["env"])
            self.assertEqual(call.kwargs["cwd"], self.root)
            self.assertTrue(call.kwargs["check"])
            self.assertEqual(call.kwargs["timeout"], 12)

    def test_bound_tss_and_probe_fixtures_survive_crlf_checkout(self):
        (self.root / ".gitattributes").write_bytes(
            (checker.ROOT / ".gitattributes").read_bytes())
        fixture_dirs = (
            "test_files/fixtures/tss_profiles",
            "test_files/fixtures/probe_region_outputs/clariom_e_mtab_14704_tp73_validation",
        )
        expected = {}
        for directory in fixture_dirs:
            for source in (checker.ROOT / directory).rglob("*"):
                if not source.is_file():
                    continue
                relative = source.relative_to(checker.ROOT)
                expected[relative] = source.read_bytes()
                destination = self.root / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(expected[relative])
        provenance_path = Path(fixture_dirs[1]) / "provenance.json"
        provenance = json.loads(expected[provenance_path])
        inputs = provenance["input_fingerprints"]
        self.assertTrue(inputs, "expected retained adapter input provenance")
        for binding in inputs:
            relative = Path(binding["path"])
            self.assertFalse(relative.is_absolute())
            self.assertNotIn("..", relative.parts)
            expected[relative] = (checker.ROOT / relative).read_bytes()
            destination = self.root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(expected[relative])
        self.assertTrue(expected)
        checker.git(self.root, "add", "--all")
        checker.git(self.root, "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "bound fixtures")
        for mode in checker.MODES:
            target = Path(self.tmp.name) / f"bound-{mode[0]}"
            checker.prepare_checkout(self.root, target, mode)
            for relative, payload in expected.items():
                with self.subTest(mode=mode[0], path=relative):
                    self.assertEqual((target / relative).read_bytes(), payload)
            for binding in inputs:
                with self.subTest(mode=mode[0], bound_input=binding["path"]):
                    self.assertEqual(
                        "sha256:" + hashlib.sha256(
                            (target / binding["path"]).read_bytes()).hexdigest(),
                        binding["sha256"],
                        "Preserve the exact input bytes with a scoped .gitattributes LF rule; "
                        "do not change provenance hashes or normalize input bytes")
            bundle = target / fixture_dirs[0]
            for line in (bundle / "SHA256SUMS").read_text().splitlines():
                digest, relative = line.split(maxsplit=1)
                self.assertEqual(hashlib.sha256((bundle / relative).read_bytes()).hexdigest(), digest)

    def test_failed_check_or_timeout_cannot_be_reported_as_success(self):
        for error in (subprocess.CalledProcessError(1, "validator"),
                      subprocess.TimeoutExpired("validator", 12)):
            for failed_call in (0, 1, 2):
                with self.subTest(error=type(error).__name__, failed_call=failed_call), \
                        patch.object(checker.subprocess, "run",
                                     side_effect=[None] * failed_call + [error]):
                    with self.assertRaisesRegex(RuntimeError, "crlf checkout failed"):
                        checker.check_checkout(Path("existing-binary"), self.root, "crlf", 12)

    def test_main_checks_both_modes_at_one_revision_without_building(self):
        binary = Path(self.tmp.name) / "existing-binary"
        binary.touch()
        arguments = ["--repo-root", str(self.root), "--binary", str(binary)]
        with patch.object(checker, "prepare_checkout") as prepare, \
                patch.object(checker, "check_checkout") as check:
            self.assertEqual(checker.main(arguments), 0)
        self.assertEqual([call.args[2][0] for call in prepare.call_args_list], ["lf", "crlf"])
        self.assertEqual(prepare.call_args_list[0].args[4], prepare.call_args_list[1].args[4])
        self.assertEqual([call.args[2] for call in check.call_args_list], ["lf", "crlf"])
        with patch.object(checker, "prepare_checkout"), \
                patch.object(checker, "check_checkout", side_effect=[None, RuntimeError("crlf failed")]):
            with self.assertRaisesRegex(RuntimeError, "crlf failed"):
                checker.main(arguments)


class RetainedTutorialReceiptTests(unittest.TestCase):
    def test_serpine1_occupancy_fixture_matches_retained_figure_provenance(self):
        root = checker.ROOT
        fixture = root / "test_files/fixtures/genomic_regions/serpine1_offline/synthetic_serpine1_cutrun.bed"
        figure = root / "docs/tutorial/generated/artifacts/portable_genomic_regions_offline/artifacts/portable_genomic_regions.locus.svg"
        sources = [element for element in ET.parse(figure).iter()
                   if element.get("data-gentle-occupancy-source") == "synthetic_serpine1_cutrun_v1"]
        self.assertEqual(len(sources), 1, "expected retained synthetic occupancy lane")
        self.assertEqual("sha256:" + hashlib.sha256(fixture.read_bytes()).hexdigest(),
                         sources[0].get("data-gentle-occupancy-source-sha256"),
                         "Preserve the synthetic BED input's LF checkout rule")

    def test_retained_bed_bytes_match_manifest_digests(self):
        root = checker.ROOT / "docs/tutorial/generated"
        manifests = sorted(root.rglob("*.bed.manifest.json"))
        self.assertTrue(manifests, "expected retained BED receipt coverage")
        for path in manifests:
            with self.subTest(manifest=path.relative_to(root)):
                manifest = json.loads(path.read_bytes())
                self.assertEqual(manifest["schema"], "gentle.genomic_region_bed_manifest.v1")
                name = manifest["bed_file_name"]
                self.assertEqual(Path(name).name, name)
                digest = "sha256:" + hashlib.sha256((path.parent / name).read_bytes()).hexdigest()
                self.assertEqual(digest, manifest["bed_sha256"],
                                 "Preserve the retained BED export's LF checkout rule")


if __name__ == "__main__":
    unittest.main()
