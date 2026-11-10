"""Synthetic fixture-only checks; opt-in real GENtle CLI covers LF/CRLF ingestion. No build."""

import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from scripts import prepare_exon_locus_rehearsal as fixture
from scripts import prepare_read_end_rehearsal as base

BIN = os.environ.get("GENTLE_READ_END_BIN_DIR")


class ExonLocusInputs(unittest.TestCase):
    def test_synthetic_both_strands_and_raw_fixture_binding_graph(self):
        hashes = []
        for ending in ("\n", "\r\n"):
            with self.subTest(ending=repr(ending)), tempfile.TemporaryDirectory() as tmp:
                self._check_fixture(Path(tmp), ending, hashes)
        self.assertNotEqual(hashes[0], hashes[1])

    def _check_fixture(self, temporary, ending, hashes):
        root = temporary.resolve() / "locus's α"
        fixture.prepare(root, ending=ending)
        request = json.loads((root / "sam-request.json").read_bytes())
        annotation = json.loads((root / "annotation.json").read_bytes())
        self.assertEqual({t["genomic_strand"] for t in annotation["transcripts"]}, {"+", "-"})
        for file in [request["sam"], request["dictionary"], request["manifest"], request["annotation_sources"][0]["annotation"]]:
            self.assertEqual(file, base.bound(Path(file["path"])))
            raw = Path(file["path"]).read_bytes()
            self.assertEqual(file["sha256"], "sha256:" + hashlib.sha256(raw).hexdigest())
            self.assertEqual(b"\r\n" in raw, ending == "\r\n")
        raw = (root / "annotation.json").read_bytes()
        hashes.append(hashlib.sha256(raw).hexdigest())
        model = json.loads((root / "models.json").read_bytes())
        self.assertEqual(model["annotation_sha256"], base.bound(root / "annotation.json")["sha256"])
        units = json.loads((root / "units.json").read_bytes())["units"]
        self.assertEqual(len(units), 3)
        self.assertTrue(all(not u["search"]["search_complete"] for u in units))
        with self.assertRaises(FileExistsError):
            fixture.prepare(root)

    def test_invalid_fixture_ending_refuses_before_creating_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "not-created"
            with self.assertRaises(ValueError):
                fixture.prepare(root, ending="\r")
            self.assertFalse(root.exists())


@unittest.skipUnless(BIN, "Set GENTLE_READ_END_BIN_DIR for same-source actual CLI acceptance")
class ExonLocusProcess(unittest.TestCase):
    def test_actual_locus_cli_models_originals_exports_and_no_clobber(self):
        for ending in ("\n", "\r\n"):
            with self.subTest(ending=repr(ending)):
                self._replay(ending)

    def _replay(self, ending):
        binary = Path(BIN) / ("gentle_cli.exe" if os.name == "nt" else "gentle_cli")
        self.assertTrue(binary.is_file())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve() / "locus's α"
            fixture.prepare(root, ending=ending)
            state = root / "not-created-project.json"

            def run(*args, success=True):
                result = subprocess.run([str(binary), "--state", str(state), *args], capture_output=True, timeout=90)
                self.assertEqual(result.returncode == 0, success, result.stderr.decode(errors="replace"))
                return result

            def command(name, value, *, success=True):
                path = root / name
                base.write(path, value)
                return run("rna-reads", "read-ends", "@" + str(path), success=success)

            run("rna-reads", "transcript-structures-from-sam", "@" + str(root / "sam-request.json"))
            build = fixture.bind(root)
            run("rna-reads", "read-ends", "@" + str(root / "locus-build.json"))
            run("rna-reads", "read-ends", "@" + str(root / "locus-build.json"), success=False)
            report = json.loads((root / "locus/artifacts/pileup.json").read_bytes())
            self.assertEqual(len(report["reads"]), 3)
            shared, = [r for r in report["reads"] if r["unit"]["unit_id"] == "synthetic_shared"]
            self.assertTrue(all(a["class"] == "uninformative_overlap" for a in shared["assignments"]))
            self.assertIsNone(shared["molecule_strand"]["strand"])
            self.assertEqual(shared["membership"], "admitted-hypothesis")
            self.assertGreater(report["work"], 0)
            print(f"exon-locus CLI work ({ending!r}): {report['work']}")
            self.assertEqual({r["base"] for e in shared["genomic_basecalls"] for v in e["variants"] for r in v["runs"]}, {"A", "T"})
            binding = base.bound(root / "locus/report.json")
            view = fixture.view()
            inspect = {"action": "locus", "command": {"action": "inspect", "report": binding, "view": view}}
            command("inspect-locus.json", inspect)
            for index, model in enumerate(report["models"]):
                projected = copy.deepcopy(inspect)
                projected["command"]["view"]["projection"] = {"kind": "transcript", "model_key": model["key"]}
                command(f"inspect-model-{index}.json", projected)
            export = {"action": "locus", "command": {"action": "export", "report": binding, "view": view,
                      "formats": ["json", "tsv", "svg", "png", "pdf"], "output_dir": str(root / "exports"),
                      "max_output_bytes": 32 * 1024 * 1024, "cancel_path": None}}
            command("export-locus.json", export)
            run("rna-reads", "read-ends", "@" + str(root / "export-locus.json"), success=False)
            first = {p.name: p.read_bytes() for p in (root / "exports/artifacts").iterdir()}
            export["command"]["output_dir"] = str(root / "repeat")
            command("repeat-locus.json", export)
            for name in ("page.json", "lanes.tsv", "pileup.svg"):
                self.assertEqual(first[name], (root / "repeat/artifacts" / name).read_bytes())
            page = json.loads(first["page.json"])
            self.assertEqual(page["total_input_units"], 3)
            self.assertEqual(len({r["read_key"] for r in page["reads"]}), 3)
            narrow = copy.deepcopy(build)
            narrow["command"]["request"]["region"]["flank_left"] = 0
            narrow["command"]["request"]["output_dir"] = str(root / "narrow")
            command("narrow-build.json", narrow)
            command("narrow-json-export.json", {"action": "locus", "command": {
                "action": "export", "report": base.bound(root / "narrow/report.json"), "view": fixture.view(),
                "formats": ["json"], "output_dir": str(root / "narrow-json"), "max_output_bytes": 32 * 1024 * 1024, "cancel_path": None}})
            narrow_page = json.loads((root / "narrow-json/artifacts/page.json").read_bytes())
            self.assertTrue(any(r["reference_anchor"] is not None and r["reference_anchor"] < 200
                                for read in narrow_page["reads"] for e in read["genomic_basecalls"]
                                for v in e["variants"] for r in v["runs"]))
            self.assertFalse(state.exists())
