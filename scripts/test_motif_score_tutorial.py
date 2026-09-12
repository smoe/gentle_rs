"""Deterministic checks for tutorial 08.13 and its shared GENtle route.

The four-column matrix and DNA words are synthetic. Set GENTLE_TUTORIAL_BIN_DIR
to exercise the exact built CLI against the existing checksum-bound plus/minus
TSS fixture. No network, private package or genome-scale input is used.
"""

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = ROOT / "docs/tutorial/08-13_motif_logo_to_promoter_trace.md"
SOURCE = ROOT / "docs/tutorial/sources/08-13_motif_logo_to_promoter_trace.json"
REPRO = ROOT / "docs/tutorial/reproducibility/motif_logo_to_promoter_trace"
FIXTURE = ROOT / "test_files/fixtures/tss_profiles"
BIN_DIR = os.environ.get("GENTLE_TUTORIAL_BIN_DIR")
COUNTS = [[8.0, 1.0, 1.0, 0.0], [1.0, 7.0, 1.0, 1.0],
          [0.0, 1.0, 8.0, 1.0], [1.0, 1.0, 1.0, 7.0]]


def score(word, family, policy):
    baseline = max(1.0, max(sum(column) for column in COUNTS))
    total = 0.0
    for base, column in zip(word, COUNTS):
        index = "ACGT".index(base)
        if policy == "scanner_a1":
            p = (column[index] + 1.0) / (sum(column) + 4.0)
        else:
            adjusted = list(column)
            if sum(adjusted) < baseline:
                padding = (baseline - sum(adjusted)) / 4.0
                adjusted = [value + padding for value in adjusted]
            epsilon = baseline * 1e-9
            p = (adjusted[index] + epsilon) / (sum(adjusted) + 4.0 * epsilon)
        if family == "llr":
            total += math.log2(p / 0.25)
        else:
            total += math.log2((p / (1.0 - p)) / (0.25 / 0.75))
    return total


class MotifScoreTutorialSourceTests(unittest.TestCase):
    def test_catalog_source_and_local_links(self):
        source = json.loads(SOURCE.read_text(encoding="utf-8"))
        self.assertEqual(source["schema"], "gentle.tutorial_source.v4")
        self.assertEqual((source["catalog"]["group"],
                          source["catalog"]["group_position"]), ("08", 13))
        self.assertEqual(source["catalog"]["status"], "manual/hybrid")
        catalog = json.loads((ROOT / "docs/tutorial/catalog.json").read_text(encoding="utf-8"))
        entry, = [row for row in catalog["entries"] if row["id"] == source["id"]]
        self.assertEqual(entry["notes"], source["catalog"]["notes"])
        self.assertEqual(entry["title"], source["title"])
        for term in ["PWM/PSSM", "JASPAR", "TFBS", "pseudocounts", "tail probability",
                     "binding affinity", "TP73", "no loaded project required"]:
            self.assertIn(term, entry["notes"])
        text = TUTORIAL.read_text(encoding="utf-8")
        self.assertEqual(text.splitlines()[0], f'# {source["title"]}')
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if "://" not in target:
                self.assertTrue((TUTORIAL.parent / target.split("#")[0]).exists(), target)

    def test_worked_table_recomputes_from_declared_synthetic_counts(self):
        with (REPRO / "worked_scores.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        expected = {
            ("jaspar-mapping educational replay", "log2_relative_risk", "ACGT"):
                score("ACGT", "llr", "scanner_a1"),
            ("jaspar-mapping educational replay", "log_odds", "ACGT"):
                score("ACGT", "lor", "scanner_a1"),
            ("GENtle policy educational replay", "llr_bits", "ACGT"):
                score("ACGT", "llr", "gentle"),
            ("GENtle policy educational replay", "true_log_odds_bits", "ACGT"):
                score("ACGT", "lor", "gentle"),
            ("jaspar-mapping educational replay", "log2_relative_risk", "ACGC"):
                score("ACGC", "llr", "scanner_a1"),
            ("GENtle policy educational replay", "llr_bits", "ACGC"):
                score("ACGC", "llr", "gentle"),
            ("jaspar-mapping educational replay", "log2_relative_risk", "TCGA"):
                score("TCGA", "llr", "scanner_a1"),
            ("GENtle policy educational replay", "llr_bits", "TCGA"):
                score("TCGA", "llr", "gentle"),
        }
        found = {(row["producer"], row["formula"], row["word"]):
                 float(row["score"]) for row in rows}
        for key, value in expected.items():
            self.assertAlmostEqual(found[key], value, places=11, msg=str(key))
        tp73 = [row for row in rows if row["matrix"] == "MA0861.2"]
        self.assertEqual(len(tp73), 2)
        self.assertAlmostEqual(float(tp73[0]["score"]), 19.543680326691, places=11)
        self.assertAlmostEqual(float(tp73[1]["score"]),
                               -math.log10(4.0 ** -16), places=11)
        for row in rows:
            expected_units = ("-log10(probability)" if row["formula"].endswith("tail_log10")
                              else "bits")
            self.assertEqual(row["units"], expected_units)

    def test_teaching_diagram_sites_and_information_match_the_pfm(self):
        svg = ET.parse(REPRO / "motif_to_score.svg").getroot()
        sites = [node.text for node in svg.find(".//*[@id='aligned-sites']")]
        self.assertEqual(len(sites), 10)
        for position, counts in enumerate(COUNTS):
            self.assertEqual([sum(word[position] == base for word in sites)
                              for base in "ACGT"], counts)
        stacks = svg.find(".//*[@id='information-stacks']")
        scale = float(stacks.attrib["data-pixels-per-bit"])
        heights = {(int(node.attrib["data-position"]), node.attrib["data-base"]):
                   float(node.attrib["height"]) for node in stacks}
        for position, counts in enumerate(COUNTS, 1):
            probabilities = [count / sum(counts) for count in counts]
            ic = 2 + sum(p * math.log2(p) for p in probabilities if p)
            for base, p in zip("ACGT", probabilities):
                self.assertAlmostEqual(heights.get((position, base), 0),
                                       scale * p * ic, delta=0.0001)

    def test_fixture_integrity_orientations_and_matrix_pin(self):
        manifest = json.loads((FIXTURE / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual({row["geometry"]["strand"] for row in manifest["records"]}, {"+", "-"})
        sums = {}
        for line in (FIXTURE / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split(maxsplit=1)
            sums[name] = digest
        for name, digest in sums.items():
            self.assertEqual(hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest(), digest)
        panel = json.loads((REPRO / "shared_across_tss_panel.json").read_text(encoding="utf-8"))
        self.assertEqual(panel["scale_mode"], "shared_across_tss")
        self.assertEqual(panel["factors"][0]["source_id"], "MA0004.1")
        fasta = (REPRO / "tiny.fa").read_bytes()
        self.assertEqual(fasta, b">synthetic_chr\nACGTACGCTCGAACGT\n")
        self.assertEqual((REPRO / "tiny.fa.fai").read_text(encoding="ascii"),
                         "synthetic_chr\t16\t15\t16\t17\n")
        pfm = (REPRO / "synthetic_4bp.pfm").read_text(encoding="ascii")
        self.assertIn(">SYNTH4.1 SYNTHETIC4", pfm)
        self.assertIn("A [ 8 1 0 1 ]", pfm)

    def test_documented_implementation_revisions_and_limits(self):
        text = TUTORIAL.read_text(encoding="utf-8")
        for term in ["8b7b6b472c0d2df3177b31901ceab69932ab90b1",
                     "0883ec719abee70bafbb2e8abfcae45e5bc9bcd9",
                     "uniform_iid_quantized_conservative_survival_v2",
                     "shared_across_tss", "float32", "0.95", "4^-16",
                     "prepared package", "binding probability"]:
            self.assertIn(term.lower(), text.lower())
        for stale in ["promoter_design_artifact_slice_offline",
                      "promoter_gene_set_ortholog_cohort_offline",
                      "gene_set_ortholog_promoter_cohorts_offline"]:
            self.assertIn(stale, text)

    def test_committed_shared_engine_output_is_receipt_bound(self):
        generated = REPRO / "generated"
        report = json.loads((generated / "report.json").read_text(encoding="utf-8"))
        receipt = json.loads((generated / "receipt.json").read_text(encoding="utf-8"))
        self.assertIn("git.8b7b6b472c0d2df3177b31901ceab69932ab90b1",
                      report["producer_revision"])
        matrix, = report["panel_resolution"]["matrices"]
        self.assertEqual(matrix["specification"]["source_id"], "MA0004.1")
        self.assertEqual(matrix["matrix_sha256"],
                         "64b5dcdb8dec059a43ad95819d9eca9d40abeb136cb0735421b4cba15badc4f9")
        self.assertEqual(report["panel_resolution"]["panel"]["scale_mode"],
                         "shared_across_tss")
        self.assertEqual(len(report["windows"]), 2)
        inventory = receipt["outputs"]
        self.assertTrue(inventory)
        for name, digest in inventory.items():
            path = generated / name
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)
        provenance = json.loads((REPRO / "provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(provenance["gentle"]["matrix_sha256"], matrix["matrix_sha256"])
        for field, name in [("panel_sha256", "shared_across_tss_panel.json")]:
            self.assertEqual(provenance["gentle"][field],
                             hashlib.sha256((REPRO / name).read_bytes()).hexdigest())
        self.assertEqual(provenance["jaspar_mapping"]["matrix_file_sha256"],
                         hashlib.sha256((REPRO / "synthetic_4bp.pfm").read_bytes()).hexdigest())


@unittest.skipUnless(BIN_DIR, "Set GENTLE_TUTORIAL_BIN_DIR for actual CLI replay")
class MotifScoreTutorialReplayTests(unittest.TestCase):
    def test_tss_profiles_execute_through_shared_cli_engine(self):
        suffix = ".exe" if os.name == "nt" else ""
        cli = Path(BIN_DIR).resolve() / f"gentle_cli{suffix}"
        with tempfile.TemporaryDirectory(prefix="gentle-motif-score-tutorial-") as tmp:
            output = Path(tmp).resolve() / "output"
            command = [str(cli), "features", "tss-tfbs-profiles",
                       "--manifest", str(FIXTURE / "manifest.json"),
                       "--fasta", str(FIXTURE / "plus.fa"),
                       "--fasta", str(FIXTURE / "minus.fa"),
                       "--panel", str(REPRO / "shared_across_tss_panel.json"),
                       "--selection", str(FIXTURE / "selection.json"),
                       "--expected-genome-id", "synthetic-genome-v1",
                       "--expected-assembly", "synthetic-assembly-v1",
                       "--expected-annotation-release", "synthetic-annotation-v1",
                       "--output-dir", str(output), "--formats", "svg"]
            result = subprocess.run(command, cwd=ROOT, capture_output=True,
                                    text=True, timeout=180)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads((output / "report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["panel_resolution"]["panel"]["scale_mode"],
                             "shared_across_tss")
            self.assertEqual(len(report["windows"]), 2)
            self.assertEqual({window["record"]["geometry"]["strand"]
                              for window in report["windows"]}, {"+", "-"})
            self.assertTrue(list(output.glob("*.svg")))


if __name__ == "__main__":
    unittest.main()
