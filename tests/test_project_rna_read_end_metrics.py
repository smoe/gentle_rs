"""Synthetic transcript-aware UTR reach projection regressions."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "project_rna_read_end_metrics.py"
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("project_rna_read_end_metrics", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

SHA = "sha256:" + "a" * 64
ANNOTATION_SHA = "sha256:" + "b" * 64


def model(transcript_id="tx1", strand="+", cds=None):
    return {
        "transcript_id": transcript_id,
        "gene_id": "gene1",
        "contig": "chr1",
        "strand": strand,
        "genomic_exons_0based_half_open": [[100, 110], [120, 140]],
        "genomic_cds_0based_half_open": cds or ([[105, 110], [120, 130]] if strand == "+" else [[105, 110], [120, 135]]),
    }


def model_document(transcripts=None):
    return {
        "schema": MODULE.MODEL_SCHEMA,
        "source_id": "annotation1",
        "assembly": "synthetic-v1",
        "reference_fasta_sha256": SHA,
        "annotation_sha256": ANNOTATION_SHA,
        "coordinate_convention": "zero_based_half_open",
        "model_source": {"format": "gtf", "name": "synthetic.gtf", "sha256": "sha256:" + "f" * 64},
        "transcripts": transcripts or [model()],
    }


def observation(identifier="obs1", strand="+", axis=25, qend=25, models=None, blocks=None):
    return {
        "observation_id": identifier,
        "on_optimal_path": True,
        "contig": "chr1",
        "query_axis_length": axis,
        "query_to_reference_strand": strand,
        "blocks": blocks or [
            {"query": [0, 10], "target": [100, 110]},
            {"query": [10, qend], "target": [120, 120 + qend - 10]},
        ],
        "comparisons": [
            {
                "source_id": "annotation1",
                "availability": "available",
                "comparison": {
                    "resolution": "resolved",
                    "models": [
                        {"transcript_id": transcript_id, "relation": "exact_partial_chain"}
                        for transcript_id in (models or ["tx1"])
                    ],
                },
            }
        ],
    }


def data_document(observations=None, units=None):
    return {
        "schema": MODULE.DATA_SCHEMA,
        "request_sha256": "sha256:" + "c" * 64,
        "intake": None,
        "annotations": [],
        "analysis": {
            "reference": {"assembly": "synthetic-v1", "fasta_sha256": SHA},
            "annotation_bindings": [
                {"source_id": "annotation1", "annotation_sha256": ANNOTATION_SHA}
            ],
            "observations": observations or [observation()],
            "units": units
            or [
                {
                    "unit": {"source_sha256": "sha256:" + "d" * 64, "unit_id": "read1"},
                    "unclipped_axis_length": 25,
                    "hypotheses": [["obs1"]],
                }
            ],
        },
    }


class ProjectionTests(unittest.TestCase):
    def execute(self, data=None, models=None):
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        data_path, model_path, output = root / "data.json", root / "models.json", root / "out"
        data_path.write_text(json.dumps(data or data_document()), encoding="utf-8")
        model_path.write_text(json.dumps(models or model_document()), encoding="utf-8")
        MODULE.run(data_path, model_path, "barcode01", output)
        with (output / "metrics.tsv").open(encoding="utf-8", newline="") as handle:
            metrics = list(csv.DictReader(handle, delimiter="\t"))
        with (output / "projections.tsv").open(encoding="utf-8", newline="") as handle:
            details = list(csv.DictReader(handle, delimiter="\t"))
        report = json.loads((output / "report.json").read_text(encoding="utf-8"))
        return directory, metrics, details, report

    def test_plus_strand_projects_read_length_and_both_utr_reaches(self):
        directory, metrics, details, report = self.execute()
        self.addCleanup(directory.cleanup)
        self.assertEqual(
            [(row["metric"], row["value_bp"], row["observation_status"]) for row in metrics],
            [
                ("read_length_bp", "25", "observed"),
                ("utr5_reach_bp", "5", "observed"),
                ("utr3_reach_bp", "5", "observed"),
            ],
        )
        self.assertEqual(details[1]["transcript_ids"], "tx1")
        self.assertEqual(report["counts"], {"missing": 0, "resolved": 3, "right_censored": 0})

    def test_clipped_transcript_facing_end_is_right_censored_not_zero(self):
        data = data_document(
            observations=[observation(axis=30, qend=25)],
            units=[
                {
                    "unit": {"source_sha256": "sha256:" + "d" * 64, "unit_id": "read1"},
                    "unclipped_axis_length": 30,
                    "hypotheses": [["obs1"]],
                }
            ],
        )
        directory, metrics, _details, report = self.execute(data=data)
        self.addCleanup(directory.cleanup)
        self.assertEqual(metrics[1]["observation_status"], "observed")
        self.assertEqual(metrics[2]["value_bp"], "5")
        self.assertEqual(metrics[2]["observation_status"], "right_censored")
        self.assertEqual(report["counts"]["right_censored"], 1)

    def test_minus_strand_uses_transcript_orientation(self):
        minus = model("tx_minus", "-")
        obs = observation(
            strand="-",
            axis=30,
            qend=30,
            models=["tx_minus"],
            blocks=[
                {"query": [20, 30], "target": [100, 110]},
                {"query": [0, 20], "target": [120, 140]},
            ],
        )
        directory, metrics, details, _report = self.execute(
            data=data_document(observations=[obs], units=[{
                "unit": {"source_sha256": "sha256:" + "d" * 64, "unit_id": "read1"},
                "unclipped_axis_length": 30,
                "hypotheses": [["obs1"]],
            }]),
            models=model_document([minus]),
        )
        self.addCleanup(directory.cleanup)
        self.assertEqual([(row["value_bp"], row["observation_status"]) for row in metrics[1:]], [("5", "observed"), ("5", "observed")])
        self.assertEqual(details[2]["transcript_ids"], "tx_minus")

    def test_transcript_geometry_disagreement_stays_missing(self):
        models = model_document([
            model("tx1"),
            model("tx2", cds=[[105, 110], [120, 128]]),
        ])
        obs = observation(models=["tx1", "tx2"])
        directory, metrics, details, _report = self.execute(data=data_document(observations=[obs]), models=models)
        self.addCleanup(directory.cleanup)
        self.assertEqual(metrics[1]["observation_status"], "observed")
        self.assertEqual(metrics[2]["observation_status"], "missing")
        self.assertEqual(details[2]["reason"], "ambiguous_transcript_end_geometry")
        self.assertEqual(details[2]["transcript_ids"], "tx1,tx2")

    def test_alignment_without_cds_boundary_anchor_is_missing(self):
        obs = observation(axis=4, qend=4, blocks=[{"query": [0, 4], "target": [131, 135]}])
        directory, metrics, details, _report = self.execute(data=data_document(observations=[obs]))
        self.addCleanup(directory.cleanup)
        self.assertEqual([row["observation_status"] for row in metrics[1:]], ["missing", "missing"])
        self.assertEqual([row["reason"] for row in details[1:]], ["resolved", "resolved"])

    def test_projection_rows_retain_gene_identity(self):
        directory, _metrics, details, _report = self.execute()
        self.addCleanup(directory.cleanup)
        self.assertEqual(details[0]["gene_ids"], "gene1")
        self.assertEqual(details[1]["gene_ids"], "gene1")
        self.assertEqual(details[2]["gene_ids"], "gene1")

    def test_model_binding_mismatch_fails_closed(self):
        models = model_document()
        models["annotation_sha256"] = "sha256:" + "e" * 64
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data_path, model_path = root / "data.json", root / "models.json"
            data_path.write_text(json.dumps(data_document()), encoding="utf-8")
            model_path.write_text(json.dumps(models), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "does not match exactly one saved annotation binding"):
                MODULE.run(data_path, model_path, "barcode01", root / "out")


if __name__ == "__main__":
    unittest.main()
