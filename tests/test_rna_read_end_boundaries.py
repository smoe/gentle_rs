"""Hand-crafted synthetic UTR boundary and actual Python-CLI regressions.

No external annotation, sequence or private data is used. Fixtures are created
inline in temporary directories; recreate with unittest discovery for
test_*rna_read_end*.py. The tests exercise the standalone projector, GTF model
preparer and unchanged comparator, not Rust, aligner or biological acceptance.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from test_project_rna_read_end_metrics import (
    ANNOTATION_SHA, ROOT, MODULE, data_document, model, model_document,
    observation,
)

import prepare_rna_read_end_transcript_models as preparation


def units(axis, hypotheses=None):
    return [{
        "unit": {"source_sha256": "sha256:" + "d" * 64, "unit_id": "read1"},
        "unclipped_axis_length": axis,
        "hypotheses": hypotheses or [["obs1"]],
    }]


class BoundaryTests(unittest.TestCase):
    def test_normal_ci_includes_complete_python_scope_on_linux_and_windows(self):
        workflow = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
        release_policy = workflow.split("\n  release-policy:", 1)[1].split("\n  select-platform:", 1)[0]
        windows = workflow.split("\n  windows:", 1)[1].split("\n  ci-summary:", 1)[0]
        scope = "-m unittest discover -s tests -p 'test_*rna_read_end*.py' -v"
        self.assertIn("python3 " + scope, release_policy)
        self.assertIn("python " + scope, windows)

    def project(self, models=None, observations=None, axis=25, hypotheses=None):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        data = data_document(observations=observations, units=units(axis, hypotheses))
        (root / "data.json").write_text(json.dumps(data), encoding="utf-8")
        (root / "models.json").write_text(json.dumps(models or model_document()), encoding="utf-8")
        try:
            MODULE.run(root / "data.json", root / "models.json", "synthetic_sample", root / "out")
        except ValueError:
            self.assertFalse((root / "out").exists())
            raise
        def rows(name):
            with (root / "out" / name).open(encoding="utf-8", newline="") as handle:
                return list(csv.DictReader(handle, delimiter="\t"))
        return rows("metrics.tsv"), rows("projections.tsv"), json.loads((root / "out" / "report.json").read_text())

    def test_missing_compatible_model_preserves_uncertainty_and_ids(self):
        for candidates in (["tx1", "without_cds"], ["without_cds"]):
            with self.subTest(candidates=candidates):
                metrics, details, _ = self.project(observations=[observation(models=candidates)])
                self.assertEqual([row["observation_status"] for row in metrics[1:]], ["missing", "missing"])
                self.assertEqual(details[2]["reason"], "candidate_model_unavailable")
                self.assertEqual(details[2]["transcript_ids"], ",".join(sorted(candidates)))
                self.assertEqual(details[2]["gene_ids"], "")

    def test_complete_agreeing_models_keep_observed_values(self):
        metrics, details, _ = self.project(
            models=model_document([model("tx1"), model("tx2")]),
            observations=[observation(models=["tx1", "tx2"])],
        )
        self.assertEqual([row["value_bp"] for row in metrics[1:]], ["5", "5"])
        self.assertEqual(details[2]["transcript_ids"], "tx1,tx2")

    def test_candidate_intersection_is_computed_before_model_availability(self):
        first = observation(models=["tx1", "without_cds"], blocks=[{"query": [0, 10], "target": [100, 110]}])
        second = observation("obs2", models=["tx1"], blocks=[{"query": [10, 25], "target": [120, 135]}])
        metrics, details, _ = self.project(observations=[first, second], hypotheses=[["obs1", "obs2"]])
        self.assertEqual([row["value_bp"] for row in metrics[1:]], ["5", "5"])
        self.assertEqual(details[2]["transcript_ids"], "tx1")

    def test_unknown_model_in_another_optimal_path_is_not_discarded(self):
        metrics, details, _ = self.project(
            observations=[observation(), observation("obs2", models=["without_cds"])],
            hypotheses=[["obs1"], ["obs2"]],
        )
        self.assertEqual(metrics[2]["observation_status"], "missing")
        self.assertEqual(details[2]["transcript_ids"], "tx1,without_cds")

    def test_utr_only_boundary_touch_is_not_a_cds_anchor_on_either_strand(self):
        cases = (("+", [100, 105], 1), ("+", [130, 135], 2),
                 ("-", [135, 140], 1), ("-", [100, 105], 2))
        for strand, target, metric in cases:
            with self.subTest(strand=strand, target=target):
                obs = observation(strand=strand, axis=5, blocks=[{"query": [0, 5], "target": target}])
                metrics, details, _ = self.project(models=model_document([model(strand=strand)]), observations=[obs], axis=5)
                self.assertEqual(metrics[metric]["observation_status"], "missing")
                self.assertEqual(details[metric]["reason"], "cds_boundary_unanchored")

    def test_coding_side_anchored_zero_remains_observed_on_either_strand(self):
        cases = (("+", [105, 110], 1), ("+", [125, 130], 2),
                 ("-", [130, 135], 1), ("-", [105, 110], 2))
        for strand, target, metric in cases:
            with self.subTest(strand=strand, target=target):
                obs = observation(strand=strand, axis=5, blocks=[{"query": [0, 5], "target": target}])
                metrics, _, _ = self.project(models=model_document([model(strand=strand)]), observations=[obs], axis=5)
                self.assertEqual((metrics[metric]["value_bp"], metrics[metric]["observation_status"]), ("0", "observed"))

    def test_disconnected_coding_block_cannot_anchor_a_utr_only_block(self):
        obs = observation(axis=10, blocks=[
            {"query": [0, 5], "target": [120, 125]},
            {"query": [5, 10], "target": [130, 135]},
        ])
        metrics, _, _ = self.project(observations=[obs], axis=10)
        self.assertEqual(metrics[2]["observation_status"], "missing")

    def test_unknown_partial_and_legacy_stop_boundaries_do_not_invent_utr3(self):
        documents = [model_document([model(stops=[])]), model_document([model(stops=[[130, 132]])])]
        legacy = model_document([model()])
        legacy["schema"] = MODULE.LEGACY_MODEL_SCHEMA
        del legacy["transcripts"][0]["genomic_stop_codon_0based_half_open"]
        documents.append(legacy)
        for document in documents:
            with self.subTest(document=document):
                metrics, details, _ = self.project(models=document)
                self.assertEqual(metrics[1]["value_bp"], "5")
                self.assertEqual((metrics[2]["value_bp"], metrics[2]["observation_status"]), ("", "missing"))
                self.assertEqual(details[2]["reason"], "stop_codon_unavailable")

    def test_missing_v2_stop_field_is_not_silently_read_as_legacy(self):
        document = model_document()
        del document["transcripts"][0]["genomic_stop_codon_0based_half_open"]
        with self.assertRaisesRegex(ValueError, "unexpected fields"):
            self.project(models=document)

    def test_separate_included_and_overlapping_stop_codons_have_same_boundary(self):
        for strand in ("+", "-"):
            for status in ("separate", "included", "overlapping"):
                with self.subTest(strand=strand, status=status):
                    if strand == "+":
                        end = {"separate": 130, "included": 133, "overlapping": 131}[status]
                        row = model(cds=[[105, 110], [120, end]], stops=[[130, 133]])
                        target = [125, 135]
                    else:
                        start = {"separate": 108, "included": 105, "overlapping": 107}[status]
                        row = model(strand="-", cds=[[start, 110], [120, 135]], stops=[[105, 108]])
                        target = [103, 110]
                    axis = target[1] - target[0]
                    obs = observation(strand=strand, axis=axis, blocks=[{"query": [0, axis], "target": target}])
                    metrics, _, report = self.project(models=model_document([row]), observations=[obs], axis=axis)
                    self.assertEqual(metrics[2]["value_bp"], "2")
                    self.assertEqual(report["transcript_models"]["stop_codon_evidence"][status], 1)
                    parsed = MODULE.prepare_model(row, 0)
                    self.assertEqual(parsed.cds, tuple(map(tuple, row["genomic_cds_0based_half_open"])))

    def test_split_stop_codons_project_in_transcript_order_on_both_strands(self):
        cases = (
            ("+", [[105, 109]], [[109, 110], [120, 122]], [
                {"query": [0, 3], "target": [107, 110]},
                {"query": [3, 7], "target": [120, 124]},
            ], 7),
            ("-", [[121, 135]], [[108, 110], [120, 121]], [
                {"query": [4, 8], "target": [106, 110]},
                {"query": [0, 4], "target": [120, 124]},
            ], 8),
        )
        for strand, cds, stops, blocks, axis in cases:
            with self.subTest(strand=strand):
                row = model(strand=strand, cds=cds, stops=stops)
                obs = observation(strand=strand, axis=axis, blocks=blocks)
                metrics, _, _ = self.project(models=model_document([row]), observations=[obs], axis=axis)
                self.assertEqual(metrics[2]["value_bp"], "2")

    def test_invalid_stop_geometry_is_refused_before_output_publication(self):
        for stops in ([[110, 113]], [[122, 125]], [[130, 134]], [[125, 126], [130, 132]]):
            with self.subTest(stops=stops), self.assertRaises(ValueError):
                self.project(models=model_document([model(stops=stops)]))

    def test_gtf_preparation_retains_raw_cds_and_split_stop_intervals(self):
        for strand, cds, stops in (("+", [[105, 109]], [[109, 110], [120, 122]]),
                                   ("-", [[121, 135]], [[108, 110], [120, 121]])):
            with self.subTest(strand=strand):
                lines = []
                for kind, intervals in (("exon", [[100, 110], [120, 140]]), ("CDS", cds), ("stop_codon", stops)):
                    for start, end in intervals:
                        lines.append(f'chr1\tsynthetic\t{kind}\t{start + 1}\t{end}\t.\t{strand}\t0\tgene_id "gene1"; transcript_id "tx1";')
                row = preparation.parse_gtf("\n".join(lines) + "\n")[0]
                self.assertEqual(row["genomic_cds_0based_half_open"], cds)
                self.assertEqual(row["genomic_stop_codon_0based_half_open"], stops)
                self.assertEqual(MODULE.prepare_model(row, 0).stop_codon_evidence, "separate")

    def test_gtf_refuses_invalid_complete_stop_before_model_publication(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf, reference, output = root / "input.gtf", root / "reference.fa", root / "model.json"
            attributes = 'gene_id "gene1"; transcript_id "tx1";'
            gtf.write_text("\n".join(
                f"chr1\tsynthetic\t{kind}\t{start}\t{end}\t.\t+\t0\t{attributes}"
                for kind, start, end in (("exon", 101, 140), ("CDS", 106, 130), ("stop_codon", 101, 103))
            ) + "\n", encoding="utf-8")
            reference.write_text(">chr1\n" + "ACGT" * 50 + "\n", encoding="utf-8")
            before = {path: path.read_bytes() for path in (gtf, reference)}
            with self.assertRaisesRegex(ValueError, "terminal CDS boundary"):
                preparation.run(gtf, reference, "synthetic-v1", "annotation1", ANNOTATION_SHA, output)
            self.assertFalse(output.exists())
            self.assertEqual({path: path.read_bytes() for path in before}, before)

    def test_actual_python_cli_chain_preserves_censoring_hashes_and_inputs(self):
        semantic_results = []
        for newline in ("\n", "\r\n"):
            with self.subTest(newline=repr(newline)), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary) / "synthetic UTR's \u00e4 paths"
                root.mkdir()
                attributes = 'gene_id "gene1"; transcript_id "tx1";'
                gtf = root / "synthetic models.gtf"
                gtf.write_bytes((newline.join(
                    f"chr1\tsynthetic\t{kind}\t{start}\t{end}\t.\t+\t0\t{attributes}"
                    for kind, start, end in (("exon", 101, 110), ("exon", 121, 140),
                                             ("CDS", 106, 110), ("CDS", 121, 130),
                                             ("stop_codon", 131, 133))
                ) + newline).encode())
                reference = root / "reference.fa"
                reference.write_bytes((">chr1" + newline + "ACGT" * 50 + newline).encode())
                model_path = root / "models.json"
                before = {path: path.read_bytes() for path in (gtf, reference)}
                result = subprocess.run([
                    sys.executable, str(ROOT / "scripts/prepare_rna_read_end_transcript_models.py"),
                    "--gtf", str(gtf), "--reference-fasta", str(reference),
                    "--assembly", "synthetic-v1", "--source-id", "annotation1",
                    "--annotation-sha256", ANNOTATION_SHA, "--output", str(model_path),
                ], capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                document = json.loads(model_path.read_text())
                self.assertEqual(document["model_source"]["sha256"], "sha256:" + hashlib.sha256(gtf.read_bytes()).hexdigest())
                self.assertEqual(document["reference_fasta_sha256"], "sha256:" + hashlib.sha256(reference.read_bytes()).hexdigest())
                before[model_path] = model_path.read_bytes()
                outputs = []
                for sample in ("sample1", "sample2"):
                    observed = [observation(identifier=f"obs{i}", axis=axis, qend=25) for i, axis in enumerate((25, 25, 30))]
                    unit_rows = [
                        {"unit": {"source_sha256": "sha256:" + "d" * 64, "unit_id": f"read{i}"},
                         "unclipped_axis_length": axis, "hypotheses": [[f"obs{i}"]]}
                        for i, axis in enumerate((25, 25, 30))
                    ]
                    data = root / f"{sample}.json"
                    payload = data_document(observed, unit_rows)
                    payload["analysis"]["reference"]["fasta_sha256"] = document["reference_fasta_sha256"]
                    data.write_bytes((json.dumps(payload) + newline).encode())
                    before[data] = data.read_bytes()
                    output = root / sample
                    result = subprocess.run([
                        sys.executable, str(ROOT / "scripts/project_rna_read_end_metrics.py"),
                        "--observed-sam-data", str(data), "--transcript-models", str(model_path),
                        "--sample-id", sample, "--output-dir", str(output),
                    ], capture_output=True, text=True, timeout=30)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    report = json.loads((output / "report.json").read_text())
                    self.assertEqual(report["source"]["sha256"], hashlib.sha256(data.read_bytes()).hexdigest())
                    self.assertEqual(report["transcript_models"]["sha256"], hashlib.sha256(model_path.read_bytes()).hexdigest())
                    with (output / "projections.tsv").open(newline="", encoding="utf-8") as handle:
                        clipped = [row for row in csv.DictReader(handle, delimiter="\t") if row["observation_status"] == "right_censored"]
                    self.assertEqual([(row["metric"], row["value_bp"]) for row in clipped], [("utr3_reach_bp", "2")])
                    outputs.append(output)
                combined = root / "combined.tsv"
                lines = []
                for index, output in enumerate(outputs):
                    rows = (output / "metrics.tsv").read_text().splitlines()
                    lines.extend(rows if index == 0 else rows[1:])
                combined.write_bytes((newline.join(lines) + newline).encode())
                result = subprocess.run([
                    sys.executable, str(ROOT / "scripts/compare_rna_read_end_metrics.py"),
                    str(combined), str(root / "comparison"),
                ], capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                report = json.loads((root / "comparison/report.json").read_text())
                for sample in ("sample1", "sample2"):
                    row = next(row for row in report["counts"] if row["sample_id"] == sample and row["metric"] == "utr3_reach_bp")
                    self.assertEqual((row["observed"], row["right_censored"]), (2, 1))
                self.assertEqual(report["source"]["sha256"], hashlib.sha256(combined.read_bytes()).hexdigest())
                self.assertEqual({path: path.read_bytes() for path in before}, before)
                semantic_results.append((report["counts"], report["pairs"]))
        self.assertEqual(semantic_results[0], semantic_results[1])


if __name__ == "__main__":
    unittest.main()
