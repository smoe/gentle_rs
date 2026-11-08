"""Synthetic shared UTR corpus against the unchanged Python reference, not Rust execution.

Origin/recreation/use: test_files/rna_read_end_core/README.md. No private data,
external tool, frozen generated-file bytes or new Python business logic.
"""

from __future__ import annotations

import copy
import json
import unittest

from test_project_rna_read_end_metrics import MODULE, ROOT


CORPUS = ROOT / "test_files/rna_read_end_core/cases.json"


def inputs(corpus, case):
    document = copy.deepcopy(corpus["binding"])
    document["transcripts"] = [copy.deepcopy(corpus["models"][name]) for name in case["models"]]
    analysis = {
        "reference": {"assembly": document["assembly"], "fasta_sha256": document["reference_fasta_sha256"]},
        "annotation_bindings": [{"source_id": document["source_id"], "annotation_sha256": document["annotation_sha256"]}],
    }
    observations = {}
    for row in case["observations"]:
        observations[row["id"]] = {
            "observation_id": row["id"], "contig": "synthetic_contig",
            "query_axis_length": row["axis"], "query_to_reference_strand": row["strand"],
            "blocks": [{"query": b[:2], "target": b[2:]} for b in row["blocks"]],
            "comparisons": [{"source_id": document["source_id"], "availability": "available",
                             "comparison": {"resolution": "resolved", "models": [
                                 {"transcript_id": identifier, "relation": row.get("relation", "exact_partial_chain")}
                                 for identifier in row["candidates"]]}}],
        }
    return document, analysis, observations


class SharedUtrCoreContractTests(unittest.TestCase):
    def test_shared_literals_match_existing_python_for_lf_and_crlf(self):
        raw = CORPUS.read_bytes()
        for newline in (b"\n", b"\r\n"):
            corpus = json.loads(raw.replace(b"\r\n", b"\n").replace(b"\n", newline))
            self.assertEqual(corpus["schema"], "gentle.synthetic_utr_core_cases.v1")
            self.assertEqual(len(corpus["cases"]), 31)
            names = [row["name"] for row in corpus["cases"]]
            self.assertEqual(len(names), len(set(names)))
            for case in corpus["cases"]:
                with self.subTest(newline=newline, case=case["name"]):
                    document, analysis, observations = inputs(corpus, case)
                    before = copy.deepcopy((document, analysis, observations))
                    source_id, models = MODULE.read_models(document, analysis)
                    values, reason, ids = MODULE.project_unit(
                        {"hypotheses": case["hypotheses"]}, observations, source_id, models)
                    if values is None:
                        values = ((None, "missing", reason),) * 2
                    genes = sorted({models[identifier].gene_id for identifier in ids}) if set(ids) <= models.keys() else []
                    actual = {"five": list(values[0]), "three": list(values[1]), "ids": ids, "genes": genes}
                    self.assertEqual(actual, case["expected"])
                    self.assertEqual((document, analysis, observations), before)

    def test_core_is_registered_without_new_operation_or_parser_route(self):
        parent = (ROOT / "src/observed_transcript_structures.rs").read_text()
        self.assertIn("pub mod read_ends;", parent)
        core = (ROOT / "src/observed_transcript_structures/read_ends.rs").read_text()
        self.assertIn("geometry::project(", core)
        for forbidden in ("std::fs", "std::process", "GENES =", "barcode15", "barcode16"):
            self.assertNotIn(forbidden, core)

    def test_clip_descriptors_do_not_change_python_censoring_under_either_encoding(self):
        raw = CORPUS.read_bytes()
        for newline in (b"\n", b"\r\n"):
            corpus = json.loads(raw.replace(b"\r\n", b"\n").replace(b"\n", newline))
            self.assertEqual(len(corpus["clip_cases"]), 6)
            for entry in corpus["clip_cases"]:
                with self.subTest(case=entry, newline=newline):
                    case = copy.deepcopy(next(c for c in corpus["cases"] if c["name"] == entry["case"]))
                    if "axis" in entry:
                        case["axis"] = entry["axis"]
                        for row in case["observations"]:
                            row["axis"] = entry["axis"]
                    document, analysis, observations = inputs(corpus, case)
                    observations["o1"]["clipping"] = None if entry["soft"] is None else {"soft": entry["soft"], "hard": entry["hard"]}
                    source, models = MODULE.read_models(document, analysis)
                    values, _reason, _ids = MODULE.project_unit({"hypotheses": case["hypotheses"]}, observations, source, models)
                    self.assertEqual(list(values[1]), entry["utr3"])


if __name__ == "__main__":
    unittest.main()
