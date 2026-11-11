"""Synthetic C1 source/checkout checks and opt-in actual CLI rehearsal.

All inputs below are invented and recreated in private temporary directories.
No external provider, real gene, private read, archive or Cargo invocation is used.
Set GENTLE_CODING_MODELS_CLI to a freshly built same-source binary for process tests.
"""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def digest(raw):
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def quote(value):
    # GENtle's shared-shell grammar, not cmd.exe, PowerShell or POSIX execution.
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def fixture(root, crlf=False):
    """Invented one-exon transcript; MKP plus a separately annotated TAA stop."""
    root.mkdir(parents=True, exist_ok=True)

    def write(name, text):
        if crlf:
            text = text.replace("\n", "\r\n")
        raw = text.encode("utf-8")
        path = root / name
        path.write_bytes(raw)
        return {"path": str(path), "bytes": len(raw), "sha256": digest(raw)}

    sequence = "CCATGAAACCCTAAGG"
    reference_hash = digest(sequence.encode("ascii"))
    attrs = ('gene_id "synthetic-gene"; transcript_id "synthetic-transcript"; '
             'transcript_version "2"; transcript_biotype "protein_coding"; '
             'tag "mRNA_end_NF";')
    gtf = ""
    for kind, start, end, phase in [
        ("transcript", 101, 116, "."), ("exon", 101, 116, "."),
        ("CDS", 103, 111, "0"), ("start_codon", 103, 105, "0"),
        ("stop_codon", 112, 114, "0"), ("five_prime_utr", 101, 102, "."),
        ("three_prime_utr", 115, 116, "."),
    ]:
        extra = ' protein_id "synthetic-protein"; protein_version "3";' if kind == "CDS" else ""
        gtf += f"chrSynthetic\tsynthetic\t{kind}\t{start}\t{end}\t.\t+\t{phase}\t{attrs}{extra}\n"
    annotation = {
        "schema": "gentle.primer_payload_annotation.v1", "reference_id": "synthetic-source",
        "assembly": "synthetic-assembly", "release": "synthetic-1",
        "reference_fasta_sha256": reference_hash, "coordinate_convention": "zero_based_half_open",
        "transcripts": [{"transcript_id": "synthetic-transcript", "gene_ids": ["synthetic-gene"],
                         "target_name": "chrSynthetic", "target_length": 1000,
                         "genomic_contig": "chrSynthetic", "genomic_strand": "+",
                         "genomic_exons": [[100, 116]]}],
    }
    return {
        "schema": "gentle.coding_transcript_models_request.v1",
        "gtf": write("coding.gtf", gtf),
        "normalized_annotation": write("annotation.json", json.dumps(annotation, indent=2) + "\n"),
        "reference": {"assembly": "synthetic-assembly", "contig": "chrSynthetic", "start_0based": 100,
                      "genomic_reference_sha256": reference_hash,
                      "locus_fasta": write("locus.fa", f">chrSynthetic\n{sequence}\n")},
        "selector": {"source_id": "synthetic-source", "gene_id": "synthetic-gene"},
        "provider": {"name": "synthetic", "release": "synthetic-1",
                     "cds_row_convention": "stop_codon_excluded_from_cds_rows"},
        "reference_kind": "nuclear", "translation_table": 1,
        "protein_fasta": write("protein.fa", ">synthetic-protein.3\nMKP*\n"), "budgets": {},
    }


class CodingModelSourceTests(unittest.TestCase):
    def test_invented_fixture_has_honest_raw_bindings_in_both_line_endings(self):
        with tempfile.TemporaryDirectory() as temp:
            for crlf in (False, True):
                request = fixture(Path(temp) / str(crlf), crlf)
                files = [request["gtf"], request["normalized_annotation"], request["protein_fasta"],
                         request["reference"]["locus_fasta"]]
                for binding in files:
                    raw = Path(binding["path"]).read_bytes()
                    self.assertEqual(binding["sha256"], digest(raw))
                    self.assertEqual(binding["bytes"], len(raw))
                    self.assertEqual(b"\r\n" in raw, crlf)

    def test_c1_routes_share_service_and_opt_in_gtf_grammar(self):
        service = (ROOT / "src/allele_hash_screen/coding_models.rs").read_text()
        model = (ROOT / "src/allele_hash_screen/coding_models/model.rs").read_text()
        parser = (ROOT / "src/observed_transcript_structures/annotation/gtf.rs").read_text()
        self.assertIn("gtf::parse_coding(line)?", service)
        self.assertIn("Some(length), false", parser)
        self.assertIn("None, true", parser)
        self.assertIn("AMINO_ACIDS.codon2aa", model)
        self.assertIn("read_ends::exon_map", model)
        self.assertNotIn("CompareCdnaAlleles", model)
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())
        rows = [row for row in glossary["commands"] if row["path"] == "rna-reads coding-models-derive"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["engine_operations"], ["DeriveCodingTranscriptModels"])


@unittest.skipUnless(os.environ.get("GENTLE_CODING_MODELS_CLI"), "same-source C1 CLI not supplied")
class CodingModelCliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="gentle-coding-models-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve() / "C1 O'Brien spaced \u03b1"
        self.binary = str(Path(os.environ["GENTLE_CODING_MODELS_CLI"]).resolve(strict=True))

    def invoke(self, request):
        path = self.root / "request's.json"
        path.write_text(json.dumps(request), encoding="utf-8")
        before = {p.name: p.read_bytes() for p in self.root.iterdir()}
        result = subprocess.run([self.binary, "--state", str(self.root / "must stay absent.gentle.json"),
                                 "shell", "rna-reads coding-models-derive " + quote("@" + str(path))],
                                cwd=self.root, capture_output=True, encoding="utf-8", timeout=60)
        self.assertEqual({p.name: p.read_bytes() for p in self.root.iterdir()}, before)
        return result

    def test_real_cli_lf_crlf_determinism_and_translation(self):
        models = []
        hashes = []
        for crlf in (False, True):
            request = fixture(self.root, crlf)
            first = self.invoke(request)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(self.invoke(request).stdout, first.stdout)
            report = json.loads(first.stdout)["result"]["coding_transcript_models"]
            model = report["models"][0]
            self.assertEqual(model["translation_status"], "complete_verified")
            self.assertEqual(model["coding"]["translated_protein_sha256"], digest(b"MKP"))
            self.assertEqual(model["coding"]["supplied_protein_comparison"], "identical")
            self.assertEqual(model["transcript_completeness_tags"], ["mRNA_end_NF"])
            models.append(model)
            hashes.append(report["sequence_sha256"])
        self.assertEqual(models[0], models[1])
        self.assertEqual(hashes[0], hashes[1])

    def test_real_cli_tampering_is_refused_without_state_or_output(self):
        request = fixture(self.root)
        request["gtf"]["sha256"] = digest(b"not the input")
        result = self.invoke(request)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("file_hash_or_change", result.stderr)
        self.assertNotIn('"coding_transcript_models":', result.stdout)


if __name__ == "__main__":
    unittest.main()
