"""C1-C4 synthetic contract checks and opt-in actual CLI execution.

Every sequence, GTF and manifest is invented and recreated in private tempdirs.
The process companion reuses the existing synthetic saved-SAM producer helper,
not its old test campaign. No private reads, network or Cargo calls.
"""

import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
BIN = os.environ.get("GENTLE_CDNA_INTERPRETATION_CLI")
SEQUENCE = "CCATGTCCAAACCCGGGAAATAAGGGGGGG"


def bound(path):
    raw = path.read_bytes()
    return {"path": str(path), "sha256": "sha256:" + hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}


def write(path, value, ending=b"\n"):
    path.write_bytes((json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode().replace(b"\n", ending))
    return bound(path)


def quote(value):
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def coding_gtf(gene, transcript):
    attrs = f'gene_id "{gene}"; transcript_id "{transcript}"; transcript_version "2"; transcript_biotype "protein_coding"; tag "mRNA_end_NF";'
    rows = []
    for feature, a, z, frame in [
        ("transcript", 100, 140, "."), ("exon", 100, 110, "."), ("exon", 120, 140, "."),
        ("CDS", 102, 110, "0"), ("CDS", 120, 130, "1"), ("start_codon", 102, 105, "0"),
        ("stop_codon", 130, 133, "0"), ("five_prime_utr", 100, 102, "."), ("three_prime_utr", 133, 140, "."),
    ]:
        protein = ' protein_id "synthetic_protein"; protein_version "3";' if feature == "CDS" else ""
        rows.append(f"ctg_a\tsynthetic\t{feature}\t{a+1}\t{z}\t.\t+\t{frame}\t{attrs}{protein}\n")
    return "".join(rows)


def write_coding_inputs(root, ending):
    root.mkdir(parents=True, exist_ok=True)
    gtf = root / "coding O'Brien.gtf"
    gtf.write_bytes(coding_gtf("invented_gene", "invented_transcript").encode().replace(b"\n", ending))
    protein = root / "protein.fa"
    protein.write_bytes(b">synthetic_protein.3\nMSKPGK\n".replace(b"\n", ending))
    return bound(gtf), bound(protein)


def python_client(state, binary, cwd):
    with patch.object(sys, "path", [str(ROOT / "integrations/python"), *sys.path]):
        from gentle_py.client import GentleClient
    # An explicit argv prefix prevents the wrapper's optional Cargo fallback.
    return GentleClient(state_path=str(state), cli_cmd=[str(binary)], cwd=str(cwd))


class InterpretationSourceTests(unittest.TestCase):
    def test_python_delegation_keeps_json_and_text_separate_without_cargo_fallback(self):
        client = python_client("absent project.json", "synthetic CLI's \u03b1", ROOT)
        self.assertEqual(client.cli_prefix, ("synthetic CLI's \u03b1",))
        line = "rna-reads allele-hash-screen --interpret-cdna " + quote("@request's \u03b1.json")
        report = {"result": {"cdna_allele_interpretation": {"schema": "gentle.cdna_allele_interpretation.v1"}}}
        with patch("gentle_py.client.subprocess.run", return_value=subprocess.CompletedProcess([], 0, json.dumps(report), "")) as execute:
            self.assertEqual(client.shell(line, expect_json=True).json, report)
            self.assertEqual(execute.call_args.args[0], ["synthetic CLI's \u03b1", "--state", "absent project.json", "shell", line])
        for fmt, text in (("tsv", "record\tgene\nsite\tsynthetic\n"), ("summary", "# Conditional cDNA interpretation\n")):
            with patch("gentle_py.client.subprocess.run", return_value=subprocess.CompletedProcess([], 0, text, "")):
                result = client.shell(line + " --format " + fmt, expect_json=False)
                self.assertEqual(result.stdout, text)
                self.assertIsNone(result.json)

    def test_synthetic_coding_inputs_bind_raw_lf_crlf_and_hand_computed_codons(self):
        self.assertEqual(len(SEQUENCE), 30)
        self.assertEqual([SEQUENCE[i:i+3] for i in range(2, 23, 3)], ["ATG", "TCC", "AAA", "CCC", "GGG", "AAA", "TAA"])
        hashes = []
        with tempfile.TemporaryDirectory() as temp:
            for ending in (b"\n", b"\r\n"):
                files = write_coding_inputs(Path(temp) / str(len(ending)), ending)
                for file in files:
                    self.assertEqual(file, bound(Path(file["path"])))
                hashes.append(files[0]["sha256"])
        self.assertNotEqual(hashes[0], hashes[1])

    def test_discovery_and_adapters_use_the_same_operation(self):
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())
        rows = [c for c in glossary["commands"] if c["path"] == "rna-reads allele-hash-screen --interpret-cdna"]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["engine_operations"], ["InterpretCdnaAlleles"])
        self.assertIn("--format json|tsv|summary", rows[0]["usage"])
        service = (ROOT / "src/allele_hash_screen/interpretation.rs").read_text()
        self.assertIn("comparison::lower_cdna_files", service)
        self.assertIn("coding_models::derive_with_work_limit", service)
        self.assertIn("core::interpret", service)
        self.assertNotIn("QNAME", service)


@unittest.skipUnless(BIN, "requires freshly built same-source GENTLE_CDNA_INTERPRETATION_CLI")
class InterpretationCliTests(unittest.TestCase):
    def invoke(self, state, *args, success=True, raw=False):
        result = subprocess.run([BIN, "--state", str(state), *args], capture_output=True, encoding="utf-8", timeout=120)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
            return result.stdout if raw else json.loads(result.stdout)
        self.assertNotEqual(result.returncode, 0)
        return result.stderr

    def prepare(self, root, ending):
        from scripts.test_cdna_allele_comparison import CdnaAlleleComparisonCliReplay
        helper = CdnaAlleleComparisonCliReplay("test_actual_cli_saved_evidence_replay_lf_crlf_and_quoted_routes")
        helper.invoke = self.invoke
        state, comparison = helper.prepare(root, ending)
        locus = Path(comparison["reference"]["locus_fasta"]["path"])
        locus.write_bytes((">ctg_a\n" + SEQUENCE[:10] + "N" * 10 + SEQUENCE[10:] + "N" * 10 + "\n").encode().replace(b"\n", ending))
        comparison["reference"]["locus_fasta"] = bound(locus)
        for index, sample in enumerate(comparison["samples"]):
            path = Path(sample["basecalls"]["path"])
            calls = json.loads(path.read_bytes())
            for call in calls["reads"]:
                call["sequence"] = SEQUENCE[:6] + "AA" + SEQUENCE[8:] if index == 0 else SEQUENCE
            sample["basecalls"] = write(path, calls, ending)
        gtf, protein = write_coding_inputs(root, ending)
        annotation = json.loads(Path(comparison["normalized_annotation"]["path"]).read_bytes())
        model_request = {
            "schema": "gentle.coding_transcript_models_request.v1", "gtf": gtf,
            "normalized_annotation": comparison["normalized_annotation"], "reference": comparison["reference"],
            "selector": {"source_id": comparison["source_id"], "gene_id": comparison["gene_id"]},
            "provider": {"name": "synthetic", "release": annotation["release"], "cds_row_convention": "stop_codon_excluded_from_cds_rows"},
            "reference_kind": "nuclear", "translation_table": 1, "protein_fasta": protein,
        }
        path = root / "model request.json"
        write(path, model_request, ending)
        models = self.invoke(state, "rna-reads", "coding-models-derive", "@" + str(path))["result"]["coding_transcript_models"]
        self.assertEqual(models["models"][0]["translation_status"], "complete_verified")
        models = write(root / "models.json", models, ending)
        path = root / "comparison request.json"
        write(path, comparison, ending)
        report = self.invoke(state, "allele-hash-screen", "--compare-cdna", "@" + str(path))
        report = write(root / "comparison.json", report, ending)  # Exercise the real CLI envelope.
        request = {"schema": "gentle.cdna_allele_interpretation_request.v1", "genes": [{"comparison_report": report, "coding_models": models}],
                   "linkage": "same_unit", "policy": {}, "budgets": {}}
        path = root / "interpretation O'Brien \u03b1.json"
        write(path, request, ending)
        return state, request, path

    def test_actual_cli_c1_c4_end_to_end_lf_crlf_formats_and_unchanged_files(self):
        for ending in (b"\n", b"\r\n"):
            with tempfile.TemporaryDirectory(prefix="gentle-interpretation-") as temp:
                root = Path(temp).resolve() / "work O'Brien \u03b1"
                state, request, path = self.prepare(root, ending)
                manifest = {"schema": "gentle.known_variant_manifest.v1", "source_name": "synthetic", "source_release": "invented-1",
                            "retrieved_utc": "2026-11-12T00:00:00Z", "license": "synthetic-redistributable", "attribution": "Invented test data only",
                            "assembly": "synthetic", "entries": [{"contig": "ctg_a", "position_1based": 107, "ref": "CC", "alt": "AA", "identifier": "synthetic-composite"}]}
                request["known_variants"] = write(root / "known.json", manifest, ending)
                write(path, request, ending)
                before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
                direct = self.invoke(state, "allele-hash-screen", "--interpret-cdna", "@" + str(path))
                shell = self.invoke(state, "shell", "rna-reads allele-hash-screen --interpret-cdna " + quote("@" + str(path)))
                self.assertEqual(direct, shell)
                client = python_client(state, BIN, root)
                line = "rna-reads allele-hash-screen --interpret-cdna " + quote("@" + str(path))
                self.assertEqual(client.shell(line, expect_json=True).json, direct)
                report = direct["result"]["cdna_allele_interpretation"]
                combos = report["genes"][0]["codon_combinations"]
                self.assertTrue(any("stop_gained" in c["consequence"]["consequence_terms"] and c["known_variants"] for c in combos))
                for fmt in ("tsv", "summary"):
                    text = self.invoke(state, "allele-hash-screen", "--interpret-cdna", "@" + str(path), "--format", fmt, raw=True)
                    repeated = self.invoke(state, "shell", "rna-reads allele-hash-screen --interpret-cdna " + quote("@" + str(path)) + " --format " + fmt, raw=True)
                    self.assertEqual(text, repeated)
                    self.assertEqual(client.shell(line + " --format " + fmt, expect_json=False).stdout, text)
                    self.assertIn("Ser2Ter", text)
                self.assertEqual({p: p.read_bytes() for p in root.rglob("*") if p.is_file()}, before)
                self.assertFalse(state.exists())

    def test_actual_cli_tamper_duplicate_gene_and_assembly_refusals(self):
        with tempfile.TemporaryDirectory(prefix="gentle-interpretation-negative-") as temp:
            root = Path(temp).resolve() / "work O'Brien"
            state, request, path = self.prepare(root, b"\n")
            for fault in ("hash", "duplicate", "assembly"):
                altered = copy.deepcopy(request)
                if fault == "hash":
                    altered["genes"][0]["coding_models"]["sha256"] = "sha256:" + "0" * 64
                elif fault == "duplicate":
                    altered["genes"] *= 2
                else:
                    manifest = {"schema": "gentle.known_variant_manifest.v1", "source_name": "synthetic", "source_release": "invented-1",
                                "retrieved_utc": "2026-11-12T00:00:00Z", "license": "synthetic-redistributable", "attribution": "Invented test data only",
                                "assembly": "wrong-assembly", "entries": []}
                    altered["known_variants"] = write(root / "wrong assembly.json", manifest)
                write(path, altered)
                before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
                self.invoke(state, "allele-hash-screen", "--interpret-cdna", "@" + str(path), success=False)
                self.assertEqual({p: p.read_bytes() for p in root.rglob("*") if p.is_file()}, before)
                self.assertFalse(state.exists())


if __name__ == "__main__":
    unittest.main()
