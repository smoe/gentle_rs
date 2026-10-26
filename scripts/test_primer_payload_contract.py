"""Source/fixture checks, plus opt-in actual CLI acceptance; no Rust build or private input.

Fixture origin/recreation: test_files/primer_payload/README.md. The optional
process test needs GENTLE_PAYLOAD_CLI pointing to a source-matched gentle_cli.
"""

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from scripts import primer_payload_rehearsal as fixture

ROOT = Path(__file__).resolve().parents[1]
BIN = os.environ.get("GENTLE_PAYLOAD_CLI")


class PrimerPayloadSourceTests(unittest.TestCase):
    def test_mask_clipping_is_shared_and_matches_independent_basewise_oracle(self):
        base = ROOT / "src/engine/analysis/rna_panel/primer_payload"
        preparation = (base / "prepare.rs").read_text()
        binding = (base / "binding.rs").read_text()
        self.assertIn("payload_geometry(row.read_length, hit, &row.hits)?", preparation)
        self.assertIn("prepare::payload_geometry(row.read_length, h, &row.hits)?", binding)
        self.assertNotIn(".then_some(if", preparation + binding)
        # Independent oracle cross-check only, not execution of Rust.
        cases = 0
        for start in range(12):
            for end in range(start + 1, 13):
                for other_start in range(12):
                    for other_end in range(other_start + 1, 13):
                        for reverse in (False, True):
                            a, z = (0, start) if reverse else (end, 12)
                            raw_positions = list(range(a, z))
                            if reverse:
                                raw_positions.reverse()
                            observed = [i for i, raw in enumerate(raw_positions)
                                        if other_start <= raw < other_end]
                            s, e = max(a, other_start), min(z, other_end)
                            clipped = [] if s >= e else list(range(z - e, z - s) if reverse else range(s - a, e - a))
                            self.assertEqual(observed, clipped)
                            cases += 1
        self.assertEqual(cases, 12168)

    def test_invented_inputs_vary_oligos_orientation_and_opaque_names(self):
        for prefix in ("toy", "renamed_arbitrary_ids"):
            with tempfile.TemporaryDirectory() as scratch:
                root = Path(scratch) / "synthetic work's directory"
                fixture.inputs(root, prefix)
                data = json.loads((root / "fixture.json").read_text())
                reads = (root / "invented reads.fastq").read_text().splitlines()
                self.assertEqual(reads[0], reads[4])
                self.assertEqual(list(map(len, data["cores"])), [10, 12])
                self.assertNotIn(44, map(len, data["cores"]))
                rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                self.assertEqual(reads[1], data["cores"][0] + data["payloads"][0])
                self.assertEqual(rc(reads[5]), data["cores"][1] + data["payloads"][1])
                for read in (reads[1], reads[5]):
                    self.assertEqual(sum(read.count(core) + read.count(rc(core)) for core in data["cores"]), 1)
                request = json.loads((root / "trace.json").read_text())
                self.assertEqual(request["files"][0]["expected_sha256"], fixture.bound(root / "invented reads.fastq")["sha256"])
                self.assertEqual(request["primer_panel"]["panel"]["content_sha256"], "")
                with self.assertRaises(FileExistsError):
                    fixture.inputs(root, prefix)

    def test_new_schemas_are_strict_and_phase_one_public_fields_unchanged(self):
        text = (ROOT / "src/engine/protocol/primer_payload.rs").read_text()
        self.assertEqual(text.count("pub struct "), text.count("#[serde(deny_unknown_fields)]"))
        old = (ROOT / "src/engine/protocol/primer_trace.rs").read_text()
        self.assertNotIn("PrimerPayload", old)
        self.assertNotIn("primer_payload", old)
        self.assertNotIn("#[serde(default", text)
        for required in ("source_mode", "anchor_selection", "internal_primer_policy", "max_chain_comparisons", "expected_targets", "accounting", "competitor_search_limitations"):
            self.assertIn("pub " + required + ":", text)

    def test_shared_routes_and_no_aligner_or_current_panel_fixations(self):
        glossary = json.loads((ROOT / "docs/glossary.json").read_text())
        commands = {c["path"]: c for c in glossary["commands"]}
        routes = dict(prepare="PreparePrimerPayloadAlignmentQueries", summarize="SummarizePrimerPayloadAlignments", show="InspectPrimerPayload", cancel="CancelPrimerPayload")
        engine = (ROOT / "src/engine.rs").read_text()
        registry = (ROOT / "crates/gentle-protocol/src/lib.rs").read_text()
        managed = (ROOT / "src/command_execution/rna_artifacts.rs").read_text()
        for route, operation in routes.items():
            self.assertEqual(commands["rna-reads primer-payload-" + route]["engine_operations"], [operation])
            self.assertIn(operation, engine)
            self.assertIn('"' + operation + '"', registry)
            self.assertIn(operation, managed)
        for path in (ROOT / "src/engine/analysis/rna_panel/primer_payload").glob("*.rs"):
            if "tests" in path.name:
                continue
            code = path.read_text()
            for forbidden in ("Command::", "reqwest", "barcode13", "barcode14", "TP73", "PATZ1", "E2F1", "HDAC2", "GRCh38"):
                self.assertNotIn(forbidden, code, str(path))

    def test_paf_and_claim_guards_are_not_truncation_or_implicit_replay(self):
        base = ROOT / "src/engine/analysis/rna_panel/primer_payload"
        prepare = (base / "prepare.rs").read_text()
        parser = (base / "paf.rs").read_text()
        summary = (base / "summarize.rs").read_text()
        self.assertIn("ReplayVerifiedFastq =>", prepare)
        self.assertIn("retained_fasta_missing_use_explicit_replay", prepare)
        for guard in ("contradictory_detail", "cigar_span_or_counts", "cigar_operation_budget", "paf_no_hit"):
            self.assertIn(guard, parser)
        for guard in ("requires_query_grouping", "accounting_conflict", "accounting_coverage"):
            self.assertIn(guard, summary)
        self.assertNotIn(".truncate(", "\n".join(p.read_text() for p in base.glob("*.rs")))

    def test_docs_cover_limits_without_new_acceptance_claims(self):
        text = (ROOT / "docs/primer_payload_evidence.md").read_text()
        for required in ("Runtime", "pending", "not a signature", "128 MB", "255", "zero", "unmasked", "one", "cg", "cs", "cancel", "No molecule"):
            self.assertIn(required.lower(), text.lower())
        guide = (ROOT / "docs/tutorial/10-04_targeted_nanopore_fastq_panel.md").read_text()
        for stage in ("inputs", "seal", "prepare", "alignments"):
            self.assertIn("primer_payload_rehearsal.py " + stage, guide)
        for newline in ("\n", "\r\n"):
            parsed = json.loads((ROOT / "docs/glossary.json").read_text().replace("\n", newline))
            self.assertEqual(len([c for c in parsed["commands"] if "primer-payload-" in c["path"]]), 4)


@unittest.skipUnless(BIN, "requires explicitly supplied source-matched GENTLE_PAYLOAD_CLI; no build")
class PrimerPayloadCliAcceptance(unittest.TestCase):
    def test_actual_cli_rehearsal_lf_crlf_quoted_paths_and_no_source_overwrite(self):
        for newline in (b"\n", b"\r\n"):
            with tempfile.TemporaryDirectory(prefix="payload CLI's ") as scratch:
                root = Path(scratch) / "fixture work's directory"
                fixture.inputs(root)
                fastq = root / "invented reads.fastq"
                fastq.write_bytes(fastq.read_bytes().replace(b"\n", newline))
                request = json.loads((root / "trace.json").read_text())
                request["files"][0]["expected_sha256"] = fixture.bound(fastq)["sha256"]
                (root / "trace.json").write_text(json.dumps(request), encoding="utf-8")
                state = root / "state must remain absent.json"

                def run(*args):
                    proc = subprocess.run([BIN, "--state", str(state), *args], capture_output=True, text=True, timeout=120)
                    self.assertEqual(proc.returncode, 0, proc.stderr)
                    return json.loads(proc.stdout)

                sealed = run("rna-reads", "primer-panel-validate", "@" + str(root / "panel-draft-binding.json"), "--seal")
                fixture.write(root / "sealed.json", sealed)
                fixture.sealed_trace(root, root / "sealed.json")
                run("rna-reads", "primer-trace", "@" + str(root / "sealed-trace.json"))
                trace = root / "phase-one/report.json"
                original = trace.read_bytes()
                raw = fastq.read_bytes()
                fixture.preparation(root, trace)
                run("rna-reads", "primer-payload-prepare", "@" + str(root / "prepare.json"))
                fixture.alignments(root, root / "queries/report.json")
                run("rna-reads", "primer-payload-summarize", "@" + str(root / "summarize.json"))
                report = root / "alignment-evidence/report.json"
                result = run("rna-reads", "primer-payload-show", str(report), "--rows", "2", "--verify-artifacts")
                result = result["result"]["primer_payload_inspection"]
                self.assertTrue(result["artifact_content_verified"])
                self.assertEqual(sorted(r["classification"] for r in result["rows"]), ["expected_dominant", "off_target_dominant"])
                self.assertEqual(trace.read_bytes(), original)
                self.assertEqual(fastq.read_bytes(), raw)
                self.assertFalse(state.exists())
                # A completed child must have closed files before publication/reopen on Windows.
                renamed = report.with_name("renamed.json")
                report.rename(renamed)
                renamed.rename(report)


if __name__ == "__main__":
    unittest.main()
