"""Opt-in actual-process cDNA comparison replay; never builds Rust or uses sample data.

Origin/recreation: every SAM, JSON, sequence and quality is invented here in a
private temporary directory. Geometry mirrors the shared observed-SAM Rust test
fixture (ctg_a:100-110,120-140), not TP53/TP73 reference biology. Run with an
explicit source-matched GENTLE_CDNA_CLI executable to exercise real CLI routes.
"""

import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


BIN = os.environ.get("GENTLE_CDNA_CLI")
REFERENCE_HASH = "sha256:" + "1" * 64


def bound(path):
    raw = path.read_bytes()
    return dict(path=str(path), sha256="sha256:" + hashlib.sha256(raw).hexdigest(), bytes=len(raw))


def write_json(path, value, ending=b"\n"):
    raw = (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    path.write_bytes(raw.replace(b"\n", ending))
    return bound(path)


def quote(value):
    # GENtle shared-shell quoting, with no platform shell or command expansion.
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


@unittest.skipUnless(BIN, "requires source-matched GENTLE_CDNA_CLI; no Rust build")
class CdnaAlleleComparisonCliReplay(unittest.TestCase):
    def invoke(self, state, *args, success=True):
        run = subprocess.run(
            [BIN, "--state", str(state), *args], capture_output=True,
            text=True, encoding="utf-8", timeout=120,
        )
        if success:
            self.assertEqual(run.returncode, 0, run.stderr)
            return json.loads(run.stdout)
        self.assertNotEqual(run.returncode, 0, run.stdout)
        return run.stderr

    def prepare(self, root, ending):
        root.mkdir()
        state = root / "project must remain absent.json"
        dictionary = write_json(root / "synthetic dictionary's α.json", dict(
            schema="gentle.reference_contig_dictionary.v1",
            reference_id="synthetic_reference", assembly="synthetic",
            reference_fasta_sha256=REFERENCE_HASH,
            contigs=[dict(name="ctg_a", length=1000), dict(name="ctg_b", length=1000)],
        ), ending)
        annotation = write_json(root / "synthetic annotation's α.json", dict(
            schema="gentle.primer_payload_annotation.v1", reference_id="synthetic_reference",
            assembly="synthetic", release="invented_cli_fixture",
            reference_fasta_sha256=REFERENCE_HASH, coordinate_convention="zero_based_half_open",
            transcripts=[dict(
                transcript_id="invented_transcript", gene_ids=["invented_gene"],
                target_name="ctg_a", target_length=1000, genomic_contig="ctg_a",
                genomic_strand="+", genomic_exons=[[100, 110], [120, 140]],
            )],
        ), ending)
        locus = root / "synthetic locus's α.fa"
        locus.write_bytes(b">ctg_a invented locus\n".replace(b"\n", ending) + b"A" * 50 + ending)
        samples = []
        for index in range(2):
            qname = f"invented_read_{index}"
            sam = root / f"invented sample {index}'s α.sam"
            raw = ("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ctg_a\tLN:1000\n"
                   "@SQ\tSN:ctg_b\tLN:1000\n@RG\tID:library\n"
                   f"{qname}\t0\tctg_a\t101\t60\t10M10N20M\t*\t0\t0\t*\t*\tRG:Z:library\n")
            sam.write_bytes(raw.encode().replace(b"\n", ending))
            sam_binding = bound(sam)
            manifest = write_json(root / f"sample {index}'s units α.json", dict(
                schema="gentle.observed_sam_query_units.v1", sam_sha256=sam_binding["sha256"],
                units=[dict(
                    unit_id=f"invented_unit_{index}", unclipped_axis_length=30,
                    rows=[dict(ordinal=1, qname=qname, read_group="library")],
                    masks=dict(coordinate_convention="zero_based_half_open_unclipped_query",
                               unclipped_axis_length=30, intervals=[], source_sha256="sha256:" + "c" * 64,
                               derivation="Invented explicit empty mask; no biological record"),
                    ends=[[], []], search=dict(
                        search_complete=True, views_complete=True, limitations=[],
                        producer="Invented exhaustive fixture", version="1", arguments=[],
                        selection_scope="Only the invented row exists",
                    ),
                )],
            ), ending)
            selector = dict(source_id="invented_source", gene_id="invented_gene")
            output = root / f"sample {index}'s saved α"
            request = dict(
                schema="gentle.observed_sam_structures_request.v1", sam=sam_binding,
                manifest=manifest, dictionary=dictionary,
                conversion=dict(kind="direct_sam", producer="Hand-crafted fixture", version="1",
                                arguments=[], selection_scope="One invented mapping"),
                reference_id="synthetic_reference",
                reference=dict(assembly="synthetic", fasta_sha256=REFERENCE_HASH),
                reference_gene=selector, comparison_genes=[],
                annotation_sources=[dict(kind="normalized", selector=selector,
                                         provider="synthetic", annotation=annotation)],
                version_policy="strict16", unpadded_reference=True, required_flank_bases=2,
                text_budgets=dict(max_input_bytes=1000000, max_line_bytes=65536, max_lines=100,
                                  max_rows=100, max_dictionary_entries=100, max_header_ids=10,
                                  max_fields_per_line=64, max_operations=10000, max_sa_entries=100,
                                  max_retained_bytes=4194304, max_work=100000),
                kernel_budgets=dict(max_queries=100, max_records=100, max_segments_per_query=100,
                                    max_operations=10000, max_skips_per_chain=100, max_chain_families=100,
                                    max_alternatives_per_record=100, max_annotation_rows=100,
                                    max_model_comparisons_per_family=100, max_total_work=1000000),
                budgets=dict(max_input_bytes=1000000, max_annotation_bytes=65536,
                             max_output_bytes=2097152, max_output_row_bytes=65536, max_units=100,
                             max_retained_bytes=8388608, max_projection_work=1000000),
                output_dir=str(output), cancel_path=None,
            )
            request_path = root / f"sample {index}'s producer request α.json"
            write_json(request_path, request, ending)
            self.invoke(state, "rna-reads", "transcript-structures-from-sam", "@" + str(request_path))
            report_path = output / "report.json"
            report = json.loads(report_path.read_bytes())
            self.assertEqual(report["status"], "completed")
            analysis = json.loads((output / "artifacts/analysis.json").read_bytes())["analysis"]
            sequence = list("A" * 30)
            if index == 0:
                sequence[3], sequence[15] = "C", "G"
            calls = write_json(root / f"sample {index}'s original basecalls α.json", dict(
                schema="gentle.cdna_allele_basecalls.v1",
                analysis_context_sha256=analysis["context_sha256"],
                orientation="original_sequencing_5_to_3",
                provenance="Invented original bases/Phred40; no sequenced cell-line data",
                reads=[dict(unit=analysis["units"][0]["unit"], read_id=qname,
                            sequence="".join(sequence), qualities_phred33="I" * 30)],
            ), ending)
            samples.append(dict(
                sample_id=f"synthetic_{index}", label=f"Invented sample {index}; not a cell line",
                analysis=dict(schema="gentle.observed_sam_structures_result_binding.v1",
                              request_sha256=report["request_sha256"],
                              report_content_sha256=report["content_sha256"],
                              report_path=str(report_path), outcome="completed"),
                basecalls=calls,
            ))
        self.assertFalse(state.exists())
        return state, dict(
            schema="gentle.cdna_allele_comparison_request.v1", gene="Invented gene",
            source_id="invented_source", gene_id="invented_gene",
            reference=dict(assembly="synthetic", contig="ctg_a", start_0based=100,
                           genomic_reference_sha256=REFERENCE_HASH, locus_fasta=bound(locus)),
            normalized_annotation=annotation, samples=samples, supplied_sites=[],
            policy=dict(minimum_alternate_units=1, minimum_haplotype_units=1),
        )

    def test_actual_cli_saved_evidence_replay_lf_crlf_and_quoted_routes(self):
        for ending in (b"\n", b"\r\n"):
            with self.subTest(ending=ending), tempfile.TemporaryDirectory(prefix="cDNA CLI's ") as scratch:
                root = Path(scratch).resolve() / "invented work's α"
                state, request = self.prepare(root, ending)
                for display in ("TP53 (synthetic route label only)", "TP73 (synthetic route label only)"):
                    request["gene"] = display
                    path = root / "comparison request's α.json"
                    write_json(path, request, ending)
                    before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
                    direct = self.invoke(state, "allele-hash-screen", "--compare-cdna", "@" + str(path))
                    shell = self.invoke(state, "shell", "rna-reads allele-hash-screen --compare-cdna "
                                        + quote("@" + str(path)))
                    a, b = (x["result"]["cdna_allele_comparison"] for x in (direct, shell))
                    self.assertEqual(a, b)
                    self.assertEqual(a["content_sha256"], b["content_sha256"])
                    self.assertEqual([s["position_0based"] for s in a["comparison"]["sites"]], [103, 125])
                    self.assertEqual(a["comparison"]["sites"][0]["transcripts"][0]["cdna_position_1based"], 4)
                    counts = a["comparison"]["sites"][0]["samples"]
                    self.assertEqual([s["total_units"] for s in counts], [1, 1])
                    self.assertEqual(counts[0]["base_counts"], {"C": 1})
                    self.assertEqual(counts[1]["base_counts"], {"A": 1})
                    self.assertEqual({p: p.read_bytes() for p in root.rglob("*") if p.is_file()}, before)
                    self.assertFalse(state.exists())

    def test_actual_cli_missing_qualities_remain_unassessed_and_reference_mismatch_is_refused(self):
        with tempfile.TemporaryDirectory(prefix="cDNA negative's ") as scratch:
            root = Path(scratch).resolve() / "invented work's α"
            state, request = self.prepare(root, b"\n")
            calls_path = Path(request["samples"][0]["basecalls"]["path"])
            calls = json.loads(calls_path.read_bytes())
            calls["reads"][0]["qualities_phred33"] = None
            request["samples"][0]["basecalls"] = write_json(calls_path, calls)
            request["supplied_sites"] = [103]
            path = root / "comparison request's α.json"
            write_json(path, request)
            report = self.invoke(state, "allele-hash-screen", "--compare-cdna", "@" + str(path))
            site = report["result"]["cdna_allele_comparison"]["comparison"]["sites"][0]
            self.assertEqual(site["candidate_alternate_bases"], [])
            self.assertEqual(site["samples"][0]["unassessed_units"], 1)
            self.assertEqual(site["samples"][0]["callable_units"], 0)
            self.assertEqual(site["samples"][1]["base_counts"], {"A": 1})
            bad = copy.deepcopy(request)
            bad["reference"]["genomic_reference_sha256"] = "sha256:" + "f" * 64
            write_json(path, bad)
            before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
            self.invoke(state, "allele-hash-screen", "--compare-cdna", "@" + str(path), success=False)
            self.assertEqual({p: p.read_bytes() for p in root.rglob("*") if p.is_file()}, before)
            self.assertFalse(state.exists())


if __name__ == "__main__":
    unittest.main()
