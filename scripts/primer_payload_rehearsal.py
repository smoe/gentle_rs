"""Invented Phase 2 inputs only; no aligner, network, private data or GENtle build.

Four explicit stages consume actual preceding GENtle reports. Refuse existing
outputs. See test_files/primer_payload/README.md for provenance and commands.
"""

import argparse
import hashlib
import json
from pathlib import Path


def digest(value):
    return "sha256:" + hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()


def bound(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256="sha256:" + hashlib.sha256(data).hexdigest())


def write(path, value):
    with Path(path).open("x", encoding="utf-8", newline="\n") as output:
        json.dump(value, output, indent=2, ensure_ascii=False)
        output.write("\n")


def budgets():
    return dict(max_input_bytes=20_000_000, max_output_bytes=20_000_000,
                max_records=100, max_queries=100, max_query_bases=10_000,
                max_line_bytes=100_000, max_alignments_per_query=100,
                max_cigar_operations=1000, max_annotation_rows=100,
                max_reference_context_bases=100_000, max_chain_comparisons=100_000)


def inputs(out, prefix="toy"):
    out = Path(out).resolve()
    out.mkdir()
    cores = ["AGTCCGATCA", "TGCATACGTGGA"]
    payloads = ["GATCCATGAGCTAGTCGACATCGTAGCTACGATGCTAGCATCG",
                "CGTACAGTCAGATCGATCGGATACCTAGCTAGCATCTAGGCA"]
    rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    primers = [dict(primer_id=f"{prefix}_primer_{i}", assay_label="Invented teaching oligo",
                    sequence=core, components=[dict(component_id="core", role="target_core", start=0, end=len(core))],
                    core_orientation="sense_5prime", expected_targets=[dict(locus_id=f"{prefix}_locus", transcript_feature_id=None)],
                    variant="synthetic", source="scripts/primer_payload_rehearsal.py", source_sha256=None)
               for i, core in enumerate(cores)]
    reads = [cores[0] + payloads[0], rc(cores[1] + payloads[1])]
    fastq = "".join(f"@repeated_name\n{read}\n+\n{'I' * len(read)}\n" for read in reads)
    path = out / "invented reads.fastq"
    with path.open("x", encoding="ascii", newline="\n") as output:
        output.write(fastq)
    write(out / "trace.json", dict(schema="gentle.rna_primer_trace_request.v1",
        primer_panel=dict(kind="inline", panel=dict(schema="gentle.rna_primer_panel.v1", content_sha256="",
            primers=primers, policy=dict(max_full_edits=0, max_core_edits_for_attribution=0))),
        files=[dict(sample_id=f"{prefix}_sample", barcode=f"{prefix}_barcode", priming="declared_synthetic",
            path=str(path), declared_complete=True, expected_sha256=bound(path)["sha256"])],
        export_payload_fasta=True, output_dir=str(out / "phase-one")))
    write(out / "panel-draft-binding.json", json.loads((out / "trace.json").read_text())["primer_panel"])
    write(out / "fixture.json", dict(prefix=prefix, cores=cores, payloads=payloads))


def preparation(out, report_path):
    out, report_path = Path(out).resolve(), Path(report_path).resolve()
    report = json.loads(report_path.read_text(encoding="utf-8"))
    trace = dict(schema="gentle.rna_primer_trace_result_binding.v1",
        submitted_request_sha256=report["submitted_request_sha256"],
        effective_request_sha256=digest(report["request"]), primer_panel_sha256=report["primer_panel"]["content_sha256"],
        reference_panel_sha256=report["reference_panel_sha256"], report_content_sha256=report["content_sha256"],
        report_path=str(report_path), outcome=report["status"])
    write(out / "prepare.json", dict(schema="gentle.primer_payload_prepare_request.v1", trace=trace,
        file_ids=[f["file_id"] for f in report["files_completed"]],
        cohort=dict(cohort_id="invented_cohort", label="All verified synthetic units", all_verified_units=True,
            selection_provenance="All verified units of this two-read invented fixture"),
        source_mode="retained_primary_fasta", anchor_selection="primary", internal_primer_policy="mask",
        min_payload_bases=12, budgets=budgets(), output_dir=str(out / "queries"), cancel_path=None))


def sealed_trace(out, panel_result):
    out = Path(out).resolve()
    panel = json.loads(Path(panel_result).read_text(encoding="utf-8"))["result"]["rna_primer_panel"]
    trace = json.loads((out / "trace.json").read_text(encoding="utf-8"))
    trace["primer_panel"] = dict(kind="inline", panel=panel)
    write(out / "sealed-trace.json", trace)


def alignments(out, manifest_path):
    out, manifest_path = Path(out).resolve(), Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    fixture = json.loads((out / "fixture.json").read_text())
    prefix, cores, payloads = fixture["prefix"], fixture["cores"], fixture["payloads"]
    targets = [f"{prefix}_expected", f"{prefix}_foreign"]
    fasta = out / "invented reference.fa"
    with fasta.open("x", encoding="ascii", newline="\n") as output:
        for target, core, payload in zip(targets, cores, payloads):
            output.write(f">{target}\n{core}{payload}\n")
    fasta_binding = bound(fasta)
    reference_id = f"{prefix}_reference"
    annotation_path = out / "annotation.json"
    write(annotation_path, dict(schema="gentle.primer_payload_annotation.v1", reference_id=reference_id,
        assembly="invented_assembly", release="invented_release", reference_fasta_sha256=fasta_binding["sha256"],
        coordinate_convention="zero_based_half_open", transcripts=[dict(transcript_id=target, gene_ids=[target + "_gene"],
            target_name=target, target_length=len(core + payload), genomic_contig=None, genomic_strand=None, genomic_exons=[])
            for target, core, payload in zip(targets, cores, payloads)]))
    paf = out / "invented query-grouped.paf"
    with paf.open("x", encoding="ascii", newline="\n") as output:
        for query in manifest["queries"]:
            i = [p["primer_id"] for p in manifest["primer_panel"]["primers"]].index(query["primer_ids"][0])
            n, start = query["query_length"], len(cores[i])
            output.write("\t".join(map(str, [query["query_id"], n, 0, n, "+", targets[i], start + n,
                start, start + n, n, n, 60, f"cg:Z:{n}M", f"cs:Z::{n}", "NM:i:0", "tp:A:P"])) + "\n")
    write(out / "summarize.json", dict(schema="gentle.primer_payload_alignment_request.v1", trace=manifest["request"]["trace"],
        queries=bound(manifest_path), query_manifest_sha256=manifest["content_sha256"],
        references=[dict(reference_id=reference_id, scope="transcriptome", assembly="invented_assembly", release="invented_release",
            fasta=fasta_binding, annotation=bound(annotation_path))],
        expected_targets=[dict(primer_id=p["primer_id"], reference_id=reference_id, phase1_locus_id=f"{prefix}_locus",
            phase1_transcript_feature_id=None, gene_ids=[targets[0] + "_gene"], transcript_ids=[targets[0]]) for p in manifest["primer_panel"]["primers"]],
        batches=[dict(batch_id="invented_batch", reference_id=reference_id, paf=bound(paf), query_manifest_sha256=manifest["content_sha256"],
            producer="literal fixture, not aligner", producer_version="1", arguments=[], producer_exit_code=0, declared_complete=True,
            competitor_search_limitations="Only literal teaching rows; no search or completeness proof",
            accounting={q["query_id"]: "mapped" for q in manifest["queries"]})],
        policy=dict(strict_alignment_detail=True, min_segment_bases=12, min_identity_per_mille=800, min_mapq=1,
            dominant_fraction_per_mille=800, min_cross_locus_segment_bases=12, breakpoint_tolerance_bases=2, terminal_seed_lengths=[6, 8, 10, 12]),
        budgets=budgets(), output_dir=str(out / "alignment-evidence"), cancel_path=None))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["inputs", "seal", "prepare", "alignments"])
    parser.add_argument("output")
    parser.add_argument("report", nargs="?")
    args = parser.parse_args()
    if args.stage == "inputs":
        inputs(args.output)
    elif not args.report:
        parser.error("this stage requires an actual preceding report path")
    elif args.stage == "seal":
        sealed_trace(args.output, args.report)
    elif args.stage == "prepare":
        preparation(args.output, args.report)
    else:
        alignments(args.output, args.report)


if __name__ == "__main__":
    main()
