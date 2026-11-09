#!/usr/bin/env python3
"""Create explicit synthetic inputs/command drafts; never run GENtle or process private data."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def bound(path: Path) -> dict:
    raw = path.read_bytes()
    return {"path": str(path.resolve()), "bytes": len(raw), "sha256": "sha256:" + hashlib.sha256(raw).hexdigest()}


def write(path: Path, value) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as target:
        json.dump(value, target, sort_keys=True, indent=2, allow_nan=False)
        target.write("\n")


def view() -> dict:
    return {"schema": "gentle.rna_read_end_pileup_view.v1", "sort": "unit_key", "descending": False,
            "filters": [], "offset": 0, "count": 20, "max_sublanes": 100, "width": 1400, "scale": "genomic",
            "channels": dict(utr5=True, utr3=True, read_length=True, clips=True, composition=True, tail=True)}


def prepare(root: Path, *, oriented: bool = False, strand: str = "+") -> None:
    if strand not in ("+", "-"):
        raise ValueError("explicit synthetic transcript strand required")
    root.mkdir(exist_ok=False)
    reference = {"assembly": "synthetic_assembly", "fasta_sha256": "sha256:" + "a" * 64}
    selector = {"source_id": "synthetic_annotation", "gene_id": "SYN_GENE"}
    sam = root / "synthetic.sam"
    with sam.open("x", encoding="ascii", newline="\n") as target:
        target.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:synthetic_contig\tLN:1000\n")
        originals = {"synthetic_exact": "A" * 12 + "C" * 70 + "T" * 18,
                     "synthetic_clipped": "A" * 12 + "C" * 85 + "T" * 18}
        for index, (name, cigar) in enumerate([("synthetic_exact", "50M50N50M"), ("synthetic_clipped", "5S50M50N50M10S")]):
            flag = 16 if oriented and ((index == 0) == (strand == "-")) else 0
            sequence = originals[name] if oriented else "*"
            if flag:
                sequence = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            target.write(f"{name}\t{flag}\tsynthetic_contig\t101\t60\t{cigar}\t*\t0\t0\t{sequence}\t*\n")
    if oriented:
        write(root / "original-sequences.json", originals)
    dictionary = root / "dictionary.json"
    write(dictionary, {"schema": "gentle.reference_contig_dictionary.v1", "reference_id": "synthetic_reference",
                       "assembly": reference["assembly"], "reference_fasta_sha256": reference["fasta_sha256"],
                       "contigs": [{"name": "synthetic_contig", "length": 1000}]})
    manifest = root / "units.json"
    search = {"search_complete": False, "views_complete": False, "limitations": ["Synthetic supplied records only; no exhaustive mapping search."],
              "producer": "GENtle synthetic rehearsal", "version": "1", "arguments": [], "selection_scope": "two declared units"}
    write(manifest, {"schema": "gentle.observed_sam_query_units.v1", "sam_sha256": bound(sam)["sha256"], "units": [
        {"unit_id": name, "unclipped_axis_length": size, "rows": [{"ordinal": i, "qname": name, "read_group": None}],
         "masks": None, "ends": [[], []], "search": search}
        for i, (name, size) in enumerate([("synthetic_exact", 100), ("synthetic_clipped", 115)], 1)]})
    annotation = root / "annotation.json"
    exons = [[100, 150], [200, 250]]
    ids = ["SYN_TX_A", "SYN_TX_B"]
    write(annotation, {"schema": "gentle.primer_payload_annotation.v1", "reference_id": "synthetic_reference",
                       "assembly": reference["assembly"], "release": "synthetic_snapshot", "reference_fasta_sha256": reference["fasta_sha256"],
                       "coordinate_convention": "zero_based_half_open", "transcripts": [
                           {"transcript_id": transcript, "gene_ids": [selector["gene_id"]], "target_name": "synthetic_contig", "target_length": 1000,
                            "genomic_contig": "synthetic_contig", "genomic_strand": strand, "genomic_exons": exons} for transcript in ids]})
    write(root / "models.json", {"schema": "gentle.rna_read_end_transcript_models.v2", "source_id": selector["source_id"],
        "assembly": reference["assembly"], "reference_fasta_sha256": reference["fasta_sha256"], "annotation_sha256": bound(annotation)["sha256"],
        "coordinate_convention": "zero_based_half_open", "model_source": {"format": "gtf", "name": "hand-crafted synthetic geometry, not a GTF retrieval", "sha256": "sha256:" + "b" * 64},
        "transcripts": [{"transcript_id": transcript, "gene_id": selector["gene_id"], "contig": "synthetic_contig", "strand": strand,
                         "genomic_exons_0based_half_open": exons,
                         "genomic_cds_0based_half_open": [[120, 150], [200, 220]] if strand == "+" else [[130, 150], [200, 230]],
                         "genomic_stop_codon_0based_half_open": [[220, 223]] if strand == "+" else [[127, 130]]} for transcript in ids]})
    text_limits = dict(max_input_bytes=1_000_000, max_line_bytes=65536, max_lines=100, max_rows=100, max_dictionary_entries=100,
                       max_header_ids=10, max_fields_per_line=64, max_operations=10000, max_sa_entries=100, max_retained_bytes=8*1024*1024, max_work=100000)
    kernel_limits = dict(max_queries=100, max_records=100, max_segments_per_query=100, max_operations=10000, max_skips_per_chain=100,
                         max_chain_families=100, max_alternatives_per_record=100, max_annotation_rows=100, max_model_comparisons_per_family=100, max_total_work=1_000_000)
    write(root / "sam-request.json", {"schema": "gentle.observed_sam_structures_request.v1", "sam": bound(sam), "dictionary": bound(dictionary), "manifest": bound(manifest),
        "conversion": {"kind": "direct_sam", "producer": "hand-crafted synthetic fixture", "version": "1", "arguments": [], "selection_scope": "two synthetic units"},
        "reference_id": "synthetic_reference", "reference": reference, "reference_gene": selector, "comparison_genes": [],
        "annotation_sources": [{"kind": "normalized", "selector": selector, "provider": "synthetic", "annotation": bound(annotation)}],
        "version_policy": "strict16", "unpadded_reference": True, "required_flank_bases": 2, "text_budgets": text_limits, "kernel_budgets": kernel_limits,
        "budgets": dict(max_input_bytes=1_000_000, max_annotation_bytes=65536, max_output_bytes=8*1024*1024, max_output_row_bytes=65536,
                        max_units=100, max_retained_bytes=16*1024*1024, max_projection_work=1_000_000),
        "output_dir": str(root / "saved-sam"), "cancel_path": None})


def bind_analysis(root: Path) -> None:
    report_path = root / "saved-sam/report.json"
    report = json.loads(report_path.read_bytes())
    if report["status"] != "completed":
        raise ValueError("completed synthetic SAM report required")
    binding = {"schema": "gentle.observed_sam_structures_result_binding.v1", "request_sha256": report["request_sha256"],
               "report_content_sha256": report["content_sha256"], "report_path": str(report_path), "outcome": report["status"]}
    write(root / "build-command.json", {"action": "build", "request": {"schema": "gentle.rna_read_end_build.v1",
        "source_id": "synthetic_annotation", "gene_id": "SYN_GENE", "transcript_id": "SYN_TX_A", "samples": [
            {"sample_id": "synthetic_sample", "analysis": binding, "model": bound(root / "models.json"), "normalized_annotation": bound(root / "annotation.json"),
             "composition": {"orientation": {"status": "undeclared"}, "trimming": {"status": "undeclared"}, "adapter_sequences": [], "primer_sequences": [],
                             "window_bp": 20, "minimum_run_bp": 10, "maximum_mismatches": 1}, "tail": None}],
        "exact_quantile_probabilities": [0, 0.5, 1], "limits": dict(max_models=100, max_intervals=1000, max_units=100, max_observations=1000, max_work=1_000_000),
        "output_dir": str(root / "read-ends"), "max_output_bytes": 8*1024*1024, "cancel_path": str(root / "cancel.marker")}})


def bind_pileup(root: Path) -> None:
    report = bound(root / "read-ends/report.json")
    write(root / "inspect-command.json", {"action": "inspect", "report": report, "view": view()})
    write(root / "export-command.json", {"action": "export", "report": report, "view": view(), "formats": ["json", "tsv", "svg", "png", "pdf"],
          "output_dir": str(root / "exports"), "max_output_bytes": 16*1024*1024, "cancel_path": None})


def bind_originals(root: Path) -> None:
    """Bind only the independently generated synthetic originals, never reconstruct from SAM SEQ."""
    data = json.loads((root / "saved-sam/artifacts/analysis.json").read_bytes())["analysis"]
    originals = json.loads((root / "original-sequences.json").read_bytes())
    rows = []
    for unit in data["units"]:
        name = unit["unit"]["unit_id"]
        sequence = originals[name]
        if len(sequence) != unit["unclipped_axis_length"] or any(row["qname"] != name for row in unit["rows"]):
            raise ValueError("synthetic original query binding mismatch")
        rows.append({"unit": unit["unit"], "read_id": name, "sequence": sequence})
    bases = root / "original-bases.json"
    write(bases, {"schema": "gentle.rna_read_end_original_bases.v1",
                  "analysis_context_sha256": data["context_sha256"], "orientation": "original_sequencing_5_to_3",
                  "provenance": "Hand-crafted A12/Cn/T18 originals; separate from strand-oriented synthetic SAM SEQ.", "reads": rows})
    command = json.loads((root / "build-command.json").read_bytes())
    command["request"]["samples"][0]["original_bases"] = bound(bases)
    command["request"]["terminal_runs"] = {"window_bp": 60, "minimum_run_bp": 10}
    write(root / "oriented-build-command.json", command)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "prepare-oriented-plus", "prepare-oriented-minus", "bind-analysis", "bind-originals", "bind-pileup"])
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    root = args.directory.resolve()
    if args.stage.startswith("prepare-oriented-"):
        prepare(root, oriented=True, strand="-" if args.stage.endswith("minus") else "+")
    else:
        {"prepare": prepare, "bind-analysis": bind_analysis, "bind-originals": bind_originals, "bind-pileup": bind_pileup}[args.stage](root)


if __name__ == "__main__":
    main()
