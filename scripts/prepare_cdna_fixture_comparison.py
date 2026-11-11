#!/usr/bin/env python3
"""Bind Glen's explicit BAM/FASTQ fixture to existing GENtle annotation/SAM/comparison contracts.

This is evidence preparation and command orchestration, not an aligner, GTF/CIGAR
geometry interpreter, SNP caller or phaser. Run with an externally supplied
Python/pysam installation and source-matched gentle_cli. Every output directory
must be fresh; source fixture bytes are never changed. Original orientation is
provided by pysam and independently checked against the paired FASTQ records.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess


GENES = {
    "TP53": ("ENSG00000141510", "ENSG00000141510.21"),
    "TP73": ("ENSG00000078900", "ENSG00000078900.18"),
}
SOURCE_ID = "ensembl116_tp53_tp73_fixture"
REFERENCE_ID = "GRCh38_Ensembl116_primary_assembly"
MAX_RECORDS = 4096
MAX_BASES = 16_000_000


def sha_bytes(value: bytes) -> str:
    return "sha256:" + hashlib.sha256(value).hexdigest()


def bound(path: Path) -> dict:
    digest, size = hashlib.sha256(), 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
            size += len(chunk)
    return {"path": str(path.resolve()), "sha256": "sha256:" + digest.hexdigest(), "bytes": size}


def write_json(path: Path, value) -> dict:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
    return bound(path)


def table(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_inventory(data: Path) -> list[dict]:
    """Verify the complete supplied raw-byte inventory, including compressed bytes."""
    seen, result = set(), []
    for line in (data / "SHA256SUMS").read_text(encoding="ascii").splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (?:\./)?(.+)", line)
        if not match:
            raise ValueError("invalid SHA256SUMS line")
        digest, relative = match.groups()
        raw = PurePosixPath(relative)
        if raw.is_absolute() or any(p in ("", ".", "..") for p in relative.split("/")) or "\\" in relative:
            raise ValueError("unsafe fixture inventory path")
        path = data / relative
        if relative in seen or path.is_symlink() or not path.resolve().is_relative_to(data.resolve()):
            raise ValueError("duplicate/symlink/foreign fixture inventory path")
        record = bound(path)
        if record["sha256"] != "sha256:" + digest:
            raise ValueError(f"fixture hash mismatch: {relative}")
        seen.add(relative)
        result.append(record)
    actual = {str(p.relative_to(data)).replace(os.sep, "/") for p in data.rglob("*") if p.is_file()}
    if actual != seen | {"SHA256SUMS"}:
        raise ValueError("fixture files do not exactly match SHA256SUMS")
    return result


def region_sequences(path: Path) -> dict[str, str]:
    records, name, pieces = {}, None, []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(pieces)
            name, pieces = line[1:].split()[0], []
            if name in records:
                raise ValueError("duplicate FASTA region")
        elif line:
            if name is None or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
                raise ValueError("invalid reference FASTA sequence")
            pieces.append(line.upper())
    if name is not None:
        if name in records:
            raise ValueError("duplicate FASTA region")
        records[name] = "".join(pieces)
    return records


def read_fastq(path: Path) -> dict[str, tuple[int, str, str]]:
    records, bases = {}, 0
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        while True:
            head = handle.readline()
            if not head:
                break
            sequence, separator, qualities = [handle.readline().rstrip("\r\n") for _ in range(3)]
            head = head.rstrip("\r\n")
            if not head.startswith("@") or not head[1:] or not separator.startswith("+"):
                raise ValueError("invalid four-line FASTQ record")
            name = head[1:].split()[0]
            bases += len(sequence)
            if (name in records or not sequence or len(sequence) != len(qualities)
                    or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence)
                    or any(not 33 <= ord(q) <= 126 for q in qualities)
                    or len(records) >= MAX_RECORDS or bases > MAX_BASES):
                raise ValueError("FASTQ identity/sequence/quality/budget mismatch")
            records[name] = (len(records) + 1, sequence, qualities)
    if not records:
        raise ValueError("empty FASTQ")
    return records


def matched_original(record, fastq: tuple[int, str, str]) -> tuple[str, str]:
    """Use library orientation/length APIs; never reconstruct hard clips or CIGAR geometry."""
    _, sequence, qualities = fastq
    if record.infer_read_length() != record.query_length:
        raise ValueError("hard-clipped or unavailable original query axis")
    forward, phred = record.get_forward_sequence(), record.get_forward_qualities()
    if forward is None or phred is None:
        raise ValueError("BAM original bases/qualities unavailable; FASTQ cannot invent them")
    encoded = "".join(chr(q + 33) for q in phred)
    if sequence != forward or qualities != encoded:
        raise ValueError("BAM/FASTQ original-orientation bases or qualities disagree")
    return sequence, qualities


def convert_bam(bam: Path, fastq: Path, sam: Path) -> tuple[list[dict], list[dict], list[dict]]:
    import pysam  # Deliberate external tool prerequisite; no repository dependency.

    originals, witnesses, names, contigs = read_fastq(fastq), [], set(), []
    with pysam.AlignmentFile(str(bam), "rb") as source, sam.open("x", encoding="ascii", newline="\n") as output:
        contigs = [{"name": name, "length": length} for name, length in zip(source.references, source.lengths)]
        output.write(str(source.header))
        for ordinal, record in enumerate(source.fetch(until_eof=True), 1):
            if ordinal > MAX_RECORDS or record.flag & 2308 or not record.query_name or record.query_name in names:
                raise ValueError("input must contain distinct primary mapped records within budget")
            if record.query_name not in originals:
                raise ValueError("BAM record missing from paired FASTQ")
            sequence, qualities = matched_original(record, originals[record.query_name])
            output.write(record.to_string() + "\n")
            names.add(record.query_name)
            witnesses.append({"ordinal": ordinal, "unit_id": f"bam_record_{ordinal:08d}",
                "qname": record.query_name, "read_group": record.get_tag("RG") if record.has_tag("RG") else None,
                "sequence": sequence, "qualities_phred33": qualities,
                "fastq_ordinal": originals[record.query_name][0], "flag": record.flag,
                "sequence_sha256": sha_bytes(sequence.encode("ascii")),
                "quality_sha256": sha_bytes(qualities.encode("ascii"))})
    if names != set(originals):
        raise ValueError("FASTQ population differs from BAM population")
    proof = [{k: v for k, v in row.items() if k not in ("sequence", "qualities_phred33")} for row in witnesses]
    return contigs, witnesses, proof


def invoke(cli: Path, root: Path, *args: str) -> dict:
    state = root / "project_must_remain_absent.json"
    run = subprocess.run([str(cli), "--state", str(state), *args], capture_output=True,
                         text=True, encoding="utf-8", timeout=600)
    if state.exists():
        raise ValueError("file-only engine route unexpectedly wrote project state")
    if run.returncode:
        raise ValueError(f"GENtle command failed ({args[0]}): {run.stderr.strip()}")
    return json.loads(run.stdout)


def annotation_request(gtf: Path, dictionary: Path, selector: dict, reference: dict, output: Path) -> dict:
    return {"schema": "gentle.observed_transcript_annotation_request.v1", "format": "gtf",
        "gtf": bound(gtf), "dictionary": bound(dictionary), "reference_id": REFERENCE_ID,
        "reference": reference, "selector": selector,
        "source": {"scope": "producer_excerpt", "provider": "Ensembl",
                   "declared_provider_release": "116", "declared_derivation": "Glen's pinned TP53/TP73 per-gene GTF excerpt",
                   "source_uri": None, "notice": "Parent genome/source hashes are declarations; only supplied files are reverified."},
        "policy": {"transcript_feature_rows": "required", "minimum_reference_gap_bases": 1},
        "budgets": {"max_input_bytes": 16_000_000, "max_dictionary_bytes": 1_000_000,
                    "max_line_bytes": 65536, "max_lines": 100000, "max_dictionary_entries": 1000,
                    "max_indexed_genes": 10, "max_indexed_transcripts": 1000, "max_exon_rows": 10000,
                    "max_retained_bytes": 128 * 1024 * 1024, "max_work": 100_000_000,
                    "max_diagnostics": 100, "max_output_bytes": 16_000_000},
        "output_dir": str(output), "cancel_path": None}


def sam_request(sam: Path, manifest: Path, dictionary: Path, annotation: dict, selector: dict,
                reference: dict, output: Path, bam: Path, version: str) -> dict:
    return {"schema": "gentle.observed_sam_structures_request.v1", "sam": bound(sam),
        "manifest": bound(manifest), "dictionary": bound(dictionary),
        "conversion": {"kind": "bam_conversion", "bam_sha256": bound(bam)["sha256"],
                       "producer": "pysam", "version": version,
                       "arguments": ["AlignmentFile(rb)", "fetch(until_eof=True)", "AlignedSegment.to_string"],
                       "selection_scope": "Every record in supplied primary-only overlapping-gene BAM; no new alignment/filtering."},
        "reference_id": REFERENCE_ID, "reference": reference, "reference_gene": selector,
        "comparison_genes": [], "annotation_sources": [{"kind": "normalized", "selector": selector,
             "provider": "Ensembl", "annotation": annotation}], "version_policy": "strict16",
        "unpadded_reference": True, "required_flank_bases": 2,
        "text_budgets": {"max_input_bytes": 32_000_000, "max_line_bytes": 1_000_000,
            "max_lines": 100000, "max_rows": MAX_RECORDS, "max_dictionary_entries": 1000,
            "max_header_ids": 1000, "max_fields_per_line": 1035, "max_operations": 1_000_000,
            "max_sa_entries": 65536, "max_retained_bytes": 128 * 1024 * 1024, "max_work": 100_000_000},
        "kernel_budgets": {"max_queries": MAX_RECORDS, "max_records": MAX_RECORDS,
            "max_segments_per_query": 128, "max_operations": 1_000_000, "max_skips_per_chain": 1000,
            "max_chain_families": MAX_RECORDS, "max_alternatives_per_record": 128,
            "max_annotation_rows": 10000, "max_model_comparisons_per_family": 10000,
            "max_total_work": 100_000_000},
        "budgets": {"max_input_bytes": 64_000_000, "max_annotation_bytes": 8_000_000,
            "max_output_bytes": 128 * 1024 * 1024, "max_output_row_bytes": 8_000_000,
            "max_units": MAX_RECORDS, "max_retained_bytes": 128 * 1024 * 1024,
            "max_projection_work": 100_000_000}, "output_dir": str(output), "cancel_path": None}


def partition_sam(sam: Path, originals: list[dict], root: Path, size: int) -> list[tuple[Path, Path, list[dict]]]:
    """Partition complete records, preserving one explicit original-BAM ordinal per unit."""
    import pysam

    if not 1 <= size <= 32 or (len(originals) + size - 1) // size > 64:
        raise ValueError("record partition size/count exceeds explicit budget")
    parent = bound(sam)
    result = []
    with pysam.AlignmentFile(str(sam), "r") as source:
        records = iter(source.fetch(until_eof=True))
        for start in range(0, len(originals), size):
            part = root / f"part_{len(result) + 1:04d}"
            part.mkdir()
            path = part / "selected.sam"
            rows = []
            with path.open("x", encoding="ascii", newline="\n") as target:
                target.write(str(source.header))
                for local_ordinal, original in enumerate(originals[start:start + size], 1):
                    record = next(records)
                    if record.query_name != original["qname"]:
                        raise ValueError("partition order differs from original record proof")
                    target.write(record.to_string() + "\n")
                    rows.append({**original, "original_bam_ordinal": original["ordinal"], "ordinal": local_ordinal})
            write_json(part / "partition-proof.json", {"parent_full_sam": parent,
                "partition_index": len(result), "records": [{k: v for k, v in row.items()
                    if k not in ("sequence", "qualities_phred33")} for row in rows], "sam": bound(path)})
            result.append((part, path, rows))
        if next(records, None) is not None:
            raise ValueError("partitioning omitted source records")
    return result


def prepare(fixture: Path, cli: Path, output: Path, commit: str, genes: list[str], compare: bool,
            records_per_partition: int = 32) -> dict:
    import pysam

    fixture, cli, output = fixture.resolve(), cli.resolve(), output.resolve()
    if (not re.fullmatch(r"[0-9a-f]{40}", commit) or not cli.is_file()
            or len(set(genes)) != len(genes) or not 1 <= records_per_partition <= 32):
        raise ValueError("explicit fixture commit, executable and unique gene selection required")
    inventory = verify_inventory(fixture / "data")
    loci = {r["gene"]: r for r in table(fixture / "data/loci.tsv")}
    samples = table(fixture / "data/samples.tsv")
    if [(r["barcode"], r["cell_line"]) for r in samples] != [("barcode15", "SK-MEL-28"), ("barcode16", "WM793")]:
        raise ValueError("unexpected fixture sample identities")
    refs = [r for r in table(fixture / "data/source_bindings.tsv") if r["role"] == "reference_fasta"]
    if len(refs) != 1 or not re.fullmatch(r"[0-9a-f]{64}", refs[0]["sha256_or_bound_digest"]):
        raise ValueError("one explicit parent reference binding required")
    reference = {"assembly": "GRCh38", "fasta_sha256": "sha256:" + refs[0]["sha256_or_bound_digest"]}
    regions = region_sequences(fixture / "data/reference/ensembl116.tp53_tp73.genomic_regions.fa")
    output.mkdir(parents=True, exist_ok=False)
    receipt = {"schema": "gentle.cdna_fixture_preparation.v1", "fixture_commit": commit,
        "fixture_files": inventory, "cli": bound(cli), "pysam_version": pysam.__version__,
        "parent_reference": refs[0], "records_per_partition": records_per_partition,
        "genes": [], "limitations": [
            "Primary mapped overlap-selected records only; secondary, supplementary and unmapped records are absent.",
            "FASTQ was derived from archived BAM/CRAM, not independently reverified raw original FASTQ.",
            "No additional query-position exclusions are applied; no target-specific RT primer was used, but adapter/mismatch cleanliness is not inferred.",
            "RNA evidence and empirical candidate templates do not establish genomic genotype, copy number or a trans-splicing mechanism."]}
    expected_contigs = None
    for gene in genes:
        selected, locus = GENES[gene], loci[gene]
        if locus["gene_id"] != selected[1] or locus["assembly"] != reference["assembly"]:
            raise ValueError("fixture locus identity mismatch")
        start, end = int(locus["start_1based"]) - 1, int(locus["end_1based"])
        region = f"{locus['chromosome']}:{start + 1}-{end}"
        sequence = regions[region]
        if len(sequence) != end - start:
            raise ValueError("region FASTA length does not match declared genomic interval")
        root = output / gene
        root.mkdir()
        fasta = root / "locus.fa"
        with fasta.open("x", encoding="ascii", newline="\n") as handle:
            handle.write(f">{locus['chromosome']}\n{sequence}\n")
        selector = {"source_id": SOURCE_ID, "gene_id": selected[0]}
        prepared = []
        for sample in samples:
            sample_root = root / sample["barcode"]
            sample_root.mkdir()
            prefix = f"{sample['barcode']}.{gene}.primary"
            bam, fastq = [fixture / "data/reads" / (prefix + ending) for ending in (".bam", ".fastq.gz")]
            sam = sample_root / "selected.sam"
            contigs, originals, proof = convert_bam(bam, fastq, sam)
            if expected_contigs is None:
                expected_contigs = contigs
                write_json(output / "dictionary.json", {"schema": "gentle.reference_contig_dictionary.v1",
                    "reference_id": REFERENCE_ID, "assembly": reference["assembly"],
                    "reference_fasta_sha256": reference["fasta_sha256"], "contigs": contigs})
            if contigs != expected_contigs:
                raise ValueError("BAM reference dictionaries disagree")
            write_json(sample_root / "record-binding-proof.json", {"bam": bound(bam), "fastq": bound(fastq),
                "sam": bound(sam), "records": proof, "both_populations_eof_verified": True})
            parts = partition_sam(sam, originals, sample_root, records_per_partition)
            prepared.append((sample, bam, parts))
        dictionary = output / "dictionary.json"
        request = annotation_request(fixture / "data/annotation/ensembl116.tp53_tp73.gtf", dictionary,
                                     selector, reference, root / "prepared-annotation")
        path = root / "annotation.request.json"
        write_json(path, request)
        invoke(cli, output, "rna-reads", "transcript-annotation-prepare", "@" + str(path))
        annotation_report_path = root / "prepared-annotation/report.json"
        annotation_report = json.loads(annotation_report_path.read_bytes())
        if annotation_report["status"] != "completed" or annotation_report["analysis"]["outcome"] != "completed":
            raise ValueError(f"engine annotation normalization unavailable: {gene}")
        annotation = bound(root / "prepared-annotation" / annotation_report["annotation"]["path"])
        comparison_samples, counts = [], []
        for sample, bam, parts in prepared:
            sample_parts = []
            count = {"sample_id": sample["barcode"], "units": 0, "partitions": len(parts),
                "units_without_complete_search_declaration": 0,
                "units_without_complete_view_declaration": 0, "units_with_unassessed_supplied_views": 0}
            for sample_root, sam, originals in parts:
                part, part_counts = prepare_partition(cli, output, sample_root, sam, bam, originals,
                    dictionary, annotation, selector, reference, receipt["limitations"], commit, pysam.__version__)
                sample_parts.append(part)
                for name, number in part_counts.items():
                    count[name] += number
            comparison_samples.append({"sample_id": sample["barcode"], "label": sample["cell_line"],
                **sample_parts[0], "additional_partitions": sample_parts[1:],
                "partition_source_sha256": bound(bam)["sha256"]})
            counts.append(count)
        comparison = {"schema": "gentle.cdna_allele_comparison_request.v1", "gene": gene,
            "source_id": SOURCE_ID, "gene_id": selected[0],
            "reference": {"assembly": reference["assembly"], "contig": locus["chromosome"],
                "start_0based": start, "genomic_reference_sha256": reference["fasta_sha256"], "locus_fasta": bound(fasta)},
            "normalized_annotation": annotation, "samples": comparison_samples, "supplied_sites": [],
            "policy": {}, "budgets": {}}
        path = root / "comparison.request.json"
        write_json(path, comparison)
        if compare:
            result = invoke(cli, output, "allele-hash-screen", "--compare-cdna", "@" + str(path))
            write_json(root / "comparison.json", result["result"]["cdna_allele_comparison"])
        receipt["genes"].append({"gene": gene, "locus": locus, "selected_exact_gene_id": selected[0],
            "request": bound(path), "annotation_report": bound(annotation_report_path), "counts": counts,
            "comparison": bound(root / "comparison.json") if compare else None})
    receipt["content_sha256"] = sha_bytes(json.dumps(receipt, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode())
    write_json(output / "preparation.json", receipt)
    return receipt


def prepare_partition(cli: Path, output: Path, sample_root: Path, sam: Path, bam: Path,
                      originals: list[dict], dictionary: Path, annotation: dict, selector: dict,
                      reference: dict, limitations: list[str], commit: str, version: str) -> tuple[dict, dict]:
    """Use the unchanged saved-SAM producer for one bounded record partition."""
    sam_sha = bound(sam)["sha256"]
    manifest = sample_root / "units.json"
    write_json(manifest, {"schema": "gentle.observed_sam_query_units.v1", "sam_sha256": sam_sha,
        "units": [{"unit_id": row["unit_id"], "unclipped_axis_length": len(row["sequence"]),
            "rows": [{"ordinal": row["ordinal"], "qname": row["qname"], "read_group": row["read_group"]}],
            "masks": {"coordinate_convention": "zero_based_half_open_unclipped_query",
                "unclipped_axis_length": len(row["sequence"]), "intervals": [], "source_sha256": sam_sha,
                "derivation": "Explicit policy: no additional query-position exclusions; regular cDNA without gene-specific RT primer. Aligned support remains engine-derived."},
            "ends": [[], []], "search": {"search_complete": False, "views_complete": False,
                "limitations": limitations[:3], "producer": "Glen primary-only fixture via pysam",
                "version": version, "arguments": ["-F", "2308", "no MAPQ floor"],
                "selection_scope": "All supplied primary rows; original mapping alternatives and raw FASTQ are unavailable."}}
            for row in originals]})
    request = sam_request(sam, manifest, dictionary, annotation, selector, reference,
                          sample_root / "saved-sam", bam, version)
    path = sample_root / "sam.request.json"
    write_json(path, request)
    invoke(cli, output, "rna-reads", "transcript-structures-from-sam", "@" + str(path))
    report_path = sample_root / "saved-sam/report.json"
    report = json.loads(report_path.read_bytes())
    analysis_file = sample_root / "saved-sam/artifacts/analysis.json"
    analysis = json.loads(analysis_file.read_bytes())["analysis"]
    if report["status"] != "completed" or not report["source_eof_verified"]:
        raise ValueError("incomplete saved-SAM producer result")
    by_id = {row["unit_id"]: row for row in originals}
    reads = []
    for unit in analysis["units"]:
        original = by_id.pop(unit["unit"]["unit_id"])
        if (unit["unclipped_axis_length"] != len(original["sequence"])
                or any(r["qname"] != original["qname"] for r in unit["rows"])):
            raise ValueError("saved unit and original-record proof disagree")
        reads.append({"unit": unit["unit"], "read_id": original["qname"],
            "sequence": original["sequence"], "qualities_phred33": original["qualities_phred33"]})
    if by_id:
        raise ValueError("producer omitted declared units")
    calls = write_json(sample_root / "basecalls.json", {"schema": "gentle.cdna_allele_basecalls.v1",
        "analysis_context_sha256": analysis["context_sha256"], "orientation": "original_sequencing_5_to_3",
        "provenance": f"Pinned fixture {commit}; paired BAM/FASTQ exact base+Phred check, pysam {version}; see parent record-binding-proof.json and local partition-proof.json. Original pre-mapping FASTQ not reverified.",
        "reads": reads})
    return {"analysis": {"schema": "gentle.observed_sam_structures_result_binding.v1",
            "request_sha256": report["request_sha256"], "report_content_sha256": report["content_sha256"],
            "report_path": str(report_path), "outcome": "completed"}, "basecalls": calls}, {"units": len(reads),
        "units_without_complete_search_declaration": sum(not u["search_declared_complete"] for u in analysis["units"]),
        "units_without_complete_view_declaration": sum(not u["views_declared_complete"] for u in analysis["units"]),
        "units_with_unassessed_supplied_views": sum(not u["all_query_views_assessed"] for u in analysis["units"])}

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--fixture-commit", required=True)
    parser.add_argument("--gentle-cli", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--genes", nargs="+", choices=GENES, default=list(GENES))
    parser.add_argument("--compare", action="store_true", help="Also execute the unchanged default comparison policy.")
    parser.add_argument("--records-per-partition", type=int, default=32, choices=range(1, 33))
    args = parser.parse_args()
    try:
        receipt = prepare(args.fixture, args.gentle_cli, args.output, args.fixture_commit, args.genes, args.compare,
                          args.records_per_partition)
    except (OSError, ValueError, KeyError, ImportError, subprocess.SubprocessError) as exc:
        parser.exit(2, f"preparation failed: {exc}\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
