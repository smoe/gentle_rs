#!/usr/bin/env python3
"""Clip Ensembl regulatory features to selected TSS-window stretches for similarity search."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

try:
    from .prepare_regulatory_region_indexes import require
    from .prepare_transcript_promoterome import validate_promoterome, window_id
    from .prepare_tp73_cutrun_promoter_candidates import extract_fasta, load_tsv
    from .tss_regulatory_report_binding import read_locus_report, validate_locus_svg
except ImportError:
    from prepare_regulatory_region_indexes import require
    from prepare_transcript_promoterome import validate_promoterome, window_id
    from prepare_tp73_cutrun_promoter_candidates import extract_fasta, load_tsv
    from tss_regulatory_report_binding import read_locus_report, validate_locus_svg


SCHEMA = "gentle.regulatory_region_comparison_sequences.v1"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def selected_window(tss: int, strand: str, upstream: int, downstream: int) -> tuple[int, int]:
    if strand == "+":
        return tss - upstream, tss + downstream
    if strand == "-":
        return tss - downstream, tss + upstream
    raise ValueError(f"unsupported transcript strand: {strand!r}")


def connected_stretches(windows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    require(len({(row.get("chromosome"), row.get("strand")) for row in windows}) <= 1,
            "TSS stretches must share a chromosome and strand")
    stretches: list[dict[str, Any]] = []
    for window in sorted(windows, key=lambda row: (row["start_1based"], row["end_1based"])):
        if not stretches or window["start_1based"] > stretches[-1]["end_1based"] + 1:
            stretches.append({
                "start_1based": window["start_1based"],
                "end_1based": window["end_1based"],
                "tss_windows": [window],
            })
        else:
            stretches[-1]["end_1based"] = max(
                stretches[-1]["end_1based"], window["end_1based"]
            )
            stretches[-1]["tss_windows"].append(window)
    return stretches


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def assembly_forward_slice(
    sequence: str,
    window: dict[str, str],
    start_1based: int,
    end_1based: int,
) -> str:
    """Slice a genomic interval from one receipt-bound transcript-oriented promoter."""
    window_start = int(window["start_0based"])
    window_end = int(window["end_0based_exclusive"])
    start_0based = start_1based - 1
    require(window_start <= start_0based < end_1based <= window_end,
            "feature intersection is outside its selected promoter window")
    if window["strand"] == "+":
        return sequence[start_0based - window_start:end_1based - window_start]
    require(window["strand"] == "-", "invalid promoter strand")
    biological_slice = sequence[window_end - end_1based:window_end - start_0based]
    return reverse_complement(biological_slice)


def validate_selected_region(row, reference, windows, transcripts, sequences) -> None:
    """Resolve all labels against one receipt-bound promoter, not independent IDs."""
    promoter_id = row["promoterome_id"]
    require(promoter_id in windows, "selected TSS promoter is absent from the prepared promoterome")
    window = windows[promoter_id]
    extraction = row["genome_extraction"]
    start, end, tss = (int(window[key]) for key in
                       ("start_0based", "end_0based_exclusive", "tss_1based"))
    strand = window["strand"]
    require(strand in {"+", "-"} and 0 <= start < end,
            "invalid receipt-bound promoter geometry")
    expected_start, expected_end = selected_window(
        tss, strand, reference["upstream_bp"], reference["downstream_bp"])
    require((start + 1, end) == (expected_start, expected_end)
            and str(window.get("boundary_clipped", "")).lower() == "false"
            and promoter_id == window_id(window["chromosome"], start, end, strand),
            "receipt-bound promoter identity, clipping or TSS geometry mismatch")
    expected = dict(genome_id=reference["genome_id"], chromosome=window["chromosome"],
                    strand=strand, tss_1based=tss, start_1based=start + 1, end_1based=end,
                    promoter_upstream_bp=reference["upstream_bp"],
                    promoter_downstream_bp=reference["downstream_bp"])
    require(all(extraction.get(key) == value for key, value in expected.items()),
            "selected TSS geometry differs from its receipt-bound promoter")
    ids = row["transcript_ids"]
    require(ids and len(set(ids)) == len(ids)
            and sorted(extraction.get("transcript_ids", [])) == sorted(ids),
            "missing, duplicate or inconsistent selected transcript membership")
    require(all(tx in transcripts and transcripts[tx]["promoter_id"] == promoter_id
                and transcripts[tx]["gene_name"] == row["gene_query"] for tx in ids),
            "selected transcript does not belong to the selected gene/promoter")
    gene_ids = {transcripts[tx]["gene_id"] for tx in ids}
    require(len(gene_ids) == 1 and next(iter(gene_ids))
            and extraction.get("gene_id") == next(iter(gene_ids))
            and extraction.get("gene_name") == row["gene_query"],
            "selected gene identity disagrees with transcript membership")
    sequence = sequences[promoter_id]
    require(len(sequence) == end - start == row["sequence_length_bp"]
            and row.get("sequence_orientation") == "biological_5prime_to_3prime"
            and row["sequence_sha256"].removeprefix("sha256:") == sha256_bytes(sequence.encode("ascii")),
            "selected sequence digest, length or orientation mismatch")


def validate_reference_assembly(reference, assembly_id: str, catalog_path: Path | None) -> None:
    """Resolve assembly identity from the receipt-bound catalog entry, not its label."""
    require(isinstance(assembly_id, str) and assembly_id and assembly_id == assembly_id.strip(),
            "assembly identifier must be nonempty and exact")
    filename = catalog_path or reference.get("catalog_path")
    require(filename, "promoterome receipt lacks a catalog path; provide --catalog")
    payload = Path(filename).read_bytes()
    require(sha256_bytes(payload) == reference.get("catalog_sha256", "").removeprefix("sha256:"),
            "promoterome catalog hash mismatch")
    catalog = json.loads(payload)
    entry = catalog.get(reference["genome_id"]) if isinstance(catalog, dict) else None
    require(isinstance(entry, dict), "prepared genome is absent from the receipt-bound catalog")
    assemblies = []
    stem = (entry.get("ensembl_template") or {}).get("file_stem")
    if isinstance(stem, str) and "." in stem:
        assemblies.append(stem.partition(".")[2])
    if entry.get("ncbi_assembly_name"):
        assemblies.append(entry["ncbi_assembly_name"])
    require(assemblies and all(value == assembly_id for value in assemblies),
            "declared assembly disagrees with the receipt-bound prepared genome")


def validate_locus_reference(report, windows, assembly_id) -> None:
    anchor = (report.get("sequence_binding") or {}).get("genome_anchor") or {}
    chromosomes = {window["chromosome"] for window in windows}
    strands = {window["strand"] for window in windows}
    require(len(chromosomes) == len(strands) == 1, "gene has mixed chromosome/strand TSS windows")
    require(anchor.get("genome_id") == assembly_id
            and anchor.get("chromosome") == next(iter(chromosomes))
            and report.get("gene_strand") == next(iter(strands))
            and report.get("isoform_evidence", {}).get("chromosome") == anchor.get("chromosome"),
            "locus report lacks a matching genome/chromosome/strand binding; re-export it")
    require(all(anchor.get("start_1based", 0) <= row["start_1based"]
                and row["end_1based"] <= anchor.get("end_1based", 0) for row in windows),
            "selected TSS window lies outside the locus sequence binding")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-tss", type=Path, required=True)
    parser.add_argument("--locus-report", type=Path, action="append", required=True)
    parser.add_argument("--locus-svg", action="append", default=[], metavar="GENE=PATH",
                        help="Declare the original SVG for each tall-report gene before comparison")
    parser.add_argument("--promoterome", type=Path, required=True)
    parser.add_argument("--catalog", type=Path,
                        help="Relocated genome catalog; must match the promoterome receipt hash")
    parser.add_argument(
        "--assembly-id", required=True,
        help="Exact independently declared assembly identifier (for example GRCh38)",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--upstream-bp", type=int, default=500)
    parser.add_argument("--downstream-bp", type=int, default=200)
    args = parser.parse_args()
    if args.upstream_bp < 0 or args.downstream_bp < 0:
        raise SystemExit("window spans must be non-negative")

    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()),
            "Output directory must be absent or empty")
    promoterome = args.promoterome.resolve(strict=True)
    reference = validate_promoterome(promoterome)
    validate_reference_assembly(reference, args.assembly_id, args.catalog)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[1], text=True
    ).strip()
    require(args.source_revision == revision,
            "--source-revision must be the current checkout's full HEAD SHA")
    selected_bytes = args.selected_tss.read_bytes()
    selected = json.loads(selected_bytes)
    require(selected.get("schema") == SCHEMA, "unsupported selected-TSS candidate schema")
    require(selected.get("source_bindings", {}).get("promoterome_receipt_sha256", "").removeprefix("sha256:")
            == sha256_file(promoterome / "receipt.json").removeprefix("sha256:"),
            "selected TSS candidates and feature comparison use different promoterome receipts")
    window_rows = load_tsv(promoterome / "promoter_windows.tsv")
    windows = {row["promoter_id"]: row for row in window_rows}
    mappings = load_tsv(promoterome / "promoter_transcripts.tsv")
    transcripts = {row["transcript_id"]: row for row in mappings}
    require(len(windows) == len(window_rows) == reference["unique_promoter_window_count"]
            and len(transcripts) == len(mappings) == reference["included_transcript_count"],
            "duplicate or inconsistent promoter/transcript inventory")
    require(selected["regions"] and len({row["region_id"] for row in selected["regions"]})
            == len(selected["regions"]), "empty or duplicate selected regions")
    selected_promoter_ids = {row["promoterome_id"] for row in selected["regions"]}
    promoter_sequences = extract_fasta(promoterome / "promoter_windows.fa", selected_promoter_ids)
    selected_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in selected["regions"]:
        extraction = row["genome_extraction"]
        validate_selected_region(row, reference, windows, transcripts, promoter_sequences)
        start, end = selected_window(
            extraction["tss_1based"], extraction["strand"],
            args.upstream_bp, args.downstream_bp,
        )
        require(extraction["start_1based"] <= start <= end <= extraction["end_1based"],
                "requested TSS window exceeds the selected promoter")
        selected_by_gene[row["gene_query"]].append({
            "tss_1based": extraction["tss_1based"],
            "strand": extraction["strand"],
            "chromosome": extraction["chromosome"],
            "transcript_ids": row["transcript_ids"],
            "promoterome_id": row["promoterome_id"],
            "start_1based": start,
            "end_1based": end,
        })

    reports = {}
    report_hashes = {}
    for path in args.locus_report:
        report, digest = read_locus_report(path)
        gene = report["gene_symbol"]
        if gene in reports:
            raise SystemExit(f"duplicate locus report for {gene}")
        require(gene in selected_by_gene, f"unexpected locus report for {gene}")
        validate_locus_reference(report, selected_by_gene[gene], args.assembly_id)
        binding = report["ensembl_regulation"]["source_binding"]
        if binding["content_identity_verified"] is not True or binding["truncated"] is not False:
            raise SystemExit(f"unverified or truncated Ensembl evidence for {gene}")
        assembly_names = {row["assembly_name"] for row in report["ensembl_regulation"]["rows"]}
        require(assembly_names == {args.assembly_id},
                f"Ensembl feature assembly disagrees with promoterome for {gene}")
        reports[gene] = (path.resolve(), report)
        report_hashes[gene] = digest

    missing = sorted(set(selected_by_gene) - set(reports))
    if missing:
        raise SystemExit(f"missing locus reports: {missing}")

    svg_hashes = {}
    for value in args.locus_svg:
        gene, separator, filename = value.partition("=")
        require(separator and filename and gene in reports and gene not in svg_hashes,
                "--locus-svg requires a unique selected GENE=PATH")
        payload = Path(filename).read_bytes()
        validate_locus_svg(payload.decode("utf-8"), reports[gene][1])
        svg_hashes[gene] = sha256_bytes(payload)

    regions = []
    sequences: dict[str, str] = {}
    stretches_out = []
    for gene in sorted(selected_by_gene):
        report_path, report = reports[gene]
        rows = report["ensembl_regulation"]["rows"]
        gene_transcripts = {transcript for window in selected_by_gene[gene]
                            for transcript in window["transcript_ids"]}
        gene_ids = {transcripts[transcript]["gene_id"] for transcript in gene_transcripts}
        require(len(gene_ids) == 1 and all(transcripts[transcript]["gene_name"] == gene
                                           for transcript in gene_transcripts),
                f"selected transcript identity disagrees with gene {gene}")
        gene_id = next(iter(gene_ids))
        for stretch_index, stretch in enumerate(connected_stretches(selected_by_gene[gene]), 1):
            chromosome = stretch["tss_windows"][0]["chromosome"]
            require(all(window["chromosome"] == chromosome for window in stretch["tss_windows"]),
                    "connected TSS stretch crosses chromosomes")
            stretch_id = f"{gene}_tss_stretch_{stretch_index}"
            stretches_out.append({"stretch_id": stretch_id, "gene": gene, **stretch})
            for feature in rows:
                feature_start = feature["core_genomic_start_1based"]
                feature_end = feature["core_genomic_end_1based"]
                start = max(stretch["start_1based"], feature_start)
                end = min(stretch["end_1based"], feature_end)
                if start > end:
                    continue
                region_id = f"{gene}_{stretch_id}_{feature['feature_id']}"
                containing = [window for window in stretch["tss_windows"]
                              if int(windows[window["promoterome_id"]]["start_0based"]) <= start - 1
                              and end <= int(windows[window["promoterome_id"]]["end_0based_exclusive"])]
                require(containing, f"no receipt-bound promoter contains {region_id}")
                extracted = {
                    assembly_forward_slice(
                        promoter_sequences[window["promoterome_id"]],
                        windows[window["promoterome_id"]], start, end,
                    )
                    for window in containing
                }
                require(len(extracted) == 1,
                        f"receipt-bound promoter sequences disagree for {region_id}")
                sequence = next(iter(extracted))
                if len(sequence) != end - start + 1:
                    raise SystemExit(f"reference extraction length mismatch for {region_id}")
                sequences[region_id] = sequence
                regions.append({
                    "region_id": region_id,
                    "input_kind": "canonical_region",
                    "comparison_class": "ensembl_regulation_tss_window_intersection",
                    "gene_query": gene,
                    "transcript_ids": sorted({transcript for window in stretch["tss_windows"]
                                               for transcript in window["transcript_ids"]}),
                    "sequence_length_bp": len(sequence),
                    "sequence_sha256": f"sha256:{sha256_bytes(sequence.encode('ascii'))}",
                    "sequence_orientation": "assembly_reference_forward",
                    "genome_extraction": {
                        "genome_id": reference["genome_id"],
                        "chromosome": chromosome,
                        "start_1based": start,
                        "end_1based": end,
                        "gene_id": gene_id,
                        "gene_name": gene,
                    },
                    "source_region": {
                        "region_id": feature["feature_id"],
                        "purpose": "regulatory_feature_similarity_within_selected_tss_window",
                        "interval": {
                            "coordinate_convention": "zero_based_half_open",
                            "reference": {
                                "assembly_name": feature["assembly_name"],
                                "assembly_accession": feature["assembly_accession"],
                                "contig_name": chromosome,
                            },
                            "start_0based": start - 1,
                            "end_0based_exclusive": end,
                            "strand": "unstranded",
                        },
                        "feature_type": feature["feature_type"],
                        "core_genomic_start_1based": feature_start,
                        "core_genomic_end_1based": feature_end,
                        "clipped_to_tss_stretch": start != feature_start or end != feature_end,
                        "stretch_id": stretch_id,
                        "tss_windows": stretch["tss_windows"],
                        "canonical_feature_url": feature["canonical_feature_url"],
                        "source_id": feature["source_id"],
                        "annotation_release": feature["annotation_release"],
                        "source_report_sha256": f"sha256:{report_hashes[gene]}",
                    },
                })

    sequence_classes: dict[str, list[str]] = defaultdict(list)
    for region_id, sequence in sequences.items():
        sequence_classes[sequence].append(region_id)
    equivalence = []
    fasta_rows = []
    for sequence, members in sorted(sequence_classes.items(), key=lambda item: sorted(item[1])):
        members = sorted(members)
        representative = members[0]
        equivalence.append({
            "sequence_equivalence_class_id": f"region_seq_eq_{sha256_bytes(sequence.encode('ascii'))[:16]}",
            "representative_region_id": representative,
            "member_region_ids": members,
            "sequence_length_bp": len(sequence),
            "sequence_sha256": f"sha256:{sha256_bytes(sequence.encode('ascii'))}",
        })
        fasta_rows.append((representative, sequence))

    payload = {
        "schema": SCHEMA,
        "dataset_id": "selected_gene_ensembl_regulation_intersections_tss_500_200_v1",
        "source_revision": revision,
        "window_policy": {
            "coordinate_basis": "transcript_oriented_tss",
            "upstream_bp": args.upstream_bp,
            "downstream_bp": args.downstream_bp,
            "feature_sequence_policy": "intersection_with_connected_selected_tss_windows",
        },
        "stretches": stretches_out,
        "regions": regions,
        "sequence_equivalence_classes": equivalence,
        "source_bindings": {
            "assembly_id": args.assembly_id,
            "selected_tss_sha256": f"sha256:{sha256_bytes(selected_bytes)}",
            "promoterome_receipt_sha256": sha256_file(promoterome / "receipt.json"),
            "promoterome_fasta_sha256": reference["artifacts"]["promoter_windows.fa"],
            "producer_sha256": sha256_file(Path(__file__)),
            "locus_reports": {
                gene: f"sha256:{digest}" for gene, digest in sorted(report_hashes.items())
            },
            "locus_svgs": {gene: f"sha256:{digest}" for gene, digest in sorted(svg_hashes.items())},
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "candidate_regions.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (output / "candidate_regions.fa").open("w", encoding="ascii") as handle:
        for representative, sequence in fasta_rows:
            handle.write(f">{representative}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")
    print(json.dumps({
        "stretches": len(stretches_out),
        "regions": len(regions),
        "shorter_than_40_bp": sum(len(sequence) < 40 for sequence in sequences.values()),
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
