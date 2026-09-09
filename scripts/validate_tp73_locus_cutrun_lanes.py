#!/usr/bin/env python3
"""Fail closed unless the three TP73 locus reports retain all known BigWig lanes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import pyBigWig


SCHEMA = "gentle.tp73_locus_cutrun_lane_validation.v1"
LOCUS_SCHEMA = "gentle.gene_locus_evidence_display.v1"
REQUEST_SCHEMA = "gentle.gene_locus_evidence_preparation_request.v1"
EXPECTED_GENES = {"CD44", "TGFB1", "SERPINE1"}
EXPECTED_LANES = 12


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def declared_files(values: list[str], option: str) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        gene, separator, filename = value.partition("=")
        require(separator == "=" and gene in EXPECTED_GENES and gene not in result,
                f"{option} requires one unique GENE=PATH for each expected gene")
        result[gene] = Path(filename).resolve(strict=True)
    require(set(result) == EXPECTED_GENES,
            f"{option} must declare exactly {sorted(EXPECTED_GENES)}")
    return result


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"expected a JSON object: {path}")
    return value


def compatible_chromosome(bigwig: pyBigWig.pyBigWig, declared: str) -> str:
    inventory = bigwig.chroms()
    candidates = [declared]
    if declared.startswith("chr"):
        candidates.append(declared[3:])
    else:
        candidates.append(f"chr{declared}")
    matches = [candidate for candidate in candidates if candidate in inventory]
    require(len(matches) == 1,
            f"BigWig must contain exactly one compatible chromosome for {declared}")
    return matches[0]


def expected_local(anchor: dict[str, Any], strand: str, start: int, end: int) -> tuple[int, int]:
    if strand == "+":
        return start - anchor["start_1based"] + 1, end - anchor["start_1based"] + 1
    require(strand == "-", "gene strand must be + or -")
    return anchor["end_1based"] - end + 1, anchor["end_1based"] - start + 1


def validate_gene(gene: str, report_path: Path, request_path: Path,
                  assembly_id: str) -> dict[str, Any]:
    report = load_json(report_path)
    request = load_json(request_path)
    require(report.get("schema") == LOCUS_SCHEMA and report.get("gene_symbol") == gene,
            f"{gene}: incompatible locus report")
    require(request.get("schema") == REQUEST_SCHEMA and request.get("gene_query") == gene,
            f"{gene}: incompatible preparation request")
    require(request.get("assembly") == assembly_id,
            f"{gene}: request assembly does not equal {assembly_id}")
    anchor = (report.get("sequence_binding") or {}).get("genome_anchor") or {}
    require(anchor.get("genome_id") == assembly_id,
            f"{gene}: locus anchor assembly does not equal {assembly_id}")
    chromosome = anchor.get("chromosome")
    anchor_start = anchor.get("start_1based")
    anchor_end = anchor.get("end_1based")
    strand = report.get("gene_strand")
    require(isinstance(chromosome, str) and chromosome
            and isinstance(anchor_start, int) and isinstance(anchor_end, int)
            and 1 <= anchor_start <= anchor_end and strand in {"+", "-"},
            f"{gene}: invalid genome anchor")

    tracks = request.get("local_tracks")
    require(isinstance(tracks, list) and len(tracks) == EXPECTED_LANES,
            f"{gene}: expected exactly {EXPECTED_LANES} requested local tracks")
    tracks_by_source: dict[str, dict[str, Any]] = {}
    for track in tracks:
        source_id = track.get("source_id")
        require(track.get("source_kind") == "big_wig" and isinstance(source_id, str)
                and source_id and source_id not in tracks_by_source,
                f"{gene}: local tracks require unique BigWig source IDs")
        require(track.get("assembly") == assembly_id,
                f"{gene}/{source_id}: track assembly does not equal {assembly_id}")
        tracks_by_source[source_id] = track

    lanes = [lane for group in report.get("occupancy_groups", [])
             for lane in group.get("lanes", [])]
    require(len(lanes) == EXPECTED_LANES,
            f"{gene}: expected exactly {EXPECTED_LANES} rendered lanes")
    lanes_by_source = {lane.get("source_id"): lane for lane in lanes}
    require(len(lanes_by_source) == EXPECTED_LANES
            and set(lanes_by_source) == set(tracks_by_source),
            f"{gene}: rendered/requested source IDs differ or are duplicated")

    lane_receipts = []
    for source_id, track in tracks_by_source.items():
        lane = lanes_by_source[source_id]
        lane_data = lane.get("lane") or {}
        intervals = lane_data.get("intervals")
        source = Path(track.get("path", "")).resolve(strict=True)
        source_digest = sha256(source)
        require(lane.get("source_sha256", "").removeprefix("sha256:") == source_digest,
                f"{gene}/{source_id}: source hash mismatch")
        require(lane.get("source_assembly") == assembly_id,
                f"{gene}/{source_id}: rendered assembly does not equal {assembly_id}")
        require(lane.get("state") == "available",
                f"{gene}/{source_id}: missing intervals are not measured zero")
        require(lane_data.get("track_name") == track.get("track_name"),
                f"{gene}/{source_id}: track-name mismatch")
        require(isinstance(intervals, list) and intervals
                and lane_data.get("interval_count") == len(intervals),
                f"{gene}/{source_id}: lane must contain counted compatible intervals")

        try:
            bigwig = pyBigWig.open(str(source))
        except RuntimeError as error:
            raise SystemExit(f"{gene}/{source_id}: cannot open BigWig: {error}") from error
        require(bigwig is not None and bigwig.isBigWig(),
                f"{gene}/{source_id}: source is not BigWig")
        try:
            native_chromosome = compatible_chromosome(bigwig, chromosome)
            native_intervals = bigwig.intervals(
                native_chromosome, anchor_start - 1, anchor_end
            ) or ()
        finally:
            bigwig.close()
        require(bool(native_intervals),
                f"{gene}/{source_id}: BigWig has no native overlap with the locus anchor")

        for interval in intervals:
            start = interval.get("genomic_start_1based")
            end = interval.get("genomic_end_1based")
            require(isinstance(start, int) and isinstance(end, int)
                    and anchor_start <= start <= end <= anchor_end,
                    f"{gene}/{source_id}: rendered interval lies outside the locus anchor")
            local_start, local_end = expected_local(anchor, strand, start, end)
            require(interval.get("local_start_1based") == local_start
                    and interval.get("local_end_1based") == local_end,
                    f"{gene}/{source_id}: local/genomic coordinates disagree")
            require(math.isfinite(float(interval.get("score"))),
                    f"{gene}/{source_id}: non-finite rendered score")
        lane_receipts.append({
            "source_id": source_id,
            "track_name": track["track_name"],
            "source_path": str(source),
            "source_sha256": f"sha256:{source_digest}",
            "source_assembly": assembly_id,
            "source_chromosome": native_chromosome,
            "source_overlapping_interval_count": len(native_intervals),
            "rendered_interval_count": len(intervals),
        })
    return {
        "gene": gene,
        "chromosome": chromosome,
        "strand": strand,
        "anchor_start_1based": anchor_start,
        "anchor_end_1based": anchor_end,
        "request_path": str(request_path),
        "request_sha256": f"sha256:{sha256(request_path)}",
        "report_path": str(report_path),
        "report_sha256": f"sha256:{sha256(report_path)}",
        "lane_count": len(lane_receipts),
        "rendered_interval_count": sum(row["rendered_interval_count"] for row in lane_receipts),
        "lanes": lane_receipts,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-id", required=True)
    parser.add_argument("--report", action="append", required=True, metavar="GENE=PATH")
    parser.add_argument("--request", action="append", required=True, metavar="GENE=PATH")
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.assembly_id.strip() == args.assembly_id and args.assembly_id,
            "--assembly-id must be a non-empty exact identifier")
    repository = Path(__file__).resolve().parents[1]
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repository, text=True
    ).strip()
    require(args.source_revision == revision,
            "--source-revision must be the current checkout's full HEAD SHA")
    reports = declared_files(args.report, "--report")
    requests = declared_files(args.request, "--request")
    require(not args.output.exists(), f"Refusing to overwrite output: {args.output}")
    genes = [validate_gene(gene, reports[gene], requests[gene], args.assembly_id)
             for gene in sorted(EXPECTED_GENES)]
    receipt = {
        "schema": SCHEMA,
        "source_revision": revision,
        "assembly_id": args.assembly_id,
        "requirements": {
            "expected_genes": sorted(EXPECTED_GENES),
            "expected_lanes_per_gene": EXPECTED_LANES,
            "missing_intervals_are_zero": False,
        },
        "genes": genes,
        "total_lane_count": sum(row["lane_count"] for row in genes),
        "total_rendered_interval_count": sum(row["rendered_interval_count"] for row in genes),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "total_lane_count": receipt["total_lane_count"],
                      "total_rendered_interval_count": receipt["total_rendered_interval_count"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
