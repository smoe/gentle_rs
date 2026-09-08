#!/usr/bin/env python3
"""Compare candidate regions to a transcript-linked promoterome background."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

try:
    from .prepare_regulatory_region_indexes import require, run_command, sha256_file, write_json
except ImportError:
    from prepare_regulatory_region_indexes import (  # type: ignore[no-redef]
        require, run_command, sha256_file, write_json,
    )


SCHEMA = "gentle.candidate_promoterome_comparison.v1"
BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].strip().split()[0])
    return ids


def expand_equivalent_query_rows(rows: list[dict[str, Any]],
                                 equivalence_classes: list[dict[str, Any]],
                                 indexed_query_ids: set[str]) -> list[dict[str, Any]]:
    members_by_representative: dict[str, list[str]] = {}
    for row in equivalence_classes:
        representative = row.get("representative_region_id")
        members = row.get("member_region_ids")
        require(isinstance(representative, str) and isinstance(members, list)
                and all(isinstance(member, str) for member in members),
                "invalid candidate sequence-equivalence class")
        members_by_representative[representative] = members
    require(indexed_query_ids == set(members_by_representative),
            "query FASTA IDs do not match sequence-equivalence representatives")
    expanded: list[dict[str, Any]] = []
    for row in rows:
        for member in members_by_representative.get(row["qseqid"], []):
            expanded.append({**row, "qseqid": member,
                             "indexed_representative_query_id": row["qseqid"]})
    return expanded


def normalized_contig(value: str) -> str:
    return value[3:] if value.lower().startswith("chr") else value


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def query_interval(row: dict[str, Any]) -> tuple[str, int, int] | None:
    if row.get("input_kind") == "canonical_region":
        interval = row.get("source_region", {}).get("interval", {})
        reference = interval.get("reference", {})
        chromosome = reference.get("contig_name")
        start = interval.get("start_0based")
        end = interval.get("end_0based_exclusive")
    else:
        extraction = row.get("genome_extraction", {})
        chromosome = extraction.get("chromosome")
        start_1based = extraction.get("start_1based")
        end = extraction.get("end_1based")
        start = start_1based - 1 if isinstance(start_1based, int) else None
    if isinstance(chromosome, str) and isinstance(start, int) and isinstance(end, int):
        return chromosome, start, end
    return None


def coverage_segments(target_intervals: dict[str, list[tuple[int, int]]],
                      genes_by_target: dict[str, set[str]],
                      transcripts_by_target: dict[str, set[str]]) -> list[dict[str, int]]:
    events: dict[int, dict[str, set[str]]] = defaultdict(lambda: {"add": set(), "remove": set()})
    for target, intervals in target_intervals.items():
        for start, end in merge_intervals(intervals):
            events[start]["add"].add(target)
            events[end]["remove"].add(target)
    active_targets: set[str] = set()
    gene_counts: Counter[str] = Counter()
    transcript_counts: Counter[str] = Counter()
    segments: list[dict[str, int]] = []
    previous: int | None = None
    for position in sorted(events):
        if previous is not None and previous < position and active_targets:
            segments.append({
                "query_start_0based": previous, "query_end_0based_exclusive": position,
                "distinct_promoter_windows": len(active_targets),
                "distinct_genes": len(gene_counts), "distinct_transcripts": len(transcript_counts),
            })
        for target in events[position]["remove"]:
            if target in active_targets:
                active_targets.remove(target)
                for gene in genes_by_target[target]:
                    gene_counts[gene] -= 1
                    if gene_counts[gene] == 0:
                        del gene_counts[gene]
                for transcript in transcripts_by_target[target]:
                    transcript_counts[transcript] -= 1
                    if transcript_counts[transcript] == 0:
                        del transcript_counts[transcript]
        for target in events[position]["add"]:
            if target not in active_targets:
                active_targets.add(target)
                gene_counts.update(genes_by_target[target])
                transcript_counts.update(transcripts_by_target[target])
        previous = position
    return segments


def parse_blast(path: Path, min_length: int, min_identity: float,
                max_evalue: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            require(len(fields) == len(BLAST_COLUMNS), f"invalid BLAST row in {path}")
            row = dict(zip(BLAST_COLUMNS, fields))
            length = int(row["length"])
            identity = float(row["pident"])
            evalue = float(row["evalue"])
            if length < min_length or identity < min_identity or evalue > max_evalue:
                continue
            rows.append({**row, "length": length, "pident": identity,
                         "qstart": int(row["qstart"]), "qend": int(row["qend"]),
                         "evalue": evalue, "bitscore": float(row["bitscore"])})
    return rows


def summarize_task(rows: list[dict[str, Any]], queries: dict[str, dict[str, Any]],
                   windows: dict[str, dict[str, str]], mappings: list[dict[str, str]],
                   matches_path: Path) -> list[dict[str, Any]]:
    genes_by_target: dict[str, set[str]] = defaultdict(set)
    transcripts_by_target: dict[str, set[str]] = defaultdict(set)
    names_by_target: dict[str, set[str]] = defaultdict(set)
    for mapping in mappings:
        genes_by_target[mapping["promoter_id"]].add(mapping["gene_id"])
        transcripts_by_target[mapping["promoter_id"]].add(mapping["transcript_id"])
        if mapping["gene_name"]:
            names_by_target[mapping["promoter_id"]].add(mapping["gene_name"])

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        require(row["qseqid"] in queries, f"unknown query id {row['qseqid']}")
        require(row["sseqid"] in windows, f"unknown promoterome target id {row['sseqid']}")
        grouped[(row["qseqid"], row["sseqid"])].append(row)

    summaries: dict[str, dict[str, Any]] = {}
    target_intervals_by_query: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(dict)
    match_rows: list[dict[str, Any]] = []
    for (query_id, target_id), hsps in sorted(grouped.items()):
        query_row = queries[query_id]
        target = windows[target_id]
        qinterval = query_interval(query_row)
        same_locus = False
        if (qinterval is not None
                and normalized_contig(qinterval[0]) == normalized_contig(target["chromosome"])):
            same_locus = qinterval[1] < int(target["end_0based_exclusive"]) and int(target["start_0based"]) < qinterval[2]
        intervals = merge_intervals([
            (min(row["qstart"], row["qend"]) - 1, max(row["qstart"], row["qend"]))
            for row in hsps
        ])
        aligned_bp = sum(end - start for start, end in intervals)
        query_length = int(query_row["sequence_length_bp"])
        match_rows.append({
            "query_id": query_id, "promoter_id": target_id,
            "same_locus_overlap": same_locus, "aligned_query_bp": aligned_bp,
            "aligned_query_fraction": aligned_bp / query_length,
            "query_intervals_0based_half_open": ";".join(f"{a}-{b}" for a, b in intervals),
            "best_identity_pct": max(row["pident"] for row in hsps),
            "best_bitscore": max(row["bitscore"] for row in hsps),
            "chromosome": target["chromosome"], "tss_1based": target["tss_1based"],
            "strand": target["strand"],
            "gene_ids": ";".join(sorted(genes_by_target[target_id])),
            "gene_names": ";".join(sorted(names_by_target[target_id])),
            "transcript_ids": ";".join(sorted(transcripts_by_target[target_id])),
        })
        if not same_locus:
            target_intervals_by_query[query_id][target_id] = intervals

    columns = ["query_id", "promoter_id", "same_locus_overlap", "aligned_query_bp",
               "aligned_query_fraction",
               "query_intervals_0based_half_open", "best_identity_pct", "best_bitscore",
               "chromosome", "tss_1based", "strand", "gene_ids", "gene_names", "transcript_ids"]
    with matches_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(match_rows)

    for query_id, query_row in sorted(queries.items()):
        all_matches = [row for row in match_rows if row["query_id"] == query_id]
        other_matches = [row for row in all_matches if not row["same_locus_overlap"]]
        other_targets = {row["promoter_id"] for row in other_matches}
        coverage_tiers: dict[str, dict[str, int]] = {}
        for threshold in (0.25, 0.5, 0.8):
            tier_targets = {
                row["promoter_id"] for row in other_matches
                if row["aligned_query_fraction"] >= threshold
            }
            coverage_tiers[f"{threshold:.2f}"] = {
                "distinct_promoter_windows": len(tier_targets),
                "distinct_genes": len({gene for target in tier_targets for gene in genes_by_target[target]}),
                "distinct_transcripts": len({tx for target in tier_targets for tx in transcripts_by_target[target]}),
            }
        summaries[query_id] = {
            "query_id": query_id,
            "comparison_class": query_row.get("comparison_class"),
            "query_length_bp": query_row.get("sequence_length_bp"),
            "all_matches": {
                "distinct_promoter_windows": len(all_matches),
                "distinct_genes": len({gene for row in all_matches for gene in row["gene_ids"].split(";") if gene}),
                "distinct_transcripts": len({tx for row in all_matches for tx in row["transcript_ids"].split(";") if tx}),
            },
            "other_promoters": {
                "distinct_promoter_windows": len(other_targets),
                "distinct_genes": len({gene for target in other_targets for gene in genes_by_target[target]}),
                "distinct_transcripts": len({tx for target in other_targets for tx in transcripts_by_target[target]}),
            },
            "other_promoters_by_min_query_coverage": coverage_tiers,
            "recurrent_query_segments": coverage_segments(
                target_intervals_by_query.get(query_id, {}), genes_by_target, transcripts_by_target
            ),
        }
    return list(summaries.values())


def compare(args: argparse.Namespace) -> None:
    candidates_json = args.candidates_json.resolve(strict=True)
    query_fasta = args.query_fasta.resolve(strict=True)
    promoterome = args.promoterome.resolve(strict=True)
    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()), "Output directory must be absent or empty")
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    logs.mkdir()
    candidate_payload = json.loads(candidates_json.read_text(encoding="utf-8"))
    queries = {row["region_id"]: row for row in candidate_payload.get("regions", [])}
    require(queries, "candidate metadata contains no regions")
    query_fasta_ids = fasta_ids(query_fasta)
    windows = {row["promoter_id"]: row for row in load_tsv(promoterome / "promoter_windows.tsv")}
    mappings = load_tsv(promoterome / "promoter_transcripts.tsv")
    require(windows and mappings, "promoterome metadata is empty")
    blast_db = promoterome / "indexes" / "promoterome"
    receipts: list[dict[str, Any]] = []
    task_summaries: dict[str, Any] = {}
    outfmt = "6 " + " ".join(BLAST_COLUMNS)
    for task in args.task:
        require(task in {"megablast", "blastn", "dc-megablast"}, f"unsupported task '{task}'")
        raw_path = output / f"hits.{task}.tsv"
        run_command([
            args.blastn, "-task", task, "-query", str(query_fasta), "-db", str(blast_db),
            "-dust", "yes", "-soft_masking", "true", "-max_target_seqs", str(args.max_target_seqs),
            "-max_hsps", str(args.max_hsps), "-outfmt", outfmt, "-out", str(raw_path),
        ], cwd=output, timeout=args.timeout, log_dir=logs,
            label=f"blast_{task}", receipts=receipts)
        filtered = parse_blast(raw_path, args.min_alignment_bp, args.min_identity_pct, args.max_evalue)
        filtered = expand_equivalent_query_rows(
            filtered, candidate_payload.get("sequence_equivalence_classes", []), query_fasta_ids
        )
        matches_path = output / f"matches.{task}.tsv"
        query_summaries = summarize_task(filtered, queries, windows, mappings, matches_path)
        for query_summary in query_summaries:
            query_summary["target_cap_reached"] = (
                query_summary["all_matches"]["distinct_promoter_windows"] >= args.max_target_seqs
            )
        task_summaries[task] = {
            "filtered_hsp_count": len(filtered),
            "queries": query_summaries,
            "raw_hits_path": raw_path.name, "matches_path": matches_path.name,
        }
    write_json(output / "comparison.json", {
        "schema": SCHEMA,
        "candidate_metadata_path": str(candidates_json),
        "candidate_metadata_sha256": sha256_file(candidates_json),
        "query_fasta_path": str(query_fasta), "query_fasta_sha256": sha256_file(query_fasta),
        "promoterome_path": str(promoterome),
        "promoterome_receipt_sha256": sha256_file(promoterome / "receipt.json"),
        "thresholds": {"min_alignment_bp": args.min_alignment_bp,
                       "min_identity_pct": args.min_identity_pct, "max_evalue": args.max_evalue,
                       "max_target_seqs": args.max_target_seqs, "max_hsps_per_target": args.max_hsps},
        "counting_policy": (
            "Frequency is reported separately for distinct genomic promoter windows, genes, and transcripts. "
            "Targets overlapping the query's own genomic locus are retained but excluded from other_promoters."
        ),
        "tasks": task_summaries, "commands": receipts,
        "non_claims": [
            "A recurrent sequence segment is not proof of promoter activity or functional interchangeability.",
            "Absence of a reported hit is conditional on the declared task, thresholds, prepared annotation, and target cap.",
            "Motif/module similarity and reporter sufficiency require separate evidence and experimental contrasts.",
        ],
    })
    print(json.dumps({"status": "ok", "queries": len(queries), "tasks": args.task,
                      "output": str(output / "comparison.json")}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates-json", type=Path, required=True)
    parser.add_argument("--query-fasta", type=Path, required=True)
    parser.add_argument("--promoterome", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--task", action="append")
    parser.add_argument("--min-alignment-bp", type=int, default=40)
    parser.add_argument("--min-identity-pct", type=float, default=80.0)
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    parser.add_argument("--max-target-seqs", type=int, default=1_000_000)
    parser.add_argument("--max-hsps", type=int, default=10)
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--timeout", type=int, default=7200)
    args = parser.parse_args()
    args.task = args.task or ["megablast", "blastn", "dc-megablast"]
    require(args.min_alignment_bp > 0 and 0 <= args.min_identity_pct <= 100
            and args.max_evalue >= 0 and args.max_target_seqs > 0 and args.max_hsps > 0
            and args.timeout > 0, "invalid comparison thresholds")
    compare(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
