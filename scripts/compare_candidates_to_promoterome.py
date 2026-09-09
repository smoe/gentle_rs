#!/usr/bin/env python3
"""Compare candidate regions to a transcript-linked promoterome background."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any

try:
    from .prepare_regulatory_region_indexes import require, run_command, sha256_file, write_json
    from .prepare_transcript_promoterome import validate_promoterome
except ImportError:
    from prepare_regulatory_region_indexes import (  # type: ignore[no-redef]
        require, run_command, sha256_file, write_json,
    )
    from prepare_transcript_promoterome import validate_promoterome


SCHEMA = "gentle.candidate_promoterome_comparison.v1"
COUNTING_POLICY = "exclude_overlapping_or_shared_gene_windows.v1"
BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def same_digest(left: str, right: str) -> bool:
    """Existing preparation and publication receipts use both SHA-256 spellings."""
    return (isinstance(left, str) and isinstance(right, str)
            and left.removeprefix("sha256:") == right.removeprefix("sha256:"))


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


def validate_candidate_sequences(payload: dict[str, Any], path: Path) -> dict[str, str]:
    """Bind every region and equivalence-class member to the supplied FASTA."""
    require(payload.get("schema") == "gentle.regulatory_region_comparison_sequences.v1",
            "unsupported candidate sequence schema")
    sequences: dict[str, str] = {}
    current = None
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            require(current not in sequences, "duplicate candidate FASTA ID")
            sequences[current] = ""
        elif line.strip():
            require(current is not None, "candidate FASTA sequence precedes its header")
            sequences[current] += line.strip()
    regions = {row["region_id"]: row for row in payload["regions"]}
    require(regions and len(regions) == len(payload["regions"]), "missing or duplicate candidate regions")
    representatives: dict[str, str] = {}
    indexed = set()
    for group in payload["sequence_equivalence_classes"]:
        representative = group["representative_region_id"]
        require(representative in sequences and representative not in indexed,
                "missing or duplicate equivalence representative")
        indexed.add(representative)
        sequence = sequences[representative]
        digest = "sha256:" + hashlib.sha256(sequence.encode("ascii")).hexdigest()
        require(sequence and set(sequence.upper()) <= set("ACGTRYSWKMBDHVN"), "invalid candidate DNA")
        require(representative in group["member_region_ids"], "representative is not a class member")
        for member in group["member_region_ids"]:
            require(member in regions and member not in representatives, "invalid or repeated equivalence member")
            representatives[member] = representative
            require(regions[member]["sequence_length_bp"] == len(sequence)
                    and regions[member]["sequence_sha256"] == digest, "candidate sequence hash/length mismatch")
        require(group["sequence_length_bp"] == len(sequence) and group["sequence_sha256"] == digest,
                "equivalence-class sequence hash/length mismatch")
    require(indexed == set(sequences) and set(representatives) == set(regions),
            "candidate FASTA/classes do not cover exactly the declared regions")
    return representatives


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


def iter_blast_rows(path: Path):
    """Read normalized HSPs; callers apply the same declared filter."""
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
            require(length > 0 and math.isfinite(identity) and 0 <= identity <= 100
                    and math.isfinite(evalue) and evalue >= 0, "invalid BLAST scores")
            yield {**row, "length": length, "pident": identity,
                         "qstart": int(row["qstart"]), "qend": int(row["qend"]),
                         "sstart": int(row["sstart"]), "send": int(row["send"]),
                         "evalue": evalue, "bitscore": float(row["bitscore"])}


def hsp_passes(row: dict[str, Any], min_length: int, min_identity: float, max_evalue: float) -> bool:
    return row["length"] >= min_length and row["pident"] >= min_identity and row["evalue"] <= max_evalue


def parse_blast(path: Path, min_length: int, min_identity: float,
                max_evalue: float, *, limits: tuple[int, int] | None = None,
                audit: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    rows = []
    raw_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in iter_blast_rows(path):
        if limits is not None:
            raw_counts[row["qseqid"]][row["sseqid"]] += 1
        if hsp_passes(row, min_length, min_identity, max_evalue):
            rows.append(row)
    if limits is not None and audit is not None:
        for query, counts in raw_counts.items():
            target_hit = len(counts) >= limits[0]
            hsp_hit = max(counts.values()) >= limits[1]
            audit[query] = {"raw_target_count": len(counts), "max_hsps_observed": max(counts.values()),
                            "target_cap_reached": target_hit, "hsp_cap_reached": hsp_hit,
                            "counts_are_lower_bounds": target_hit or hsp_hit}
    return rows


def query_gene_ids(query: dict[str, Any], mappings: list[dict[str, str]]) -> set[str]:
    extraction = query.get("genome_extraction", {})
    explicit = extraction.get("gene_id")
    ids = {explicit} if explicit else set()
    transcripts = set(query.get("transcript_ids", extraction.get("transcript_ids", [])))
    matched = [row for row in mappings if row["transcript_id"] in transcripts]
    require({row["transcript_id"] for row in matched} == transcripts,
            "query transcripts cannot all be resolved in the promoterome")
    transcript_genes = {row["gene_id"] for row in matched if row["gene_id"]}
    require(not explicit or explicit in {row["gene_id"] for row in mappings},
            "query gene ID is absent from the promoterome")
    require(not explicit or not transcripts or explicit in transcript_genes,
            "query gene ID disagrees with its transcripts")
    ids.update(transcript_genes)
    if not ids and query.get("gene_query"):
        ids.update(row["gene_id"] for row in mappings
                   if row["gene_name"] == query["gene_query"] and row["gene_id"])
        require(ids, "query gene cannot be resolved in the promoterome")
    return ids


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
    query_genes = {key: query_gene_ids(query, mappings) for key, query in queries.items()}

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
        require(all(0 <= start < end <= query_length for start, end in intervals),
                "BLAST query interval exceeds candidate length")
        same_gene = bool(query_genes[query_id] & genes_by_target[target_id])
        match_rows.append({
            "query_id": query_id, "promoter_id": target_id,
            "same_locus_overlap": same_locus, "same_gene": same_gene, "aligned_query_bp": aligned_bp,
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
        if not same_locus and not same_gene:
            target_intervals_by_query[query_id][target_id] = intervals

    columns = ["query_id", "promoter_id", "same_locus_overlap", "same_gene", "aligned_query_bp",
               "aligned_query_fraction",
               "query_intervals_0based_half_open", "best_identity_pct", "best_bitscore",
               "chromosome", "tss_1based", "strand", "gene_ids", "gene_names", "transcript_ids"]
    with matches_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(match_rows)

    for query_id, query_row in sorted(queries.items()):
        all_matches = [row for row in match_rows if row["query_id"] == query_id]
        other_matches = [row for row in all_matches if not row["same_locus_overlap"] and not row["same_gene"]]
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
            "query_gene_ids": sorted(query_genes[query_id]),
            "other_gene_exclusion_verified": bool(query_genes[query_id]),
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
    reference = validate_promoterome(promoterome, include_blast=True)
    candidate_digest = sha256_file(candidates_json)
    fasta_digest = sha256_file(query_fasta)
    reference_digest = sha256_file(promoterome / "receipt.json")
    candidate_payload = json.loads(candidates_json.read_text(encoding="utf-8"))
    representatives = validate_candidate_sequences(candidate_payload, query_fasta)
    bound_reference = candidate_payload.get("source_bindings", {}).get("promoterome_receipt_sha256")
    require(bound_reference is None or same_digest(bound_reference, reference_digest),
            "candidate was prepared against a different promoterome receipt")
    queries = {row["region_id"]: row for row in candidate_payload.get("regions", [])}
    require(queries, "candidate metadata contains no regions")
    for query in queries.values():
        genome = query.get("genome_extraction", {}).get("genome_id")
        require(genome is None or genome == reference["genome_id"], "candidate/reference genome mismatch")
    query_fasta_ids = set(representatives.values())
    windows = {row["promoter_id"]: row for row in load_tsv(promoterome / "promoter_windows.tsv")}
    mappings = load_tsv(promoterome / "promoter_transcripts.tsv")
    require(windows and mappings, "promoterome metadata is empty")
    require(len(windows) == reference["unique_promoter_window_count"]
            and len(mappings) == reference["included_transcript_count"], "reference inventory count mismatch")
    for query in queries.values():
        query_gene_ids(query, mappings)
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    logs.mkdir()
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
        audit: dict[str, Any] = {}
        filtered = parse_blast(raw_path, args.min_alignment_bp, args.min_identity_pct, args.max_evalue,
                               limits=(args.max_target_seqs, args.max_hsps), audit=audit)
        require(set(audit) <= query_fasta_ids, "BLAST output contains an unknown query")
        filtered = expand_equivalent_query_rows(
            filtered, candidate_payload.get("sequence_equivalence_classes", []), query_fasta_ids
        )
        matches_path = output / f"matches.{task}.tsv"
        query_summaries = summarize_task(filtered, queries, windows, mappings, matches_path)
        for query_summary in query_summaries:
            query_summary.update(audit.get(representatives[query_summary["query_id"]], {
                "raw_target_count": 0, "max_hsps_observed": 0,
                "target_cap_reached": False, "hsp_cap_reached": False, "counts_are_lower_bounds": False,
            }))
        task_summaries[task] = {
            "filtered_hsp_count": len(filtered),
            "queries": query_summaries,
            "raw_hits_path": raw_path.name, "matches_path": matches_path.name,
            "raw_hits_sha256": sha256_file(raw_path), "matches_sha256": sha256_file(matches_path),
        }
    write_json(output / "comparison.json", {
        "schema": SCHEMA,
        "candidate_metadata_path": str(candidates_json),
        "candidate_metadata_sha256": candidate_digest,
        "query_fasta_path": str(query_fasta), "query_fasta_sha256": fasta_digest,
        "promoterome_path": str(promoterome),
        "promoterome_receipt_sha256": reference_digest,
        "background": {key: reference[key] for key in ("genome_id", "upstream_bp", "downstream_bp",
                                                        "unique_promoter_window_count", "sequence_orientation")},
        "producer_sha256": sha256_file(Path(__file__)),
        "thresholds": {"min_alignment_bp": args.min_alignment_bp,
                       "min_identity_pct": args.min_identity_pct, "max_evalue": args.max_evalue,
                       "max_target_seqs": args.max_target_seqs, "max_hsps_per_target": args.max_hsps},
        "search_policy": {"dust": "yes", "soft_masking": True},
        "counting_policy_id": COUNTING_POLICY,
        "counting_policy": (
            "Frequency is reported separately for distinct genomic promoter windows, genes, and transcripts. "
            "Overlapping windows and windows sharing any resolved query gene ID are excluded together "
            "from other_promoters, coverage tiers, frequency strips and ranked other-gene rows. "
            "Without a resolved query gene, only non-overlap is established; this is explicitly labelled."
        ),
        "completeness_policy": (
            "Target/HSP cap saturation is assessed on raw output before filters. Saturation makes counts "
            "lower bounds; no saturation means only that no output cap was observed, not exhaustive homology."
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
