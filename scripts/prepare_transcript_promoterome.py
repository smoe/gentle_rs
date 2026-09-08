#!/usr/bin/env python3
"""Prepare a transcript-linked promoterome background from a GENtle genome.

GENtle resolves the exact prepared genome and transcript index. BEDTools performs
strand-aware bulk extraction. Exact genomic windows are indexed once, while a
separate mapping retains every transcript and gene represented by each window.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

try:
    from .prepare_regulatory_region_indexes import (
        checked_id, require, run_command, sha256_file, write_json,
    )
except ImportError:
    from prepare_regulatory_region_indexes import (  # type: ignore[no-redef]
        checked_id, require, run_command, sha256_file, write_json,
    )


SCHEMA = "gentle.transcript_promoterome_preparation.v1"
RECEIPT_SCHEMA = "gentle.transcript_promoterome_preparation_receipt.v1"


def promoter_interval(transcript: dict[str, Any], upstream_bp: int,
                      downstream_bp: int, contig_length: int) -> tuple[int, int, int, bool]:
    strand = transcript.get("strand")
    require(strand in {"+", "-"}, "transcript strand must be '+' or '-'")
    start = transcript.get("transcript_start_1based")
    end = transcript.get("transcript_end_1based")
    require(isinstance(start, int) and isinstance(end, int) and 1 <= start <= end,
            "transcript coordinates must be positive 1-based inclusive")
    tss = end if strand == "-" else start
    if strand == "+":
        start_1based = max(1, tss - upstream_bp)
        end_1based = min(contig_length, tss + downstream_bp)
    else:
        start_1based = max(1, tss - downstream_bp)
        end_1based = min(contig_length, tss + upstream_bp)
    clipped = (start_1based != tss - (downstream_bp if strand == "-" else upstream_bp)
               or end_1based != tss + (upstream_bp if strand == "-" else downstream_bp))
    return start_1based - 1, end_1based, tss, clipped


def window_id(chromosome: str, start_0based: int, end_0based_exclusive: int,
              strand: str) -> str:
    identity = f"{chromosome}\t{start_0based}\t{end_0based_exclusive}\t{strand}".encode()
    return "promoter_" + hashlib.sha256(identity).hexdigest()[:20]


def build_windows(transcripts: list[dict[str, Any]], contig_lengths: dict[str, int],
                  upstream_bp: int, downstream_bp: int) -> tuple[list[dict[str, Any]], list[dict[str, str]], dict[str, int]]:
    windows: dict[tuple[str, int, int, str], dict[str, Any]] = {}
    mappings: list[dict[str, str]] = []
    excluded = {"invalid_record": 0, "contig_absent": 0}
    for transcript in transcripts:
        chromosome = transcript.get("chromosome")
        if not isinstance(chromosome, str) or chromosome not in contig_lengths:
            excluded["contig_absent"] += 1
            continue
        try:
            start_0based, end_exclusive, tss, clipped = promoter_interval(
                transcript, upstream_bp, downstream_bp, contig_lengths[chromosome]
            )
        except RuntimeError:
            excluded["invalid_record"] += 1
            continue
        strand = transcript["strand"]
        key = (chromosome, start_0based, end_exclusive, strand)
        row = windows.get(key)
        if row is None:
            row = {
                "promoter_id": window_id(*key),
                "chromosome": chromosome,
                "start_0based": start_0based,
                "end_0based_exclusive": end_exclusive,
                "strand": strand,
                "tss_1based": tss,
                "boundary_clipped": clipped,
                "gene_ids": set(),
                "transcript_count": 0,
            }
            windows[key] = row
        gene_id = str(transcript.get("gene_id") or "")
        gene_name = str(transcript.get("gene_name") or "")
        transcript_id = str(transcript.get("transcript_id") or "")
        require(transcript_id, "valid transcript record lacks transcript_id")
        row["gene_ids"].add(gene_id)
        row["transcript_count"] += 1
        mappings.append({
            "promoter_id": row["promoter_id"], "gene_id": gene_id,
            "gene_name": gene_name, "transcript_id": transcript_id,
        })
    ordered = sorted(windows.values(), key=lambda row: (
        row["chromosome"], row["start_0based"], row["end_0based_exclusive"], row["strand"]
    ))
    for row in ordered:
        row["gene_count"] = len(row.pop("gene_ids"))
    mappings.sort(key=lambda row: (row["promoter_id"], row["gene_id"], row["transcript_id"]))
    return ordered, mappings, excluded


def read_fai(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            require(len(fields) >= 2, f"invalid FAI row in {path}")
            lengths[fields[0]] = int(fields[1])
    return lengths


def normalize_bedtools_fasta(path: Path, expected_ids: set[str]) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    seen: set[str] = set()
    normalized: list[str] = []
    for line in lines:
        if line.startswith(">"):
            record_id = line[1:].split("::", 1)[0]
            if record_id.endswith("(+)") or record_id.endswith("(-)"):
                record_id = record_id[:-3]
            require(record_id in expected_ids, f"unexpected BEDTools FASTA id '{record_id}'")
            require(record_id not in seen, f"duplicate BEDTools FASTA id '{record_id}'")
            seen.add(record_id)
            normalized.append(">" + record_id)
        else:
            normalized.append(line.upper())
    require(seen == expected_ids, "BEDTools FASTA records do not match promoter windows")
    path.write_text("\n".join(normalized) + "\n", encoding="utf-8")


def prepare(args: argparse.Namespace) -> None:
    require(args.upstream_bp >= 0 and args.downstream_bp >= 0
            and args.upstream_bp + args.downstream_bp > 0, "invalid promoter span")
    dataset_id = checked_id(args.dataset_id, "dataset_id")
    repo_root = args.repo_root.resolve(strict=True)
    gentle = args.gentle.resolve(strict=True)
    catalog = args.catalog.resolve(strict=True)
    cache_dir = args.cache_dir.resolve(strict=True)
    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()), "Output directory must be absent or empty")
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    indexes = output / "indexes"
    logs.mkdir()
    indexes.mkdir()
    receipts: list[dict[str, Any]] = []

    version = run_command([str(gentle), "--version"], cwd=repo_root, timeout=30,
                          log_dir=logs, label="gentle_version", receipts=receipts)
    status_result = run_command([
        str(gentle), "genomes", "status", args.genome_id,
        "--catalog", str(catalog), "--cache-dir", str(cache_dir),
    ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
        label="genome_status", receipts=receipts)
    status = json.loads(status_result.stdout)
    require(status.get("component_ready") is True, "GENtle reports that the genome is not prepared")
    components = status.get("components", {})
    require(components.get("transcript_index_ready") is True, "transcript index is not ready")
    require(components.get("fasta_index_ready") is True, "FASTA index is not ready")
    transcript_path = Path(components["transcript_index_path"]).resolve(strict=True)
    sequence_path = Path(components["sequence_path"]).resolve(strict=True)
    fai_path = Path(components["fasta_index_path"]).resolve(strict=True)

    transcripts = json.loads(transcript_path.read_text(encoding="utf-8"))
    require(isinstance(transcripts, list), "GENtle transcript index must be an array")
    windows, mappings, excluded = build_windows(
        transcripts, read_fai(fai_path), args.upstream_bp, args.downstream_bp
    )
    require(windows, "no promoter windows were resolved")

    bed_path = output / "promoter_windows.bed"
    with bed_path.open("w", encoding="utf-8", newline="") as handle:
        for row in windows:
            handle.write(
                f"{row['chromosome']}\t{row['start_0based']}\t{row['end_0based_exclusive']}\t"
                f"{row['promoter_id']}\t0\t{row['strand']}\n"
            )
    windows_path = output / "promoter_windows.tsv"
    with windows_path.open("w", encoding="utf-8", newline="") as handle:
        columns = ["promoter_id", "chromosome", "start_0based", "end_0based_exclusive",
                   "strand", "tss_1based", "boundary_clipped", "gene_count", "transcript_count"]
        handle.write("\t".join(columns) + "\n")
        for row in windows:
            handle.write("\t".join(str(row[column]).lower() if isinstance(row[column], bool)
                                   else str(row[column]) for column in columns) + "\n")
    mappings_path = output / "promoter_transcripts.tsv"
    with mappings_path.open("w", encoding="utf-8", newline="") as handle:
        columns = ["promoter_id", "gene_id", "gene_name", "transcript_id"]
        handle.write("\t".join(columns) + "\n")
        for row in mappings:
            handle.write("\t".join(row[column].replace("\t", " ") for column in columns) + "\n")

    fasta_path = output / "promoter_windows.fa"
    run_command([
        args.bedtools, "getfasta", "-fi", str(sequence_path), "-bed", str(bed_path),
        "-fo", str(fasta_path), "-nameOnly", "-s",
    ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
        label="bedtools_getfasta", receipts=receipts)
    normalize_bedtools_fasta(fasta_path, {row["promoter_id"] for row in windows})

    blast_prefix = indexes / "promoterome"
    run_command([args.makeblastdb, "-version"], cwd=repo_root, timeout=30,
                log_dir=logs, label="makeblastdb_version", receipts=receipts)
    run_command([
        args.makeblastdb, "-in", str(fasta_path), "-dbtype", "nucl", "-parse_seqids",
        "-blastdb_version", "5", "-title", dataset_id, "-out", str(blast_prefix),
    ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
        label="makeblastdb", receipts=receipts)
    if not args.skip_minimap2:
        run_command([args.minimap2, "--version"], cwd=repo_root, timeout=30,
                    log_dir=logs, label="minimap2_version", receipts=receipts)
        run_command([args.minimap2, "-d", str(indexes / "promoterome.mmi"), str(fasta_path)],
                    cwd=repo_root, timeout=args.timeout, log_dir=logs,
                    label="minimap2_index", receipts=receipts)

    artifact_hashes = {
        path.relative_to(output).as_posix(): sha256_file(path)
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name != "receipt.json"
    }
    write_json(output / "receipt.json", {
        "schema": RECEIPT_SCHEMA, "bundle_schema": SCHEMA, "dataset_id": dataset_id,
        "genome_id": args.genome_id, "upstream_bp": args.upstream_bp,
        "downstream_bp": args.downstream_bp,
        "coordinate_policy": "GENtle extract-promoter-compatible 1-based inclusive TSS geometry; BED is 0-based half-open",
        "sequence_orientation": "transcript_5prime_to_3prime_via_bedtools_strand",
        "input_transcript_count": len(transcripts), "included_transcript_count": len(mappings),
        "unique_promoter_window_count": len(windows), "excluded": excluded,
        "catalog_path": str(catalog), "catalog_sha256": sha256_file(catalog),
        "transcript_index_path": str(transcript_path), "transcript_index_sha256": sha256_file(transcript_path),
        "sequence_path": str(sequence_path), "sequence_fai_sha256": sha256_file(fai_path),
        "gentle_binary": str(gentle), "gentle_binary_sha256": sha256_file(gentle),
        "gentle_version_stdout": version.stdout.decode(errors="replace").strip(),
        "commands": receipts, "artifacts": artifact_hashes,
        "counting_policy": (
            "One indexed record per exact genomic window/strand. Gene and transcript frequencies "
            "must be derived through promoter_transcripts.tsv, not inferred from FASTA record count."
        ),
        "non_claims": [
            "A transcript-TSS window is a search background, not a claim of promoter activity.",
            "Sequence recurrence does not establish functional equivalence or reporter sufficiency.",
        ],
    })
    print(json.dumps({"status": "ok", "dataset_id": dataset_id,
                      "transcripts": len(mappings), "promoter_windows": len(windows),
                      "output": str(output)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--upstream-bp", type=int, default=2000)
    parser.add_argument("--downstream-bp", type=int, default=200)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--gentle", type=Path, default=Path(__file__).resolve().parents[1] / "target/debug/gentle_cli")
    parser.add_argument("--catalog", type=Path, default=Path(__file__).resolve().parents[1] / "assets/genomes.json")
    parser.add_argument("--cache-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data/genomes")
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--skip-minimap2", action="store_true")
    parser.add_argument("--timeout", type=int, default=7200)
    args = parser.parse_args()
    require(args.timeout > 0, "--timeout must be positive")
    prepare(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
