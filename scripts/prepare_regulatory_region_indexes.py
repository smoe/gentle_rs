#!/usr/bin/env python3
"""Prepare reproducible upstream-region comparison resources through GENtle.

The script delegates sequence extraction to GENtle, exports each extracted
sequence through the engine SaveFile operation, and then prepares two distinct
comparison resources:

* a BLAST database and all-vs-all reports for local nucleotide homology; and
* strand-neutral k-mer signatures for deliberately lower-stringency candidate
  discovery.

The k-mer output is not a uniqueness, orthology, binding, or regulatory-
function claim. Whole-genome uniqueness remains a separate GENtle homology
screen against explicitly prepared genomic-DNA indexes.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
from typing import Any


SCHEMA = "gentle.regulatory_region_index_preparation.v1"
RECEIPT_SCHEMA = "gentle.regulatory_region_index_preparation_receipt.v1"
REGION_SET_SCHEMA = "gentle.regulatory_region_comparison_sequences.v1"
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
DNA_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return "sha256:" + hashlib.file_digest(handle, "sha256").hexdigest()


def canonical_json_bytes(payload: Any) -> bytes:
    return (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def checked_id(value: Any, label: str) -> str:
    require(isinstance(value, str) and SAFE_ID.fullmatch(value) is not None,
            f"{label} must match {SAFE_ID.pattern}")
    return value


def run_command(command: list[str], *, cwd: Path, timeout: int,
                log_dir: Path, label: str, receipts: list[dict[str, Any]]) -> subprocess.CompletedProcess[bytes]:
    stdout_path = log_dir / f"{label}.stdout"
    stderr_path = log_dir / f"{label}.stderr"
    started = time.time()
    completed = subprocess.run(command, cwd=cwd, capture_output=True, timeout=timeout)
    stdout_path.write_bytes(completed.stdout)
    stderr_path.write_bytes(completed.stderr)
    receipts.append({
        "label": label,
        "command": command,
        "exit_code": completed.returncode,
        "elapsed_seconds": round(time.time() - started, 6),
        "stdout_path": stdout_path.relative_to(log_dir.parent).as_posix(),
        "stdout_sha256": sha256_file(stdout_path),
        "stdout_bytes": stdout_path.stat().st_size,
        "stderr_path": stderr_path.relative_to(log_dir.parent).as_posix(),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
    })
    require(completed.returncode == 0,
            f"Command '{label}' failed; inspect {stdout_path} and {stderr_path}")
    return completed


def parse_single_fasta(path: Path) -> tuple[str, str]:
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                require(header is None, f"Expected one FASTA record in {path}")
                header = line[1:]
            else:
                require(header is not None, f"Sequence appeared before FASTA header in {path}")
                sequence.append(line)
    require(header is not None and sequence, f"Empty FASTA record in {path}")
    joined = "".join(sequence).upper()
    require(set(joined) <= set("ACGTNRYKMSWBDHV"), f"Unsupported nucleotide symbols in {path}")
    return header, joined


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def canonical_kmers(sequence: str, k: int) -> set[str]:
    kmers: set[str] = set()
    for offset in range(0, max(0, len(sequence) - k + 1)):
        word = sequence[offset:offset + k]
        if set(word) <= set("ACGT"):
            kmers.add(min(word, reverse_complement(word)))
    return kmers


def write_kmer_resources(records: list[tuple[str, str]], k_values: list[int], output: Path) -> dict[str, str]:
    signatures: dict[str, dict[str, list[str]]] = {}
    sets: dict[tuple[str, int], set[str]] = {}
    for record_id, sequence in records:
        signatures[record_id] = {}
        for k in k_values:
            words = canonical_kmers(sequence, k)
            sets[(record_id, k)] = words
            signatures[record_id][str(k)] = sorted(words)

    signature_path = output / "kmer_signatures.json.gz"
    payload = {
        "schema": "gentle.regulatory_region_kmer_signatures.v1",
        "strand_policy": "canonical_forward_or_reverse_complement",
        "k_values": k_values,
        "interpretation": (
            "Candidate-generation evidence only. Shared short words can reveal weak local resemblance, "
            "but do not establish alignment, genomic uniqueness, orthology, TF binding, or regulation."
        ),
        "signatures": signatures,
    }
    with signature_path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as handle:
            handle.write(canonical_json_bytes(payload))

    table_path = output / "kmer_pairwise.tsv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("left_id\tright_id\tk\tleft_unique\tright_unique\tshared\tjaccard\tleft_containment\tright_containment\n")
        for left_index, (left_id, _) in enumerate(records):
            for right_id, _ in records[left_index + 1:]:
                for k in k_values:
                    left = sets[(left_id, k)]
                    right = sets[(right_id, k)]
                    shared = len(left & right)
                    union = len(left | right)
                    handle.write(
                        f"{left_id}\t{right_id}\t{k}\t{len(left)}\t{len(right)}\t{shared}\t"
                        f"{shared / union if union else 0:.8f}\t"
                        f"{shared / len(left) if left else 0:.8f}\t"
                        f"{shared / len(right) if right else 0:.8f}\n"
                    )
    return {
        signature_path.name: sha256_file(signature_path),
        table_path.name: sha256_file(table_path),
    }


def collapse_equivalent_records(records: list[tuple[str, str]]) -> tuple[list[tuple[str, str]], list[dict[str, Any]], dict[str, str]]:
    region_ids = [region_id for region_id, _ in records]
    require(len(region_ids) == len(set(region_ids)), "region_id values must be unique")
    equivalence_by_digest: dict[str, list[str]] = {}
    sequence_by_region = dict(records)
    for region_id, sequence in records:
        equivalence_by_digest.setdefault(sha256_bytes(sequence.encode()), []).append(region_id)
    equivalence_classes: list[dict[str, Any]] = []
    index_records: list[tuple[str, str]] = []
    class_by_region: dict[str, str] = {}
    for sequence_digest, member_ids in equivalence_by_digest.items():
        class_id = "region_seq_eq_" + sequence_digest.removeprefix("sha256:")[:16]
        representative_id = member_ids[0]
        equivalence_classes.append({
            "class_id": class_id,
            "sequence_sha256": sequence_digest,
            "representative_region_id": representative_id,
            "member_region_ids": member_ids,
        })
        index_records.append((representative_id, sequence_by_region[representative_id]))
        for member_id in member_ids:
            class_by_region[member_id] = class_id
    return index_records, equivalence_classes, class_by_region


def latest_genome_extraction(state: dict[str, Any], seq_id: str) -> dict[str, Any]:
    extractions = state.get("metadata", {}).get("provenance", {}).get("genome_extractions", [])
    rows = [row for row in extractions if row.get("seq_id") == seq_id]
    require(rows, f"GENtle state has no genome-extraction provenance for '{seq_id}'")
    return rows[-1]


def resolve_repo_path(repo_root: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (repo_root / path).resolve()


def prepare(args: argparse.Namespace) -> None:
    manifest_path = args.manifest.resolve(strict=True)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(manifest.get("schema") == SCHEMA, f"Expected schema '{SCHEMA}'")
    dataset_id = checked_id(manifest.get("dataset_id"), "dataset_id")
    genome = manifest.get("genome")
    require(isinstance(genome, dict), "manifest.genome must be an object")
    genome_id = genome.get("genome_id")
    require(isinstance(genome_id, str) and genome_id.strip(), "genome.genome_id is required")
    repo_root = args.repo_root.resolve(strict=True)
    catalog = (args.catalog.resolve(strict=True) if args.catalog else
               resolve_repo_path(repo_root, genome.get("catalog_path", "assets/genomes.json")))
    cache_dir = (args.cache_dir.resolve() if args.cache_dir else
                 resolve_repo_path(repo_root, genome.get("cache_dir", "data/genomes")))
    gentle = args.gentle.resolve(strict=True)
    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()), "Output directory must be absent or empty")
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    sequences_dir = output / "sequences"
    indexes_dir = output / "indexes"
    logs.mkdir()
    sequences_dir.mkdir()
    indexes_dir.mkdir()
    state_path = output / "regions.project.gentle.json"
    command_receipts: list[dict[str, Any]] = []

    version = run_command([str(gentle), "--version"], cwd=repo_root, timeout=30,
                          log_dir=logs, label="gentle_version", receipts=command_receipts)
    if args.prepare_genome:
        run_command([
            str(gentle), "genomes", "prepare", genome_id,
            "--catalog", str(catalog), "--cache-dir", str(cache_dir),
            "--timeout-secs", str(args.timeout),
        ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
            label="prepare_genome", receipts=command_receipts)

    region_specs = manifest.get("regions")
    require(isinstance(region_specs, list) and region_specs, "manifest.regions must be a non-empty array")
    declared_region_ids = [
        checked_id(spec.get("region_id"), f"regions[{index}].region_id")
        if isinstance(spec, dict) else ""
        for index, spec in enumerate(region_specs)
    ]
    require(all(declared_region_ids), "every regions[] entry must be an object")
    require(len(declared_region_ids) == len(set(declared_region_ids)), "region_id values must be unique")
    records: list[tuple[str, str]] = []
    region_rows: list[dict[str, Any]] = []
    for index, spec in enumerate(region_specs, start=1):
        require(isinstance(spec, dict), f"regions[{index - 1}] must be an object")
        region_id = checked_id(spec.get("region_id"), f"regions[{index - 1}].region_id")
        gene_query = spec.get("gene_query")
        require(isinstance(gene_query, str) and gene_query.strip(), f"{region_id}: gene_query is required")
        upstream = int(spec.get("upstream_bp", 5000))
        downstream = int(spec.get("downstream_bp", 1000))
        require(upstream >= 0 and downstream >= 0 and upstream + downstream > 0,
                f"{region_id}: invalid upstream/downstream span")
        annotation_scope = spec.get("annotation_scope", "core")
        require(annotation_scope in {"none", "core", "full"}, f"{region_id}: invalid annotation_scope")
        command = [
            str(gentle), "--state", str(state_path), "genomes", "extract-promoter",
            genome_id, gene_query, "--output-id", region_id,
            "--upstream-bp", str(upstream), "--downstream-bp", str(downstream),
            "--annotation-scope", annotation_scope,
            "--catalog", str(catalog), "--cache-dir", str(cache_dir),
        ]
        if spec.get("transcript_id"):
            command.extend(["--transcript-id", str(spec["transcript_id"])])
        if spec.get("occurrence") is not None:
            command.extend(["--occurrence", str(int(spec["occurrence"]))])
        run_command(command, cwd=repo_root, timeout=args.timeout, log_dir=logs,
                    label=f"{index:03d}_{region_id}_extract", receipts=command_receipts)

        fasta_path = sequences_dir / f"{region_id}.fa"
        operation_path = output / f"{region_id}.save_fasta.operation.json"
        operation = {"SaveFile": {"seq_id": region_id, "path": str(fasta_path), "format": "Fasta"}}
        write_json(operation_path, operation)
        run_command([str(gentle), "--state", str(state_path), "op", str(operation_path)],
                    cwd=repo_root, timeout=args.timeout, log_dir=logs,
                    label=f"{index:03d}_{region_id}_export", receipts=command_receipts)
        _, sequence = parse_single_fasta(fasta_path)
        records.append((region_id, sequence))

        state = json.loads(state_path.read_text(encoding="utf-8"))
        extraction = latest_genome_extraction(state, region_id)
        region_rows.append({
            "region_id": region_id,
            "gene_query": gene_query,
            "transcript_id": spec.get("transcript_id"),
            "upstream_bp": upstream,
            "downstream_bp": downstream,
            "sequence_length_bp": len(sequence),
            "sequence_sha256": sha256_bytes(sequence.encode()),
            "genome_extraction": extraction,
        })

    index_records, equivalence_classes, class_by_region = collapse_equivalent_records(records)
    for row in region_rows:
        row["sequence_equivalence_class_id"] = class_by_region[row["region_id"]]

    combined_fasta = output / "candidate_regions.fa"
    with combined_fasta.open("w", encoding="utf-8", newline="") as handle:
        for region_id, sequence in index_records:
            handle.write(f">{region_id}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")

    region_report = {
        "schema": REGION_SET_SCHEMA,
        "dataset_id": dataset_id,
        "genome_id": genome_id,
        "catalog_path": str(catalog),
        "catalog_sha256": sha256_file(catalog),
        "cache_dir": str(cache_dir),
        "coordinate_interpretation": (
            "Coordinates and strand are copied from GENtle genome-extraction provenance. "
            "Each FASTA sequence is exported in GENtle's biological 5-prime-to-3-prime orientation."
        ),
        "input_region_count": len(records),
        "unique_sequence_count": len(index_records),
        "sequence_equivalence_classes": equivalence_classes,
        "regions": region_rows,
    }
    write_json(output / "candidate_regions.json", region_report)

    k_values = manifest.get("indexing", {}).get("canonical_kmer_lengths", [7, 11, 15, 21])
    require(isinstance(k_values, list) and k_values, "canonical_kmer_lengths must be non-empty")
    k_values = sorted({int(value) for value in k_values})
    require(all(3 <= value <= 63 for value in k_values), "canonical k-mer lengths must be 3..63")
    write_kmer_resources(index_records, k_values, output)

    blast_enabled = not args.skip_blast and manifest.get("indexing", {}).get("blast", True)
    if blast_enabled:
        run_command([args.makeblastdb, "-version"], cwd=repo_root, timeout=30,
                    log_dir=logs, label="makeblastdb_version", receipts=command_receipts)
        run_command([args.blastn, "-version"], cwd=repo_root, timeout=30,
                    log_dir=logs, label="blastn_version", receipts=command_receipts)
        blast_prefix = indexes_dir / "candidate_regions"
        run_command([
            args.makeblastdb, "-in", str(combined_fasta), "-dbtype", "nucl",
            "-parse_seqids", "-blastdb_version", "5", "-title", dataset_id,
            "-out", str(blast_prefix),
        ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
            label="makeblastdb", receipts=command_receipts)
        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        for task in manifest.get("indexing", {}).get("blast_tasks", ["megablast", "blastn"]):
            require(task in {"megablast", "blastn", "dc-megablast"}, f"Unsupported BLAST task '{task}'")
            result_path = output / f"all_vs_all.{task}.tsv"
            run_command([
                args.blastn, "-task", task, "-query", str(combined_fasta),
                "-db", str(blast_prefix), "-dust", "yes", "-soft_masking", "true",
                "-outfmt", outfmt, "-out", str(result_path),
            ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
                label=f"all_vs_all_{task}", receipts=command_receipts)

    manifest_sha = sha256_file(manifest_path)
    artifact_hashes = {
        path.relative_to(output).as_posix(): sha256_file(path)
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name not in {"receipt.json", "checksums.sha256"}
    }
    receipt = {
        "schema": RECEIPT_SCHEMA,
        "dataset_id": dataset_id,
        "manifest_path": str(manifest_path),
        "manifest_sha256": manifest_sha,
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "producer": "Preparation/orchestration script; scientific sequence and coordinate authority remains GENtle.",
        "gentle_binary": str(gentle),
        "gentle_binary_sha256": sha256_file(gentle),
        "gentle_version_stdout": version.stdout.decode(errors="replace").strip(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "prepare_genome_authorized": bool(args.prepare_genome),
        "commands": command_receipts,
        "artifacts": artifact_hashes,
        "non_claims": [
            "The candidate-region BLAST database is not a whole genome and cannot establish genomic uniqueness.",
            "Short-word similarity is candidate-generation evidence, not an alignment, orthology, TF-binding, or regulatory-function claim.",
            "Cross-species orthology and same-genome repetition require separate declared GENtle homology targets and validated genomic indexes.",
        ],
    }
    write_json(output / "receipt.json", receipt)
    with (output / "checksums.sha256").open("w", encoding="utf-8", newline="") as handle:
        for relative, digest in sorted({**artifact_hashes, "receipt.json": sha256_file(output / "receipt.json")}.items()):
            handle.write(f"{digest.removeprefix('sha256:')}  {relative}\n")
    print(json.dumps({
        "status": "ok", "dataset_id": dataset_id, "regions": len(records),
        "unique_sequences": len(index_records),
        "output": str(output), "receipt": str(output / "receipt.json"),
    }, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--gentle", type=Path, default=Path(__file__).resolve().parents[1] / "target/debug/gentle_cli")
    parser.add_argument("--catalog", type=Path,
                        help="Override the manifest catalog path without changing its scientific genome_id")
    parser.add_argument("--cache-dir", type=Path,
                        help="Override the host-local prepared-cache root")
    parser.add_argument("--prepare-genome", action="store_true",
                        help="Explicitly authorize GENtle to prepare/download the declared genome before extraction")
    parser.add_argument("--skip-blast", action="store_true",
                        help="Prepare FASTA and k-mer evidence without running BLAST tools")
    parser.add_argument("--makeblastdb", default="makeblastdb")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    require(args.timeout > 0, "--timeout must be positive")
    prepare(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
