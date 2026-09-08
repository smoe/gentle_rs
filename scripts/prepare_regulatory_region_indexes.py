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
GENOMIC_REGION_SET_SCHEMA = "gentle.genomic_region_set.v1"
GENOMIC_REGION_SCHEMA = "gentle.genomic_region_of_interest.v1"
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
    extraction = dict(rows[-1])
    # Wall-clock recording remains in the retained GENtle project and command
    # receipt. It is not part of the scientific coordinate/source binding and
    # would otherwise make candidate_regions.json differ across identical runs.
    extraction.pop("recorded_at_unix_ms", None)
    return extraction


def resolve_repo_path(repo_root: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (repo_root / path).resolve()


def matches_optional_filter(value: Any, allowed: set[str]) -> bool:
    return not allowed or value in allowed


def canonical_region_tasks(manifest: dict[str, Any], manifest_path: Path,
                           genome: dict[str, Any]) -> list[dict[str, Any]]:
    """Expand canonical region-set sources into extraction tasks.

    Region-set paths are resolved relative to the manifest so a comparison
    bundle can move as a unit. Assembly/taxon checks happen before GENtle is
    invoked; a human-readable genome_id alone is not treated as an assembly
    proof.
    """
    sources = manifest.get("region_sets", [])
    require(isinstance(sources, list), "manifest.region_sets must be an array")
    if not sources:
        return []

    expected = genome.get("expected_reference")
    require(isinstance(expected, dict),
            "genome.expected_reference is required when region_sets are used")
    assembly_name_values = expected.get("assembly_names", [])
    assembly_accession_values = expected.get("assembly_accessions", [])
    require(isinstance(assembly_name_values, list)
            and all(isinstance(value, str) for value in assembly_name_values),
            "genome.expected_reference.assembly_names must be an array of strings")
    require(isinstance(assembly_accession_values, list)
            and all(isinstance(value, str) for value in assembly_accession_values),
            "genome.expected_reference.assembly_accessions must be an array of strings")
    assembly_names = set(assembly_name_values)
    assembly_accessions = set(assembly_accession_values)
    expected_taxon = expected.get("taxon_id")
    require(assembly_names or assembly_accessions,
            "genome.expected_reference must declare assembly_names or assembly_accessions")

    tasks: list[dict[str, Any]] = []
    for source_index, source in enumerate(sources):
        require(isinstance(source, dict), f"region_sets[{source_index}] must be an object")
        comparison_class = checked_id(
            source.get("comparison_class"),
            f"region_sets[{source_index}].comparison_class",
        )
        path_value = source.get("path")
        require(isinstance(path_value, str) and path_value.strip(),
                f"region_sets[{source_index}].path is required")
        source_path = Path(path_value)
        if not source_path.is_absolute():
            source_path = manifest_path.parent / source_path
        source_path = source_path.resolve(strict=True)
        payload = json.loads(source_path.read_text(encoding="utf-8"))
        require(payload.get("schema") == GENOMIC_REGION_SET_SCHEMA,
                f"{source_path}: expected schema '{GENOMIC_REGION_SET_SCHEMA}'")
        regions = payload.get("regions")
        require(isinstance(regions, list), f"{source_path}: regions must be an array")

        prefix = source.get("id_prefix", "")
        require(isinstance(prefix, str), f"region_sets[{source_index}].id_prefix must be a string")
        if prefix:
            checked_id(prefix.rstrip("_.-"), f"region_sets[{source_index}].id_prefix")
        filter_values: list[list[Any]] = []
        for filter_name in ("region_ids", "purposes", "selection_methods"):
            values = source.get(filter_name, [])
            require(isinstance(values, list),
                    f"region_sets[{source_index}].{filter_name} must be an array")
            require(all(isinstance(value, str) for value in values),
                    f"region_sets[{source_index}].{filter_name} must contain strings")
            filter_values.append(values)
        region_ids, purposes, methods = (set(values) for values in filter_values)
        annotation_scope = source.get("annotation_scope", "none")
        require(annotation_scope in {"none", "core", "full"},
                f"region_sets[{source_index}]: invalid annotation_scope")

        selected = 0
        for region_index, region in enumerate(regions):
            require(isinstance(region, dict), f"{source_path}: regions[{region_index}] must be an object")
            require(region.get("schema") == GENOMIC_REGION_SCHEMA,
                    f"{source_path}: region must use schema '{GENOMIC_REGION_SCHEMA}'")
            source_region_id = checked_id(region.get("region_id"),
                                          f"{source_path}: regions[{region_index}].region_id")
            if not matches_optional_filter(source_region_id, region_ids):
                continue
            if not matches_optional_filter(region.get("purpose"), purposes):
                continue
            if not matches_optional_filter(region.get("selection_method"), methods):
                continue

            interval = region.get("interval")
            require(isinstance(interval, dict), f"{source_region_id}: interval is required")
            require(interval.get("coordinate_convention") == "zero_based_half_open",
                    f"{source_region_id}: expected zero_based_half_open coordinates")
            reference = interval.get("reference")
            require(isinstance(reference, dict), f"{source_region_id}: interval.reference is required")
            assembly_name = reference.get("assembly_name")
            assembly_accession = reference.get("assembly_accession")
            name_ok = bool(assembly_names and assembly_name in assembly_names)
            accession_ok = bool(assembly_accessions and assembly_accession in assembly_accessions)
            require(name_ok or accession_ok,
                    f"{source_region_id}: reference assembly is not allowed by genome.expected_reference")
            if expected_taxon is not None:
                require(reference.get("taxon_id") == expected_taxon,
                        f"{source_region_id}: taxon_id does not match genome.expected_reference")
            contig = reference.get("contig_name")
            require(isinstance(contig, str) and contig.strip(), f"{source_region_id}: contig_name is required")
            start_0based = interval.get("start_0based")
            end_exclusive = interval.get("end_0based_exclusive")
            require(isinstance(start_0based, int) and not isinstance(start_0based, bool) and start_0based >= 0,
                    f"{source_region_id}: invalid start_0based")
            require(isinstance(end_exclusive, int) and not isinstance(end_exclusive, bool)
                    and end_exclusive > start_0based,
                    f"{source_region_id}: invalid end_0based_exclusive")
            strand = interval.get("strand", "unstranded")
            require(strand in {"plus", "minus", "unstranded"},
                    f"{source_region_id}: invalid strand")

            region_id = checked_id(prefix + source_region_id,
                                   f"region_sets[{source_index}] output region_id")
            tasks.append({
                "kind": "canonical_region",
                "region_id": region_id,
                "comparison_class": comparison_class,
                "annotation_scope": annotation_scope,
                "contig_name": contig,
                # GENtle's CLI boundary is 1-based inclusive. For a canonical
                # [start,end) interval this conversion intentionally leaves
                # the numeric end unchanged.
                "start_1based": start_0based + 1,
                "end_1based": end_exclusive,
                "source_region": region,
                "source_region_set": {
                    "path": str(source_path),
                    "file_sha256": sha256_file(source_path),
                    "source_index": source_index,
                    "set_id": payload.get("set_id"),
                    "declared_content_sha256": payload.get("content_sha256"),
                    "source_region_id": source_region_id,
                },
            })
            selected += 1
        require(selected > 0, f"region_sets[{source_index}] filters selected no regions")
    return tasks


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

    region_specs = manifest.get("regions", [])
    require(isinstance(region_specs, list), "manifest.regions must be an array")
    tasks: list[dict[str, Any]] = []
    for index, spec in enumerate(region_specs):
        require(isinstance(spec, dict), f"regions[{index}] must be an object")
        task = dict(spec)
        task["kind"] = "transcript_promoter"
        task["region_id"] = checked_id(spec.get("region_id"), f"regions[{index}].region_id")
        task["comparison_class"] = checked_id(
            spec.get("comparison_class", "transcript_tss_window"),
            f"regions[{index}].comparison_class",
        )
        tasks.append(task)
    tasks.extend(canonical_region_tasks(manifest, manifest_path, genome))
    require(tasks, "manifest must declare at least one regions[] or region_sets[] input")
    declared_region_ids = [task["region_id"] for task in tasks]
    require(len(declared_region_ids) == len(set(declared_region_ids)),
            "region_id values must be unique after region-set prefixes are applied")
    records: list[tuple[str, str]] = []
    region_rows: list[dict[str, Any]] = []
    imported_sources: set[str] = set()
    for task in tasks:
        if task["kind"] != "canonical_region":
            continue
        source = task["source_region_set"]
        source_path = source["path"]
        if source_path in imported_sources:
            continue
        imported_sources.add(source_path)
        source_index = source["source_index"]
        import_request_path = output / f"region_set_{source_index:03d}.import_request.json"
        imported_set_id = checked_id(
            f"{dataset_id}.input{source_index:03d}", "canonical region import set_id"
        )
        write_json(import_request_path, {
            "path": source_path,
            "format": "json",
            "set_id_override": imported_set_id,
            "collision_policy": "reject",
            "max_bytes": 10_485_760,
            "max_rows": 100_000,
        })
        run_command([
            str(gentle), "--state", str(state_path), "shell",
            f"regions import @{import_request_path}",
        ], cwd=repo_root, timeout=args.timeout, log_dir=logs,
            label=f"region_set_{source_index:03d}_validate_import", receipts=command_receipts)
        for matching_task in tasks:
            if (matching_task["kind"] == "canonical_region"
                    and matching_task["source_region_set"]["path"] == source_path):
                matching_task["source_region_set"]["validated_import_set_id"] = imported_set_id

    for index, task in enumerate(tasks, start=1):
        region_id = task["region_id"]
        annotation_scope = task.get("annotation_scope", "core")
        require(annotation_scope in {"none", "core", "full"}, f"{region_id}: invalid annotation_scope")
        if task["kind"] == "transcript_promoter":
            gene_query = task.get("gene_query")
            require(isinstance(gene_query, str) and gene_query.strip(),
                    f"{region_id}: gene_query is required")
            upstream = int(task.get("upstream_bp", 5000))
            downstream = int(task.get("downstream_bp", 1000))
            require(upstream >= 0 and downstream >= 0 and upstream + downstream > 0,
                    f"{region_id}: invalid upstream/downstream span")
            command = [
                str(gentle), "--state", str(state_path), "genomes", "extract-promoter",
                genome_id, gene_query, "--output-id", region_id,
                "--upstream-bp", str(upstream), "--downstream-bp", str(downstream),
                "--annotation-scope", annotation_scope,
                "--catalog", str(catalog), "--cache-dir", str(cache_dir),
            ]
            if task.get("transcript_id"):
                command.extend(["--transcript-id", str(task["transcript_id"])])
            if task.get("occurrence") is not None:
                command.extend(["--occurrence", str(int(task["occurrence"]))])
        else:
            command = [
                str(gentle), "--state", str(state_path), "genomes", "extract-region",
                genome_id, task["contig_name"], str(task["start_1based"]),
                str(task["end_1based"]), "--output-id", region_id,
                "--annotation-scope", annotation_scope,
                "--catalog", str(catalog), "--cache-dir", str(cache_dir),
            ]
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
        row = {
            "region_id": region_id,
            "input_kind": task["kind"],
            "comparison_class": task["comparison_class"],
            "sequence_length_bp": len(sequence),
            "sequence_sha256": sha256_bytes(sequence.encode()),
            "genome_extraction": extraction,
        }
        if task["kind"] == "transcript_promoter":
            row.update({
                "gene_query": gene_query,
                "transcript_id": task.get("transcript_id"),
                "upstream_bp": upstream,
                "downstream_bp": downstream,
                "sequence_orientation": "biological_5prime_to_3prime",
            })
        else:
            row.update({
                "sequence_orientation": "assembly_reference_forward",
                "source_region_set": task["source_region_set"],
                "source_region": task["source_region"],
            })
        region_rows.append(row)

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
            "Coordinates and strand are copied from GENtle genome-extraction provenance and canonical "
            "region-set inputs. Transcript/TSS windows are exported in biological 5-prime-to-3-prime "
            "orientation; explicit genomic intervals are exported in assembly-reference-forward "
            "orientation. BLAST searches both strands and k-mers are strand-neutral."
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
