#!/usr/bin/env python3
"""Prepare TP73 CUT&RUN-supported TSS-window queries for promoterome comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import pyBigWig


TRACKS = {
    "SAOS-2": {
        "GFP": "tp73_saos2_GFP_R1.bigWig",
        "TAp73alpha": "tp73_saos2_TA_R1.bigWig",
        "DNp73beta": "tp73_saos2_DN_R1.bigWig",
    },
    "SK-MEL-29-2": {
        "GFP": "tp73_skmel29_2_GFP_R1.bigWig",
        "TAp73alpha": "tp73_skmel29_2_TA_R1.bigWig",
        "DNp73beta": "tp73_skmel29_2_DN_R1.bigWig",
    },
}

SELECTED_TRANSCRIPTS = {
    "CD44": {"ENST00000263398", "ENST00000278386", "ENST00000428726"},
    "TGFB1": {"ENST00000221930", "ENST00001090430", "ENST00001114525"},
    "SERPINE1": {"ENST00000223095", "ENST00000870828", "ENST00000950058"},
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mean_zero_filled(bw: pyBigWig.pyBigWig, chromosome: str, start: int, end: int) -> float:
    values = bw.values(chromosome, start, end)
    return sum(0.0 if math.isnan(value) else value for value in values) / len(values)


def extract_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    found: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                current = name if name in wanted else None
                if current is not None:
                    found[current] = []
            elif current is not None:
                found[current].append(line.strip())
    missing = wanted - set(found)
    if missing:
        raise RuntimeError(f"missing promoter FASTA records: {sorted(missing)}")
    return {name: "".join(chunks) for name, chunks in found.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--promoterome", type=Path, required=True)
    parser.add_argument("--track-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    promoterome = args.promoterome.resolve()
    track_root = args.track_root.resolve()

    windows = {row["promoter_id"]: row for row in load_tsv(promoterome / "promoter_windows.tsv")}
    mappings = load_tsv(promoterome / "promoter_transcripts.tsv")
    selected_ids = set().union(*SELECTED_TRANSCRIPTS.values())
    selected_mappings = [row for row in mappings if row["transcript_id"] in selected_ids]
    promoter_ids = sorted({row["promoter_id"] for row in selected_mappings})
    sequences = extract_fasta(promoterome / "promoter_windows.fa", set(promoter_ids))

    opened = {
        cell: {condition: pyBigWig.open(str(track_root / filename))
               for condition, filename in conditions.items()}
        for cell, conditions in TRACKS.items()
    }
    try:
        regions = []
        for promoter_id in promoter_ids:
            window = windows[promoter_id]
            chromosome = window["chromosome"]
            start = int(window["start_0based"])
            end = int(window["end_0based_exclusive"])
            members = [row for row in selected_mappings if row["promoter_id"] == promoter_id]
            evidence = []
            supported = False
            for cell in sorted(opened):
                means = {condition: mean_zero_filled(bw, chromosome, start, end)
                         for condition, bw in opened[cell].items()}
                deltas = {condition: means[condition] - means["GFP"]
                          for condition in ("TAp73alpha", "DNp73beta")}
                if max(deltas.values()) > 0:
                    supported = True
                evidence.append({
                    "cell_line": cell,
                    "mean_signal_per_base": means,
                    "experimental_minus_matched_gfp": deltas,
                    "source_ids": {
                        condition: f"E-MTAB-15709:tp73:{cell}:{condition}:R1"
                        for condition in means
                    },
                })
            if not supported:
                continue
            sequence = sequences[promoter_id]
            gene_names = sorted({row["gene_name"] for row in members})
            transcript_ids = sorted({row["transcript_id"] for row in members})
            region_id = f"{gene_names[0]}_{promoter_id}_tp73_cutrun_supported"
            regions.append({
                "region_id": region_id,
                "input_kind": "transcript_promoter",
                "comparison_class": "tp73_cutrun_supported_tss_window",
                "sequence_length_bp": len(sequence),
                "sequence_sha256": f"sha256:{sha256_bytes(sequence.encode('ascii'))}",
                "sequence_orientation": "biological_5prime_to_3prime",
                "gene_query": gene_names[0],
                "transcript_ids": transcript_ids,
                "cutrun_support": {
                    "criterion": (
                        "At least one TP73 experimental BigWig mean across the exact "
                        "TSS window exceeds the matched GFP-control mean in the same cell line"
                    ),
                    "dataset": "E-MTAB-15709",
                    "factor": "TP73",
                    "evidence": evidence,
                    "non_claim": "Signal support is occupancy evidence, not proof of direct binding or promoter activity",
                },
                "genome_extraction": {
                    "genome_id": "Human GRCh38 Ensembl 116",
                    "chromosome": chromosome,
                    "start_1based": start + 1,
                    "end_1based": end,
                    "strand": window["strand"],
                    "tss_1based": int(window["tss_1based"]),
                    "promoter_upstream_bp": 2000,
                    "promoter_downstream_bp": 200,
                    "gene_id": members[0]["gene_id"],
                    "gene_name": gene_names[0],
                    "transcript_ids": transcript_ids,
                    "anchor_verified": True,
                },
                "promoterome_id": promoter_id,
            })

        sequence_to_members: dict[str, list[str]] = {}
        sequence_by_region: dict[str, str] = {}
        for region in regions:
            sequence = sequences[region["promoterome_id"]]
            sequence_by_region[region["region_id"]] = sequence
            sequence_to_members.setdefault(sequence, []).append(region["region_id"])
        equivalence_classes = []
        representatives = []
        for sequence, member_ids in sorted(sequence_to_members.items(), key=lambda item: sorted(item[1])):
            member_ids = sorted(member_ids)
            representative = member_ids[0]
            representatives.append((representative, sequence))
            equivalence_classes.append({
                "sequence_equivalence_class_id": f"region_seq_eq_{sha256_bytes(sequence.encode('ascii'))[:16]}",
                "representative_region_id": representative,
                "member_region_ids": member_ids,
                "sequence_length_bp": len(sequence),
                "sequence_sha256": f"sha256:{sha256_bytes(sequence.encode('ascii'))}",
            })

        payload = {
            "schema": "gentle.regulatory_region_comparison_sequences.v1",
            "dataset_id": "tp73_cutrun_supported_selected_gene_tss_windows_grch38_ensembl116_v1",
            "source_revision": args.source_revision,
            "regions": regions,
            "sequence_equivalence_classes": equivalence_classes,
            "selection_policy": {
                "genes": sorted(SELECTED_TRANSCRIPTS),
                "transcripts": sorted(selected_ids),
                "window": {"upstream_bp": 2000, "downstream_bp": 200},
                "cutrun_support_rule": (
                    "At least one experimental TP73 mean exceeds its matched GFP-control mean "
                    "within the exact promoter window"
                ),
            },
            "source_bindings": {
                "promoterome_receipt_sha256": sha256_file(promoterome / "receipt.json"),
                "track_sha256": {
                    f"{cell}:{condition}": sha256_file(track_root / filename)
                    for cell, conditions in TRACKS.items()
                    for condition, filename in conditions.items()
                },
            },
        }
        (output / "candidate_regions.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (output / "candidate_regions.fa").open("w", encoding="ascii") as handle:
            for representative, sequence in sorted(representatives):
                handle.write(f">{representative}\n")
                for offset in range(0, len(sequence), 80):
                    handle.write(sequence[offset:offset + 80] + "\n")
        print(json.dumps({
            "regions": len(regions),
            "sequence_equivalence_classes": len(equivalence_classes),
            "candidate_json": str(output / "candidate_regions.json"),
            "candidate_fasta": str(output / "candidate_regions.fa"),
        }, indent=2))
    finally:
        for conditions in opened.values():
            for bw in conditions.values():
                bw.close()


if __name__ == "__main__":
    main()
