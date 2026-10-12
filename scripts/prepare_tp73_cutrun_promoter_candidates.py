#!/usr/bin/env python3
"""Prepare TP73 CUT&RUN-supported TSS-window queries for promoterome comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

try:
    from .prepare_transcript_promoterome import validate_promoterome, window_id
    from .prepare_regulatory_region_indexes import require
except ImportError:
    from prepare_transcript_promoterome import validate_promoterome, window_id
    from prepare_regulatory_region_indexes import require


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


def mean_zero_filled(bw, chromosome: str, start: int, end: int) -> float:
    require(end > start >= 0, "invalid BigWig interval")
    values = bw.values(chromosome, start, end)
    require(len(values) == end - start, "BigWig did not return the complete window")
    require(all(math.isfinite(value) or math.isnan(value) for value in values),
            "BigWig signal contains infinite values")
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
                    require(current not in found, f"duplicate promoter FASTA record: {current}")
                    found[current] = []
            elif current is not None:
                found[current].append(line.strip())
    missing = wanted - set(found)
    if missing:
        raise RuntimeError(f"missing promoter FASTA records: {sorted(missing)}")
    return {name: "".join(chunks) for name, chunks in found.items()}


def selected_inputs(promoterome: Path):
    receipt = validate_promoterome(promoterome)
    require(receipt["genome_id"] == "Human GRCh38 Ensembl 116",
            "TP73 track selection requires Human GRCh38 Ensembl 116")
    require((receipt["upstream_bp"], receipt["downstream_bp"]) == (2000, 200),
            "TP73 selection requires the declared -2000/+200 promoterome")
    window_rows = load_tsv(promoterome / "promoter_windows.tsv")
    windows = {row["promoter_id"]: row for row in window_rows}
    mappings = load_tsv(promoterome / "promoter_transcripts.tsv")
    require(len(windows) == len(window_rows) == receipt["unique_promoter_window_count"],
            "promoter window inventory disagrees with receipt")
    require(len(mappings) == receipt["included_transcript_count"],
            "transcript inventory disagrees with receipt")
    expected = {transcript: gene for gene, transcripts in SELECTED_TRANSCRIPTS.items()
                for transcript in transcripts}
    selected = [row for row in mappings if row["transcript_id"] in expected]
    require(len(selected) == len(expected)
            and {row["transcript_id"] for row in selected} == set(expected),
            "selected transcripts are missing or duplicated in the promoterome")
    require(all(row["gene_name"] == expected[row["transcript_id"]] and row["gene_id"]
                and row["promoter_id"] in windows for row in selected),
            "selected transcript gene/window mapping disagrees with selection")
    promoter_ids = {row["promoter_id"] for row in selected}
    sequences = extract_fasta(promoterome / "promoter_windows.fa", promoter_ids)
    for promoter_id in sorted(promoter_ids):
        row = windows[promoter_id]
        start, end, tss = (int(row[key]) for key in ("start_0based", "end_0based_exclusive", "tss_1based"))
        strand = row["strand"]
        require(strand in {"+", "-"} and start >= 0 and end > start, "invalid selected window geometry")
        upstream, downstream = receipt["upstream_bp"], receipt["downstream_bp"]
        expected_span = ((tss - upstream - 1, tss + downstream) if strand == "+"
                         else (tss - downstream - 1, tss + upstream))
        require((start, end) == expected_span
                and str(row.get("boundary_clipped", "")).lower() == "false",
                "selected window is clipped or disagrees with its strand/TSS")
        require(promoter_id == window_id(row["chromosome"], start, end, strand),
                "selected promoter identity disagrees with geometry")
        require(len(sequences[promoter_id]) == end - start
                and set(sequences[promoter_id].upper()) <= set("ACGTRYSWKMBDHVN"),
                "selected promoter sequence disagrees with window length/alphabet")
    return receipt, windows, selected, sequences


def prepare(args, *, open_bigwig=None) -> None:
    output = args.output.resolve()
    require(not output.exists() or not any(output.iterdir()), "Output directory must be absent or empty")
    promoterome = args.promoterome.resolve()
    track_root = args.track_root.resolve()
    receipt, windows, selected_mappings, sequences = selected_inputs(promoterome)
    selected_ids = set().union(*SELECTED_TRANSCRIPTS.values())
    promoter_ids = sorted(sequences)
    repo = Path(__file__).resolve().parents[1]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    require(args.source_revision == revision, "--source-revision must be the current checkout's full HEAD SHA")
    source_bindings = {
        "promoterome_receipt_sha256": sha256_file(promoterome / "receipt.json"),
        "producer_sha256": sha256_file(Path(__file__)),
        "track_sha256": {f"{cell}:{condition}": sha256_file(track_root / filename)
                         for cell, conditions in TRACKS.items() for condition, filename in conditions.items()},
    }
    if open_bigwig is None:
        import pyBigWig
        open_bigwig = pyBigWig.open
    opened = {}
    try:
        for cell, conditions in TRACKS.items():
            opened[cell] = {}
            for condition, filename in conditions.items():
                opened[cell][condition] = open_bigwig(str(track_root / filename))
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
                    "non_claim": "Positive window-mean difference only; no peak, significance, direct-binding or activity claim",
                },
                "genome_extraction": {
                    "genome_id": receipt["genome_id"],
                    "chromosome": chromosome,
                    "start_1based": start + 1,
                    "end_1based": end,
                    "strand": window["strand"],
                    "tss_1based": int(window["tss_1based"]),
                    "promoter_upstream_bp": receipt["upstream_bp"],
                    "promoter_downstream_bp": receipt["downstream_bp"],
                    "gene_id": members[0]["gene_id"],
                    "gene_name": gene_names[0],
                    "transcript_ids": transcript_ids,
                    "anchor_verified": True,
                    "anchor_verification_basis": "receipt-bound window, transcript membership, TSS geometry and FASTA length",
                },
                "promoterome_id": promoter_id,
            })

        require(regions, "no selected TSS window passes the declared signal-difference rule")
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
            "source_revision": revision,
            "regions": regions,
            "sequence_equivalence_classes": equivalence_classes,
            "selection_policy": {
                "genes": sorted(SELECTED_TRANSCRIPTS),
                "transcripts": sorted(selected_ids),
                "window": {"upstream_bp": receipt["upstream_bp"], "downstream_bp": receipt["downstream_bp"]},
                "cutrun_support_rule": (
                    "At least one experimental TP73 mean exceeds its matched GFP-control mean "
                    "within the exact promoter window"
                ),
            },
            "source_bindings": source_bindings,
        }
        output.mkdir(parents=True, exist_ok=True)
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--promoterome", type=Path, required=True)
    parser.add_argument("--track-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
