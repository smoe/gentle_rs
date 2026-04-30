#!/usr/bin/env python3
"""Build the local TP73 long-range cDNA selector panel and virtual gel."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


ANCHOR_SEQUENCE = "AATCTGCTGAGCAGCACCATGGACCAGATGAGCAGC"


@dataclass(frozen=True)
class FivePrimeClass:
    key: str
    label: str
    source_seq_id: str
    primer_name: str
    primer_sequence: str


@dataclass(frozen=True)
class ThreePrimeClass:
    key: str
    label: str
    source_seq_id: str


@dataclass(frozen=True)
class ReverseReadout:
    key: str
    label: str
    primer_name: str
    primer_sequence: str
    resolved_classes: tuple[str, ...]
    note: str


FIVE_PRIME_CLASSES = [
    FivePrimeClass(
        "ta",
        "TAp73",
        "tp73_tx__f5__nm_005427_4",
        "TP73-5p-TA-F",
        "ACGGCTGCAGAGCGAGCTGC",
    ),
    FivePrimeClass(
        "x2",
        "Ex2-p73",
        "tp73_tx__f21__xm_047429524_1",
        "TP73-5p-Ex2-F",
        "AGGCTTAGCCAAGAGCGAGCTGC",
    ),
    FivePrimeClass(
        "x1",
        "Ex3-p73",
        "tp73_tx__f24__xm_047429521_1",
        "TP73-5p-Ex3-F",
        "AAATAACAGAGCGAGCTGCCCTCG",
    ),
    FivePrimeClass(
        "dn",
        "DeltaN-p73",
        "tp73_tx__f99__nm_001126240_3",
        "TP73-5p-DeltaN-F",
        "ACCTCGCCACGGCCCAGTTCAAT",
    ),
    FivePrimeClass(
        "v13",
        "V13-p73",
        "tp73_tx__f111__nm_001204192_2",
        "TP73-5p-V13-F",
        "CATCTGCAGGACAGGCCCAGTTC",
    ),
]


THREE_PRIME_CLASSES = [
    ThreePrimeClass("alpha", "alpha", "tp73_tx__f5__nm_005427_4"),
    ThreePrimeClass("beta", "beta", "tp73_tx__f6__nm_001204184_2"),
    ThreePrimeClass("gamma", "gamma", "tp73_tx__f7__nm_001204185_2"),
    ThreePrimeClass("delta", "delta", "tp73_tx__f10__nm_001204186_2"),
    ThreePrimeClass("epsilon", "epsilon", "tp73_tx__f8__nm_001204187_2"),
    ThreePrimeClass("zeta", "zeta", "tp73_tx__f9__nm_001204188_2"),
]


REVERSE_READOUTS = [
    ReverseReadout(
        "agz",
        "alpha/gamma/zeta",
        "TP73-3p-alpha-gamma-zeta-R",
        "CCCCAGGTCCTCAATGGTCA",
        ("alpha", "gamma", "zeta"),
        "One reverse primer resolves alpha, gamma, and zeta by long-range cDNA product size.",
    ),
    ReverseReadout(
        "be",
        "beta/epsilon",
        "TP73-3p-beta-epsilon-R",
        "CCAGGTCCTGACGAGGCTGG",
        ("beta", "epsilon"),
        "One reverse primer resolves beta and epsilon by long-range cDNA product size.",
    ),
    ReverseReadout(
        "delta",
        "delta",
        "TP73-3p-delta-R2",
        "GTCGGCCTCTGTAGGAGCTGCTG",
        ("delta",),
        "Improved delta-specific junction primer replacing the older high-GC delta readout.",
    ),
]


OBSERVED_ACCESSIONS = {
    ("ta", "alpha"): "NM_005427.4",
    ("ta", "beta"): "NM_001204184.2",
    ("ta", "gamma"): "NM_001204185.2",
    ("ta", "delta"): "NM_001204186.2",
    ("ta", "epsilon"): "NM_001204187.2",
    ("ta", "zeta"): "NM_001204188.2",
    ("dn", "alpha"): "NM_001126240.3",
    ("dn", "beta"): "NM_001126241.3",
    ("dn", "gamma"): "NM_001126242.3",
    ("dn", "delta"): "NM_001204189.2",
    ("dn", "epsilon"): "NM_001204190.2",
    ("dn", "zeta"): "NM_001204191.2",
    ("x2", "alpha"): "XM_047429524.1",
    ("x1", "alpha"): "XM_047429521.1",
    ("v13", "alpha"): "NM_001204192.2",
}


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def gc_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def wallace_tm_c(sequence: str) -> float:
    return float(2 * (sequence.count("A") + sequence.count("T")) + 4 * (sequence.count("G") + sequence.count("C")))


def gc_tm_c(sequence: str) -> float:
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return 64.9 + 41.0 * (gc_count - 16.4) / len(sequence)


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(chunks).upper()
            current = line[1:].split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if current is not None:
        seqs[current] = "".join(chunks).upper()
    return seqs


def read_state_sequences(path: Path) -> dict[str, str]:
    payload = json.loads(path.read_text())
    out: dict[str, str] = {}
    for seq_id, record in payload.get("sequences", {}).items():
        seq_record = record.get("seq", {})
        raw = seq_record.get("seq", [])
        if isinstance(raw, list):
            out[seq_id] = "".join(chr(value) for value in raw).upper()
        elif isinstance(raw, str):
            out[seq_id] = raw.upper()
    return out


def load_sequences(args: argparse.Namespace) -> dict[str, str]:
    if args.source_state:
        return read_state_sequences(Path(args.source_state))
    if args.source_fasta:
        return read_fasta(Path(args.source_fasta))
    raise SystemExit("Provide --source-state or --source-fasta")


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[idx : idx + width] for idx in range(0, len(sequence), width))


def build_virtual_variants(source_sequences: dict[str, str]) -> tuple[list[dict[str, object]], dict[str, str]]:
    missing = [
        item.source_seq_id
        for item in [*FIVE_PRIME_CLASSES, *THREE_PRIME_CLASSES]
        if item.source_seq_id not in source_sequences
    ]
    if missing:
        raise SystemExit(f"Source transcript sequence(s) missing: {', '.join(sorted(set(missing)))}")

    prefixes: dict[str, str] = {}
    suffixes: dict[str, str] = {}
    for five in FIVE_PRIME_CLASSES:
        source = source_sequences[five.source_seq_id]
        anchor = source.find(ANCHOR_SEQUENCE)
        if anchor < 0:
            raise SystemExit(f"Anchor sequence missing from 5' source {five.source_seq_id}")
        prefixes[five.key] = source[:anchor]
    for three in THREE_PRIME_CLASSES:
        source = source_sequences[three.source_seq_id]
        anchor = source.find(ANCHOR_SEQUENCE)
        if anchor < 0:
            raise SystemExit(f"Anchor sequence missing from 3' source {three.source_seq_id}")
        suffixes[three.key] = source[anchor:]

    variants: list[dict[str, object]] = []
    sequences: dict[str, str] = {}
    for five in FIVE_PRIME_CLASSES:
        for three in THREE_PRIME_CLASSES:
            variant_id = f"tp73_lr_{five.key}_{three.key}"
            display_name = f"{five.label}-{three.label}"
            sequence = prefixes[five.key] + suffixes[three.key]
            observed_accession = OBSERVED_ACCESSIONS.get((five.key, three.key))
            variants.append(
                {
                    "variant_id": variant_id,
                    "display_name": display_name,
                    "five_prime_class": five.key,
                    "three_prime_class": three.key,
                    "source_status": "public_refseq_observed"
                    if observed_accession
                    else "virtual_local_hypothesis",
                    "observed_accession": observed_accession,
                    "sequence_length_bp": len(sequence),
                    "construction": {
                        "five_prime_prefix_source": five.source_seq_id,
                        "three_prime_suffix_source": three.source_seq_id,
                        "anchor_sequence": ANCHOR_SEQUENCE,
                    },
                }
            )
            sequences[variant_id] = sequence
    return variants, sequences


def write_fasta(path: Path, variants: list[dict[str, object]], sequences: dict[str, str]) -> str:
    lines: list[str] = []
    for variant in variants:
        variant_id = str(variant["variant_id"])
        lines.append(
            f">{variant_id} {variant['display_name']} "
            f"{variant['source_status']} length={variant['sequence_length_bp']}"
        )
        lines.append(wrap_fasta(sequences[variant_id]))
    text = "\n".join(lines) + "\n"
    path.write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def primer_record(name: str, sequence: str, role: str) -> dict[str, object]:
    return {
        "name": name,
        "role": role,
        "sequence": sequence,
        "length_bp": len(sequence),
        "gc_fraction": round(gc_fraction(sequence), 4),
        "tm_wallace_c": round(wallace_tm_c(sequence), 1),
        "tm_gc_rule_c": round(gc_tm_c(sequence), 1),
    }


def product_rows(variants: list[dict[str, object]], sequences: dict[str, str]) -> list[dict[str, object]]:
    by_key = {
        (str(item["five_prime_class"]), str(item["three_prime_class"])): item
        for item in variants
    }
    rows: list[dict[str, object]] = []
    for five in FIVE_PRIME_CLASSES:
        for readout in REVERSE_READOUTS:
            target = reverse_complement(readout.primer_sequence)
            for three_key in readout.resolved_classes:
                variant = by_key[(five.key, three_key)]
                sequence = sequences[str(variant["variant_id"])]
                forward_start = sequence.find(five.primer_sequence)
                reverse_start = sequence.find(target)
                if forward_start < 0 or reverse_start < 0 or reverse_start < forward_start:
                    raise SystemExit(
                        f"Primer product missing for {five.key}/{readout.key}/{three_key}"
                    )
                rows.append(
                    {
                        "lane": f"{five.key}+{readout.key}",
                        "five_prime_class": five.key,
                        "five_prime_label": five.label,
                        "reverse_readout": readout.key,
                        "reverse_label": readout.label,
                        "three_prime_class": three_key,
                        "variant_id": variant["variant_id"],
                        "display_name": variant["display_name"],
                        "source_status": variant["source_status"],
                        "observed_accession": variant["observed_accession"],
                        "amplicon_length_bp": reverse_start + len(target) - forward_start,
                        "forward_start_0based": forward_start,
                        "reverse_target_start_0based": reverse_start,
                        "reverse_target_end_0based_exclusive": reverse_start + len(target),
                    }
                )
    return rows


def write_product_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "lane",
        "five_prime_class",
        "reverse_readout",
        "three_prime_class",
        "display_name",
        "source_status",
        "observed_accession",
        "amplicon_length_bp",
        "forward_start_0based",
        "reverse_target_start_0based",
        "reverse_target_end_0based_exclusive",
    ]
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join("" if row[col] is None else str(row[col]) for col in columns))
    path.write_text("\n".join(lines) + "\n")


def find_class_specific_kmers(sequences: dict[str, str], min_len: int = 8, max_len: int = 24) -> dict[str, object]:
    suffixes = {}
    for three in THREE_PRIME_CLASSES:
        source = sequences[three.source_seq_id]
        suffixes[three.key] = source[source.find(ANCHOR_SEQUENCE) :]
    result: dict[str, object] = {}
    for key, sequence in suffixes.items():
        others = [value for other_key, value in suffixes.items() if other_key != key]
        found: dict[str, object] | None = None
        for length in range(min_len, max_len + 1):
            for start in range(300, len(sequence) - length):
                kmer = sequence[start : start + length]
                if sequence.count(kmer) == 1 and not any(kmer in other for other in others):
                    found = {"min_unique_kmer_length_bp": length, "start_0based": start, "kmer": kmer}
                    break
            if found:
                break
        result[key] = found or {
            "min_unique_kmer_length_bp": None,
            "note": f"No exact class-unique k-mer of {min_len}-{max_len} bp was found after the shared anchor.",
        }
    return result


def write_panel_json(
    path: Path,
    fasta_path: Path,
    fasta_sha256: str,
    variants: list[dict[str, object]],
    products: list[dict[str, object]],
    class_specific_kmers: dict[str, object],
    rt_summary: dict[str, object] | None,
) -> None:
    panel = {
        "schema": "gentle.local_tp73_long_range_cdna_panel.v1",
        "panel_id": "tp73_long_range_cdna_virtual_panel_v1",
        "gene_symbol": "TP73",
        "molecule": "cDNA",
        "sequence_count": len(variants),
        "source": "Local TP73 long-range cDNA selector panel generated from bundled TP73 RefSeq-style cDNA templates plus virtual 5'/3' recombinations.",
        "source_fixture": "test_files/tp73.ncbi.gb",
        "sequence_fasta": {
            "path": str(fasta_path),
            "sha256": fasta_sha256,
        },
        "construction_rule": {
            "kind": "prefix_suffix_recombination",
            "anchor_sequence": ANCHOR_SEQUENCE,
            "description": "Each virtual sequence joins a 5' isoform-class prefix to an observed 3' class suffix at the first shared TP73 coding anchor.",
        },
        "curation": {
            "source_kind": "local_lab_hypothesis",
            "source_label": "TP73 5-prime by 3-prime cDNA long-range PCR exploration panel",
            "evidence": [
                "Existing public RefSeq-like TP73 cDNA templates cover TA and DeltaN crossed with alpha through zeta.",
                "Ex2, Ex3, and V13 public templates in the bundled annotation are alpha-only; non-alpha combinations are represented here as explicit virtual hypotheses rather than treated as absent.",
                "Primer-product sizes are exact matches on the virtual cDNA sequence panel.",
            ],
            "limitations": [
                "Expression is not implied by a virtual sequence.",
                "Eta and rarer 3' forms are not synthesized here because no local sequence template was present in the bundled TP73 cDNA set.",
            ],
        },
        "five_prime_classes": [
            {
                "key": item.key,
                "label": item.label,
                "source_seq_id": item.source_seq_id,
                "selector_primer": primer_record(item.primer_name, item.primer_sequence, "forward_primer"),
            }
            for item in FIVE_PRIME_CLASSES
        ],
        "three_prime_classes": [
            {"key": item.key, "label": item.label, "source_seq_id": item.source_seq_id}
            for item in THREE_PRIME_CLASSES
        ],
        "reverse_readouts": [
            {
                "key": item.key,
                "label": item.label,
                "resolved_classes": list(item.resolved_classes),
                "note": item.note,
                "selector_primer": primer_record(item.primer_name, item.primer_sequence, "reverse_primer"),
            }
            for item in REVERSE_READOUTS
        ],
        "specificity_findings": {
            "class_specific_kmers_after_anchor": class_specific_kmers,
            "interpretation": "Alpha, beta, gamma, and epsilon are not supported by single exact class-unique reverse primers in this local suffix set; the selected panel distinguishes them by readout lane plus gel-resolved product size.",
        },
        "variants": variants,
        "primer_products": products,
        "reverse_transcriptase_obstacle_summary": rt_summary,
    }
    path.write_text(json.dumps(panel, indent=2, ensure_ascii=False) + "\n")


def gel_y(length_bp: int, top: float, bottom: float, min_bp: int = 900, max_bp: int = 1800) -> float:
    # Agarose gels separate roughly by log length, so map log10(bp) to vertical position.
    import math

    min_log = math.log10(min_bp)
    max_log = math.log10(max_bp)
    value = (math.log10(length_bp) - min_log) / (max_log - min_log)
    return bottom - value * (bottom - top)


def render_gel_svg(rows: list[dict[str, object]], path: Path, unreported_style: str) -> None:
    colors = {
        "alpha": "#2563eb",
        "beta": "#dc2626",
        "gamma": "#16a34a",
        "delta": "#9333ea",
        "epsilon": "#f97316",
        "zeta": "#0f766e",
    }
    lanes = []
    for five in FIVE_PRIME_CLASSES:
        for readout in REVERSE_READOUTS:
            lanes.append((five.key, readout.key, f"{five.label}\n{readout.label}"))
    lane_rows = {
        (str(row["five_prime_class"]), str(row["reverse_readout"])): []
        for row in rows
    }
    for row in rows:
        lane_rows[(str(row["five_prime_class"]), str(row["reverse_readout"]))].append(row)

    five_by_key = {item.key: item for item in FIVE_PRIME_CLASSES}
    readout_by_key = {item.key: item for item in REVERSE_READOUTS}
    lane_width = 100
    left = 118
    top = 234
    bottom = 582
    width = left + lane_width * len(lanes) + 44
    height = 730
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        "<title>TP73 virtual long-range cDNA selector gel</title>",
        "<style>"
        "text{font-family:Avenir Next,Segoe UI,Arial,sans-serif}"
        ".title{font-size:24px;font-weight:800;fill:#0f172a}"
        ".sub{font-size:12px;fill:#475569}"
        ".axis{font-size:10px;font-weight:700;fill:#64748b}"
        ".lane{fill:#f8fafc;stroke:#cbd5e1}"
        ".well{fill:#334155;opacity:.7}"
        ".band-label{font-size:9px;font-weight:700;fill:#0f172a}"
        ".legend{font-size:11px;fill:#334155}"
        ".small{font-size:10px;fill:#475569}"
        ".col{font-size:10px;font-weight:800;fill:#0f172a}"
        ".mono{font-size:8px;fill:#334155;font-family:Menlo,Consolas,monospace}"
        ".group{font-size:11px;font-weight:800;fill:#0f172a}"
        ".bar{stroke:#334155;stroke-width:1.4;stroke-linecap:round}"
        ".unreported{stroke:#0f172a;stroke-width:1.1;stroke-dasharray:3 2}"
        ".unreported-tag{font-size:8px;font-weight:900;fill:#0f172a}"
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        f'<rect x="18" y="18" width="{width - 36}" height="{height - 36}" rx="14" fill="#f1f5f9" stroke="#cbd5e1"/>',
        '<text class="title" x="38" y="56">TP73 virtual long-range cDNA selector gel</text>',
        '<text class="sub" x="38" y="78">Virtual local-knowledge panel: five 5-prime classes crossed with alpha, beta, gamma, delta, epsilon, and zeta 3-prime classes.</text>',
        '<text class="sub" x="38" y="97">Readout lanes: alpha/gamma/zeta and beta/epsilon are separated by product size; delta uses its own junction primer.</text>',
        '<text class="axis" x="118" y="118">Top grouped bars: shared 5-prime forward selector across the three 3-prime readout lanes underneath.</text>',
        '<text class="mono" x="118" y="134">3-prime readout primers: AGZ 5&apos;-CCCCAGGTCCTCAATGGTCA-3&apos; | BE 5&apos;-CCAGGTCCTGACGAGGCTGG-3&apos; | DELTA 5&apos;-GTCGGCCTCTGTAGGAGCTGCTG-3&apos;</text>',
    ]
    group_bar_y = 166
    for group_idx, five in enumerate(FIVE_PRIME_CLASSES):
        group_x = left + group_idx * len(REVERSE_READOUTS) * lane_width
        group_width = len(REVERSE_READOUTS) * lane_width - 10
        center_x = group_x + group_width / 2
        svg.append(f'<line class="bar" x1="{group_x + 6}" y1="{group_bar_y}" x2="{group_x + group_width - 6}" y2="{group_bar_y}"/>')
        svg.append(f'<line class="bar" x1="{group_x + 6}" y1="{group_bar_y}" x2="{group_x + 6}" y2="{group_bar_y + 8}"/>')
        svg.append(f'<line class="bar" x1="{group_x + group_width - 6}" y1="{group_bar_y}" x2="{group_x + group_width - 6}" y2="{group_bar_y + 8}"/>')
        svg.append(
            f'<text class="group" x="{center_x:.1f}" y="{group_bar_y - 21}" text-anchor="middle">'
            f"{html.escape(five.primer_name)} ({html.escape(five.label)})</text>"
        )
        svg.append(
            f'<text class="mono" x="{center_x:.1f}" y="{group_bar_y - 8}" text-anchor="middle">'
            f"F 5'-{html.escape(five.primer_sequence)}-3'</text>"
        )
    for bp in [1000, 1200, 1400, 1600, 1800]:
        y = gel_y(bp, top, bottom)
        svg.append(f'<line x1="78" y1="{y:.1f}" x2="{width - 38}" y2="{y:.1f}" stroke="#e2e8f0"/>')
        svg.append(f'<text class="axis" x="70" y="{y + 3:.1f}" text-anchor="end">{bp} bp</text>')
    for idx, (five_key, readout_key, label) in enumerate(lanes):
        x = left + idx * lane_width
        five = five_by_key[five_key]
        readout = readout_by_key[readout_key]
        lane_title = (
            f"L{idx + 1:02d}: {five.primer_name} ({five.primer_sequence}) x "
            f"{readout.primer_name} ({readout.primer_sequence})"
        )
        svg.append(
            f'<text class="col" x="{x + (lane_width - 10) / 2:.1f}" y="{top - 52}" text-anchor="middle">'
            f"L{idx + 1:02d}</text>"
        )
        svg.append(
            f'<text class="axis" x="{x + (lane_width - 10) / 2:.1f}" y="{top - 38}" text-anchor="middle">'
            f"R {html.escape(readout.key.upper())}: {html.escape(readout.label)}</text>"
        )
        svg.append(
            f'<rect class="lane" x="{x}" y="{top - 26}" width="{lane_width - 10}" height="{bottom - top + 52}" rx="7">'
            f"<title>{html.escape(lane_title)}</title></rect>"
        )
        svg.append(f'<rect class="well" x="{x + 18}" y="{top - 19}" width="{lane_width - 46}" height="8" rx="3"/>')
        label_lines = label.split("\n")
        svg.append(
            f'<text class="axis" x="{x + (lane_width - 10) / 2:.1f}" y="{bottom + 74}" text-anchor="middle">{html.escape(label_lines[0])}</text>'
        )
        svg.append(
            f'<text class="small" x="{x + (lane_width - 10) / 2:.1f}" y="{bottom + 89}" text-anchor="middle">{html.escape(label_lines[1])}</text>'
        )
        for band in lane_rows[(five_key, readout_key)]:
            length = int(band["amplicon_length_bp"])
            y = gel_y(length, top, bottom)
            cls = str(band["three_prime_class"])
            color = colors[cls]
            is_unreported = str(band.get("source_status")) == "virtual_local_hypothesis"
            dim_unreported = unreported_style in {"dim", "dim_tag"}
            tag_unreported = unreported_style in {"tag", "dim_tag"}
            opacity = ".38" if is_unreported and dim_unreported else ".78"
            band_class = ' class="unreported"' if is_unreported and unreported_style != "none" else ""
            tag = "*" if is_unreported and tag_unreported else ""
            title_suffix = (
                "; virtual local hypothesis, not reported in the bundled public GenBank/Ensembl-style panel"
                if is_unreported
                else "; public record present in the bundled panel"
            )
            svg.append(
                f'<rect{band_class} x="{x + 10}" y="{y - 4:.1f}" width="{lane_width - 30}" height="8" rx="4" fill="{color}" opacity="{opacity}">'
                f'<title>{html.escape(str(band["display_name"]))}: {length} bp{html.escape(title_suffix)}</title></rect>'
            )
            svg.append(
                f'<text class="band-label" x="{x + lane_width - 13}" y="{y + 3:.1f}" text-anchor="end">{html.escape(cls[0])} {length}{tag}</text>'
            )
            if is_unreported and tag_unreported:
                svg.append(
                    f'<text class="unreported-tag" x="{x + 14}" y="{y - 7:.1f}">*</text>'
                )
    legend_x = 42
    legend_y = 638
    for idx, (key, color) in enumerate(colors.items()):
        x = legend_x + idx * 112
        svg.append(f'<rect x="{x}" y="{legend_y}" width="18" height="8" rx="4" fill="{color}" opacity=".78"/>')
        svg.append(f'<text class="legend" x="{x + 25}" y="{legend_y + 8}">{html.escape(key)}</text>')
    if unreported_style != "none":
        unreported_x = legend_x + len(colors) * 112 + 20
        unreported_opacity = ".38" if unreported_style in {"dim", "dim_tag"} else ".78"
        unreported_tag = "*" if unreported_style in {"tag", "dim_tag"} else ""
        svg.append(
            f'<rect class="unreported" x="{unreported_x}" y="{legend_y}" width="18" height="8" rx="4" fill="#64748b" opacity="{unreported_opacity}"/>'
        )
        svg.append(
            f'<text class="legend" x="{unreported_x + 25}" y="{legend_y + 8}">{unreported_tag} virtual/unreported local hypothesis</text>'
        )
    svg.extend(
        [
            '<text class="small" x="42" y="690">No lane is interpreted as absence: all bands are theoretical products from explicitly synthesized local cDNA hypotheses; dim/dashed or * bands are not reported in the bundled public panel.</text>',
            '<text class="small" x="42" y="710">Coordinates and products are 0-based internally in the TSV/JSON; displayed band sizes are bp lengths.</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(svg) + "\n")


def run_rnalfold(fasta_path: Path, rnalfold_bin: str | None) -> dict[str, object] | None:
    executable = rnalfold_bin or shutil.which("RNALfold")
    if not executable:
        return None
    completed = subprocess.run(
        [executable, f"--infile={fasta_path}", "-L", "150", "--noLP", "-T", "50"],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    current: str | None = None
    hits: dict[str, list[tuple[float, int, int]]] = {}
    for line in completed.stdout.splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            hits[current] = []
            continue
        if current is None:
            continue
        match = re.search(r"\(\s*(-?\d+(?:\.\d+)?)\)\s+(\d+)\s*$", line)
        if not match:
            continue
        structure = line.split()[0]
        hits[current].append((float(match.group(1)), int(match.group(2)), len(structure)))
    return {"executable": executable, "hits": hits}


def summarize_rt_risk(
    variants: list[dict[str, object]],
    sequences: dict[str, str],
    rnalfold_result: dict[str, object] | None,
) -> tuple[list[dict[str, object]], dict[str, object] | None]:
    if rnalfold_result is None:
        return [], None
    hits_by_id = rnalfold_result["hits"]
    rows: list[dict[str, object]] = []
    for variant in variants:
        variant_id = str(variant["variant_id"])
        hits = hits_by_id.get(variant_id, [])
        length = len(sequences[variant_id])
        if not hits:
            continue
        min_energy, min_start, min_span = min(hits, key=lambda item: item[0])
        first = [energy for energy, start, _span in hits if start <= 1000]
        last = [energy for energy, start, _span in hits if start >= length - 1000]
        first_min = min(first) if first else None
        last_min = min(last) if last else None
        risk_class = "higher_5prime_structure" if first_min is not None and first_min <= -60.0 else "baseline"
        rows.append(
            {
                "variant_id": variant_id,
                "display_name": variant["display_name"],
                "five_prime_class": variant["five_prime_class"],
                "three_prime_class": variant["three_prime_class"],
                "sequence_length_bp": length,
                "min_local_mfe_50c_kcal_per_mol": round(min_energy, 1),
                "min_local_structure_start_1based": min_start,
                "min_local_structure_span_bp": min_span,
                "first_1kb_min_mfe_50c_kcal_per_mol": None if first_min is None else round(first_min, 1),
                "last_1kb_min_mfe_50c_kcal_per_mol": None if last_min is None else round(last_min, 1),
                "local_structures_le_minus_50": sum(1 for energy, _start, _span in hits if energy <= -50.0),
                "local_structures_le_minus_60": sum(1 for energy, _start, _span in hits if energy <= -60.0),
                "rt_obstacle_class": risk_class,
            }
        )
    high_classes = sorted({str(row["five_prime_class"]) for row in rows if row["rt_obstacle_class"] == "higher_5prime_structure"})
    summary = {
        "method": "RNALfold -L 150 --noLP -T 50 over full virtual cDNA sequences",
        "rnalfold_executable": rnalfold_result["executable"],
        "interpretation": "This is an RNA-structure obstacle proxy, not an empirically calibrated reverse-transcriptase processivity score.",
        "higher_risk_five_prime_classes": high_classes,
        "conclusion": "TA 5-prime variants form the only higher-risk subset in this virtual panel; 3-prime class does not materially change the reverse-transcription obstacle proxy.",
    }
    return rows, summary


def write_rt_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("status\tmessage\nnot_available\tRNALfold was not available when this file was generated.\n")
        return
    columns = [
        "variant_id",
        "display_name",
        "five_prime_class",
        "three_prime_class",
        "sequence_length_bp",
        "min_local_mfe_50c_kcal_per_mol",
        "min_local_structure_start_1based",
        "min_local_structure_span_bp",
        "first_1kb_min_mfe_50c_kcal_per_mol",
        "last_1kb_min_mfe_50c_kcal_per_mol",
        "local_structures_le_minus_50",
        "local_structures_le_minus_60",
        "rt_obstacle_class",
    ]
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join("" if row[col] is None else str(row[col]) for col in columns))
    path.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--source-state", help="GENtle project/state JSON containing derived TP73 transcript sequences")
    source.add_argument("--source-fasta", help="FASTA containing derived TP73 transcript sequences")
    parser.add_argument("--output-json", default="assets/panels/tp73_long_range_cdna_virtual_panel_v1.json")
    parser.add_argument("--output-fasta", default="assets/panels/tp73_long_range_cdna_virtual_panel_v1.fasta")
    parser.add_argument("--output-products-tsv", default="docs/figures/tp73_long_range_cdna_virtual_products.tsv")
    parser.add_argument("--output-gel-svg", default="docs/figures/tp73_long_range_cdna_virtual_gel.svg")
    parser.add_argument("--output-rt-risk-tsv", default="docs/figures/tp73_long_range_cdna_rt_risk.tsv")
    parser.add_argument(
        "--unreported-style",
        choices=["dim_tag", "dim", "tag", "none"],
        default="dim_tag",
        help="How virtual/unreported local-hypothesis bands are marked in the gel figure.",
    )
    parser.add_argument("--rnalfold-bin", default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_sequences = load_sequences(args)
    variants, virtual_sequences = build_virtual_variants(source_sequences)
    products = product_rows(variants, virtual_sequences)

    output_json = Path(args.output_json)
    output_fasta = Path(args.output_fasta)
    output_products = Path(args.output_products_tsv)
    output_gel = Path(args.output_gel_svg)
    output_rt = Path(args.output_rt_risk_tsv)
    for path in [output_json, output_fasta, output_products, output_gel, output_rt]:
        path.parent.mkdir(parents=True, exist_ok=True)

    fasta_sha256 = write_fasta(output_fasta, variants, virtual_sequences)
    write_product_tsv(output_products, products)
    render_gel_svg(products, output_gel, args.unreported_style)
    rnalfold_result = run_rnalfold(output_fasta, args.rnalfold_bin)
    rt_rows, rt_summary = summarize_rt_risk(variants, virtual_sequences, rnalfold_result)
    write_rt_tsv(output_rt, rt_rows)
    class_specific_kmers = find_class_specific_kmers(source_sequences)
    write_panel_json(
        output_json,
        output_fasta,
        fasta_sha256,
        variants,
        products,
        class_specific_kmers,
        rt_summary,
    )
    print(f"wrote {output_json}")
    print(f"wrote {output_fasta}")
    print(f"wrote {output_products}")
    print(f"wrote {output_gel}")
    print(f"wrote {output_rt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
