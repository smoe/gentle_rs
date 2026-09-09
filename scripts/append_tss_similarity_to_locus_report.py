#!/usr/bin/env python3
"""Append TSS-local regulatory-feature promoter similarity to a GENtle locus SVG."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from html import escape
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any


FEATURE_COLOURS = {
    "promoter": "#5b8ff9",
    "enhancer": "#61d9a8",
    "emar": "#c7a0e8",
    "ctcf": "#f6bd5b",
    "open_chromatin_region": "#78d3f8",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def frequency_segments(rows: list[dict[str, str]]) -> list[tuple[int, int, int]]:
    events: dict[int, dict[str, set[str]]] = defaultdict(lambda: {"add": set(), "remove": set()})
    genes_by_target: dict[str, set[str]] = {}
    for row in rows:
        target = row["promoter_id"]
        genes_by_target[target] = {value for value in row["gene_ids"].split(";") if value}
        for value in row["query_intervals_0based_half_open"].split(";"):
            start, end = map(int, value.split("-"))
            events[start]["add"].add(target)
            events[end]["remove"].add(target)
    active: set[str] = set()
    gene_counts: Counter[str] = Counter()
    segments = []
    previous = None
    for position in sorted(events):
        if previous is not None and previous < position and active:
            segments.append((previous, position, len(gene_counts)))
        for target in events[position]["remove"]:
            if target in active:
                active.remove(target)
                for gene in genes_by_target[target]:
                    gene_counts[gene] -= 1
                    if gene_counts[gene] == 0:
                        del gene_counts[gene]
        for target in events[position]["add"]:
            if target not in active:
                active.add(target)
                gene_counts.update(genes_by_target[target])
        previous = position
    return segments


def blue(identity: float) -> str:
    fraction = min(1.0, max(0.0, (identity - 80.0) / 20.0))
    return "#{:02x}{:02x}{:02x}".format(
        round(210 - 170 * fraction), round(230 - 120 * fraction), round(250 - 35 * fraction)
    )


def is_order_break(
    previous_query_mid: float | None,
    previous_orientation: int | None,
    query_mid: float,
    orientation: int,
) -> bool:
    """Return whether target-order traversal breaks one collinear query chain."""
    if previous_query_mid is None:
        return False
    if orientation != previous_orientation:
        return True
    if orientation == 1:
        return query_mid < previous_query_mid
    return query_mid > previous_query_mid


def svg_text(x: float, y: float, text: str, *, size: float = 9,
             fill: str = "#475569", weight: str | None = None,
             family: str = "sans-serif", anchor: str | None = None) -> str:
    attrs = [f'x="{x:.2f}"', f'y="{y:.2f}"', f'font-size="{size}"',
             f'fill="{fill}"', f'font-family="{family}"']
    if weight:
        attrs.append(f'font-weight="{weight}"')
    if anchor:
        attrs.append(f'text-anchor="{anchor}"')
    return f"<text {' '.join(attrs)}>{escape(text)}</text>"


def x_for(genomic: int, stretch: dict[str, Any], strand: str) -> float:
    start, end = stretch["start_1based"], stretch["end_1based"]
    fraction = (genomic - start) / max(1, end - start)
    if strand == "-":
        fraction = 1.0 - fraction
    return 255 + fraction * 795


def append_section(
    base_svg: str,
    gene: str,
    candidates: dict[str, Any],
    matches_by_query: dict[str, list[dict[str, str]]],
    hsps: dict[tuple[str, str], list[dict[str, str]]],
    complete: dict[str, bool],
) -> tuple[str, int]:
    root = re.search(r"<svg[^>]*height=\"(\d+)\"[^>]*viewBox=\"0 0 1400 (\d+)\"[^>]*>", base_svg)
    if not root or root.group(1) != root.group(2):
        raise ValueError("unsupported base SVG geometry")
    old_height = int(root.group(1))
    stretches = [row for row in candidates["stretches"] if row["gene"] == gene]
    regions_by_stretch = {
        stretch["stretch_id"]: [
            row for row in candidates["regions"]
            if row["source_region"]["stretch_id"] == stretch["stretch_id"]
        ]
        for stretch in stretches
    }
    section_height = 104
    for stretch in stretches:
        section_height += 76 + sum(
            54 if region["sequence_length_bp"] < 40
            else 80 if matches_by_query.get(region["region_id"])
            else 58
            for region in regions_by_stretch[stretch["stretch_id"]]
        )
    new_height = old_height + section_height
    opening = root.group(0).replace(f'height="{old_height}"', f'height="{new_height}"')
    opening = opening.replace(f'viewBox="0 0 1400 {old_height}"', f'viewBox="0 0 1400 {new_height}"')
    base_svg = base_svg.replace(root.group(0), opening, 1)
    base_svg = base_svg.replace(
        f'<rect fill="#ffffff" height="{old_height}" width="1400" x="0" y="0"/>',
        f'<rect fill="#ffffff" height="{new_height}" width="1400" x="0" y="0"/>', 1,
    )
    body = [f'<g data-gentle-panel="tss-local-promoter-similarity" data-gentle-gene="{escape(gene)}">']
    y = old_height + 30
    body.append(svg_text(34, y, "TSS-local Ensembl-feature promoterome similarity",
                         size=14, fill="#1f2937", weight="bold"))
    y += 20
    body.append(svg_text(
        34, y,
        "Transcript-oriented −500/+200 bp; features clipped to the displayed stretch; "
        "BLASTN ≥40 bp, ≥80% identity, E≤1e−5; self-locus and same-gene targets excluded.",
        size=9, family="monospace", fill="#64748b",
    ))
    y += 18
    body.append(svg_text(
        34, y,
        "Blue frequency strips count distinct other genes by query position. Numbered blocks give target-promoter order; red outlines mark order/orientation breaks.",
        size=9, fill="#64748b",
    ))
    y += 26
    for stretch in stretches:
        strand = stretch["tss_windows"][0]["strand"]
        regions = regions_by_stretch[stretch["stretch_id"]]
        body.append(svg_text(34, y, stretch["stretch_id"], size=11,
                             fill="#334155", weight="bold"))
        body.append(svg_text(
            255, y,
            f"GRCh38 {stretch['tss_windows'][0]['chromosome']}:{stretch['start_1based']:,}–{stretch['end_1based']:,} ({strand})",
            size=9, family="monospace", fill="#64748b",
        ))
        y += 22
        body.append(f'<line x1="255" x2="1050" y1="{y}" y2="{y}" stroke="#64748b" stroke-width="1"/>')
        for window in stretch["tss_windows"]:
            x = x_for(window["tss_1based"], stretch, strand)
            body.append(f'<line x1="{x:.2f}" x2="{x:.2f}" y1="{y-5}" y2="{y+8}" stroke="#111827" stroke-width="1"/>')
            body.append(svg_text(x, y - 7, f"TSS {window['tss_1based']:,}", size=8,
                                 fill="#334155", family="monospace", anchor="middle"))
        y += 20
        for region in regions:
            query_id = region["region_id"]
            source = region["source_region"]
            length = region["sequence_length_bp"]
            rows = [row for row in matches_by_query.get(query_id, [])
                    if row["same_locus_overlap"] == "False"
                    and gene not in set(row["gene_names"].split(";"))]
            genes = {value for row in rows for value in row["gene_ids"].split(";") if value}
            genes_25 = {value for row in rows if float(row["aligned_query_fraction"]) >= 0.25
                        for value in row["gene_ids"].split(";") if value}
            genes_50 = {value for row in rows if float(row["aligned_query_fraction"]) >= 0.50
                        for value in row["gene_ids"].split(";") if value}
            feature_start = source["interval"]["start_0based"] + 1
            feature_end = source["interval"]["end_0based_exclusive"]
            fx0 = min(x_for(feature_start, stretch, strand), x_for(feature_end, stretch, strand))
            fw = abs(x_for(feature_end, stretch, strand) - x_for(feature_start, stretch, strand))
            body.append(svg_text(
                34, y + 8,
                f"{source['region_id']} · {source['feature_type']} · {length} bp",
                size=9, family="monospace", fill="#374151",
            ))
            body.append(f'<rect x="255" y="{y-3}" width="795" height="16" fill="#f8fafc" stroke="#cbd5e1" stroke-width="0.7"/>')
            body.append(f'<rect x="{fx0:.2f}" y="{y-3}" width="{max(2,fw):.2f}" height="16" fill="{FEATURE_COLOURS.get(source["feature_type"], "#aaaaaa")}" fill-opacity="0.34" stroke="#475569" stroke-width="0.7"/>')
            if length < 40:
                body.append(svg_text(1070, y + 8, "motif-scale; below 40-bp search gate",
                                     size=8, family="monospace", fill="#64748b"))
                y += 54
                continue
            status = "complete" if complete.get(query_id, True) else "lower bound"
            body.append(svg_text(
                1070, y + 8,
                f"other genes {len(genes):,}; ≥25% {len(genes_25):,}; ≥50% {len(genes_50):,}; {status}",
                size=8, family="monospace", fill="#475569",
            ))
            y += 22
            segments = frequency_segments(rows)
            maximum = max([count for _, _, count in segments] or [0])
            body.append(f'<rect x="255" y="{y-6}" width="795" height="12" fill="#f8fafc" stroke="#cbd5e1" stroke-width="0.5"/>')
            for start, end, count in segments:
                fraction = 0 if maximum == 0 else math.log1p(count) / math.log1p(maximum)
                colour = "#{:02x}{:02x}{:02x}".format(
                    round(235 - 190 * fraction), round(242 - 130 * fraction),
                    round(252 - 35 * fraction),
                )
                sx = fx0 + start / length * max(2, fw)
                sw = (end - start) / length * max(2, fw)
                body.append(f'<rect x="{sx:.2f}" y="{y-6}" width="{max(0.7,sw):.2f}" height="12" fill="{colour}"/>')
            body.append(svg_text(1070, y + 4, f"position frequency; max {maximum:,} genes",
                                 size=8, family="monospace", fill="#64748b"))
            y += 20
            top = sorted(rows, key=lambda row: (float(row["aligned_query_fraction"]),
                                                float(row["best_bitscore"])), reverse=True)[:2]
            for rank, target in enumerate(top, 1):
                target_names = [value for value in target["gene_names"].split(";") if value]
                target_label = ",".join(target_names[:2]) or target["gene_ids"].split(";")[0]
                body.append(svg_text(34, y + 5,
                                     f"{rank}. {target_label} · {float(target['aligned_query_fraction']):.1%}",
                                     size=8, family="monospace", fill="#475569"))
                body.append(f'<rect x="255" y="{y-6}" width="795" height="14" fill="#f8fafc" stroke="#e2e8f0" stroke-width="0.5"/>')
                ordered = sorted(hsps.get((query_id, target["promoter_id"]), []),
                                 key=lambda row: min(int(row["sstart"]), int(row["send"])))
                previous_mid = None
                previous_orientation = None
                for number, hsp in enumerate(ordered, 1):
                    q0 = min(int(hsp["qstart"]), int(hsp["qend"])) - 1
                    q1 = max(int(hsp["qstart"]), int(hsp["qend"]))
                    orientation = 1 if int(hsp["send"]) >= int(hsp["sstart"]) else -1
                    midpoint = (q0 + q1) / 2
                    broken = is_order_break(
                        previous_mid, previous_orientation, midpoint, orientation
                    )
                    bx = fx0 + q0 / length * max(2, fw)
                    bw = (q1 - q0) / length * max(2, fw)
                    body.append(f'<rect x="{bx:.2f}" y="{y-5}" width="{max(1,bw):.2f}" height="12" fill="{blue(float(hsp["pident"]))}" stroke="{"#dc2626" if broken else "#475569"}" stroke-width="{1.4 if broken else 0.5}"/>')
                    if bw >= 10:
                        body.append(svg_text(bx + bw / 2, y + 4, str(number), size=7,
                                             fill="#111827", anchor="middle"))
                    previous_mid = midpoint
                    previous_orientation = orientation
                y += 19
            if not top:
                body.append(svg_text(652, y + 4, "No qualifying other-gene promoter match",
                                     size=8, fill="#64748b", anchor="middle"))
                y += 18
            y += 16
        y += 18
    body.append(svg_text(
        34, new_height - 18,
        "Interpretation boundary: sequence recurrence is structural evidence; Ensembl classes, predicted TFBS and CUT&RUN enrichment do not prove reporter activity or sufficiency.",
        size=9, fill="#64748b",
    ))
    body.append("</g>")
    return base_svg.rsplit("</svg>", 1)[0] + "\n" + "\n".join(body) + "\n</svg>\n", new_height


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-svg", type=Path, required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--candidates-json", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--matches", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    args = parser.parse_args()

    candidates = json.loads(args.candidates_json.read_text())
    comparison = json.loads(args.comparison.read_text())
    regions = {row["region_id"]: row for row in candidates["regions"]}
    matches_by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in load_tsv(args.matches):
        if regions[row["query_id"]]["gene_query"] == args.gene:
            matches_by_query[row["query_id"]].append(row)
    selected_pairs = set()
    for query_id, rows in matches_by_query.items():
        eligible = [row for row in rows if row["same_locus_overlap"] == "False"
                    and args.gene not in set(row["gene_names"].split(";"))]
        for row in sorted(eligible, key=lambda item: (
                float(item["aligned_query_fraction"]), float(item["best_bitscore"])),
                reverse=True)[:2]:
            selected_pairs.add((query_id, row["promoter_id"]))
    hsps: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    with args.hits.open(encoding="utf-8") as handle:
        for line in handle:
            row = dict(zip(columns, line.rstrip("\n").split("\t")))
            pair = (row["qseqid"], row["sseqid"])
            if pair in selected_pairs and int(row["length"]) >= 40 \
                    and float(row["pident"]) >= 80 and float(row["evalue"]) <= 1e-5:
                hsps[pair].append(row)
    complete = {row["query_id"]: not row["target_cap_reached"]
                for row in comparison["tasks"]["blastn"]["queries"]}
    output, height = append_section(
        args.base_svg.read_text(), args.gene, candidates, matches_by_query, hsps, complete
    )
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(output, encoding="utf-8")
    print(json.dumps({
        "output_svg": str(args.output_svg),
        "height": height,
        "base_svg_sha256": f"sha256:{sha256(args.base_svg)}",
        "candidate_sha256": f"sha256:{sha256(args.candidates_json)}",
    }, indent=2))


if __name__ == "__main__":
    main()
