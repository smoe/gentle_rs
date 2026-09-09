#!/usr/bin/env python3
"""Render the first TP73-supported promoterome comparison as an auditable SVG."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from html import escape
import hashlib
import json
import math
from pathlib import Path


ROOT = Path.cwd()
RUN = ROOT / "blastn-40bp-80pct"
TOP_N = 5
PLOT_X = 500
PLOT_W = 1000
ROW_H = 17
PANEL_H = 142


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def blue(identity: float) -> str:
    fraction = min(1.0, max(0.0, (identity - 80.0) / 20.0))
    r = round(210 - 170 * fraction)
    g = round(230 - 120 * fraction)
    b = round(250 - 35 * fraction)
    return f"rgb({r},{g},{b})"


def frequency_blue(count: int, maximum: int) -> str:
    if count <= 0 or maximum <= 0:
        return "#f5f5f5"
    fraction = math.log1p(count) / math.log1p(maximum)
    r = round(235 - 195 * fraction)
    g = round(244 - 134 * fraction)
    b = round(252 - 37 * fraction)
    return f"rgb({r},{g},{b})"


def main() -> None:
    global ROOT, RUN, TOP_N
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True,
                        help="Directory containing candidate_regions.json and the comparison task directory")
    parser.add_argument("--task-directory", default="blastn-40bp-80pct")
    parser.add_argument("--top", type=int, default=5)
    args = parser.parse_args()
    ROOT = args.root.resolve()
    RUN = ROOT / args.task_directory
    TOP_N = args.top

    candidates = json.loads((ROOT / "candidate_regions.json").read_text())
    comparison = json.loads((RUN / "comparison.json").read_text())
    regions = {row["region_id"]: row for row in candidates["regions"]}
    summaries = {row["query_id"]: row for row in comparison["tasks"]["blastn"]["queries"]}

    top: dict[str, list[dict[str, str]]] = {query_id: [] for query_id in regions}
    with (RUN / "matches.blastn.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            query_gene = regions[row["query_id"]]["gene_query"]
            if row["same_locus_overlap"] == "True" or query_gene in row["gene_names"].split(";"):
                continue
            top[row["query_id"]].append(row)
    for query_id in top:
        top[query_id] = sorted(
            top[query_id],
            key=lambda row: (float(row["aligned_query_fraction"]), float(row["best_bitscore"])),
            reverse=True,
        )[:TOP_N]

    selected_pairs = {
        (query_id, row["promoter_id"])
        for query_id, rows in top.items()
        for row in rows
    }
    hsps: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    with (RUN / "hits.blastn.tsv").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(columns, fields))
            pair = (row["qseqid"], row["sseqid"])
            if pair not in selected_pairs:
                continue
            if int(row["length"]) >= 40 and float(row["pident"]) >= 80.0 and float(row["evalue"]) <= 1e-5:
                hsps[pair].append(row)

    ordered_queries = sorted(
        regions,
        key=lambda query_id: (
            regions[query_id]["gene_query"],
            regions[query_id]["genome_extraction"]["tss_1based"],
        ),
    )
    height = 230 + PANEL_H * len(ordered_queries) + 100
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="{height}" viewBox="0 0 1600 {height}">',
        '<rect width="1600" height="100%" fill="white"/>',
        '<style>text{font-family:DejaVu Sans,Arial,sans-serif;fill:#202124}.title{font-size:25px;font-weight:700}.sub{font-size:14px}.head{font-size:15px;font-weight:700}.small{font-size:11px}.tiny{font-size:9px}</style>',
        '<text x="50" y="40" class="title">TP73 CUT&amp;RUN-supported TSS windows vs the human promoterome</text>',
        '<text x="50" y="67" class="sub">GRCh38 / Ensembl 116 · transcript-oriented −2,000 to +200 bp · BLASTN ≥40 bp, ≥80% identity, E ≤ 1e−5</text>',
        '<text x="50" y="91" class="sub">Support: at least one TP73 experimental-window mean exceeds its matched GFP-control mean in the same cell line.</text>',
        '<text x="50" y="115" class="sub">Frequency strip: distinct other genes per query position (log scale). Rows: top other-gene promoter windows by aggregate query coverage.</text>',
        '<rect x="50" y="140" width="18" height="12" fill="#d2e6fa"/><text x="75" y="151" class="small">80% identity</text>',
        '<rect x="180" y="140" width="18" height="12" fill="#2878d7"/><text x="205" y="151" class="small">100% identity</text>',
        '<rect x="330" y="138" width="24" height="16" fill="none" stroke="#d62728" stroke-width="2"/><text x="362" y="151" class="small">order/orientation break</text>',
        f'<line x1="{PLOT_X}" y1="174" x2="{PLOT_X + PLOT_W}" y2="174" stroke="#555"/>',
    ]
    for position, label in [(0, "−2000"), (1000, "−1000"), (2000, "TSS"), (2200, "+200")]:
        x = PLOT_X + PLOT_W * position / 2201
        svg.append(f'<line x1="{x:.1f}" y1="168" x2="{x:.1f}" y2="180" stroke="#555"/>')
        svg.append(f'<text x="{x:.1f}" y="164" text-anchor="middle" class="small">{label}</text>')

    summary_rows = []
    top_rows = []
    y = 205
    for query_id in ordered_queries:
        region = regions[query_id]
        summary = summaries[query_id]
        gene = region["gene_query"]
        tss = region["genome_extraction"]["tss_1based"]
        transcripts = region["transcript_ids"]
        deltas = {
            evidence["cell_line"]: max(evidence["experimental_minus_matched_gfp"].values())
            for evidence in region["cutrun_support"]["evidence"]
        }
        max_frequency = max(
            [segment["distinct_genes"] for segment in summary["recurrent_query_segments"]] or [0]
        )
        counts = summary["other_promoters"]
        coverage = summary["other_promoters_by_min_query_coverage"]
        title = f"{gene} · TSS {tss:,} · {', '.join(transcripts)}"
        support = f"Δmean max: SAOS-2 {deltas['SAOS-2']:+.3f}; SK-MEL-29-2 {deltas['SK-MEL-29-2']:+.3f}"
        svg.append(f'<text x="50" y="{y}" class="head">{escape(title)}</text>')
        svg.append(f'<text x="50" y="{y + 17}" class="small">{escape(support)}</text>')
        svg.append(f'<text x="50" y="{y + 34}" class="small">Any tract: {counts["distinct_genes"]:,} genes · ≥25%: {coverage["0.25"]["distinct_genes"]:,} · ≥50%: {coverage["0.50"]["distinct_genes"]:,} · ≥80%: {coverage["0.80"]["distinct_genes"]:,}</text>')
        strip_y = y + 9
        svg.append(f'<rect x="{PLOT_X}" y="{strip_y}" width="{PLOT_W}" height="12" fill="#f5f5f5" stroke="#bbb"/>')
        for segment in summary["recurrent_query_segments"]:
            x = PLOT_X + PLOT_W * segment["query_start_0based"] / 2201
            width = max(0.5, PLOT_W * (segment["query_end_0based_exclusive"] - segment["query_start_0based"]) / 2201)
            colour = frequency_blue(segment["distinct_genes"], max_frequency)
            svg.append(f'<rect x="{x:.2f}" y="{strip_y}" width="{width:.2f}" height="12" fill="{colour}"/>')
        svg.append(f'<line x1="{PLOT_X + PLOT_W * 2000 / 2201:.1f}" y1="{strip_y - 2}" x2="{PLOT_X + PLOT_W * 2000 / 2201:.1f}" y2="{strip_y + 14}" stroke="#111"/>')

        row_y = y + 48
        for rank, target in enumerate(top[query_id], 1):
            names = [name for name in target["gene_names"].split(";") if name]
            genes = names or [name for name in target["gene_ids"].split(";") if name]
            target_label = ",".join(genes[:2])
            if len(genes) > 2:
                target_label += f" +{len(genes) - 2}"
            transcript_count = len([value for value in target["transcript_ids"].split(";") if value])
            label = f"{rank}. {target_label} · chr{target['chromosome']}:{int(target['tss_1based']):,} · {transcript_count} tx · {float(target['aligned_query_fraction']):.1%}"
            svg.append(f'<text x="50" y="{row_y + 11}" class="tiny">{escape(label)}</text>')
            svg.append(f'<rect x="{PLOT_X}" y="{row_y}" width="{PLOT_W}" height="13" fill="#fafafa" stroke="#ddd"/>')
            blocks = hsps.get((query_id, target["promoter_id"]), [])
            ordered = sorted(blocks, key=lambda row: min(int(row["sstart"]), int(row["send"])))
            previous_query_mid = None
            previous_orientation = None
            for block_number, hsp in enumerate(ordered, 1):
                q0 = min(int(hsp["qstart"]), int(hsp["qend"])) - 1
                q1 = max(int(hsp["qstart"]), int(hsp["qend"]))
                x = PLOT_X + PLOT_W * q0 / 2201
                width = max(2.0, PLOT_W * (q1 - q0) / 2201)
                orientation = 1 if int(hsp["send"]) >= int(hsp["sstart"]) else -1
                query_mid = (q0 + q1) / 2
                broken = previous_query_mid is not None and (
                    orientation != previous_orientation or query_mid < previous_query_mid
                )
                stroke = "#d62728" if broken else "#4a6572"
                stroke_width = 2 if broken else 0.6
                svg.append(f'<rect x="{x:.2f}" y="{row_y + 1}" width="{width:.2f}" height="11" fill="{blue(float(hsp["pident"]))}" stroke="{stroke}" stroke-width="{stroke_width}"/>')
                if width >= 12:
                    svg.append(f'<text x="{x + width / 2:.2f}" y="{row_y + 10}" text-anchor="middle" class="tiny">{block_number}</text>')
                previous_query_mid = query_mid
                previous_orientation = orientation
            top_rows.append({
                "query_id": query_id,
                "rank": rank,
                **target,
            })
            row_y += ROW_H
        if not top[query_id]:
            svg.append(f'<text x="{PLOT_X}" y="{row_y + 11}" class="small">No other-gene promoter hit passed the declared threshold.</text>')
        svg.append(f'<line x1="50" y1="{y + PANEL_H - 10}" x2="1550" y2="{y + PANEL_H - 10}" stroke="#ddd"/>')
        summary_rows.append({
            "gene": gene,
            "tss_1based": tss,
            "transcript_ids": ";".join(transcripts),
            "max_saos2_delta": deltas["SAOS-2"],
            "max_skmel29_2_delta": deltas["SK-MEL-29-2"],
            "other_promoter_windows_any": counts["distinct_promoter_windows"],
            "other_genes_any": counts["distinct_genes"],
            "other_transcripts_any": counts["distinct_transcripts"],
            "other_genes_25pct": coverage["0.25"]["distinct_genes"],
            "other_genes_50pct": coverage["0.50"]["distinct_genes"],
            "other_genes_80pct": coverage["0.80"]["distinct_genes"],
            "target_cap_reached": summary["target_cap_reached"],
        })
        y += PANEL_H

    svg.extend([
        f'<text x="50" y="{height - 62}" class="small">Numbers are block order in the target promoter’s 5′→3′ sequence. Red outlines begin a changed-order or reversed-orientation block.</text>',
        f'<text x="50" y="{height - 42}" class="small">Similarity and CUT&amp;RUN enrichment motivate reporter contrasts; neither establishes direct TP73 binding, promoter activity, or fragment sufficiency.</text>',
        f'<text x="50" y="{height - 22}" class="tiny">GENtle {escape(candidates["source_revision"])} · candidate SHA-256 {sha256(ROOT / "candidate_regions.json")} · comparison SHA-256 {sha256(RUN / "comparison.json")}</text>',
        '</svg>',
    ])
    (ROOT / "first_comparison.svg").write_text("\n".join(svg) + "\n")

    with (ROOT / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    with (ROOT / "top_matches.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(top_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(top_rows)

    report = [
        "# First TP73 CUT&RUN-supported promoterome comparison",
        "",
        "This run compares every distinct selected-gene TSS window that passes the declared matched-control TP73 BigWig rule against 389,722 GRCh38/Ensembl-116 transcript-linked promoter windows.",
        "",
        "## Parameters",
        "",
        "- Window: transcript-oriented −2,000/+200 bp around each annotated TSS.",
        "- Support: at least one TAp73α or DNp73β mean exceeds the matched GFP-control mean in the same cell line across that exact window.",
        "- Similarity: BLASTN, ≥40 aligned bp, ≥80% identity, E ≤ 1e−5; dust and soft masking enabled.",
        "- Counting: self-locus overlaps excluded; promoter windows, genes, and transcripts counted separately.",
        "",
        "## First observations",
        "",
        "- The SERPINE1 ENST00000950058 window contains a highly recurrent component: its strongest MAN2A2 promoter matches cover about 47% of the proposed window.",
        "- The two promoter-proximal SERPINE1 windows have many short/repeat-rich hits, but only tens of other genes reach 25% aggregate coverage and none reaches 50%.",
        "- CD44 and TGFB1 have only sparse short matches under this threshold; no other gene reaches 25% aggregate window coverage.",
        "- The CD44 ENST00000428726 window has no qualifying other-gene promoter hit in this BLASTN pass.",
        "",
        "These are sequence-recurrence observations, not reporter recommendations. The next decision layer should intersect the recurrent blocks with localized CUT&RUN peaks, Ensembl Regulation features, motif/module evidence, and planned deletion/split constructs.",
    ]
    (ROOT / "report.md").write_text("\n".join(report) + "\n")
    receipt = {
        "schema": "gentle.tp73_cutrun_supported_promoter_comparison_receipt.v1",
        "source_revision": candidates["source_revision"],
        "inputs": {
            "candidate_regions.json": sha256(ROOT / "candidate_regions.json"),
            "candidate_regions.fa": sha256(ROOT / "candidate_regions.fa"),
            "comparison.json": sha256(RUN / "comparison.json"),
        },
        "outputs": {
            "first_comparison.svg": sha256(ROOT / "first_comparison.svg"),
            "summary.tsv": sha256(ROOT / "summary.tsv"),
            "top_matches.tsv": sha256(ROOT / "top_matches.tsv"),
            "report.md": sha256(ROOT / "report.md"),
        },
        "target_cap_reached": any(row["target_cap_reached"] for row in summary_rows),
    }
    (ROOT / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"svg": str(ROOT / "first_comparison.svg"), "queries": len(ordered_queries)}, indent=2))


if __name__ == "__main__":
    main()
