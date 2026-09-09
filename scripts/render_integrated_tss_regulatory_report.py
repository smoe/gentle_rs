#!/usr/bin/env python3
"""Render TSS-local Ensembl features, TFBS, CUT&RUN and promoterome similarity."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter


FEATURE_COLOURS = {
    "promoter": "#5b8ff9",
    "enhancer": "#61d9a8",
    "emar": "#c7a0e8",
    "ctcf": "#f6bd5b",
    "open_chromatin_region": "#78d3f8",
}
FACTOR_COLOURS = [
    "#2f5597", "#c00000", "#7030a0", "#548235", "#bf9000",
    "#00b0f0", "#c55a11", "#4472c4", "#a5a5a5", "#70ad47",
    "#5b9bd5", "#ed7d31", "#8064a2", "#9e480e", "#264478",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


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


def blue(identity: float) -> tuple[float, float, float]:
    fraction = min(1.0, max(0.0, (identity - 80.0) / 20.0))
    return ((210 - 170 * fraction) / 255, (230 - 120 * fraction) / 255,
            (250 - 35 * fraction) / 255)


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


def abbreviate_lane(label: str) -> str:
    label = label.replace(" (", "|").split("|", 1)[0]
    return (label.replace("SAOS-2 ", "SAOS ")
            .replace("SK-MEL-29-2 ", "SKMEL ")
            .replace("TAp73alpha", "TA")
            .replace("DNp73beta", "DN"))


def context_axis(ax: Any, stretch: dict[str, Any], report: dict[str, Any],
                 regions: list[dict[str, Any]]) -> None:
    start, end = stretch["start_1based"], stretch["end_1based"]
    ax.set_xlim((end, start) if report["gene_strand"] == "-" else (start, end))
    ax.set_ylim(-0.3, len(regions) + 1.5)
    ax.set_yticks([])
    ax.set_title("Selected TSS windows and all overlapping Ensembl Regulation features",
                 loc="left", fontsize=9, fontweight="bold")
    for window in stretch["tss_windows"]:
        ax.axvspan(window["start_1based"], window["end_1based"], color="#e7e9ed", alpha=0.45)
        ax.axvline(window["tss_1based"], color="#111", lw=1.0)
        transcript_label = ",".join(value.replace("ENST0000", "ENST…")
                                    for value in window["transcript_ids"])
        ax.text(window["tss_1based"], len(regions) + 1.08,
                f"TSS {window['tss_1based']:,}\n{transcript_label}", fontsize=6,
                ha="center", va="bottom")
    for index, region in enumerate(regions):
        source = region["source_region"]
        interval = source["interval"]
        y = len(regions) - index - 0.1
        x0 = interval["start_0based"] + 1
        width = interval["end_0based_exclusive"] - interval["start_0based"]
        colour = FEATURE_COLOURS.get(source["feature_type"], "#aaaaaa")
        ax.add_patch(Rectangle((x0, y - 0.27), width, 0.54,
                               facecolor=colour, edgecolor="#333", lw=0.5))
        clipped = " clipped" if source["clipped_to_tss_stretch"] else ""
        ax.text(start if report["gene_strand"] == "+" else end, y,
                f"{source['region_id']} · {source['feature_type']}{clipped} · {width} bp",
                fontsize=6.2, va="center",
                ha="right" if report["gene_strand"] == "+" else "left",
                bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none", "pad": 0.3})
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.set_xlabel(f"GRCh38 chromosome {stretch['tss_windows'][0]['chromosome']} genomic coordinate",
                  fontsize=7)


def tfbs_axis(ax: Any, stretch: dict[str, Any], report: dict[str, Any]) -> None:
    start, end = stretch["start_1based"], stretch["end_1based"]
    ax.set_xlim((end, start) if report["gene_strand"] == "-" else (start, end))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_title("Candidate/cofactor TFBS predictions — one shared layer for this genomic stretch",
                 loc="left", fontsize=8, fontweight="bold")
    shown = []
    for index, track in enumerate(report["regulatory_score_tracks"]):
        factor = track["label"].split()[0]
        colour = FACTOR_COLOURS[index % len(FACTOR_COLOURS)]
        for site in track["sites"]:
            position = (site["genomic_start_1based"] + site["genomic_end_1based"]) / 2
            if start <= position <= end:
                shown.append((position, factor, site["strand"], colour))
    if not shown:
        ax.text((start + end) / 2, 0.5,
                "No retained top-three JASPAR sites from the 15 prior score tracks fall in this stretch.",
                ha="center", va="center", fontsize=7, color="#555")
    else:
        for rank, (position, factor, strand, colour) in enumerate(sorted(shown)):
            y = 0.3 + 0.22 * (rank % 3)
            marker = ">" if strand == "+" else "<"
            ax.scatter([position], [y], marker=marker, s=34, color=colour, zorder=3)
            ax.text(position, y + 0.13, factor, fontsize=5.5, ha="center", color=colour)
    ax.spines[:].set_visible(False)
    ax.set_xticks([])


def cutrun_axis(ax: Any, stretch: dict[str, Any], report: dict[str, Any]) -> None:
    start, end = stretch["start_1based"], stretch["end_1based"]
    groups = report["occupancy_groups"]
    lanes = [(group, item["lane"]) for group in groups for item in group["lanes"]]
    ax.set_xlim((end, start) if report["gene_strand"] == "-" else (start, end))
    ax.set_ylim(0, len(lanes))
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_title("CUT&RUN BigWig lanes from the prior report (group-normalized display)",
                 loc="left", fontsize=8, fontweight="bold")
    for lane_index, (group, lane) in enumerate(lanes):
        y = len(lanes) - lane_index - 0.8
        is_tp73 = group["group_id"].startswith("tp73")
        colour = "#4f81bd" if is_tp73 else "#d9a441"
        scale = max(float(group.get("group_abs_max_score", 0)), 1.0)
        ax.hlines(y, start, end, color="#dddddd", lw=0.4)
        for interval in lane["intervals"]:
            a = max(start, interval["genomic_start_1based"])
            b = min(end, interval["genomic_end_1based"])
            if a > b:
                continue
            height = 0.62 * min(1.0, abs(float(interval["score"])) / scale)
            ax.add_patch(Rectangle((a, y), b - a + 1, height,
                                   facecolor=colour, edgecolor="none", alpha=0.9))
        label_x = end if report["gene_strand"] == "+" else start
        ax.text(label_x, y + 0.08, abbreviate_lane(lane["display_label"]),
                ha="right" if report["gene_strand"] == "+" else "left",
                va="bottom", fontsize=5.2,
                bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 0.2})
    ax.spines[:].set_visible(False)


def similarity_axis(ax: Any, stretch: dict[str, Any], report: dict[str, Any],
                    regions: list[dict[str, Any]], matches_by_query: dict[str, list[dict[str, str]]],
                    hsps: dict[tuple[str, str], list[dict[str, str]]], complete: dict[str, bool]) -> None:
    start, end = stretch["start_1based"], stretch["end_1based"]
    row_count = sum(2 if region["sequence_length_bp"] < 40 else 4 for region in regions)
    ax.set_xlim((end, start) if report["gene_strand"] == "-" else (start, end))
    ax.set_ylim(0, row_count + 0.5)
    ax.set_yticks([])
    ax.set_title("Promoterome recurrence by Ensembl feature intersection",
                 loc="left", fontsize=8.5, fontweight="bold")
    visual_left = start if report["gene_strand"] == "+" else end
    y = row_count - 0.4
    query_gene = stretch["gene"]
    for region in regions:
        query_id = region["region_id"]
        source = region["source_region"]
        x0 = source["interval"]["start_0based"] + 1
        length = region["sequence_length_bp"]
        rows = [row for row in matches_by_query.get(query_id, [])
                if row["same_locus_overlap"] == "False"
                and query_gene not in {value for value in row["gene_names"].split(";") if value}]
        label = f"{source['region_id']} · {source['feature_type']} · {length} bp"
        if length < 40:
            ax.add_patch(Rectangle((x0, y - 0.22), length, 0.44,
                                   facecolor=FEATURE_COLOURS.get(source["feature_type"], "#aaa"),
                                   edgecolor="#333", lw=0.5))
            ax.text(visual_left, y,
                    f"{label}: motif-scale; not tested by ≥40-bp fragment search",
                    fontsize=5.9, va="center", ha="left",
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 0.25})
            y -= 2
            continue
        targets = {row["promoter_id"] for row in rows}
        genes = {gene for row in rows for gene in row["gene_ids"].split(";") if gene}
        tier_counts = {}
        for threshold in (0.25, 0.50, 0.80):
            tier_counts[threshold] = len({gene for row in rows
                                          if float(row["aligned_query_fraction"]) >= threshold
                                          for gene in row["gene_ids"].split(";") if gene})
        suffix = "complete" if complete.get(query_id, True) else "lower bound"
        ax.text(visual_left, y + 0.18,
                f"{label} · other genes any {len(genes):,}; ≥25% {tier_counts[0.25]:,}; "
                f"≥50% {tier_counts[0.50]:,}; ≥80% {tier_counts[0.80]:,} · {suffix}",
                fontsize=5.7, va="bottom", ha="left",
                bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none", "pad": 0.2})
        ax.add_patch(Rectangle((x0, y - 0.18), length, 0.32,
                               facecolor="#f5f5f5", edgecolor="#aaa", lw=0.45))
        segments = frequency_segments(rows)
        maximum = max([count for _, _, count in segments] or [0])
        for a, b, count in segments:
            fraction = 0 if maximum == 0 else math.log1p(count) / math.log1p(maximum)
            colour = (0.92 - 0.72 * fraction, 0.95 - 0.50 * fraction, 0.99 - 0.12 * fraction)
            ax.add_patch(Rectangle((x0 + a, y - 0.18), b - a, 0.32,
                                   facecolor=colour, edgecolor="none"))
        y -= 0.8
        top = sorted(rows, key=lambda row: (float(row["aligned_query_fraction"]),
                                            float(row["best_bitscore"])), reverse=True)[:2]
        for rank, target in enumerate(top, 1):
            target_genes = [value for value in target["gene_names"].split(";") if value]
            target_label = ",".join(target_genes[:2]) or target["gene_ids"].split(";")[0]
            if len(target_genes) > 2:
                target_label += f" +{len(target_genes) - 2}"
            ax.add_patch(Rectangle((x0, y - 0.16), length, 0.32,
                                   facecolor="#fafafa", edgecolor="#dddddd", lw=0.4))
            ordered = sorted(hsps.get((query_id, target["promoter_id"]), []),
                             key=lambda row: min(int(row["sstart"]), int(row["send"])))
            previous_mid = None
            previous_orientation = None
            for block_number, hsp in enumerate(ordered, 1):
                q0 = min(int(hsp["qstart"]), int(hsp["qend"])) - 1
                q1 = max(int(hsp["qstart"]), int(hsp["qend"]))
                orientation = 1 if int(hsp["send"]) >= int(hsp["sstart"]) else -1
                midpoint = (q0 + q1) / 2
                broken = is_order_break(
                    previous_mid, previous_orientation, midpoint, orientation
                )
                ax.add_patch(Rectangle((x0 + q0, y - 0.14), q1 - q0, 0.28,
                                       facecolor=blue(float(hsp["pident"])),
                                       edgecolor="#d62728" if broken else "#4a6572",
                                       lw=1.0 if broken else 0.35))
                if q1 - q0 >= 12:
                    ax.text(x0 + midpoint, y, str(block_number), fontsize=4.8,
                            ha="center", va="center")
                previous_mid = midpoint
                previous_orientation = orientation
            ax.text(visual_left, y,
                    f"{rank}. {target_label} · {float(target['aligned_query_fraction']):.1%}",
                    fontsize=5.3, va="center", ha="left",
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 0.2})
            y -= 0.7
        if not top:
            ax.text((start + end) / 2, y, "No other-gene hit",
                    fontsize=5.5, ha="center", va="center", color="#555")
            y -= 0.7
        y -= 0.25
    ax.spines[:].set_visible(False)
    ax.set_xticks([])


def render_page(gene: str, stretch: dict[str, Any], report: dict[str, Any],
                regions: list[dict[str, Any]], matches_by_query: dict[str, list[dict[str, str]]],
                hsps: dict[tuple[str, str], list[dict[str, str]]], complete: dict[str, bool]) -> Any:
    fig = plt.figure(figsize=(11.69, 8.27))
    grid = fig.add_gridspec(
        4, 1, height_ratios=[1.45, 0.7, 1.9, 3.15],
        left=0.23, right=0.985, top=0.86, bottom=0.055, hspace=0.52,
    )
    fig.text(
        0.5, 0.975,
        f"{gene}: TP73/CUT&RUN-supported TSS-local regulatory comparison\n"
        f"{stretch['stretch_id']} · transcript-oriented −500/+200 bp windows",
        fontsize=13, fontweight="bold", ha="center", va="top",
    )
    context_axis(fig.add_subplot(grid[0]), stretch, report, regions)
    tfbs_axis(fig.add_subplot(grid[1]), stretch, report)
    cutrun_axis(fig.add_subplot(grid[2]), stretch, report)
    similarity_axis(fig.add_subplot(grid[3]), stretch, report, regions,
                    matches_by_query, hsps, complete)
    fig.text(
        0.01, 0.004,
        "BLASTN ≥40 bp, ≥80% identity, E≤1e−5; self-locus and same-gene targets excluded. "
        "Numbers are target-promoter block order; red outlines mark order/orientation breaks. "
        "TFBS are JASPAR predictions; CUT&RUN is occupancy/enrichment evidence. None proves reporter sufficiency.",
        fontsize=5.8,
    )
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates-json", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--matches", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--locus-report", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    candidates = json.loads(args.candidates_json.read_text())
    comparison = json.loads(args.comparison.read_text())
    regions_by_id = {row["region_id"]: row for row in candidates["regions"]}
    reports = {}
    for path in args.locus_report:
        report = json.loads(path.read_text())
        reports[report["gene_symbol"]] = report
    matches_by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in load_tsv(args.matches):
        matches_by_query[row["query_id"]].append(row)

    selected_pairs = set()
    for query_id, rows in matches_by_query.items():
        query_gene = regions_by_id[query_id]["gene_query"]
        eligible = [row for row in rows if row["same_locus_overlap"] == "False"
                    and query_gene not in set(row["gene_names"].split(";"))]
        for row in sorted(eligible, key=lambda item: (
                float(item["aligned_query_fraction"]), float(item["best_bitscore"])),
                reverse=True)[:2]:
            selected_pairs.add((query_id, row["promoter_id"]))
    hsps: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    with args.hits.open(encoding="utf-8") as handle:
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(columns, values))
            pair = (row["qseqid"], row["sseqid"])
            if pair in selected_pairs and int(row["length"]) >= 40 \
                    and float(row["pident"]) >= 80 and float(row["evalue"]) <= 1e-5:
                hsps[pair].append(row)
    complete = {
        row["query_id"]: not row["target_cap_reached"]
        for row in comparison["tasks"]["blastn"]["queries"]
    }

    summary_rows = []
    top_rows = []
    for query_id, region in sorted(regions_by_id.items()):
        source = region["source_region"]
        query_gene = region["gene_query"]
        rows = [row for row in matches_by_query.get(query_id, [])
                if row["same_locus_overlap"] == "False"
                and query_gene not in set(row["gene_names"].split(";"))]
        genes_any = {gene for row in rows for gene in row["gene_ids"].split(";") if gene}
        summary = {
            "gene": query_gene,
            "stretch_id": source["stretch_id"],
            "feature_id": source["region_id"],
            "feature_type": source["feature_type"],
            "start_1based": source["interval"]["start_0based"] + 1,
            "end_1based": source["interval"]["end_0based_exclusive"],
            "length_bp": region["sequence_length_bp"],
            "clipped_to_tss_stretch": source["clipped_to_tss_stretch"],
            "other_genes_any": len(genes_any),
            "other_genes_25pct": len({gene for row in rows
                                       if float(row["aligned_query_fraction"]) >= 0.25
                                       for gene in row["gene_ids"].split(";") if gene}),
            "other_genes_50pct": len({gene for row in rows
                                       if float(row["aligned_query_fraction"]) >= 0.50
                                       for gene in row["gene_ids"].split(";") if gene}),
            "other_genes_80pct": len({gene for row in rows
                                       if float(row["aligned_query_fraction"]) >= 0.80
                                       for gene in row["gene_ids"].split(";") if gene}),
            "frequency_complete": complete.get(query_id, True),
        }
        summary_rows.append(summary)
        for rank, row in enumerate(sorted(rows, key=lambda item: (
                float(item["aligned_query_fraction"]), float(item["best_bitscore"])),
                reverse=True)[:5], 1):
            top_rows.append({
                "gene": query_gene,
                "feature_id": source["region_id"],
                "rank": rank,
                "target_gene_names": row["gene_names"],
                "target_gene_ids": row["gene_ids"],
                "target_transcript_ids": row["transcript_ids"],
                "target_promoter_id": row["promoter_id"],
                "aligned_query_fraction": row["aligned_query_fraction"],
                "query_intervals_0based_half_open": row["query_intervals_0based_half_open"],
                "best_identity_pct": row["best_identity_pct"],
            })
    with (args.output / "feature_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    with (args.output / "top_matches.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(top_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(top_rows)

    combined_path = args.output / "CD44_TGFB1_SERPINE1_tss_local_integrated.pdf"
    gene_pdfs = {gene: PdfPages(args.output / f"{gene}_tss_local_integrated.pdf")
                 for gene in sorted(reports)}
    with PdfPages(combined_path) as combined:
        for stretch in candidates["stretches"]:
            gene = stretch["gene"]
            regions = [row for row in candidates["regions"]
                       if row["source_region"]["stretch_id"] == stretch["stretch_id"]]
            fig = render_page(gene, stretch, reports[gene], regions,
                              matches_by_query, hsps, complete)
            combined.savefig(fig, dpi=180)
            gene_pdfs[gene].savefig(fig, dpi=180)
            fig.savefig(args.output / f"{stretch['stretch_id']}.png", dpi=180)
            plt.close(fig)
    for pdf in gene_pdfs.values():
        pdf.close()
    receipt = {
        "schema": "gentle.tss_local_integrated_regulatory_report_receipt.v1",
        "source_revision": candidates["source_revision"],
        "inputs": {
            "candidates": f"sha256:{sha256(args.candidates_json)}",
            "comparison": f"sha256:{sha256(args.comparison)}",
            "matches": f"sha256:{sha256(args.matches)}",
            "hits": f"sha256:{sha256(args.hits)}",
            "locus_reports": {path.name: f"sha256:{sha256(path)}" for path in args.locus_report},
        },
        "outputs": {
            path.name: f"sha256:{sha256(path)}"
            for path in sorted(list(args.output.glob("*.pdf")) + list(args.output.glob("*.png"))
                               + [args.output / "feature_summary.tsv", args.output / "top_matches.tsv"])
        },
        "non_claim": "Similarity, predicted TFBS and CUT&RUN enrichment do not prove reporter sufficiency or direct regulation.",
    }
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"combined_pdf": str(combined_path), "pages": len(candidates["stretches"])}, indent=2))


if __name__ == "__main__":
    main()
