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
import subprocess
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

try:
    from . import compare_candidates_to_promoterome as comparison_tools
    from .render_tp73_cutrun_promoter_comparison import ordered_blocks
except ImportError:
    import compare_candidates_to_promoterome as comparison_tools
    from render_tp73_cutrun_promoter_comparison import ordered_blocks

require = comparison_tools.require


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


def load_bound_comparison(
    candidates_path: Path,
    comparison_path: Path,
    matches_path: Path,
    hits_path: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, dict[str, Any]],
           dict[str, list[dict[str, str]]], dict[tuple[str, str], list[dict[str, Any]]],
           dict[str, dict[str, Any]], dict[str, float]]:
    """Load only corrected, hash-bound comparison evidence and its retained HSPs."""
    candidates = json.loads(candidates_path.read_text(encoding="utf-8"))
    comparison = json.loads(comparison_path.read_text(encoding="utf-8"))
    require(comparison.get("schema") == comparison_tools.SCHEMA
            and comparison.get("counting_policy_id") == comparison_tools.COUNTING_POLICY,
            "comparison lacks corrected counting/provenance; rerun it")
    fasta_path = candidates_path.with_suffix(".fa")
    require(fasta_path.is_file(), "candidate FASTA is missing beside candidate JSON")
    representatives = comparison_tools.validate_candidate_sequences(candidates, fasta_path)
    require(comparison_tools.same_digest(sha256(candidates_path), comparison["candidate_metadata_sha256"])
            and comparison_tools.same_digest(sha256(fasta_path), comparison["query_fasta_sha256"]),
            "comparison candidate hash mismatch")
    require(comparison_tools.same_digest(
        candidates["source_bindings"]["promoterome_receipt_sha256"],
        comparison["promoterome_receipt_sha256"]), "candidate/reference receipt mismatch")

    regions = {row["region_id"]: row for row in candidates["regions"]}
    require(regions and len(regions) == len(candidates["regions"]),
            "missing or duplicate candidate regions")
    require(set(comparison["tasks"]) == {"blastn"}, "expected one corrected BLASTN task")
    task = comparison["tasks"]["blastn"]
    summaries = {row["query_id"]: row for row in task["queries"]}
    require(set(summaries) == set(regions) and len(summaries) == len(task["queries"]),
            "comparison query mismatch")
    run_directory = comparison_path.parent.resolve()
    declared_matches = (run_directory / task["matches_path"]).resolve(strict=True)
    declared_hits = (run_directory / task["raw_hits_path"]).resolve(strict=True)
    require(matches_path.resolve(strict=True) == declared_matches
            and hits_path.resolve(strict=True) == declared_hits,
            "supplied tables are not the comparison-bound artifacts")
    require(comparison_tools.same_digest(sha256(declared_matches), task["matches_sha256"])
            and comparison_tools.same_digest(sha256(declared_hits), task["raw_hits_sha256"]),
            "comparison table hash mismatch")

    for query_id, summary in summaries.items():
        require(summary["query_length_bp"] == regions[query_id]["sequence_length_bp"],
                "candidate/comparison length mismatch")
        require(summary.get("other_gene_exclusion_verified") is True,
                "query gene identity is unresolved")
        require(all(type(summary.get(field)) is bool for field in
                    ("target_cap_reached", "hsp_cap_reached", "counts_are_lower_bounds")),
                "missing search-limit audit")
        require(summary["counts_are_lower_bounds"] ==
                (summary["target_cap_reached"] or summary["hsp_cap_reached"]),
                "inconsistent search-limit audit")

    matches_by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
    with declared_matches.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames and "same_gene" in reader.fieldnames,
                "matches lack corrected same-gene policy")
        for row in reader:
            require(row["query_id"] in regions, "unknown match query")
            require(row["same_locus_overlap"] in {"True", "False"}
                    and row["same_gene"] in {"True", "False"}, "invalid exclusion flags")
            if row["same_locus_overlap"] == "False" and row["same_gene"] == "False":
                matches_by_query[row["query_id"]].append(row)

    selected_pairs = {
        (query_id, row["promoter_id"])
        for query_id, rows in matches_by_query.items()
        for row in sorted(rows, key=lambda item: (
            -float(item["aligned_query_fraction"]), -float(item["best_bitscore"]),
            item["promoter_id"]))[:5]
    }
    members: dict[str, list[str]] = defaultdict(list)
    for member, representative in representatives.items():
        members[representative].append(member)
    thresholds = comparison["thresholds"]
    minimum_bp = int(thresholds["min_alignment_bp"])
    minimum_identity = float(thresholds["min_identity_pct"])
    maximum_evalue = float(thresholds["max_evalue"])
    require(minimum_bp > 0 and 0 <= minimum_identity <= 100 and maximum_evalue >= 0,
            "invalid comparison thresholds")
    hsps: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in comparison_tools.iter_blast_rows(declared_hits):
        require(row["qseqid"] in members, "unknown HSP representative")
        if comparison_tools.hsp_passes(row, minimum_bp, minimum_identity, maximum_evalue):
            for member in members[row["qseqid"]]:
                pair = (member, row["sseqid"])
                if pair in selected_pairs:
                    length = regions[member]["sequence_length_bp"]
                    require(1 <= min(row["qstart"], row["qend"])
                            <= max(row["qstart"], row["qend"]) <= length,
                            "HSP query coordinate exceeds candidate")
                    hsps[pair].append(row)
    return (candidates, comparison, regions, matches_by_query, hsps, summaries,
            {"minimum_bp": minimum_bp, "minimum_identity": minimum_identity,
             "maximum_evalue": maximum_evalue})


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


def blue(identity: float, minimum_identity: float) -> tuple[float, float, float]:
    fraction = min(1.0, max(0.0, (identity - minimum_identity)
                            / max(1.0, 100.0 - minimum_identity)))
    return ((210 - 170 * fraction) / 255, (230 - 120 * fraction) / 255,
            (250 - 35 * fraction) / 255)


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
                    hsps: dict[tuple[str, str], list[dict[str, Any]]],
                    summaries: dict[str, dict[str, Any]], minimum_bp: int,
                    minimum_identity: float) -> None:
    start, end = stretch["start_1based"], stretch["end_1based"]
    row_count = sum(2 if region["sequence_length_bp"] < minimum_bp else 4 for region in regions)
    ax.set_xlim((end, start) if report["gene_strand"] == "-" else (start, end))
    ax.set_ylim(0, row_count + 0.5)
    ax.set_yticks([])
    ax.set_title("Promoterome recurrence by Ensembl feature intersection",
                 loc="left", fontsize=8.5, fontweight="bold")
    visual_left = start if report["gene_strand"] == "+" else end
    y = row_count - 0.4
    for region in regions:
        query_id = region["region_id"]
        source = region["source_region"]
        x0 = source["interval"]["start_0based"] + 1
        length = region["sequence_length_bp"]
        rows = matches_by_query.get(query_id, [])
        summary = summaries[query_id]
        label = f"{source['region_id']} · {source['feature_type']} · {length} bp"
        if length < minimum_bp:
            ax.add_patch(Rectangle((x0, y - 0.22), length, 0.44,
                                   facecolor=FEATURE_COLOURS.get(source["feature_type"], "#aaa"),
                                   edgecolor="#333", lw=0.5))
            ax.text(visual_left, y,
                    f"{label}: motif-scale; not tested by ≥{minimum_bp:g}-bp fragment search",
                    fontsize=5.9, va="center", ha="left",
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 0.25})
            y -= 2
            continue
        genes = summary["other_promoters"]["distinct_genes"]
        tier_counts = {
            threshold: summary["other_promoters_by_min_query_coverage"][f"{threshold:.2f}"]["distinct_genes"]
            for threshold in (0.25, 0.50, 0.80)
        }
        suffix = ("LOWER BOUNDS" if summary["counts_are_lower_bounds"]
                  else "observed; no cap saturation")
        ax.text(visual_left, y + 0.18,
                f"{label} · other genes any {genes:,}; ≥25% {tier_counts[0.25]:,}; "
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
            blocks = hsps.get((query_id, target["promoter_id"]), [])
            require(blocks, "ranked match lacks qualifying bound HSPs")
            for block_number, (hsp, broken) in enumerate(ordered_blocks(blocks), 1):
                q0 = min(int(hsp["qstart"]), int(hsp["qend"])) - 1
                q1 = max(int(hsp["qstart"]), int(hsp["qend"]))
                midpoint = (q0 + q1) / 2
                ax.add_patch(Rectangle((x0 + q0, y - 0.14), q1 - q0, 0.28,
                                       facecolor=blue(float(hsp["pident"]), minimum_identity),
                                       edgecolor="#d62728" if broken else "#4a6572",
                                       lw=1.0 if broken else 0.35))
                if q1 - q0 >= 12:
                    ax.text(x0 + midpoint, y, str(block_number), fontsize=4.8,
                            ha="center", va="center")
            ax.text(visual_left, y,
                    f"{rank}. {target_label} · {float(target['aligned_query_fraction']):.1%}",
                    fontsize=5.3, va="center", ha="left",
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 0.2})
            y -= 0.7
        if not top:
            ax.text((start + end) / 2, y,
                    "No qualifying other-gene hit reported; not proof of uniqueness",
                    fontsize=5.5, ha="center", va="center", color="#555")
            y -= 0.7
        y -= 0.25
    ax.spines[:].set_visible(False)
    ax.set_xticks([])


def render_page(gene: str, stretch: dict[str, Any], report: dict[str, Any],
                regions: list[dict[str, Any]], matches_by_query: dict[str, list[dict[str, str]]],
                hsps: dict[tuple[str, str], list[dict[str, Any]]],
                summaries: dict[str, dict[str, Any]], thresholds: dict[str, float]) -> Any:
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
                    matches_by_query, hsps, summaries, int(thresholds["minimum_bp"]),
                    thresholds["minimum_identity"])
    fig.text(
        0.01, 0.004,
        f"BLASTN ≥{thresholds['minimum_bp']:g} bp, ≥{thresholds['minimum_identity']:g}% identity, "
        f"E≤{thresholds['maximum_evalue']:g}; self-locus and same-gene targets excluded. "
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
    require(not args.output.exists() or not any(args.output.iterdir()),
            "Output directory must be absent or empty")
    (candidates, comparison, regions_by_id, matches_by_query, hsps, summaries,
     thresholds) = load_bound_comparison(
        args.candidates_json, args.comparison, args.matches, args.hits)
    reports = {}
    for path in args.locus_report:
        report = json.loads(path.read_text())
        reports[report["gene_symbol"]] = report
    summary_rows = []
    top_rows = []
    for query_id, region in sorted(regions_by_id.items()):
        source = region["source_region"]
        query_gene = region["gene_query"]
        rows = matches_by_query.get(query_id, [])
        evidence = summaries[query_id]
        counts = evidence["other_promoters"]
        tiers = evidence["other_promoters_by_min_query_coverage"]
        summary = {
            "gene": query_gene,
            "stretch_id": source["stretch_id"],
            "feature_id": source["region_id"],
            "feature_type": source["feature_type"],
            "start_1based": source["interval"]["start_0based"] + 1,
            "end_1based": source["interval"]["end_0based_exclusive"],
            "length_bp": region["sequence_length_bp"],
            "clipped_to_tss_stretch": source["clipped_to_tss_stretch"],
            "other_promoter_windows_any": counts["distinct_promoter_windows"],
            "other_genes_any": counts["distinct_genes"],
            "other_transcripts_any": counts["distinct_transcripts"],
            "other_genes_25pct": tiers["0.25"]["distinct_genes"],
            "other_genes_50pct": tiers["0.50"]["distinct_genes"],
            "other_genes_80pct": tiers["0.80"]["distinct_genes"],
            "target_cap_reached": evidence["target_cap_reached"],
            "hsp_cap_reached": evidence["hsp_cap_reached"],
            "counts_are_lower_bounds": evidence["counts_are_lower_bounds"],
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
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "feature_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    with (args.output / "top_matches.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "gene", "feature_id", "rank", "target_gene_names", "target_gene_ids",
            "target_transcript_ids", "target_promoter_id", "aligned_query_fraction",
            "query_intervals_0based_half_open", "best_identity_pct",
        ], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(top_rows)
    top_by_feature: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in top_rows:
        top_by_feature[row["feature_id"]].append(row)
    interpretation = [
        "# TSS-local regulatory-feature promoterome comparison",
        "",
        (f"Parameters: BLASTN >={thresholds['minimum_bp']:g} bp, "
         f">={thresholds['minimum_identity']:g}% identity, "
         f"E<={thresholds['maximum_evalue']:g}; same-locus and same-gene windows excluded."),
        "",
        "Counts marked lower bounds reached a configured target or per-target HSP cap. "
        "No reported hit is not proof of uniqueness.",
        "",
        "## Derived observations",
        "",
    ]
    for row in summary_rows:
        qualifier = "lower bounds" if row["counts_are_lower_bounds"] else "observed"
        if row["length_bp"] < thresholds["minimum_bp"]:
            interpretation.append(
                f"- {row['gene']} {row['feature_id']} ({row['feature_type']}, "
                f"{row['length_bp']} bp): below the {thresholds['minimum_bp']:g}-bp "
                "fragment-search gate; displayed as motif-scale and not tested."
            )
            continue
        statement = (
            f"- {row['gene']} {row['feature_id']} ({row['feature_type']}, "
            f"{row['length_bp']} bp): {qualifier} — {row['other_genes_any']} other genes; "
            f">=25% / >=50% / >=80% coverage: {row['other_genes_25pct']} / "
            f"{row['other_genes_50pct']} / {row['other_genes_80pct']}."
        )
        best = top_by_feature.get(row["feature_id"], [])
        if best:
            first = best[0]
            statement += (f" Highest ranked observed window: {first['target_promoter_id']} "
                          f"({first['target_gene_names'] or first['target_gene_ids']}), "
                          f"{float(first['aligned_query_fraction']):.1%} aggregate coverage.")
        else:
            statement += " No qualifying other-gene hit was reported under this search."
        interpretation.append(statement)
    interpretation.extend([
        "",
        "## Interpretation boundary",
        "",
        "Sequence recurrence is structural evidence. Ensembl feature classes, predicted TFBS "
        "and CUT&RUN enrichment do not establish direct binding, reporter activity, or sufficiency.",
    ])
    interpretation_path = args.output / "interpretation.md"
    interpretation_path.write_text("\n".join(interpretation) + "\n", encoding="utf-8")

    combined_path = args.output / "CD44_TGFB1_SERPINE1_tss_local_integrated.pdf"
    gene_pdfs = {gene: PdfPages(args.output / f"{gene}_tss_local_integrated.pdf")
                 for gene in sorted(reports)}
    with PdfPages(combined_path) as combined:
        for stretch in candidates["stretches"]:
            gene = stretch["gene"]
            regions = [row for row in candidates["regions"]
                       if row["source_region"]["stretch_id"] == stretch["stretch_id"]]
            fig = render_page(gene, stretch, reports[gene], regions,
                              matches_by_query, hsps, summaries, thresholds)
            combined.savefig(fig, dpi=180)
            gene_pdfs[gene].savefig(fig, dpi=180)
            fig.savefig(args.output / f"{stretch['stretch_id']}.png", dpi=180)
            plt.close(fig)
    for pdf in gene_pdfs.values():
        pdf.close()
    renderer_revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[1], text=True
    ).strip()
    receipt = {
        "schema": "gentle.tss_local_integrated_regulatory_report_receipt.v1",
        "source_revision": candidates["source_revision"],
        "candidate_source_revision": candidates["source_revision"],
        "renderer_revision": renderer_revision,
        "renderer_sha256": f"sha256:{sha256(Path(__file__))}",
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
                               + [args.output / "feature_summary.tsv", args.output / "top_matches.tsv",
                                  interpretation_path])
        },
        "non_claim": "Similarity, predicted TFBS and CUT&RUN enrichment do not prove reporter sufficiency or direct regulation.",
        "counting_policy_id": comparison["counting_policy_id"],
        "thresholds": comparison["thresholds"],
        "search_limit_status": {
            query_id: {field: summaries[query_id][field] for field in
                       ("target_cap_reached", "hsp_cap_reached", "counts_are_lower_bounds")}
            for query_id in sorted(summaries)
        },
    }
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"combined_pdf": str(combined_path), "pages": len(candidates["stretches"])}, indent=2))


if __name__ == "__main__":
    main()
