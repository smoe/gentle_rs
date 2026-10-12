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
import shutil
import subprocess

try:
    from . import compare_candidates_to_promoterome as comparison_tools
except ImportError:
    import compare_candidates_to_promoterome as comparison_tools

require = comparison_tools.require

PLOT_X = 500
PLOT_W = 1000
ROW_H = 17


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def blue(identity: float, minimum: float) -> str:
    fraction = min(1.0, max(0.0, (identity - minimum) / max(1.0, 100.0 - minimum)))
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


def ordered_blocks(blocks):
    """Mark changes relative to alignment direction, not absolute query order."""
    previous_mid = previous_orientation = None
    for row in sorted(blocks, key=lambda row: (
            min(row["sstart"], row["send"]), max(row["sstart"], row["send"]),
            row["qstart"], row["qend"], -row["bitscore"])):
        orientation = (1 if row["send"] >= row["sstart"] else -1) * (1 if row["qend"] >= row["qstart"] else -1)
        mid = (row["qstart"] + row["qend"]) / 2
        broken = previous_mid is not None and (
            orientation != previous_orientation or orientation * (mid - previous_mid) < 0)
        yield row, broken
        previous_mid, previous_orientation = mid, orientation


def render_png(output: Path, renderer: str) -> dict:
    executable = shutil.which(renderer)
    require(executable is not None, "PNG renderer unavailable; provide the path to rsvg-convert")
    executable = str(Path(executable).resolve())
    version = subprocess.run([executable, "--version"], check=True, capture_output=True, timeout=30)
    command = [executable, "--format=png", "--output=first_comparison.png", "first_comparison.svg"]
    result = subprocess.run(command, cwd=output, check=True, capture_output=True, timeout=120)
    png = output / "first_comparison.png"
    with png.open("rb") as handle:
        require(handle.read(8) == b"\x89PNG\r\n\x1a\n" and png.stat().st_size > 8, "renderer did not produce PNG")
    return {"executable": executable, "executable_sha256": sha256(Path(executable)),
            "version": version.stdout.decode(errors="replace").strip(), "command": command,
            "stdout_sha256": hashlib.sha256(result.stdout).hexdigest(),
            "stderr_sha256": hashlib.sha256(result.stderr).hexdigest()}


def render(root: Path, output: Path, task_directory: str = "blastn-40bp-80pct",
           top_n: int = 5, png_renderer: str | None = None) -> dict:
    require(1 <= top_n <= 20, "--top must be between 1 and 20")
    require(not output.exists() or not any(output.iterdir()), "Output directory must be absent or empty")
    ROOT = root.resolve()
    RUN = ROOT / task_directory
    PANEL_H = 65 + ROW_H * top_n

    candidates = json.loads((ROOT / "candidate_regions.json").read_text(encoding="utf-8"))
    comparison = json.loads((RUN / "comparison.json").read_text(encoding="utf-8"))
    require(comparison.get("schema") == comparison_tools.SCHEMA
            and comparison.get("counting_policy_id") == comparison_tools.COUNTING_POLICY,
            "comparison lacks corrected counting/provenance; rerun the comparison, do not relabel historical evidence")
    require(comparison_tools.same_digest(sha256(ROOT / "candidate_regions.json"), comparison["candidate_metadata_sha256"])
            and comparison_tools.same_digest(sha256(ROOT / "candidate_regions.fa"), comparison["query_fasta_sha256"]),
            "comparison candidate input hash mismatch")
    representatives = comparison_tools.validate_candidate_sequences(candidates, ROOT / "candidate_regions.fa")
    require(comparison_tools.same_digest(candidates["source_bindings"]["promoterome_receipt_sha256"], comparison["promoterome_receipt_sha256"]),
            "candidate/reference receipt mismatch")
    regions = {row["region_id"]: row for row in candidates["regions"]}
    task = comparison["tasks"]["blastn"]
    summaries = {row["query_id"]: row for row in task["queries"]}
    require(set(summaries) == set(regions) and len(summaries) == len(task["queries"]), "comparison query mismatch")
    inputs = {"candidate_regions.json": comparison["candidate_metadata_sha256"],
              "candidate_regions.fa": comparison["query_fasta_sha256"],
              "comparison.json": sha256(RUN / "comparison.json")}
    paths = {}
    for kind in ("raw_hits", "matches"):
        path = (RUN / task[kind + "_path"]).resolve(strict=True)
        require(path.is_relative_to(RUN.resolve()) and comparison_tools.same_digest(sha256(path), task[kind + "_sha256"]),
                f"comparison {kind} hash/path mismatch")
        paths[kind] = path
        inputs[kind] = task[kind + "_sha256"]
    background = comparison["background"]
    upstream, downstream = background["upstream_bp"], background["downstream_bp"]
    length = upstream + downstream + 1
    require(type(upstream) is int and type(downstream) is int and upstream >= 0 and downstream >= 0
            and length > 1, "invalid promoter span")
    thresholds = comparison["thresholds"]
    minimum_bp, minimum_identity, maximum_evalue = (
        thresholds[key] for key in ("min_alignment_bp", "min_identity_pct", "max_evalue"))
    require(minimum_bp > 0 and 0 <= minimum_identity <= 100 and maximum_evalue >= 0,
            "invalid comparison thresholds")
    for key, region in regions.items():
        extraction = region["genome_extraction"]
        require(region["sequence_length_bp"] == summaries[key]["query_length_bp"] == length
                and extraction["genome_id"] == background["genome_id"]
                and extraction["promoter_upstream_bp"] == upstream
                and extraction["promoter_downstream_bp"] == downstream, "candidate geometry/reference mismatch")
        require(summaries[key].get("other_gene_exclusion_verified") is True, "query gene identity is unresolved")
        require(all(type(summaries[key].get(field)) is bool for field in
                    ("target_cap_reached", "hsp_cap_reached", "counts_are_lower_bounds")), "missing search-limit audit")
        require(summaries[key]["counts_are_lower_bounds"] == (
            summaries[key]["target_cap_reached"] or summaries[key]["hsp_cap_reached"]), "inconsistent search-limit audit")

    top: dict[str, list[dict[str, str]]] = {query_id: [] for query_id in regions}
    with paths["matches"].open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        match_columns = reader.fieldnames
        require(match_columns and "same_gene" in match_columns, "matches lack shared gene-exclusion policy")
        for row in reader:
            require(row["query_id"] in regions, "unknown match query")
            require(row["same_gene"] in {"True", "False"} and row["same_locus_overlap"] in {"True", "False"},
                    "invalid match exclusion flags")
            if row["same_locus_overlap"] == "True" or row["same_gene"] == "True":
                continue
            top[row["query_id"]].append(row)
    for query_id in top:
        top[query_id] = sorted(
            top[query_id],
            key=lambda row: (-float(row["aligned_query_fraction"]), -float(row["best_bitscore"]), row["promoter_id"]),
        )[:top_n]

    selected_pairs = {
        (query_id, row["promoter_id"])
        for query_id, rows in top.items()
        for row in rows
    }
    hsps: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    members = defaultdict(list)
    for member, representative in representatives.items():
        members[representative].append(member)
    for row in comparison_tools.iter_blast_rows(paths["raw_hits"]):
        require(row["qseqid"] in members, "unknown HSP representative")
        if comparison_tools.hsp_passes(row, minimum_bp, minimum_identity, maximum_evalue):
            for member in members[row["qseqid"]]:
                pair = (member, row["sseqid"])
                if pair in selected_pairs:
                    require(1 <= min(row["qstart"], row["qend"]) <= max(row["qstart"], row["qend"]) <= length,
                            "HSP query coordinate exceeds candidate")
                    hsps[pair].append(row)

    ordered_queries = sorted(
        regions,
        key=lambda query_id: (
            regions[query_id]["gene_query"],
            regions[query_id]["genome_extraction"]["tss_1based"],
        ),
    )
    height = 230 + PANEL_H * len(ordered_queries) + 100
    parameters = (f'{background["genome_id"]} · transcript-oriented -{upstream}/+{downstream} bp · '
                  f'BLASTN ≥{minimum_bp} bp, ≥{minimum_identity:g}% identity, E ≤ {maximum_evalue:g}')
    limited = any(row["counts_are_lower_bounds"] for row in summaries.values())
    limit_note = ("LOWER BOUNDS: a target or HSP cap was reached; absence and coverage are not exhaustive."
                  if limited else "No observed target/HSP cap saturation; absence remains conditional on this search and reference.")
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="{height}" viewBox="0 0 1600 {height}">',
        '<rect width="1600" height="100%" fill="white"/>',
        '<style>text{font-family:DejaVu Sans,Arial,sans-serif;fill:#202124}.title{font-size:25px;font-weight:700}.sub{font-size:14px}.head{font-size:15px;font-weight:700}.small{font-size:11px}.tiny{font-size:9px}</style>',
        '<text x="50" y="40" class="title">TP73 CUT&amp;RUN-supported TSS windows vs the human promoterome</text>',
        f'<text x="50" y="67" class="sub">{escape(parameters)}</text>',
        '<text x="50" y="91" class="sub">Support: at least one TP73 experimental-window mean exceeds its matched GFP-control mean in the same cell line.</text>',
        '<text x="50" y="115" class="sub">Frequency strip: distinct other genes per query position (log scale). Rows: top other-gene promoter windows by aggregate query coverage.</text>',
        f'<rect x="50" y="140" width="18" height="12" fill="#d2e6fa"/><text x="75" y="151" class="small">{minimum_identity:g}% identity</text>',
        '<rect x="180" y="140" width="18" height="12" fill="#2878d7"/><text x="205" y="151" class="small">100% identity</text>',
        '<rect x="330" y="138" width="24" height="16" fill="none" stroke="#d62728" stroke-width="2"/><text x="362" y="151" class="small">order/orientation break</text>',
        f'<line x1="{PLOT_X}" y1="174" x2="{PLOT_X + PLOT_W}" y2="174" stroke="#555"/>',
    ]
    for position, label in [(0, f"-{upstream}"), (upstream // 2, f"-{upstream - upstream // 2}"),
                            (upstream, "TSS"), (length - 1, f"+{downstream}")]:
        x = PLOT_X + PLOT_W * position / length
        svg.append(f'<line x1="{x:.1f}" y1="168" x2="{x:.1f}" y2="180" stroke="#555"/>')
        svg.append(f'<text x="{x:.1f}" y="164" text-anchor="middle" class="small">{label}</text>')

    summary_rows = []
    top_rows = []
    observations = []
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
        qualifier = "LOWER BOUNDS: " if summary["counts_are_lower_bounds"] else "Observed: "
        svg.append(f'<text x="50" y="{y + 34}" class="small">{qualifier}{counts["distinct_genes"]:,} genes · ≥25%: {coverage["0.25"]["distinct_genes"]:,} · ≥50%: {coverage["0.50"]["distinct_genes"]:,} · ≥80%: {coverage["0.80"]["distinct_genes"]:,}</text>')
        strip_y = y + 9
        svg.append(f'<rect x="{PLOT_X}" y="{strip_y}" width="{PLOT_W}" height="12" fill="#f5f5f5" stroke="#bbb"/>')
        for segment in summary["recurrent_query_segments"]:
            require(0 <= segment["query_start_0based"] < segment["query_end_0based_exclusive"] <= length,
                    "frequency segment exceeds candidate")
            x = PLOT_X + PLOT_W * segment["query_start_0based"] / length
            width = max(0.5, PLOT_W * (segment["query_end_0based_exclusive"] - segment["query_start_0based"]) / length)
            colour = frequency_blue(segment["distinct_genes"], max_frequency)
            svg.append(f'<rect x="{x:.2f}" y="{strip_y}" width="{width:.2f}" height="12" fill="{colour}"/>')
        svg.append(f'<line x1="{PLOT_X + PLOT_W * upstream / length:.1f}" y1="{strip_y - 2}" x2="{PLOT_X + PLOT_W * upstream / length:.1f}" y2="{strip_y + 14}" stroke="#111"/>')

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
            require(blocks, "ranked match lacks qualifying bound HSPs")
            for block_number, (hsp, broken) in enumerate(ordered_blocks(blocks), 1):
                q0 = min(int(hsp["qstart"]), int(hsp["qend"])) - 1
                q1 = max(int(hsp["qstart"]), int(hsp["qend"]))
                x = PLOT_X + PLOT_W * q0 / length
                width = max(2.0, PLOT_W * (q1 - q0) / length)
                stroke = "#d62728" if broken else "#4a6572"
                stroke_width = 2 if broken else 0.6
                svg.append(f'<rect x="{x:.2f}" y="{row_y + 1}" width="{width:.2f}" height="11" fill="{blue(float(hsp["pident"]), minimum_identity)}" stroke="{stroke}" stroke-width="{stroke_width}"/>')
                if width >= 12:
                    svg.append(f'<text x="{x + width / 2:.2f}" y="{row_y + 10}" text-anchor="middle" class="tiny">{block_number}</text>')
            top_rows.append({
                "query_id": query_id,
                "rank": rank,
                **target,
                "counts_are_lower_bounds": summary["counts_are_lower_bounds"],
            })
            row_y += ROW_H
        if not top[query_id]:
            svg.append(f'<text x="{PLOT_X}" y="{row_y + 11}" class="small">No qualifying other-gene hit was reported; this is not proof of uniqueness.</text>')
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
            "hsp_cap_reached": summary["hsp_cap_reached"],
            "counts_are_lower_bounds": summary["counts_are_lower_bounds"],
        })
        observation = (f'- {gene}, TSS {tss:,} ({", ".join(transcripts)}): {qualifier.lower()}'
                       f'{counts["distinct_promoter_windows"]:,} other-gene windows, '
                       f'{counts["distinct_genes"]:,} genes and {counts["distinct_transcripts"]:,} transcripts. '
                       f'Genes at >=25% / >=50% / >=80% aggregate coverage: '
                       f'{coverage["0.25"]["distinct_genes"]:,} / {coverage["0.50"]["distinct_genes"]:,} / '
                       f'{coverage["0.80"]["distinct_genes"]:,}.')
        if top[query_id]:
            best = top[query_id][0]
            observation += (f' Highest ranked observed window: {best["promoter_id"]} '
                            f'({best["gene_names"] or best["gene_ids"]}), '
                            f'{float(best["aligned_query_fraction"]):.1%} aggregate coverage.')
        else:
            observation += ' No qualifying other-gene hit was reported; this does not establish uniqueness.'
        observations.append(observation)
        y += PANEL_H

    svg.extend([
        f'<text x="50" y="{height - 82}" class="small">{escape(limit_note)}</text>',
        f'<text x="50" y="{height - 62}" class="small">Numbers are block order in the target promoter’s 5′→3′ sequence. Red outlines begin a changed-order or reversed-orientation block.</text>',
        f'<text x="50" y="{height - 42}" class="small">Window-mean differences and similarity do not prove peaks, direct binding, rearrangements, activity or sufficiency.</text>',
        f'<text x="50" y="{height - 22}" class="tiny">GENtle {escape(candidates["source_revision"])} · candidate SHA-256 {sha256(ROOT / "candidate_regions.json")} · comparison SHA-256 {sha256(RUN / "comparison.json")}</text>',
        '</svg>',
    ])
    output.mkdir(parents=True, exist_ok=True)
    (output / "first_comparison.svg").write_text("\n".join(svg) + "\n", encoding="utf-8")

    with (output / "summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    with (output / "top_matches.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rank", *match_columns, "counts_are_lower_bounds"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(top_rows)

    report = [
        "# First TP73 CUT&RUN-supported promoterome comparison",
        "",
        f'This run contains {len(regions)} selected TSS windows and a background of '
        f'{background["unique_promoter_window_count"]:,} transcript-linked promoter windows.',
        "",
        "## Parameters",
        "",
        f"- {parameters}",
        "- Support: at least one TAp73α or DNp73β mean exceeds the matched GFP-control mean in the same cell line across that exact window.",
        f'- Search: dust={comparison["search_policy"]["dust"]}; soft masking={comparison["search_policy"]["soft_masking"]}; '
        f'target cap={thresholds["max_target_seqs"]}; HSP cap per target={thresholds["max_hsps_per_target"]}.',
        "- Counting: overlapping windows and windows sharing any query gene ID are excluded consistently from counts and detail rows.",
        f"- {limit_note}",
        "",
        "## First observations",
        "",
        *observations,
        "",
        "These are sequence-recurrence observations, not reporter recommendations. The next decision layer should intersect the recurrent blocks with localized CUT&RUN peaks, Ensembl Regulation features, motif/module evidence, and planned deletion/split constructs.",
        "Positive window-mean differences do not establish peaks, significance or activity. Block-order markers are descriptive, not proof of rearrangement.",
    ]
    (output / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    png_provenance = render_png(output, png_renderer) if png_renderer else {"status": "not_requested"}
    output_names = ["first_comparison.svg", "summary.tsv", "top_matches.tsv", "report.md"]
    if png_renderer:
        output_names.append("first_comparison.png")
    receipt = {
        "schema": "gentle.tp73_cutrun_supported_promoter_comparison_receipt.v1",
        "source_revision": candidates["source_revision"],
        "renderer_sha256": sha256(Path(__file__)),
        "comparison_helpers_sha256": sha256(Path(comparison_tools.__file__)),
        "inputs": inputs,
        "outputs": {name: sha256(output / name) for name in output_names},
        "png_rendering": png_provenance,
        "render_parameters": {"top": top_n, "task": "blastn", "width_px": 1600},
        "counting_policy_id": comparison_tools.COUNTING_POLICY,
        "target_cap_reached": any(row["target_cap_reached"] for row in summary_rows),
        "hsp_cap_reached": any(row["hsp_cap_reached"] for row in summary_rows),
        "counts_are_lower_bounds": limited,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True,
                        help="Directory containing candidates and a corrected comparison task directory")
    parser.add_argument("--output", type=Path, required=True, help="New/empty directory; never overwrites retained evidence")
    parser.add_argument("--task-directory", default="blastn-40bp-80pct")
    parser.add_argument("--top", type=int, default=5)
    parser.add_argument("--png-renderer", help="Optional rsvg-convert executable; its version, digest and invocation are recorded")
    args = parser.parse_args()
    receipt = render(args.root, args.output.resolve(), args.task_directory, args.top, args.png_renderer)
    print(json.dumps({"output": str(args.output), "counts_are_lower_bounds": receipt["counts_are_lower_bounds"]}, indent=2))


if __name__ == "__main__":
    main()
