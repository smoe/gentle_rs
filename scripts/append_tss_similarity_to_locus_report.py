#!/usr/bin/env python3
"""Append TSS-local regulatory-feature promoter similarity to a GENtle locus SVG."""

from __future__ import annotations

import argparse
from html import escape
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any
import xml.etree.ElementTree as ET

try:
    from .render_integrated_tss_regulatory_report import (
        comparison_tools, frequency_segments, load_bound_comparison,
    )
    from .render_tp73_cutrun_promoter_comparison import ordered_blocks
    from .tss_regulatory_report_binding import load_bound_locus_report, load_bound_locus_svg
except ImportError:
    from render_integrated_tss_regulatory_report import (
        comparison_tools, frequency_segments, load_bound_comparison,
    )
    from render_tp73_cutrun_promoter_comparison import ordered_blocks
    from tss_regulatory_report_binding import load_bound_locus_report, load_bound_locus_svg


FEATURE_COLOURS = {
    "promoter": "#5b8ff9",
    "enhancer": "#61d9a8",
    "emar": "#c7a0e8",
    "ctcf": "#f6bd5b",
    "open_chromatin_region": "#78d3f8",
}
STRETCH_COLOURS = ("#0f766e", "#b45309", "#1d4ed8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def render_derivatives(svg: Path, pdf: Path, png: Path, renderer: str) -> dict[str, Any]:
    executable = shutil.which(renderer)
    if executable is None:
        raise RuntimeError("SVG renderer unavailable; provide --renderer rsvg-convert path")
    executable = str(Path(executable).resolve())
    version = subprocess.run(
        [executable, "--version"], check=True, capture_output=True, timeout=30,
    )
    commands = []
    for output_format, output in (("pdf", pdf), ("png", png)):
        command = [executable, f"--format={output_format}", f"--output={output}", str(svg)]
        result = subprocess.run(command, check=True, capture_output=True, timeout=180)
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(f"renderer did not produce {output_format.upper()}")
        commands.append({
            "command": command,
            "stdout_sha256": hashlib.sha256(result.stdout).hexdigest(),
            "stderr_sha256": hashlib.sha256(result.stderr).hexdigest(),
        })
    with png.open("rb") as handle:
        if handle.read(8) != b"\x89PNG\r\n\x1a\n":
            raise RuntimeError("renderer output is not a PNG")
    return {
        "executable": executable,
        "executable_sha256": f"sha256:{sha256(Path(executable))}",
        "version": version.stdout.decode(errors="replace").strip(),
        "commands": commands,
    }


def blue(identity: float, minimum_identity: float) -> str:
    fraction = min(1.0, max(0.0, (identity - minimum_identity)
                            / max(1.0, 100.0 - minimum_identity)))
    return "#{:02x}{:02x}{:02x}".format(
        round(210 - 170 * fraction), round(230 - 120 * fraction), round(250 - 35 * fraction)
    )
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


def query_interval_pixels(start: int, end: int, length: int,
                          left: float, width: float, strand: str) -> tuple[float, float]:
    """Project assembly-forward, half-open query offsets onto the gene-oriented axis."""
    if strand not in {"+", "-"} or not 0 <= start < end <= length:
        raise ValueError("invalid query interval or gene strand")
    offset = start if strand == "+" else length - end
    return left + offset / length * width, (end - start) / length * width


def require_new_outputs(paths: list[Path]) -> None:
    resolved = [path.resolve() for path in paths]
    if len(set(resolved)) != len(resolved) or any(path.exists() or path.is_symlink() for path in paths):
        raise ValueError("output paths must be distinct and absent; do not overwrite bound evidence")


def add_stretch_overview(svg: str, gene: str, candidates: dict[str, Any],
                         report: dict[str, Any]) -> tuple[str, int]:
    """Insert linked stretch labels and background bands on the bound genomic axis."""
    root = ET.fromstring(svg)
    lines = [item for item in root.iter() if item.get("data-gentle-transcript")]
    heading = re.search(r'<text\b[^>]*>\s*Transcript models and annotation-derived metrics\s*</text>', svg)
    stretches = [row for row in candidates["stretches"] if row["gene"] == gene]
    if not heading or not lines or not stretches:
        raise ValueError("base SVG lacks a transcript-model axis for TSS-stretch references")
    frames = {(float(item.get("x1")), float(item.get("x2"))) for item in lines}
    if len(frames) != 1:
        raise ValueError("base SVG has inconsistent transcript-model axes")
    x0, x1 = next(iter(frames))
    left, right = report["axis_left_genomic_1based"], report["axis_right_genomic_1based"]
    if left == right or x0 >= x1:
        raise ValueError("invalid genomic display axis")
    y = float(ET.fromstring(heading.group()).get("y"))
    start_y = y
    similarity = next((item for item in root.iter()
                       if item.get("data-gentle-panel") == "tss-local-promoter-similarity"), None)
    similarity_heading = next((item for item in ([] if similarity is None else similarity)
                               if item.tag.rsplit("}", 1)[-1] == "text"), None)
    if similarity_heading is None:
        raise ValueError("base SVG lacks the similarity boundary for TSS background bands")
    band_top = start_y - 14
    band_bottom = float(similarity_heading.get("y")) - 20
    if not all(math.isfinite(value) for value in (band_top, band_bottom)) or band_bottom <= band_top:
        raise ValueError("invalid upper genomic overview height for TSS background bands")
    backgrounds = ['<g data-gentle-tss-stretch-backgrounds="true" pointer-events="none">']
    body = ['<g data-gentle-tss-stretch-overview="true">',
            svg_text(34, y, "TSS stretches: reference for similarity below", size=12, weight="bold")]
    y += 25
    for index, stretch in enumerate(stretches):
        start, end = stretch["start_1based"], stretch["end_1based"]
        if not min(left, right) <= start <= end <= max(left, right):
            raise ValueError("TSS stretch is outside the original genomic display axis")
        a, b = sorted(x0 + (position - left) / (right - left) * (x1 - x0)
                      for position in (start, end))
        name = escape(stretch["stretch_id"], quote=True)
        colour = STRETCH_COLOURS[index % len(STRETCH_COLOURS)]
        backgrounds.append(
            f'<rect data-gentle-tss-stretch-band="{name}" data-gentle-genomic-start="{start}" '
            f'data-gentle-genomic-end="{end}" x="{a:.2f}" y="{band_top:.2f}" '
            f'width="{max(1, b-a):.2f}" height="{band_bottom-band_top:.2f}" '
            f'fill="{colour}" fill-opacity="0.12"/>'
        )
        body.extend([
            f'<a id="overview-{name}" href="#similarity-{name}">',
            svg_text(34, y + 4, stretch["stretch_id"], size=9, family="monospace", fill=colour),
            f'<line x1="{x0:.2f}" x2="{x1:.2f}" y1="{y}" y2="{y}" stroke="#cbd5e1"/>',
            f'<rect data-gentle-tss-stretch="{name}" data-gentle-genomic-start="{start}" '
            f'data-gentle-genomic-end="{end}" x="{a:.2f}" y="{y-6:.2f}" '
            f'width="{max(1, b-a):.2f}" height="12" fill="{colour}"/>',
            svg_text(x1 + 20, y + 4, f"{start:,}..{end:,}", size=8, family="monospace"),
            '</a>',
        ])
        y += 26
    body.append('</g>')
    backgrounds.append('</g>')
    shift = math.ceil(y - start_y + 18)
    old_height = int(root.get("height"))
    height = old_height + shift
    prefix = svg[:heading.start()]
    prefix = prefix.replace(f'height="{old_height}"', f'height="{height}"')
    prefix = prefix.replace(f'viewBox="0 0 1400 {old_height}"', f'viewBox="0 0 1400 {height}"')
    tail = svg[heading.start():].rsplit('</svg>', 1)[0]
    # Paint bands first, behind the untouched scientific lanes, never over their glyphs.
    return (prefix + '\n'.join(body) + f'<g data-gentle-shifted-locus="true" transform="translate(0 {shift})">'
            + '\n'.join(backgrounds) + tail + '</g></svg>\n', height)


def append_section(
    base_svg: str,
    gene: str,
    candidates: dict[str, Any],
    matches_by_query: dict[str, list[dict[str, str]]],
    hsps: dict[tuple[str, str], list[dict[str, Any]]],
    summaries: dict[str, dict[str, Any]],
    thresholds: dict[str, float],
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
    boundary = re.search(
        r'<text[^>]*data-gentle-overlay-non-claims="true"[^>]*\by="([0-9.]+)"[^>]*>\s*\n?'
        r'Reporter interpretation boundaries', base_svg,
    )
    if not boundary:
        raise ValueError("base SVG lacks the interpretation/provenance insertion boundary")
    boundary_y = float(boundary.group(1))
    boundary_start = boundary.start()
    body = [f'<g data-gentle-panel="tss-local-promoter-similarity" data-gentle-gene="{escape(gene)}">']
    y = boundary_y
    body.append(svg_text(34, y, "TSS-local Ensembl-feature promoterome similarity",
                         size=14, fill="#1f2937", weight="bold"))
    y += 20
    body.append(svg_text(
        34, y,
        "Transcript-oriented −500/+200 bp; features clipped to the displayed stretch; "
        f"BLASTN ≥{thresholds['minimum_bp']:g} bp, ≥{thresholds['minimum_identity']:g}% identity, "
        f"E≤{thresholds['maximum_evalue']:g}; self-locus and same-gene targets excluded.",
        size=9, family="monospace", fill="#64748b",
    ))
    y += 18
    body.append(svg_text(
        34, y,
        "Blue frequency strips count distinct other genes by query position. Numbered blocks give target-promoter order; red outlines mark order/orientation breaks.",
        size=9, fill="#64748b",
    ))
    y += 26
    for index, stretch in enumerate(stretches):
        strand = stretch["tss_windows"][0]["strand"]
        regions = regions_by_stretch[stretch["stretch_id"]]
        name = escape(stretch["stretch_id"], quote=True)
        body.append(f'<a id="similarity-{name}" href="#overview-{name}">')
        body.append(svg_text(34, y, stretch["stretch_id"], size=11,
                             fill=STRETCH_COLOURS[index % len(STRETCH_COLOURS)], weight="bold"))
        body.append('</a>')
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
            rows = matches_by_query.get(query_id, [])
            summary = summaries[query_id]
            counts = summary["other_promoters"]
            tiers = summary["other_promoters_by_min_query_coverage"]
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
            if length < thresholds["minimum_bp"]:
                body.append(svg_text(1070, y + 8,
                                     f"motif-scale; below {thresholds['minimum_bp']:g}-bp search gate",
                                     size=8, family="monospace", fill="#64748b"))
                y += 54
                continue
            status = ("LOWER BOUNDS" if summary["counts_are_lower_bounds"]
                      else "observed; no cap saturation")
            body.append(svg_text(
                1070, y + 8,
                f"other genes {counts['distinct_genes']:,}; "
                f"≥25% {tiers['0.25']['distinct_genes']:,}; "
                f"≥50% {tiers['0.50']['distinct_genes']:,}; {status}",
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
                sx, sw = query_interval_pixels(start, end, length, fx0, max(2, fw), strand)
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
                blocks = hsps.get((query_id, target["promoter_id"]), [])
                if not blocks:
                    raise ValueError("ranked match lacks qualifying bound HSPs")
                for number, (hsp, broken) in enumerate(ordered_blocks(blocks), 1):
                    q0 = min(int(hsp["qstart"]), int(hsp["qend"])) - 1
                    q1 = max(int(hsp["qstart"]), int(hsp["qend"]))
                    bx, bw = query_interval_pixels(q0, q1, length, fx0, max(2, fw), strand)
                    body.append(f'<rect x="{bx:.2f}" y="{y-5}" width="{max(1,bw):.2f}" height="12" fill="{blue(float(hsp["pident"]), thresholds["minimum_identity"])}" stroke="{"#dc2626" if broken else "#475569"}" stroke-width="{1.4 if broken else 0.5}"/>')
                    if bw >= 10:
                        body.append(svg_text(bx + bw / 2, y + 4, str(number), size=7,
                                             fill="#111827", anchor="middle"))
                y += 19
            if not top:
                body.append(svg_text(652, y + 4, "No qualifying other-gene promoter match",
                                     size=8, fill="#64748b", anchor="middle"))
                y += 18
            y += 16
        y += 18
    body.append(svg_text(
        34, y,
        "Similarity is structural evidence; Ensembl classes, predicted TFBS and CUT&RUN enrichment do not prove reporter activity or sufficiency.",
        size=9, fill="#64748b",
    ))
    body.append("</g>")
    # Reserve the height actually drawn, including all ranked matches and empty rows.
    section_height = math.ceil(y - boundary_y + 32)
    new_height = old_height + section_height
    prefix = base_svg[:boundary_start]
    opening = root.group(0).replace(f'height="{old_height}"', f'height="{new_height}"')
    opening = opening.replace(f'viewBox="0 0 1400 {old_height}"', f'viewBox="0 0 1400 {new_height}"')
    prefix = prefix.replace(root.group(0), opening, 1)
    prefix = prefix.replace(
        f'<rect fill="#ffffff" height="{old_height}" width="1400" x="0" y="0"/>',
        f'<rect fill="#ffffff" height="{new_height}" width="1400" x="0" y="0"/>', 1,
    )
    tail = base_svg[boundary_start:].rsplit("</svg>", 1)[0]
    translated_tail = f'<g data-gentle-shifted-footer="true" transform="translate(0 {section_height})">\n{tail}\n</g>'
    return prefix + "\n" + "\n".join(body) + "\n" + translated_tail + "\n</svg>\n", new_height


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-svg", type=Path, required=True)
    parser.add_argument("--locus-report", type=Path, required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--candidates-json", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--matches", type=Path, required=True)
    parser.add_argument("--hits", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--output-png", type=Path, required=True)
    parser.add_argument("--renderer", default="rsvg-convert")
    args = parser.parse_args()
    receipt_path = args.output_svg.with_suffix(".receipt.json")
    require_new_outputs([args.output_svg, args.output_pdf, args.output_png, receipt_path])

    candidates, _, regions, matches_by_query, hsps, summaries, thresholds = \
        load_bound_comparison(
            args.candidates_json, args.comparison, args.matches, args.hits)
    require_gene = {row["gene_query"] for row in regions.values() if row["gene_query"] == args.gene}
    if require_gene != {args.gene}:
        raise ValueError("requested gene has no candidate regions")
    report, report_digest = load_bound_locus_report(candidates, args.locus_report, args.gene)
    base_svg, svg_digest = load_bound_locus_svg(candidates, args.base_svg, report)
    matches_by_query = {
        query_id: rows for query_id, rows in matches_by_query.items()
        if regions[query_id]["gene_query"] == args.gene
    }
    output, height = append_section(
        base_svg, args.gene, candidates, matches_by_query, hsps,
        summaries, thresholds
    )
    output, height = add_stretch_overview(output, args.gene, candidates, report)
    for path in (args.output_svg, args.output_pdf, args.output_png):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(output, encoding="utf-8")
    renderer = render_derivatives(
        args.output_svg.resolve(), args.output_pdf.resolve(), args.output_png.resolve(),
        args.renderer,
    )
    renderer_revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[1], text=True
    ).strip()
    receipt = {
        "schema": "gentle.tss_local_similarity_locus_report_receipt.v1",
        "gene": args.gene,
        "source_revision": candidates["source_revision"],
        "candidate_source_revision": candidates["source_revision"],
        "renderer_revision": renderer_revision,
        "producer_sha256": f"sha256:{sha256(Path(__file__))}",
        "counting_policy_id": comparison_tools.COUNTING_POLICY,
        "inputs": {
            "base_svg": f"sha256:{svg_digest}",
            "locus_report": f"sha256:{report_digest}",
            "candidates": f"sha256:{sha256(args.candidates_json)}",
            "comparison": f"sha256:{sha256(args.comparison)}",
            "matches": f"sha256:{sha256(args.matches)}",
            "hits": f"sha256:{sha256(args.hits)}",
        },
        "outputs": {
            args.output_svg.name: f"sha256:{sha256(args.output_svg)}",
            args.output_pdf.name: f"sha256:{sha256(args.output_pdf)}",
            args.output_png.name: f"sha256:{sha256(args.output_png)}",
        },
        "renderer": renderer,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    print(json.dumps({
        "output_svg": str(args.output_svg),
        "height": height,
        "base_svg_sha256": f"sha256:{svg_digest}",
        "candidate_sha256": f"sha256:{sha256(args.candidates_json)}",
        "receipt": str(receipt_path),
    }, indent=2))


if __name__ == "__main__":
    main()
