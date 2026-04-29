#!/usr/bin/env python3
"""Render the README TP73 cDNA isoform-selector PCR matrix from GENtle reports."""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path


FORWARD_SELECTORS = [
    ("TA", "TA / full N-terminus", "ACGGCTGCAGAGCGAGCTGC"),
    ("X2", "X2 alternative 5' cap", "AGGCTTAGCCAAGAGCGAGCTGC"),
    ("X1", "X1 alternative 5' cap", "AAATAACAGAGCGAGCTGCCCTCG"),
    ("DN", "DeltaN / internal promoter", "ACCTCGCCACGGCCCAGTTCAAT"),
    ("V13", "variant-13 5' exon", "CATCTGCAGGACAGGCCCAGTTC"),
]

REVERSE_READOUTS = [
    (
        "R-alpha-beta",
        "alpha/beta junction",
        "GAGCATCCCGGGGCCCAC",
        "77485..77633 -> 78412..78550",
    ),
    (
        "R-beta-epsilon",
        "beta/epsilon junction",
        "CCAGGTCCTGACGAGGCTGG",
        "78412..78550 -> 80232..83686",
    ),
    (
        "R-gamma-epsilon",
        "gamma/epsilon junction",
        "ATCCCGGGGCGGCCTCTGTAG",
        "76812..76933 -> 78412..78550",
    ),
    (
        "R-delta",
        "delta junction",
        "CCAGGTCGGCCTCTGTAGGAGCT",
        "76812..76933 -> 80232..83686",
    ),
    (
        "R-zeta",
        "zeta junction",
        "TCCTGTTAAAAAAGGCCTCTGTA",
        "76812..76933 -> 78948..79041",
    ),
    (
        "R-alpha-gamma-zeta",
        "alpha/gamma/zeta junction",
        "CCCCAGGTCCTCAATGGTCA",
        "78948..79041 -> 80232..83686",
    ),
    (
        "R-common",
        "common terminal exon",
        "AGTCGTGGCCCTGCTTCAGGTCCT",
        "common terminal exon",
    ),
]

ISOFORM_LABELS = {
    "NM_005427.4": "TAp73α",
    "NM_001204184.2": "TAp73β",
    "NM_001204185.2": "TAp73γ",
    "NM_001204186.2": "TAp73δ",
    "NM_001204187.2": "TAp73ε",
    "NM_001204188.2": "TAp73ζ",
    "NM_001126240.3": "ΔNp73α",
    "NM_001126241.3": "ΔNp73β",
    "NM_001126242.3": "ΔNp73γ",
    "NM_001204189.2": "ΔNp73δ",
    "NM_001204190.2": "ΔNp73ε",
    "NM_001204191.2": "ΔNp73ζ",
    "XM_047429524.1": "X2-p73α",
    "XM_047429521.1": "X1-p73α",
    "NM_001204192.2": "V13-p73α",
}


def svg_escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def wrap_text(text: str, max_chars: int = 24, max_lines: int = 3) -> list[str]:
    lines: list[str] = []
    current = ""
    for word in text.split(" "):
        candidate = f"{current} {word}".strip()
        if current and len(candidate) > max_chars:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines[:max_lines]


def report_cell(report_dir: Path, forward_key: str, reverse_key: str) -> dict[str, object]:
    path = report_dir / f"{forward_key}_{reverse_key}.json"
    report = json.loads(path.read_text())
    hits: list[tuple[str, str, str]] = []
    for row in report["transcript_results"]:
        products = row.get("products", [])
        if not products:
            continue
        lengths = "/".join(
            str(length)
            for length in sorted({product["amplicon_length_bp"] for product in products})
        )
        transcript_id = row["transcript_id"]
        hits.append((transcript_id, ISOFORM_LABELS.get(transcript_id, transcript_id), lengths))
    return {
        "status": report["overall_status"],
        "count": report["product_count"],
        "lengths": ",".join(map(str, report["construct_lengths"]["detected_amplicon_lengths_bp"]))
        or "none",
        "hits": hits,
        "qc": report["oligo_qc"]["status"],
        "gdna": report["genomic_carryover_risk"]["risk_level"],
    }


def render_matrix(report_dir: Path) -> str:
    cell_width = 185
    cell_height = 98
    left_width = 210
    top_height = 175
    width = left_width + cell_width * len(REVERSE_READOUTS) + 50
    height = top_height + cell_height * len(FORWARD_SELECTORS) + 235
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">',
        "<title>TP73 cDNA isoform-selector PCR matrix</title>",
        "<style>"
        "text{font-family:Avenir Next,Segoe UI,Arial,sans-serif}"
        ".title{font-size:25px;font-weight:800;fill:#0f172a}"
        ".sub{font-size:13px;fill:#475569}"
        ".hdr{font-size:12px;font-weight:800;fill:#0f172a}"
        ".small{font-size:10px;fill:#475569}"
        ".seq{font-size:9px;fill:#334155;font-family:Menlo,Consolas,monospace}"
        ".celltext{font-size:11px;font-weight:700;fill:#0f172a}"
        ".len{font-size:10px;fill:#334155}"
        ".none{fill:#f1f5f9;stroke:#cbd5e1}"
        ".single{fill:#dcfce7;stroke:#16a34a}"
        ".multi{fill:#fef3c7;stroke:#d97706}"
        ".fail{stroke:#dc2626;stroke-width:2}"
        ".warn{stroke:#f97316;stroke-width:1.4}"
        ".pass{stroke-width:1}"
        "</style>",
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        f'<rect x="18" y="18" width="{width - 36}" height="{height - 36}" '
        'rx="12" fill="#f8fafc" stroke="#cbd5e1"/>',
        '<text class="title" x="36" y="55">TP73 cDNA isoform-selector PCR matrix</text>',
        '<text class="sub" x="36" y="78">Five 5-prime selector primers crossed '
        "with seven 3-prime splice/readout primers; products were tested with "
        "GENtle primers test-cdna-pcr on the bundled TP73 RefSeq locus.</text>",
        '<text class="sub" x="36" y="99">Readout logic: α=α/β + α/γ/ζ, '
        "β=α/β + β/ε, γ=γ/ε + α/γ/ζ, δ=δ, ε=γ/ε + β/ε, ζ=ζ + α/γ/ζ; "
        "η has no separate local RefSeq cDNA row here.</text>",
    ]

    for column_idx, (_key, label, seq, junction) in enumerate(REVERSE_READOUTS):
        x = left_width + column_idx * cell_width + 8
        svg.extend(
            [
                f'<text class="hdr" x="{x}" y="128">{svg_escape(label)}</text>',
                f'<text class="small" x="{x}" y="145">{svg_escape(junction)}</text>',
                f'<text class="seq" x="{x}" y="160">{svg_escape(seq[:28])}</text>',
            ]
        )
    svg.append('<text class="hdr" x="36" y="145">5-prime selector</text>')

    for row_idx, (forward_key, forward_label, forward_seq) in enumerate(FORWARD_SELECTORS):
        y = top_height + row_idx * cell_height
        svg.extend(
            [
                f'<text class="hdr" x="36" y="{y + 23}">{svg_escape(forward_key)}: '
                f"{svg_escape(forward_label)}</text>",
                f'<text class="seq" x="36" y="{y + 41}">{svg_escape(forward_seq[:28])}</text>',
            ]
        )
        for column_idx, (reverse_key, _label, _seq, _junction) in enumerate(REVERSE_READOUTS):
            x = left_width + column_idx * cell_width
            cell = report_cell(report_dir, forward_key, reverse_key)
            count = int(cell["count"])
            css_class = "none" if count == 0 else "single" if count == 1 else "multi"
            qc = str(cell["qc"])
            if qc == "fail":
                css_class += " fail"
            elif qc == "warning":
                css_class += " warn"
            else:
                css_class += " pass"
            svg.append(
                f'<rect class="{css_class}" x="{x + 5}" y="{y + 8}" '
                f'width="{cell_width - 10}" height="{cell_height - 14}" rx="8"/>'
            )
            if count == 0:
                svg.append(f'<text class="celltext" x="{x + 15}" y="{y + 36}">no product</text>')
            else:
                isoforms = ", ".join(sorted({hit[1] for hit in cell["hits"]}))
                text_y = y + 28
                for line in wrap_text(isoforms):
                    svg.append(
                        f'<text class="celltext" x="{x + 15}" y="{text_y}">'
                        f"{svg_escape(line)}</text>"
                    )
                    text_y += 14
                svg.append(
                    f'<text class="len" x="{x + 15}" y="{y + 80}">'
                    f'{svg_escape(cell["lengths"])} bp; QC {svg_escape(cell["qc"])}; '
                    f'gDNA {svg_escape(cell["gdna"])}</text>'
                )

    legend_y = top_height + cell_height * len(FORWARD_SELECTORS) + 34
    svg.extend(
        [
            f'<rect class="single pass" x="36" y="{legend_y - 12}" width="18" height="12" rx="3"/>'
            f'<text class="small" x="62" y="{legend_y - 2}">single product</text>',
            f'<rect class="multi pass" x="170" y="{legend_y - 12}" width="18" height="12" rx="3"/>'
            f'<text class="small" x="196" y="{legend_y - 2}">multiple products; interpret by band size and companion readouts</text>',
            f'<rect class="none pass" x="505" y="{legend_y - 12}" width="18" height="12" rx="3"/>'
            f'<text class="small" x="531" y="{legend_y - 2}">not detected</text>',
            f'<rect fill="none" stroke="#dc2626" stroke-width="2" x="650" y="{legend_y - 12}" width="18" height="12" rx="3"/>'
            f'<text class="small" x="676" y="{legend_y - 2}">red border: oligo QC fail; orange border: QC warning</text>',
            f'<text class="small" x="36" y="{legend_y + 24}">Full per-pair JSON/SVG/PNG maps are generated from the same reports. The SVG rows use cDNA templates derived from test_files/tp73.ncbi.gb feature 4 on TargetGroupTargetStrand.</text>',
            "</svg>",
        ]
    )
    return "\n".join(svg)


def write_markdown_index(report_dir: Path) -> None:
    index_path = report_dir / "index.md"
    with index_path.open("w") as out:
        out.write("# TP73 full cDNA isoform-selector PCR matrix\n\n")
        out.write(
            "Generated with `primers test-cdna-pcr` against the TP73 splicing "
            "group in `test_files/tp73.ncbi.gb`.\n\n"
        )
        out.write(
            "| Pair | Status | Products | Amplicons bp | Detected isoforms | "
            "Oligo QC | gDNA carryover | PNG | SVG |\n"
        )
        out.write("|---|---:|---:|---|---|---|---|---|---|\n")
        for forward_key, _forward_label, _forward_seq in FORWARD_SELECTORS:
            for reverse_key, _reverse_label, _reverse_seq, _junction in REVERSE_READOUTS:
                cell = report_cell(report_dir, forward_key, reverse_key)
                pair = f"{forward_key}_{reverse_key}"
                detected = ", ".join(f"{label}={lengths}" for _tid, label, lengths in cell["hits"])
                if not detected:
                    detected = "none"
                out.write(
                    f"| `{pair}` | `{cell['status']}` | {cell['count']} | "
                    f"{cell['lengths']} | {detected} | `{cell['qc']}` | `{cell['gdna']}` | "
                    f"[PNG]({report_dir / (pair + '.png')}) | "
                    f"[SVG]({report_dir / (pair + '.svg')}) |\n"
                )


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(
            "usage: tp73_isoform_selector_matrix.py REPORT_DIR OUTPUT.svg",
            file=sys.stderr,
        )
        return 2
    report_dir = Path(argv[1])
    output = Path(argv[2])
    output.write_text(render_matrix(report_dir))
    write_markdown_index(report_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
