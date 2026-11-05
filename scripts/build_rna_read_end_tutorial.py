#!/usr/bin/env python3
"""Regenerate the synthetic RNA read-end tutorial fixture and figure."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

import rank_rna_read_end_genes as ranking


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "docs/tutorial/reproducibility/rna_read_end_gene_ranking"
FIGURE = ROOT / "docs/figures/rna_read_end_gene_ranking.synthetic.svg"
GEOMETRY = ROOT / "docs/figures/rna_read_end_projection_geometry.svg"
GENES = ("MDM2", "MDM4", "FUS", "ACTB")


def series(center: int) -> list[int]:
    return [center - 4, center - 2, center, center + 2, center + 4]


def main() -> None:
    if FIXTURE.exists():
        shutil.rmtree(FIXTURE)
    FIXTURE.mkdir(parents=True)
    projections = []
    for sample_number in range(1, 7):
        sample = f"sample{sample_number}"
        path = FIXTURE / f"{sample}.projections.tsv"
        projections.append(path)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(ranking.HEADER)
            unit = 0
            for gene_index, gene in enumerate(GENES):
                read_center = 1000 + 30 * gene_index + 4 * sample_number
                five_center = 70 + 5 * gene_index + sample_number
                ordinary_three = 180 + 30 * gene_index + 3 * sample_number
                focus_shift = {"MDM2": 90, "MDM4": -65, "FUS": 35, "ACTB": 0}[gene]
                three_center = ordinary_three + (focus_shift if sample_number == 2 else 0)
                for metric, values in (
                    ("read_length_bp", series(read_center)),
                    ("utr5_reach_bp", series(five_center)),
                    ("utr3_reach_bp", series(three_center)),
                ):
                    for value in values:
                        unit += 1
                        writer.writerow([
                            sample, "sha256:" + f"{sample_number:x}" * 64,
                            f"synthetic-{sample}-{unit}", metric, value, "observed",
                            "resolved", f"synthetic-{gene}-tx", gene,
                        ])
    output = FIXTURE / "ranked"
    ranking.run(projections, "sample2", output, 4, 0.25, 12)
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(output / "top_genes.svg", FIGURE)
    GEOMETRY.write_text("""<svg xmlns="http://www.w3.org/2000/svg" width="1180" height="380" viewBox="0 0 1180 380">
<rect width="1180" height="380" fill="white"/>
<style>text{font-family:sans-serif}.h{font-size:21px;font-weight:700}.l{font-size:14px}.s{font-size:12px;fill:#444}</style>
<text x="28" y="35" class="h">What GENtle measures at transcript ends</text>
<text x="28" y="62" class="s">Coordinates are projected through exons into transcript orientation; the CDS boundaries are explicit anchors.</text>
<text x="28" y="112" class="l">frozen transcript model</text>
<rect x="220" y="92" width="150" height="28" rx="4" fill="#6bb7a8"/><text x="295" y="111" text-anchor="middle" class="l">5′ UTR</text>
<rect x="370" y="92" width="430" height="28" fill="#44546a"/><text x="585" y="111" text-anchor="middle" class="l" fill="white">CDS</text>
<rect x="800" y="92" width="260" height="28" rx="4" fill="#dc7783"/><text x="930" y="111" text-anchor="middle" class="l">3′ UTR</text>
<line x1="370" x2="370" y1="82" y2="236" stroke="#222" stroke-dasharray="4 3"/><line x1="800" x2="800" y1="82" y2="236" stroke="#222" stroke-dasharray="4 3"/>
<text x="370" y="76" text-anchor="middle" class="s">CDS start</text><text x="800" y="76" text-anchor="middle" class="s">CDS end</text>
<text x="28" y="174" class="l">one aligned read</text><line x1="286" x2="955" y1="166" y2="166" stroke="#263238" stroke-width="8" stroke-linecap="round"/>
<line x1="286" x2="370" y1="192" y2="192" stroke="#2a9d8f" stroke-width="4"/><path d="M286 192l10 -6v12zM370 192l-10 -6v12z" fill="#2a9d8f"/><text x="328" y="215" text-anchor="middle" class="s">5′ reach</text>
<line x1="800" x2="955" y1="192" y2="192" stroke="#d1495b" stroke-width="4"/><path d="M800 192l10 -6v12zM955 192l-10 -6v12z" fill="#d1495b"/><text x="878" y="215" text-anchor="middle" class="s">3′ reach</text>
<text x="28" y="275" class="l">clipped 3′ end</text><line x1="286" x2="900" y1="267" y2="267" stroke="#263238" stroke-width="8" stroke-linecap="round"/><line x1="900" x2="1005" y1="267" y2="267" stroke="#d1495b" stroke-width="4" stroke-dasharray="7 5"/><text x="952" y="251" text-anchor="middle" class="s">unknown continuation</text>
<rect x="781" y="292" width="260" height="42" rx="6" fill="#fff3cd" stroke="#b58b00"/><text x="911" y="309" text-anchor="middle" class="s">right-censored lower bound</text><text x="911" y="325" text-anchor="middle" class="s">never converted to an exact zero</text>
<text x="28" y="358" class="s">If the alignment never crosses a CDS boundary, that UTR metric is missing. Alternative transcripts must agree on value and status.</text>
</svg>
""", encoding="utf-8")


if __name__ == "__main__":
    main()
