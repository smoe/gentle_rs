#!/usr/bin/env python3
"""Derive a calmer README hero SVG from a technical serial gel export.

This keeps the rendered gel bands and lane geometry, but replaces the dense
header/axis/lane labels with a tighter presentation-oriented layout.
"""

from __future__ import annotations

import math
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def qname(tag: str) -> str:
    return f"{{{SVG_NS}}}{tag}"


def parse_rect(element: ET.Element) -> tuple[float, float, float, float]:
    return (
        float(element.get("x", "0")),
        float(element.get("y", "0")),
        float(element.get("width", "0")),
        float(element.get("height", "0")),
    )


def main() -> int:
    if len(sys.argv) < 4:
        print(
            "usage: render_serial_gel_hero.py INPUT.svg OUTPUT.svg LABEL [LABEL ...]",
            file=sys.stderr,
        )
        return 2

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    lane_labels = sys.argv[3:]

    tree = ET.parse(input_path)
    root = tree.getroot()

    background = None
    gel_panel = None
    lane_rects: list[ET.Element] = []
    band_rects: list[ET.Element] = []
    scale_labels: list[tuple[int, float]] = []

    for child in list(root):
        local = child.tag.split("}")[-1]
        if local == "rect":
            x, y, w, h = parse_rect(child)
            fill = child.get("fill", "")
            if w > 900 and h > 700:
                background = child
            elif fill == "#111315" and w > 500 and h > 300:
                gel_panel = child
            elif fill in {"#1a2028", "#1f252e"} and w < 100 and h > 400:
                lane_rects.append(child)
            elif fill == "#e5e7eb" or fill == "#f59e0b":
                band_rects.append(child)
        elif local == "text":
            text = "".join(child.itertext()).strip()
            match = re.match(r"^(\d+)\s*bp$", text)
            if match:
                scale_labels.append((int(match.group(1)), float(child.get("y", "0"))))

    for child in list(root):
        if child.tag == qname("text") or child.tag == qname("line"):
            root.remove(child)
    if background is not None:
        root.remove(background)

    defs = ET.Element(qname("defs"))
    panel_shadow = ET.SubElement(defs, qname("filter"), {
        "id": "gel-shadow",
        "x": "-15%",
        "y": "-15%",
        "width": "130%",
        "height": "160%",
    })
    ET.SubElement(panel_shadow, qname("feDropShadow"), {
        "dx": "0.0",
        "dy": "8.0",
        "stdDeviation": "8.0",
        "flood-color": "#0f172a",
        "flood-opacity": "0.16",
    })
    panel_grad = ET.SubElement(defs, qname("linearGradient"), {
        "id": "gel-panel-grad",
        "x1": "0%",
        "y1": "0%",
        "x2": "0%",
        "y2": "100%",
    })
    ET.SubElement(panel_grad, qname("stop"), {"offset": "0%", "stop-color": "#171a1f"})
    ET.SubElement(panel_grad, qname("stop"), {"offset": "100%", "stop-color": "#0f1216"})
    root.insert(0, defs)

    if gel_panel is None:
        raise SystemExit("Could not locate gel panel in input SVG")

    gel_panel.set("fill", "url(#gel-panel-grad)")
    gel_panel.set("filter", "url(#gel-shadow)")

    lane_rects_sorted = sorted(lane_rects, key=lambda el: float(el.get("x", "0")))
    for idx, lane in enumerate(lane_rects_sorted):
        lane.set("rx", "8")
        lane.set("ry", "8")
        if idx in {0, len(lane_rects_sorted) - 1}:
            lane.set("fill", "#1d2530")
        else:
            lane.set("fill", "#202835")

    x, y, w, h = parse_rect(gel_panel)
    title = ET.Element(qname("text"), {
        "x": f"{x:.1f}",
        "y": f"{(y - 18):.1f}",
        "font-family": "monospace",
        "font-size": "18",
        "font-weight": "700",
        "fill": "#173042",
    })
    title.text = "Gibson gel readout"
    root.append(title)

    subtitle = ET.Element(qname("text"), {
        "x": f"{x:.1f}",
        "y": f"{(y - 2):.1f}",
        "font-family": "monospace",
        "font-size": "12",
        "fill": "#475569",
    })
    subtitle.text = "left ladder: 100 bp | right ladder: 1 kb"
    root.append(subtitle)

    desired_scale_values = [5000, 3000, 1000, 500, 100]
    scale_lookup = {bp: y_pos for bp, y_pos in scale_labels}
    for bp in desired_scale_values:
        y_pos = scale_lookup.get(bp)
        if y_pos is None:
            continue
        if bp >= 1000:
            label = f"{bp // 1000} kb" if bp % 1000 == 0 else f"{bp / 1000:.1f} kb"
        else:
            label = f"{bp} bp"
        for text_x, anchor in ((x - 14.0, "end"), (x + w + 14.0, "start")):
            text = ET.Element(qname("text"), {
                "x": f"{text_x:.1f}",
                "y": f"{y_pos:.1f}",
                "text-anchor": anchor,
                "font-family": "monospace",
                "font-size": "11",
                "fill": "#64748b",
            })
            text.text = label
            root.append(text)

    labels_y = y + h + 26.0
    for lane, label in zip(lane_rects_sorted, lane_labels):
        lane_x, _, lane_w, _ = parse_rect(lane)
        text = ET.Element(qname("text"), {
            "x": f"{lane_x + lane_w * 0.5:.1f}",
            "y": f"{labels_y:.1f}",
            "text-anchor": "middle",
            "font-family": "monospace",
            "font-size": "13",
            "fill": "#334155",
        })
        text.text = label
        root.append(text)

    min_x = math.inf
    min_y = math.inf
    max_x = -math.inf
    max_y = -math.inf

    for child in list(root):
        local = child.tag.split("}")[-1]
        if local == "rect":
            cx, cy, cw, ch = parse_rect(child)
            min_x = min(min_x, cx)
            min_y = min(min_y, cy)
            max_x = max(max_x, cx + cw)
            max_y = max(max_y, cy + ch)
        elif local == "text":
            tx = float(child.get("x", "0"))
            ty = float(child.get("y", "0"))
            min_x = min(min_x, tx - 80)
            min_y = min(min_y, ty - 20)
            max_x = max(max_x, tx + 180)
            max_y = max(max_y, ty + 12)

    margin = 16.0
    min_x -= margin
    min_y -= margin
    max_x += margin
    max_y += margin
    width = max_x - min_x
    height = max_y - min_y
    root.set("viewBox", f"{min_x:.2f} {min_y:.2f} {width:.2f} {height:.2f}")
    root.set("width", f"{width:.1f}")
    root.set("height", f"{height:.1f}")

    output_path.write_text(
        ET.tostring(root, encoding="unicode", short_empty_elements=True)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
