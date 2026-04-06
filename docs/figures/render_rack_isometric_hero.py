#!/usr/bin/env python3
"""Derive a calmer README hero SVG from a technical rack isometric export.

The input is the deterministic SVG produced by `racks isometric-svg`. This
script keeps the real rack geometry and occupied positions, but removes most
technical annotations and applies softer gradients/shadows so the resulting
figure reads more like an object hero than a diagram.
"""

from __future__ import annotations

import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def qname(tag: str) -> str:
    return f"{{{SVG_NS}}}{tag}"


def parse_points(points: str) -> list[tuple[float, float]]:
    result: list[tuple[float, float]] = []
    for token in points.split():
        x_str, y_str = token.split(",")
        result.append((float(x_str), float(y_str)))
    return result


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgb_to_hex(rgb: tuple[int, int, int]) -> str:
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def mix(color: str, target: str, amount: float) -> str:
    base = hex_to_rgb(color)
    tgt = hex_to_rgb(target)
    mixed = tuple(
        max(0, min(255, round(base[i] * (1.0 - amount) + tgt[i] * amount)))
        for i in range(3)
    )
    return rgb_to_hex(mixed)


def set_points(element: ET.Element, points: list[tuple[float, float]]) -> None:
    element.set(
        "points",
        " ".join(f"{x:.2f},{y:.2f}" for x, y in points),
    )


def main() -> int:
    if len(sys.argv) != 3:
        print(
            "usage: render_rack_isometric_hero.py INPUT.svg OUTPUT.svg",
            file=sys.stderr,
        )
        return 2

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    tree = ET.parse(input_path)
    root = tree.getroot()
    children = list(root)

    background = None
    polygons: list[ET.Element] = []
    rects: list[ET.Element] = []
    ellipses: list[ET.Element] = []
    paths: list[ET.Element] = []
    title_texts: list[ET.Element] = []
    footer_texts: list[ET.Element] = []

    for child in children:
        local = child.tag.split("}")[-1]
        if local == "rect":
            rects.append(child)
            if child.get("width") == "100%" and child.get("height") == "100%":
                background = child
        elif local == "polygon":
            polygons.append(child)
        elif local == "ellipse":
            ellipses.append(child)
        elif local == "path":
            paths.append(child)
        elif local == "text":
            text = "".join(child.itertext()).strip()
            if text == "GENtle rack isometric sketch":
                title_texts.append(child)
            elif text.startswith("template="):
                footer_texts.append(child)

    if background is not None:
        root.remove(background)
    for element in title_texts + footer_texts:
        root.remove(element)

    legend_start = None
    legend_end = None
    for idx, child in enumerate(list(root)):
        if child.tag == qname("text") and "".join(child.itertext()).strip() == "Occupancy":
            legend_start = idx
            break
    if legend_start is not None:
        legend_end = len(list(root))
        for idx in range(legend_start, len(list(root))):
            child = list(root)[idx]
            if child.tag == qname("text"):
                text = "".join(child.itertext()).strip()
                if text.startswith("template="):
                    legend_end = idx
                    break
        for child in list(root)[legend_start:legend_end]:
            root.remove(child)

    # Remove coordinate labels and front-strip status text. Keep only geometry.
    for child in list(root):
        if child.tag != qname("text"):
            continue
        text = "".join(child.itertext()).strip()
        if (
            text in {str(n) for n in range(1, 97)}
            or len(text) == 1
            and text.isalpha()
            or text.startswith("rack-")
            or "occupied" in text
        ):
            root.remove(child)

    # Geometry assumptions follow the deterministic export layout:
    # polygons[0..2] are top/front/side, rects[1] is the front strip.
    rack_top, rack_front, rack_side = polygons[:3]
    strip_rect = None
    for rect in list(root):
        if rect.tag != qname("rect"):
            continue
        if rect.get("rx") == "1.2" and rect.get("ry") == "1.2":
            strip_rect = rect
            break

    top_fill = rack_top.get("fill", "#e7f6f2")
    front_fill = rack_front.get("fill", "#b8dfd5")
    side_fill = rack_side.get("fill", "#8ecaba")
    edge_stroke = rack_top.get("stroke", "#115e59")
    strip_fill = strip_rect.get("fill", "#0f766e") if strip_rect is not None else "#0f766e"

    defs = ET.Element(qname("defs"))

    rack_shadow = ET.SubElement(defs, qname("filter"), {
        "id": "rack-shadow",
        "x": "-20%",
        "y": "-20%",
        "width": "160%",
        "height": "180%",
    })
    ET.SubElement(rack_shadow, qname("feDropShadow"), {
        "dx": "0.0",
        "dy": "2.4",
        "stdDeviation": "2.2",
        "flood-color": "#0f172a",
        "flood-opacity": "0.18",
    })

    top_grad = ET.SubElement(defs, qname("linearGradient"), {
        "id": "rack-top-grad",
        "x1": "0%",
        "y1": "0%",
        "x2": "100%",
        "y2": "100%",
    })
    ET.SubElement(top_grad, qname("stop"), {"offset": "0%", "stop-color": mix(top_fill, "#ffffff", 0.30)})
    ET.SubElement(top_grad, qname("stop"), {"offset": "100%", "stop-color": mix(top_fill, "#0f172a", 0.08)})

    front_grad = ET.SubElement(defs, qname("linearGradient"), {
        "id": "rack-front-grad",
        "x1": "0%",
        "y1": "0%",
        "x2": "0%",
        "y2": "100%",
    })
    ET.SubElement(front_grad, qname("stop"), {"offset": "0%", "stop-color": mix(front_fill, "#ffffff", 0.22)})
    ET.SubElement(front_grad, qname("stop"), {"offset": "100%", "stop-color": mix(front_fill, "#0f172a", 0.10)})

    side_grad = ET.SubElement(defs, qname("linearGradient"), {
        "id": "rack-side-grad",
        "x1": "0%",
        "y1": "0%",
        "x2": "0%",
        "y2": "100%",
    })
    ET.SubElement(side_grad, qname("stop"), {"offset": "0%", "stop-color": mix(side_fill, "#ffffff", 0.18)})
    ET.SubElement(side_grad, qname("stop"), {"offset": "100%", "stop-color": mix(side_fill, "#0f172a", 0.12)})

    strip_grad = ET.SubElement(defs, qname("linearGradient"), {
        "id": "rack-strip-grad",
        "x1": "0%",
        "y1": "0%",
        "x2": "100%",
        "y2": "0%",
    })
    ET.SubElement(strip_grad, qname("stop"), {"offset": "0%", "stop-color": mix(strip_fill, "#ffffff", 0.16)})
    ET.SubElement(strip_grad, qname("stop"), {"offset": "100%", "stop-color": mix(strip_fill, "#0f172a", 0.10)})

    hole_grad = ET.SubElement(defs, qname("radialGradient"), {
        "id": "hole-grad",
        "cx": "50%",
        "cy": "45%",
        "r": "70%",
    })
    ET.SubElement(hole_grad, qname("stop"), {"offset": "0%", "stop-color": "#ffffff"})
    ET.SubElement(hole_grad, qname("stop"), {"offset": "55%", "stop-color": "#f1f5f9"})
    ET.SubElement(hole_grad, qname("stop"), {"offset": "100%", "stop-color": "#cbd5e1"})

    shadow_grad = ET.SubElement(defs, qname("radialGradient"), {
        "id": "ground-shadow-grad",
        "cx": "50%",
        "cy": "50%",
        "r": "50%",
    })
    ET.SubElement(shadow_grad, qname("stop"), {"offset": "0%", "stop-color": "#0f172a", "stop-opacity": "0.16"})
    ET.SubElement(shadow_grad, qname("stop"), {"offset": "100%", "stop-color": "#0f172a", "stop-opacity": "0"})

    root.insert(0, defs)

    rack_top.set("fill", "url(#rack-top-grad)")
    rack_front.set("fill", "url(#rack-front-grad)")
    rack_side.set("fill", "url(#rack-side-grad)")
    for surface in (rack_top, rack_front, rack_side):
        surface.set("filter", "url(#rack-shadow)")
        surface.set("stroke-width", "0.45")

    if strip_rect is not None:
        strip_rect.set("fill", "url(#rack-strip-grad)")
        strip_rect.set("fill-opacity", "0.98")
        strip_rect.set("stroke-width", "0.22")

    holes = [
        ellipse
        for ellipse in list(root)
        if ellipse.tag == qname("ellipse")
        and ellipse.get("fill") == "#ffffff"
        and ellipse.get("stroke") == "#475569"
    ]
    for hole in holes:
        hole.set("fill", "url(#hole-grad)")
        hole.set("stroke", mix("#475569", "#94a3b8", 0.45))
        hole.set("stroke-width", "0.28")

    # Slightly deepen tube side walls and keep highlights.
    for rect in list(root):
        if rect.tag != qname("rect") or rect is strip_rect:
            continue
        fill = rect.get("fill", "")
        if fill.startswith("#"):
            rect.set("fill", mix(fill, "#0f172a", 0.10))
            rect.set("fill-opacity", "0.95")

    for path in list(root):
        if path.tag != qname("path"):
            continue
        if path.get("data-rack-tube-shell") == "1":
            fill = path.get("fill", "")
            if fill.startswith("#"):
                path.set("fill", mix(fill, "#0f172a", 0.10))
                path.set("fill-opacity", "0.95")
        elif path.get("data-rack-tube-interface") == "1":
            stroke = path.get("stroke", "")
            if stroke.startswith("#"):
                path.set("stroke", mix(stroke, "#0f172a", 0.12))
                path.set("stroke-opacity", "0.42")
                path.set("stroke-width", "0.34")

    for ellipse in list(root):
        if ellipse.tag != qname("ellipse"):
            continue
        fill = ellipse.get("fill", "")
        if fill.startswith("#") and fill not in {"#ffffff", "#0f172a"}:
            ellipse.set("fill", mix(fill, "#ffffff", 0.08))
            ellipse.set("stroke-width", "0.20")

    top_points = parse_points(rack_top.get("points", ""))
    front_points = parse_points(rack_front.get("points", ""))
    side_points = parse_points(rack_side.get("points", ""))
    top_front_left = top_points[0]
    top_front_right = top_points[1]
    top_back_right = top_points[2]
    top_back_left = top_points[3]
    front_bottom_right = front_points[2]
    front_bottom_left = front_points[3]

    shadow = ET.Element(qname("ellipse"), {
        "cx": f"{(front_bottom_left[0] + top_back_right[0]) * 0.5:.2f}",
        "cy": f"{front_bottom_left[1] + 4.8:.2f}",
        "rx": f"{(top_back_right[0] - top_front_left[0]) * 0.44:.2f}",
        "ry": f"{(front_bottom_left[1] - top_back_left[1]) * 0.18:.2f}",
        "fill": "url(#ground-shadow-grad)",
    })
    root.insert(1, shadow)

    # Add one soft top-face highlight.
    highlight = ET.Element(qname("path"), {
        "d": (
            f"M {top_front_left[0] + 6.0:.2f} {top_front_left[1] + 1.7:.2f} "
            f"C {top_front_left[0] + 24.0:.2f} {top_front_left[1] - 2.2:.2f}, "
            f"{top_front_right[0] - 24.0:.2f} {top_front_right[1] - 2.0:.2f}, "
            f"{top_back_right[0] - 8.0:.2f} {top_back_right[1] + 1.2:.2f}"
        ),
        "fill": "none",
        "stroke": "#ffffff",
        "stroke-opacity": "0.35",
        "stroke-width": "0.9",
        "stroke-linecap": "round",
    })
    root.insert(2, highlight)

    # Tight crop around the rack object.
    min_x = math.inf
    min_y = math.inf
    max_x = -math.inf
    max_y = -math.inf

    for child in list(root):
        local = child.tag.split("}")[-1]
        if local == "polygon":
            for x, y in parse_points(child.get("points", "")):
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                max_x = max(max_x, x)
                max_y = max(max_y, y)
        elif local == "rect":
            x = float(child.get("x", "0"))
            y = float(child.get("y", "0"))
            w = float(child.get("width", "0"))
            h = float(child.get("height", "0"))
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            max_x = max(max_x, x + w)
            max_y = max(max_y, y + h)
        elif local == "ellipse":
            cx = float(child.get("cx", "0"))
            cy = float(child.get("cy", "0"))
            rx = float(child.get("rx", "0"))
            ry = float(child.get("ry", "0"))
            min_x = min(min_x, cx - rx)
            min_y = min(min_y, cy - ry)
            max_x = max(max_x, cx + rx)
            max_y = max(max_y, cy + ry)

    margin = 6.0
    min_x -= margin
    min_y -= margin
    max_x += margin
    max_y += margin
    width = max_x - min_x
    height = max_y - min_y
    root.set("viewBox", f"{min_x:.2f} {min_y:.2f} {width:.2f} {height:.2f}")
    root.set("width", f"{width:.1f}mm")
    root.set("height", f"{height:.1f}mm")

    output_path.write_text(
        ET.tostring(root, encoding="unicode", short_empty_elements=True)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
