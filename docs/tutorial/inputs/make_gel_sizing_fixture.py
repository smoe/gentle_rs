#!/usr/bin/env python3
"""Draw fixed synthetic bands and requests for the measured-gel tutorial.

Uses original 160x220 pixel coordinates, not GENtle's virtual-gel model.
The uncompressed PNG deflate block has no timestamps or compressor heuristics.
No image processing, sizing algorithm, private input or third-party dependency.
"""

import argparse
import hashlib
import json
from pathlib import Path
import struct
import zlib


def png_bytes():
    pixels = [bytearray([12] * 160) for _ in range(220)]
    for x, y in [(20, 20), (20, 80), (20, 140), (20, 200), (90, 110), (90, 210)]:
        for row in range(y - 1, y + 2):
            pixels[row][x - 8:x + 9] = bytes([210] * 17)
    raw = b"".join(b"\0" + row for row in pixels)
    # One stored DEFLATE block makes exact recreation independent of zlib tuning.
    stream = (b"\x78\x01\x01" + struct.pack("<HH", len(raw), len(raw) ^ 0xFFFF)
              + raw + struct.pack(">I", zlib.adler32(raw)))

    def chunk(kind, data):
        return (struct.pack(">I", len(data)) + kind + data
                + struct.pack(">I", zlib.crc32(kind + data)))

    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", 160, 220, 8, 0, 0, 0, 0))
            + chunk(b"IDAT", stream) + chunk(b"IEND", b""))


def request(digest, protein=False):
    return {
        "report_id": "synthetic_sds_1" if protein else "synthetic_dna_1",
        "image_id": "gel_teaching", "image_sha256": digest,
        "size_kind": "sds_protein_kda" if protein else "linear_dna_bp",
        "migration": "down",
        "lanes": [
            {"id": "marker", "label": "Synthetic ladder",
             "min": {"x": 0, "y": 0}, "max": {"x": 40, "y": 219}},
            {"id": "sample", "label": "Synthetic sample",
             "min": {"x": 60, "y": 0}, "max": {"x": 130, "y": 219}},
        ],
        "ladder": {
            "lane_id": "marker", "label": "Invented calibration, not a vendor ladder",
            "source": "synthetic: explicit sizes halve every 60 pixels",
            "prestained": protein,
            "gel_system": "synthetic Tris-glycine SDS exercise" if protein else None,
            "bands": [
                {"id": f"r{i + 1}", "center": {"x": 20, "y": y},
                 "size": size / 20 if protein else size}
                for i, (y, size) in enumerate([(20, 2000), (80, 1000), (140, 500), (200, 250)])
            ],
        },
        "sample_bands": [
            {"id": "inside", "lane_id": "sample", "label": "Bracketed band",
             "center": {"x": 90, "y": 110}, "position_half_width_px": 2},
            {"id": "outside", "lane_id": "sample", "label": "Outside calibration",
             "center": {"x": 90, "y": 210}},
        ],
    }


def fixture_files():
    image = png_bytes()
    digest = hashlib.sha256(image).hexdigest()
    files = {"gel_sizing_synthetic.png": image}
    for name, protein in [("dna-analysis.json", False), ("protein-analysis.json", True)]:
        files[name] = (json.dumps(request(digest, protein), indent=2, sort_keys=True) + "\n").encode()
    files["SHA256SUMS"] = "".join(
        f"{hashlib.sha256(data).hexdigest()}  {name}\n" for name, data in files.items()
    ).encode()
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Fixture directory; never overwrites members")
    parser.add_argument("--check", action="store_true", help="Compare exact bytes without writing")
    args = parser.parse_args()
    files = fixture_files()
    if args.check:
        for name, expected in files.items():
            if (args.output / name).read_bytes() != expected:
                raise SystemExit(f"Fixture mismatch: {name}")
        print("Synthetic PNG, requests and checksums match; no GENtle execution claimed.")
    else:
        if any((args.output / name).exists() for name in files):
            raise SystemExit("Refusing to overwrite an existing fixture member")
        args.output.mkdir(parents=True, exist_ok=True)
        for name, data in files.items():
            with (args.output / name).open("xb") as handle:
                handle.write(data)


if __name__ == "__main__":
    main()
