# Measured-Gel Teaching Inputs

Origin: hand-drawn synthetic 160 x 220, 8-bit grayscale PNG, created on
2026-10-12. No photograph, biological sample, vendor image or GENtle virtual-gel
migration algorithm is involved. Rectangles were chosen in original pixel
coordinates: background 12, bands 210, centers (20,20), (20,80), (20,140),
(20,200), (90,110) and (90,210). Each band spans x +/-8 and y +/-1 inclusive.

Used by tutorial 10.03, `scripts/test_tutorial_gap_lessons.py`, and the
LF/CRLF/binary-preservation gate in `scripts/test_tutorial_checkouts.py`.
Recreate all four generated members in a fresh directory:

```sh
python3 docs/tutorial/inputs/make_gel_sizing_fixture.py /tmp/new-gel-teaching
python3 docs/tutorial/inputs/make_gel_sizing_fixture.py --check docs/tutorial/inputs/gel_sizing
```

`SHA256SUMS` binds exact PNG and request bytes; it does not bind itself. A stored
uncompressed DEFLATE block makes the PNG independent of compressor heuristics.
Requests/checksums use UTF-8/LF. The fixture has no EXIF rotation or other image
transformation. Image coordinates are zero-based original pixel centers.

`dna-analysis.json` uses invented 2000/1000/500/250 bp reference sizes at
y=20/80/140/200. At y=110, log-size interpolation gives sqrt(1000*500), about
707.106781 bp. y=210 is outside calibration and must not receive an estimate.
`protein-analysis.json` reuses the image with sizes divided by 20 and explicitly
synthetic SDS/prestained-system labels: the inner band is about 35.355339 apparent
kDa. These assignments are not calibrated commercial ladder values and never
measure band identity, abundance or concentration.

JSON files are input requests, not engine results, approvals or execution
receipts. `image_sha256` is derived from the generated original PNG. GENtle's
import must return this exact digest before analysis. Runtime reports, exports
and project states belong in a fresh disposable directory. Native GUI and
scientific acceptance remain with Glen/the owner.
