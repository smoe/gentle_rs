# Tube visualization reference provenance

These are static reference dimensions, not physical inventory or fit approvals.
The runtime table lives in `crates/gentle-protocol/src/labware.rs`; the shared
OpenSCAD renderer lives in `crates/gentle-render/src/labware.rs`. Runtime use is
offline and does not download a drawing. No manufacturer PDF is redistributed.

## Sources

Retrieved and visually read on 2026-10-17 from Eppendorf's English technical
drawings, page 1, dimensions in millimetres. Hashes bind the raw PDF bytes, not
an extracted text representation:

| Profile | Manufacturer drawing | Raw PDF SHA-256 |
| --- | --- | --- |
| `eppendorf_safe_lock_0_5ml` | [0.5 ml, 0017010409-04/042016](https://www.eppendorf.com/product-media/doc/en/140020_Technical-Data/Eppendorf_Consumables_Technical-data_Safe-Lock-Tube-05-mL_Safe-Lock-05-mL-technical-drawing.pdf) | `256d42cddac90d2276f0552cffa7723ef607288461e05469914c270759cb928f` |
| `eppendorf_safe_lock_1_5ml` | [1.5 ml, 0017010417-05/032024](https://www.eppendorf.com/product-media/doc/en/140027_Technical-Data/Eppendorf_Consumables_Technical-data_Safe-Lock-Tube-15-mL_Safe-Lock-15-mL-technical-drawing.pdf) | `f62ad456464047ea8108fb5117638e03036be7c08880748bf9060ebdecb7fa58` |
| `eppendorf_safe_lock_2_0ml` | [2.0 ml, 0017010425-05/032024](https://www.eppendorf.com/product-media/doc/en/140034_Technical-Data/Eppendorf_Consumables_Technical-data_Safe-Lock-Tube-20-mL_Safe-Lock-20-mL-technical-drawing.pdf) | `03a3adcfe777f85662550b0d85c1a64b715ba812afc087848c33588dc09ff5c3` |

To reproduce the reference check, download each linked PDF, compute SHA-256 of
its original bytes, and inspect the callouts on page 1. A changed hash is a new
source revision requiring review, not an automatic profile update.

| Callout / interpretation | 0.5 ml | 1.5 ml | 2.0 ml |
| --- | --- | --- | --- |
| Body outer diameter | 6.7 | 9.85 | 9.85 |
| Body inner diameter | 6.1 | 8.7 | 8.8 |
| Collar outer diameter | 7.8 | 10.8 | 10.8 |
| Body length, cap open | 30 | 38.9 | 40 |
| Cylindrical section | 16.5 | 20 | 32.3 |
| Approximate closed height | 31 | 41 | 42 |
| Approximate closed side-view width | 18 | 20 | 20 |
| Approximate height with cap at 90 degrees | 48 | 59 | 60 |
| Approximate side-view width with cap at 90 degrees | 16 | 19 | 19 |

The last four fields retain the drawing's approximate (`ca.`) status. Side-view
width is **not** a cap/hinge outline viewed from above. That outline remains
unknown. Ambiguous diameter callouts 7.3/10.7/10.5 are not used. Safe-Lock 0.5 ml
is not a 0.5 ml PCR tube; these sources do not dimension the QV product.

## Model assumptions and nonclaims

The model is a solid outer body and collar, not a hollow printable tube. Cap,
hinge and cavity are omitted. Collar thickness (1 mm), conical tip radius
(0.6 mm), the conical transition and the 2.0 ml ellipsoidal bottom are explicitly
illustrative; they are not manufacturing tolerances or fully resolved drawing
geometry. The table preserves inner diameters without inventing wall thickness
through the cone or bottom. Both the portable profile and SCAD header disclose
these boundaries.

No physical fit, closed/open-cap clearance, fabrication, material, sterility,
freezer, autoclave or cell-culture qualification follows. Existing PCR-rack and
plate templates are unchanged and are not manufacturer-backed compatibility
evidence. A 96-position layout is not proof that these tubes fit an ordinary
96-well plate or an ANSI/SLAS carrier. New carriers, support mode, boxes,
inserts and lids remain deferred.

## Verification boundary

Rust unit tests check the table, portable JSON, parameterized SCAD and shared
engine/shell/GUI paths. `python3 -m unittest scripts.test_labware_models -v`
independently exercises the actual SCAD body template with literal drawing
parameters through OpenSCAD, if installed. It checks finite, nonempty mesh
geometry and outer bounds; it does not run the Rust renderer or prove fit.
Synthetic test racks and numeric cases are hand-crafted in those tests; no
private sample data is used.

Legacy carrier byte-exact baselines still need to be captured by Glen/CI from
the pre-change source `1c5be3175323eec41359b980b5a41e22c7edaec1`, for each of
the three existing templates on one fixed rack (SCAD, simulation JSON and spec
JSON). Compare them with this change before accepting legacy byte identity.
If those baselines are committed, include provenance, scoped LF checkout policy
and LF/CRLF checkout coverage. No legacy geometry formula or template identity
is intentionally changed here; source inspection is not a byte-exact pass.
