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
96-well plate or an ANSI/SLAS carrier. Safe-Lock carriers now offer the explicit
support modes below. Boxes, inserts, lids and denser unchecked spacing stay deferred.

## Safe-Lock 1.5/2.0 ml carriers

`safe_lock_tube_rack_collar_seated` suspends the tube at its collar;
`safe_lock_tube_rack_floor_seated` rests its tip on the floor. Both derive
worst-case dimensions over both manufacturer profiles, not nominal volume.
They are `tube_rack`, format `eppendorf_safe_lock_1_5ml_or_2_0ml`, with
`geometry_basis=drawing_derived_unverified`. They do not modify the three older
templates. New optional spec fields are omitted for old templates.

All constants below are **engineering assumptions**, not manufacturer
tolerances, fabrication guidance or physical qualification:

| Assumption (mm) | Value / derivation |
| --- | --- |
| `RADIAL_CLEARANCE_COLLAR` | 0.15; opening 9.85 + 2x0.15 = 10.15 |
| `RADIAL_CLEARANCE_FLOOR` | 0.30; opening 9.85 + 2x0.30 = 10.45 |
| `MIN_COLLAR_LEDGE` | 0.20 radially; opening must not exceed 10.8 - 0.4 |
| `TIP_CLEARANCE` | 2.0; suspended bore depth 40 + 2 = 42 |
| `CYLINDER_ENGAGEMENT` | strictly more than 3.0; longest lower part 38.9 - 20 = 18.9 |
| Floor bore rounding | strictly next whole mm above 18.9 + 3 = 21.9, hence 22 |
| `COLLAR_TOP_MARGIN` | strictly more than 5.0 above a floor-seated bore; depth < 38.9 - 5 |
| Floor thickness | 2.0; total heights 44.0 / 24.0 |
| Closed-cap pitch | 20.0 both axes; conservative approximate side-view width, **not** measured cap outline |
| Edge margin | (pitch - opening)/2; 4.925 / 4.775 |
| Front clearance / label-strip depth / recess | 10.0 / 6.0 / 0.8 |
| Corner-radius parameter | 2.0; existing straight-edged SCAD block remains unchanged |

The existing overall-size formula is retained: a 4x9 grid is 180.0 mm wide
and 85.075 / 85.225 mm deep. A1 is (10, 15.075 / 15.225); D9 is
(170, 75.075 / 75.225). Closed-cap clearance remains unverified. Straight
bores preserve the floor and open top; no material, ice/freezer/autoclave,
sterility or commercial-rack compatibility claim follows.

Shared derivation rejects absent/nonfinite/unproven dimensions and unknown
side-view envelopes. `TubeCarrierGeometry::validate` also checks nominated
bores against synthetic profiles: `collar_ledge_insufficient`,
`bore_engages_cone_only`, `collar_not_above_top`, `cap_envelope_unknown` retain
stable prefixes. Cone-only regression validates a proposed 22-mm bore against
a synthetic longer cone; derivation otherwise raises the bore when possible.
Every new-template SVG/SCAD/JSON includes unverified geometry disclosures.

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
The new carrier Rust tests are authored but not locally executed. Native
GUI controls, fresh source-matched exports and physical fit remain unverified.
