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
support modes below and the explicit compact override. Physical box/insert/lid
profiles still need their own evidence; no material qualification follows.

## Safe-Lock 1.5/2.0 ml carriers

`safe_lock_tube_rack_collar_seated` suspends the tube at its collar;
`safe_lock_tube_rack_floor_seated` rests its tip on the floor. Both derive
worst-case dimensions over both manufacturer profiles, not nominal volume.
They are `tube_rack`, format `eppendorf_safe_lock_1_5ml_or_2_0ml`, with
`geometry_basis=drawing_derived_unverified`. They do not modify the three older
templates. New optional spec fields are omitted for old templates.

### Explicit Compact Spacing

`Rack.carrier_pitch_override_um` is optional, integer micrometres. Without it,
all existing template numbers and serialized rack/export bytes stay unchanged.
Use `racks set-carrier-pitch RACK_ID 14` or `--clear`; grid dimensions are a
separate `racks set-custom-profile` choice, not a box-format assumption.

The pure shared derivation retains conservative outer edge margins (4.925 mm
collar-seated; 4.775 mm floor-seated) and changes only inter-hole spacing.
`MIN_WEB=0.8 mm` is a per-side allowance: the admitted pitch must exceed the
opening by **two** such allowances (1.6 mm total). `COLLAR_GAP=0.5 mm` is the
separate adjacent-collar allowance. Both are engineering assumptions, not
manufacturing tolerances or measured mechanical strength.

| Admission order | Reason prefix | Limit |
| --- | --- | --- |
| 1 | `compact_pitch_invalid` | finite and positive |
| 2 | `compact_pitch_not_compact` | at most conservative 20 mm |
| 3 | `compact_pitch_below_web` | opening + 2 x MIN_WEB: 11.75/12.05 mm collar/floor |
| 4 | `compact_pitch_collar_overlap` | maximum collar diameter + COLLAR_GAP: 11.3 mm |

The setter admits collar-compatible pitch because the rack has no support-mode
identity; every export rechecks its selected mode. Thus 11.75 mm is admitted on
the rack but a floor-seated export refuses it. Other templates, including PCR96,
refuse `carrier_pitch_override_unsupported` before any write. Clearing restores
the old output. Row/column changes retain the saved override.

Worked synthetic geometry: a 9x9 floor-seated layout at 14 mm has width
`2 x 4.775 + 10.45 + 8 x 14 = 132 mm`; A1 centre is `(10, 15.225)` mm and
I9 is `(122, 127.225)` mm. These are geometry assertions only, not box dimensions.
Compact exports use `drawing_derived_unverified_compact` and disclose possible
closed-cap/hinge collisions. Keeping outer margins does not qualify internal
cap clearance, physical fit, material, freezing, ice, sterility or autoclaving.

### Different Storage-Box Formats

The owner requires different formats, not a single 9x9 convention. Treat a
commercial paper/cardboard box as a possible destination: independent row/column
grid, usable interior width/depth, divider layout/thickness, interior height,
closed-cap height and lid constraints. Square/rectangular custom grids and
position labels already exist; they are not a physical box profile. Do not infer
fit from nominal tube volume, a chosen grid, material name or generated pitch.
Source/measured box profiles and compatibility remain separate future work.

Owner clarification, 2026-10-20: regular 1.5/2.0 ml tubes rest on the box bottom.
The **openings are square in top view**, not a claim about divider cross-section.
Tubes should fit snugly against the divider walls; adjacent cells are separated
only by wall thickness. Literal tube-to-tube contact is not required. Any
eventual model must distinguish that desired fit from sourced dimensions and
explicit handling/manufacturing tolerances.

The owner supplied [Semadeni's cardboard cryobox example](https://eshop.semadeni.com/de/PRODUKTE/Life-Science_-Proben-nehmen-und-Analysieren/Kryo-Zubehoer/Kryo-Boxen---Karton-9290.html)
(consulted 2026-10-20) as a storage concept, not a purchase or fit endorsement.
The page describes an integrated grid and cardboard lid, with 81/100-position
variants; article 6345 lists 9x9, **exterior** 133x133x50 mm and nominal 1.2/2 ml
tubes. It does not supply clear cell width, divider thickness or usable internal
height/lid clearance. Do not divide exterior width by nine to invent cell size,
or interpret that volume label as evidence of Safe-Lock 1.5/2.0 ml compatibility.
No supplier temperature statement is adopted as GENtle qualification, no vendor
profile or runtime dependency is added, and no model dimensions change.

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

## Derived Storage Box

`safe_lock_storage_box_derived` uses the existing rack's admitted row/column
grid (up to 1,536 positions), `family=storage_box`, floor support and both
1.5/2.0 ml profiles. It is a **GENtle-derived design**, not a model of any
commercial box. `geometry_basis=drawing_derived_unverified_storage_box` and
the optional `storage_box` record preserve its dimensions/assumptions in the
spec and simulation JSON. The field is absent from every older template.

| Engineering assumption (mm) | Value |
| --- | --- |
| `RADIAL_CLEARANCE_SNUG` | 0.15 per side |
| `BOX_DIVIDER_THICKNESS` | 1.2; chosen above the 1.15 collar-gap minimum |
| `BOX_WALL_THICKNESS` | 2.0, walls and floor |
| `BOX_LID_CLEARANCE` | 3.0 above approximate closed height |
| `LID_FIT_CLEARANCE` | 0.3 per side |
| `LID_TOP_THICKNESS` / `LID_SKIRT_HEIGHT` | 2.0 / 8.0 |

| Derived dimension | Rule / result for both profiles |
| --- | --- |
| Square cell | max body diameter + 2 radial clearances = 10.15 mm |
| Pitch | cell + divider = 11.35 mm; collar diameter + 0.5-mm gap is 11.3 mm |
| Divider height | existing validated floor-seated bore depth = 22.0 mm |
| Interior height | max approximate (`ca.`) closed height + 3 = 45 mm |
| Interior width/depth | n cells + (n - 1) dividers |
| Exterior width/depth | interior + two 2-mm walls |
| Exterior height | 2-mm floor + interior height = 47 mm |
| Lid inner width/depth | box exterior + twice 0.3-mm clearance |
| Lid exterior width/depth | lid interior + two 2-mm walls |
| Lid height | 8-mm skirt + 2-mm top = 10 mm |

Worked 9x9 example: interior 100.95 x 100.95 mm, exterior 104.95 x 104.95 x
47 mm; A1 centre (7.075, 7.075), I9 (97.875, 97.875) mm. A 4x12 grid instead
has interior 135 x 44.2 mm and exterior 139 x 48.2 x 47 mm; 1x1 is also valid.
These are derived design dimensions, not measured fits. A tube rests on the
box floor and dividers end below the collar; snugness is a target, not evidence.
Closed caps/hinges may collide at 11.35-mm pitch because the approximate closed
side-view envelope is 20 mm. No cap clearance, freezer/ice/thermal, cardboard,
material or lid-retention suitability is established.

Use the same `racks ... --template safe_lock_storage_box_derived` export routes
and GUI physical-template selector. Clear any saved compact-carrier pitch first;
otherwise every export refuses `carrier_pitch_override_unsupported` before
writing. Grid and source-geometry errors propagate; `box_collar_overlap` and
`cap_height_unknown` are explicit refusal prefixes. Extremely long synthetic
cones can make the reused carrier derivation refuse `collar_not_above_top`;
do not rename that existing diagnostic as cone-only engagement.

SCAD provides `gentle_storage_box`, `gentle_storage_box_insert` and
`gentle_storage_box_lid`; the lid is displayed alongside the box, not fitted
on it. All grid positions remain physical cells, including blocked positions.
Optional OpenSCAD `part="box"`, `"insert"` or `"lid"` isolates a component;
the default `"all"` includes all parts. The ten-mm display separation is not
a fit dimension. SVG views retain round markers with the visible line
**Square cells drawn as round markers**; use SCAD for the square-cell design.
Source-matched Rust exports/native controls and measured fit remain Glen/CI
or physical validation work. The Python harness checks the actual SCAD template
with hand-derived dimensions, not execution of the Rust derivation.

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

## PCR 0.2 ml drawing and 96-position carrier

On 2026-10-20, Eppendorf's [product page](https://www.eppendorf.com/de-en/Products/PCR/PCR-Consumables/PCR-Tubes-PCR-Tube-Strips-Cap-Strips-p-PF-8634)
confirms individual 0.2 ml tubes (0030124332). Its Technical Drawings tab links
[Technical Data Sheet 0012 595 670-00/2026-05](https://www.eppendorf.com/product-media/doc/en/17413036/Consumables_Technical-data_PCR-Tubes_PCR-Tube-Strips_Cap-Strips.pdf).
Raw PDF SHA-256: `67b84dea22bcb2cbea26907b8f1116c506210f1917744feb64b2fafb6745b40d`.
The PDF is not redistributed. Re-download, hash the original bytes and visually
inspect **page 4** (individual 0.2 ml tube); page 7 depicts the 0.2 ml strip.
No 0.1 ml or 0.5 ml PCR tube is added.

Page 4 callouts: body length 20.8, outer conical transition at 11.1 from the
top, main cylindrical OD 5.8, upper OD 6.0, upper inner diameter 5.2, all mm.
The profile conservatively uses **6.0 mm maximum OD**, not just the 5.8-mm
lower cylinder. The rim does not document a 0.20-mm radial seating ledge over
the main body; `collar_outer_diameter` remains `not_documented`. Closed height,
closed-cap width, 90-degree cap geometry and cap plan outline are also unknown.
The **22-mm open-cap side-view width is not a closed-cap envelope**. The
reference profile is inspectable; the current body-and-collar visualization
refuses this profile with `tube_model_dimension_unavailable` rather than
inventing a collar. Missing/nonfinite dimensions in future profiles now refuse
before file writes instead of panicking.

`pcr_tube_96_position_carrier` is floor-seated only, `family=tube_rack`, format
`eppendorf_pcr_tube_0_2ml`. It requires exactly 8 rows and 12 columns, including
custom grids; 8x11, 12x8 and 16x24 refuse `carrier_requires_8x12` before writes.
The same floor-support constraints produce an opening of 6.0 + 2x0.30 = 6.6 mm,
a 13-mm bore (strictly above 20.8 - 11.1 + 3), and 2-mm floor / 15-mm total
height. Using the maximum upper OD is a deliberate conservative interpretation
of the drawing, not a changed clearance constant.

Nominal positions were checked against official SLAS standard text:
[ANSI/SLAS 1-2004 (R2012), section 4.1.1, page 4](https://www.slas.org/SLAS/assets/File/public/standards/ANSI_SLAS_1-2004_FootprintDimensions.pdf)
and [ANSI/SLAS 4-2004 (R2012), sections 4.1.1-4.1.3, pages 4-5](https://www.slas.org/SLAS/assets/File/public/standards/ANSI_SLAS_4-2004_WellPositions.pdf).
The spec directly sets 127.76 x 85.48 mm, 9-mm pitch, A1 (14.38, 11.24) and
H12 (113.38, 74.24). `edge_margin=11.08` and `front_top_clearance=7.94`
retain those centers in the existing coordinate generator; the old symmetric
size formula is not reused. Outer-wall parameter 2.0, label depth 6.0,
recess 0.8 and corner parameter 2.0 mm are engineering assumptions. Full
standard compliance, flange/holder fit, thermal contact and thermocycler use
are **not** established. Nine-mm cap/strip compatibility is an explicit,
unverified format assumption, not a measured single-tube closed-cap clearance.

`scripts.test_labware_models` additionally exercises the actual rack SCAD
format strings with literal expected dimensions for both Safe-Lock modes and
PCR96. It checks finite nonempty meshes and outer bounds, with no Rust build.
This does not run the Rust derivation/exporter and cannot replace Glen's
source-matched export, legacy-byte or native/physical checks.

## Cut PCR strips, 0.2 ml per tube

`brand_pcr_strip_0_2ml_781320` models **one through eight contiguous tubes**,
default eight, following the owner's 2026-10-22 request to represent cut strips.
It is distinct from the individual Eppendorf PCR profile and from rack carriers.
Counts do not assign sample inventory or claim a new purchasable manufacturer
article. The model omits separate caps, cavity, numbering and end tab/cut stubs;
only N-1 internal bridges are retained, including zero bridges for one tube.

Source: [BRAND's product page](https://shop.brand.de/de/pcr-strips-pp-bio-certr-pcr-quality-p545.html)
and its linked [German technical datasheet](https://shop.brand.de/media/import/1/27/32406/42485/53752/56329/55264/Data_sheet_PCR_strips_DE.pdf),
**page 6, article 781320**, retrieved and visually inspected 2026-10-22.
No revision/date is printed on that page. Raw PDF SHA-256:
`8ba466736960b6933c642568d03ac5901c09010db176b7278d6e2989ffbd3c13`.
Re-download the PDF, hash its original bytes and inspect page 6 to reproduce
the check; the PDF is not redistributed or fetched at runtime.

| Drawing callout (mm) | Value |
| --- | --- |
| Adjacent tube centers | 9 |
| First-to-last centers, eight tubes | 63 |
| Overall uncapped height | 20.75 |
| Main body OD | 6 |
| Collar OD | 8 |
| Upper rim OD | 7.05 |
| Lower OD before rounded tip | 2.75 |

The catalogue records 0.2 ml **per tube** and a single-bridge standard strip;
the 0.15-ml low-profile, 12-tube and attached-cap variants are not added.
Cylindrical-section length, collar elevation and closed-cap envelope are not
dimensioned sufficiently for this model and remain unknown. The dedicated
strip renderer uses the known envelope with these **illustrative assumptions**:
taper transition at 55% of height, spherical tip radius 2.75/2, collar bottom
4.5 mm below the rim with 0.9-mm thickness, rim thickness 0.6 mm, and rectangular
bridges 1.0 mm wide by 0.7 mm high at the collar. Do not read the drawing's
end-tab dimensions as bridge dimensions. Cut ends are idealized with no stubs.
This is a solid visualization, not a printable vessel, measured cavity,
seating model, cap-clearance test or thermal/material/sterility qualification.
Other manufacturer profiles and their missing-dimension refusals are unchanged.

Rust tests authored for this addition cover all counts, invalid-count refusal
before overwrite, serde compatibility, deterministic exports and the shared
GUI/shell/engine path. They require Glen/CI under the no-local-Rust instruction:
`cargo test -q -p gentle-protocol labware`,
`cargo test -q -p gentle-render labware`, `cargo test -q --lib labware`, then
`cargo check -q` and `cargo fmt --check`. Native GUI interaction remains separate.
The existing Python/OpenSCAD harness reads the actual strip template with
literal source dimensions and declared assumptions; it checks all eight cut
lengths, finite bounds, connected geometry and rejection of zero/nine tubes.
These synthetic cases live in `scripts/test_labware_models.py` and reproduce
with `python3 -m unittest scripts.test_labware_models -v`. They are not execution
of the Rust exporter or evidence of physical fit.

On 2026-10-22, owner-approved OpenSCAD PNG previews of one, three and eight tubes
were visually inspected: complete framing, the expected tube counts and only
internal bridges. These previews use `strip_template` from the Python harness,
not an executed Rust export. The solid tops are outer-envelope surfaces, not
caps; cavities are intentionally absent. Preview PNG/SCAD files remain outside
the repository and establish neither native GUI acceptance nor physical fit.
