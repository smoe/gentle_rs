# Documentation Figures

`gentle_system_overview.svg` is the hand-authored homepage schematic used near
the top of the README. Unlike the cloning and analysis showcases below, it is
not rendered from an engine operation; it is a maintained explanatory figure
that summarizes how interfaces, imported context, the shared engine, and
provenance interact. The Mermaid block in `README.md` mirrors the same
structure in text-native form.

`gibson_two_fragment_protocol_cartoon.svg` is a deterministic render of the
built-in protocol cartoon `gibson.two_fragment`.

Regenerate these rack/gel companion figures from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  gibson.two_fragment \
  docs/figures/gibson_two_fragment_protocol_cartoon.svg
```

`gibson_single_insert_protocol_cartoon.svg` is a deterministic render of the
built-in protocol cartoon `gibson.single_insert_dual_junction`. It is the
canonical single-insert Gibson mechanism strip used in the README to show both
destination-insert junctions explicitly.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  gibson.single_insert_dual_junction \
  docs/figures/gibson_single_insert_protocol_cartoon.svg
```

`pcr_pair_protocol_cartoon.svg` is the deterministic success-path render of the
built-in protocol cartoon `pcr.assay.pair`. It keeps the PCR story
selection-first: context + highlighted ROI first, then explicit forward/reverse
constraint windows (red/blue) around the ROI (green), followed by readable
primer placement and accepted ROI-bounded amplicon outcome.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  pcr.assay.pair \
  docs/figures/pcr_pair_protocol_cartoon.svg
```

`pcr_pair_no_product_protocol_cartoon.svg` is the deterministic report-only
render of `pcr.assay.pair.no_product`. It shows the same ROI-centered workflow
when no accepted primer pair yields a product.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  pcr.assay.pair.no_product \
  docs/figures/pcr_pair_no_product_protocol_cartoon.svg
```

`pcr_pair_with_tail_protocol_cartoon.svg` is the deterministic tailed-PCR strip
for `pcr.assay.pair.with_tail`. It is insertion-first:

- panel 1: extension sequences + requested insertion anchors
- panel 2: anchor-adjacent primer windows
- panel 3: tailed primer placement with shift-compensation semantics
- panel 4: amplified product carrying inserted extensions

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  pcr.assay.pair.with_tail \
  docs/figures/pcr_pair_with_tail_protocol_cartoon.svg
```

`pcr_overlap_extension_substitution_fig1_style.svg` is a deterministic render
of built-in protocol cartoon `pcr.oe.substitution`
(template source:
`docs/examples/protocol_cartoon/oe_substitution_figure1_template.json`).
It currently captures Figures 1-6 in the same renderer style used by built-in
Gibson cartoons:

- Figure 1: starting dsDNA segments (target + insert), substitution boundaries,
  and outer secondary primers (`a`,`f`).
- Figure 2: primer-design assignment of all six primers (`a`..`f`), including
  chimeric overhang intent for `b`,`c`,`d`,`e`.
- Figure 3: three independent first-step PCR products: AB (`a+b`), CD
  (`c+d`), and EF (`e+f`).
- Figure 4: denaturation into single strands for AB/CD/EF products.
- Figure 5: overlap annealing (b↔c, d↔e) into one annealed intermediate.
- Figure 6: polymerase fills annealed gaps to form the continuous duplex
  template for second-step amplification.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  pcr.oe.substitution \
  docs/figures/pcr_overlap_extension_substitution_fig1_style.svg
```

`qpcr_assay_protocol_cartoon.svg` is the deterministic probe-bearing qPCR strip
for `pcr.assay.qpcr`. It reuses the same PCR family layout while adding an
explicit third probe window (separate from ROI semantics), keeps panel 2
constraint-only, and shows panels 3/4 on amplified amplicon geometry (no
genomic flanks) with the probe constrained inside its dedicated window.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  pcr.assay.qpcr \
  docs/figures/qpcr_assay_protocol_cartoon.svg
```

`pcr_blunt_vector_ligation_template.json` is the small custom protocol-cartoon
template for the blunt-end cloning README hero. It is intentionally more
minimal than the Gibson strips: one blunt PCR product, one intact circular
pGEX vector with a separately colored MCS, one `SmaI` opening event that leaves
tiny yellow blunt vector ends, and two possible blunt ligation products.

Regenerate the hero SVG from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon template-validate \
  docs/examples/protocol_cartoon/pcr_blunt_vector_ligation_template.json

cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-template-svg \
  docs/examples/protocol_cartoon/pcr_blunt_vector_ligation_template.json \
  docs/figures/pgex_blunt_pcr_ligation_hero.svg
```

`pgex_blunt_pcr_ligation.workflow.json` is the corresponding replayable
off-line workflow behind that README section. It deliberately splits the same
cloning story across two stable local source files:

- `test_files/pGEX_3X.fa`
  - used as the linear PCR template so the shared `Pcr` route can derive the
    250 bp amplicon with primers `ATCCGCTCATGAGACAATAA` and
    `AAACGTTCTTCGGGGCGAAA`
- `test_files/pGEX-3X.gb`
  - used as the circular pGEX vector so `Digest` with `SmaI` yields one opened
    blunt backbone instead of two linear fragments

The workflow then runs blunt ligation with the PCR product and writes
`pgex_blunt_pcr_ligation_lineage.svg`, which shows the exact two orientation
products expected for a non-directional blunt insertion.

Regenerate that workflow figure from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/pgex_blunt_pcr_ligation.state.json \
  workflow @docs/figures/pgex_blunt_pcr_ligation.workflow.json
```

`gibson_single_insert_readme.plan.json` is the deterministic pGEX + insert demo
Gibson plan used to generate the README lineage figure. It is not hand-drawn;
it reuses the same tutorial destination/insert pair and the same `SmaI`
cutpoint (`941..941`) used in the Gibson specialist walkthrough.

`gibson_single_insert_lineage.svg` and `gibson_single_insert_lineage.png` are
the corresponding project-level lineage graph exports after applying that
single-insert Gibson plan. They show one Gibson operation with two inputs and
three concrete outputs (left primer, right primer, assembled product). This is
the same visible DALG/lineage graph users can export in the GUI with
`File -> Export DALG SVG...` or the graph-canvas context menu item
`Save Graph as SVG...`; it is not a screenshot.

Regenerate the lineage assets from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_readme_lineage.state.json \
  workflow @docs/examples/workflows/gibson_specialist_testing_baseline.json

cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_readme_lineage.state.json \
  gibson apply @docs/figures/gibson_single_insert_readme.plan.json

cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_readme_lineage.state.json \
  render-lineage-svg docs/figures/gibson_single_insert_lineage.svg

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/gibson_single_insert_lineage.svg \
  docs/figures/gibson_single_insert_lineage.png \
  --scale 2
```

`gibson_single_insert_rack_isometric.svg` is the physical-rack companion
figure for the same single-insert Gibson tutorial baseline. It is exported
from the saved linked rack draft, not drawn by hand, and shows how the
arrangement/vector/insert/product story projects into one pseudo-3D wet-lab
carrier using the shared `pipetting_pcr_tube_rack` physical template.

`gibson_single_insert_rack_hero.svg` is the README-facing presentation variant
derived from that same technical isometric export. It keeps the real rack
geometry and occupied positions, but removes most annotation and applies a
softer, more natural visual treatment for top-level documentation.

`gibson_single_insert_arrangement_gel.svg` is the matching README gel export
for the same single-insert Gibson baseline. For presentation, it intentionally
reorders the same saved samples as `insert -> vector -> product` so the insert
sits immediately beside the finer ladder and the vector/product comparison
reads more directly. It also uses an asymmetric ladder pair so the small insert
is easier to size while the vector/product bands still have a broader
reference:

- left ladder: `GeneRuler 100bp DNA Ladder Plus`
- right ladder: `Plasmid Factory 1kb DNA Ladder`

`gibson_single_insert_arrangement_gel_hero.svg` is the README-facing
presentation variant of that gel. It keeps the real lane/band geometry from the
shared serial-gel export, but simplifies the header and lane labels so the rack
and gel story can be read together at a glance.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_rack_hero.state.json \
  workflow @docs/examples/workflows/gibson_arrangements_baseline.json

cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_rack_hero.state.json \
  racks isometric-svg \
  rack-1 \
  docs/figures/gibson_single_insert_rack_isometric.svg \
  --template pipetting_pcr_tube_rack

python3 docs/figures/render_rack_isometric_hero.py \
  docs/figures/gibson_single_insert_rack_isometric.svg \
  docs/figures/gibson_single_insert_rack_hero.svg

cargo run --quiet --bin gentle_cli -- \
  --state /tmp/gibson_rack_hero.state.json \
  render-pool-gel-svg \
  - \
  docs/figures/gibson_single_insert_arrangement_gel.svg \
  --containers container-2,container-1,container-5 \
  --ladders "GeneRuler 100bp DNA Ladder Plus,Plasmid Factory 1kb DNA Ladder"

python3 docs/figures/render_serial_gel_hero.py \
  docs/figures/gibson_single_insert_arrangement_gel.svg \
  docs/figures/gibson_single_insert_arrangement_gel_hero.svg \
  "100 bp ladder" insert vector product "1 kb ladder"
```

`tp53_ensembl116_panel_source.gb` is a synthetic TP53 locus slice whose
feature geometry comes from the Ensembl 116 GRCh38 TP53 annotation. Its
sequence bases are placeholder `N`s because the isoform-architecture render
consumes transcript/CDS coordinates plus `cds_ranges_1based`, not genomic base
content.

`tp53_isoform_architecture.svg` is rendered from that source plus the curated
panel resource `assets/panels/tp53_isoforms_v1.json`.

Regenerate it from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  workflow @docs/figures/tp53_isoform_architecture.workflow.json
```

`tp73_cdna_genomic_dotplot.workflow.json` renders an offline TP73 cDNA-vs-
genomic dotplot from the local `test_files/tp73.ncbi.gb` fixture. The workflow
loads the TP73 genomic locus, derives transcript feature `n-100`
(`feature_ids=[99]`, `NM_001126241.3`), computes a pair-forward dotplot against
the same genomic source, and writes the auditable SVG intermediate
`tp73_cdna_genomic_dotplot.svg`.

`tp73_cdna_genomic_dotplot.png` is the README-friendly raster derivative. It is
rendered from the SVG with `resvg` through `gentle_examples_docs svg-png` while
dropping the export metadata/title text and keeping the basepair axis labels.

Regenerate the TP73 dotplot assets from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/tp73_readme_dotplot.state.json \
  workflow @docs/figures/tp73_cdna_genomic_dotplot.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/tp73_cdna_genomic_dotplot.svg \
  docs/figures/tp73_cdna_genomic_dotplot.png \
  --drop-dotplot-metadata
```

`tp73_multi_isoform_overlay_dotplot.workflow.json` renders an offline
shared-reference TP73 multi-isoform dotplot overlay from the same local
`test_files/tp73.ncbi.gb` fixture. The workflow loads the TP73 genomic locus,
derives four curated transcript variants (`NM_005427.4`, `NM_001126240.3`,
`NM_001126241.3`, `NM_001204189.2`), computes one
`ComputeDotplotOverlay` payload against the shared genomic reference, and
writes the combined SVG `tp73_multi_isoform_overlay_dotplot.svg`. The stored
query specs now also include their original TP73 transcript feature ids so the
same overlay payload remains eligible for shared-exon-anchor exports.

`tp73_multi_isoform_overlay_dotplot.png` is the README-friendly raster
derivative. It is rendered from the SVG with `resvg` through
`gentle_examples_docs svg-png` while dropping the dotplot export metadata/title
text and keeping the basepair axis labels.

Regenerate the TP73 multi-isoform overlay dotplot assets from the repository
root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/tp73_overlay_dotplot.state.json \
  workflow @docs/figures/tp73_multi_isoform_overlay_dotplot.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/tp73_multi_isoform_overlay_dotplot.svg \
  docs/figures/tp73_multi_isoform_overlay_dotplot.png \
  --drop-dotplot-metadata
```

`tp73_multi_isoform_anchor_dotplot.workflow.json` renders the README-facing
TP73 shared-exon-anchor variant of that overlay. It uses the same four curated
TP73 transcripts, restricts the genomic reference span to the shared downstream
region (`54001..83686` local TP73 slice coordinates), and renders the overlay
with `overlay_x_axis_mode=shared_exon_anchor` around the common exon
`55034..55276`. The output SVG is
`tp73_multi_isoform_anchor_dotplot.svg`.

`tp73_multi_isoform_anchor_dotplot.png` is the README-friendly raster
derivative of that anchored export.

Regenerate the TP73 shared-exon-anchor dotplot assets from the repository root
with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/tp73_anchor_readme.state.json \
  workflow @docs/figures/tp73_multi_isoform_anchor_dotplot.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/tp73_multi_isoform_anchor_dotplot.svg \
  docs/figures/tp73_multi_isoform_anchor_dotplot.png \
  --drop-dotplot-metadata
```

`tp73_upstream_tfbs_score_tracks.workflow.json` renders a README-facing stacked
TFBS/PSSM score-track figure for a TP73 promoter-proximal window. The workflow
loads the local `test_files/tp73.ncbi.gb` locus, keeps the continuous score
tracks for `TP53`, `TP63`, `TP73`, `PATZ1`, `SP1`, `BACH2`, and `REST` over the
internal TP73 transcript-start neighborhood `15564..16764` (`1000 bp`
upstream plus `200 bp` into the transcribed region around
`XM_047429524.1`), and writes the shared SVG
`tp73_upstream_tfbs_score_tracks.svg` through `RenderTfbsScoreTracksSvg`.

This figure intentionally uses `clip_negative=false`. In this local TP73 slice,
`SP1` is the clearest positive-support factor, while the remaining requested
motifs are still useful to show as continuous background scores instead of
flattening them to zero. When a transcription start falls inside or directly at
the edge of the plotted span, the renderer also adds a short hooked-arrow TSS
cue above the score tracks.

`tp73_upstream_tfbs_score_tracks.png` is the README-friendly raster derivative
of that shared export.

Regenerate the TP73 TFBS score-track assets from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/tp73_upstream_tfbs_score_tracks.state.json \
  workflow @docs/figures/tp73_upstream_tfbs_score_tracks.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/tp73_upstream_tfbs_score_tracks.svg \
  docs/figures/tp73_upstream_tfbs_score_tracks.png
```

`toy_shared_exon_anchor_source.gb` is a hand-crafted 150 bp synthetic teaching
locus with four short exons and three transcripts. All three transcripts share
exon `85..108`, but one starts later (`delta 5'`) and one ends earlier
(`delta 3'`), so the anchored x-axis behavior can be read immediately without
TP73-scale biological density. It is intended as the compact companion figure
for README/ClawBio materials when the visual grammar matters more than the
real-locus complexity.

`toy_shared_exon_anchor_dotplot.workflow.json` loads that source, derives the
three transcript sequences from mRNA features `n-7..n-9`, computes one
reference-centered overlay against the same toy genomic source, and renders the
shared-exon-anchor view around exon `85..108`.

`toy_shared_exon_anchor_dotplot.svg` and
`toy_shared_exon_anchor_dotplot.png` are the corresponding export artifacts.

Regenerate the toy shared-exon-anchor teaching figure from the repository root
with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/toy_anchor.state.json \
  workflow @docs/figures/toy_shared_exon_anchor_dotplot.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/toy_shared_exon_anchor_dotplot.svg \
  docs/figures/toy_shared_exon_anchor_dotplot.png \
  --drop-dotplot-metadata
```

`p53_family_tp73_longest.fasta`, `p53_family_tp63_canonical.fasta`, and
`p53_family_tp53_canonical.fasta` are a small curated cross-family demo set:

- `p53_family_tp73_longest.fasta`
  - copied from record 1 of
    `test_files/fixtures/mapping/ensembl_human_tp73_all.fasta`
    (`ENST00000378295.9`, protein-coding, longest bundled human TP73 transcript)
- `p53_family_tp63_canonical.fasta`
  - copied from record 2 of
    `test_files/fixtures/mapping/ensembl_human_tp63_all.fasta`
    (`ENST00000264731.8`, protein-coding TP63 transcript)
- `p53_family_tp53_canonical.fasta`
  - copied from record 1 of
    `test_files/fixtures/mapping/ensembl_human_tp53_all.fasta`
    (`ENST00000413465.6`, protein-coding TP53 transcript)

The copied headers retain their Ensembl transcript provenance and add the exact
source fixture path plus selected record index so the demo stays reviewable and
re-creatable from the committed mapping fixture pack.

`p53_family_query_anchor_dotplot.workflow.json` renders the first cross-family
anchored overlay demo. It loads those three curated transcripts, keeps TP73 as
the shared reference/Y-axis, restricts the reference span to the conserved
TP73 core region (`701..1700` in local cDNA coordinates), and overlays TP63
plus TP53 using explicit per-query anchors at the conserved motif
`CATGTGTAACAG`. The workflow exercises the new `query_anchor_bp` overlay
x-axis mode rather than the same-locus `shared_exon_anchor` mode.

`p53_family_query_anchor_dotplot.svg` and
`p53_family_query_anchor_dotplot.png` are the corresponding cross-family demo
artifacts.

Regenerate the anchored p53-family demo from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/p53_family_anchor.state.json \
  workflow @docs/figures/p53_family_query_anchor_dotplot.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/p53_family_query_anchor_dotplot.svg \
  docs/figures/p53_family_query_anchor_dotplot.png \
  --drop-dotplot-metadata
```

`vkorc1_rs9923231_context_map.workflow.json` is the first concrete
ClawBio-facing genomic-context asset. It resolves dbSNP `rs9923231` against
prepared `Human GRCh38 Ensembl 116`, extracts `+/- 3000 bp`, keeps full locus
annotation, overlays a focused JASPAR motif panel (`SP1`, `CTCF`, `HNF4A`),
hides non-context tracks such as restriction enzymes/GC/ORFs, and writes the
map-first SVG `vkorc1_rs9923231_context_map.svg`.

For legibility, the exported figure is an mRNA-focused asymmetric crop of that
larger retrieval rather than a full 6 kb viewport:

- right edge: `500 bp` to the right of `rs9923231`
- left edge: `500 bp` into the first gene body encountered on the left side of
  the SNP in this prepared Ensembl 116 context
- transcript (`mRNA`) features remain visible, while gene-level lanes are
  hidden to reduce redundancy in the exported figure
- export-time label pruning now prefers gene-style names over raw accession
  labels where possible and compacts nearby repeats so the context figure reads
  as a communication asset instead of an annotation dump
- the SNP itself is rendered as a baseline marker on the DNA rather than as a
  floating feature block
- the TFBS overlay is intentionally figure-focused rather than exhaustive: it
  uses a small JASPAR PSSM panel with stricter per-motif thresholds and
  quantile-based display filtering so the promoter context stays legible

So the retrieval stays tutorial/replay-friendly while the figure itself stays
focused on the local `VKORC1` regulatory neighborhood.

The shared `RenderSequenceSvg` path now honors that stored linear viewport, so
the figure workflow's asymmetric crop is applied in the exported SVG itself
rather than only in the workflow intent.

Preflight the target GENtle instance first so the workflow uses the exact
catalog/cache you intend:

```sh
cargo run --quiet --bin gentle_cli -- \
  genomes status "Human GRCh38 Ensembl 116" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes
```

If `prepared` is `false`, prepare the selected cache first:

```sh
cargo run --quiet --bin gentle_cli -- \
  genomes prepare "Human GRCh38 Ensembl 116" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes
```

Then regenerate the context-map asset from the repository root with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/vkorc1_rs9923231_context.state.json \
  workflow @docs/figures/vkorc1_rs9923231_context_map.workflow.json
```

`vkorc1_rs9923231_luciferase_hero.svg` is the complementary community-facing
figure for the same story. It is now a hybrid asset:

- the left panel is a hand-authored editorial explainer for the upstream
  promoter-selection logic
- the right panel is a real circular DNA-window export from
  `vkorc1_rs9923231_luciferase_construct.svg`, now embedded inline so the hero
  SVG stays self-contained in viewers that do not resolve external image refs

Its honest job is not to claim that the direct molecular effect of
`rs9923231` is already known. Its job is to show one plausible
question-to-assay bridge:

- the biological question can be framed as: how do we learn more about the
  direct molecular effect of the warfarin-response SNP `rs9923231`?
- because `rs9923231` lies upstream of `VKORC1`, one plausible hypothesis is a
  direct effect on expression of that gene
- one plausible ClawBio/agent-suggested answer is therefore to compare matched
  `VKORC1` promoter->indicator-gene reporters
- the study fragment should be a TSS-centered `VKORC1` promoter window that
  keeps `rs9923231` inside the insert rather than on its edge
- that fragment is then placed upstream of luciferase for the follow-up assay

ClawBio's role in this story is not to claim that it already knows the direct
molecular effect. Its value is that it can help answer a complicated question
with a tractable experimental next step:

- ClawBio carries the pharmacogenomic alert context and helps propose the
  promoter->luciferase follow-up
- GENtle provides the explicit deterministic operations, figure exports,
  tutorials, workflows, and online help that let that assay design be reviewed
  and executed rather than hidden inside free-form prompting
- put differently: GENtle describes its own capabilities clearly enough that
  ClawBio can bridge into them instead of improvising the cloning logic

For slide decks and similar communication assets, prefer the reduced left
cartoon derived from the hero figure together with the standalone construct
export:

- left panel: `docs/figures/vkorc1_rs9923231_luciferase_left_cartoon.svg`
  (and matching `.png` / `.pdf`), extracted from the hand-authored left hero
  panel so the locus cartoon can stand alone without the surrounding hero-card
  frame
- right panel: `docs/figures/vkorc1_rs9923231_luciferase_construct.png`
  or `.svg` / `.pdf` directly, without any surrounding hero card
  The standalone construct export now carries larger `Insert` and
  `Luciferase` callouts so it reads more clearly in slide contexts.

Keep the raw source figures as supporting evidence assets:

- `docs/figures/vkorc1_rs9923231_context_map.png`
- `docs/figures/vkorc1_rs9923231_luciferase_construct.png`

The raw context export is still useful for documentation and provenance. For
talks, the extracted left cartoon reads more clearly, while the right
standalone construct export keeps the biology explicit without the extra
hero-frame.

Keep the existing TP73 tutorial for GUI/CLI parity walkthroughs and detailed
stepwise testing.

For ongoing work on a ClawBio integration, also treat this asset as another
hero image:

- SVG source:
  `docs/figures/vkorc1_rs9923231_luciferase_hero.svg`
- raster derivative for slides/docs:
  `docs/figures/vkorc1_rs9923231_luciferase_hero.png`
- role in the story:
  show the bridge from a complicated pharmacogenomic question into a
  GENtle-guided follow-up assay design, with the left panel explaining the
  reverse-strand `VKORC1` promoter interval and the right panel showing the
  luciferase construct context

`vkorc1_rs9923231_luciferase_construct.gb` is the synthetic source used to get
that real circular export style without requiring a live luciferase-cloning
workflow during figure generation. It is intentionally illustrative rather than
wet-lab-ready. Regenerate the construct-map assets with:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/vkorc1_rs9923231_luciferase_construct.state.json \
  workflow @docs/figures/vkorc1_rs9923231_luciferase_construct.workflow.json

cargo run --quiet --bin gentle_examples_docs -- \
  svg-png \
  docs/figures/vkorc1_rs9923231_luciferase_construct.svg \
  docs/figures/vkorc1_rs9923231_luciferase_construct.png
```

The current construct source deliberately keeps `rs9923231` inside the cloned
insert rather than on its exact edge, and flips `bla` onto the opposite strand
so the circular export can demonstrate the transcription-start/direction
markers clearly.

## Planned VKORC1/rs9923231 PGx-Alert-to-Construct Showcase

The next README/community-facing showcase should be a
**VKORC1/rs9923231 question-to-assay-design** story aligned with the
ClawBio/OpenClaw handoff model:

- ClawBio starts from the complicated warfarin/`rs9923231` question and helps
  bridge it to one plausible follow-up assay idea
- GENtle retrieves the relevant `VKORC1` context, derives allele-specific
  upstream regulatory fragments, and plans luciferase reporter constructs
- the run yields both explanation figures and a reproducibility bundle
- the baseline claim is not that the molecular mechanism is already known, but
  that allele-specific promoter activity is a reasonable first functional
  follow-up to test; any warfarin treatment arm is a later extension
- one later extension can borrow the warfarin/Factor X theme more directly:
  keep a `factor Xa` recognition site in the construct, but place it between
  reporter/indicator modules rather than in a purification-first GST layout, so
  cleavage efficiency becomes a ratio readout instead of a protein-extraction
  problem
- if that extension is taken into human cells or an adenoviral delivery setup,
  the same assay logic should move into an appropriate mammalian/viral
  expression backbone rather than stay in `pGEX`
- the opening visual should look like a real ClawBio handoff:
  a hard question first, then a concrete assay-design explanation

The intended asset set is:

1. `vkorc1_rs9923231_pgx_alert_panel.*`
   - upstream warfarin/PGx alert panel from the ClawBio side
   - likely assembled from the ClawBio skill/report layer rather than from a
     pure GENtle render
   - should carry the actual alert semantics:
     drug, genotype call, `VKORC1`, and `rs9923231`
2. `vkorc1_rs9923231_context_map.*`
   - first concrete GENtle-side asset is now implemented as
     `docs/figures/vkorc1_rs9923231_context_map.workflow.json` plus
     `docs/figures/vkorc1_rs9923231_context_map.svg`
   - shows the local `rs9923231` / `VKORC1` / `LOC124903680` neighborhood
   - should make the assembly/build and reverse-strand promoter orientation
     explicit
3. `vkorc1_rs9923231_luciferase_hero.svg`
   - now implemented as the preferred luciferase showcase opener
   - hybrid explainer figure showing the exact reverse-strand promoter
     fragment taken from the `VKORC1 5'` boundary up to `rs9923231`
   - right panel now uses a real circular DNA-window export rather than a
     pure construct cartoon
   - still explains the reporter design, with the reference and variant
     constructs differing only at the SNP allele
   - intended to replace TP73 as the community-facing luciferase opener while
     keeping TP73 as the parity/tutorial walkthrough
4. `vkorc1_rs9923231_luciferase_construct.svg`
   - now implemented as the actual circular DNA-window export used inside the
     right panel of the hero figure
   - generated from `docs/figures/vkorc1_rs9923231_luciferase_construct.gb`
     through `docs/figures/vkorc1_rs9923231_luciferase_construct.workflow.json`
   - intentionally illustrative: it provides real GENtle map styling for the
     communication asset before the full luciferase-cloning workflow is wired
     end to end
5. `vkorc1_rs9923231_luciferase_protocol_cartoon.svg`
   - planned VKORC1-targeted promoter->luciferase cloning mechanism strip
   - should come from the same protocol-cartoon engine family already used for
     Gibson README figures, but presented as a reporter-assay build rather than
     as a generic cloning strip
6. `vkorc1_rs9923231_luciferase_lineage.svg`
   - planned lineage/provenance export showing allele-specific inserts,
     primers, and assembled reporter construct(s) from the same project state
7. `vkorc1_rs9923231_clawbio_bundle_panel.*`
   - compact visual summary of `report.md`, `result.json`, and reproducibility
     bundle contents from the ClawBio/OpenClaw wrapper

The existing repository assets that should anchor the first implementation are:

- `docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json`
  - current VKORC1/warfarin promoter->luciferase planning skeleton already
    covering dbSNP-driven locus retrieval, TSS-centered fragment extraction,
    vector import, assembly preview, and junction-PCR reporting
- `docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json`
  - older promoter->luciferase planning skeleton retained as historical
    precedent rather than the canonical community-facing tutorial path
- `docs/examples/workflows/tp63_extend_anchor_online.json`
  - anchored-region extension baseline for genome-context expansion
- `docs/examples/workflows/digest_ligation_extract_region_minimal.json`
  - cloning/export baseline for deterministic construct materialization
- `assets/helper_genomes.json`
  - includes the luciferase destination helper entry used for existing reporter
    planning examples

Recommended execution order for implementing the full showcase:

1. Add one genomic-context figure showing `rs9923231` in the local `VKORC1` /
   `LOC124903680` neighborhood on the correct assembly/build.
2. Add one allele-paired promoter-fragment figure that answers:
   "what exact regulatory sequence are we taking forward into the reporter
   assay?"
3. Add one `VKORC1`-targeted luciferase construct workflow that produces:
   - a promoter->luciferase mechanism/cartoon export
   - a lineage/provenance export
   - primer/qPCR reports for the reporter construct
4. Add one ClawBio-wrapper demo/report snapshot tying the whole run back to the
   pharmacogenomics story and reproducibility bundle.

Current status:

- the first shared-engine handoff artifact now exists:
  `docs/examples/workflows/vkorc1_rs9923231_context_map_online.json`
- the first figure-workflow asset now exists:
  `docs/figures/vkorc1_rs9923231_context_map.workflow.json`
- the corresponding genomic-context SVG is now part of the repository:
  `docs/figures/vkorc1_rs9923231_context_map.svg`
- the community-facing luciferase hero figure now also exists:
  `docs/figures/vkorc1_rs9923231_luciferase_hero.svg`
- the hero figure now also has a real circular construct-map companion export:
  `docs/figures/vkorc1_rs9923231_luciferase_construct.svg`
- protocol-cartoon exports, lineage/provenance exports, and the ClawBio bundle
  panel remain planned
