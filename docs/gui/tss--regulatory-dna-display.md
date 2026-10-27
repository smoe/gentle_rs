# TSS / Regulatory DNA Display

<!-- manual-counterparts:start -->
Related CLI reference: [Promoter Design](../cli/promoters.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Open a GENtle **annotated TSS EMBL or GenBank** sequence in the DNA viewer,
then choose **TSS / Regulatory** beside Standard map. The same sequence stays
loaded; switching display does not extract a new molecule or change its bases.
The title becomes gene, chromosome and genomic TSS position after validation.
The file's original name and provenance remain unchanged.

The native display separates transcript exon/coding context, each supplied
CUT&RUN/chromatin source, retained motif peaks per matrix and score kind, and
other annotations. Category checkboxes and a lane-name/source-ID filter limit
clutter. The ordinary zoom/pan controls apply to its shared local axis. Rulers
show local bases, genomic bases and signed TSS-relative offsets; a negative
genomic strand decreases genomic coordinates without reversing the local DNA.
Side summaries show interval counts and displayed coordinate envelopes (not
continuous coverage). Hover retains original feature/source notes. Click a
feature to select its exact local DNA span, then **Inspect in standard DNA map**
to see that annotation and sequence. For report-only hits or curve footprints,
**Show DNA span in standard map** selects the bases and clears any previous
annotation highlight; it does not create an annotation. Wrapped/clipped original spans remain in
details; unsupported compound locations are explicitly reported, not flattened.

Signal intervals preserve gaps and supplied source scale; no-interval and
not-prepared sources are not measured zeros. Motif triangles above/below the
baseline distinguish local plus/minus strands, not positive/negative values.
Their height uses the original score kind, never relabels tail probabilities
as LLR bits, and hides values below zero without deleting data. These annotation lanes are stored
peaks, **not full continuous score traces or experimentally measured binding**.
Display clipping, filtering and selection never rescore motifs.

Recognition is cached; hash/geometry validation and grouping run off the paint
path. Only visible rows are painted, with an explicit 10,000-feature drawing
limit per lane/viewport; zoom in if that limit is reached. The source data are
not truncated. Editing/replacing the sequence invalidates the presentation.
Mismatched bases, TSS markers, or motif-strand metadata fail closed. Sequence
hash validation is not external reference authentication or a receipt audit.

### Attach the Quantitative Report

Choose **Load TSS profile report...** and select `report.json` from the TSS
SVG/PDF export bundle belonging to this window. Loading and validation run in
the background, independently of the TFBS panel and its caches. No database
query or scoring is performed. The file is bounded to 256 MiB and passes the
same report validation as export (including matrix and imported-report hashes,
array geometry and the ten-million-score-cell limit). Genome ID, assembly,
annotation release, promoter ID, exact window/TSS geometry and sequence digest
must match. Equal bases at another locus are not sufficient. Failed attachment
leaves the previous view intact; replacing/editing the sequence discards the
attachment. **Detach report** removes report-only lanes, retaining annotations
and any independently computed local curves.

- **Report TFBS curves** shows complete stored forward/reverse arrays at motif
  window starts. Local `+` is solid blue; local `-` is dashed rose. Both use the
  same increasing local axis, also on negative genomic strands. Hover reports
  raw scores and both coordinate/strand conventions. Clicking an evaluable
  start selects its full motif footprint without creating a DNA annotation.
- Each curve has numeric Y ticks and its own explicit score kind. Scaling and
  negative-score clipping follow the saved panel policy, not current TFBS-panel
  settings or an export-only scale override. A background-tail score is not
  relabelled LLR bits. Amber areas and line gaps mean unavailable windows:
  the upper half represents local `+`, the lower half local `-`, relative to
  the displayed DNA, including on negative genomic strands. The
  grey terminal areas cannot hold a full motif window. A genuine zero remains
  an evaluable value. An empty/all-zero range uses a labelled display fallback.
- **DuckDB peaks** shows the report's attached sparse-query results, not a live
  query. Amber triangles encode footprint and local strand; height uses the
  original source/report/matrix score range, including signed raw values, never
  the local-curve scale. Side summaries retain coverage and truncation states;
  hover includes original genomic spans, retention/density policy and hashes.
  A valid empty subset with no declared score mode says **score units unavailable**,
  rather than inventing units or displaying an empty label.
  Missing or partially queried evidence does not establish absence of binding.

The original annotated structure, CUT&RUN/chromatin and stored-peak lanes stay
separate. Report annotations are not silently substituted for edited flat-file
features. Attachment is session-local and not a receipt audit or independent
reference authentication. Display preparation is cached; visible curve drawing
is limited to 10,000 window starts per row. Above that budget the view asks for
zooming rather than silently downsampling away narrow peaks or gaps.

### Compute Local TSS Scores

Open **Local scoring**, enter exact matrix accessions (comma/space-separated),
choose `score_kind` and clipping, then **Compute local scores**. These settings,
worker, cancellation and cached result belong to this TSS view, not the separate
TFBS panel. Scoring uses the shared engine's linear `InlineSequence` target,
bound to the displayed, transcript-oriented DNA. It never fetches sequence,
queries DuckDB, reconstructs an old report, or modifies sequence annotations.

**Locally computed curves** is a separate lane class and visibility toggle.
Local `+` curves are solid violet, distinct from blue report curves and green
stored peaks; local `-` curves remain dashed rose. The native view and its SVG
export share the colour and local-strand legend, including unavailable halves.
Each lane shows its actual score kind, exact matrix hash and evaluability count;
its own range does not normalize imported/report scores. Zero is a scored
value, while ambiguous/unscorable windows remain gaps via the engine validity
masks. These are model predictions, not experimental binding or a prediction
that a mutation will reduce luciferase activity. Stored report curves retain
their original settings even when the local settings change.

Admission requires 1..32 unique full-PFM IDs, at most 64 columns per matrix,
50,000 window bases and 1,000,000 matrix/base combinations. Aliases, `ALL` and
consensus fallback are refused here; existing headless routes remain available.
Resolution/calibration/scoring run on a background worker with stage/matrix
progress, not an ETA. **Cancel scoring** keeps the last complete result.
Replaced documents or changed requests cannot accept late results. Panning
does not recompute; changing settings leaves previous curves explicitly labelled
until an explicit new run. Clearing local scores preserves attached reports,
and detaching reports preserves local scores. Jobs never share a result cache.
The cache receipt binds reference/geometry/sequence, effective local request,
actual matrix hashes, background policy and scorer; a fresh run always resolves
the current registry anew. A matrix mismatch between admission and scoring
rejects the result.

### Agent and Headless Routes

Every step of the native view is reachable without the mouse. `ui open|focus
tss-view --report REPORT_JSON` attaches a report, and `ui open|focus tss-view
--local-score ACCESSION[,ACCESSION...] [--score-kind KIND] [--keep-negative]`
computes local curves; the viewer adopts those settings into its own controls so
the displayed panel always describes the curves that were actually computed.
Both need a GUI host, and headless shells record `applied=false`.

`promoters tss-view-svg SEQ_ID OUTPUT.svg` needs no GUI host at all: it decodes
the annotated window, optionally attaches a report with `--report`, optionally
computes local curves with `--motif`, and uses the native view's renderer for
those attachments. It reports the lane count, exported local span and SVG SHA-256, and
refuses a sequence without GENtle annotated-TSS metadata. `--score-kind` is not
inherited from the viewer; repeat it per export.

With `--motif --span START..END`, the headless operation scores only complete
motif footprints inside that local, 1-based inclusive span, on both strands.
It excludes footprints crossing either edge. Counts, local axis ranges and
empirical quantiles describe the scanned span, whose exact coordinates and
sequence hash are retained alongside the full annotated-window binding.
The denominator is `2 * max(0, span_length - motif_length + 1)`; ambiguous
windows are unavailable, not zero, and the trailing grey band means no complete
footprint **inside the scored span**. Admission still applies to the full
annotated window, even if the requested span is small.

This differs from panning and **Export View SVG**, which do not rescore and
keep the existing full-window count, calibration and scale. Attached `--report`
curves likewise remain untouched. A cropped existing plot and a newly scored
subspan therefore need not have the same edge coverage, scale or empirical
quantiles. Same-window, same-scoring-scope GUI/headless parity remains required.

Each full-curve lane's side summary lists up to three highest raw-score window
starts in the displayed span, with signed TSS offsets and local strands. Ranking
is strictly within that matrix lane, using raw scores even when negative values
are clipped for drawing. Unavailable starts are excluded and valid zero is not.
Ties use start coordinate, then local `+` before `-`; these positions are not a
binding-likelihood ranking and do not rank one factor against another.

For all validated collection members, the headless/agent route is
`promoters tss-view-svg --collection COLLECTION_ID OUTPUT_DIR` with the same
attachment/scoring/span flags. At most 32 members and 32 MiB per complete bundle
are allowed; every member is checked before scoring. Ordered SVG pages,
`index.html` links and a typed `receipt.json` are staged before publishing to a
new directory. Per-window/per-lane scale policy preserves supplied report
scaling, with no new cross-window calibration. This command does not change
approvals, select more starts or open a GUI. The single-window route also returns
the typed `tss_view_svg_export` artifact receipt.

See the [TP73 ΔNp73 factor-curve walkthrough](../tutorial/08-17_tp73_dnp73_factor_curves.md)
for all three routes on a real locus.

### Export the Native TSS View

**Export View SVG** exports the current horizontal span and enabled/filter-matched
TSS lanes, including lanes vertically scrolled out of sight. It exports the
TSS evidence plot, not a pixel screenshot or the separate sequence-text panel.
Wide-context and A3 variants use the existing expanded-span/print profiles.
Rendering and atomic file publication run in a separate worker on the selected
immutable presentation; later navigation does not alter that export.

The headless renderer preserves local/genomic/TSS coordinates on both strands,
separate source scales, stored peak footprints, missing-window gaps, terminal
unavailable regions, hover provenance and coverage/truncation warnings. SVG
metadata binds the presentation digest, sequence, reference, attached report
file hash, local-scoring provenance/report hash and exact span/lane options.
This is not a fresh scientific analysis or a full receipt audit. It never
substitutes the standard map. A loading or
stale presentation cannot be exported; excessive work is refused (10,000 items
per lane, 100,000 overall, 128 lanes and 32 MiB), with a request to narrow the
span or lane selection, rather than silently dropping evidence.

SVG score-window hover annotations are capped at 128 per lane and local strand,
retaining strongest and evenly sampled valid starts. Each lane discloses
emitted/total and omitted-title counts. Unsampled starts lose their own
coordinate/raw-score/footprint tooltip; every curve point and validity band
remains. Native GUI hovering is not capped. The same top-three within-lane
raw-score summary appears in the exported SVG, outside its clipped plot area.

GENtle's annotated TSS export grammar is still required, including existing
hashed record names. FASTA alone, arbitrary locus features and generic prose
labels cannot establish its TSS geometry. New DuckDB
queries from imported windows and reporter rows remain
follow-ups. Existing receipt-bound TSS report exports remain the quantitative
publication path; a native-view SVG records the selected presentation instead.
The existing single-sequence file loader selects the last EMBL/GenBank record;
multi-record import/selection is a separate follow-up, not provided by this view.

Shared GUI/agent navigation: focus the intended DNA window, then use
`ui open tss-view` or `ui focus tss-view`; `ui close tss-view` restores Standard
map. No arbitrary project sequence is selected when the active DNA context is
missing. Headless execution records the intent with `applied=false`.
To attach the same report as the native file-picker path, use
`ui open tss-view --report REPORT_JSON` (or `ui focus ...`). The active annotated
TSS viewer performs the normal binding validation asynchronously; the intent
does not select a sequence implicitly, rescore motifs or query DuckDB. Tutorial
[08.16](../tutorial/08-16_tss_regulatory_view_gui.md) demonstrates GUI, Agent Assistant
and inner-agent wording side by side. Agent Assistant preserves the explicitly
focused DNA viewer as its launch subject while it has focus; stale or closed
subjects are rejected rather than replaced from project inventory.

For an annotated anchored locus, **TFBS scan > Transcript starts / TSS windows...**
now provides a background inventory, explicit row selection and approval, then
collection-targeted opening. Truncated or uncertain 5-prime starts and unlinked
transcripts appear as nonselectable diagnostic rows, not candidate TSSs; unlinked
records are labelled locus-level unassigned rather than exclusions from the gene.
**Inspect stored collection** validates the entered collection ID in the
background and shows a concise **Validated collection** summary. Expand
**Technical details** for the member table, full fingerprint on hover and
copyable JSON. Stale or
legacy metadata produces the shared engine diagnostic. **Forget registry
entry...** requires a second, ID-bound confirmation and retains all sequences,
windows and lineage; use a new ID for re-derivation rather than overwriting them.
Feature
cache refresh does not mark an intact collection edited; changed bases,
annotations or genomic anchors still prevent reuse. `ui open|focus|close tss-view --collection ID`
targets that persisted collection independently of the active viewer. See the
[TSS workspace walkthrough and limitations](../tss_workspace.md).

Collection-window requests originating in a child DNA viewer wake the root
application poller before validation and window creation continue. This keeps
the asynchronous request bounded while ensuring an otherwise idle root does
not leave it queued indefinitely.
New sequence windows use a bounded 48-point cascade (six positions); existing
windows retain their positions. Collection opening reuses pending/open members
and requests focus for the first member after new windows are registered.
Use the native **Window** menu to raise another named member. Native focus
remains subject to the window manager; a queued request is not proof of focus.

**Refresh collections** lists stored collection IDs, gene queries, source loci
and window counts. Select a row, then explicitly **Inspect stored collection**;
the browser's **not checked** label is never a validity claim. Invalid/legacy
records remain visible. Selecting another ID clears old validation and pending
forget confirmation; project edits invalidate the displayed validation. The
workspace scrolls on smaller screens. The [synthetic TSS tutorial](../tutorial/08-15_tss_collection_gui.md)
has ten learner steps from preview through explicit approval, validation and
opening members. Repeated opening has focused application regressions;
Forget/Undo, save/restart and stale-member editing remain separate contributor
checks. Its retained Linux/X11 evidence names the exact source revision and
does not establish merged-candidate, Windows/macOS or live TP73 acceptance.
Prepared GenBank/XML references do not supply a tabular transcript index:
extraction reports unavailable transcript enrichment, not zero transcripts.
The workspace can still use transcript features already imported into the
loaded locus; see the walkthrough's annotation-source distinction.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: TFBS Annotation (Engine Ops)](tfbs-annotation-engine-ops.md) | [Next: RNA Read Panel](rna-read-panel.md)
<!-- manual-page-nav:end -->
