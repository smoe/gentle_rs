# Evidence-Guided Reporter Boundaries

`promoters fragment-candidates` proposes inserts from a bound gene-locus report.
It does not replace fixed-window exports, save regions, alter approved plans,
design primers, or create constructs. Use it before the exact-ROI contrast
planner when the question is **which sequence should we test?**

## First Proposal

Load the annotated source locus into the project and export its composed
`gentle.gene_locus_evidence_display.v1` report (a reporter-comparison envelope
containing that report is also accepted). Its sequence/anchor binding must
match the loaded locus. Specify the exact file SHA-256, panel ID, annotation
release and transcript IDs, not a gene-symbol guess:

```json
{
  "schema": "gentle.reporter_fragment_selection_request.v1",
  "locus": {
    "path": "locus.json",
    "sha256": "sha256:REPLACE_WITH_EXACT_FILE_SHA256",
    "panel_id": "REPLACE_WITH_PANEL_ID",
    "annotation_release": "REPLACE_WITH_RECORDED_RELEASE"
  },
  "anchors": [{ "transcript_id": "REPLACE_WITH_EXACT_TRANSCRIPT_ID" }]
}
```

```bash
cargo run --quiet --bin gentle_cli -- --state STATE.json \
  promoters fragment-candidates @selection_request.json --path candidates.json
```

These are placeholders for your own bound inputs, not an executable biological
example. The command uses the existing shell parser; the identical request is
available to GUI Shell, MCP `op`, workflows and operation-based scripting as
`PlanEvidenceGuidedFragmentCandidates { request, path? }`. The optional file
contains the exact report returned in `OpResult.reporter_fragment_selection`.

## What Determines The Insert

The default search envelope is TSS -700/+300, **1,001 bp including the TSS**.
It is a discovery envelope, not the insert. Available Ensembl regulatory core
intervals, already-scored TFBS site footprints and explicitly supplied,
provenance-bound called peaks may seed proposals. Raw
CUT&RUN coverage alone does not constitute a peak or a regulatory boundary.
To limit alternatives, supply `seed_evidence_ids` from the returned inventory.
`required_evidence_ids` specifies context that no eligible alternative may lose.

For `purpose: "endogenous_promoter"` (default), a promoter annotation covering
the TSS is retained. Without one, the configured provisional -100/+50 TSS
context is retained and explicitly labelled unproven. For a response element
to be tested with a separately supplied minimal promoter, choose
`purpose: "response_element"`; it does not silently inherit that TSS context.
Current transcript exons/strand are checked against the report. The shared
CDS-start/5'-UTR audit reports start-codon, upstream-ATG/uORF and UTR-intron
context when assessable; absent CDS evidence is not a clean translation audit.
Shared TSSs retain transcript membership, while opposite-strand TSSs are flagged
rather than merged into a single promoter identity.

The initial insert contains the seed and required context, with 20 bp padding.
The preferred 700 bp length is **soft**: padding may be removed, but required
context is not discarded to fit. Longer alternatives stay visible. The default
hard maximum is 5,000 bp. All values are configurable and recorded in `request`.
The source must cover the proposed geometry; no sequence is downloaded or
silently invented beyond the loaded locus.

Ranking is lexicographic and exposed in `ranking`: preserve required context,
declared descriptive support, retained annotation, retained model site, retained
called peak, fewer
bisected biological features, compactness, length, stable ID. It is not a
weighted biological confidence score. Signal-bin boundaries do not influence
the bisection rank. `blockers` controls eligibility; a high-ranked row with a
hard-limit blocker must not be treated as an accepted insert.
These are presence flags, not independent votes. Several overlapping peaks do
not contribute more votes, and a called peak is not independent of its coverage.

## Revising Either End

Add an `adjustments` entry to the same request, using the returned seed ID:

```json
{
  "adjustment_id": "inspect-additional-upstream-context",
  "transcript_id": "REPLACE_WITH_EXACT_TRANSCRIPT_ID",
  "seed_id": "REPLACE_WITH_SEED_ID",
  "upstream_delta_bp": 60,
  "downstream_delta_bp": -10,
  "reason": "tfbs_context",
  "explanation": "Retain the additional model-site context inspected in the track.",
  "evidence_ids": ["REPLACE_WITH_EVIDENCE_ID"]
}
```

Positive deltas extend outward; negative deltas shorten. Directions are in
**transcript orientation**, including a reverse-oriented loaded sequence.
Reasons are `human_review`, `tfbs_context`, `cutrun_context`, or
`restriction_site`. Each adjustment retains the original candidate, has its
own identity and links to its parent. Removing required context produces a
visible blocked alternative, not an implicit waiver.

Both the per-end shift limit and extension beyond the initial search envelope
default to 200 bp and are independently configurable. For an existing human
selection, place its complete, digest-valid portable ROI in the anchor's
`seed_region` instead of `seed_evidence_ids`. The old ROI remains immutable;
required context is checked again. Adjustments that cannot resolve to an
anchor/seed, or exceed a declared shift budget, fail explicitly.

## Promega MCS And Cloning Practicality

For pGL4.10[luc2], use the existing exact helper-catalog identity, not an
arbitrary sequence labelled "luciferase":

```json
{
  "seq_id": "YOUR_LOADED_VECTOR_SEQUENCE_ID",
  "catalog_id": "Reporter Promega Luciferase AY738222 (online)",
  "suggest_restriction_adjustments": true
}
```

Put this object in `vector`. The catalog describes Promega E6651,
AY738222.1, 4,242 bp, circular, with explicit MCS/luc2 feature expectations.
Its commercial sequence is not bundled. Retrieve and validate it separately;
the candidate operation never fetches it. A custom `helper_catalog_path` is
supported, but the loaded vector must pass that catalog's validation and have
an identifiable MCS. The proposal records vector and catalog hashes, validation,
MCS coordinates, source recognition footprints and cleavage positions.

For each insert, GENtle checks MCS enzyme-pair alternatives and internal insert
sites using its existing restriction engine. Automatic suggestions remove
internal recognition footprints in dispensable boundary context only, bounded
by the shift limit and preserving every required span. Full site counts and
blocked pairs remain inspectable. No compatible pair produces an explicit
alternative-assembly suggestion, not a claim of successful Gibson assembly.

**PCR-added restriction sites and native genomic sites are different.** The
MCS check assumes sites added in primer tails; the genomic insert need not end
at a native site. Native sites are included as context for manual extensions
or shortening. Automatic native-fragment excision, enzyme-specific end-padding
requirements, methylation sensitivity, primer-tail design and reaction
simulation are not performed. Extending a fragment can introduce another
internal site, so each alternative is checked anew. Experimental activity and
cloning feasibility remain separate questions.

## CUT&RUN And Historical Selection

### Explicit Called Peaks

Optionally add `called_peak_sources` to the selection request. Each source is
a **locus-scoped BED** file of existing peak calls, with 0-based half-open
coordinates, a file SHA-256, assembly/chromosome and peak-caller provenance:

```json
{
  "source_id": "reviewed_peak_subset",
  "path": "locus_peaks.bed",
  "sha256": "sha256:REPLACE_WITH_PEAK_FILE_HASH",
  "assembly": "GRCh38",
  "chromosome": "chr7",
  "caller": "YOUR_PEAK_CALLER",
  "caller_version": "EXACT_VERSION",
  "parameters": {"threshold": "EXACT_SETTING", "subset": "EXACT_LOCUS_SELECTION"},
  "cell_line": "EXACT_LOCUS_CELL_LINE",
  "sample": {
    "lane_id": "EXACT_SAMPLE_LANE",
    "source_sha256": "sha256:EXACT_SAMPLE_SOURCE_HASH",
    "replicate_id": "DECLARED_SAMPLE_REPLICATE"
  },
  "control": {"mode": "matched", "control": {
    "lane_id": "EXACT_CONTROL_LANE",
    "source_sha256": "sha256:EXACT_CONTROL_SOURCE_HASH",
    "replicate_id": "DECLARED_CONTROL_REPLICATE"
  }}
}
```

GENtle requires exactly one available source lane per binding and checks its
hash, assembly and cell-line label against the bound locus report. Matched
sample/control sources must differ. A control-free call must explicitly use
`{"mode":"not_used","reason":"YOUR_REASON"}`; it does not acquire an
enrichment claim. Caller version, nonempty parameter map and replicate labels
are mandatory declarations, not independently verified experimental QA. No
peak caller is run here. Use the hash of the lane's actual retained source,
not a BAM hash substituted for a lane sourced from BigWig.

Missing/stale files, invalid intervals/strands/non-finite scores, another
chromosome or any interval outside the loaded locus fail before a proposal is
returned. Calls are not silently clipped. BED3 through extended BED files are
accepted; column 5 is retained as a caller-defined score, never interpreted as
calibrated confidence. Limits are 16 sources, 8 MiB per file and the configured
combined evidence budget. Duplicate source IDs or file hashes are rejected.

The report preserves the full request and exposes `called_peak` evidence with
stable `peak:SOURCE_ID:LINE_NUMBER` IDs and `may_seed_boundary=true`.
Exported ROI references retain the peak caller's version, not the locus's
annotation release, and link the complete selection request by its digest.
Raw coverage, a `Strong` label and `overlapping_peak_count > 0` do **not** enter
this tier automatically. Peak retention cannot establish direct binding,
independent replication or reporter activity. Legacy requests with no
`called_peak_sources` retain coverage-only behavior.

### Fixed-Envelope Comparisons

Without `signal_comparisons`, enrichment is `not_evaluated`. An optional
comparison names sample/control lane IDs, cell line, replicate IDs, compatible
units, a `minimum_mean_difference`, and `missing_policy` (`require_complete`
or explicitly `missing_as_zero`). Source hashes must differ, cell-line labels
must match and overlapping signal bins are rejected. Replicate identity and
units are recorded **caller declarations**, not independently verified facts.
The mean difference is descriptive, not statistical significance.

Each comparison is measured once on the fixed 1,001 bp envelope (or configured
replacement), never optimized over trial insert boundaries. Any passing
declared comparison supplies the descriptive-support ranking flag; individual
results remain available. Truncated/missing coverage stays unassessed. A broad
signal difference prioritizes a neighbourhood; it does not locate the bases
responsible for reporter activity.

Separately, normalized TSS exports now carry historical `selection_window`
metadata when supplied completely: upstream/downstream extent, length and
sequence hash. That window may legitimately differ from the detailed display
window. JSON, figure/legend text, annotated GenBank/EMBL comments and PDF
composition bindings preserve the distinction. Legacy missing geometry stays
unavailable, not assumed equal to the display. This metadata never reselects
TSSs or substitutes historical signal for a new insert assessment.
New binaries read legacy reports with this field absent. Older binaries using
strict unknown-field rejection may reject newly enriched reports; this is not
a promise of forward compatibility with those binaries.

## Review And Continue

Inspect `evidence`, `anchors`, `fixed_comparisons`, each candidate's retained,
excluded and bisected IDs, `blockers`, and optional `cloning`. The report binds
the normalized request, locus file and resulting proposal. IDs and hashes
identify content; they are not scientific approval.

`proposed_region_set` contains eligible alternatives only. After choosing
regions, write that complete set to `regions.json`, then run
`regions import '{"path":"regions.json","format":"json"}'` through the
existing shell operation. Inspect them in Saved Genomic Regions / the locus view and use
the exact-ROI contrast planner. Import is a separate explicit mutation. The
selector has no dedicated graphical editor yet; GUI Shell exposes the same
request without creating a second boundary algorithm.

Conservation, motif scores, annotation and occupancy motivate reporter tests;
none proves autonomous promoter activity or predicts a measured decrease in
luciferase activity. Existing fixed-window, manual ROI and exact-panel routes
remain available as alternatives.

## Fixed-Window Similarity Geometry

Fixed similarity-search windows are not adaptive reporter inserts. The shared,
stateless `ComputeTssWindowGeometry { request }` operation owns strand-aware
window bounds, touching/overlapping unions and regulatory-feature intersections.
Its `gentle.tss_window_geometry_request.v1` request declares an assembly label,
upstream/downstream lengths and groups of same-chromosome/same-strand anchors.
Each anchor binds its existing `TssGeometry` source window. Coordinates are
1-based inclusive; -500/+200 is 701 bp. Gaps are never bridged. Each intersection
must fit a single prepared source, whose anchor IDs remain in the result.
The pure operation checks geometry, not FASTA hashes or assembly authenticity;
the importing/preparation adapter must verify those source bindings separately.

For example, save this operation as `geometry_op.json`:

```json
{"ComputeTssWindowGeometry":{"request":{
  "schema":"gentle.tss_window_geometry_request.v1","assembly":"synthetic",
  "upstream_bp":500,"downstream_bp":200,"groups":[{
    "group_id":"toy","chromosome":"1","strand":"-",
    "anchors":[{"anchor_id":"T1","source":{
      "chromosome":"1","strand":"-","tss_1based":1000,
      "start_1based":200,"end_1based":1800,"upstream_bp":800,"downstream_bp":800
    }}],"features":[{"feature_id":"E1","start_1based":990,"end_1based":1010}]
  }]
}}}
```

Run `gentle_cli shell 'op @geometry_op.json'` or `op @geometry_op.json` in GUI
Shell; MCP `op` and operation-based scripting use the identical contract. The
result in `result.tss_window_geometry` is `gentle.tss_window_geometry.v1`, with
normalized request/hash, per-anchor windows, connected stretches and feature
intersections. This example returns window 800..1500, not 500..1200.
Duplicate IDs, unavailable source span, mixed strand/chromosome or excessive
work fail explicitly. No project is required or changed, and no files are written.

`scripts/prepare_tss_regulatory_similarity_candidates.py --gentle PATH ...`
now consumes that operation via a built `gentle_cli`. It retains
`window_geometry.json`, its hash, the request hash and the engine binary hash.
It still verifies/extracts from the bound promoterome FASTA; there is no Python
geometry fallback and no silent switch to evidence-guided reporter selection.
