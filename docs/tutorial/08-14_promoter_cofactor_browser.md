# Explore promoter cofactors: from TA/DN associations to a saved motif region

**Question:** which sequence motifs occur near TP73 anchors, and which specific
sites would be worth inspecting before a reporter experiment?

Start with an empty GENtle project. No genome, network request, SQL knowledge or
agent login is needed to inspect a supplied package. DuckDB is a local table
reader, not a molecular database: GENtle uses it to query already-delivered
Parquet files. This does not run a new motif scan or fit new statistics.

This registered guide has a reproducible synthetic CLI companion. Native GUI
interaction and human scientific approval remain separate acceptance steps;
no live screenshots or biological validation are claimed here.

## 1. Prepare The Small Teaching Package

With a built GENtle CLI and an explicitly installed DuckDB executable, run from
the repository root (use a **new** output directory):

```sh
python3 scripts/promoter_cofactor_tutorial.py \
  --cli target/debug/gentle_cli --duckdb /absolute/path/to/duckdb \
  --output /tmp/gentle-cofactor-tutorial
```

On Windows the binary is `gentle_cli.exe`; choose a new local output directory.
The script creates `package/`, exact request/report JSONs, a saved project,
`regions.json`, a tiny local reference/catalog and `replay.json`. It executes only repository-authored fixture
SQL to prepare data; every biological query and region capture goes through
GENtle's shared CLI with an explicit temporary project, never the checkout's
current state. It refuses to overwrite an existing destination. Parquet
bytes can differ between DuckDB versions; fresh inventories bind the bytes
actually produced, not a claimed cross-version binary hash. The DNA attachment
exercise prepares only a 240-base synthetic reference; it downloads nothing.
Although its teaching assembly label matches the report, its bases are not
human GRCh38 and do not biologically validate the imported motif score.

The [fixture and origin](../../test_files/fixtures/promoter_cofactors/README.md)
are artificial: two anchors, two shared promoters, a positional motif
`MA9000.1` (SYNTHETIC-A), and an overview-only `MA9001.1` (SYNTHETIC-B).
These are fabricated accessions, not new JASPAR matrices or experimental data.
Do not interpret the toy effect sizes as biological findings.

Open **File > Promoter Cofactors...** or the **Promoter Cofactors** command-palette
entry. Choose `/tmp/gentle-cofactor-tutorial/package`, then **Inspect package**.
Inspection verifies hashes, completion and declared assembly without DuckDB or
a loaded sequence. Subsequent queries need DuckDB; use **Query limits > DuckDB
executable** if the desktop app cannot find the terminal's executable.

The coverage panel is your first checkpoint. Check GRCh38, chromosomes 1 and 2,
one positional motif, retention policy and score settings. The missing requested
factor `SYNTHETIC-MISSING` is unavailable, not a negative experimental result.

## 2. Compare Cohort Associations Before Choosing A Site

Choose distance band `gap_6_20`, then **Rank motifs**. Read the columns together:

| SYNTHETIC-A quantity | Toy value | Interpretation |
|---|---:|---|
| TA adjusted odds ratio [95% CI] | 2.5 [1.1, 4.0] | Estimated cohort association, with uncertainty |
| DN adjusted odds ratio [95% CI] | 0.8 [0.2, 1.8] | Its interval crosses 1; do not overstate depletion |
| TA/DN ratio of adjusted odds ratios | 3.125 | `2.5 / 0.8`, not a binding-probability ratio |
| Difference CI / original BH q | [1.4, 5.0] / 0.03 | Original model diagnostics, not recomputed here |
| Frequency at score >= 0 | 50% (1/2) | Positive anchors / eligible cohort anchors in this band |

The denominator is included TP73 anchors overlapping extended regulatory
promoters, not genes, reads or all genomic bases. Two promoters owning the same
anchor do not make two independent observations. Ranking by effect size is not
a significance filter: select **BH q <=** separately when appropriate.

SYNTHETIC-B can have overview statistics without positional evidence. Its
overview-only status does not mean there were zero sites. Likewise, restricting
the anchor search to one gene does **not** refit cohort statistics for that gene.

## 3. Find An Anchor And Interpret The Retained Match

Set chromosome `1`, start `100`, end `116`, then **Find anchors**. These are BED
0-based half-open coordinates: `[100,116)` is 16 bases; its human-readable
1-based inclusive equivalent is `101..116`. Select the returned anchor,
choose `MA9000.1`, and open **Cofactor detail**. Anchor IDs are package-local;
use coordinates when changing packages.

The retained match occupies `[133,148)`, has raw score `4.25` and strand `+`.
Its gap is **`133 - 116 = 17 bp`**, putting it in `gap_6_20`. A BED end is not
an additional occupied base. This example preserves the whole hit rather than
clipping it to the anchor or to a promoter boundary.

| Exclusive band | Interval distance |
|---|---|
| `overlap` | Negative; shared occupied bases |
| `adjacent_0_5` | 0 through 5 bp; zero means abutting, not overlapping |
| `gap_6_20` | 6 through 20 bp |
| `gap_21_50` | 21 through 50 bp |
| `gap_51_100` | 51 through 100 bp |
| `gap_101_150` | 101 through 150 bp |

Distance is `max(starts) - min(ends)`. Genomic left/right does not automatically
mean transcript upstream/downstream. An anchor may have matches in several
bands: **do not add band frequencies** to obtain an overall frequency.

The package retains one strongest physical locus per anchor/motif/band, not
every motif match. Plus and minus scores remain separate; `.` indicates a
strongest-strand tie, not an inferred direction. A missing strand score is
censored below the retained floor, not a numeric zero.

In this fixture the anchor has SAOS2-TA TP73 support with depth `7`. Control
fields are absent, so control support/depth is **unavailable**, not zero. Other
samples are also unavailable. Read sample and negative-control columns
separately; a depth maximum is not a normalized fold enrichment. H3K4me3 effects,
raw coverage and SK-MEL-29_1 observations are unavailable in this contract.

## 4. Change Presence Without Rewriting The Evidence

Set **Query limits > Presence score >=** from `0` to `5`. Before rerunning,
the form-changed warning means the visible results still belong to the earlier
request. **Copy exact request** and region saving remain bound to that displayed
report, not to newly typed controls. The displayed-result line states the
threshold actually used.

Rerun cofactor detail. The retained `4.25` score no longer meets threshold 5,
but neither its span nor its score disappears. In this fixture both original
counts remain 1:

- `n_source_loci` counts loci retained at score >= -1.
- `n_score_zero_loci` counts loci at the original analytical threshold >= 0.

The changed presence flag neither estimates how many sites pass threshold 5
nor recalculates the cohort odds ratios, q-values or frequencies. A missing
positional row can yield zero only inside verified anchor/motif/band coverage.
Outside that scope GENtle reports unavailability rather than inventing absence.

These imported scores use the producer's declared **raw log2-relative-risk**
configuration. They are not GENtle's `-log10(Ptail)` promoter traces, normalized
affinities, or biochemical binding probabilities. See
[PWM/PSSM scores and promoter traces](08-13_motif_logo_to_promoter_trace.md)
before comparing different score families or matrices.

## 5. Preserve The Site And Navigate Deliberately

Return to threshold 0 if desired and rerun. The anchor, retained-hit and promoter
rows offer **Save evidence region** and **Copy region request**. Saving creates
an explicit record in the existing `promoter_cofactors` region set. It does not
download DNA, edit a sequence or create a CUT&RUN track. Existing records are
not silently replaced; a duplicate save is rejected.

The saved hit retains genomic `[133,148)`, strand `+`, exact raw scores, promoter
membership and gene links, source hashes, query parameters and separate support
summaries in its digest-bound `evidence[].source_record`. The evidence hover in
the existing **Regions...** view summarizes the selected raw score and links.
`promoter-a` has both GENE-A and GENE-B links; keep both instead of picking an
invented single owner. Anchors/promoters without a declared strand stay
unstranded. Neither a TP73 anchor nor a motif hit establishes cofactor occupancy.

**Open hit region...** (or the anchor/promoter **Open region...**) prepares the
existing reference-retrieval dialog. It shows explicit BED-to-1-based conversion;
the hit becomes `134..148`. Choose the matching prepared reference and extract
explicitly to open DNA. This navigation alone does **not** attach the saved
evidence as DNA features.

For the smaller evidence-view handoff, optionally enter an already-loaded
sequence ID before saving. GENtle requires a verified genome anchor, exact
assembly/contig identity and full span containment; it records sequence hash
and a strand-aware local projection. Reverse-oriented sequence views reverse
local strand/coordinates, never the original genomic record. A catalog label
that merely contains `GRCh38` is not an identity proof: unmatched labels are
rejected, without guessing aliases or silently switching assembly. Leave the
field empty to save a portable genomic region without any loaded sequence.

### Explicitly Attach A Stranded Motif To DNA

For a prepared GUI practice project, open the companion's
`demo_plus.preview.state.json` through **File > Open Project...**, then reopen
**Promoter Cofactors...** and click **Reload saved motifs**. It contains the
saved hit and verified synthetic DNA but no attached feature. The companion's
final `tutorial.state.json` already contains the feature and deliberately rejects
duplicate attachment.

Saving remains separate from annotation. After saving a retained hit, open
**Saved motif: DNA annotation**, choose a loaded, genome-anchored sequence and
**Preview DNA annotation**. For the companion's `demo_plus` sequence, genomic
`[133,148)` projects to local `[33,48)` because the sequence starts at genomic
100. Review the reference, interval and strand, tick the approval checkbox,
then **Attach DNA annotation**. Normal project undo removes the new feature.
After reopening the project, or editing a saved region elsewhere, use
**Reload saved motifs** and select the saved hit before making a fresh preview.

The engine checks structured catalog assembly/taxon metadata, the exact contig,
full containment and all anchored DNA bases against the locally prepared
reference. A reference display label is not evidence. It does not download or
substitute another genome. Changing DNA, catalog, saved evidence or existing
annotations invalidates approval and requires a fresh preview. Repeating
attachment of the same saved region is rejected. Keep the original report;
the ordinary feature carries its ROI/source-report JSON, score and gene links.
`gentle_roi_json_base64` encodes that JSON so GenBank line wrapping cannot break
strings inside it; strip wrapping whitespace before base64 decoding.

The offline CLI companion exercises plus-strand attachment end to end; Rust
tests also exercise a verified reverse-oriented `[100,200)` view, where the
same hit is local `[52,67)` on the minus strand. Do not use a newly reverse-
complemented unanchored sequence as proof of genomic orientation. Anchors,
promoters and tied/unstranded hits remain saved regions because the simple
feature editor represents forward/reverse, not unknown orientation. These
annotations indicate sequence association, never protein occupancy or causality.
See [the reference-binding decision](../decisions.md#dec-045-portable-genomic-regions-are-assembly-bound-evidence-ledgers).

## 6. Reproduce Through CLI Or MCP

The preparation script writes these concrete requests and checks their results:

```sh
gentle_cli features promoter-cofactors @/tmp/gentle-cofactor-tutorial/inspect.request.json
gentle_cli features promoter-cofactors @/tmp/gentle-cofactor-tutorial/rankings.request.json
gentle_cli features promoter-cofactors @/tmp/gentle-cofactor-tutorial/detail.request.json
gentle_cli --state /tmp/cofactor-replay.state.json shell 'regions capture @/tmp/gentle-cofactor-tutorial/capture.request.json'
gentle_cli --state /tmp/cofactor-replay.state.json shell 'regions inspect @/tmp/gentle-cofactor-tutorial/inspect_region.request.json'
```

Use a fresh state for this manual replay to avoid intentional duplicate rejection.
`capture.request.json` embeds the complete **displayed** report and target
`{"kind":"hit","anchor_id":1,"motif_id":"MA9000.1","distance_band":"gap_6_20"}`.
Here 1 is the fixture's returned ID; the script obtains it from coordinates.

An MCP client can call the existing `op` tool with
`{"confirm":true,"operation":{"QueryPromoterCofactors":{"request":REQUEST}}}` or
`{"confirm":true,"operation":{"CaptureGenomicRegion":{"request":CAPTURE_REQUEST}}}`,
substituting the JSON objects from the files, not their filenames. Querying is
read-only; capture modifies the project. The generic MCP `op` adapter requires
explicit confirmation for both calls, including the read-only query; do not
set `confirm` without that consent. No SQL or nested agent call is needed. Both routes execute the
same engine operations as the GUI.

The script also writes `demo_plus.preview.request.json`, `demo_plus.preview.json`,
`demo_plus.apply.request.json` and `demo_plus.applied.json`. The apply request
copies the preview's approval digest; it is not reusable after annotation changes.
Use `regions preview-feature` and `regions materialize-feature` with these
request shapes in a fresh matching project. MCP clients use
`PreviewGenomicRegionFeature` and `MaterializeGenomicRegionFeature` through the
same confirmed `op` tool. A successful preview alone changes no annotations.

## Optional: Inspect The Local IRF9 Handoff

This is not a committed dataset or a prerequisite for the tutorial. If you have
the supplied `glen_promoter_cofactors_score0_20260909_v1/package/`, inspect its
manifest and hashes first. Use assembly GRCh38 and **Find anchors** at chr1
`[1778943,1778959)`. Select the exact matching coordinates, then request IRF9
`MA0653.1` detail using the returned ID. Do not reuse the toy ID or an ID from
another package.

The supplied walkthrough describes a forward IRF9 locus `[1778976,1778991)`,
gap `1778976 - 1778959 = 17 bp`, raw score approximately `4.30329`, and
SAOS2-TA TP73 maximum depth `9`. Treat these as expectations to verify against
the actual hash-bound package, not values to repair into a report. Record any
difference and keep sample/control fields separate. Private candidate lists
and the production package must remain outside the repository.

## Acceptance And Next Experiment

The automated companion checks inspection, ranking, coordinate-selected detail,
threshold invariants, lossless shared-engine region capture/export and exact
local-reference preview/attachment on synthetic DNA. Run:

```sh
GENTLE_TEST_DUCKDB=/absolute/path/to/duckdb \
GENTLE_TUTORIAL_BIN_DIR=/absolute/path/to/gentle/bin \
  python3 scripts/test_promoter_cofactor_tutorial.py
```

For Glen's live acceptance: repeat steps 1-5 in the native GUI, edit a threshold
without rerunning, verify the old-result warning and copied request, save a hit,
inspect its evidence, preview/confirm a motif annotation, undo it, and test
matching plus/minus anchors and a mismatched assembly. Record exact revision/package hashes. Screenshot capture requires
the repository's explicit consent policy; a CLI replay is not a screenshot or
human scientific sign-off.

Use these observations to nominate a controlled reporter contrast, not to
declare that the cofactor binds or causes regulation. Retained maxima cannot
reconstruct raw signal or arbitrary multi-fragment architecture. Demonstrating
occupancy or causality requires additional evidence and experiment.
