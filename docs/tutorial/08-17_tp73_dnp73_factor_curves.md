# Compare E2F1, PATZ1 and TP73 motif curves at the human ΔNp73 TSS

This tutorial answers one concrete regulatory question on a real human locus:

> At the internal ΔNp73 transcription start of human *TP73*, where in the
> −500..+200 bp window do the E2F1, PATZ1 and TP73 sequence models score, and
> which windows could not be scored at all?

Everything runs offline from the retained public RefSeq record already in this
repository. No network request, no prepared genome, no private sample and no
motif database query is involved.

Sequence models locate *candidate* sites. Nothing below shows measured
occupancy, and a high model score is not evidence that a factor binds, that it
acts at this start, or that a mutation would change reporter activity.

## Why this promoter and these three factors

*TP73* carries two promoters. The upstream P1 promoter produces full-length
TAp73; an internal promoter in intron 3 produces the N-terminally truncated
ΔNp73, which antagonises TAp73 and p53. The two are usually discussed together
because TAp73 induces ΔNp73, so the internal promoter is a documented point of
negative feedback, and E2F1 is a documented activator of the *TP73* locus.

That makes three motifs worth showing side by side in one window:

| Factor | Exact accession | Matrix length | Why it is in this view |
| --- | --- | --- | --- |
| TP73 | `MA0861.2` | 16 bp | p73 family motif; the autoregulatory candidate |
| E2F1 | `MA0024.3` | 12 bp | documented activator at this locus |
| PATZ1 | `MA1961.2` | 11 bp | additional candidate factor under test |

The accessions are exact and must stay exact: GENtle refuses factor-name
aliases, `ALL` and consensus fallbacks on this route, because a curve is only
interpretable next to the matrix that produced it.

Choosing three motifs is a display decision. It is not a claim that these are
the only, or the strongest, candidates in this window.

## What you will accomplish

In about 15 minutes you will:

1. open the retained public *TP73* locus and let GENtle recover its genome anchor;
2. preview every annotated *TP73* start and see which ones are unusable;
3. approve exactly one start, the ΔNp73 window, as a TSS window;
4. inspect it in the native **TSS / Regulatory** view;
5. compute the three local motif curves and read them without over-claiming;
6. export the same view as a deterministic SVG.

Every step below is available three ways — GUI, Agent Assistant / shared
operation, and a request you can hand to the inner agent — and all three drive
the same engine contract. Steps 5 and 6 additionally have a fully headless
route that needs no GUI host at all.

## Before you start

- Start GENtle from a repository checkout so the retained record is present.
- Keep the repository root as your working directory.
- Before any `ui ...` command, focus the intended DNA window and open Agent
  Assistant from there. GENtle binds that explicit launch subject and never
  falls back to some other project sequence.
- Retained public input: [`test_files/tp73.ncbi.gb`](../../test_files/tp73.ncbi.gb)
  — NCBI RefSeq `NC_000001.11 REGION: 3652516..3736201`, GRCh38.p14, 83,686 bp.
  Its origin and recreation steps are recorded in
  [`test_files/README.md`](../../test_files/README.md).

The CLI examples use `target/debug/gentle_cli`; a release or packaged binary
behaves identically.

## 1. Open the TP73 locus

**GUI**

**File → Open Sequence…**, then `test_files/tp73.ncbi.gb`.

**Agent Assistant / shared operation**

```text
/open file test_files/tp73.ncbi.gb --id tp73_locus
ui open sequence-window tp73_locus
```

**Ask the inner agent**

> Open the retained public TP73 RefSeq record `test_files/tp73.ncbi.gb` as
> `tp73_locus` and open its sequence window. Report the genome anchor GENtle
> detected and any anchor verification warning verbatim.

GENtle reads the anchor out of the record's own `ACCESSION ... REGION:` header:

```text
Detected GenBank genome anchor for 'tp73_locus': 1:3652516-3736201 (GRCh38.p14, strand +, verification n/a)
```

You will also see a warning that the anchor could not be verified against a
local catalog, because GRCh38.p14 is not prepared here. That is expected and
honest: the coordinates come from the file, and nothing has re-derived them from
a reference. The warning is not a defect and does not block this tutorial.

## 2. Preview every annotated start before choosing one

Save this request as `tp73-inventory.json`:

```json
{
  "seq_id": "tp73_locus",
  "gene_query": "TP73",
  "collection_id": "tp73_tss",
  "upstream_bp": 500,
  "downstream_bp": 200
}
```

**GUI**

**TFBS scan → Transcript starts / TSS windows…**, enter the same gene query,
collection ID and flank sizes, then **Refresh inventory**.

**Agent Assistant / shared operation**

```text
promoters tss-inventory @tp73-inventory.json
```

**Ask the inner agent**

> Preview the annotated TP73 transcript starts in `tp73_locus` with 500 bp
> upstream and 200 bp downstream into collection `tp73_tss`. List each start's
> genomic position, transcript count and availability, and tell me which ones I
> cannot use. Do not approve anything yet.

Five annotated starts are reported:

| Genomic start (1-based) | Transcripts | Availability |
| --- | --- | --- |
| 1:3,652,516 | 6 | `missing_flanks` |
| 1:3,669,080 | 1 | `available` |
| 1:3,675,326 | 1 | `available` |
| **1:3,690,672** | **6** | **`available`** |
| 1:3,698,046 | 1 | `available` |

Read the first row carefully. `1:3,652,516` is the P1 / TAp73 start, and it is
the very first base of this record — there is no upstream sequence here, so a
−500 bp window cannot be built. GENtle reports `missing_flanks` instead of
silently returning a shorter or padded window.

That is a presentation limit of *this excerpt*, not a statement about the
promoter. P1 is real and well studied; it is simply not inspectable in this
file. To work on it you would extend the anchor upstream from a prepared genome
first — a separate task, deliberately not hidden inside this tutorial.

`1:3,690,672` is the internal ΔNp73 start, it carries six transcript models,
and it has the full flank. That is the window used below.

The preview also prints an `approval_sha256`. It binds this exact inventory,
and the next step refuses to run without it.

## 3. Approve exactly one start

Take the `approval_sha256` and the `tss_id` of the `1:3,690,672` row from your
own preview output — do not copy an ID from this page, and never invent one.

```json
{
  "inventory": {
    "seq_id": "tp73_locus",
    "gene_query": "TP73",
    "collection_id": "tp73_tss",
    "upstream_bp": 500,
    "downstream_bp": 200
  },
  "expected_approval_sha256": "sha256:PASTE_FROM_YOUR_PREVIEW",
  "selected_tss_ids": ["tss_PASTE_FROM_YOUR_PREVIEW"]
}
```

**GUI**

Select only the `1:3,690,672` row, then **Approve selected starts**.

**Agent Assistant / shared operation**

```text
promoters tss-materialize @tp73-materialize.json
```

**Ask the inner agent**

> Using the preview you just showed me, approve only the 1:3,690,672 start into
> collection `tp73_tss`. Use that preview's own approval hash and TSS ID. Show
> me the request before you run it.

One window sequence is created. Its ID is `tp73_tss_` followed by the TSS ID;
the commands below call it `TP73_WINDOW`. Re-running the same approval reuses
the existing window rather than creating a duplicate.

The window is annotation-derived. Its `Origin=project_annotation_derivation`
comment states that it is not a verified external report bundle.

## 4. Inspect the window natively

**GUI**

Focus the new window, then choose **TSS / Regulatory** beside Standard map.

**Agent Assistant / shared operation**

```text
ui open sequence-window TP73_WINDOW
ui open tss-view --collection tp73_tss
```

**Ask the inner agent**

> Open the approved `tp73_tss` window in the native TSS / Regulatory view and
> describe the lanes that exist before any scoring.

The ruler shows three aligned coordinate systems for the same base: local
1..701, genomic 3,690,172..3,690,872, and signed TSS-relative −500..+200. Local
base 501 is the annotated start itself.

The title reads `TP73 | TSS 1:3690672 (+) | GRCh38.p14`.

At this point the view contains only what the record supplied: the gene lane,
the TSS marker, and one exon/CDS context lane per transcript model that starts
here — `NM_001126240.3`, `NM_001126241.3`, `NM_001126242.3`, `NM_001204189.2`,
`NM_001204190.2` and `NM_001204191.2`. There are no score curves yet, and their
absence is not a statement about the sequence.

## 5. Compute the three factor curves

**GUI**

Open **Local scoring**, enter `MA0861.2, MA0024.3, MA1961.2`, choose score kind
`llr_background_tail_log10`, then **Compute local scores**.

**Agent Assistant / shared operation**

```text
ui open tss-view --local-score MA0861.2,MA0024.3,MA1961.2 --score-kind llr_background_tail_log10
```

**Ask the inner agent**

> In the active TP73 TSS view, compute local curves for exactly MA0861.2,
> MA0024.3 and MA1961.2 using `llr_background_tail_log10`. Then tell me, per
> factor, how many strand-windows were evaluable and where the strongest window
> starts lie relative to the TSS.

**Headless, without a GUI host** — the same scoring, straight to a file:

```text
promoters tss-view-svg TP73_WINDOW tp73-dnp73-factors.svg --motif MA0861.2,MA0024.3,MA1961.2 --score-kind llr_background_tail_log10
```

Three **Locally computed** lanes appear, each with its own scale, its own
printed score kind and its own matrix SHA-256. For this window all three are
fully evaluable:

| Lane | Evaluated strand-windows |
| --- | --- |
| `Locally computed \| TP73 \| MA0861.2` | 1372 / 1372 |
| `Locally computed \| E2F1 \| MA0024.3` | 1380 / 1380 |
| `Locally computed \| PATZ1 \| MA1961.2` | 1382 / 1382 |

The counts differ only because the matrices differ in length: a 16 bp matrix
has fewer complete window starts in 701 bp than an 11 bp matrix, and each start
is scored on both strands. Every lane ends in a grey band marking the trailing
positions where no complete motif window fits.

### Reading these lanes honestly

- Each x position is a motif-window **start** in the displayed orientation, for
  forward (solid) and reverse (dashed) matches. It is not a base-wise
  occupancy profile.
- The three lanes keep **separate scales** and are not cross-calibrated. A
  taller curve in one lane does not mean that factor is more likely to bind.
  Comparing across matrices is exactly what these axes do not support.
- Amber marks windows that could not be scored — in the upper half for local
  `+`, the lower half for local `−`. Amber is not a low score, and the trailing
  grey band is not a score of zero either. This window happens to contain no
  ambiguous bases, so you should see no amber here; on windows with `N` runs you
  will.
- With `llr_background_tail_log10` and clipping on, zero is a real displayed
  value. Do not read a flat stretch at zero as missing data — that is what the
  amber and grey bands are for.
- The curves are computed from this window's bases by the shared engine scorer.
  They are independent of any attached report and of imported DuckDB evidence,
  and they never overwrite either.

Locally computed curves are model output on genomic DNA. They do not show
occupancy, they do not establish that TAp73 autoregulates through a site in this
window, and they do not predict reporter behaviour.

## 6. Export the view

**GUI**

**Export View SVG**. This exports the enabled and filter-matched lanes over the
current horizontal span, including lanes scrolled out of sight — not a
screenshot.

**Agent Assistant / shared operation, and headless**

```text
promoters tss-view-svg TP73_WINDOW tp73-dnp73-factors.svg --motif MA0861.2,MA0024.3,MA1961.2 --score-kind llr_background_tail_log10
```

Narrow the export to the proximal promoter with a local 1-based span:

```text
promoters tss-view-svg TP73_WINDOW tp73-dnp73-proximal.svg --motif MA0861.2,MA0024.3,MA1961.2 --score-kind llr_background_tail_log10 --span 401..701 --width 1600
```

Repeat `--score-kind` on every export. Omitting it does not inherit the score
kind you used in the view — it falls back to the route default `llr_bits`, and
two exports that differ only in that flag are not comparable.

**Ask the inner agent**

> Export the native TP73 TSS view with those three motif curves to
> `tp73-dnp73-factors.svg`, then export local 401..701 separately. Report the
> lane count and each file's SVG SHA-256.

The operation reports the lane count, the exported local span and the SVG's
SHA-256, so the same command can be shown to have produced the same figure. The
export writes atomically and refuses to render non-finite scores.

The headless route accepts `--report REPORT_JSON` as well, so an agent can
attach a validated TSS profile report and locally computed curves in one
export. Report-backed curves, imported DuckDB hits and locally computed curves
then stay in separate lane classes with separate scales.

## 7. State the result conservatively

A defensible summary of this session:

> In the −500..+200 bp window around the annotated internal ΔNp73 start of human
> *TP73* (GRCh38.p14, 1:3,690,672, plus strand), the E2F1 `MA0024.3`, PATZ1
> `MA1961.2` and TP73 `MA0861.2` matrices were scored on both strands with
> `llr_background_tail_log10` over all complete window starts. All windows were
> evaluable. Per-matrix score positions are reported on separate,
> non-cross-calibrated scales. The upstream P1 / TAp73 start was not assessed,
> because the retained record provides no upstream flank for it.

What this does **not** support: that any of the three factors occupies this
promoter; that ΔNp73 is autoregulated through a site shown here; that one factor
outranks another; or that the unassessed P1 promoter lacks such sites.

Occupancy needs CUT&RUN or ChIP at this locus; functional relevance needs
perturbation and a reporter or endogenous readout.

## How this fits the neighbouring tutorials

- [08-15](08-15_tss_collection_gui.md) covers the TSS collection workflow
  itself — inventory, approval, stored collections and their registry.
- [08-16](08-16_tss_regulatory_view_gui.md) teaches the same native viewer on a
  synthetic **minus-strand** fixture and, unlike this page, walks through
  attaching a hash-bound profile report with imported DuckDB evidence. Use it
  for the minus-strand coordinate convention and the report-attachment path.
- [08-06](08-06_promoter_motif_control_comparison.md) takes these same factors
  across several promoters against a matched control set, which is the right
  tool for asking whether a motif is enriched.
- [08-13](08-13_motif_logo_to_promoter_trace.md) connects a single matrix logo
  to its promoter trace.

This page deliberately stays at one window and one locus.

## Provenance of the retained teaching input

`test_files/tp73.ncbi.gb` is a public NCBI RefSeq excerpt
(`NC_000001.11 REGION: 3652516..3736201`, GRCh38.p14, retrieved 26-AUG-2024),
already used by other GENtle tests and runbooks; see
[`test_files/README.md`](../../test_files/README.md). No new fixture is added by
this tutorial, and no private or experimental data is used anywhere in it.

Genome coordinates, transcript models and the two-promoter architecture come
from that record's own annotation. The motif matrices come from GENtle's
bundled JASPAR registry, identified by exact accession and by matrix SHA-256 in
every lane and export.
