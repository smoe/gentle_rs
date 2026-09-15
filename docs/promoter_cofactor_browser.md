# Explore a promoter-cofactor package

This browser joins existing cohort statistics, strongest motif matches and
promoter ownership. It does not scan DNA, refit a model or download data.
Start with the registered [empty-project synthetic tutorial](tutorial/08-14_promoter_cofactor_browser.md)
for GUI steps, worked numbers, CLI replay and an optional coordinate-selected IRF9 example.
DuckDB and the package are optional. Only explicit queries need DuckDB:
`duckdb` on PATH, `GENTLE_DUCKDB_BIN`, or request `duckdb_executable`.
The browser's **Query limits > DuckDB executable** field supplies an explicit
path when the desktop application's PATH differs from the terminal's.

## Wet-lab walkthrough

1. Open **File > Promoter Cofactors...**, or find **Promoter Cofactors** in the
   command palette. Select the delivered `package/` directory, not the ZIP or
   its parent. **Inspect package** verifies completion, assembly and SHA-256
   file inventory without DuckDB. Requested factors without a matrix remain
   unavailable. Private groups come from your package, not the public catalog.
2. **Rank motifs** compares TA/DN adjusted odds ratios, their confidence
   intervals, original BH q-values and the ratio of odds ratios. Choose TA/DN
   enrichment, depletion or strongest isoform difference. Effect-size ranking
   is not a significance filter; enable **BH q <=** deliberately. Non-estimable
   rows remain visible after estimable rows unless filtered. Top-N output flags
   additional matches. Filtering never triggers a new multiple-testing correction.
3. Frequency is a **percentage with positive/eligible anchor counts** for
   TP73 anchors overlapping Ensembl extended regulatory promoters in the named
   distance band. It is not gene expression, binding probability or genome-wide
   frequency. Source species describes matrix provenance, not binding exclusivity.
4. Select an exact motif accession. **Find anchors** accepts a chromosome and
   BED interval or an exact Ensembl **promoter gene ID**. The motif/TF name
   filter concerns cofactor rankings; promoter gene ownership is a different
   question. Shared promoters never multiply physical anchor observations.
5. **Cofactor detail** shows the strongest physical match per band, fractional
   or negative scores, strand, distance, fixed-cutoff counts and context
   statistics. TP73/control support and depth stay separate for SAOS2 and
   SK-MEL-29_2. Missing fields mean unavailable. SK-MEL-29_1 is excluded;
   H3K4me3 model effects and raw CUT&RUN coverage tracks are not included.
6. **Open region...** prepares the existing retrieval window using the source
   genome ID and explicit BED-to-1-based conversion. Choose the matching
   prepared reference and explicitly extract to reach the DNA viewer. Nothing
   is downloaded automatically; the browser works without that reference.
7. **Copy exact request** copies the request that produced the displayed
   result, even after form edits. **Copy report JSON** retains full diagnostics,
   source hashes and interpretation limits for collaborators or agents.
   Form edits now show a stale-result warning; the displayed-result line retains
   the actual motif, band and presence threshold instead of adopting edited values.
8. **Save evidence region** on an anchor, retained-hit or promoter row uses
   `regions capture` to save an assembly-bound region in set `promoter_cofactors`.
   **Copy region request** copies that same report-bound request. This works
   without a reference; it preserves evidence, not DNA features. Use **Open hit
   region...** for explicit reference retrieval of the retained motif span.

### Portable Evidence Handoff

The existing `CaptureGenomicRegion` request accepts
`source: {"source_kind":"promoter_cofactor","report":REPORT,"target":TARGET}`.
`REPORT` is the entire displayed query report. Targets are
`{"kind":"anchor","anchor_id":ID}`,
`{"kind":"hit","anchor_id":ID,"motif_id":"ACCESSION","distance_band":"BAND"}`,
or `{"kind":"promoter","regulatory_feature_id":"ID"}`. Use the existing
`regions capture @request.json` route or its typed operation; no new query or
SQL implementation is introduced.

Saved evidence retains the original report as `source_record`, bound by SHA-256
of its compact key-sorted JSON representation. The report contains manifest,
completion and Parquet digests, request, raw strand scores, source coordinates,
promoter/gene links and separate support summaries. Capture validates coverage,
selected geometry and score/strand consistency; it does not reopen the package
or independently authenticate producer claims. Source records are capped at
8 MiB per capture. A report checksum is not a new experimental observation.
The generic signal-intensity field remains empty: raw motif scores are not depth.

Optional source `seq_id` requests existing checked local projection. It requires
an exact verified anchor, matching assembly/contig and full containment; the
sequence digest and strand transform are retained, and genomic coordinates are
unchanged. Catalog labels merely containing an assembly token are not silently
treated as aliases. Leave the GUI sequence-ID field empty to save portable
regions when that exact reference binding is unavailable. Native feature
attachment needs a further checked reference identity/materialization handoff;
see [DEC-045](decisions.md#dec-045-portable-genomic-regions-are-assembly-bound-evidence-ledgers).

This can nominate sites for testing whether a defined sequence change reduces
luciferase activity in an otherwise identical reporter. A cohort association
or motif score does not predict that measured response, prove occupancy or
identify the protein responsible for an observed effect.

## Shared queries

Example `cofactors.request.json`:

```json
{
  "package_path": "/absolute/path/to/package",
  "assembly": "GRCh38",
  "query": "rankings",
  "ranking": "ta_enriched",
  "distance_band": "gap_6_20",
  "max_q_value": 0.05,
  "max_rows": 20,
  "timeout_seconds": 30
}
```

```sh
gentle_cli features promoter-cofactors @cofactors.request.json
gentle_cli shell 'features promoter-cofactors @cofactors.request.json'
```

Workflows and MCP `op` accept
`{"QueryPromoterCofactors":{"request":{...}}}`. GUI Shell, JS/Lua/Python
shared-shell wrappers and inner-agent introspection use the same contract.

| query | Additional input | Result |
|---|---|---|
| `inspect` / `candidates` | none | Integrity, coverage, requested-factor availability |
| `rankings` | optional `motif`, `source_species`, `distance_band`, `ranking`, `max_q_value` | Original cohort statistics, never local-gene refits |
| `anchors` / `promoters` | `region`, `gene_id` or `anchor_id` | Anchors, membership edges and separately owned promoters |
| `anchor_detail` | exact `motif`, package-local `anchor_id` | Zero-complete bands, context statistics and promoter/gene links |

`region` is `{"chromosome":"1","start_0based":1000000,
"end_0based_exclusive":1100000}`. This example does not assert an anchor exists
there. Obtain IDs using `anchors`; bind them to the package hashes.
Cross-package identity uses assembly/chromosome/start/end, not anchor ID alone.

`presence_threshold` defaults to 0 and accepts finite values >= -1. It changes
only the presence flag, not published statistics or the fixed counts
`n_source_loci` (>= -1) and `n_score_zero_loci` (>= 0). Both retained strand
scores survive; a strongest-strand tie is `.`. A missing orientation is censored
below the source floor, not zero. Maxima/counts cannot reconstruct all sites,
arbitrary-threshold counts or pair architecture. Frequencies across exclusive
bands cannot be summed; any future pooled frequency must use anchor-level presence.

Bands: overlap, `adjacent_0_5`, `gap_6_20`, `gap_21_50`, `gap_51_100`,
`gap_101_150`. Distance is `max(starts)-min(ends)`: overlaps negative,
abutment zero. Genomic left/right is not transcript upstream/downstream.
Cofactor spans are not clipped to a query or promoter boundary.

## Safety and limits

The report is `gentle.promoter_cofactor_query.v1`. Inspect `availability`:
missing package/runtime, invalid integrity, assembly mismatch, unsupported
coverage, query failure and excessive detail output are explicit states. No
failure publishes partial biological rows. Sparse zero counts are justified
only inside verified anchor/motif/band scope; overview-only motifs and uncovered
regions are unavailable, not biologically absent.

Requests allow 1-2,000 combined rows and 1-120 seconds including integrity work.
DuckDB uses two threads, 512 MiB memory, no disk spill, no initialization scripts
or automatic extensions, and an 8 MiB output cap. GENtle constructs SQL over
exact inventoried Parquet paths. Package SQL, helpers and DuckDB views are never
executed. No working-directory changes or source-package writes occur. Files
are streamed through SHA-256 checks, not loaded as whole tables. Rechecks reject
replacements during queries. Checksums establish consistency, not independent
authentication of the producer's biological conclusions.

This slice does not materialize annotations, propose sequence edits, export
GenBank/EMBL or reconstruct >=2-fragment blocks from summarized CUT&RUN depth.
Those remain separate explicit operations requiring appropriate source evidence.

## Verification

The [synthetic fixture](../test_files/fixtures/promoter_cofactors/README.md)
contains no private genes/groups. Require real Parquet tests with:

```sh
GENTLE_TEST_DUCKDB=/absolute/path/to/duckdb \
  cargo test --locked --lib promoter_cofactors -- --test-threads=1 --nocapture
```

Without the variable, the Parquet test skips explicitly; pure contract tests
still run. A configured but broken runtime fails. Desktop interaction and real
producer-data acceptance are separate from these synthetic tests.
