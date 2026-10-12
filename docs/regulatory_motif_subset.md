# Explore Regulatory/TSS Motif Hits Offline

This walkthrough answers: **which retained motif predictions overlap both a
transcript-start window and an annotated regulatory region?** It starts from an
empty project. It neither retrieves DNA nor calls new peaks, changes annotations,
or measures binding. The same shared operation serves CLI, GUI Shell and agents.

This package is different from the [TP73 promoter-cofactor browser](promoter_cofactor_browser.md):
it is not conditioned on a TP73 hit and does not reduce hits to a strongest
representative. It is also not an unfiltered genome scan. Do not substitute one
package's manifest for another's.

## Prepare The Small Example

You need a built `gentle_cli`, Python 3 and a local DuckDB CLI. No genome, JASPAR
download, network connection or private sample is needed. Run from the GENtle
checkout root and choose a fresh output directory:

```sh
python3 test_files/fixtures/regulatory_motif_subset/make_fixture.py /tmp/regulatory-toy --duckdb duckdb
```

On Windows use `python` and a writable absolute directory instead of `/tmp`.
On a Mac where DuckDB is not on PATH, pass `/opt/homebrew/bin/duckdb` and also
use `--duckdb /opt/homebrew/bin/duckdb` in the queries below. GENtle never installs
DuckDB implicitly. The generator refuses an existing destination.

The [fixture README](../test_files/fixtures/regulatory_motif_subset/README.md)
documents the synthetic origin and recreation. `TOY`, its sequences/coordinates,
scores and transcript owners are invented software-test data. The two matrix
accessions are identifier examples, **not scores calculated from those PFMs**.
This is not a model of PATZ1, TP73 or any real gene.

## 1. Discover Matrices, Not Target Genes

```sh
gentle_cli features genomic-motif-evidence --package /tmp/regulatory-toy --inspect --catalog-limit 20 --path /tmp/regulatory-catalog.json
gentle_cli features genomic-motif-evidence --package /tmp/regulatory-toy --inspect --search PATZ --catalog-limit 20
```

The first report contains two matrix records; the second finds `MA1961.2`.
Inspect `regulatory_subset.motif_metadata`, `catalog_total_motifs`,
`catalog_matched_motifs` and `catalog_has_more`. For larger catalogs, repeat with
`--catalog-offset N`; pages contain at most 1,000 records. Search matches
accession or factor name, case-insensitively, not TSS ownership. Original taxa
metadata stays visible; no human-only filter is implicit. A production catalog
of 2,633 matrices is not a list of 2,633 human target genes.

Inspection verifies bounded package metadata but opens no hit Parquets. An
empty search means no matching catalog label, not no genomic binding site.

## 2. Query The Gene's Physical TSS Windows

```sh
gentle_cli features genomic-motif-evidence --package /tmp/regulatory-toy --genome-id fixture-genome --gene TOY --motifs MA0861.2,MA1961.2 --path /tmp/regulatory-all.json
```

Expect **7 orientation-specific hit records**, two physical TSS windows, and
four transcript-owner records. Two transcripts of `TOY` and one of `NEIGHBOR`
share `tss1`; another `TOY` transcript owns the minus-strand `tss2`. A query for
`TOY` retains the other owner of its shared physical start, rather than duplicating
the hits or concealing that owner. `--gene` matches an exact, case-sensitive
annotated gene ID or name. It does not search TF names in the motif catalog.

Useful report fields:

| Field | What to look for |
| --- | --- |
| `hits` | Full BED spans, strand, native score and exact motif version |
| `motif_coverage` | Per-motif retained-score floors, explicit truncation |
| `regulatory_subset.tss_windows` | Physical windows and TSS positions |
| `regulatory_subset.transcript_owners` | Separate TSS-to-transcript/gene links |
| `regulatory_subset.regulatory_gene_links` | Native regulatory annotation links, not nearest-gene guesses |
| `regulatory_subset.hit_annotations` | Tags and all overlapping query/window IDs, without multiplying physical hits |
| `regulatory_subset.verified_files` | Exact metadata and selected-file hashes checked for this query |

For an individual start, use `--tss-id tss1` instead of `--gene TOY` (six hits).
For a genomic interval, no loaded sequence is necessary:

```sh
gentle_cli features genomic-motif-evidence --package /tmp/regulatory-toy --genome-id fixture-genome --region toy=1:300..2101 --motifs MA0861.2,MA1961.2
```

BED starts are 0-based and ends exclusive. Abutting an interval is not overlap;
a hit extending across the boundary retains its full footprint. Both genomic
orientations remain separate. Conversion to a loaded reverse-complemented
sequence happens through the shared checked projection, not by negating scores.

The producer's windows include the TSS base: plus-strand base `t` gives
`[t-700,t+301)`; minus-strand base `t` gives `[t-300,t+701)`, clipped at chromosome
ends. Unclipped windows are **1,001 bp**, not 1,000. GENtle reads their stored
bounds. Selection requires overlap with an actual regulatory/TSS intersection;
a motif bridging two separated annotations is not sufficient.

## 3. Compare A Score Filter Without Losing Negative Values

```sh
gentle_cli features genomic-motif-evidence --package /tmp/regulatory-toy --gene TOY --motifs MA0861.2,MA1961.2 --min-score 0 --path /tmp/regulatory-nonnegative.json
gentle_cli features genomic-motif-evidence --package /tmp/regulatory-toy --gene TOY --motifs MA0861.2,MA1961.2 --min-score -6 --path /tmp/regulatory-below-floor.json
```

The first retains **4** records, including two genuine scores of zero. Without
`--min-score`, the original seven include TP73-identifier scores -4 and -3.5,
and a PATZ1-identifier score -0.5. Floors in this fixture are -5 and -1.
The second request cannot recover scores discarded upstream: it returns retained
records but marks `incomplete_below_storage_floor` and `query_complete=false`.
A source floor is not an empirical operating threshold or a biological cutoff.

Package-native scores are not GENtle-local background-tail significance scores,
occupancy, affinity or a probability of binding. A missing pseudocount is
`regulatory_subset.score_pseudocount=null`, never an assertion that zero was used.

## 4. Distinguish Empty, Missing And Partial

Repeat the gene query with `--gene EMPTY`: both verified chromosome-2 hit
files contain zero rows. Repeat with `--gene MITO`: the MT inventory explicitly
declares `known_empty_intersection` and has no hit file to open. Neither result
means no motif matches exist elsewhere on those chromosomes.

Other outcomes stay separate:

- `--region chr1:300..1301` is unsupported in this fixture, which names its
  chromosome `1`. This reader makes no contig-alias or liftover guesses.
- An accession absent from the catalog is `motif_not_in_package`, not zero hits.
- A missing/corrupt selected file is `invalid_package`, not an empty answer.
- `--max-rows 2` returns an explicitly truncated, incomplete list.
- Missing DuckDB/package, incompatible assembly and execution-budget failures
  are typed states. Inspect `availability` before counting hits.

Completeness always means the declared **regulatory/TSS subset under its source
floors**, never the complete genome or experimental binding. Source assembly
and hashes establish internal consistency, not independent reference validation.

## 5. Save, Inspect And Reuse

`--path` saves a standalone `gentle.genomic_motif_evidence.v1` report. It retains
the original values, scope, annotation releases, ownership, limits and selected
file hashes. It does not change the empty project. To inspect it without DuckDB:

```sh
python3 -m json.tool /tmp/regulatory-all.json
```

For real, **exactly matching** TSS reports, use the existing attachment route:

```sh
gentle_cli features tss-tfbs-profiles-export --report /path/to/original/report.json --genomic-motif-evidence /path/to/query.json --output-dir /fresh/enriched-tss --formats svg
```

Do not attach the synthetic example to a real gene. Reference, chromosome,
matrix width and geometry checks must pass. Imported hits stay in separate
labelled/scaled lanes, never spliced into local score arrays. Open the matching
annotated TSS DNA in GENtle and attach the enriched `report.json` via
[Attach the quantitative report](gui.md#attach-the-quantitative-report).
The regulatory-subset label and coverage survive report replay, even without
the package. Re-export that enriched report without new attachment arguments
into another fresh directory; it needs neither DuckDB nor the original query
file. Historical results are not overwritten.

## GUI And Inner-Agent Parity

Every `features ...` line above works unchanged in the GUI Shell; omit the
`gentle_cli` prefix. An agent can propose exactly those commands or call the
typed `QueryGenomicMotifEvidence` operation. For example:

> Inspect this local regulatory/TSS package. Find the exact PATZ1 matrix
> accession, then query TOY's annotated starts for it and MA0861.2. Keep
> negative retained scores, show ownership and incomplete coverage, and save
> the report. Do not annotate DNA or download anything.

From a genome-anchored DNA viewer, **TFBS scan > Precomputed hits** also uses
this reader when `GENTLE_JASPAR_GENOME_SCAN_PACKAGE` points to this package.
The inspector shows subset scope, physical windows, transcript owners and
native regulatory links. Opening the package never creates sequence features.
Flat-file TSS imports need not have an engine genome anchor: use the exact BED
interval/shared Shell path and report attachment rather than inventing one.
This release does not add direct package-to-DNA annotation approval or a
dedicated package-selection wizard.

## Verification And Limits

The real-DuckDB tests recreate these Parquets offline and exercise query/export
inventory selection, saved JSON, relocation, negative scores, strands, zero/MT
coverage, integrity failures and bounded metadata. Run explicitly:

```sh
cargo test --locked --lib genomic_motif_evidence::regulatory::tests -- --include-ignored --test-threads=1
cargo test --locked --lib regulatory_subset_scope_survives
```

Ordinary test runs skip the four real-DuckDB tests when their external tools
are not requested. They are not mocked integration tests. Metadata validation
also tests the 65,825-entry production inventory shape without creating a large
genome dataset. Limits include 64 motifs, 256 physical windows, 20,000 ownership
or regulatory-context rows, 1 GiB selected hit files, 64 MiB subprocess output
and one overall timeout (default 30 seconds). Catalog inspection is paginated;
queries never glob or hash all genome payloads. Large queries must be split.

The remote production package and live native GUI still require Glen's
exact-revision acceptance. This synthetic example does not certify that export
as complete or any real biological conclusion.
