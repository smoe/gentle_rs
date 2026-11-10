# Transcript-oriented RNA read-end metrics

GENtle's read-end utilities separate three questions that must not be collapsed:

1. **whole-read length** — a library and sequencing control;
2. **5′-UTR reach** — a useful truncation/control comparison; and
3. **3′-UTR reach** — the annotation-relative measurement of interest.

The utilities are deliberately downstream of the existing observed-SAM engine.
They do not align reads, choose a transcript silently, or infer a poly(A) site.

## Shared Rust Core: Step 1

### Integrated Saved-Service Contract

`rna-reads read-ends` / `RnaReadEnds` now share the typed commands in
`src/rna_read_ends.rs`. They are experimental; scoped execution returns and
the remaining native/platform boundaries are recorded in the
[current return and remaining checks](glen_next_verification_handoff.md#glen-oriented-pileup-return-2026-11-09).
The build schema is `gentle.rna_read_end_build.v1`; explicit source/gene/
transcript and ordered sample bindings are mandatory. One analysis per sample/
gene is admitted. Saved report status, result binding and artifact hashes are
verified before the pure core is called; the supplied normalized annotation is
authenticated against the original S2.3 annotation binding, then checked against
V2 model geometry. Do not substitute the GTF hash for the normalized file hash.

Output completion `gentle.rna_read_end_completion.v1` binds `pileup.json` and
`coverage.json`. Inspection takes the completion's **raw file hash and size**
as `report`, plus a `gentle.rna_read_end_pileup_view.v1` spec. Export uses the
same page and accepts explicit `json`, `tsv`, `svg`, `png`, `pdf` formats and
a fresh output directory. It writes artifacts with closed/synced handles and
publishes completion last through the existing owned S2.3 publication helper.
JSON/TSV/SVG are deterministic for identical saved input and view; PNG/PDF use
the same SVG with the established converters, not a cross-platform pixel promise.
Raster/vector conversion has an explicit 32-million-pixel page ceiling.
PNG/PDF requests check actual dimensions from the shared SVG layout, including
wrapped headers and lanes, with overflow-safe multiplication before conversion.
An oversized mixed-format request refuses before publishing any artifacts;
JSON/TSV/SVG-only export remains available under its normal output-byte limits.
Long labels/statuses wrap without dropping text; SVG and virtualized native
lanes use the renderer's deterministic line heights. Smaller explicit pages may
be needed when wrapping makes an export exceed the unchanged pixel/output cap.
JSON/TSV remain complete and their scientific values are unaffected by layout.

Coverage reports declared/assessed, gene-assigned and selected-transcript-compatible
counts, exclusions, exact-only quantiles and **both** censoring fractions.
Filters/page ranges do not recompute denominators. The genomic frame stays
bound to the selected model plus all gene-assigned
saved blocks, independent of the visible page/filter. Without the optional
original-base import below, A/T remains unavailable (`no_bound_basecalls`),
not an observed zero; legacy composition parameters/provenance remain explicit.
Separate tail input binds an imported file and raw-signal manifest plus supported
caller/version/model/kit/basecaller/configuration. Unsupported/QC-failing rows
never acquire accepted values. These channels never change UTR censoring.
Status/value sorts and filters use conservative unit-level collapsed metrics.
Selected-transcript metrics under each hypothesis remain visible separately;
the viewer does not pick one hypothesis to manufacture a sortable value.

Indexed ingress binds configured executable/version, local alignment/index,
dictionary, interval and declared filters. Full alignment SHA verification is
optional and explicitly labelled; index pairing is a producer declaration plus
hash/size/mtime and tool-warning checks, not an exhaustive index proof. CRAM
requires a bound uncompressed FASTA and adjacent FAI, matching dictionary and
private REF_PATH/REF_CACHE. Alternate header UR locations are refused before
decode, conservatively including otherwise valid external-reference CRAMs.
Non-UTF-8 `samtools --version` diagnostic bytes are explicitly ASCII-escaped
without loss; valid UTF-8 is unchanged. This does not relax strict SAM/header
validation or interpret binary alignments as text inside GENtle.
See [samtools view](https://www.htslib.org/doc/samtools-view.html) and the
[HTSlib reference search order](https://www.htslib.org/doc/reference_seqs.html).
This is a trusted direct-child supervisor, not a process-tree security sandbox
or an operating-system guarantee against stalled filesystem calls.

### Original Bases And Oriented Single-Read Lanes

New builds write `gentle.rna_read_end_pileup.v2`. Existing v1 reports remain
readable without manufacturing basecalls, new summary counters or changing their
content hashes. No new operation, shell grammar or registry path is introduced.

Optionally supply each sample's `original_bases` as a raw path/size/SHA-256
`PrimerPayloadFile`. Its strict JSON document has schema
`gentle.rna_read_end_original_bases.v1`, `analysis_context_sha256` equal to the
saved analysis, `orientation: "original_sequencing_5_to_3"`, a nonempty
`provenance` description, and `reads`. Every read explicitly supplies `unit`
(`source_sha256`, `unit_id`), `read_id` matching its saved QNAME witnesses, and
the full original `sequence`, including soft/hard-clipped or unaligned bases.
Only ASCII DNA/IUPAC letters are admitted; length must equal the saved unclipped
query axis. Missing originals stay unavailable. Do not fill missing hard clips
with reference bases or pass reverse-complemented/trimmed SAM SEQ as an original.
This file is an explicit input/declaration, not an implicit FASTQ join or proof
of the producer's physical orientation/trimming claim. Authentication checks
the supplied bytes and saved association, not biological truth.

`terminal_runs`, optionally on the build request, defaults to
`{"window_bp":60,"minimum_run_bp":10}`. The window must be 1..1,000,000 and
the minimum 1..window. Preserve each maximal exact A/T run at least that long,
with its literal base, full length and zero-based half-open **original-query**
coordinates. A run intersecting the first/last window has `terminal=true`;
internal runs remain visible with `terminal=false` and do not enter terminal
counters. Overlapping windows name both read ends on one run, not two copies.
No mismatches, clipping adjustment or reverse-complement base renaming occurs.
Legacy `composition` parameters remain retained provenance; their
`maximum_mismatches` does not alter this separate exact-run scanner.

Every hypothesis retains an `oriented` record: stable `read_key` over sample
and declared unit, same/opposite/alternatives direction relative to the selected
mRNA, basecall availability, original-file hash/provenance, policy, and all
sequence-hash-qualified variants/runs. Distinct original sequences for one unit
stay `alternatives`; identical duplicates refuse. No more than 64 variants/unit
or 100,000 run records (also bounded by the request observation budget) are
admitted. Resource failures refuse, never truncate.

`read_ends` refer to original-query sides. `transcript_ends` project those sides
only with a known direction; they do **not** establish reaching an annotated
RNA end. Each run keeps its aligned query/reference intervals and
`aligned_to_reference`, true only when saved blocks cover the entire run.
Partial alignments stay false with their covered intervals retained. Off-axis
glyphs attach to the nearest saved block; they do not extrapolate an alignment
through clips, insertions or unavailable sequence.

The x-axis is left-to-right mRNA 5'-to-3' on both genomic strands, with genomic
scale/introns retained. Same-direction lanes are above the divider, opposite
below. Conflicting hypotheses/sequence variants are shown in **Alternatives**,
not resolved or chosen. This is not pore travel direction. Literal A is magenta,
T turquoise, with labels such as `12A` / `18T`; hollow means not fully aligned.
The native virtualized lane and SVG/PNG/PDF use the same layout, including
wrapped labels, group dividers and glyph heights in the actual pixel guard.

JSON and TSV keep the same saved evidence, both UTR statuses, model/annotation
bindings and full unfiltered `orientation_summary`. Counters are per declared
unit, not per hypothesis/variant. Same/opposite/alternative direction counts
partition units; original-base alternatives are a separate counter. A/T counts
can overlap and do not make exclusive read assignments. Unknown bases are not
zero. Filtering/paging never changes these counters or either original coverage
denominator. There is no physical poly(A)-tail length, polyadenylation-site,
complete-transcript or causality verdict.

For executable synthetic input, use the existing preparer's new
`prepare-oriented-plus` or `prepare-oriented-minus` stage instead of `prepare`.
Run the saved-SAM command and `bind-analysis` as in tutorial 07.01, then run
`bind-originals`. This creates `original-bases.json` from independently generated
synthetic originals and a fresh `oriented-build-command.json`; execute that
instead of the legacy build command. Continue with `bind-pileup`, inspect and
export. All stages require fresh paths. The source-matched process companion is
`tests/test_read_end_rehearsal.py`; skipped process cases are not acceptance.

Glen reports the focused `7fa6b15a` return green: renderer 6/6, root library/
actual JS/Lua/MCP 47/47 and rehearsal 8/8, including both real CLI cases without
skips, plus both locked checks. Native GUI is not_run; final platform/release
acceptance is separate. This is the **transcript-bound** contract. A later
exon-locus/common-genome-track/scored-Assignments extension below is authored
separately and is not certified by those results.

### Exon-Locus And Coupled Transcript Projection

The existing operation accepts a new nested command:
`{"action":"locus","command":{"action":"build","request":{...}}}`.
Build schema `gentle.rna_read_end_locus_build.v1` supplies `source_id`,
`region` (contig, contig_length, zero-based half-open exon, flank_left,
flank_right, optional anchor_transcript_id), ordered `samples`, `terminal_runs`,
core `limits`, fresh `output_dir`, `max_output_bytes` and optional `cancel_path`.
The report/page records effective limits, capped at 100,000 models/units,
1,000,000 intervals/observations and 100,000,000 work units; candidate paths
and copied run records have additional 100,000 ceilings. Bounds refuse instead
of truncating scientific output; they are not a claim of an allocator RSS cap.
Each sample has sample_id, the completed observed-SAM `analysis` binding, and
explicit optional raw file bindings `annotation`, `model` and `original_bases`.
No annotation is chosen implicitly. Missing annotation leaves the fixed genome
view usable; missing CDS gives missing UTR facts, not zero reach. A present file
with wrong hash/size still refuses. A new annotation binding recomputes scores;
historical comparisons are neither upgraded nor reused as assignment authority.
Omitted files retain the existing unavailable diagnostics. An explicitly declared
but absent file is tolerated separately as
`annotation_declared_not_found:<sha256>` or
`cds_model_declared_not_found:<sha256>`. Therefore identical request hashes can
produce different content when these optional files appear/disappear; retain the
completion/content hash and actual diagnostics, not just the request hash.
An explicitly bound `original_bases` file still refuses when absent or changed.

`inspect` and `export` under the same nested command accept the completion
`report` binding and `gentle.rna_read_end_locus_view.v1` view. Default projection
is `{"kind":"genome"}`; `{"kind":"transcript","model_key":"EXPLICIT_KEY"}`
uses one key from saved `models`. Other fields are sample_ids (empty means all),
optional class, offset, count (1..2000), max_sublanes (1..20000), width (900..8000).
Genome coordinates stay fixed low-to-high. Transcript coordinates are mature
mRNA 5'-to-3', with introns omitted; full genomic blocks/junctions remain in the
same read record. Both-strand model tracks name genes, models, annotation release
and hashes. Hypothesis sublanes, clipping, outside-window continuation and
supplementary/chimeric hints stay explicit; each read is counted/displayed once.

The locus report `gentle.rna_read_end_locus.v1` ranks all overlapping models
using exported `gentle.locus_assignment.lexicographic.v1` rules. Fewer structural
conflicts precede matched junctions/full chain and discriminating exon bases;
endpoints/proximity only order weak evidence and cannot break structural ties.
An overhang remains a structural conflict, with no new tolerance. A unique
containment-only preferred model is not resolved RNA origin: `model_supported`
requires at least one matched junction on **every** preferred assignment. Pure
exon containment stays `unresolved` with
`exon_containment_only_no_junction_support`; the per-model score, preferred flag
and `antisense_compatible` class are descriptive and otherwise unchanged.
Uncalibrated literal A/T has zero strand weight. Shared sense/antisense sequence
without a discriminator remains `uninformative_overlap`; tied assignments stay
alternatives. Exact/compatible, novel-compatible, antisense-compatible, chimeric
and unresolved states retain reasons/conflicts and per-hypothesis scores.
As in the existing `ExactFullChain` comparison relation, `score.full_chain`
means the full annotated junction chain with at least one intron and no
structural conflict. Internal deletion gaps and short endpoints do not break
that chain; endpoint matches/distances remain separate. `exact_transcript_match`
uses this full-chain condition (subject to the existing higher-priority
antisense/ambiguity classes), not base-exact exon coverage or complete RNA ends.
Unknown splice motifs are not called canonical. Alignment direction, model
strand and conditional `model_supported` molecule-strand evidence are distinct:
double-stranded cDNA BAM direction alone does not establish RNA origin.
Saved search/view completeness declarations remain per read; conditional model
support does not exclude other loci when the original mapping search was limited.

Input/locus/visible denominators, conditional per-model UTR statuses and
nonadditive A/T/alignment summaries are separate. Model selection never changes
genomic evidence, stable read IDs or input counts. Duplicate source units across
samples refuse. Original-query A/T letters, full/partial/no aligned intervals
and missing original-base diagnostics keep the existing semantics. No physical
tail, RNA completeness, APA or causal claim follows.
`Read.membership` distinguishes `admitted-hypothesis` from
`non-admitted-observation-only`; the latter remains in locus counts but has a
separate per-sample `non_admitted_observation_only_units` count. Units without
overlap retain `no_saved_block_in_locus`, except upstream-unavailable units,
which use `upstream_unit_unavailable:` plus the upstream serde snake_case reasons.
Read-level unavailable names use that same vocabulary, not Rust Debug names.
Runs are counted once per unit/query interval/base, not once per original-base
variant; differing alignment states remain `alternatives`. Disagreeing
per-hypothesis UTR metrics also count as `alternatives`, not `missing`.
When CDS exists but projection conflicts, UTR reason is
`model_projection_incompatible`, not `candidate_model_unavailable`.

Genome tracks are unique by annotation hash/transcript ID with all sample model
keys retained in details. A selected transcript track is always above the reads,
including for a minus-strand model. An original run attaching outside the fixed
genome window remains unchanged in JSON/TSV and explicit off-axis text; it is
neither an invalid glyph nor a glyph clamped to an invented edge position.
Masked query bases are not removed from locus base counts. A supplementary
alignment overlapping the primary on the query never joins a hypothesis;
classification uses admitted paths, with non-admitted observations retained.

`Report.work` and `Page.work` record consumed deterministic charge units (as
`Projection.work`), not wall time, CPU time, reads/second or allocator memory.
The union of competitor exons is cached once per sample/model, respecting both
chain and strand; all remaining charge sites and limits are unchanged. The meter
is shared across samples and original-base scanning still charges per base.
The synthetic 30-model/12-exon/20-read deletion-split test is authored, not run
locally; real capacity remains unmeasured. S2 intake still caps one analysis at
4096 units. No throughput promise or cap increase follows.

JSON/TSV retain all candidates, rule/version bindings and unavailable facts;
SVG/PNG/PDF share the same checked scene, with annotation tracks separate from
read rows. Existing 32-Mpixel, fresh-parent/no-clobber and request-bound cancellation
guards still apply. Requesting oversized raster formats refuses publication;
request SVG/JSON/TSV separately to retain those non-raster artifacts.
Use `tests/test_rna_read_end_locus.py` for the general synthetic CLI companion.
Source authoring/Python checks are not Rust, VM, Windows or native GUI acceptance.

`observed_transcript_structures::read_ends::project_read_ends` now provides the
same V2 UTR semantics as a pure Rust-level building block. It consumes a strict
V2 `ModelDocument`, existing typed `ObservedSamAnalysis` and the parsed normalized
annotation authenticated by the caller. It verifies reference/binding equality
and each used model's exact contig/strand/exons/gene agreement, then returns length and 5-prime/
3-prime reach. It reuses the shared genomic-to-query block transform for the
exon-to-transcript coordinates. It never reruns alignment or chooses a model.

`Metric::Observed` carries a value, `RightCensored` a lower bound and `Missing`
only a reason. Derived boundaries expose complete/included/separate/overlapping
or partial/unknown stop evidence without editing original CDS. Missing candidates
remain visible and withhold a complete gene assignment. Results do not promote
SAM query-unit identity into verified original reads or independent molecules.

Each optimal hypothesis retains its observation IDs and every compatible
candidate's saved relations, transcript intervals, metrics and additive clips.
The collapsed unit metrics derive from these records, not another calculation.
`annotation_gene_ids` retain saved match identities even without a CDS model;
`model_gene_ids` match Python's conservative complete-model assignment. Clip
lengths are unaligned query extents, with soft/hard/mixed/none/unknown kinds
from the terminal saved observation where determinable. They never change a
censored value into an observed end. A model/annotation mismatch refuses with
`rna_read_end_model_annotation_geometry`.

This is **not a new CLI/MCP command, file importer or frozen report schema**.
Callers must first authenticate saved artifacts; in-memory hash equality does
not verify actual files, reference bases or exhaustive search. Explicit model,
interval, observation, unit and work budgets fail without partial output.
Malformed identity/path/axis links refuse; incompatible geometry remains missing.
Existing standalone scripts/outputs are unchanged, including their conservative
legacy-v1 reader; this first Rust boundary accepts V2 only.

The [shared synthetic corpus](../test_files/rna_read_end_core/README.md) checks
31 literal cases in both implementations. Python's LF/CRLF corpus checks are
not Rust execution; see the [Glen checklist](glen_next_verification_handoff.md#shared-utr-core-step-1-2026-11-08)
for its historical compilation/tests. Later service/coverage/pileup and explicit
basecall-import boundaries are documented above, not accepted by that checkpoint. No real-data
analysis, hard-wired gene/barcode selection or biological acceptance is added.

## 1. Freeze the transcript model

The integration's step 2 also defines pure coverage/evidence building blocks;
see [the plan and authorized defaults](rna_read_end_integration_plan.md).
Gene-assigned and transcript-compatible denominators are distinct. Coverage
is conditional on the selected model, with min/max bounds over hypotheses,
not an isoform allocation. Exact-only quantiles retain both censoring fractions
alongside missing reasons. The shared service is described above.
Terminal composition stays `no_bound_basecalls` without the explicit original import;
configured tail-import rows preserve their own QC/provenance and cannot change
UTR metrics. No physical tail estimate is inferred from A/T or clips.

`scripts/project_rna_read_end_metrics.py` requires a CDS-aware model document
with schema `gentle.rna_read_end_transcript_models.v2`. Coordinates are genomic,
zero-based and half-open. Each transcript has stable transcript and gene IDs,
contig, strand, exons, original CDS and explicit stop-codon intervals. The document also binds the assembly,
reference FASTA SHA-256 and annotation SHA-256 already present in the completed
`gentle.observed_sam_structures_data.v1` artifact.

This extra binding is necessary because the observed-SAM report retains exon
comparisons but does not itself define CDS boundaries. A mismatched model is
refused.

The 3′-UTR begins **after the annotated stop codon**, not automatically after
the last raw GTF `CDS` interval. The preparer preserves `stop_codon` features,
including codons split across exons. Complete three-base stops that are separate
from, included in or partly overlap the CDS give the same derived boundary.
The raw CDS is never extended or rewritten. Missing/partial stop evidence leaves
3′ reach `missing`; there is no automatic `+3`. Legacy v1 models remain readable
but have unknown stop evidence and cannot supply a known 3′ boundary. The
[projection decision](decisions.md#standalone-transcript-read-end-projection)
records the source-format and compatibility limits. Reprepare a new model from
explicit evidence rather than relabelling an old v1 document as v2.

Prepare the model from an explicit GTF and the exact reference FASTA. Supply
the annotation digest recorded in the completed observed-SAM analysis; it is
kept separate from the raw GTF digest recorded as model provenance.

```sh
python3 scripts/prepare_rna_read_end_transcript_models.py \
  --gtf annotation.gtf.gz \
  --reference-fasta reference.fa \
  --assembly GRCh38 \
  --source-id ensembl-116 \
  --annotation-sha256 sha256:... \
  --output frozen-transcript-models.json
```

## 2. Project one sample

![Transcript-oriented UTR reach and right-censoring](figures/rna_read_end_projection_geometry.svg)

The solid read segment is projected through annotated exons. The 5′ and 3′
measurements begin at the corresponding CDS boundary. A clipped end is a lower
bound; it is not an observed short UTR.

```sh
python3 scripts/project_rna_read_end_metrics.py \
  --observed-sam-data barcode16.observed-sam.data.json \
  --transcript-models frozen-transcript-models.json \
  --sample-id barcode16 \
  --output-dir barcode16-read-end-metrics
```

The fresh output directory contains:

- `metrics.tsv`: the four-column input expected by
  `compare_rna_read_end_metrics.py`, with blank values for censored/missing rows;
- `projections.tsv`: per-read metrics plus transcript and gene identities,
  retaining numeric censored lower bounds separately from exact observations; and
- `report.json`: v2 projection policy, exact input hashes, model schema,
  stop-evidence counts and nonclaims.

Both gene-ranking utilities validate every detailed `projections.tsv` row
before counting it or excluding an unavailable/ambiguous gene assignment.
`observed` and `right_censored` require a nonnegative integer `value_bp`,
including valid zero; `missing` requires an exactly empty value. A missing
censored bound is invalid evidence, not an extra censored count. Refusals
occur before creating ranking outputs and do not rewrite inputs. This is the
**detailed-table** contract: the four-column comparator `metrics.tsv` still
requires blank censored values, with numeric lower bounds kept only in details.

Only saved optimal hypotheses and resolved `exact_full_chain`,
`exact_partial_chain` or `within_exon` transcript relations are eligible. A UTR
value is retained across transcript ambiguity only when all eligible paths and
models give the same value and status.
All compatible candidates are considered before model lookup. A remaining
candidate without a CDS model (for example, a noncoding alternative) keeps
both UTR metrics missing with `candidate_model_unavailable`; it does not
silently disappear. Candidate IDs remain visible; an incomplete model set
does not supply a supposedly complete gene assignment for ranking.

An aligned end that is clipped is **right-censored**: its value is a lower
bound, not an exact length. An alignment that does not cross the relevant CDS
boundary with coding-side support in the same connected projected interval is
**missing**, not zero. Merely touching the boundary from the UTR side is not
an anchor. Conversely, a coding-anchored read ending exactly at the boundary
may legitimately have observed reach zero. Unknown stop evidence affects only
the 3′ metric, not independently available 5′ or whole-read measurements.

## 3. Compare all six samples globally

Concatenate the six `metrics.tsv` files (retaining one header) and run:

```sh
python3 scripts/compare_rna_read_end_metrics.py \
  six-sample-metrics.tsv six-sample-qq
```

This produces the 15 sample-pair comparisons for each metric. Read the three
Q–Q matrices together: a 3′-UTR displacement is more interesting when the
whole-read and 5′-UTR controls do not show the same displacement.
Each sample/metric needs at least two exact observed values. An entirely
missing/censored metric cannot produce a Q–Q distribution; the comparator
refuses rather than inventing zero observations or relaxing its input contract.

## 4. Rank genes with explicit like-library contrasts

The biological samples are SK-MEL-28 (barcodes 13, 15 and 17) and WM793
(barcodes 14, 16 and 18), but the three barcodes per cell line are not
interchangeable replicates. Barcode 16 belongs to the regular-cDNA pair 15/16.
Pairs 13/14 and 17/18 used enrichment and therefore remain separate
directional counterchecks.

Create `contrasts.tsv`:

```text
contrast_id	focus_sample	comparison_sample	library_preparation	evidence_role
regular_cdna	barcode16	barcode15	regular_cdna	primary
enrichment_13_14	barcode14	barcode13	enrichment_13_14	directional_countercheck
enrichment_17_18	barcode18	barcode17	enrichment_17_18	directional_countercheck
```

Then use the detailed projection files:

```sh
python3 scripts/compare_rna_read_end_gene_contrasts.py \
  --projection barcode13-read-end-metrics/projections.tsv \
  --projection barcode14-read-end-metrics/projections.tsv \
  --projection barcode15-read-end-metrics/projections.tsv \
  --projection barcode16-read-end-metrics/projections.tsv \
  --projection barcode17-read-end-metrics/projections.tsv \
  --projection barcode18-read-end-metrics/projections.tsv \
  --contrast-sheet contrasts.tsv \
  --focus-label 'WM793 vs SK-MEL-28' \
  --minimum-observed 10 \
  --maximum-censored-fraction 0.25 \
  --output-dir wm793-gene-ranking
```

The ranking uses the regular-cDNA barcode16-minus-barcode15 median difference.
The enriched pairs contribute separate directional deltas and concordance
counts; they are never pooled into the primary score. Individual reads are not
treated as biological replicates. The report retains censored and missing
counts and excludes a pair when its evidence is too sparse or too censored.

For each gene and sample, the censoring fraction is
`right_censored / (observed + right_censored)`; missing rows are reported
separately and are not in that denominator. The minimum-observed and maximum-
censoring requirements must pass in both samples of a contrast. A percentage
of censored rows among **all** rows is a different summary, not this gate.
Aggregate barcode percentages cannot substitute for per-gene eligibility.
If no gene passes the prespecified thresholds, report that outcome rather
than loosening the gate to force a ranking.

The resulting candidates can motivate follow-up around the TP73–MDM2/MDM4
hypothesis. They do not demonstrate that TP73, MDM2 or MDM4 caused a change.
Library preparation, read truncation, mapping, transcript selection and the
frozen annotation remain alternative explanations.

## Interpretation picture

![Synthetic gene-wise regular-cDNA comparison with enriched counterchecks for WM793 versus SK-MEL-28](figures/rna_read_end_gene_ranking.synthetic.svg)

The committed figure is synthetic. Regenerate its six-sample fixture and SVG:

```sh
python3 scripts/build_rna_read_end_tutorial.py
```

The solid points are the regular-cDNA pair; small ticks are enriched-library
counterchecks. The three panels deliberately use separate scales. Compare direction and
candidate ordering, not the apparent horizontal distance between panels.
