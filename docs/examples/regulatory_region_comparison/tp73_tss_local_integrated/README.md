# TSS-local regulatory-feature similarity in the tall locus reports

> **Corrected replacement.** These retained reports were regenerated from exact
> producer revision `ff44ebdec0ff46e3ea485cd3a36ef3d4c403b6f9`, above the
> reviewed TSS-background renderer change `c6e34de3`. They replace, but do not
> rewrite, the historical reports. The withdrawn `c91bf912` PDFs remain
> non-mergeable because their chromosome-7-only adapter emptied the CD44 and
> TGFB1 lanes.

The replacement uses the chromosome-general converter below and a separate
validation receipt. All 12 requested CUT&RUN/H3K4me3 lanes are nonempty for
each of CD44, TGFB1, and SERPINE1; their native BigWig chromosome overlaps,
source hashes, assembly identifier, interval counts, and strand-aware local
coordinates were checked. Missing or incompatible intervals are never treated
as measured zero signal. The assembly check accepts only the independently
declared `GRCh38` identifier, not words or numbers borrowed from the promoterome
catalog label.

Follow-up review found that the original validation could accept a consistently
mislabelled locus against a different promoterome, or fabricated interval values
within an otherwise overlapping BigWig. Preparation now checks the assembly in
the receipt-bound genome catalog entry, and lane validation compares the complete
interval/score multiset after requested filtering and locus clipping. It follows
the importer's six-decimal score storage and the sequence anchor's orientation,
which need not equal the gene's strand. The retained `ff44ebde` validation receipt
predates these stronger checks and must be refreshed against the original
request/report JSON; do not relabel it as a new validation. Independent local
review did match all 36 lanes / 11,797 SVG interval coordinates, source hashes
and plotted heights to the local BigWigs. The figures need not be redrawn unless
the stronger full-input check identifies a mismatch.

This bundle extends the existing tall `CD44`, `TGFB1`, and `SERPINE1`
promoter–reporter architecture reports. Their transcript models, proposed
reporters, Ensembl Regulation annotations, TP73 CUT&RUN/H3K4me3 lanes, and
15-factor JASPAR score tracks remain in place. A new section immediately before
the interpretation/provenance footer compares every Ensembl regulatory feature
intersecting a transcript-oriented **−500/+200 bp** TSS window with the prepared
human promoterome.

The TFBS tracks are deliberately shown once per gene report. The similarity
section reuses the same genomic stretch and adds, per overlapping regulatory
feature:

- a position-frequency strip counting distinct other genes;
- top other-gene promoter windows with their genes and transcript mappings;
- identity intensity and target-promoter block numbers;
- red boundaries only when alignment-relative order or orientation changes;
- an explicit lower-bound label when either the target cap or per-target HSP
  cap was reached.

New tall reports include a labelled TSS-stretch track immediately above the
transcript models, on the same genomic axis. Its exact stretch IDs and colours
reappear in the similarity section, with reciprocal links in SVG viewers that
support them. Each stretch occupies a separate row to avoid label collisions.
Translucent bands in the same colours extend behind the upper genomic lanes,
from the transcript-model heading to just before the similarity section. They
follow the bound genomic axis on either strand, do not intercept mouse input,
and do not extend into the independently scaled similarity plots or footer.
The bands identify comparison windows, not extra occupancy or functional evidence.
The original transcript, reporter and evidence lanes are translated intact,
not recalculated. Similarity-strip query offsets are assembly-forward; they
are mirrored when displayed on a minus-strand gene axis.

The comparison uses the corrected `exclude_overlapping_or_shared_gene_windows.v1`
policy. Gene-name text is not used to exclude targets. A run without observed
cap saturation is not described as complete, and an absent hit is not evidence
of uniqueness.

## Retained outputs

- `*_with_TSS_similarity.{svg,pdf,png}`: the three tall reports in the requested
  presentation grammar;
- matching `*.receipt.json`: input/output hashes plus exact SVG renderer command,
  version, and executable hash;
- `candidate_regions.{json,fa}`: 15 clipped Ensembl-feature intersections;
- `selected_tss_candidate_regions.{json,fa}`: the seven TP73-supported source
  TSS windows;
- `cutrun_lane_validation.json`: exact request/report/source bindings and the
  nonempty per-lane interval counts for all 36 known input lanes;
- `comparison.json`: corrected counting, cap audit, thresholds, and SHA-256
  bindings for the regenerable raw tables;
- `feature_summary.tsv`, `top_matches.tsv`, and `interpretation.md`: conclusions
  derived from validated outputs rather than fixed prose;
- `CD44_TGFB1_SERPINE1_tss_local_integrated.pdf`: compact four-page supplement.

The raw BLAST HSP and match tables are about 109 MiB and are intentionally not
versioned. Their exact hashes and relative paths are retained in
`comparison.json`; the commands below regenerate them.

## Reproduce

For a corrected run, use a new output directory and record the current GENtle
revision. Replaying the historical revision reproduces the historical defects.
The promoterome receipt, selected transcript membership, chromosome, strand,
TSS geometry, and selected sequence digest are validated before extraction.

Re-export each original locus request with the repository's chromosome-general
BigWig converter. It implements the two-path `bigWigToBedGraph` interface and
emits every chromosome in the source file; it must not be replaced by the
historical chromosome-7 analysis helper.

Both the converter and lane validator require `pyBigWig` in the Python
environment used to launch them. The converter's executable wrapper uses the
`python3` selected by the environment.

```bash
export GENTLE_BIGWIG_TO_BEDGRAPH_BIN="$PWD/scripts/bigwig_to_bedgraph.py"
gentle_cli /path/to/original.state.json \
  gene-locus prepare @/path/to/fresh/request.json
```

Before feature preparation, validate the three re-exported reports against
their exact requests and source BigWigs. This gate requires 12 available,
nonempty lanes per gene, verifies native BigWig overlap on the bound chromosome,
checks source hashes and assembly identity, and checks rendered local/genomic
coordinates on both strands. Missing intervals are rejected rather than treated
as measured zero signal.

```bash
python3 scripts/validate_tp73_locus_cutrun_lanes.py \
  --assembly-id GRCh38 \
  --request CD44=/path/to/fresh/CD44.request.json \
  --request TGFB1=/path/to/fresh/TGFB1.request.json \
  --request SERPINE1=/path/to/fresh/SERPINE1.request.json \
  --report CD44=/path/to/fresh/CD44.report.json \
  --report TGFB1=/path/to/fresh/TGFB1.report.json \
  --report SERPINE1=/path/to/fresh/SERPINE1.report.json \
  --source-revision "$(git rev-parse HEAD)" \
  --output /path/to/run/cutrun_lane_validation.json
```

The locus JSON must carry a matching `sequence_binding.genome_anchor`; re-export
legacy reports that lack it. Both renderers require the exact preparation-bound
JSON bytes, including each candidate's source-report digest. For tall reports,
preparation also records the declared original SVG's hash in
`source_bindings.locus_svgs`, keyed by gene, and verifies its schema/panel ID
against that JSON. This binds the supplied pair; it does not prove the biological
truth of a drawing. The appender requires both original files and rejects an
already-extended SVG or existing output paths. Legacy candidate files without
SVG bindings must be prepared and compared again, not edited to bypass checks.

```bash
python3 scripts/prepare_tp73_cutrun_promoter_candidates.py \
  --promoterome /path/to/human-grch38-ensembl116-promoterome \
  --track-root /path/to/cutandrun_20250602_noDuplicates \
  --output /path/to/run/selected-tss \
  --source-revision "$(git rev-parse HEAD)"

python3 scripts/prepare_tss_regulatory_similarity_candidates.py \
  --selected-tss /path/to/run/selected-tss/candidate_regions.json \
  --promoterome /path/to/human-grch38-ensembl116-promoterome \
  --catalog /path/to/exact/genomes.json \
  --assembly-id GRCh38 \
  --locus-report /path/to/CD44_luciferase_planning_EnsemblReg_15TF.report.json \
  --locus-report /path/to/TGFB1_luciferase_planning_EnsemblReg_15TF.report.json \
  --locus-report /path/to/SERPINE1_luciferase_planning_EnsemblReg_15TF.report.json \
  --locus-svg CD44=/path/to/CD44_luciferase_planning_EnsemblReg_15TF.svg \
  --locus-svg TGFB1=/path/to/TGFB1_luciferase_planning_EnsemblReg_15TF.svg \
  --locus-svg SERPINE1=/path/to/SERPINE1_luciferase_planning_EnsemblReg_15TF.svg \
  --output /path/to/run/feature-candidates \
  --source-revision "$(git rev-parse HEAD)" \
  --upstream-bp 500 --downstream-bp 200

python3 scripts/compare_candidates_to_promoterome.py \
  --candidates-json /path/to/run/feature-candidates/candidate_regions.json \
  --query-fasta /path/to/run/feature-candidates/candidate_regions.fa \
  --promoterome /path/to/human-grch38-ensembl116-promoterome \
  --output /path/to/run/feature-candidates/blastn-40bp-80pct \
  --task blastn --min-alignment-bp 40 --min-identity-pct 80 \
  --max-evalue 1e-5 --max-target-seqs 1000000 --max-hsps 10
```

`--catalog` may relocate the catalog recorded in the promoterome receipt; its
SHA-256 must still match. Without the option, preparation uses the recorded
path. Assembly identity comes from the selected entry's Ensembl `file_stem`
or NCBI `ncbi_assembly_name`, never an arbitrary word in the catalog label.
Missing, conflicting or nonmatching assembly metadata is rejected.

Render a tall report and its bound PDF/PNG derivatives (repeat for each gene):

```bash
python3 scripts/append_tss_similarity_to_locus_report.py \
  --base-svg /path/to/SERPINE1_luciferase_planning_EnsemblReg_15TF.svg \
  --locus-report /path/to/SERPINE1_luciferase_planning_EnsemblReg_15TF.report.json \
  --gene SERPINE1 \
  --candidates-json /path/to/run/feature-candidates/candidate_regions.json \
  --comparison /path/to/run/feature-candidates/blastn-40bp-80pct/comparison.json \
  --matches /path/to/run/feature-candidates/blastn-40bp-80pct/matches.blastn.tsv \
  --hits /path/to/run/feature-candidates/blastn-40bp-80pct/hits.blastn.tsv \
  --output-svg /path/to/run/reports/SERPINE1_with_TSS_similarity.svg \
  --output-pdf /path/to/run/reports/SERPINE1_with_TSS_similarity.pdf \
  --output-png /path/to/run/reports/SERPINE1_with_TSS_similarity.png \
  --renderer rsvg-convert
```

Preparation and the tall-SVG path do not require Matplotlib; only the compact
PDF renderer does. Focused offline regression tests, including synthetic source
mismatches, both axis directions and dense-footer geometry, run with:

```bash
python3 -m unittest \
  scripts.test_tss_regulatory_integrated_report \
  scripts.test_tp73_locus_cutrun_lane_validation
```

## Interpretation boundary

Sequence recurrence is structural evidence useful for reporter contrasts.
Ensembl feature classes, predicted TFBS, and CUT&RUN enrichment do not by
themselves establish direct binding, reporter activity, dispensability, or
sufficiency. See `interpretation.md` for observations generated from this exact
comparison rather than a hand-written biological conclusion.
