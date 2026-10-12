# CUT&RUN Motif-Pilot Region Selection Runbook

`scripts/cutandrun_select_regions.py` turns a candidate BED/bedGraph and one
BigWig signal track into a deterministic set of centered genomic sequences for
motif-discovery pilots. It is a standalone analysis helper, not a GENtle engine
operation. It does not change project state or invent candidate intervals from
the BigWig.

## Reproducibility Contract

All coordinates are **0-based half-open** at every stage: candidate intervals,
BigWig requests, merged islands, centered windows, output BED files, FASTA
slices, and FASTA headers. A header such as `>2:100-300` therefore describes
bases in `[100, 300)`. No 1-based conversion is performed.

The pipeline order is fixed:

1. Filter candidate segments by chromosome and score.
2. Merge same-chromosome segments whose gap is at most `--merge-distance`.
3. Score each merged island from the BigWig.
4. Rank islands in a deterministic total order.
5. Select, center, clip, and extract sequences.

BigWig scoring uses both maximum and mean signal through
`pyBigWig.stats(chrom, start, end, type=..., exact=True)`. The `exact=True`
argument is mandatory: pyBigWig's default may use zoom-level summaries and is
not suitable for this reproducibility contract. A `[None]` result, meaning no
coverage in the interval, becomes `0.0`.

The total-order ranking key is:

```text
bw_max descending,
bw_mean descending,
support descending,
length descending,
chrom ascending,
start ascending,
end ascending
```

Here, `support` is the sum of `candidate_score * candidate_length` across the
segments merged into an island, and `length` is `end - start`. Support and
length are tie-break fields; the exact BigWig statistics rank first.

## Chromosome And Assembly Policy

Chromosome names are normalized for allowlist matching by removing an optional
`chr` prefix and upper-casing the remainder. `M`, `MT`, and `chrM` all normalize
to mitochondrial `MT`. The concrete chromosome spelling is not rewritten.

The default `--include-chroms canonical` policy for this pilot keeps
chromosomes `2` through `22`, `X`, and `Y`. It intentionally excludes:

- chromosome 1, because the present TP73 pilot may be confounded there by
  adenoviral-vector artifact risk;
- mitochondrial sequence;
- alternative and unplaced contigs.

Use a comma list to select another normalized allowlist. `--exclude-chroms` is
applied after inclusion. Both the resolved policy and the concrete chromosome
names actually used are recorded in provenance.

Every candidate chromosome that passes the allowlist must occur under the same
concrete name and with the same chromosome length in both the BigWig header and
FASTA `.fai`. Retained intervals must also fit within that shared length. A
candidate named `chr2` therefore does not silently match a BigWig and FASTA that
call it `2`. Conflicting lengths or out-of-bounds intervals always stop the run.
For an exploratory run, `--allow-missing-chroms` warns and skips chromosomes
that are absent from either source; it does not reinterpret names or waive
conflicting assembly lengths.

The genome FASTA must use the same assembly and naming as the signal and
candidate tracks. For the example below, use Ensembl GRCh38 primary assembly
with bare chromosome names. The existence guard catches obvious name/header
mismatches, but it cannot prove biological assembly identity by itself.

## Threshold, Merge, And Window Definitions

Candidate score columns are 1-based at the command line. By default, the tool
uses column 4 for `.bedGraph`/`.bdg`, column 5 for BED with at least five
columns, and column 4 for BED4. Override this with `--score-column` when the
producer uses another layout.

`--min-score` applies an absolute candidate-score threshold.
`--score-percentile P` instead computes the nearest-rank percentile over the
scores of chromosome-filtered segments:

```text
sorted_scores[max(1, ceil(P / 100 * N)) - 1]
```

Segments at or above the resulting threshold are retained. Omitting both
threshold options retains all chromosome-filtered segments.

Filtered segments are sorted by chromosome and start. Overlapping,
book-ended, or nearby segments merge whenever the half-open gap is at most the
configured distance (default 50 bp).

For a merged island and requested width `W`, the centered window is:

```text
center = (start + end) // 2
half = W // 2
win_start = center - half
win_end = win_start + W
```

The result is then intersected with `[0, chromosome_length)`. A chromosome-end
window is shorter than `W`; it is never shifted or padded to restore the width.

## Requirements

Install pyBigWig in the Python environment used to run the script:

```bash
python3 -m pip install pyBigWig
```

The import is lazy, so `--help` and the pure selection helpers remain usable
without pyBigWig. The default FASTA extractor is an internal `.fai`-indexed
reader that accounts for FASTA newline bytes. `--fasta-tool bedtools` uses
`bedtools getfasta` as an alternative and therefore requires the validated index
at the conventional adjacent `<genome-fasta>.fai` path; bedtools cannot consume
a separate `--fai` path. Neither extractor changes sequence case, so soft
masking is preserved.

## TP73 SK-MEL-29 Pilot

Use the SK-MEL-29 `_2` subline for this example because the corresponding
`tp73_skmel29_2_TA_R1.clipped.clean.bed` candidate file is available locally.
Do not substitute the `_1` track while retaining the `_2` candidates.

The tracks under `data/resources/` are local analysis inputs and are not
committed:

```bash
python3 scripts/cutandrun_select_regions.py \
  --bigwig data/resources/cutandrun_20250602_noDuplicates/tp73_skmel29_2_TA_R1.bigWig \
  --candidates data/resources/cutandrun_20250602_noDuplicates/tp73_skmel29_2_TA_R1.clipped.clean.bed \
  --genome-fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --fai /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
  --output-dir out/tp73_skmel29_2_TA_R1_motif_pilot \
  --include-chroms canonical --window 200 --merge-distance 50 \
  --score-percentile 90 --top-n 500
```

This pilot ranks a **single BigWig signal track**. It performs no input,
negative-control, or matched-condition subtraction and must not be interpreted
as differential occupancy. Use the selected sequences as motif-discovery input,
then return to the underlying controlled experiment for biological inference.

## Outputs And Methods Record

With the default prefix `cutandrun_regions`, the output directory contains:

- `cutandrun_regions.filtered_segments.bed`: chromosome- and
  threshold-filtered candidates;
- `cutandrun_regions.merged_islands.bed`: merged candidate islands;
- `cutandrun_regions.scored_ranked.bed`: exact-scored islands in rank order;
- `cutandrun_regions.scored_ranked.columns.txt`: the nine ranked-file columns;
- `cutandrun_regions.topN.centered_wW.bed`: clipped centered windows;
- `cutandrun_regions.topN.fasta`: extracted sequences;
- `cutandrun_regions.provenance.json` and `.tsv`: machine-readable and flat
  methods records.

`N` in output names is the number actually selected, which can be lower than a
requested `--top-n` when fewer islands survive. BED files have no header row.
The provenance records every input's absolute path, byte size, and SHA-256;
the FASTA and `.fai`; chromosome policy; score/threshold, merge, window, and
selection settings; pyBigWig and optional bedtools versions; UTC timestamp;
best-effort Git commit; output paths; and exact invocation arguments.

For a custom provenance location, pass `--provenance PATH.json`; the flat TSV
is written next to it with a `.tsv` suffix.
