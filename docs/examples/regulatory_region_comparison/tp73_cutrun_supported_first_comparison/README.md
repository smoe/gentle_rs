# TP73 CUT&RUN-supported promoterome comparison

**Historical evidence, not a corrected result.** The candidate bundle, report,
tables, SVG/PNG and receipt below are preserved byte-for-byte from `b4ed17d4`.
Review found unchecked reference assumptions, hardcoded report conclusions,
incorrect reverse-block break markers, inconsistent query-gene exclusion and
missing HSP-cap accounting in that pipeline. Do not use its counts, absence
statements or red outlines as validated design evidence. The corrected scripts
require a new run; passing the retained-file hash tests is not scientific
acceptance. The historical PNG's conversion environment was not recorded.

This retained first comparison binds TP73 CUT&RUN BigWig signal to exact
GRCh38/Ensembl-116 transcript promoter windows for the selected `CD44`,
`TGFB1`, and `SERPINE1` transcript models. Seven distinct −2,000/+200 bp TSS
windows passed the declared matched-control rule and were compared with the
389,722-window human promoterome prepared by
`scripts/prepare_transcript_promoterome.py`.

The support rule is deliberately inspectable: at least one TAp73alpha or
DNp73beta mean signal across the exact window must exceed the matched GFP mean
from the same cell line. The retained candidate JSON records every lane mean,
delta, source identifier, and BigWig digest. This is a positive window-mean
difference screen, not a peak/significance test or proof of direct binding or
promoter activity. Assembly and compatible signal normalization must be
established for the six supplied tracks; filenames and file hashes alone
cannot establish either biological property.

## Produce A Corrected Run

Use Python with `pyBigWig`, the original six signal tracks, and the complete
prepared promoterome including its receipt. The scripts do not download data.
Preparation verifies the reference schema, genome, orientation and input hashes,
requires every requested transcript/gene mapping, and checks selected sequence
lengths and strand/TSS geometry. `anchor_verified` means consistent with that
receipt-bound reference, not independently revalidated annotation. Empty or
partial transcript selections and clipped windows are errors, not silent drops.
The current checkout SHA and producer-file digest are recorded separately.

Prepare the candidate bundle from the promoterome and the six local
E-MTAB-15709 BigWigs:

```bash
python3 scripts/prepare_tp73_cutrun_promoter_candidates.py \
  --promoterome /path/to/human-grch38-ensembl116-promoterome \
  --track-root /path/to/cutandrun_20250602_noDuplicates \
  --output /path/to/run \
  --source-revision "$(git rev-parse HEAD)"
```

Run the declared local-alignment comparison:

```bash
python3 scripts/compare_candidates_to_promoterome.py \
  --candidates-json /path/to/run/candidate_regions.json \
  --query-fasta /path/to/run/candidate_regions.fa \
  --promoterome /path/to/human-grch38-ensembl116-promoterome \
  --output /path/to/run/blastn-40bp-80pct \
  --task blastn --min-alignment-bp 40 --min-identity-pct 80 \
  --max-evalue 1e-5 --max-target-seqs 1000000 --max-hsps 10
```

Then render into a separate, new/empty directory. PNG output is explicit and
uses [rsvg-convert](https://manpages.debian.org/trixie/librsvg2-bin/rsvg-convert.1.en.html).
Its executable digest, version, invocation and PNG digest are recorded in the
new receipt. Pin the renderer environment, including fonts and libraries, when
comparing PNG bytes across hosts. Omit `--png-renderer` for SVG/TSV/Markdown
only; the receipt explicitly records that PNG was not requested.

```bash
python3 scripts/render_tp73_cutrun_promoter_comparison.py \
  --root /path/to/run --task-directory blastn-40bp-80pct --top 5 \
  --output /path/to/corrected-figure \
  --png-renderer rsvg-convert
```

The historical 376 MiB raw HSP table and 98 MiB target table are intentionally
not versioned. They are regenerable from the bound inputs. The repository
retains the candidates, compact summary/top-five table, report, SVG/PNG figure,
and receipt. The historical `receipt.json` binds the omitted full comparison
report by digest. New receipts additionally bind both raw HSP and target tables.
The renderer rejects mixed inputs and old comparison policies; merely adding a
new policy label to an old report does not recreate the missing analysis.

For Glen's offline regression gate (no BLAST/BigWig dependency required):

```bash
python3 -B -m unittest scripts.test_tp73_cutrun_promoter_comparison \
  scripts.test_compare_candidates_to_promoterome \
  scripts.test_prepare_transcript_promoterome
```

These tests generate hand-crafted temporary DNA, mappings, mock signals and
HSPs to exercise the scripts. They are not a replay of E-MTAB-15709. Glen still
needs to regenerate the real comparison and inspect the new figure before
publishing a replacement evidence bundle.

## Interpretation

In corrected runs, blue frequency strips show how many distinct *other genes* have a qualifying
hit at each candidate coordinate. Detailed rows are distinct genomic promoter
windows, not transcripts; all genes and transcripts sharing a TSS remain
attached to that one occurrence. Windows overlapping the query or sharing any
resolved query gene ID are excluded consistently from strips, counts, coverage
tiers and detailed rows, including shared windows containing that gene and
another gene. Numbers record block order in the target promoter. Decreasing
query coordinates are expected for a consistently reverse-oriented alignment;
red outlines mark changes relative to alignment direction, not proof of a
biological rearrangement.

Target and HSP caps are assessed on raw output before filtering. Reaching
either cap marks the figure, report, tables and receipt as **lower bounds**.
No cap saturation means only that the configured output limit was not observed;
it does not prove exhaustive homology or whole-genome uniqueness. A zero-hit
query still produces a report and an empty, header-bearing top-match table.

The historical first report suggested a broadly recurrent component in the distal
`SERPINE1` TSS window, including MAN2A2 promoter matches covering about 47% of
the 2.2 kb candidate. `CD44` and `TGFB1` show only sparse short recurrence at
this threshold. Those statements require the corrected replay above before use
in reporter decisions; they do not establish fragment sufficiency or regulatory
function. New report observations and parameters are generated from the actual
supplied comparison rather than repeating this historical narrative.
