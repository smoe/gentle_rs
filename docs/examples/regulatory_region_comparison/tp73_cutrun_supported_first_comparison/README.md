# TP73 CUT&RUN-supported promoterome comparison

This retained first comparison binds TP73 CUT&RUN BigWig signal to exact
GRCh38/Ensembl-116 transcript promoter windows for the selected `CD44`,
`TGFB1`, and `SERPINE1` transcript models. Seven distinct −2,000/+200 bp TSS
windows passed the declared matched-control rule and were compared with the
389,722-window human promoterome prepared by
`scripts/prepare_transcript_promoterome.py`.

The support rule is deliberately inspectable: at least one TAp73alpha or
DNp73beta mean signal across the exact window must exceed the matched GFP mean
from the same cell line. The retained candidate JSON records every lane mean,
delta, source identifier, and BigWig digest. This is occupancy/enrichment
support, not proof of direct binding or promoter activity.

## Reproduce

Prepare the candidate bundle from the promoterome and the six local
E-MTAB-15709 BigWigs:

```bash
python3 scripts/prepare_tp73_cutrun_promoter_candidates.py \
  --promoterome /path/to/human-grch38-ensembl116-promoterome \
  --track-root /path/to/cutandrun_20250602_noDuplicates \
  --output /path/to/run \
  --source-revision c23ea8110cfbad1cd0dc09ffcc638ddf12609d07
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

Then render the compact figure and tables:

```bash
python3 scripts/render_tp73_cutrun_promoter_comparison.py \
  --root /path/to/run --task-directory blastn-40bp-80pct --top 5
```

The 376 MiB raw HSP table and 98 MiB complete target table are intentionally
not versioned. They are regenerable from the bound inputs. The repository
retains the candidates, compact summary/top-five table, report, SVG/PNG figure,
and receipt. `receipt.json` binds the omitted full comparison report by digest.

## Interpretation

Blue frequency strips show how many distinct *other genes* have a qualifying
hit at each candidate coordinate. Detailed rows are distinct genomic promoter
windows, not transcripts; all genes and transcripts sharing a TSS remain
attached to that one occurrence. Numbers record block order in the target
promoter, and red outlines mark an order/orientation break.

The first result identifies a broadly recurrent component in the distal
`SERPINE1` TSS window, including MAN2A2 promoter matches covering about 47% of
the 2.2 kb candidate. `CD44` and `TGFB1` show only sparse short recurrence at
this threshold. These observations motivate reporter split/deletion contrasts;
they do not establish fragment sufficiency or regulatory function.
