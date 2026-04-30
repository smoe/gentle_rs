# TFBS Similarity Ranking Tutorial

This is a short sign-off tutorial for the new DNA-window `TFBS similarity`
path.

The goal is practical, not theoretical:

- verify that the GUI can rank candidate TF motifs against one anchor motif
  over the same DNA span
- verify that the cached ranked table is inspectable in-window
- verify that the same shared report can be exported as JSON
- verify that the same biology can be replayed through one small offline
  workflow and one ClawBio wrapper request

This tutorial intentionally reuses the same tiny synthetic FASTA from the
stateless direct-inspection walkthrough so the sign-off path stays local and
fast.

## What You Will Test

By the end of this tutorial, you should have verified all of these:

- the DNA-window toolbar `TFBS similarity` action runs without needing a
  promoter-design specialist window
- the `TFBS similarity ranking` subpanel in `TFBS annotation` shows the cached
  ranked report
- the ranked table includes the expected anchor, metric, and candidate count
- `Export cached TFBS similarity JSON...` replays the same shared engine route
- the same ranking can be reproduced by one offline workflow example
- the same workflow can be wrapped for ClawBio without first building a GENtle
  state by hand

## Synthetic Input

Local GUI input:

- [`docs/tutorial/inputs/inline_sequence_inspection_demo.fa`](./inputs/inline_sequence_inspection_demo.fa)

Exact sequence:

```text
GAATTCCCGGGATCCGGGCGGGGCGCATGTGTAACAGGGGCGGGGC
```

Why this sequence is still good enough for the ranking tutorial:

- it is tiny, so repeated GUI reruns stay fast
- it contains one GC-rich `SP1`-like block plus one p53-family-like teaching
  block
- the ranking route still has a continuous signal to compare even when the
  sequence is much shorter than a real promoter

## GUI Walkthrough

1. Open the FASTA:
   [`docs/tutorial/inputs/inline_sequence_inspection_demo.fa`](./inputs/inline_sequence_inspection_demo.fa)
2. In the DNA window, open the `TFBS annotation (log-likelihood ratio)` panel.
3. Configure the motif set.
   - Easiest robust route:
     - use the picker list under `JASPAR filter`
     - add these motifs with the `+` button:
       - `SP1`
       - `TP53`
       - `TP63`
       - `TP73`
       - `REST`
       - `CTCF`
   - If those exact names are already accepted directly by your local motif
     snapshot, the quick text route is:
     - `Selected motifs = SP1,TP53,TP63,TP73,REST,CTCF`
4. In the score-track subsection, use:
   - `value kind = llr_background_tail_log10`
   - `clip negatives = off`
5. In the new `TFBS similarity ranking (shared report)` subsection, use:
   - `anchor motif = SP1`
   - `all JASPAR motifs = off`
   - `candidate motifs = TP53,TP63,TP73,REST,CTCF`
   - `metric = Smoothed Spearman rho`
   - `limit = 10`
   - leave `species filters` empty for the first run
   - leave `include cached remote metadata` off for the first run
6. Use the toolbar menu `TFBS similarity -> Rank similarity in whole sequence`.
7. Return to the `TFBS annotation` panel.
   - Expected result:
     - the `TFBS similarity ranking` block is no longer empty
     - the cached table shows `anchor SP1`
     - the table reports `5` requested candidates or fewer returned rows if
       your local motif snapshot resolves some names differently
     - `SP1` itself must not appear as a candidate row
8. Press `Export cached TFBS similarity JSON...`.
   - Expected result:
     - a JSON file is written through the shared
       `SummarizeTfbsTrackSimilarity` route
9. Optional metadata/sign-off branch:
   - if you already have a cached JASPAR remote-metadata snapshot locally,
     rerun once with:
     - `include cached remote metadata = on`
     - `species filters = Homo sapiens`
   - Expected result:
     - the table still loads
     - metadata cells become richer when cached rows exist
     - any filter caveat is shown explicitly as a warning or status note

What this should prove:

- the first GUI slice for TFBS similarity ranking is real and inspectable
- the GUI is still thin: it delegates ranking to the shared engine operation
- pre-rank species filtering is visible to the user instead of being hidden

## Offline Workflow Replay

Canonical offline workflow example:

- [`docs/examples/workflows/tfbs_track_similarity_stateless_offline.json`](../examples/workflows/tfbs_track_similarity_stateless_offline.json)

From the repository root:

```sh
cargo run --quiet --bin gentle_cli -- \
  --state /tmp/tfbs_track_similarity_demo.state.json \
  workflow @docs/examples/workflows/tfbs_track_similarity_stateless_offline.json
```

Expected workflow artifacts:

- `artifacts/tfbs_track_similarity_demo.score_tracks.json`
- `artifacts/tfbs_track_similarity_demo.score_tracks.svg`
- `artifacts/tfbs_track_similarity_demo.similarity.json`

This path stays fully offline and state-optional at the biology layer: the
workflow uses the same inline sequence letters directly rather than requiring
you to create or save a project first.

## ClawBio Replay

Matching ClawBio request:

- [`integrations/clawbio/skills/gentle-cloning/examples/request_workflow_tfbs_track_similarity_stateless.json`](../../integrations/clawbio/skills/gentle-cloning/examples/request_workflow_tfbs_track_similarity_stateless.json)

Example wrapper invocation from a ClawBio checkout:

```sh
python clawbio.py run gentle-cloning \
  --input skills/gentle-cloning/examples/request_workflow_tfbs_track_similarity_stateless.json \
  --output /tmp/gentle_clawbio_tfbs_similarity
```

Expected result:

- the wrapper writes the same three workflow artifacts into the output bundle
- no hand-prepared GENtle state is required before the request

## What To Mark As Successful

Mark this tutorial successful if all of these are true:

- the GUI `TFBS similarity` menu runs for the whole sequence without errors
- the `TFBS similarity ranking` subpanel shows a populated cached ranked table
- the cached table names `SP1` as the anchor and does not include `SP1` as a
  candidate row
- JSON export works from the cached report
- the offline workflow writes the three expected artifacts
- the optional species-filter rerun either:
  - works with cached remote metadata, or
  - fails gracefully with explicit metadata/filter warnings rather than hidden
    behavior
