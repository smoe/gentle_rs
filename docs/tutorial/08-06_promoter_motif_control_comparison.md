# Promoter Motif Control-Set Comparison Tutorial

> Type: `GUI/CLI walkthrough`
> Status: `manual/hybrid`
> Goal: inspect promoter regions for shared motifs, then compare a foreground
> cohort against a matched control set to find candidate over- and
> under-represented motif signals without turning sequence evidence into a
> biological verdict.

This tutorial is for the question that comes after a single promoter looks
interesting:

- Which transcription-factor motifs are common across the foreground promoters?
- Which motifs or motif pairs appear stronger in the foreground than in a
  matched control set?
- Which expected motifs are weaker or absent in the foreground compared with
  controls?
- Are the promoter score traces consistent with a declared relationship such
  as co-regulation or anti-co-regulation?

GENtle can help structure that review. It does not prove regulation from motif
matches alone. Treat every result below as promoter sequence evidence that
still needs expression, occupancy, perturbation, or reporter validation.

## Where This Fits

Use the executable TP73 promoter artifact slice first if you want a completely
offline warm-up:

- [`docs/tutorial/generated/chapters/08-03_promoter_design_artifact_slice_offline.md`](./generated/chapters/08-03_promoter_design_artifact_slice_offline.md)

That chapter teaches single-locus promoter windows, expression evidence,
TFBS score tracks, and similarity rankings. This page teaches the next layer:
foreground-versus-control interpretation across multiple promoters.

## Required Inputs

For a real cohort run you need:

- a prepared genome catalog entry, for example `Human GRCh38 Ensembl 116`
- one foreground gene list
- one matched control gene list
- one motif panel or motif group
- matching promoter window sizes for foreground and control

Control sets should be chosen before looking at the answer. Good controls are
usually matched by organism, gene-model source, promoter-window size, and a
biological property relevant to the question, such as expression range or cell
type. Random genes can be useful for a rough screen, but they are not a
substitute for a defensible matched control set.

Example foreground for a TP73/PATZ1 regulatory question:

```text
TP73
PATZ1
TP53
TP63
```

Example control set for a first technical pass:

```text
E2F1
SP1
TERT
FUS
TARDBP
IL6
IL10
```

That example control set is intentionally not a final biological control. It
is a compact tutorial set with familiar genes and mixed promoter architectures.
For publication-grade work, replace it with matched controls defined in your
study design.

## Step 1: Inspect One Promoter Region In The GUI

Start with one target gene or variant context.

1. Open or retrieve the promoter-containing genomic sequence.
2. Open the DNA window and select the gene, transcript, variant, or promoter
   annotation that anchors the question.
3. Open `Promoter design`.
4. Run:
   - `Annotate promoter windows`
   - `Build evidence matrix`
   - `Show TF score tracks`
   - `Show TFBS similarity ranking`
5. Inspect the score-track graphic and evidence matrix before running a
   cohort. Confirm that the TSS, strand, promoter window, and motif panel are
   plausible for the biological question.

What to look for:

- a TSS marker in the expected place
- strand-aware upstream/downstream orientation
- motifs whose score traces have visible positive support
- motifs whose traces are flat or negative-only
- CUT&RUN, repeat, variant, or expression rows that explain why a promoter
  should be treated cautiously

The GUI is a thin adapter here. The same artifacts can be replayed from the
shared shell commands shown below.

## Step 2: Export A Foreground Promoter TFBS Summary

Run the same motif panel across the foreground genes. Use exactly the same
window size you intend to use for the control cohort.

```bash
mkdir -p artifacts/promoter_motif_controls

cargo run --quiet --bin gentle_cli -- shell \
  'genomes promoter-tfbs-summary "Human GRCh38 Ensembl 116" --gene TP73 --gene PATZ1 --gene TP53 --gene TP63 --motif SP1 --motif TP53 --motif TP63 --motif TP73 --motif PATZ1 --motif E2F1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 --path artifacts/promoter_motif_controls/foreground.promoter_tfbs.json'
```

Optional SVG overview:

```bash
cargo run --quiet --bin gentle_cli -- shell \
  'genomes promoter-tfbs-svg "Human GRCh38 Ensembl 116" --gene TP73 --gene PATZ1 --gene TP53 --gene TP63 --motif SP1 --motif TP53 --motif TP63 --motif TP73 --motif PATZ1 --motif E2F1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 artifacts/promoter_motif_controls/foreground.promoter_tfbs.svg'
```

The JSON report is `gentle.multi_gene_promoter_tfbs.v1`. The most useful first
table is `summary_rows[]`, with one row per resolved gene and motif:

- `gene_label`
- `transcript_id`
- `tf_id` and `tf_name`
- `max_score`
- `positive_fraction`
- peak position relative to the promoter/TSS when available

## Step 3: Export The Matched Control Summary

Run the control set with the same motif panel, score kind, and promoter window.

```bash
cargo run --quiet --bin gentle_cli -- shell \
  'genomes promoter-tfbs-summary "Human GRCh38 Ensembl 116" --gene E2F1 --gene SP1 --gene TERT --gene FUS --gene TARDBP --gene IL6 --gene IL10 --motif SP1 --motif TP53 --motif TP63 --motif TP73 --motif PATZ1 --motif E2F1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 --path artifacts/promoter_motif_controls/control.promoter_tfbs.json'
```

Now compare the same fields between foreground and control.

Candidate over-representation means:

- the foreground has more promoters with positive support for a motif, or
- the foreground has higher `max_score`/`positive_fraction` for that motif, or
- the motif appears in shared peaks across more foreground promoters than in
  controls.

Candidate under-representation means the inverse:

- a motif expected from the biology is weaker, absent, or less frequent in the
  foreground than in controls.

This is a screen, not a statistical test. Use it to decide what to inspect and
validate next.

## Step 4: Make The Comparison Table Human-Readable

The structured reports are JSON so they can be re-used by GUI, CLI, agents, or
notebooks. For a quick terminal table:

```bash
jq -r '
  .summary_rows[]
  | [.tf_id, .gene_label, .max_score, .positive_fraction,
     (.peak_position_promoter_relative_bp // "")]
  | @tsv
' artifacts/promoter_motif_controls/foreground.promoter_tfbs.json \
  > artifacts/promoter_motif_controls/foreground.promoter_tfbs.tsv

jq -r '
  .summary_rows[]
  | [.tf_id, .gene_label, .max_score, .positive_fraction,
     (.peak_position_promoter_relative_bp // "")]
  | @tsv
' artifacts/promoter_motif_controls/control.promoter_tfbs.json \
  > artifacts/promoter_motif_controls/control.promoter_tfbs.tsv
```

For a compact motif-level view:

```bash
jq -r '
  .summary_rows
  | group_by(.tf_id)[]
  | [.[0].tf_id,
     length,
     (map(select(.positive_fraction > 0)) | length),
     (map(.max_score) | max),
     ((map(.positive_fraction) | add) / length)]
  | @tsv
' artifacts/promoter_motif_controls/foreground.promoter_tfbs.json
```

Read the columns as:

1. motif id
2. resolved promoter count for the motif
3. promoters with positive support
4. maximum score in the cohort
5. mean positive-support fraction

Run the same command on the control report, then compare the rows. Motifs with
foreground values above the control are candidate foreground-enriched signals.
Motifs below the control are candidate depleted or under-represented signals.

## Step 5: Review Motif Combinations

GENtle does not yet ship a formal motif-combination enrichment statistic. The
safe first-pass tutorial definition of a motif combination is therefore:

> two or more motifs with positive support in the same promoter window, reviewed
> under the same score kind and promoter span.

Create a simple co-occurrence table from the foreground report:

```bash
jq -r '
  .summary_rows
  | group_by(.gene_label)[]
  | {
      gene: .[0].gene_label,
      positive: (map(select(.positive_fraction > 0) | .tf_id) | unique)
    }
  | select(.positive | length >= 2)
  | [.gene, (.positive | join("+"))]
  | @tsv
' artifacts/promoter_motif_controls/foreground.promoter_tfbs.json
```

Repeat it for the control report. Candidate over-represented combinations are
motif pairs or groups that appear in more foreground promoters than control
promoters. Candidate under-represented combinations are expected pairs that
appear less often in the foreground than in controls.

For real analysis, refine this coarse screen with:

- a fixed positive-support threshold chosen before inspection
- distance between motif peaks
- motif orientation if the hypothesis depends on orientation
- occupancy or expression evidence
- a larger, matched control set

## Step 6: Use Relationship-Aware Cohort Comparison

When the foreground genes have a declared relationship, ask GENtle to carry
that expectation explicitly. This is not a verdict; it changes which review
cues GENtle surfaces.

For a co-regulated hypothesis:

```bash
cargo run --quiet --bin gentle_cli -- shell \
  'genomes promoter-cohort-comparison "Human GRCh38 Ensembl 116" --cohort-label tp73_axis_coregulated --cohort-kind co_regulated --gene TP73 --gene PATZ1 --gene TP53 --gene TP63 --motif SP1 --motif TP53 --motif TP63 --motif TP73 --motif PATZ1 --motif E2F1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 --path artifacts/promoter_motif_controls/foreground.coregulated_cohort.json'
```

For an anti-co-regulated hypothesis:

```bash
cargo run --quiet --bin gentle_cli -- shell \
  'genomes promoter-cohort-comparison "Human GRCh38 Ensembl 116" --cohort-label tp73_axis_anti_coregulated --cohort-kind anti_co_regulated --gene TP73 --gene PATZ1 --gene TP53 --gene TP63 --motif SP1 --motif TP53 --motif TP63 --motif TP73 --motif PATZ1 --motif E2F1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 --path artifacts/promoter_motif_controls/foreground.anti_coregulated_cohort.json'
```

Inspect these report fields:

- `shared_tfbs_peaks[]`
  - motifs with support across several foreground promoters
- `cohort_specific_tfbs_peaks[]`
  - motifs that look stronger or more specific in a subset
- `pairwise_similarity[]`
  - motif-score similarity between promoter pairs
- `relationship_flags[]`
  - review cues for unexpected divergence in co-regulated cohorts or
    unexpected concordance in anti-co-regulated cohorts
- `warnings[]`
  - unresolved genes, missing transcripts, or cautious interpretation notes

For `manual` or unspecified cohorts, GENtle does not add relationship-based
flags.

## Step 7: Add Occupancy Or Expression Evidence When Available

Promoter motif comparisons are stronger when paired with independent evidence.
Use the existing promoter design and cohort routes to keep these evidence
layers side by side:

- expression rows with `features promoter-expression-evidence`
- CUT&RUN dataset/read-report references with
  `genomes promoter-cohort-comparison --cutrun-dataset-id ...`
- single-locus promoter evidence with
  `features promoter-evidence-matrix`

Keep the wording conservative:

- a motif plus occupancy peak is support for follow-up, not proof
- expression association is not causality
- absence of a motif in one model does not rule out indirect regulation
- under-representation can be biologically meaningful, especially for
  anti-co-regulated or repressed promoter sets

## What To Mark As Successful

Mark this tutorial successful if all of these are true:

- one foreground promoter TFBS summary is written
- one control promoter TFBS summary is written with the same motif panel and
  promoter window
- you can name at least one candidate motif that is stronger in foreground
  than controls
- you can name at least one candidate motif that is weaker in foreground than
  controls, or state that none were observed under the chosen threshold
- you can list motif combinations that co-occur in foreground promoters and
  compare their presence against the control set
- relationship-aware cohort comparison emits no verdict, only shared/specific
  peaks, pairwise similarity, optional relationship flags, and warnings

## Common Pitfalls

- Do not compare foreground and control reports generated with different
  window sizes or score kinds.
- Do not treat a motif match as direct transcription-factor binding.
- Do not use a random control set as a publication-grade negative control
  unless the analysis plan justifies it.
- Do not hide under-represented motifs; depletion can be just as informative
  as enrichment.
- Do not mix plus- and minus-strand promoter windows by raw genomic direction.
  GENtle's promoter routes score in transcription-aligned orientation; preserve
  that convention when exporting rows.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context copied from GENtle Help -> Tutorial -> Copy Feedback Context.

- Tutorial title:
- Tutorial/chapter id:
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
