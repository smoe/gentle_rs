# AI Cloning Support Examples (GENtle)

Purpose: few-shot examples for local AI systems so users get actionable,
GENtle-native responses.

## Example 1: candidate generation between anchors

User prompt:

```text
Objective:
Find 20 bp candidate windows between TP53 start and TP53 end and rank by GC and distance from CDS.

Inputs:
- seq_id: grch38_tp53
- anchor A: TP53 gene start
- anchor B: TP53 gene end

Constraints:
- GC between 40 and 80%
- output top 25

Output wanted:
Exact gentle_cli shell commands and a validation checklist.

Execution policy:
ask-before-run
```

Good assistant response shape:

1. Assumptions
2. Commands:

```bash
gentle_cli shell "candidates generate-between-anchors tp53_seed grch38_tp53 --length 20 --anchor-a-json '{\"FeatureBoundary\":{\"feature_kind\":\"gene\",\"feature_label\":\"TP53\",\"boundary\":\"Start\",\"occurrence\":0}}' --anchor-b-json '{\"FeatureBoundary\":{\"feature_kind\":\"gene\",\"feature_label\":\"TP53\",\"boundary\":\"End\",\"occurrence\":0}}' --step 1 --limit 10000"
gentle_cli shell "candidates score tp53_seed gc_balance '100*(gc_fraction-at_fraction)'"
gentle_cli shell "candidates score-distance tp53_seed dist_cds --feature-kind CDS"
gentle_cli shell "candidates score-weighted tp53_seed priority --term gc_fraction:0.6:max --term dist_cds:0.4:min --normalize"
gentle_cli shell "candidates filter tp53_seed tp53_gc_ok --metric gc_fraction --min 0.40 --max 0.80"
gentle_cli shell "candidates top-k tp53_gc_ok tp53_top25 --metric priority --k 25 --direction max"
gentle_cli shell "candidates show tp53_top25 --limit 25 --offset 0"
```

3. Validation checks
4. Stop point for user confirmation

## Example 2: BLAST check against prepared reference genome

User prompt:

```text
Check whether this candidate has too many off-targets in human GRCh38.
Sequence: ACGT...
Return exact command and how to interpret top hits.
```

Good assistant response excerpt:

```bash
gentle_cli genomes blast homo_sapiens_grch38_ensembl_116 ACGT... --max-hits 20 --task blastn-short
```

Interpretation checklist:

- number of high-identity hits
- alignment length vs full query length
- chromosomal spread of hits
- mismatch/gap profile in top hits

## Example 3: import BED peaks and intersect with candidate set

User prompt:

```text
Map peaks from /data/chip/peaks.bed.gz to anchored sequence tp53_region, then keep only candidates near these peaks.
```

Good assistant response excerpt:

```bash
gentle_cli shell "tracks import-bed tp53_region /data/chip/peaks.bed.gz"
gentle_cli shell "candidates score-distance tp53_seed dist_track --feature-kind TRACK:peaks"
gentle_cli shell "candidates filter tp53_seed tp53_near_peaks --metric dist_track --max 200"
```

## Example 4: use agent in chat-only mode

User prompt:

```text
Give me a plan first. Do not execute any command. Ask me for missing values.
Task: design candidate windows around gene X for primer preselection.
```

Good assistant behavior:

- asks for exact sequence/gene IDs
- asks for window length and GC constraints
- returns command skeletons but no auto-execution intent

## Example 5: template-driven repetition

User prompt:

```text
Create a reusable candidate template for 20 bp windows with GC 40-80% and top-k 50.
```

Good assistant response excerpt:

```bash
gentle_cli shell "candidates template-put default_20bp --file ./templates/default_20bp.gsh --description '20bp default optimization' --param seq_id --param out_prefix=run"
gentle_cli shell "candidates template-run default_20bp --bind seq_id=grch38_tp53 --bind out_prefix=tp53 --transactional"
```

## Example 6: local model setup prompt

User prompt:

```text
My local endpoint is http://localhost:11964. Discover models and pick one available model; if none are found, ask me to confirm the endpoint.
```

Good assistant behavior:

- reminds user to set base URL in Agent Assistant
- requires explicit model selection when model is `unspecified`
- does not assume hidden fallback endpoints
