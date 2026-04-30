# AI Task Playbooks for GENtle

Purpose: ready-to-run patterns for local AI assistants. Each playbook specifies
inputs, expected outputs, and command families.

## 1) Prepare reference genome

Goal: install sequence + annotation + index once.

Required inputs:

- `genome_id`
- optional `catalog_path`
- optional `cache_dir`

Expected outputs:

- prepared genome status
- index/build progress summary
- clear error if annotation path/source invalid

Command pattern:

```bash
gentle_cli genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH]
gentle_cli genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]
```

## 2) Retrieve sequence by gene

Goal: extract one gene-centered region from prepared genome.

Required inputs:

- `genome_id`
- `query` (gene symbol/id)

Expected outputs:

- created sequence ID
- anchor/provenance metadata

Command pattern:

```bash
gentle_cli genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID]
```

## 3) Retrieve sequence by coordinates

Goal: extract exact genomic interval.

Required inputs:

- `genome_id`, `chromosome`, `start`, `end`

Expected outputs:

- created sequence ID and length check

Command pattern:

```bash
gentle_cli genomes extract-region GENOME_ID CHR START END [--output-id ID]
```

## 4) Generate candidate windows between two anchors

Goal: enumerate windows constrained to biologically relevant local boundaries.

Required inputs:

- `set_name`, `seq_id`
- window `length`
- anchor A and B (position or feature boundary)

Expected outputs:

- candidate count
- page preview and metric fields

Command pattern:

```bash
gentle_cli shell "candidates generate-between-anchors SET SEQ --length 20 --anchor-a-json '...' --anchor-b-json '...' --step 1 --limit 10000"
gentle_cli shell "candidates show SET --limit 50 --offset 0"
```

## 5) Add scoring metrics and optimize

Goal: transform raw windows into ranked subsets.

Required inputs:

- input set
- scoring expression/distance objective
- optimization strategy (`top-k`, `pareto`, quantile filter)

Expected outputs:

- derived metric(s)
- explicit result sets (`*_top`, `*_frontier`, `*_filtered`)

Command pattern:

```bash
gentle_cli shell "candidates score SET gc_balance '100*(gc_fraction-at_fraction)'"
gentle_cli shell "candidates score-distance SET dist_cds --feature-kind CDS"
gentle_cli shell "candidates score-weighted SET priority --term gc_fraction:0.7:max --term dist_cds:0.3:min --normalize"
gentle_cli shell "candidates top-k SET SET_top --metric priority --k 25 --direction max"
```

## 6) BLAST specificity check

Goal: verify candidate sequence uniqueness/off-target profile.

Required inputs:

- target catalog (`genomes` or `helpers`)
- target id (`genome_id` or `helper_id`)
- query sequence

Expected outputs:

- hit table summary
- top hit identity/alignment stats

Command pattern:

```bash
gentle_cli genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn]
gentle_cli helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn]
```

## 7) Import external experiment tracks

Goal: project BED/BigWig/VCF evidence onto anchored sequence.

Required inputs:

- `seq_id` (anchored to matching assembly)
- file path(s)

Expected outputs:

- new `TRACK:*` feature groups
- import summary with mapped/unmapped records

Command pattern:

```bash
gentle_cli shell "tracks import-bed SEQ_ID PATH.bed[.gz]"
gentle_cli shell "tracks import-bigwig SEQ_ID PATH.bw"
gentle_cli shell "tracks import-vcf SEQ_ID PATH.vcf[.gz]"
```

## 8) Build reproducible macro/template

Goal: package a repeated optimization script.

Required inputs:

- template name
- script text
- optional parameter names/defaults

Expected outputs:

- template stored in project metadata
- callable with bound values

Command pattern:

```bash
gentle_cli shell "candidates template-put NAME --script @script.gsh --description '...' --param seq_id --param min_gc=0.40"
gentle_cli shell "candidates template-run NAME --bind seq_id=grch38_tp53 --bind min_gc=0.45 --transactional"
```

## 9) Agent-assisted planning (chat/ask/auto)

Goal: use AI for plan + commands while preserving deterministic execution.

Required inputs:

- `system_id`
- prompt with objective + constraints + safety mode

Expected outputs:

- assistant explanation
- follow-up questions
- suggested commands with execution intent

Command pattern:

```bash
gentle_cli agents ask SYSTEM_ID --prompt "..." [--catalog PATH] [--base-url URL] [--model MODEL] [--no-state-summary]
```

## 10) Recommended answer format for AI

When running these playbooks, ask the AI to answer in this structure:

1. `Assumptions`
2. `Missing inputs`
3. `Plan`
4. `Commands`
5. `Validation checks`
6. `Stop point requiring user confirmation`
