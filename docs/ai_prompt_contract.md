# AI Prompt Contract for GENtle Agent Systems

Purpose: define a robust request/response style for local or remote AI systems
used from GENtle.

## User-side guidance (what to send)

Always include these blocks in your prompt:

1. `Objective`
2. `Context`
3. `Inputs`
4. `Constraints`
5. `Output wanted`
6. `Execution policy`

Minimal template:

```text
Objective:
<what decision or artifact you need>

Context:
<I am working on sequence/genome ID(s)...>

Inputs:
- seq_id/genome_id/helper_id: ...
- anchors or coordinates: ...
- feature kinds/labels: ...

Constraints:
- length: ...
- GC range: ...
- motifs/restriction sites to require/avoid: ...
- strand or boundary rules: ...

Output wanted:
- plan + exact gentle_cli commands
- include verification checks

Execution policy:
- chat-only | ask-before-run | allow-auto-exec
```

## Agent-side behavior rules

The assistant should:

- prefer deterministic GENtle shell/CLI commands over free-form advice
- explicitly list assumptions
- ask questions when required fields are missing
- avoid claiming biology facts not derivable from user/project data
- never imply wet-lab success; frame as in-silico candidate generation
- keep risky actions in `ask` mode unless user explicitly allows auto-exec

The assistant should not:

- confuse local sequence anchors with genome provenance anchoring
- silently choose unavailable systems/models/endpoints
- output destructive commands without explicit rationale and confirmation

## Required response shape for usefulness

Use this human-readable structure:

1. `Assumptions`
2. `Questions` (if missing inputs)
3. `Plan`
4. `Commands`
5. `Expected results`
6. `Validation`
7. `Next decision`

## Copy/paste request patterns

### A) Design candidate windows near a feature

```text
Objective:
Generate candidate DNA windows near a feature and rank them.

Context:
Project sequence: grch38_tp53

Inputs:
- seq_id: grch38_tp53
- anchor A: TP53 gene start
- anchor B: TP53 gene end

Constraints:
- length 20 bp
- GC 40-80%
- maximize distance from CDS boundaries
- return top 25

Output wanted:
- exact gentle_cli shell commands
- include top-k and validation commands

Execution policy:
ask-before-run
```

### B) Check specificity with BLAST

```text
Objective:
Check off-target risk of one candidate sequence.

Inputs:
- genome_id: homo_sapiens_grch38_ensembl_116
- query_sequence: ACGT...

Constraints:
- max hits 20
- task blastn-short

Output wanted:
- exact genomes blast command
- interpretation checklist for top hits

Execution policy:
chat-only
```

### C) Import and use external peak evidence

```text
Objective:
Overlay BED peaks on an anchored sequence and prioritize candidates near peaks.

Inputs:
- seq_id: my_anchored_region
- bed_path: /path/to/peaks.bed.gz

Constraints:
- keep separate TRACK group naming
- do not alter original sequence

Output wanted:
- import commands
- candidate generation + filter commands near imported TRACK features

Execution policy:
ask-before-run
```

## Small-model optimization tips

- keep prompts under ~250 lines
- avoid long prose; prefer bullet lists
- provide exact IDs and paths
- request explicit commands and explicit stop points
- iterate in short turns: plan first, execute second
