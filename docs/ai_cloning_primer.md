# AI Cloning Primer for GENtle Users

Purpose: give local/general LLMs enough biology context to reason about common
GENtle cloning and sequence-design tasks.

Scope: in-silico planning and sequence analysis. This is not a wet-lab protocol
or medical guidance.

## Core objects

- `sequence`: one DNA/RNA string in a project.
- `feature`: annotated region on a sequence (`CDS`, `gene`, `mRNA`, `promoter`,
  custom kinds, `TRACK:*` kinds).
- `container`: a logical collection that can hold multiple sequence candidates.
- `candidate set`: explicit windows/subsequences with metrics used for scoring
  and optimization.
- `local sequence anchor`: an in-sequence coordinate reference (`position`,
  `feature boundary`, `start/end/middle`).
- `genome anchoring`: provenance link to a reference assembly coordinate; this
  is different from local anchors.

## Cloning vocabulary (minimal)

- `vector` / `backbone`: carrier DNA receiving an insert.
- `insert`: sequence fragment to be integrated.
- `CDS`: coding sequence; expected to preserve reading frame.
- `ORF`: open reading frame.
- `promoter`: regulatory sequence driving transcription.
- `RBS`/`Kozak`: translation-initiation context (host-dependent).
- `MCS`: multiple cloning site.
- `restriction site`: enzyme-recognition motif used for digestion/ligation.
- `overhang`: short single-stranded ends after digestion or primer design.
- `ligation`: joining compatible DNA ends.
- `Gibson`: overlap-based assembly method.
- `Golden Gate`: Type IIS-based modular assembly.
- `primer`: short oligo to amplify/validate sequence.

## Coordinate and strand conventions

- GENtle uses zero-based, half-open intervals internally in many workflows.
- `start` is inclusive, `end` is exclusive unless otherwise stated by a command.
- `5'` and `3'` are strand-relative concepts.
- For minus-strand contexts, biological upstream/downstream can be opposite to
  increasing numeric coordinate.
- Always ask/confirm:
  - target sequence ID
  - coordinate convention expected by the user
  - strand assumptions

## Typical in-silico goals

- retrieve genomic region by gene or coordinates
- generate candidate windows between two anchors
- compute scores (`gc_fraction`, distance metrics, weighted objectives)
- optimize with `top-k`, `pareto`, quantile filters, and set operations
- run BLAST checks against reference/helper genomes
- import tracks (BED/BigWig/VCF) and reason near peaks/features

## Common failure modes to avoid

- confusing local sequence anchors with genome provenance anchoring
- assuming promoter/CDS orientation without checking strand
- assuming same coordinate conventions as BED/GTF without verification
- proposing unavailable model/system IDs without discovery first
- generating unsafe auto-execution suggestions instead of explicit ask/confirm

## What strong agent output looks like

- states assumptions explicitly
- asks for missing required inputs before executing
- proposes deterministic GENtle shell/CLI commands
- separates:
  - objective
  - plan
  - commands
  - expected outputs
  - verification checks

## Minimum request fields for cloning-support tasks

When users ask the agent for help, request or infer these fields:

- `objective`: what decision or artifact is needed
- `sequence_ids` and/or `genome_id`
- `constraints`: length, GC range, motif/site constraints, strand constraints
- `anchor definitions`: feature boundaries or absolute positions
- `optimization policy`: weighted score, pareto, top-k, quantiles
- `safety mode`: chat-only, ask-before-run, or allow-auto-exec

## Example compact user brief

"Design 20 bp candidate windows between TP53 gene start and end on sequence
`grch38_tp53`, keep GC 40-80%, maximize distance from CDS boundaries, return top
25 and show the exact `candidates ...` commands. Ask before running anything."
