# AI Task Playbooks for GENtle

Last updated: 2026-06-21

This document is a compact operating guide for AI assistants that suggest GENtle
workflows. It is intentionally practical: the assistant should use these
playbooks to ask for missing inputs, suggest valid GENtle commands, and stop at
the point where a human must confirm data retrieval or mutation.

## Ground Rules

- Return GENtle commands, not operating-system shell commands.
- Prefer reviewed execution: use `execution="ask"` for runnable suggestions.
- Do not invent file discovery commands. If a local path is needed and unknown,
  ask the user to provide or pick it.
- Do not use Ollama REPL commands such as `/set`, `/show`, `/load`, `/save`,
  `/clear`, or bare `/path/to/file` attachments.
- Use `/open file PATH` or `/import file PATH` only when the user supplied an
  exact sequence-file path.
- Ask before public database or network retrieval.
- In an empty project, do not refer to existing `seq_id` values.

Good first controls:

```text
help
/help
/list
state-summary
capabilities
```

`/list` means GENtle project state and loaded sequence records. It is not a
filesystem directory listing.

## Answer Shape

For planning replies, use this structure:

1. `Assumptions`
2. `Missing inputs`
3. `Plan`
4. `Suggested commands`
5. `Validation checks`
6. `Stop point requiring confirmation`

For GENtle Agent Assistant JSON replies, use `gentle.agent_response.v1` and put
commands in `suggested_commands[]`.

## Playbook 1: Orient An Empty Project

Use this when the user has just opened GENtle or asks what the assistant can do.

Ask or suggest:

```text
state-summary
capabilities
/list
help
```

Expected result:

- The assistant learns whether any sequence/project state exists.
- The user sees valid GENtle controls rather than generic shell habits.
- No network or file operation happens.

## Playbook 2: Load A Local Sequence File

Use this only when the user supplied a concrete path.

Required input:

- exact local file path
- optional desired sequence ID

Suggested commands:

```text
/open file PATH
/open file PATH --id SEQ_ID
/import file PATH
/import file PATH --id SEQ_ID
```

Stop and ask when:

- the path is missing,
- the user pasted an Ollama-style bare `/path/to/file`,
- the file type or intended sequence ID is unclear.

Do not suggest filesystem search commands.

## Playbook 3: Paste A Short Sequence

Use this when the user provides the sequence text directly.

Required input:

- DNA sequence text
- optional sequence ID

Suggested command:

```text
/paste sequence --sequence-text ACGTACGT --id demo_seq
```

Validation checks:

- sequence is not empty,
- sequence alphabet matches the intended molecule,
- chosen ID does not hide what the sequence represents.

## Playbook 4: Retrieve A Public Gene Or Sequence

Use this when the user asks for a gene, accession, protein, or public database
record. Always require confirmation before network retrieval.

Required input:

- source or database if known,
- organism/species,
- gene symbol, accession, protein query, or coordinates,
- desired output ID.

Useful commands:

```text
/fetch ensembl FUS --species homo_sapiens --id fus_live
/fetch genbank ACCESSION --id record_id
/fetch uniprot QUERY --id protein_id
/fetch ensembl-region homo_sapiens CHR START END --id region_id
```

For human genes, prefer `homo_sapiens` over informal labels such as `HUMAN`.

Validation checks:

- confirm species,
- confirm whether the user wants genomic DNA, transcript/CDS, or protein,
- confirm whether isoform annotations are required,
- avoid claiming a `seq_id` exists until the import succeeds.

## Playbook 5: Use Prepared Genome Resources

Use this when a local reference genome catalog/cache is available or the user
wants reproducible extraction from prepared resources.

Required input:

- `genome_id`,
- optional catalog/cache paths,
- gene symbol or coordinates.

Command families:

```text
genomes status GENOME_ID --catalog PATH --cache-dir PATH
genomes prepare GENOME_ID --catalog PATH --cache-dir PATH
genomes genes GENOME_ID --filter "^FUS$" --limit 20 --catalog PATH --cache-dir PATH
genomes extract-gene GENOME_ID FUS --occurrence 1 --output-id fus_grch38 --catalog PATH --cache-dir PATH
genomes extract-region GENOME_ID CHR START END --output-id region_id --catalog PATH --cache-dir PATH
```

Stop and ask when:

- multiple gene hits are possible,
- the requested assembly is not specified,
- preparation would download or create large local resources.

## Playbook 6: Inspect Or Annotate A Loaded Sequence

Use this only after a sequence ID exists.

Required input:

- `seq_id`,
- analysis goal.

Command examples:

```text
features restriction-scan SEQ_ID --enzyme EcoRI
features restriction-scan SEQ_ID --enzyme NotI
```

Validation checks:

- verify the sequence ID exists with `state-summary` or `/list`,
- avoid changing sequence content unless the user explicitly asks.

## Playbook 7: Candidate Windows And Scoring

Use this for design tasks that enumerate and rank local candidate regions.

Required input:

- sequence ID,
- anchor definitions,
- candidate length/step,
- scoring or filtering objective.

Command pattern:

```text
candidates generate-between-anchors SET SEQ_ID --length 20 --anchor-a-json '...' --anchor-b-json '...' --step 1 --limit 10000
candidates show SET --limit 50 --offset 0
candidates score SET gc_balance '100*(gc_fraction-at_fraction)'
candidates top-k SET SET_top --metric gc_balance --k 25 --direction max
```

Stop and ask when:

- anchors are ambiguous,
- the scoring objective is underspecified,
- the command would create many derived records.

## Playbook 8: BLAST Specificity Check

Use this when the user asks whether a candidate is unique or specific.

Required input:

- target catalog type: `genomes` or `helpers`,
- target ID,
- query sequence,
- max hits/task if known.

Command examples:

```text
genomes blast GENOME_ID QUERY_SEQUENCE --max-hits 20 --task blastn-short
helpers blast HELPER_ID QUERY_SEQUENCE --max-hits 20 --task blastn-short
```

Validation checks:

- confirm the target resource is prepared,
- report top-hit identity and alignment summary,
- avoid overinterpreting specificity without the selected database context.

## Playbook 9: Import Evidence Tracks

Use this for BED, BigWig, or VCF evidence mapped onto an anchored sequence.

Required input:

- loaded/anchored sequence ID,
- exact local track path,
- expected assembly/coordinate compatibility.

Command examples:

```text
tracks import-bed SEQ_ID PATH.bed
tracks import-bigwig SEQ_ID PATH.bw
tracks import-vcf SEQ_ID PATH.vcf
```

Stop and ask when:

- the sequence assembly is unknown,
- the track path is missing,
- coordinate compatibility is not established.

## Playbook 10: Ask An Agent But Keep Execution Deterministic

Use this when a configured agent should propose commands but GENtle should keep
execution under review.

Command pattern:

```text
gentle_cli agents ask SYSTEM_ID --prompt "OBJECTIVE" --no-state-summary
```

Use state summary only when project context matters and the model can handle the
extra tokens. For small local models, prefer short prompts and explicit command
examples.

Expected result:

- assistant message,
- optional questions,
- reviewed `suggested_commands[]`,
- deterministic execution through GENtle after user confirmation.

## Tiny Local-Model Primer

For fast local models such as Gemma through Ollama, prepend a small control
card rather than a long manual:

```text
You are inside GENtle, not the Ollama REPL.
Return strict gentle.agent_response.v1 JSON only.
Use GENtle controls: help, /help, /list, state-summary, capabilities.
/list means loaded GENtle records, not filesystem files.
Do not use /set, /show, /load, /save, /clear, or bare /path/to/file.
Use /open file PATH only when the user supplied an exact sequence-file path.
Ask before public database retrieval.
Mark runnable suggestions execution="ask".
```

This primer should be generated from GENtle's real command registry in the
future; until then, keep it short, explicit, and consistent with the parser.
