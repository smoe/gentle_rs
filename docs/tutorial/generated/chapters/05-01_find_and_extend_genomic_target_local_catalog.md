---
chapter_id: "find_and_extend_genomic_target_local_catalog"
title: "Find and extend the right genomic target (local catalog)"
tier: "core"
example_id: "prepare_extract_extend_localproject_gene"
source_example: "docs/examples/workflows/prepare_extract_extend_localproject_gene.json"
example_test_mode: "always"
executed_during_generation: true
automated_status: "passing"
review_status: "unreviewed"
codex_reviewed_at: null
human_reviewed_at: null
generated_artifact_dir: "docs/tutorial/generated/artifacts/find_and_extend_genomic_target_local_catalog"
---

# Find and extend the right genomic target (local catalog)

Use assets/genomes.json + annotation-driven gene targeting, then extend the anchored sequence directly in the DNA window.

Many investigations start with "which exact genomic interval should I look at" rather than "which cloning operation should I run". This chapter introduces the practical sequence-selection loop: choose a prepared genome from `assets/genomes.json`, narrow candidates with annotation-backed gene filtering, extract the target interval, and then widen context with 5'/3' anchor extension as needed for promoter/splice/regulatory interpretation. The same workflow also clarifies where GenBank/EMBL imports fit: exported GenBank/EMBL records are loadable directly, and GenBank entries carrying NCBI `ACCESSION ... REGION:` metadata can be mapped to genome coordinates and extended by the same anchor controls.

**Prerequisites:** Read [Chapter 1: Load FASTA, branch, and reverse-complement](./02-01_load_branch_reverse_complement_pgex_fasta.md) first.

## Parameters That Matter

- `PrepareGenome.genome_id / catalog_path / cache_dir` (where used: operation 1)
  - Why it matters: These fields define which reference source is prepared and where deterministic cache/index files are written.
  - How to derive it: Start with `assets/genomes.json`; use `LocalProject` for offline walkthroughs, then switch to your organism/build.
- `Retrieve dialog Gene filter (regex) and Top matches` (where used: GUI extraction step)
  - Why it matters: Filtering keeps candidate lists interpretable on large annotations and avoids selecting the wrong locus.
  - How to derive it: Use exact regex (`^GENE$`) when known; broaden (`^GENE_PREFIX`) when exploring families.
  - Omit when: Only omit filtering for tiny local annotations where the candidate list is already small.
- `ExtendGenomeAnchor.side / length_bp / output_id` (where used: operation 3 and DNA-window anchor controls)
  - Why it matters: Defines which biological flank is added and by how much; output IDs keep downstream comparisons clear.
  - How to derive it: Choose `5p` for upstream-context questions and `3p` for downstream-context questions; start with 200-2000 bp and iterate.

## When This Routine Is Useful

- You need to start from a biologically relevant gene locus instead of an arbitrary coordinate slice.
- You imported a GenBank/EMBL record and want to recover surrounding genomic context around that entry.
- You need to add upstream/downstream bases iteratively while preserving deterministic provenance.

## What You Learn

- Understand that `assets/genomes.json` defines selectable genome sources for preparation/retrieval.
- Use annotation-backed gene filtering to find target regions reproducibly.
- Extend anchored sequences in GUI/CLI without losing provenance or sequence lineage.

## Applied Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Deterministic Workflows** (`deterministic_workflows`): Operation chains should produce stable IDs and comparable outputs across repeated runs.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.

## GUI First

CLI snippets use GENtle's default `.gentle_state.json` state unless they say otherwise. Add `--state PATH` or `--project PATH` when you want an explicit sandboxed state file for copied commands.

### Step 1: Open File -> Prepare Reference Genome..., set catalog to assets/genomes.json,...

GUI: Open `File -> Prepare Reference Genome...`, set catalog to `assets/genomes.json`, choose `LocalProject`, and prepare it.

CLI:

```bash
cargo run --bin gentle_cli -- genomes prepare LocalProject --catalog assets/genomes.json --cache-dir cache/localproject
```

> Expected: The offline `LocalProject` reference is prepared into `cache/localproject` without network access.

### Step 2: Open File -> Retrieve Genome Sequence..., use Gene filter (regex) to narrow a...

GUI: Open `File -> Retrieve Genome Sequence...`, use `Gene filter` (regex) to narrow annotation hits, pick one gene match, and extract it.

CLI:

```bash
cargo run --bin gentle_cli -- genomes genes LocalProject --catalog assets/genomes.json --cache-dir cache/localproject --filter '^etp' --limit 20
cargo run --bin gentle_cli -- genomes extract-gene LocalProject etpC --occurrence 1 --output-id local_etpc --catalog assets/genomes.json --cache-dir cache/localproject
```

> Expected: The gene listing narrows annotation candidates, and extraction creates the anchored sequence id `local_etpc`.

### Step 3: In the resulting DNA window, use the Extend 5' / Extend 3' anchor controls (n...

GUI: In the resulting DNA window, use the `Extend 5'` / `Extend 3'` anchor controls (next to `Genome anchor`) to add flanking context.

CLI:

```bash
cargo run --bin gentle_cli -- genomes extend-anchor local_etpc 5p 250 --output-id local_etpc_ext5 --catalog assets/genomes.json --cache-dir cache/localproject
```

> Expected: Anchor extension creates `local_etpc_ext5` with widened interval provenance rather than a disconnected sequence copy.


## Follow-up Commands

```bash
cargo run --bin gentle_cli -- genomes list --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes genes LocalProject --catalog assets/genomes.json --cache-dir data/genomes --filter '^etp' --limit 20
cargo run --bin gentle_cli -- genomes extend-anchor local_etpc 5p 500 --output-id local_etpc_ext5_more --catalog assets/genomes.json --cache-dir data/genomes
```

## Checkpoints

- PrepareGenome succeeds for `LocalProject` without network access.
- ExtractGenomeGene produces `local_etpc` from annotation-backed lookup.
- ExtendGenomeAnchor produces `local_etpc_ext5` with widened genomic interval provenance.

## Tutorial Provenance

- Chapter id: `find_and_extend_genomic_target_local_catalog`
- Tier: `core`
- Example id: `prepare_extract_extend_localproject_gene`
- Tutorial source JSON: `docs/tutorial/sources/05-01_find_and_extend_genomic_target_local_catalog.json`
- Workflow file: `docs/examples/workflows/prepare_extract_extend_localproject_gene.json`
- Generated artifact dir: `docs/tutorial/generated/artifacts/find_and_extend_genomic_target_local_catalog`
- Example test_mode: `always`
- Executed during generation: `yes`
- Automated status: `passing`
- Review status: `unreviewed`
- Codex reviewed at: `not recorded`
- Human reviewed at: `not recorded`
- Inspect the source JSON when you need full option-level detail.

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context below.

- Tutorial title: `Find and extend the right genomic target (local catalog)`
- Tutorial/chapter id: `find_and_extend_genomic_target_local_catalog`
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
