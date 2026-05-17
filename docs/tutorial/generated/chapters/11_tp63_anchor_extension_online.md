---
chapter_id: "tp63_anchor_extension_online"
title: "Retrieve TP63 and extend the displayed region by +/-2 kb (online)"
tier: "online"
example_id: "tp63_extend_anchor_online"
source_example: "docs/examples/workflows/tp63_extend_anchor_online.json"
example_test_mode: "online"
executed_during_generation: false
---

# Retrieve TP63 and extend the displayed region by +/-2 kb (online)

Prepare Human GRCh38, inspect TP63 genomic coordinates, extract TP63, and extend the anchored sequence by 2000 bp on each side.

This chapter focuses on day-to-day genome-anchored sequence inspection: identify TP63 in GRCh38, verify coordinates before extraction, and then widen the visible anchored region directly in the DNA sequence window. The key point is that extension remains deterministic and provenance-preserving, whether triggered in the GUI or via shell/CLI operations.

## When This Routine Is Useful

- You want to inspect promoter-proximal and downstream context around TP63 without manually typing coordinates.
- You need to confirm annotated TP63 coordinates before creating an anchored sequence.
- You want a reproducible +/-2 kb extension workflow that can be replayed by GUI and CLI users.

## What You Learn

- Use annotation-backed gene retrieval to inspect coordinates before extracting a sequence.
- Apply anchored-region extension from the DNA sequence viewer without losing genome-anchor provenance.
- Map GUI actions to equivalent deterministic `ExtendGenomeAnchor` operations.

## Concepts

- **Shared Engine Contract** (`shared_engine_contract`): GUI, CLI, shell, and scripting interfaces execute the same operation semantics.
- **Genome Catalog Targeting** (`genome_catalog_targeting`): Prepared genome catalogs, annotation-based gene filters, and anchor extension connect imported entries to genomic context.
- **Sequence Lineage** (`sequence_lineage`): Derived sequences are explicit products linked to upstream inputs and operations.
- **Online Opt-in Execution** (`online_opt_in`): Network-dependent chapters remain explicit opt-in and do not break offline default CI.

## GUI First

1. Open `File -> Prepare Reference Genome...` and prepare `Human GRCh38 Ensembl 116` from `assets/genomes.json`.
2. Open `File -> Retrieve Genome Sequence...`, set gene query/filter to `TP63`, review the displayed TP63 coordinate hit, and extract the first TP63 entry.
3. In the resulting DNA sequence window (`grch38_tp63`), use `Extend 5'` with `2000 bp`, then `Extend 3'` with `2000 bp` to produce the +/-2 kb context sequence.

## Command Equivalent (After GUI)

Run the same routine non-interactively once the GUI flow is clear:

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/tp63_extend_anchor_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/tp63_extend_anchor_online.json'
```

## Parameters That Matter

- `ExtractGenomeGene.gene_query / occurrence` (where used: operation 2 and Retrieve Genome Sequence dialog)
  - Why it matters: Controls which TP63 transcript/gene match is selected when multiple annotation entries exist.
  - How to derive it: Start with `gene_query=TP63` and `occurrence=1`; inspect coordinate listing first and adjust occurrence only if you need a non-primary hit.
- `ExtendGenomeAnchor.side / length_bp` (where used: operations 3 and 4 + DNA window Extend controls)
  - Why it matters: Defines biological flank direction and exact context length added around the anchored region.
  - How to derive it: Use `five_prime,2000` then `three_prime,2000` for symmetric +/-2 kb context expansion.

## Follow-up Commands

```bash
cargo run --bin gentle_cli -- genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP63$" --limit 20
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP63 --occurrence 1 --output-id grch38_tp63 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extend-anchor grch38_tp63 5p 2000 --output-id grch38_tp63_ext5_2kb --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extend-anchor grch38_tp63_ext5_2kb 3p 2000 --output-id grch38_tp63_ext5_ext3_2kb --catalog assets/genomes.json --cache-dir data/genomes
```

## Checkpoints

- TP63 coordinate candidates are visible before extraction in Retrieve Genome Sequence workflows.
- ExtractGenomeGene produces `grch38_tp63` as anchored sequence context.
- Sequential 5'/3' extension by 2000 bp yields `grch38_tp63_ext5_ext3_2kb` with preserved anchor provenance.

## Canonical Source

- Chapter id: `tp63_anchor_extension_online`
- Tier: `online`
- Example id: `tp63_extend_anchor_online`
- Workflow file: `docs/examples/workflows/tp63_extend_anchor_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.
- Inspect this JSON file directly when you need full option-level detail.
