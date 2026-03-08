# Showcase Page: TP73 Promoter -> Luciferase Assay Planning in GENtle

Last updated: 2026-03-08

This page is a narrative, web-ready demonstration of how GENtle can be used to
plan cloning of TP73 promoter fragments into a luciferase reporter context,
with deterministic replay through GUI, CLI, and agent-facing interfaces.

Primary audience: collaborators evaluating GENtle as an external skill/tool
integration target (for example ClawBio/OpenClaw ecosystems).

## Why this showcase

Promoter-reporter assay planning often gets split across many disconnected
tools. This showcase demonstrates one reproducible flow in GENtle:

1. retrieve genomic context from a prepared reference,
2. extract promoter candidates upstream of TP73,
3. import luciferase reporter vector sequence,
4. preview assembly strategy,
5. design primer/qPCR verification sets,
6. keep the complete chain exportable and replayable.

## Biological framing

TP73 has alternative promoter usage and isoform complexity. For planning
purposes, we treat at least two promoter contexts as first-class candidates:

- **TAp73-associated upstream promoter context**
- **DeltaN/DeltaNp73-associated alternative promoter context**

This page focuses on deterministic in-silico planning, not wet-lab claims.

## Canonical assets used

- GUI walkthrough:
  - `/Users/u005069/GitHub/gentle_rs/docs/tutorial/tp73_promoter_luciferase_gui.md`
- Canonical workflow skeleton:
  - `/Users/u005069/GitHub/gentle_rs/docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json`
- Vector accession used in the flow:
  - NCBI AY738222 (Promega luciferase reference context)

## Narrative flow (operator view)

### Step A: Prepare genome and pull TP73 locus

In GENtle:

- Prepare reference genome (`Human GRCh38 Ensembl 116`)
- Retrieve TP73 gene entry
- Extend anchored context 5' to capture promoter space

Deterministic parity:

```bash
cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP73 --occurrence 1 --output-id grch38_tp73 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extend-anchor grch38_tp73 5p 2500 --output-id grch38_tp73_promoter_context --catalog assets/genomes.json --cache-dir data/genomes
```

### Step B: Generate promoter candidates (TA/DeltaN planning set)

Use anchored extraction around TP73-relevant feature boundaries to generate a
candidate panel that can include both canonical and alternative promoter
contexts (TAp73 and DeltaN-associated).

In GUI, this is handled through Engine Ops / anchored extraction parameters. In
CLI, the canonical workflow JSON captures the same payload.

### Step C: Import luciferase vector and preview construct strategy

- Fetch or load AY738222 luciferase vector context.
- Use routine/macro-assisted assembly preview to generate construct candidates.

### Step D: Design validation assays

- Design standard primer pairs for construct verification.
- Design qPCR assays for junction/expression-focused checks.
- Export reports for review and reproducibility.

### Step E: Keep results traceable

- Persist lineage states in project files.
- Export machine-readable reports and SVGs.
- Re-run deterministically through CLI or agent wrappers.

## Interface parity (what this demonstrates externally)

The same planning logic can be run from:

- GUI operations
- `gentle_cli` direct commands
- shared shell (`gentle_cli shell '...'`)
- Python adapter (`integrations/python/gentle_py`)
- ClawBio skill wrapper (`integrations/clawbio/skills/gentle-cloning`)
- MCP tool-calling (`gentle_mcp`)

This is the core interoperability point for external partners: GENtle is
exposed as deterministic tooling, not opaque prompt-only behavior.

## Suggested demo script for partner conversation

1. Open this page and the GUI tutorial side-by-side.
2. Run the canonical workflow skeleton once in CLI.
3. Show the same endpoint through:
   - Python `GentleClient`
   - ClawBio `gentle-cloning` skill wrapper
4. Compare resulting report/provenance artifacts.

## Related docs

- `/Users/u005069/GitHub/gentle_rs/docs/tutorial/tp73_promoter_luciferase_gui.md`
- `/Users/u005069/GitHub/gentle_rs/docs/agent_interfaces_tutorial.md`
- `/Users/u005069/GitHub/gentle_rs/docs/protocol.md`
- `/Users/u005069/GitHub/gentle_rs/docs/release_notes_v0.1.0-internal.2.md`
