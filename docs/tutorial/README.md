# GENtle Tutorial Guide

This is the canonical entry page for GENtle tutorials.

Machine-readable catalog source:

- [`docs/tutorial/catalog.json`](./catalog.json)
- generated from:
  - [`docs/tutorial/sources/catalog_meta.json`](./sources/catalog_meta.json)
  - [`docs/tutorial/sources/*.json`](./sources/)

GENtle currently has two tutorial families:

- executable tutorials generated from versioned workflow examples
- hand-written walkthroughs and reference guides for GUI-first and agent-facing use

This page keeps those together so users do not need to know the repository
layout first.

## Start Here

If you are new to GENtle and want the most reproducible path first:

- Start with the executable tutorial hub:
  [`docs/tutorial/generated/README.md`](./generated/README.md)

If you are GUI-first and want one concrete cloning workflow:

- Use the TP73 promoter/luciferase walkthrough:
  [`docs/tutorial/tp73_promoter_luciferase_gui.md`](./tp73_promoter_luciferase_gui.md)

If you want a sequence-analysis example with screenshots:

- Use the TP73 cDNA-vs-genomic dotplot walkthrough:
  [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md)

If you want to use GENtle with agents, MCP, or the command line:

- Use the agent interfaces tutorial:
  [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md)

## Tutorial Catalog

| Tutorial | Type | Status | Best for | Notes |
| --- | --- | --- | --- | --- |
| [`docs/tutorial/generated/README.md`](./generated/README.md) | Executable tutorial collection | `generated+checked` | Reproducible learning paths, CLI parity, CI-backed examples | Generated from `docs/tutorial/manifest.json` and executable workflows; validated by `tutorial-check`. |
| [`docs/tutorial/tp73_promoter_luciferase_gui.md`](./tp73_promoter_luciferase_gui.md) | GUI walkthrough + CLI mapping | `manual/hybrid` | GUI-first cloning planning | Hand-written narrative, but intentionally mapped to engine/CLI operations and linked to executable PCR material. |
| [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md) | GUI walkthrough | `manual` | Screenshot-backed visual analysis tutorial | Good for interface learning; more exposed to UI drift than generated chapters. |
| [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md) | Operational reference tutorial | `manual/reference` | CLI, MCP, in-app agent assistant, external coding agents | Conceptual and operational guide rather than an executable biology walkthrough. |

## Type and Status Labels

Use the labels above as trust/maintenance signals:

- `generated+checked`
  - produced from executable examples and validated by
    `cargo run --bin gentle_examples_docs -- tutorial-check`
  - best choice when reproducibility matters most
- `manual`
  - hand-written walkthrough
  - useful for GUI and screenshot-driven learning
  - higher drift risk when controls or layout change
- `manual/hybrid`
  - hand-written tutorial with explicit links to executable workflows or shared
    engine command paths
  - good bridge between GUI learning and reproducible automation
- `manual/reference`
  - explanatory or operational guide
  - intended for concepts, interfaces, and mental models rather than
    end-to-end biological execution

## Suggested Learning Paths

### Path A: New to GENtle

1. Read [`docs/tutorial/generated/README.md`](./generated/README.md)
2. Run one `core` chapter
3. Open the matching state/artifact outputs

### Path B: GUI-first cloning planning

1. Read [`docs/tutorial/tp73_promoter_luciferase_gui.md`](./tp73_promoter_luciferase_gui.md)
2. Use the executable PCR chapter it references
3. Save the resulting project state for later CLI/agent replay

### Path C: Sequence-analysis and visualization

1. Read [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md)
2. Compare its screenshots to the live GUI state
3. Re-run the corresponding workflow example if you want a deterministic audit
   path

### Path D: Agents and automation

1. Read [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md)
2. Verify capabilities with `gentle_cli help` or MCP `tools/list`
3. Use executable tutorials as the reproducible backing layer

## Contributor Notes

When adding a tutorial, decide and document all three of these explicitly:

1. Tutorial type: walkthrough, executable chapter, or reference guide
2. Status: `generated+checked`, `manual`, `manual/hybrid`, or
   `manual/reference`
3. Primary audience: GUI-first, CLI-first, agent-facing, contributor, or
   biology-specific

Contributors should link new tutorials from this page in the same change so
discovery remains centralized.

Machine-readable discovery should also be updated in the same change:

- [`docs/tutorial/catalog.json`](./catalog.json)
- [`docs/tutorial/sources/catalog_meta.json`](./sources/catalog_meta.json)
- relevant [`docs/tutorial/sources/*.json`](./sources/)

## Maintenance Direction

Current chosen direction:

- keep executable tutorials grounded in `docs/tutorial/manifest.json` and
  workflow examples
- keep discovery metadata grounded in:
  - generated catalog: `docs/tutorial/catalog.json`
  - per-tutorial source units under `docs/tutorial/sources/`
- keep manual/tutorial-reference pages, but surface them through this one
  canonical landing page
- keep the distinction between tutorial, recipe, and reference material
  explicit in page metadata and catalog text

Planned next documentation evolution:

- converge toward one unified tutorial catalog
- eventually move from one large manifest plus scattered manual pages toward
  smaller per-tutorial authoring units with generated indexing
