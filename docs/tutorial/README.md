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

If you want to test the destination-first Gibson specialist end-to-end:

- Use the Gibson specialist testing tutorial:
  [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md)
- fastest setup path:
  [`docs/tutorial/generated/chapters/15_gibson_specialist_testing_baseline.md`](./generated/chapters/15_gibson_specialist_testing_baseline.md)
  - in the app, `File -> Open Tutorial Project...` -> `Gibson Specialist Starter Project (offline)`
    now opens the starter project and the matching Help/Tutorial guide together
- current guardrail:
  multi-insert Gibson execution currently requires a defined destination
  opening; `existing_termini` remains the single-fragment handoff path

If you want to inspect the arrangement that Gibson apply creates and export one
reusable gel lane set from it:

- Use the Gibson arrangements tutorial:
  [`docs/tutorial/gibson_arrangements_gui.md`](./gibson_arrangements_gui.md)
- fastest setup path:
  [`docs/tutorial/generated/chapters/16_gibson_arrangements_baseline.md`](./generated/chapters/16_gibson_arrangements_baseline.md)
  - in the app, `File -> Open Tutorial Project...` -> `Gibson Arrangements Starter Project (offline)`
    now opens the same stable Gibson starter state but lands directly on the
    arrangement guide in Help/Tutorial
- this starts from the same offline Gibson starter project, then continues into
  `Containers`, `Arrangements`, graph inspection, and arrangement-level gel
  export

If you want a sequence-analysis example with screenshots:

- Use the TP73 cDNA-vs-genomic dotplot walkthrough:
  [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md)

If you want to use GENtle with agents, MCP, or the command line:

- Use the agent interfaces tutorial:
  [`docs/agent_interfaces_tutorial.md`](../agent_interfaces_tutorial.md)

## Confidence Map

Treat the tutorial catalog as three confidence bands:

### Recommended now

- The generated executable tutorial hub:
  [`docs/tutorial/generated/README.md`](./generated/README.md)
- The single-insert Gibson specialist testing path, especially when paired
  with the offline starter project and `gibson preview` parity check.
- Tutorials that point back to shared exported artifacts or deterministic
  workflows rather than screenshots alone.

### Works with caveats

- Manual/hybrid tutorials that are intentionally mapped to engine or CLI
  routes, but still depend on live GUI wording and layout.
- Multi-insert Gibson walkthrough work for preview/review, with the current
  apply guardrail on destination opening still in force.

### Exploratory / drift-prone

- Purely manual GUI walkthroughs where screenshots and control labels can drift
  faster than generated runtime chapters.
- Topics whose underlying product areas are still explicitly marked as open
  roadmap gaps rather than stable daily-use paths.

## Tutorial Catalog

| Tutorial | Type | Status | Best for | Notes |
| --- | --- | --- | --- | --- |
| [`docs/tutorial/generated/README.md`](./generated/README.md) | Executable tutorial collection | `generated+checked` | Reproducible learning paths, CLI parity, CI-backed examples | Generated from `docs/tutorial/sources/` and executable workflows through the runtime manifest; validated by `tutorial-check`. |
| [`docs/tutorial/tp73_promoter_luciferase_gui.md`](./tp73_promoter_luciferase_gui.md) | GUI walkthrough + CLI mapping | `manual/hybrid` | GUI-first cloning planning | Hand-written narrative, but intentionally mapped to engine/CLI operations and linked to executable PCR material. |
| [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md) | GUI walkthrough + CLI parity | `manual/hybrid` | Gibson specialist testing, preview/export parity, contributor verification | Hand-written end-to-end test script for `Patterns -> Gibson...` using local inputs plus `gibson preview`; documents the current multi-insert `defined opening` guardrail. |
| [`docs/tutorial/gibson_arrangements_gui.md`](./gibson_arrangements_gui.md) | GUI walkthrough + CLI parity | `manual/hybrid` | Arrangement reuse, Gibson output inspection, gel-lane planning | Hand-written walkthrough for the arrangement that Gibson apply creates automatically, including singleton output containers and arrangement-level gel export. |
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

### Path C: Gibson specialist testing

1. Read [`docs/tutorial/gibson_specialist_testing_gui.md`](./gibson_specialist_testing_gui.md)
2. Run the GUI preview/export steps with local inputs
3. Replay the exported plan through `gibson preview` for parity checking

### Path D: Gibson arrangements and gel export

1. Read [`docs/tutorial/gibson_arrangements_gui.md`](./gibson_arrangements_gui.md)
2. Optionally begin from [`docs/tutorial/generated/chapters/16_gibson_arrangements_baseline.md`](./generated/chapters/16_gibson_arrangements_baseline.md)
3. Apply one single-insert Gibson plan
4. Inspect the resulting singleton output containers and reusable arrangement
5. Export one arrangement-level gel with ladders flanking the sample lanes

### Path E: Sequence-analysis and visualization

1. Read [`docs/tutorial/two_sequence_dotplot_gui.md`](./two_sequence_dotplot_gui.md)
2. Compare its screenshots to the live GUI state
3. Re-run the corresponding workflow example if you want a deterministic audit
   path

### Path F: Agents and automation

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

- keep tutorial authoring grounded in per-tutorial source units under
  `docs/tutorial/sources/`
- generate both:
  - discovery catalog: `docs/tutorial/catalog.json`
  - executable runtime manifest: `docs/tutorial/manifest.json`
- keep workflow execution grounded in `docs/examples/workflows/*.json`
- keep manual/tutorial-reference pages, but surface them through this one
  canonical landing page
- keep the distinction between tutorial, recipe, and reference material
  explicit in page metadata and catalog text

Planned next documentation evolution:

- keep reducing drift between source units and hand-written narrative pages
- explore stronger validation for manual tutorial bodies and page-local metadata
