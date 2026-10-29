# GENtle CLI Manual

Inspect a sequence project, repeat an analysis, or connect GENtle to a script
or agent. The CLI uses the same engine as the desktop application; you do not
need to start with raw operation JSON.

**Quick navigation:** [Start here](#start-here) |
[Find a command or workflow](#commands) | [Interfaces](#overview) |
[Build from source](#build-and-run) | [Agent integration](agent_interface.md).

## Start Here

With `gentle_cli` installed and on your `PATH`, try these discovery commands:

```sh
gentle_cli --version
gentle_cli help
gentle_cli help primers design
```

They report the build, available commands and the syntax for one command;
they do not modify a project or start a biological analysis. Use
`gentle_cli help COMMAND ...` when you need an argument or option explained.
For machine-readable discovery, use `gentle_cli capabilities` or
`gentle_cli help --format json`.

If the executable is not available yet, see [installation](../INSTALL.md).
From a source checkout, replace `gentle_cli` with
`cargo run --locked --bin gentle_cli --`, running from the repository root.
For example: `cargo run --locked --bin gentle_cli -- help primers design`.

For a worked first example, follow
[Load FASTA, branch, and reverse-complement](tutorial/generated/chapters/02-01_load_branch_reverse_complement_pgex_fasta.md).
It explains input records, derived sequences and the matching GUI/CLI steps.
The [task index](#commands) leads to primer, evidence, design and storage routes.

## Work With A Saved Project

Replace the filename below with an existing saved project or a disposable
working copy. This command inspects its contents:

```sh
gentle_cli --state "my project.gentle.json" state-summary
```

Put `--state PATH` before the command and use the same path for subsequent
steps. Commands that change the project save their results there. Without
`--state`, the default is `.gentle_state.json` in the working directory. A
separate CLI process does not inherit the open GUI project's unsaved state.

For commands shared with the GUI Shell, use `gentle_cli shell 'COMMAND ...'`
in the terminal, but enter only `COMMAND ...` in the GUI Shell. The `op` route
accepts one typed JSON operation; `workflow` runs a composed workflow. Those
are alternatives for automation, not prerequisites for using named commands.

The sections below are the detailed reference. Designs and computational
scores still require review; they are not experimental validation.


Operand metavariable conventions:

- Uppercase words in glossary `usage` strings are placeholders, not literal
  values. Replace `SEQ_ID`, `GENOME_ID`, `ENTRY_ID`, `REPORT_ID`, `PANEL_ID`,
  and similar operands with identifiers already present in the active GENtle
  project, catalog, or command result.
- `ID` flags such as `--id`, `--entry-id`, `--output-id`, and `--as-id` usually
  name the local GENtle record to create or reuse. They are not external
  database accessions unless the command explicitly says so.
- `QUERY` is intentionally command-specific. For fetch/search routes it may be
  a gene symbol, accession, stable database identifier, or text query depending
  on the route. For example, `ensembl-gene fetch QUERY` accepts an Ensembl gene
  lookup query such as a gene symbol or stable Ensembl gene identifier, while
  `--entry-id ID` names the local GENtle metadata entry that will store the
  fetched result.
- `SPECIES` and `SPECIES_DIR` should use the identifier expected by the target
  service or catalog. For Ensembl routes prefer species directory names such as
  `homo_sapiens`; avoid informal aliases such as `HUMAN` unless a command
  explicitly documents that alias.
- `PATH`, `OUTPUT`, and suffix-bearing operands such as `OUTPUT.svg` are local
  filesystem paths. GENtle can open a user-selected file or an exact path, but
  glossary commands do not search the operating system from vague descriptions.
- `START`, `END`, `CHR`, `STRAND`, coordinate ranges, and sequence strings keep
  the semantics documented by the specific command family. Do not infer strand,
  coordinate system, or biological provenance from the placeholder name alone.
- If an operand is not clear from the command summary, usage row, or surrounding
  manual section, ask the user for the missing identifier instead of inventing a
  value. This rule is especially important for AI helpers and local LLMs.


### Commands

Looking for a starting point rather than every option? Choose a task below,
or return to the [first-session guide](#start-here). Use
`gentle_cli help COMMAND ...` for the current syntax.

| I Want To... | Start With |
| --- | --- |
| Load FASTA and keep source/derived sequences separate | [02.01: load, branch and reverse-complement](tutorial/generated/chapters/02-01_load_branch_reverse_complement_pgex_fasta.md) |
| Generate and review PCR primer pairs | `primers design`; [04.02: worked offline primer-pair review](tutorial/generated/chapters/04-02_pcr_selection_batch_primer_pairs_offline.md) |
| Assess primer specificity against prepared references | [04.09: reference-bound specificity](tutorial/04-09_reference_bound_primer_specificity.md); [command details](#standalone-reference-bound-primer-specificity) |
| Inspect regulatory evidence in a genomic locus | [08.16: TSS evidence and its CLI companion](tutorial/08-16_tss_regulatory_view_gui.md) |
| Preview changes to an explicitly synthetic coding insert | `sequence-design plan`; [06.07: experimental synthetic redesign](tutorial/06-07_synthetic_sequence_design.md) |
| Draft a protein-expression request without contacting a provider | [09.03: provider-neutral review draft](tutorial/09-03_protein_expression_handoff_cli.md) |
| Project a gene set into a sample arrangement and rack | [03.10: explicit storage projection](tutorial/03-10_gene_set_storage_projection.md) |
| Connect an outer agent or compare interfaces | [01.01: agent interfaces](tutorial/01-01_agent_interfaces.md); [MCP and capability discovery](agent_interface.md) |

**Reference examples, not a script:** the commands below are independent.
Many require previously loaded sequences, reports or local tools. Replace
example IDs and filenames with your own, and use a disposable project when
trying commands that change state. Repository fixture paths assume the
repository root as the working directory.

## Detailed Reference

Choose a topic below for complete commands, examples and limitations.

<!-- manual-catalog:start -->
## Topics

### Getting started

- [Interfaces and Command Discovery](cli/interfaces.md)
- [Build and Run](cli/build.md)
- [GUI Launcher](cli/gui-launcher.md)
- [CLI Options and Progress](cli/cli-adapter.md)
- [File Formats](cli/file-formats.md)
- [Shared Shell](cli/shell.md)
- [Introspection, Facts and Readiness](cli/discovery.md)
- [Protocol Examples and Tutorials](cli/examples.md)
- [Independent Command Examples](cli/command-examples.md)
- [Operation JSON Examples](cli/operation-examples.md)
- [Operation Limitations](cli/limitations.md)
- [Error Behavior](cli/errors.md)
### Projects and navigation

- [Project and Sequence State](cli/projects.md)
- [Collections](cli/collections.md)
- [Operation History](cli/history.md)
- [Window and Display Intents](cli/display.md)
- [Batch Workflows](cli/batch.md)
- [Macros and Routines](cli/macros-routines.md)
### Sequences and annotation

- [Feature Annotation](cli/features.md)
- [Dotplots, Flexibility and Alignment](cli/dotplot.md)
- [Figure and Layout Exports](cli/rendering.md)
- [Construct Reasoning](cli/construct-reasoning.md)
- [Variant and Reporter Evidence](cli/variants-reporters.md)
- [Synthetic Sequence Design](cli/sequence-design.md)
- [Read-Only SNV Intent](cli/snv-intent.md)
### References genes and transcripts

- [Reference Genomes](cli/genomes.md)
- [Helper Sequences](cli/helpers.md)
- [Host Profiles](cli/hosts.md)
- [Gene Groups and Gene Sets](cli/gene-sets.md)
- [Source-Coherent Transcripts](cli/transcripts.md)
- [Cryptic Splicing](cli/splicing.md)
- [Isoform Architecture and Genomic Regions](cli/isoform-architecture.md)
- [Protein and Structure Evidence](cli/protein-evidence.md)
### Primers cloning and experiment planning

- [Primer Design and Handoff](cli/primers.md)
- [Reference-Bound Primer Specificity](cli/primer-specificity.md)
- [Transcript Assay Panels](cli/transcript-assays.md)
- [Restriction Cloning](cli/restriction-cloning.md)
- [Experiment Planning](cli/planning.md)
- [Candidate Sets](cli/candidates.md)
- [Guide Design](cli/guides.md)
- [Sequencing Confirmation](cli/sequencing-confirmation.md)
### Promoters and regulatory evidence

- [Promoter Design](cli/promoters.md)
- [Promoter Cofactors](cli/promoter-cofactors.md)
- [Accession-Pinned TSS Profiles](cli/tss-profiles.md)
- [TFBS and JASPAR Scans](cli/tfbs.md)
- [Repeats and Regulatory Annotations](cli/repeats-regulation.md)
- [Tracks and Arrays](cli/tracks-arrays.md)
- [Multi-Gene Publication Bundles](cli/publication-bundle.md)
### RNA and sequencing reads

- [RNA Read Panel](cli/rna-read-panel.md)
- [RNA Read Interpretation](cli/rna-reads.md)
- [Primer-Anchored Read Trace](cli/primer-trace.md)
- [Primer-Free Payload Evidence](cli/primer-payload.md)
- [Observed Transcript Structures And Junctions](cli/observed-transcripts.md)
- [Target Rescue Screen](cli/rescue-screen.md)
- [RNA Secondary Structure](cli/rna-structure.md)
- [MicroRNA Target Scans](cli/microrna.md)
### Gels pools and physical racks

- [Measured Gel Images](cli/gel-images.md)
- [Racks, Pools and Ladders](cli/racks-pools.md)
- [Panel Workflows](cli/panels.md)
### Configuration services and agents

- [Services and Resource Catalogs](cli/catalogs.md)
- [Engine Settings](cli/settings.md)
- [Latency Diagnostics](cli/diagnostics.md)
- [Screenshot Capture Policy](cli/screenshot-capture.md)
- [Agent Bridge](cli/agents.md)
- [ClawBio and OpenClaw](cli/clawbio.md)
- [JavaScript Shell](cli/javascript.md)
- [Lua Shell](cli/lua.md)
- [Python Wrapper](cli/python.md)
- [MCP Server](cli/mcp.md)
<!-- manual-catalog:end -->

<!-- manual-legacy:start -->
## Earlier Section Links

These entry points retain older bookmarks. Follow the link to the full topic.

## Standalone Reference-Bound Primer Specificity

[Open this topic](cli/primer-specificity.md#standalone-reference-bound-primer-specificity).

### Inspect Retained Panel Specificity

[Open this topic](cli/primer-specificity.md#inspect-retained-panel-specificity).

## Read-Only SNV Intent

[Open this topic](cli/snv-intent.md#read-only-snv-intent-1).

## Synthetic Sequence Design

[Open this topic](cli/sequence-design.md#synthetic-sequence-design-1).

### Explicit Conflict Search

[Open this topic](cli/sequence-design.md#explicit-conflict-search).

### CDS-Wide GC Walkthrough

[Open this topic](cli/sequence-design.md#cds-wide-gc-walkthrough).

### Sliding-Window GC Walkthrough

[Open this topic](cli/sequence-design.md#sliding-window-gc-walkthrough).

## Source-Coherent Transcripts

[Open this topic](cli/transcripts.md#source-coherent-transcripts-1).

## Promoter Cofactors

[Open this topic](cli/promoter-cofactors.md#promoter-cofactors-1).

## Accession-Pinned TSS Profiles

[Open this topic](cli/tss-profiles.md#accession-pinned-tss-profiles-1).

## Measured Gel Images (`.11` Development)

[Open this topic](cli/gel-images.md#measured-gel-images-11-development).

## Multi-gene publication bundle

[Open this topic](cli/publication-bundle.md#multi-gene-publication-bundle).

## Overview

[Open this topic](cli/interfaces.md#overview).

## ClawBio/OpenClaw skill scaffold

[Open this topic](cli/clawbio.md#clawbioopenclaw-skill-scaffold).

## Python module wrapper (`gentle_py`)

[Open this topic](cli/python.md#python-module-wrapper-gentle_py).

## Build and run

[Open this topic](cli/build.md#build-and-run-1).

## Protocol-first example source

[Open this topic](cli/examples.md#protocol-first-example-source).

## `gentle` (GUI launcher)

[Open this topic](cli/gui-launcher.md#gentle-gui-launcher).

## `gentle_js` (JavaScript shell)

[Open this topic](cli/javascript.md#gentle_js-javascript-shell).

### JavaScript shell functions

[Open this topic](cli/javascript.md#javascript-shell-functions).

### JavaScript example

[Open this topic](cli/javascript.md#javascript-example).

## `gentle_mcp` (MCP stdio server)

[Open this topic](cli/mcp.md#gentle_mcp-mcp-stdio-server).

## `gentle_lua` (Lua shell)

[Open this topic](cli/lua.md#gentle_lua-lua-shell).

### Lua shell functions

[Open this topic](cli/lua.md#lua-shell-functions).

### Lua example

[Open this topic](cli/lua.md#lua-example).

## File format expectations

[Open this topic](cli/file-formats.md#file-format-expectations).

## `gentle_cli` (automation/agent CLI)

[Open this topic](cli/cli-adapter.md#gentle_cli-automationagent-cli).

### RNA-seq target rescue screen

[Open this topic](cli/rescue-screen.md#rna-seq-target-rescue-screen).

### Example operations

[Open this topic](cli/operation-examples.md#example-operations).

### Current limitations in the new operation layer

[Open this topic](cli/limitations.md#current-limitations-in-the-new-operation-layer).

## Error behavior

[Open this topic](cli/errors.md#error-behavior-1).
<!-- manual-legacy:end -->
