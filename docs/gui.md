# GENtle GUI Manual

Open a sequence, explore its annotations, or review the steps of an experiment.
GENtle keeps sequences, derived products and analysis reports together in a
project. You can start with a guided example; no shell commands are required.

**Quick navigation:** [Start here](#start-here) |
[Find your task](#find-your-task) | [What to trust today](#what-to-trust-today) |
[Window controls](#open-windows-and-focus) | [Detailed reference](#detailed-reference) |
[CLI manual](cli.md).

## Start Here

1. **Launch GENtle.** Open the installed desktop app. See
   [installation](../INSTALL.md) if needed, or the
   [source-launch reference](#start-the-gui) if you are working from a checkout.
2. **Try a guided project first.** Choose `File -> Open Tutorial Project...`
   and start with [Load FASTA, branch, and reverse-complement](tutorial/generated/chapters/02-01_load_branch_reverse_complement_pgex_fasta.md).
   The walkthrough explains original inputs, derived sequences and lineage.
3. **Bring your own data when ready.** `File -> Open Sequence...` imports
   sequence files such as FASTA or GenBank into the current project and opens
   their sequence windows. `File -> Open Project...` instead opens a saved
   GENtle project. These are different starting points.
4. **Choose a task.** Focus the sequence you want to work on, then use its
   controls or open the command palette with `Cmd/Ctrl+K`. If an action is
   dimmed, check its details for missing prerequisites.
5. **Review and save.** Review a tool's preview or report before applying its
   proposed changes. Use `File -> Save Project...` to retain your work, using a
   separate filename for a tutorial or working copy. An exported image or report
   is not a replacement for the saved project.

The main project window shows the project and its lineage; sequence windows
show maps and bases. Use the `Windows` menu to return to an open window.
`Edit -> Operation History...` shows available undo/redo steps for this session;
reopening a saved project is not the same as restoring a session's undo stack.

## Find Your Task

| I Want To... | Start With |
| --- | --- |
| Open or create a DNA sequence | `File -> Open Sequence...` or `File -> New Sequence...`; [sequence loading](#loading-and-creating-sequences) |
| Inspect a map or edit annotations | [Map interactions](#map-interactions) and the [Feature Editor walkthrough](tutorial/02-05_feature_editor_gui_cli.md) |
| Design PCR, qPCR or RT primers | Focus the intended DNA, then open PCR Designer; [primer tools](#primer-qpcr-and-sequence-specific-rt-design) |
| Plan a Gibson assembly | `Patterns -> Gibson...`; [03.05: single-insert walkthrough](tutorial/03-05_gibson_specialist_testing_gui.md) |
| Retrieve a genomic region or add reference context | `Genome -> Prepared References...`; [reference-genome tools](#reference-genomes-main-menu) |
| Inspect TSS and regulatory evidence | [08.16: TSS evidence walkthrough](tutorial/08-16_tss_regulatory_view_gui.md) and the [regulatory-view reference](#tss--regulatory-dna-display) |
| Size bands in a gel photograph | `Patterns -> Gel Image Analysis...`; [measured gel images](#measured-gel-images-11-development) |
| Arrange samples or inspect a rack draft | `Open Rack` for an arrangement; [rack view and exports](#rack-view) |
| Review a protein-expression draft | `Services > Protein Expression Handoff...`; [read-only inspector](#protein-expression-handoff) |
| Read the explanation for the current window | Its ordinary `Help` button; [manual navigation](gui/help-manuals.md) |
| Get help with the visible window | `Agent help`; [preview-before-send screenshot help](#screenshot-capture-policy) and [Agent Assistant](#agent-assistant) |

For more worked examples, use `Help -> Tutorials`. For an action's name or
location, try the command palette before browsing the full reference.
The [confidence map](#what-to-trust-today) distinguishes mature and experimental
paths; a computed design or a software test is not biological validation.

## Detailed Reference

Choose a topic below. Every chapter retains the full specialist reference.

<!-- manual-catalog:start -->
## Topics

### Getting started

- [Detailed Reference](gui/detailed-reference.md)
- [What To Trust Today](gui/what-to-trust-today.md)
- [Start the GUI](gui/start-the-gui.md)
- [Documentation automation status](gui/documentation-automation-status.md)
- [Notes and limitations](gui/notes-and-limitations.md)
### Projects and navigation

- [Main window layout](gui/main-window-layout.md)
- [Project Overview And Collections](gui/project-overview-and-collections.md)
- [Physical Arrangements And Project Layout](gui/physical-arrangements-and-project-layout.md)
- [Project And Sequence Window Examples](gui/project-and-sequence-window-examples.md)
- [Status History And Background Jobs](gui/status-history-and-background-jobs.md)
- [Lineage Node Actions](gui/lineage-node-actions.md)
- [Lineage Node Groups](gui/lineage-node-groups.md)
- [Open windows and focus](gui/open-windows-and-focus.md)
- [Help manuals](gui/help-manuals.md)
- [Command Palette and History/Jobs Panels](gui/command-palette-and-historyjobs-panels.md)
- [Engine Ops State Persistence](gui/engine-ops-state-persistence.md)
- [Loading and creating sequences](gui/loading-and-creating-sequences.md)
### Sequences and annotation

- [DNA Window Panels](gui/dna-window-panels.md)
- [Strand Placement](gui/strand-placement.md)
- [Map Modes And Dotplot Workspace](gui/map-modes-and-dotplot-workspace.md)
- [DNA Letters And Zoom](gui/dna-letters-and-zoom.md)
- [Feature Tree And Selection](gui/feature-tree-and-selection.md)
- [Restriction-Site Inspection](gui/restriction-site-inspection.md)
- [Circular Map Labels](gui/circular-map-labels.md)
- [Toolbar buttons](gui/toolbar-buttons.md)
- [Map interactions](gui/map-interactions.md)
- [Linear map conventions](gui/linear-map-conventions.md)
- [Circular map conventions](gui/circular-map-conventions.md)
- [Anchored Region Extraction (Engine Ops)](gui/anchored-region-extraction-engine-ops.md)
### References genes and transcripts

- [Gene Set Inspector And Collection Operations](gui/gene-set-inspector-and-collection-operations.md)
- [Splicing Expert Workspace](gui/splicing-expert-workspace.md)
- [Splicing Support And Transcript Actions](gui/splicing-support-and-transcript-actions.md)
- [Splicing Structure And Evidence Views](gui/splicing-structure-and-evidence-views.md)
- [Locus Figure And Evidence Tracks](gui/locus-figure-and-evidence-tracks.md)
- [Saved Genomic Regions](gui/saved-genomic-regions.md)
- [Conservation Workspace](gui/conservation-workspace.md)
- [Cryptic Splicing Screen](gui/cryptic-splicing-screen.md)
- [Isoform Architecture Panels (Engine Ops)](gui/isoform-architecture-panels-engine-ops.md)
- [Reference Genomes (Main Menu)](gui/reference-genomes-main-menu.md)
### Primers cloning and experiment planning

- [Protein Expression Handoff](gui/protein-expression-handoff.md)
- [Patterns And gRNA Setup](gui/patterns-and-grna-setup.md)
- [Routine Assistant](gui/routine-assistant.md)
- [Gibson Workspace](gui/gibson-workspace.md)
- [Sequencing Confirmation Workspace](gui/sequencing-confirmation-workspace.md)
- [Experiment Planning](gui/experiment-planning.md)
- [Candidate-Set Workflow (Engine Ops + GUI Shell)](gui/candidate-set-workflow-engine-ops--gui-shell.md)
- [Design Constraints Filter (Engine Ops)](gui/design-constraints-filter-engine-ops.md)
- [Primer, qPCR, and Sequence-Specific RT Design](gui/primer-qpcr-and-sequence-specific-rt-design.md)
- [Sequencing Confirmation](gui/sequencing-confirmation.md)
### Promoters and regulatory evidence

- [Promoter Cofactors](gui/promoter-cofactors.md)
- [Accession-Pinned TSS Profiles](gui/accession-pinned-tss-profiles.md)
- [TATA-Box Evidence In The DNA Viewer](gui/tata-box-evidence-in-the-dna-viewer.md)
- [Array Tracks In DNA Viewer](gui/array-tracks-in-dna-viewer.md)
- [TFBS Annotation (Engine Ops)](gui/tfbs-annotation-engine-ops.md)
- [TSS / Regulatory DNA Display](gui/tss--regulatory-dna-display.md)
### RNA and sequencing reads

- [RNA Read Panel](gui/rna-read-panel.md)
- [RNA-Binding Protein Evidence](gui/rna-binding-protein-evidence.md)
- [RNA-Read Mapping Workspace](gui/rna-read-mapping-workspace.md)
- [RNA Structure (ssRNA)](gui/rna-structure-ssrna.md)
- [microRNA Target Scan](gui/microrna-target-scan.md)
### Gels pools and physical racks

- [Measured Gel Images (`.11` Development)](gui/measured-gel-images-11-development.md)
- [Pool Distribution (Engine Ops)](gui/pool-distribution-engine-ops.md)
- [Serial Gel (Ladder View)](gui/serial-gel-ladder-view.md)
- [Rack View](gui/rack-view.md)
### Configuration services and agents

- [Catalog Inspection And Provider Imports](gui/catalog-inspection-and-provider-imports.md)
- [Screenshot Capture Policy](gui/screenshot-capture-policy.md)
- [External Services Window](gui/external-services-window.md)
- [Configuration Window](gui/configuration-window.md)
- [Workspace Latency Diagnostics](gui/workspace-latency-diagnostics.md)
- [GENtle Shell (GUI)](gui/gentle-shell-gui.md)
- [Agent Assistant](gui/agent-assistant.md)
- [About GENtle](gui/about-gentle.md)
- [Engine Settings (Engine Ops)](gui/engine-settings-engine-ops.md)
<!-- manual-catalog:end -->

<!-- manual-legacy:start -->
## Earlier Section Links

These entry points retain older bookmarks. Follow the link to the full topic.

### Promoter Cofactors

[Open this topic](gui/promoter-cofactors.md#promoter-cofactors).

### Accession-Pinned TSS Profiles

[Open this topic](gui/accession-pinned-tss-profiles.md#accession-pinned-tss-profiles).

### Measured Gel Images (`.11` Development)

[Open this topic](gui/measured-gel-images-11-development.md#measured-gel-images-11-development).

## RNA Read Panel

[Open this topic](gui/rna-read-panel.md#rna-read-panel).

## Catalog Inspection And Provider Imports

[Open this topic](gui/catalog-inspection-and-provider-imports.md#catalog-inspection-and-provider-imports).

## What To Trust Today

[Open this topic](gui/what-to-trust-today.md#what-to-trust-today).

### Recommended now

[Open this topic](gui/what-to-trust-today.md#recommended-now).

### Works with caveats

[Open this topic](gui/what-to-trust-today.md#works-with-caveats).

### Exploratory / not yet first choice

[Open this topic](gui/what-to-trust-today.md#exploratory--not-yet-first-choice).

## Start the GUI

[Open this topic](gui/start-the-gui.md#start-the-gui).

### TATA-Box Evidence In The DNA Viewer

[Open this topic](gui/tata-box-evidence-in-the-dna-viewer.md#tata-box-evidence-in-the-dna-viewer).

### Screenshot Capture Policy

[Open this topic](gui/screenshot-capture-policy.md#screenshot-capture-policy).

## Protein Expression Handoff

[Open this topic](gui/protein-expression-handoff.md#protein-expression-handoff).

## External Services Window

[Open this topic](gui/external-services-window.md#external-services-window).

## Configuration Window

[Open this topic](gui/configuration-window.md#configuration-window).

## Main window layout

[Open this topic](gui/main-window-layout.md#main-window-layout).

### DNA Window Panels

[Open this topic](gui/dna-window-panels.md#dna-window-panels).

### Project Overview And Collections

[Open this topic](gui/project-overview-and-collections.md#project-overview-and-collections).

### Gene Set Inspector And Collection Operations

[Open this topic](gui/gene-set-inspector-and-collection-operations.md#gene-set-inspector-and-collection-operations).

### Physical Arrangements And Project Layout

[Open this topic](gui/physical-arrangements-and-project-layout.md#physical-arrangements-and-project-layout).

### Project And Sequence Window Examples

[Open this topic](gui/project-and-sequence-window-examples.md#project-and-sequence-window-examples).

### Strand Placement

[Open this topic](gui/strand-placement.md#strand-placement).

### Map Modes And Dotplot Workspace

[Open this topic](gui/map-modes-and-dotplot-workspace.md#map-modes-and-dotplot-workspace).

### DNA Letters And Zoom

[Open this topic](gui/dna-letters-and-zoom.md#dna-letters-and-zoom).

### Feature Tree And Selection

[Open this topic](gui/feature-tree-and-selection.md#feature-tree-and-selection).

### Restriction-Site Inspection

[Open this topic](gui/restriction-site-inspection.md#restriction-site-inspection).

### Splicing Expert Workspace

[Open this topic](gui/splicing-expert-workspace.md#splicing-expert-workspace).

### Workspace Latency Diagnostics

[Open this topic](gui/workspace-latency-diagnostics.md#workspace-latency-diagnostics).

### Splicing Support And Transcript Actions

[Open this topic](gui/splicing-support-and-transcript-actions.md#splicing-support-and-transcript-actions).

### Splicing Structure And Evidence Views

[Open this topic](gui/splicing-structure-and-evidence-views.md#splicing-structure-and-evidence-views).

### Locus Figure And Evidence Tracks

[Open this topic](gui/locus-figure-and-evidence-tracks.md#locus-figure-and-evidence-tracks).

### Saved Genomic Regions

[Open this topic](gui/saved-genomic-regions.md#saved-genomic-regions).

### Conservation Workspace

[Open this topic](gui/conservation-workspace.md#conservation-workspace).

### Cryptic Splicing Screen

[Open this topic](gui/cryptic-splicing-screen.md#cryptic-splicing-screen).

### RNA-Binding Protein Evidence

[Open this topic](gui/rna-binding-protein-evidence.md#rna-binding-protein-evidence).

### RNA-Read Mapping Workspace

[Open this topic](gui/rna-read-mapping-workspace.md#rna-read-mapping-workspace).

### Circular Map Labels

[Open this topic](gui/circular-map-labels.md#circular-map-labels).

### Status History And Background Jobs

[Open this topic](gui/status-history-and-background-jobs.md#status-history-and-background-jobs).

### Patterns And gRNA Setup

[Open this topic](gui/patterns-and-grna-setup.md#patterns-and-grna-setup).

### Routine Assistant

[Open this topic](gui/routine-assistant.md#routine-assistant).

### Gibson Workspace

[Open this topic](gui/gibson-workspace.md#gibson-workspace).

### Sequencing Confirmation Workspace

[Open this topic](gui/sequencing-confirmation-workspace.md#sequencing-confirmation-workspace).

### Experiment Planning

[Open this topic](gui/experiment-planning.md#experiment-planning).

### Lineage Node Actions

[Open this topic](gui/lineage-node-actions.md#lineage-node-actions).

### Lineage Node Groups

[Open this topic](gui/lineage-node-groups.md#lineage-node-groups).

## Toolbar buttons

[Open this topic](gui/toolbar-buttons.md#toolbar-buttons).

## RNA Structure (ssRNA)

[Open this topic](gui/rna-structure-ssrna.md#rna-structure-ssrna).

## GENtle Shell (GUI)

[Open this topic](gui/gentle-shell-gui.md#gentle-shell-gui).

## Agent Assistant

[Open this topic](gui/agent-assistant.md#agent-assistant).

### Background Commands

[Open this topic](gui/agent-assistant.md#background-commands).

## Candidate-Set Workflow (Engine Ops + GUI Shell)

[Open this topic](gui/candidate-set-workflow-engine-ops--gui-shell.md#candidate-set-workflow-engine-ops--gui-shell).

## Documentation automation status

[Open this topic](gui/documentation-automation-status.md#documentation-automation-status).

## About GENtle

[Open this topic](gui/about-gentle.md#about-gentle).

## Open windows and focus

[Open this topic](gui/open-windows-and-focus.md#open-windows-and-focus).

## Help manuals

[Open this topic](gui/help-manuals.md#help-manuals).

## Array Tracks In DNA Viewer

[Open this topic](gui/array-tracks-in-dna-viewer.md#array-tracks-in-dna-viewer).

## microRNA Target Scan

[Open this topic](gui/microrna-target-scan.md#microrna-target-scan).

## Map interactions

[Open this topic](gui/map-interactions.md#map-interactions).

### Unified scroll/zoom policy

[Open this topic](gui/map-interactions.md#unified-scrollzoom-policy).

### Linear DNA map: zoom and pan (mouse/touchpad)

[Open this topic](gui/map-interactions.md#linear-dna-map-zoom-and-pan-mousetouchpad).

### Lineage graph: zoom and pan (mouse/touchpad)

[Open this topic](gui/map-interactions.md#lineage-graph-zoom-and-pan-mousetouchpad).

## Command Palette and History/Jobs Panels

[Open this topic](gui/command-palette-and-historyjobs-panels.md#command-palette-and-historyjobs-panels).

## Linear map conventions

[Open this topic](gui/linear-map-conventions.md#linear-map-conventions).

## Circular map conventions

[Open this topic](gui/circular-map-conventions.md#circular-map-conventions).

## Pool Distribution (Engine Ops)

[Open this topic](gui/pool-distribution-engine-ops.md#pool-distribution-engine-ops).

## Serial Gel (Ladder View)

[Open this topic](gui/serial-gel-ladder-view.md#serial-gel-ladder-view).

## Rack View

[Open this topic](gui/rack-view.md#rack-view).

## Engine Settings (Engine Ops)

[Open this topic](gui/engine-settings-engine-ops.md#engine-settings-engine-ops).

## Design Constraints Filter (Engine Ops)

[Open this topic](gui/design-constraints-filter-engine-ops.md#design-constraints-filter-engine-ops).

## Primer, qPCR, and Sequence-Specific RT Design

[Open this topic](gui/primer-qpcr-and-sequence-specific-rt-design.md#primer-qpcr-and-sequence-specific-rt-design).

## Sequencing Confirmation

[Open this topic](gui/sequencing-confirmation.md#sequencing-confirmation).

## Anchored Region Extraction (Engine Ops)

[Open this topic](gui/anchored-region-extraction-engine-ops.md#anchored-region-extraction-engine-ops).

## TFBS Annotation (Engine Ops)

[Open this topic](gui/tfbs-annotation-engine-ops.md#tfbs-annotation-engine-ops).

## Isoform Architecture Panels (Engine Ops)

[Open this topic](gui/isoform-architecture-panels-engine-ops.md#isoform-architecture-panels-engine-ops).

## Reference Genomes (Main Menu)

[Open this topic](gui/reference-genomes-main-menu.md#reference-genomes-main-menu).

## Engine Ops State Persistence

[Open this topic](gui/engine-ops-state-persistence.md#engine-ops-state-persistence).

## Loading and creating sequences

[Open this topic](gui/loading-and-creating-sequences.md#loading-and-creating-sequences).

## Notes and limitations

[Open this topic](gui/notes-and-limitations.md#notes-and-limitations).

# TSS / Regulatory DNA Display

[Open this topic](gui/tss--regulatory-dna-display.md#tss--regulatory-dna-display).

### Attach the Quantitative Report

[Open this topic](gui/tss--regulatory-dna-display.md#attach-the-quantitative-report).

### Compute Local TSS Scores

[Open this topic](gui/tss--regulatory-dna-display.md#compute-local-tss-scores).

### Agent and Headless Routes

[Open this topic](gui/tss--regulatory-dna-display.md#agent-and-headless-routes).

### Export the Native TSS View

[Open this topic](gui/tss--regulatory-dna-display.md#export-the-native-tss-view).
<!-- manual-legacy:end -->
