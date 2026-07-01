# Release Notes: `v0.1.0-internal.5` (draft)

This draft internal release is centered on protein-aware planning and
interpretation, topology- and exon-aware dotplot workflows, the first real
construct-reasoning layer, and a more automation-friendly runtime story for
ClawBio and containerized deployments. Compared with
`v0.1.0-internal.4`, this cut is less about establishing one new confirmation
track and more about connecting previously separate tracks: transcript/protein
inspection, reverse translation, helper/host construct interpretation,
planning-oriented PCR seeding, and agent-friendly execution surfaces.

A second major theme is that more of GENtle's specialist behavior now feels
source-neutral and cross-interface rather than tied to one GUI path:
- Ensembl protein entry intake, UniProt projection, protein comparison,
  reverse translation, and protein-to-DNA handoff now share more explicit
  engine/shell/CLI contracts
- anchored and shared-exon dotplots moved from one-off views toward reusable
  analysis/export surfaces
- helper/host catalogs and construct reasoning now exist as reusable protocol
  records rather than only ad hoc UI state
- container/runtime work now distinguishes headless automation use from GUI use
  more clearly, including first-class Apptainer/Singularity support for the
  ClawBio environment

As with the earlier internal notes, this draft is intentionally written for
internal readers who need both a feature summary and a realistic sense of what
has become stable enough to review.

## Highlights

- Protein-first and transcript-first workflows are now much more substantial:
  - Ensembl protein entry import landed
  - transcript-first protein evidence became inspectable in the GUI and shared
    shell/CLI
  - protein comparison and UniProt projection paths are now source-neutral and
    exportable
  - reverse translation is now treated as a first-class computational artifact
    with persisted reports and codon-choice diagnostics
  - the first explicit protein-to-DNA handoff workflow is in place
- Construct reasoning is now real rather than only planned:
  - reasoning graphs, facts, decision nodes, and “Why/Reasoning” surfaces were
    introduced
  - helper/host construct interpretation moved into reusable shared protocol
    layers
  - planning-oriented helper/host catalogs are now visible in the GUI and
    automation-facing adapters
- Dotplot workflows broadened again:
  - exon-anchored multi-isoform dotplots landed
  - shared exon anchors are accepted in CLI dotplot exports
  - multi-dotplot variants and new anchored tutorial figures/documentation were
    added
- PCR planning is more direct and less hidden:
  - the simple PCR path is now wired and documented
  - current selection can be promoted directly to PCR ROI
  - PCR ROI actions are routed into the dedicated PCR Designer
- ClawBio/container integration is much more practical:
  - Apptainer/Singularity launcher support landed
  - container publication was split into `cli` and `gui` images
  - `singularity run` / `apptainer run` now behave more sensibly in headless
    contexts
- The GUI received another round of stabilization/polish:
  - hosted windows are now safer around screen edges
  - help/configuration/tutorial handling continued to harden
  - linear DNA windows now start with the sequence text panel hidden by default

## Notable Changes by Area

### 1) Protein, Transcript, and Reverse-Translation Workflows

- The first live Ensembl protein path landed and then broadened:
  - transcript-first protein evidence is now available in shared forms
  - protein entries are project-persisted and can be reopened/reused
  - source-neutral adapter layers reduce dependence on one upstream source shape
- UniProt-facing specialist behavior expanded:
  - UniProt features can now be mapped back to exons
  - UniProt specialist SVG export landed
  - direct shell/CLI protein-comparison and projection paths were added
- Reverse translation became a durable computational object rather than a
  transient action:
  - persisted reverse-translation reports now carry codon-choice diagnostics
  - GUI provenance around translation speed / codon trust was tightened
  - reverse translation now has proper shared-shell and direct-CLI parity
- Protein sequences are treated more explicitly as first-class citizens in the
  broader codebase, setting up future work where protein-centric workflows are
  not special cases.

### 2) Construct Reasoning, Helper/Host Interpretation, and Catalogs

- The first construct-reasoning slices landed and then expanded substantially:
  - read-only reasoning layers became available for sequences
  - decision nodes and explicit hard-rule facts were added
  - GUI “Reasoning” / “Why” surfaces now exist on top of the evidence graph
- Helper and host concepts are no longer only documentation-level ideas:
  - helper constructs gained initial semantics
  - host/helper interpretation artifacts moved into protocol/shared layers
  - GUI inspection surfaces now exist for host-profile catalogs
- This matters because planning and evidence are now starting to meet in one
  place: the same project state can carry both sequence edits and reasoning
  records explaining why a route is attractive or blocked.

### 3) Dotplots, Anchors, and Comparative Visualization

- Multi-sequence dotplots advanced into anchored specialist analysis:
  - exon-anchored multi-isoform dotplot mode landed
  - shared exon anchors are now accepted in CLI export paths
  - multiple multi-dotplot variants were added to better support comparison use
    cases
- Documentation/figure work now reflects these modes:
  - anchored dotplot hero assets were added
  - TP73 multi-isoform anchor/overlay figures were added
  - tutorial material now points to more realistic anchored comparison stories
- This makes dotplots more usable for transcript/protein interpretation and not
  only for one generic cDNA-vs-genomic comparison.

### 4) PCR and Primer Workflow Usability

- The simple PCR path is now operational and documented.
- Current selection can be promoted directly to PCR ROI, which makes the
  selection-first workflow much less awkward.
- PCR ROI actions are now routed into the dedicated PCR Designer instead of
  remaining hidden behind less focused entry points.
- A GUI tutorial for simple PCR selection was added, reinforcing the new user
  path.

### 5) Sequencing Confirmation and Specialist Follow-Through

- Sequencing confirmation continued to improve beyond the `internal.4` baseline:
  - GUI inspection surfaces were refined
  - coverage-gap navigation was added to the overview
  - sequencing summaries now include primer guidance
- The track is still focused on confirmation/report review rather than becoming
  a generic chromatogram workbench, but it now feels more like a coherent
  specialist feature than a first landing slice.

### 6) ClawBio, Containers, and Headless Runtime Story

- The ClawBio integration matured from Docker-only assumptions toward more
  realistic HPC/shared-cluster usage:
  - Apptainer/Singularity CLI launcher support landed
  - examples and skill packaging were expanded
  - headless invocation guidance was clarified
- Container publication now has a clearer split:
  - one `cli` image for automation/agent use
  - one `gui` image for browser/VNC-style interactive use
- Entry-point behavior was made friendlier under Singularity/Apptainer so bare
  `run` in headless contexts does not simply fall into the GUI-default path.
- This release therefore tells a stronger deployment story for ClawBio and
  agent-driven environments than `internal.4` did.

### 7) UI/Window Reliability and Ergonomics

- Hosted internal windows were refined further:
  - safer default insets within the root workspace
  - less interference from menu bar / dock edges during move/resize
  - continued hardening of help/configuration/project window behavior
- Sequence-window behavior was adjusted in practical ways:
  - linear sequences now hide the text panel by default
  - PCR/RNA specialist launchers continued to move toward clearer direct entry
    points
- Several GitHub CI/test-drift fixes were landed across the same period,
  keeping the internal baseline moving instead of stalling on accumulated drift.

### 8) Wet-Lab Realism, Gels, and Physicalization

- Gel presentation and realism work continued:
  - first realism baseline
  - further gel-readout refinements
  - lane-side comparison hints
- Container/vial and physical rack surface work also advanced, continuing the
  broader effort to connect abstract planning to bench-facing output.
- These changes are not the headline of the release, but they continue to make
  GENtle feel less like a purely virtual sequence editor.

### 9) Architecture and Internal Engineering Direction

- The workspace/protocol split continued incrementally:
  - lineage and rack state moved further into `gentle-protocol`
  - more source-neutral protocol artifacts now back GUI/CLI/adapter features
- The shell/CLI command surface kept expanding while staying structured:
  - more glossary-backed commands
  - more direct parity for protein, projection, PCR, and BED export paths
- This release continues the same architecture trend as `internal.4`: more
  shared records, fewer one-off adapter-only interpretations.

## Interim Release Readiness

- `v0.1.0-internal.5` is a stronger internal-planning-and-interpretation cut
  than `internal.4`.
- The clearest new narrative thread is now:
  1. fetch/import genomic or protein evidence,
  2. inspect transcript/protein relationships,
  3. reason about constructs/helpers/hosts,
  4. derive or reverse-translate candidate coding DNA,
  5. compare/anchor with dotplots,
  6. seed PCR or related validation workflows,
  7. export/share artifacts through GUI, CLI, MCP, or ClawBio-facing runtimes.
- The codebase is also more adaptable to automation than before:
  the container split and Apptainer support make it easier to use GENtle in
  agent/HPC environments without pretending the GUI is the only runtime.

## Release-Facing Known Limitations

- Construct reasoning is still an early slice:
  - it provides explicit facts/graphs/explanations,
  - but it does not yet rank practical cloning routes by cost, time, reagent
    availability, or local feasibility.
- Protein/transcript workflows are substantially richer now, but they are still
  focused on curated/shared-source interpretation rather than a full
  protein-analysis workbench.
- Reverse translation is more inspectable and provenance-complete, but codon
  optimization and host-specific expression heuristics remain early.
- Dotplot workflows are broader and more powerful, but some anchored/multi-mode
  paths still need more polished tutorial coverage than the baseline now offers.
- Container/runtime ergonomics improved, but native arm64 container publication
  remains intentionally conservative; the current public path is still centered
  on stable amd64 publication and explicit headless use where appropriate.

## Install / Package Notes

- `Cargo.toml` is expected to be bumped to `0.1.0-internal.5` before final
  tagging.
- Container/runtime story for this cycle should be communicated clearly:
  - `ghcr.io/smoe/gentle_rs:cli` for headless/agent/ClawBio use
  - `ghcr.io/smoe/gentle_rs:gui` or `:latest` for interactive GUI container use
- Apptainer/Singularity usage should be documented as a first-class option for
  ClawBio/HPC-style execution, especially for `gentle_cli`.

## Release-Shaped Smoke Check Expectations

- Before final tagging, verify at least:
  - `cargo check -q`
  - `cargo test -q`
  - release/container automation behavior under the split `cli/gui` image model
- Because this cycle touched both GUI defaults and headless/container behavior,
  one short manual smoke pass is also warranted for:
  - opening a fresh linear sequence window
  - opening a circular sequence window
  - one Apptainer/Singularity headless `gentle_cli` invocation
