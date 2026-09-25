# GENtle `v0.1.0-internal.11` Release Notes

Status: published prerelease on 2026-09-24 at
`ffe5c637d5dddaf2f48daa3fa529d6e9bafd4c7b`, with no attached packages.
Exact-candidate package, scientific and GUI acceptance remain pending.

This release centers on **gene-informed primer-pair studies**: inspect a gene's
transcript models and supporting evidence, design and compare assay panels, and
carry explicit unresolved questions into review rather than hiding them behind
a ranked primer list. GENtle remains experimental alpha software; interfaces
and stored formats may change. A designed or assessed assay is not a
laboratory-validated assay.

These notes describe changes since `.10` at
`84f34a9e479d0dee5d8476aba29c379d370f57e1`, published on 2026-09-18 without
recorded exact-candidate acceptance. Its [historical ledger](release_notes_v0.1.0-internal.10.md#exact-candidate-gate-ledger)
remains pending; later checks do not retroactively certify it.

## Gene-Informed Primer Pairs

- A native [primer-pair study workspace](../gene_assay_study_gui_plan.md) joins
  the PCR Designer and Splicing Expert workflow. It exposes shared-engine
  feasibility, design, coverage and discrimination results rather than adding
  a separate GUI design algorithm.
- Large transcript/assay panels now expose their complete paginated contents
  and coverage scope. Review covered mature-cDNA classes, unresolved class
  distinctions and unassessed records separately; partial coverage is not
  presented as a complete isoform-discrimination panel.
- The inner Agent Assistant can discover actual feature IDs and request
  subject-bound Splicing Expert open/focus/close actions. **Use reviewed result
  in next prompt** explicitly stages local structured query/design results for
  the next conversation turn. It does not send automatically or approve
  execution. Full command/framing text counts toward the UTF-8 draft limit;
  oversized results fail without truncating evidence.
- The authentic public [PATZ1 walkthrough](../tutorial/04-08_gene_assay_study_gui.md)
  combines 13 Ensembl transcripts with four RefSeq comparison records on a
  negative-strand locus. A retained Primer3 2.6.1 development replay selected
  seven pairs covering all 13 Ensembl mature-cDNA classes, with nine class pairs
  still unresolved. The four RefSeq records are comparison context, not
  silently added assay targets. These are fixture-bound development results,
  not a release-candidate or wet-lab verdict.

This primer-PAIR workflow is distinct from single-primer Nanopore transcript
capture. Complete-oligo thermodynamics, reference-bound specificity and order
approval remain explicit checks; this release does not make capture candidates
order-ready or reconstruct missing laboratory inputs.

## Source-Aware Transcript Inspection

Ensembl, RefSeq, shared full-exon-chain and source-only views are available in
the locus inspector, linear DNA map and Splicing Expert Structure tab. Shared
projection preserves original accessions, CDS/phase alternatives, strand and
source provenance, including reverse-complemented negative-strand loci.

Comparison filters and highlights do not modify loaded sequence features, PCR
selection or the design universe. Matching exon chains do not imply identical
CDS annotations. Missing source data stays unknown. Invalid comparisons of
assembly labels with genome-catalog keys have been removed; independent
anchor-to-assembly authority remains a [documented limitation](../architecture.md#source-coherent-locus-transcript-presentation).

## TSS Collections And Reports

- TSS inventory and extraction share compound reverse-strand endpoint handling.
  Gene lookup resolves unambiguous symbol/ID links, groups starts by gene
  identity with source/strand separation, and reports ambiguous links and
  unassigned transcripts without silently guessing ownership.
- Saved TSS collections can be inspected and opened; confirmed registry
  removal does not delete their sequences. The synthetic [TSS GUI tutorial](../tutorial/08-15_tss_collection_gui.md)
  covers preview, approval, materialization, repeated opening, Forget and Undo.
  Native-window sizing, toolbar bounds, worker stack allocation and root
  repaint notification fix specific blockers in this lifecycle, not general
  long-running-command responsiveness.
- TSS exports gain opt-in `vector_pdf` with selectable text and audited
  multipage vector composition. Existing `pdf` remains the raster-compatible
  format; SVG peers retain hover details. Font provenance, page geometry,
  receipt binding and failed-export cleanup stay part of the export contract.
- Annotated GenBank/EMBL/FASTA companions remain available through the existing
  [integrated locus/TSS export](../integrated_locus_tss_profiles.md) path.
  This update repairs Windows directory validation, writable-handle PDF
  synchronization and raw traversal checks at export and both receipt readers.
  It does not claim that new real-data companion files have been regenerated
  or independently accepted.

**Existing TSS studies:** the gene-identity corrections can change `tss_id`
and `output_seq_id` for ID-bearing rows. Preview again under a **new collection
ID** and obtain fresh approval. Historical collections remain intact; old
preview hashes and receipts must not authorize the changed grouping.

## Packages And Portability

The current native packaging recipe stages the same five GENtle entrypoints and tracked
resources on Windows, macOS and Linux, with revision/version/checksum records.
It includes the GUI, CLI, MCP, examples/docs and publication-report tools.
JavaScript/Lua remain optional source interfaces, not packaged binaries.
The `ad0338a7` installer run used Cargo's default optimized release profile,
not fat LTO or a single codegen unit. The subsequent installer repair disables
all native LTO (`lto="off"`), retaining other defaults and one build job, and
retains compiler/resource diagnostics. This changed recipe
requires fresh exact-candidate package checks; older artifacts and acceptance
records retain their original build settings.
The workflow tests actual extracted ZIP/DMG/tarball contents outside the
checkout, including CLI, MCP and tutorial-manifest entrypoints. These checks
do not establish GUI usability, signing or external-tool availability.

The container now compiles genuinely headless tools with `--no-default-features`:
no GENtle GUI, embedded JavaScript/Lua, Xvfb or VNC/noVNC. Required embedded
icon inputs are copied into the builder, but no GUI executable is distributed.
Native desktop packages retain the GUI. Container availability is
separate from native downloads and is not claimed until its own build succeeds.

Windows repairs also cover quoted shell paths, canonical path identity and
byte-bound tutorial/adapter inputs. LF/CRLF checkout tests verify retained
input fingerprints without changing scientific hashing semantics. The shared
shell escape grammar is preserved; missing path-like JSON arguments now
produce file-read diagnostics instead of misleading JSON syntax errors.

See the [release guide](../release.md) for package layout, build-only runs and
the explicit boundary between building and publishing.

## Acceptance Status

The `.11` GitHub prerelease was first published at 07:11:06 UTC on 2026-09-24
while the tag pointed to `ad0338a7e5bd14442ca50722a16eb1677cb82790`.
The owner subsequently moved the tag to
`ffe5c637d5dddaf2f48daa3fa529d6e9bafd4c7b`, whose native release profile sets
`lto="off"`, and republished the prerelease at 14:11:29 UTC. The tag remains at
`ffe5c637`; this history is recorded as fact, not as authorization to move or
re-create it. The current GitHub Release is a published prerelease with no
assets.

| Evidence | Status and scope |
| --- | --- |
| Original installer publication run | [Run 35968255595](https://github.com/smoe/gentle_rs/actions/runs/35968255595) built `ad0338a7`. Windows job `107531770703` passed in 2h12m53s. macOS job `107531770716` ended after the root-library `rustc` received SIGKILL and Cargo returned 101. Ubuntu job `107531770722` ended with runner shutdown/exit 143. Package validation and publication jobs were skipped. The retained Windows artifact is not a package for the final tag revision. |
| Final-tag installer publication run | [Run 36010984447](https://github.com/smoe/gentle_rs/actions/runs/36010984447) built `ffe5c637`. Windows job `107671439337` passed in 1h32m45s. macOS job `107671438533` again failed while compiling the root GENtle library (SIGKILL/Cargo 101); Ubuntu job `107671438797` again ended with exit 143. Package validation and publication jobs were skipped, so no native package was attached. Neither Unix result establishes OOM. |
| Headless container | Tag-push [run 35995175749](https://github.com/smoe/gentle_rs/actions/runs/35995175749) passed at `ffe5c637`, including the no-network RNAPKIN SVG/PNG smoke. Release-event runs [35968255639](https://github.com/smoe/gentle_rs/actions/runs/35968255639) at `ad0338a7` and [36010984241](https://github.com/smoe/gentle_rs/actions/runs/36010984241) at `ffe5c637` both published the headless image to GHCR. Container publication is not native-package or GUI acceptance. |
| Push CI | [Run 35982007068](https://github.com/smoe/gentle_rs/actions/runs/35982007068) passed at `ffe5c637`; [run 35994722502](https://github.com/smoe/gentle_rs/actions/runs/35994722502) passed at `0e3fd09a`. These runs do not replace failed Unix package builds. |
| Glen's scientific/GUI checks | Exact-candidate PATZ1 primer-panel, TP73/DeltaNp73, annotated export, introductory tutorials and inner-agent handoff acceptance pending. Historical synthetic runs do not replace them. |

For historical test counts and their platform/SHA boundaries, see the
[changelog](../CHANGELOG.md) and [reconciled acceptance plan](../internal_11_plan.md).
Local metadata checks alone are not native package or experimental acceptance.

## Known Limits And Next Work

- Authentic whole-reference primer specificity, full-oligo/co-present-pair
  thermodynamics and order approval must be bound to the actual assay and
  laboratory inputs; the public PATZ1 replay is not an approved order sheet.
- Gene-informed primer-study publication still needs its original typed
  dossier/source bundle. Static agent guidance is not proof of a successful
  live-provider conversation or human scientific approval.
- General managed-command migration and bounded GUI latency remain incomplete.
  Glen's new [PATZ1 benchmark evidence](../gui_usability_acceptance_20260920.md)
  separates fast steady CPU painting from slower native resize latency; X11
  snapshot-instrumented timings are not release-binary frame times.
- For `.12`, Glen owns release-like native profiling without snapshot writing,
  startup/first-open measurements, viewport/event-wakeup analysis and a fast
  prebuilt benchmark runner. Runtime optimization follows a measured hotspot,
  not speculative feature caching. Calibrated cofactor prioritization and
  automatic gel-peak detection are not promised in `.11`.
