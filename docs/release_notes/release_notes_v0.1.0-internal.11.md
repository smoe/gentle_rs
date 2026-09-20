# GENtle `v0.1.0-internal.11` Release Notes

Status: alpha candidate, prepared 2026-09-20. Final candidate SHA, package
acceptance and release-owner publication approval remain pending.

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

Native packaging now stages the same seven GENtle entrypoints and tracked
resources on Windows, macOS and Linux, with revision/version/checksum records.
The workflow tests actual extracted ZIP/DMG/tarball contents outside the
checkout, including CLI, MCP and tutorial-manifest entrypoints. These checks
do not establish GUI usability, signing or external-tool availability.

The container now compiles genuinely headless tools with `--no-default-features`:
no GENtle GUI, embedded JavaScript/Lua, Xvfb or VNC/noVNC. Required embedded
icon inputs are copied into the builder, but no GUI executable is distributed.
Native desktop packages retain GUI and scripting. Container availability is
separate from native downloads and is not claimed until its own build succeeds.

Windows repairs also cover quoted shell paths, canonical path identity and
byte-bound tutorial/adapter inputs. LF/CRLF checkout tests verify retained
input fingerprints without changing scientific hashing semantics. The shared
shell escape grammar is preserved; missing path-like JSON arguments now
produce file-read diagnostics instead of misleading JSON syntax errors.

See the [release guide](../release.md) for package layout, build-only runs and
the explicit boundary between building and publishing.

## Acceptance Status

The existing `.11` tag points to `e1c7dfb2`, which still declares `.10` in
Cargo metadata. [Container run 35507064138](https://github.com/smoe/gentle_rs/actions/runs/35507064138)
correctly rejected that mismatch before compilation. The version repair is a
later source change: rerunning that old tag cannot pick it up. Tag selection
and any tag correction remain owner-managed; these notes do not authorize them.

| Evidence | Status and scope |
| --- | --- |
| Merged source | `gentle_rs_2_main` and `main` both resolve to `e1c7dfb2` before the version/notes repair; no additional branch changes are missing. |
| Previous push CI | [35500645373](https://github.com/smoe/gentle_rs/actions/runs/35500645373) succeeded at `e1c7dfb2`: macOS full suite/all-feature checks, Linux headless check and release-policy tests. Windows and Linux desktop jobs were skipped. This is not all-platform acceptance of the new candidate. |
| Local version/notes repair | macOS, `e1c7dfb2` plus this update: four release-version tests, 48 release/package/container policy tests, locked offline Cargo check/metadata and formatting passed. The rebuilt GUI reports `.11`; debug linking produced a nonfatal large-unwind-table warning. No full-suite or native-package acceptance is inferred. |
| Corrected candidate | Final SHA and fresh native Windows/macOS/Linux verdicts pending. Record each actual SHA/run, not a nearby revision's pass. |
| Desktop packages / container | Exact-candidate extracted-package checks, headless container build and receipts pending. Failed, cancelled, skipped or absent artifacts are not available downloads. |
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
