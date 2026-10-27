## Detailed Reference

<!-- manual-counterparts:start -->
Related CLI reference: [Interfaces and Command Discovery](../cli/interfaces.md) | [Protocol Examples and Tutorials](../cli/examples.md) | [Synthetic Sequence Design](../cli/sequence-design.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The sections below retain the full interface reference, specialist details
and limitations. Use the task index above when you need a starting point.

The existing GUI Shell / Agent Assistant accepts the read-only command
`primers transcript-assay-specificity-matrix PANEL_REPORT_ID`. It inspects
retained reference-selection receipts and binding diagnostics without probing
databases or changing readiness. Policies and historical designs stay separate;
matching standalone oligos do not establish panel association. There is no new
matrix window or native GUI acceptance claim. See
[tutorial 04.09](../tutorial/04-09_reference_bound_primer_specificity.md#inspect-the-panels-retained-evidence-matrix).

Experimental synthetic sequence redesign is available through the existing
GUI Shell / Agent Assistant, not a new editor: `sequence-design plan @request.json`
previews synonymous motif removal, and `sequence-design apply @preview.json
--approve DIGEST` requires explicit approval of exact feasible output. The source
DNA is unchanged; the derived sequence is undoable and can be opened through
`ui open sequence-window OUTPUT_ID`. No input, CDS, protein, protection or digest
may be inferred from a gene label. This is not full DNA Chisel or a `.12` gate.
See [the exact request and non-claims](../cli.md#synthetic-sequence-design) and
[tutorial 06.07](../tutorial/06-07_synthetic_sequence_design.md) for the offline
joint motif/GC preview, separately approved apply and same-session undo/redo.
Plan/S1 inspections preserve the live project, mutation/structural revisions and
undo/redo, including when a report is requested or a preview is refused.
Managed generic Op/Workflow still publishes its journal through the existing
owner/stale-checked commit; it is not an independent observation lane. A concurrent
journal append may make another job stale. Existing explicit report-write approval
and cancellation remain; a written report may survive a cancelled/failed command.
Plan's report writer now requires a fresh filename in an existing directory;
it refuses existing destinations and the current project filename/derived
save-owned sidecar tree. There is no overwrite flag or implicit folder creation.
GUI load/Save As/new-project updates the session binding; a detached result
cannot restore an older binding. Same-directory non-ASCII protected-name
comparisons are conservatively refused: use another existing output folder.
A later explicit Save As can still replace a report; parent-directory races,
crash durability and native acceptance are not established by staging/no-clobber.
No new window, payload resolver or observation lane is introduced.
The same request optionally selects `search_strategy=conflict_directed`; omission
retains enumeration. GUI Shell and Agent Assistant share the engine policy and
need a freshly reviewed approval after a strategy change. This adds no separate
GUI solver, automatic recoding or experimentally validated design claim.
Optional `gc_content={min_basis_points,max_basis_points}` constrains the complete
declared CDS, including frozen start/stop, using inclusive integer GC counts.
Bounds are explicit, not organism presets; flanks are excluded. Empty motifs are
allowed only with explicit GC bounds. GUI Shell/Agent Assistant use the same
preview, integer facts and renewed approval as CLI/MCP; no standalone GC editor.
See the [synthetic GC walkthrough](../cli.md#cds-wide-gc-walkthrough).
Optional `gc_window={window_bp,min_basis_points,max_basis_points}` checks every
complete overlapping CDS window at a fixed one-base step. Flanks/partial edge
windows are excluded; length and bounds must be explicit. The same shared
preview shows every interval/count and apply validates them afresh. See the
[windowed synthetic example](../cli.md#sliding-window-gc-walkthrough); this remains
GUI Shell/Agent Assistant integration, not a new optimization editor.

Source-enriched locus SVGs and detailed TSS figures can show shared physical
exons/CDS, separate structure-reference rows and colour-labelled Ensembl/RefSeq
annotated starts. Hover retains all source identifiers without crowding print.
`=` is exact coordinate agreement, not biological preference. Source inputs are
configured through existing JSON requests, not a new GUI import wizard; the
native annotated-flat-file viewer is unchanged. See
[source-coherent presentation and interactive page export](../transcript_source_presentation.md).

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Next: What To Trust Today](what-to-trust-today.md)
<!-- manual-page-nav:end -->
