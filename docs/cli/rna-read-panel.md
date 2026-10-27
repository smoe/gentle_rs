# RNA Read Panel

<!-- manual-counterparts:start -->
Related GUI reference: [RNA Read Panel](../gui/rna-read-panel.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

For targeted multi-gene cDNA FASTQ, use [tutorial 10.04](../tutorial/10-04_targeted_nanopore_fastq_panel.md):
`rna-reads panel-build --list-sources [GENE_SET_ID]` discovers loaded groups and stored gene-set bindings without writes,
using the same inventory as the native editor. Review blank provenance; do not
infer assembly/release or stable gene IDs from names. To extend an existing
panel, add optional `base_panel_path` to the existing build request, list only
new loci in `sources`, and choose a fresh `output_path`. Original loaded sources
must still match their hashes. Recompose explicitly with all desired sources
and no base to edit/remove retained loci. Old panels/screens remain unchanged;
new panels need new screening and cannot be merged into old-panel summaries.
Reuse an existing resolved gene set with `gene_set: {resolution_id,
resolution_sha256}` from the inventory; the optional ID requests exact loaded
member mappings. Preserve the discovered `expected_source_record_sha256` values.
Unresolved, missing or competing stable IDs are not zero-transcript results.
Set `register_in_project=true` to retain the same panel/set/parent graph snapshot
as the native editor. Omit it for legacy file-only behavior; project undo does
not remove external files. The panel embeds the gene-set snapshot and, on
extension, `parent_panel_sha256` as part of its content identity.

The other shared routes are
`rna-reads panel-demo OUTPUT_DIR`, `rna-reads panel-build REQUEST_JSON_OR_@FILE`,
`rna-reads panel-screen REQUEST_JSON_OR_@FILE`,
`rna-reads panel-show PATH [--file-id ID] [--rows N] [--verify-artifacts]` and
`rna-reads panel-cancel SCREEN_REQUEST_JSON_OR_@FILE`, and
`rna-reads panel-export-summary SUMMARY_REQUEST_JSON_OR_@FILE` through `gentle_cli shell`.
Compressed FASTQ is parsed directly; expected MD5/SHA-256 verifies during that
same pass, not upfront. Retain partial seed co-nomination, not whole-read fit or exclusive isoform
counts. Requests bind completed individual files and explicit sample/barcode/
priming. Fresh external outputs require approval; no input overwrite or resource
download. [Protocol and limits](../protocol.md#multi-locus-fastq-screening).
Summary requests bind explicit `report_paths` and a fresh `output_dir`; export
sample/transcript support and class-bound co-nomination sets without reopening FASTQ.
These are reads, not molecule copies or exclusive isoform counts. Reopen a
saved summary with `rna-reads panel-show PATH`; writes require approval.
Add `--verify-artifacts` only to a summary to verify its three table receipts
within a 256 MB budget; default inspection does not read tables. Named missing
inputs and `unassessed` sample rows never mean zero expression. Legacy v1
class/member associations stay unresolved. Progress reports compressed bytes
with decoder read-ahead, not checksum success. Use fresh request-bound cancel
paths; matching pre-start cancellation works, stale markers fail explicitly.
Before startup, marker errors refuse the run. Mid-run marker errors return a
`control_error` report and warning with verified file units preserved, not a
completed or authorized-cancel analysis. Inspect its diagnostic and missing
inputs via `panel-show`; summaries retain the diagnostic and stay partial.
An ordinary process return/report publication alone is not analysis success.
Pin the reviewed panel with optional `expected_panel_sha256` in the screen
request; a different valid panel at the same path is refused before streaming.
New results distinguish submitted-request hashes from effective requests after
manifest resolution. A retained `gentle.rna_read_panel_result_binding.v1` can be
supplied as `panel-show PATH --binding BINDING_JSON_OR_@FILE` for exact reopening;
do not substitute an arbitrary existing `report.json` after a failed run. This
rechecks metadata bindings, not the original reads or a now-changed manifest.

`rna-reads panel-screen --request REQUEST_JSON_OR_@FILE [--save NEW_FILE]`
only loads/round-trips typed draft metadata. It does not inspect the panel,
read a manifest, open FASTQ or approve execution. `--save` atomically writes a
fresh JSON file and requires agent confirmation; it never overwrites a draft.
This shared-shell flag also serves MCP `agent_execute_plan`, scripts and the
native request editor. The scientific raw operation remains `ScreenRnaReadPanel`.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Multi-Gene Publication Bundles](publication-bundle.md) | [Next: RNA Read Interpretation](rna-reads.md)
<!-- manual-page-nav:end -->
