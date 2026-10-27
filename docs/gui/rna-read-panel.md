## RNA Read Panel

<!-- manual-counterparts:start -->
Related CLI reference: [RNA Read Panel](../cli/rna-read-panel.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Prepare targeted RNA-read panels in **Genome > RNA Read Panel...**, also available
through the command palette and `ui open|focus|close rna-read-panel`. Refresh
loaded loci and stored gene sets, bind the target set to independently annotated gene groups, review exact identities
and provenance, and retain all annotated transcripts unless deliberately restricted.
Use Retrieve gene... for the existing prepared-reference import window; it does
not automatically download or select genes. Unknown provenance remains blank.
Use **Gene sets... > Resolve new set > Explicit members** to reuse the existing
gene-set model. Select that stored resolution and **Bind loaded gene-set members**;
exact stable IDs are required, unresolved or competing sources refuse binding.
Manual source choices remain available to resolve competing loaded views, not
to redefine the gene set by labels. Expanding the targets requires a new resolved
set, then a new panel version. The source fields describe the reference binding.
Build new panel saves the existing shared build request to a fresh panel file.
Open panel retains its original source rows; adding loci extends it after source
hash checks. Recompose all sources explicitly unlocks retained rows for editing
or removal. Neither action overwrites earlier panels or reinterprets old screens.
Workers use history-free project snapshots; other windows remain usable. If the
project changes during a successful save, the file is disclosed as an earlier
snapshot, not promoted as the current panel. Closing the editor is not cancellation
of an approved file write. These controls do not launch FASTQ screening.
Shared-agent screens/summary exports retain their bound external reports
without merging a private journal into the active project. Unrelated edits do
not invalidate published file-only evidence. A cancelled/control-error report
remains inspectable, not successful screening; late host cancellation does not
rewrite an engine-published outcome. Metadata inspection and marker cancellation
have independent bounded capacity.

The **Screen** page accepts explicit completed FASTQ/gzip files, per-file
sample/barcode/priming and checksum bindings, an inspected panel, policy and
fresh output/control paths. Selecting files does not process them. Review the
captured JSON before **Run reviewed screen (writes files)**; later edits require
new review. **Load request...** admits a metadata worker through independent
observation capacity without reading the file on the GUI thread or locking the
live project. It does not consume an execution slot or open referenced reads;
loaded panel pins require fresh inspection. Save typed requests and copy the
matching CLI command without screening. The **Runs** page shows cached shared-service receipts, file ordinal,
current-file reads and compressed decoder read-ahead, not a checksum verdict or
ETA. **Request cancellation** writes the shared request-bound marker; request
failure is distinct from worker cancellation. Closing the window does not cancel
approved work; project edits/replacement do not discard file-only outcomes.
Session-local receipts are not restart recovery. Reopen saved reports explicitly
after restarting. Native acceptance remains with Glen/CI; see the
[active implementation plan](../rna_panel_gui_processing_plan.md).

Agent Assistant **Stop** stops model generation, not a submitted screen.
Its pending-command **Cancel** requests same-process cooperative cancellation,
even for a raw request without `cancel_path`. RNA **Request cancellation** and
`panel-cancel` use the declared marker instead (required by native reviewed
runs); marker-write failures stay visible, with no silent service-cancel fallback.
Neither stopping the model nor replacing a project implicitly cancels an
artifact-only screen. A terminal published outcome cannot be relabelled later.

**Keep panel in project graph** retains a validated snapshot through the same
build operation. The main lineage graph and its SVG export show the target gene
set, panel gene/transcript counts, loaded source DNA relationships and any retained
parent-panel version. These are reference definitions, not read-support results.
Right-click the panel to **Open RNA Read Panel**, or double-click it; gene-set
nodes open the existing inspector. Save/reload preserves the snapshot even if
its external file moves; extension still requires that file and current sources.
Undo removes project registration, not the external panel file. Clearing the
checkbox makes a file-only build; CLI/agents control the same choice with
`register_in_project`.

The native **Inspect** page projects the shared panel, screen or summary record.
New panel previews warn about exact mature sequences shared across loci;
inspection lists their sequence hashes and all locus/member IDs. Check source
provenance: shared RNA can reflect duplicate annotation views or distinct genes.
This diagnostic does not block or merge otherwise independent loci. Missing
historical diagnostics are labelled unrecorded; only selected templates are
compared. Existing overlapping-source admission remains separate.
Published runs reopen with their exact result binding; saved-file inspection
after restart is explicit and does not invent a submitted-service receipt.
The view retains verified/missing file units, control diagnostics, sample
coverage, all class/member associations and non-claims. Nomination prefixes are
bounded and explicitly unverified. New screen inspection also projects optional
index/per-file work counters; missing historical profiles are unrecorded, not
zero. Summary sets show their optional sequence-sharing evidence, not an origin
assignment or explanation of read ambiguity. **Verify summary tables** is separate
bounded background work, not a FASTQ reread. Long count tables use virtualized rows.
Panel and summary inspection also show the engine's optional exon-based source-
overlap assessment, annotated antisense intervals/member IDs and conservative
mapping fallbacks. Legacy absence keeps the historical record-extent rule;
no adapter infers genomic strands or assigns a read's origin.
**Add screen to summary** or **Export > Add screen reports...** selects explicit
inputs. Review the summary request before **Export reviewed summary (writes
files)**; compatibility/deduplication remains engine-owned. Reopen the bound
summary and review its retained source hashes and partial coverage. Nothing is
automatically attached to the project graph.

**Prepare synthetic rehearsal... (writes files)** creates only a fresh synthetic
directory through the existing shared demo operation; it does not load a project
or run reads. The starter includes a resolved toy gene set for native preparation.
Use File > Open Project explicitly and follow the complete no-shell rehearsal
in tutorial 10.04. Capture/runtime/scientific acceptance remains Glen/CI work.

[Tutorial 10.04](../tutorial/10-04_targeted_nanopore_fastq_panel.md)
provides `panel-demo`, `panel-build`, `panel-screen`, `panel-show` and
`panel-cancel` routes for every stage. Export read-support summaries with
`rna-reads panel-export-summary REQUEST_JSON_OR_@FILE`, then reopen them using
`panel-show`. The tables retain class-bound partial co-nomination sets; their counts
are supporting reads, not molecule copies. Export requires write approval.
Use `panel-show SUMMARY_PATH --verify-artifacts` explicitly to verify the three
summary tables, not FASTQ/nomination sidecars. Default inspection is metadata-only.
Missing input/sample ledgers distinguish unassessed from zero. Byte progress
includes decoder read-ahead, not checksum success; no whole-read fit is inferred.
If a marker fails mid-run, the result warns and says `control_error`. Inspect
the diagnostic and verified/missing input ledger with `panel-show`; counts
from the unfinished unit are discarded. The agent must not describe that as
completed analysis or your authorized cancellation. A summary of preserved
units remains partial; fresh control/output paths require an explicit retry.
Use a disposable project before loading
the synthetic starter. Review writes and explicit sample/file bindings; live
model/GUI responsiveness and real-library calibration are not established by
the authored headless tests. For an initial long file, the terminal CLI is an
alternative that does not occupy a GUI command session.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: TSS / Regulatory DNA Display](tss--regulatory-dna-display.md) | [Next: RNA-Binding Protein Evidence](rna-binding-protein-evidence.md)
<!-- manual-page-nav:end -->
