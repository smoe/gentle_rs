# Observed Transcript Structures And Junctions

## Saved Junction Assessment (J1)

For source-relative near offsets, optional motif context and explicit declared
sample recurrence, assess one exact saved S1 report without altering its bytes:

```text
rna-reads splice-junctions-assess REQUEST_JSON_OR_@FILE
rna-reads splice-junctions-show REPORT [--table TABLE.tsv] [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]
rna-reads splice-junctions-cancel REQUEST_JSON_OR_@FILE
```

Supply a strict `gentle.observed_splice_junction_assessment_request.v1` request
with the semantic/raw S1 bindings, near window and all eight budgets explicitly.
Optional FASTA authority and sample hierarchy remain separately declared.
Read the [J1 request and interpretation contract](../observed_splice_junction_assessment.md)
before choosing limits or interpreting support. No whole-RNA/protein/expression
claim follows. Missing evidence is not zero. Generic MCP still requires approval.

## Structural Intake (S1)

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

For saved genomic Phase 2 evidence, [observed transcript structures (S1)](../observed_transcript_structures.md)
accept an explicit stable reference gene and bound annotation sources. No gene,
organism or canonical transcript is assumed; this route does not run an aligner:

```text
rna-reads transcript-structures REQUEST_JSON_OR_@FILE
rna-reads transcript-structures-show REPORT [--table TABLE.tsv] [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]
rna-reads transcript-structures-cancel REQUEST_JSON_OR_@FILE
```

## Prepare A Reference Annotation (S2.1c)

To prepare reference models for one gene, provide a reviewed strict GTF snapshot
and an ordered contig dictionary, each with its actual byte count and SHA-256.
The [versioned preparation request](../protocol.md#s21b-annotation-preparation-direct-engine-only)
uses an exact source-specific stable gene ID, not a symbol or a default gene.
Choose the declared source scope, transcript-row policy and minimum exon gap
explicitly. GENtle does not fetch missing annotations or infer a provider release.

Use clean raw paths in the request, for example `input.gtf` and a fresh
`prepared_gene` output directory with an existing parent. `./input.gtf` and
`out/` are deliberately refused as `observed_annotation_raw_path`; spelling is
checked before normalization. The strict GTF subset rejects empty attribute
values, including `transcript_id "";` on a gene row, before publishing a receipt.
Gene/transcript version attributes must also be non-whitespace identifiers of
at most 512 bytes. Receipts distinguish absent, declared and inconsistent
versions; GENtle does not infer them from transcript names.

Budget admission is not a throughput promise. With the maximum 128-MiB
retained-memory budget, conservative JSON accounting admits at most 2,096,128
dictionary bytes even if the raw dictionary cap is 8 MiB. Distinct GTF row
fingerprints alone limit the input to at most 246,723 distinct data rows;
shared dictionary/model/scratch storage lowers that bound in practice. A
1,000,000-line admission cap does not promise that many rows will fit. See the
[accounting contract](../protocol.md#s21b-annotation-preparation-direct-engine-only)
before choosing a supplied excerpt and explicit budgets.

```text
rna-reads transcript-annotation-prepare REQUEST_JSON_OR_@FILE
rna-reads transcript-annotation-cancel REQUEST_JSON_OR_@FILE
```

Both routes require write approval, even if preparation will be refused or
cancelled. They preserve project content and do not autosave the project.
Shared-shell CLI/MCP/scripting callers use the same request and engine boundary;
generic MCP `op` still requires `confirm:true`. Request files with spaces or
apostrophes must be quoted through the shared shell, as for other RNA commands.

Preparation publishes a fresh `report.json` and, only for a completed gene model,
`artifacts/annotation.json`. A fully inspected but absent/conflicting gene yields
receipt-only evidence, not an empty annotation or a measured biological zero.
Cancellation/control failure before full analysis has no gene counts/outcome;
a late interruption retains complete evidence with the interrupted status.
MCP and managed jobs return a compact result binding with separate dispatch
`outcome` and `gene_outcome`, plus the saved report's raw file binding. An
`unavailable` dispatch result describes a source refusal, not failed file writing.

Cancel writes only a marker for the **identical** request, including its output,
policies and optional `cancel_path`. It does not inspect source bytes, remove
outputs or confirm worker termination. Current paths/aliases are still checked;
source files may be gone if their parent paths remain resolvable. A matching
existing marker is idempotent; a stale/malformed marker is never overwritten.
Managed standalone preparation retains its artifact-only result without a live
project/journal merge. Workflows retain their conservative project-bound host
policy. No SAM ingress, new alignment, transcript interpretation or GUI is
introduced by reference preparation; SAM interpretation is the separate S2.3
route below.

## Interpret Supplied SAM (S2.3)

Use a SAM text export you prepared externally, an exact ordered reference
dictionary, and a manifest assigning every alignment row to a declared query
unit. GENtle does not open BAM or run a converter/aligner. Supply exact file
hashes/byte counts, a source-bound gene ID, annotation bindings and all budgets;
the [SAM request contract](../protocol.md#s23-bound-sam-interpretation) explains
the required fields and nonclaims.

```text
rna-reads transcript-structures-from-sam REQUEST_JSON_OR_@FILE
rna-reads transcript-structures-from-sam-show REPORT [--table TABLE.tsv] [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]
rna-reads transcript-structures-from-sam-cancel REQUEST_JSON_OR_@FILE
```

For example, with your **reviewed and fully bound** request file:

```sh
gentle_cli rna-reads transcript-structures-from-sam '@my SAM request.json'
gentle_cli rna-reads transcript-structures-from-sam-show 'my evidence/report.json' --table junction_catalog.tsv --rows 20 --verify-artifacts
```

These file-only routes do not autosave or repair the active project. Writes
require approval in agent/managed contexts; generic MCP `op` also requires
`confirm:true`. Inspection is read-only, works when original sources are offline,
and defaults to metadata only. Use the returned exact result binding to reject
a different report. Table names are `query_structures.tsv`,
`alignment_observations.tsv`, `structure_catalog.tsv`,
`reference_compatibility.tsv`, `junction_catalog.tsv`, `unavailable.tsv`.

`version_policy` defaults to `strict16`. Choose `legacy10_to16` explicitly for
older version labels, retaining the same strict supported-feature subset.
Do not re-header inputs silently. Use exact S2.1 receipt selectors, preserving
version spellings; a selector mismatch is an error, not a missing gene.
Prepared sources retain refusals; normalized sources require real nonempty
geometry. Search/view completeness and mask qualification are separate inputs.

Read counts as **declared query units**, not molecules or independent samples.
All admitted observations are inspectable, but nonoptimal observations do not
contribute catalog support. Mapping certainty and flank qualification remain
separate. Successful publication is not proof of transcript identity, novel RNA,
exhaustive mapping, verified reference bytes or protein consequences. Cancelled
or control-error analysis has absent counts, not zeros. A completed run can
still disclose missing/refused annotation; inspect those source outcomes.

Cancel requires the identical request and its explicit `cancel_path`, writes
only that request's marker and does not confirm termination. Late cancellation
cannot invalidate an already published completed report. Existing S1 and
Phase 2 report formats remain unchanged.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Primer-Free Payload Evidence](primer-payload.md) | [Next: Target Rescue Screen](rescue-screen.md)
<!-- manual-page-nav:end -->
