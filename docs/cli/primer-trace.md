# Primer-Anchored Read Trace

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->


For a separate primer carry-through investigation, use the experimental
[primer-trace contract](../protocol.md#primer-anchored-read-trace-phase-01) and
[tutorial 10.04 addendum](../tutorial/10-04_targeted_nanopore_fastq_panel.md#primer-trace-addendum-synthetic-only).
These commands work through `gentle_cli shell` (and the direct command router):

```text
rna-reads primer-panel-validate PANEL_BINDING_JSON_OR_@FILE [--seal]
rna-reads primer-trace REQUEST_JSON_OR_@FILE
rna-reads primer-trace REQUEST_JSON_OR_@FILE --cancel
rna-reads primer-trace-show REPORT [--file-id ID --rows 0..1000] [--binding JSON_OR_@FILE]
```

Validation takes an inline/file **binding**, not a guessed primer name. File
bindings require the expected semantic SHA-256. An explicit `--seal` preview returns a digest
for an inline draft with an empty digest, without writing files or accepting a
mismatched pre-existing hash. It is the shared inner-agent route for computing
panel identity, not permission to invent oligos. Trace requests select exact
completed FASTQ files/checksums, reviewed oligos/policy and a fresh output
directory; optionally bind a saved RNA reference panel by digest. No library
files or kit identities are discovered automatically. Payload FASTA is off by
default. Show is metadata plus an explicitly unverified bounded TSV prefix,
not a revalidation of original reads or a second analysis.

`introspect capabilities --kind operation` exposes the corresponding
`request_template` and `operation_template`. They deliberately contain empty
oligos/digests, unset checksums/cancel paths and `declared_complete=false`;
they are typed examples, not ready-to-run requests. Replace them from the reviewed order,
reference and completed-transfer evidence. The example orientation is not
inferred from a gene or kit. Inspect uses the operation template directly.
Defaults are uncalibrated, and `--seal` cannot repair an incomplete panel.
Cancellation reuses the captured run request and its reviewed marker path;
the example's unset path cannot create a placeholder marker.

The inner agent can explain the report and draft these same shell commands.
It must obtain explicit approval for trace artifact writes and cancellation,
even if a suggested plan labels itself read-only. Generic MCP `op` takes
`TracePrimerReads`, `ValidateRnaPrimerPanel`, `InspectPrimerTrace` or
`CancelPrimerTrace` with the same typed payloads and retains its confirmation
gate. These named routes preserve the MCP project file. No GUI is added and
the existing RNA-panel screen/interpretation remains unchanged. Primer
carry-through, local template concordance and unresolved evidence are distinct;
none establishes target enrichment, an exclusive isoform or a mechanism.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: RNA Read Interpretation](rna-reads.md) | [Next: Primer-Free Payload Evidence](primer-payload.md)
<!-- manual-page-nav:end -->
