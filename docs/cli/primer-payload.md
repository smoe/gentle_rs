# Primer-Free Payload Evidence

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->


For experimental external primer-free payload alignments, see the
[Phase 2 command/claim reference](../primer_payload_evidence.md) and the synthetic
[tutorial 10.04 rehearsal](../tutorial/10-04_targeted_nanopore_fastq_panel.md#phase-2-primer-free-payload-evidence-synthetic-only).
These routes share the engine; they do not invoke an aligner or implicitly
reopen FASTQ:

```text
rna-reads primer-payload-prepare REQUEST_JSON_OR_@FILE
rna-reads primer-payload-summarize REQUEST_JSON_OR_@FILE
rna-reads primer-payload-show REPORT [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]
rna-reads primer-payload-cancel CONTROL_JSON_OR_@FILE
```

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Primer-Anchored Read Trace](primer-trace.md) | [Next: Observed Transcript Structures And Junctions](observed-transcripts.md)
<!-- manual-page-nav:end -->
