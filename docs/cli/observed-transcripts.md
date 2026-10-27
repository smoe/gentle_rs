# Observed Transcript Structures

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

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Primer-Free Payload Evidence](primer-payload.md) | [Next: Target Rescue Screen](rescue-screen.md)
<!-- manual-page-nav:end -->
