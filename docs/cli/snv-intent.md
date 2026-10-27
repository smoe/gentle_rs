# Read-Only SNV Intent

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

## Read-Only SNV Intent

The shared `EvaluateGenomeEditIntent` operation validates one proposed base
replacement without changing the project. It does **not** choose an editor or
guide, or establish that an experiment will work. Use the existing generic
`gentle_cli op OPERATION_JSON_OR_@FILE` route, also available in GUI Shell,
scripting and agent execution. No new `editing ...` command is registered.

This complete operation is a synthetic six-base example, not a gene or guide:

```json
{
  "EvaluateGenomeEditIntent": {
    "request": {
      "schema": "gentle.genome_edit_intent_request.v1",
      "target": {"kind": "inline_sequence", "sequence": "ACGTAC"},
      "reference": {"kind": "synthetic", "namespace": "synthetic:snv-tests", "reference_id": "six-base-window"},
      "anchor": {"start_1based": 1001, "end_1based": 1006, "strand": "plus"},
      "expected_sequence_sha256": "sha256:e76797dcfb9f3a4639fdf013e8a72d0bd943a105bfc397387335c6540c496559",
      "design_haplotype": "reference",
      "intent": {"start_0based": 1001, "end_0based_exclusive": 1002, "deleted": "C", "inserted": "T"}
    }
  }
}
```

Pass that complete operation as the `op` argument or explicit `@FILE` input.
The CLI's returned `genome_edit_intent.desired_genomic_forward_sequence` is
`ATGTAC`, and its local changed interval is `[1,2)`. Shared Shell wraps the same
report under `result.genome_edit_intent`. The anchor is 1-based
inclusive; the intent is genomic-forward 0-based half-open. For the same
reference in a minus-oriented local window, supply `GTACGT`, its own expected
hash and `strand: minus`; the desired local sequence is then `GTACAT` at `[4,5)`.
Do not reverse the genomic alleles a second time.

Requests require an exact uppercase-DNA semantic hash and explicit reference
and haplotype labels. Case-only DNA normalization is allowed; ambiguity and
whitespace are not. The [protocol contract](../protocol.md#read-only-snv-intent)
details strict identity, 100000-base limits, optional annotation evidence,
loaded-source provenance checks and refusals. Missing annotations remain
unavailable. A design haplotype label is not a measured cell-line genotype.
Generic MCP `op` uses the same payload with its existing `confirm=true` gate,
without saving or creating a project file for this read-only operation.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Synthetic Sequence Design](sequence-design.md) | [Next: Reference Genomes](genomes.md)
<!-- manual-page-nav:end -->
