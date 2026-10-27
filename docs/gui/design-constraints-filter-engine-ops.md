## Design Constraints Filter (Engine Ops)

<!-- manual-counterparts:start -->
Related CLI reference: [Candidate Sets](../cli/candidates.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The core operations panel includes `Filter Design` for practical design-constraint
screening.

Inputs:

- `Design inputs`: comma-separated sequence IDs
- `GC min` / `GC max`
- `max homopoly`
- `Reject ambiguous`
- `Avoid U6 TTTT`
- `Forbidden motifs` (comma-separated IUPAC motifs)
- `Unique`
- `prefix`

Defaults in the GUI form:

- `GC min=0.30`, `GC max=0.70`
- `max homopoly=4`
- `Reject ambiguous=true`
- `Avoid U6 TTTT=true`

Execution calls engine operation `FilterByDesignConstraints` and creates filtered
in-silico selection outputs.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Candidate-Set Workflow (Engine Ops + GUI Shell)](candidate-set-workflow-engine-ops--gui-shell.md) | [Next: Primer, qPCR, and Sequence-Specific RT Design](primer-qpcr-and-sequence-specific-rt-design.md)
<!-- manual-page-nav:end -->
