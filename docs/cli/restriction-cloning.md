# Restriction Cloning

<!-- manual-counterparts:start -->
Related GUI reference: [Restriction-Site Inspection](../gui/restriction-site-inspection.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - Restriction-cloning handoff notes (`primers prepare-restriction-cloning`):
      - expects an `Operation` payload whose root variant is
        `PrepareRestrictionCloningPcrHandoff`
      - creates extended forward/reverse primer artifacts plus one predicted
        tailed amplicon artifact and one per-handoff container
      - persists report schema
        `gentle.restriction_cloning_pcr_handoff.v1`
      - shell/direct output includes the saved handoff report in the command
        response `output`
    - Restriction-cloning handoff seed helper
      (`primers seed-restriction-cloning-handoff PRIMER_REPORT_ID VECTOR_SEQ_ID ...`):
      - non-mutating helper returning schema
        `gentle.restriction_cloning_pcr_handoff_seed.v1`
      - resolves one saved primer-pair report + pair rank into a ready-to-run
        `PrepareRestrictionCloningPcrHandoff` operation payload
      - if enzymes are omitted, uses the top recommended single-site or
        directed-pair suggestion for the destination vector
      - if directed-pair enzymes are provided explicitly, validates that their
        order matches the same MCS/unique-cut ordering used by vector
        suggestions
    - Restriction-cloning vector suggestion notes
      (`primers restriction-cloning-vector-suggestions SEQ_ID`):
      - non-mutating helper for command-line/agent reasoning against one
        destination vector
      - returns schema
        `gentle.restriction_cloning_vector_enzyme_suggestions.v1`
      - reports `selected_mcs`, `other_unique`, `missing_mcs`,
        `recommended_single_site`, and `recommended_directed_pairs`
      - mirrors the GUI preference order:
        annotated `mcs_expected_sites` first, then other unique cutters
    - Restriction-cloning saved-report helpers:
      - `primers list-restriction-cloning-handoffs`
      - `primers show-restriction-cloning-handoff REPORT_ID`
      - `primers export-restriction-cloning-handoff REPORT_ID OUTPUT.json`

- Shared shell:
  - Commands:
    - Restriction-site scan helper notes (`features restriction-scan`):
      - non-mutating structured result schema:
        `gentle.restriction_site_scan.v1`
      - accepts either stored `SEQ_ID` or inline ASCII DNA via
        `--sequence-text`
      - `--range` / `--start` / `--end` restrict the scan to one local span on
        the chosen operand
      - when no `--enzyme` is supplied, the shared preferred restriction-enzyme
        list is used; if that list is empty, GENtle falls back to the default
        preferred enzyme set
      - the report includes both local scan coordinates and source-sequence
        coordinates, plus optional cleavage geometry unless
        `--no-cut-geometry` is set

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Transcript Assay Panels](transcript-assays.md) | [Next: Experiment Planning](planning.md)
<!-- manual-page-nav:end -->
