# Window and Display Intents

<!-- manual-counterparts:start -->
Related GUI reference: [Main window layout](../gui/main-window-layout.md) | [Open windows and focus](../gui/open-windows-and-focus.md) | [DNA Window Panels](../gui/dna-window-panels.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Shared GUI/display control routes:

- `gentle_cli shell 'ui open sequence-window SEQ_ID'`
  - records a non-mutating request to open/focus a DNA sequence viewer for an
    already-loaded sequence record when a GUI host applies the intent.
- `gentle_cli shell 'ui selection sequence-window SEQ_ID --range START..END'`
  - records a non-mutating request to set the DNA-window selection in 0-based,
    end-exclusive coordinates.
- `gentle_cli shell 'ui selection sequence-window SEQ_ID'`
  - records a non-mutating request to inspect the current GUI selection for
    that sequence window.
- `gentle_cli shell 'display show tfbs'` and
  `gentle_cli shell 'display hide restriction-enzymes'`
  - update shared engine display state for sequence-window feature/display
    layers.


- Shared shell:
  - Commands:
    - `facts graph [--evidence SCAN.json ...]`
    - `facts eval FACT_EXPR_JSON_OR_@FILE [--evidence SCAN.json ...]`
    - `ui open TARGET [--genome-id GENOME_ID] [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]`
    - `ui focus TARGET [--genome-id GENOME_ID] [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]`
    - `ui open recent-project ITEM_ID`
    - `ui open tutorial-project CHAPTER_ID`
    - `ui open tutorial-guide TUTORIAL_ID`
    - `ui open configuration [external-applications|agent-systems|microarrays|graphics|language]`
    - `ui focus configuration [external-applications|agent-systems|microarrays|graphics|language]`
    - `ui close configuration`
    - `ui open sequence-window SEQ_ID`
    - `ui open|focus|close splicing-expert SEQ_ID FEATURE_ID`
      - Requires the loaded annotated sequence and zero-based feature id from
        `features query`; opens the same group as
        `inspect-feature-expert SEQ_ID splicing FEATURE_ID` in a GUI host.
        Repeated open/focus reuses the inspection; close retains sequences and
        reports. CLI/MCP only records `applied=false`, not remote GUI success.
    - `ui focus sequence-window SEQ_ID`

- Shared shell:
  - Commands:
    - `ui open|focus|close tss-view --collection COLLECTION_ID` (bounded, persisted TSS collection; headless reports `applied=false`)
    - `ui open|focus tss-view --report REPORT_JSON` (attaches a validated TSS
      profile report to the active annotated TSS viewer; no rescoring or database query)
    - `ui open|focus tss-view --local-score ACCESSION[,ACCESSION...] [--score-kind KIND] [--keep-negative]`
      (locally scored curves for exact full-PFM accessions on the displayed
      window; aliases, `ALL` and consensus fallback are refused. The viewer
      adopts these settings into its own controls and computes on its own
      worker; locally computed lanes stay separate from an attached report and
      from imported evidence. Headless records `applied=false`.)

- Shared shell:
  - Commands:
    - `ui close TARGET`
    - `ui close sequence-window SEQ_ID`
    - `ui prepared-genomes [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]`
    - `ui latest-prepared SPECIES [--helpers] [--catalog PATH] [--cache-dir PATH]`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Operation History](history.md) | [Next: Batch Workflows](batch.md)
<!-- manual-page-nav:end -->
