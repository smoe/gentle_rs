# Guide Design

<!-- manual-counterparts:start -->
Related GUI reference: [Patterns And gRNA Setup](../gui/patterns-and-grna-setup.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Guide-design capability status:

- `gentle_cli`: supported as first-class `guides ...` commands and shared-shell `guides ...` commands, backed by shared engine operations (`UpsertGuideSet`, `DeleteGuideSet`, `FilterGuidesPractical`, `GenerateGuideOligos`, `ExportGuideOligos`, `ExportGuideProtocolText`)
- `gentle_js`: supported via `apply_operation` with the same guide-design operations
- `gentle_lua`: supported via `apply_operation` with the same guide-design operations


- Shared shell:
  - Commands:
    - `guides list`
    - `guides show GUIDE_SET_ID [--limit N] [--offset N]`
    - `guides put GUIDE_SET_ID (--json JSON|@FILE|--file PATH)`
    - `guides delete GUIDE_SET_ID`
    - `guides filter GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]`
    - `guides filter-show GUIDE_SET_ID`
    - `guides oligos-generate GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]`
    - `guides oligos-list [--guide-set GUIDE_SET_ID]`
    - `guides oligos-show OLIGO_SET_ID`
    - `guides oligos-export GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]`
    - `guides protocol-export GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]`

Guide-design commands (`gentle_cli guides ...` and `gentle_cli shell 'guides ...'`):

- `guides list`
  - Lists persisted guide sets from guide-design metadata.
- `guides show GUIDE_SET_ID [--limit N] [--offset N]`
  - Pages guide rows from one guide set.
- `guides put GUIDE_SET_ID (--json JSON|@FILE|--file PATH)`
  - Creates/replaces one guide set from JSON array payload.
  - Input payload is `Vec<GuideCandidate>` objects.
- `guides delete GUIDE_SET_ID`
  - Deletes one guide set.
- `guides filter GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]`
  - Applies practical guide constraints and persists a filter report.
  - Optional `--output-set` writes passed guides into a new set.
- `guides filter-show GUIDE_SET_ID`
  - Returns persisted practical-filter report (reasons, warnings, metrics).
- `guides oligos-generate GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]`
  - Generates guide oligos from selected template.
  - Built-ins include `lenti_bsmbi_u6_default` and `plain_forward_reverse`.
- `guides oligos-list [--guide-set GUIDE_SET_ID]`
  - Lists persisted oligo sets (optionally filtered by source guide set).
- `guides oligos-show OLIGO_SET_ID`
  - Shows one oligo set with generated forward/reverse oligos.
- `guides oligos-export GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]`
  - Exports oligo rows for ordering or plate layouts.
- `guides protocol-export GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]`
  - Exports a human-readable wet-lab protocol text for generated oligos.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Candidate Sets](candidates.md) | [Next: Sequencing Confirmation](sequencing-confirmation.md)
<!-- manual-page-nav:end -->
