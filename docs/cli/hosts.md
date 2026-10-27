# Host Profiles

<!-- manual-counterparts:start -->
Related GUI reference: [Protein Expression Handoff](../gui/protein-expression-handoff.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Host convenience commands:

- `hosts list [--catalog PATH] [--filter TEXT]`
  - Lists structured host/strain profile rows used by construct-reasoning host-fit logic.
  - Default catalog is the bundled starter catalog at `assets/host_profiles.json`.
  - Output now also includes the same optional normalized `interpretation`
    record for helper rows that carry structured helper semantics.
- `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Same behavior as `genomes genes`, with helper-catalog default.
- `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
  - Same behavior as `genomes prepare`, with helper-catalog default.
  - `--timeout-secs N`: optional prepare-job timebox.
- `helpers remove-prepared HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes remove-prepared`, with helper-catalog default.
- `helpers extract-promoter HELPER_ID QUERY [--occurrence N] [--transcript-id ID] [--output-id ID] [--upstream-bp N] [--downstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-promoter`, with helper-catalog default.
- `helpers promoter-tfbs-summary HELPER_ID --gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] [--gene ...|--gene-json JSON] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--cache-dir PATH] [--path FILE.json]`
  - Same behavior as `genomes promoter-tfbs-summary`, with helper-catalog default.
- `helpers promoter-cohort-comparison HELPER_ID --cohort-label LABEL --cohort-kind manual|co_regulated|anti_co_regulated --gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] [--gene ...|--gene-json JSON] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--source-seq-id SEQ_ID ...] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--expression-json JSON] [--expression-source-label LABEL] [--cutrun-dataset ID ...] [--cutrun-read-report ID ...] [--catalog PATH] [--cache-dir PATH] [--path FILE.json]`
  - Same behavior as `genomes promoter-cohort-comparison`, with helper-catalog default.
- `helpers promoter-tfbs-svg HELPER_ID --gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] [--gene ...|--gene-json JSON] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--cache-dir PATH] OUTPUT.svg`
  - Same behavior as `genomes promoter-tfbs-svg`, with helper-catalog default.
- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
  - Lists prepared installs and orphaned remnants under the selected cache
    roots with artifact-group byte totals.
  - Default roots are `data/genomes` for references and
    `data/helper_genomes` for helpers.
  - Runtime defaults can be redirected for sibling worktrees via:
    - `GENTLE_REFERENCE_CACHE_DIR`
    - `GENTLE_HELPER_CACHE_DIR`
- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache [--references|--helpers|--both] [--cache-dir PATH ...] [--prepared-id ID ...] [--prepared-path PATH ...] [--include-orphans]`
  - Conservatively deletes derived cache artifacts or prepared installs inside
    the selected roots only.
  - Partial modes keep cached FASTA/annotation sources and manifests intact so
    reindex-from-cached-files remains possible.
  - Catalog JSON, `.gentle_state.json`, MCP/runtime files, backdrop/runtime
    caches, and `target/` are not treated as cache.
- `racks create-from-arrangement ARR_ID [--rack-id ID] [--name TEXT] [--profile small_tube_4x6|plate_6|plate_96|plate_384]`
  - Same shared rack-draft creation path used by GUI `Open Rack`.
- `racks place-arrangement ARR_ID --rack RACK_ID`
  - Same shared rack-append path used by GUI `Place on Existing Rack...`.
- `racks move RACK_ID --from A1 --to B1 [--block]`
  - Same shared rack move/reorder contract used by the rack editor.
- `racks move-samples RACK_ID --from A1 [--from A2 ...] --to B1`
  - Same shared multi-sample rack move contract used by rack-editor
    `Command`/`Ctrl`-selected sample groups.
- `racks move-blocks RACK_ID --arrangement ARR_ID [--arrangement ARR_ID ...] --to B1`
  - Same shared multi-block rack move contract used by rack-editor
    Command/Ctrl-selected arrangement chips.
- `racks show RACK_ID`
  - Emits structured rack placement JSON (`gentle.rack_state.v1`).
  - Rack profile payload includes current `fill_direction` and
    `blocked_coordinates`.
- `racks labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--preset compact_cards|print_a4|wide_cards]`
  - Exports deterministic rack or arrangement-scoped label SVG through the
    same engine path as GUI `Labels SVG`.
- `racks fabrication-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Exports one deterministic top-view fabrication/planning SVG through the
    same engine path as GUI `Fabrication SVG...`.
- `racks isometric-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Exports one deterministic pseudo-3D/isometric rack SVG through the same
    engine path as GUI `Isometric SVG...`.
- `racks hero-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Exports one README-facing top-down rack/plate SVG with labeled
    empty slots, saved-arrangement rings, and an upper-left orientation cut.
- `racks openscad RACK_ID OUTPUT.scad [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Exports one deterministic parameterized OpenSCAD file through the same
    engine path as GUI `OpenSCAD...`.
- `racks set-profile RACK_ID small_tube_4x6|plate_6|plate_96|plate_384`
  - Reflows one saved rack onto another built-in profile.
- `racks set-custom-profile RACK_ID ROWS COLUMNS`
  - Reflows one saved rack onto a custom A1-style geometry while preserving
    arrangement order.
- `helpers remove-catalog-entry HELPER_ID [--catalog PATH] [--output-catalog PATH]`
  - Same behavior as `genomes remove-catalog-entry`, with helper-catalog default.
- `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast`, with helper-catalog default.
- `helpers blast-start HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast-start`, with helper-catalog default.
- `helpers blast-status JOB_ID [--with-report]`
  - Same behavior as `genomes blast-status`, scoped to helper jobs.
- `helpers blast-cancel JOB_ID`
  - Same behavior as `genomes blast-cancel`, scoped to helper jobs.
- `helpers blast-list`
  - Same behavior as `genomes blast-list`, scoped to helper jobs.
- `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-region`, with helper-catalog default.
  - For helper IDs containing `pUC18`/`pUC19`, GENtle auto-attaches a
    canonical MCS `misc_feature` when no MCS annotation is present in source
    annotation and exactly one canonical MCS motif is found.
  - MCS features expose `mcs_expected_sites` with REBASE-normalized enzyme names
    when recognizable from source/fallback annotation text.
- `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-gene`, with helper-catalog default.
  - pUC18/pUC19 helper extractions apply the same automatic MCS fallback
    annotation behavior when applicable (non-unique motif matches are warned and
    skipped).
- `helpers extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Same behavior as `genomes extend-anchor`, with helper-catalog default.
- `helpers verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Same behavior as `genomes verify-anchor`, with helper-catalog default.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Helper Sequences](helpers.md) | [Next: Gene Groups and Gene Sets](gene-sets.md)
<!-- manual-page-nav:end -->
