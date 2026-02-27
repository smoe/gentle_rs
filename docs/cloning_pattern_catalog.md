# Cloning Pattern Catalog (File-Hierarchy Baseline)

This document defines the current technical catalog layout for workflow macros.

## Catalog root

- `assets/cloning_patterns_catalog/`
- one JSON file per pattern (`gentle.cloning_pattern_template.v1`)
- directory structure defines category/subcategory and is used directly by the
  GUI `Patterns` menu hierarchy

## Current catalog tree

- `restriction/digest_ligation/digest_ligate_extract_sticky.json`
- `sequence/transform/branch_reverse_complement.json`
- `pcr/site_insertion/pcr_site_insertion_then_digest.json`
- `crispr/guides/candidate_scans/grna_candidate_priority_scan.json`
- `crispr/guides/candidate_scans/grna_anchor_window_scan.json`
- `crispr/guides/oligos/grna_practical_filter_and_oligos.json`

## Supported import surfaces

- GUI:
  - `Patterns -> Import Pattern File...`
  - `Patterns -> Import Pattern Folder...`
  - `Patterns -> Import Full Catalog`
  - direct hierarchical selection in `Patterns` submenu (folder tree mirrors
    file tree)
- CLI/GUI shell:
  - `macros template-import PATH`
  - `PATH` can be:
    - single template file (`gentle.cloning_pattern_template.v1`)
    - pack file (`gentle.cloning_patterns.v1`)
    - directory tree (recursive `*.json` import)

## Legacy compatibility

- `assets/cloning_patterns.json` remains supported as the legacy pack format
  (`gentle.cloning_patterns.v1`).
- New catalog files are intentionally additive and do not break existing
  template-import commands.
