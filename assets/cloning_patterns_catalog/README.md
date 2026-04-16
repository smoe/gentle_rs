# Cloning Pattern Catalog (Hierarchical)

This directory is the hierarchical source catalog for workflow macro patterns.

Rules:

- One pattern per file (`*.json`).
- File schema: `gentle.cloning_pattern_template.v1`.
- Directory hierarchy defines GUI menu hierarchy (`Patterns` menu).
- `macros template-import PATH` accepts:
  - one template file (`gentle.cloning_pattern_template.v1`),
  - one pack file (`gentle.cloning_patterns.v1`), or
  - a directory tree (recursive `*.json` import).

The legacy pack file `assets/cloning_patterns.json` is still supported for
backward compatibility.

Current family examples in this hierarchical catalog:

- `restriction/.../digest_ligate_extract_sticky.json`
- `gibson/.../gibson_two_fragment_overlap_preview.json`
- `golden_gate/.../golden_gate_single_insert.json`
- `gateway/.../gateway_bp_single_insert.json`
- `topo/.../topo_ta_single_insert.json`
- `ta_gc/.../ta_clone_single_insert.json`
- `infusion/.../infusion_two_fragment_overlap.json`
- `nebuilder_hifi/.../nebuilder_two_fragment_overlap.json`
- `pcr/.../pcr_site_insertion_then_digest.json`
- `crispr/...` guide candidate/oligo workflows
- `sequence/.../branch_reverse_complement.json`
- `reporter/.../allele_paired_promoter_luciferase_reporter.json`
