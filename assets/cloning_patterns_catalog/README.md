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
