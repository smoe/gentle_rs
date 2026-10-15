# Cloning Pattern Catalog

Source review: 2026-10-07. The maintained inventory is the JSON directory,
not the original six-pattern starter list.

## Catalog Root And Inventory

`assets/cloning_patterns_catalog/` contains one
`gentle.cloning_pattern_template.v1` JSON file per template. The directory
hierarchy supplies the GUI `Patterns` menu categories. Current inventory:

| Directory family | Templates |
|---|---:|
| crispr | 3 |
| gateway | 3 |
| gibson | 1 |
| golden_gate | 3 |
| infusion | 2 |
| nebuilder_hifi | 2 |
| pcr | 1 |
| reporter | 1 |
| restriction | 1 |
| sequence | 1 |
| ta_gc | 2 |
| topo | 3 |

These are 23 templates, not 23 independently validated laboratory protocols.
The 22-entry `assets/cloning_routines.json` catalog references the cloning/
sequence/guide families; the promoter-luciferase reporter template is separate.
See [routine status](cloning_routine_catalog_plan.md).

From a source checkout, inspect the actual files with
`rg --files assets/cloning_patterns_catalog -g '*.json'`.
Use `routines list` / `routines explain ROUTINE_ID` for declared routine
semantics, then `macros template-list` for templates imported into the project.

## Supported Import Surfaces

- GUI: `Patterns -> Import Pattern File...`, `Import Pattern Folder...`,
  `Import Full Catalog`, or a template in the hierarchical menu.
- CLI/GUI shell: `macros template-import PATH`. The path may be one template,
  a `gentle.cloning_patterns.v1` pack, or a recursively imported directory.

`assets/cloning_patterns.json` remains a supported legacy pack, not the
authoritative list of all current templates. Import/execution semantics and
typed ports are documented in [CLI](cli.md) and [protocol](protocol.md).
The fast documentation tests compare the table to the actual JSON inventory;
they do not regenerate scientific outputs or assert wet-lab acceptance.
