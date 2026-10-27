# Interfaces and Command Discovery

<!-- manual-counterparts:start -->
Related GUI reference: [Detailed Reference](../gui/detailed-reference.md) | [GENtle Shell (GUI)](../gui/gentle-shell-gui.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## Overview

GENtle currently provides seven binaries:

- `gentle`: graphical desktop app (GUI)
- `gentle_cli`: JSON operation/workflow CLI for automation and AI tools
- `gentle_js`: interactive JavaScript shell (optional build target; requires
  feature `js-interface`)
- `gentle_lua`: interactive Lua shell (optional build target; requires feature
  `lua-interface`)
- `gentle_examples_docs`: generates adapter snippets and tutorial artifacts from canonical protocol examples;
  `--version` reports the shared package/build/source identity for candidate-bound acceptance
- `gentle_mcp`: MCP stdio server (guarded mutating + UI-intent parity baseline;
  includes standardized capability discovery via `tools/list`,
  `capabilities`, and `help`)
- `gentle_publication_report`: manifest-driven responsive HTML and printable
  gene-set report bundle generator; `--version`/`-V` reports the shared
  package/build/source identity without loading inputs or creating output

In addition, the GUI includes an embedded `Shell` panel that uses the same
shared shell parser/executor as `gentle_cli shell`.

Structured command glossary:

- Built-in per-command help (`help ...`) and help catalog export
  (`help --format json|markdown`) are rendered from the shared protocol
  `capability_registry`. Glossary command descriptors carry the
  glossary-derived usage, interface, alias, mutation-class, adapter-surfacing,
  and schema metadata together. During the migration, `docs/glossary.json`
  remains the compile-time seed for those command descriptors.
- `gentle_cli capabilities` includes the shared protocol
  `capability_registry`, so CLI/MCP/JavaScript/Lua discovery rows use the same
  names, descriptions, mutation classes, and adapter projections.
- Shell-only dry-run catalog preview routes are available through the shared
  shell, for example:
  - `gentle_cli shell 'genomes preview-ensembl-specs --catalog assets/genomes.json'`
  - `gentle_cli shell 'helpers preview-ensembl-specs --catalog assets/helper_genomes.json'`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Next: Build and Run](build.md)
<!-- manual-page-nav:end -->
