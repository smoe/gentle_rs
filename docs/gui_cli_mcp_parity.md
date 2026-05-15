# GUI / CLI / MCP Parity Matrix

Generated automatically from the protocol capability registry.

## Method

This file is a projection of `gentle_protocol::capability_registry()`. The cell vocabulary is the architecture-doc parity policy:

- `prominent`: `AdapterSurfacing::Prominent` for this adapter
- `shell-only`: `AdapterSurfacing::ShellPassthrough`
- `n/a`: `AdapterSurfacing::NotApplicable`, with a one-line Notes justification
- `gap`: operation should be reachable but currently is not

Only `gap` signals implementation work. Human-readable Notes are populated from each descriptor's `surfacing_justifications` map.

## Findings

| Adapter | prominent | shell-only | gap |
|---|---:|---:|---:|
| GUI | 20 | 529 | 25 |
| gentle_cli | 315 | 234 | 25 |
| MCP | 41 | 206 | 188 |
| JS | 39 | 206 | 181 |
| Lua | 38 | 206 | 182 |
| ClawBio | 0 | 0 | 0 |

## Glossary Commands

| Capability | Source | GUI | gentle_cli | MCP | JS | Lua | ClawBio | Notes |
|---|---|---|---|---|---|---|---|---|
| help | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| capabilities | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| state-summary | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| history status | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| history undo | glossary-command | prominent | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| history redo | glossary-command | prominent | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| load-project | glossary-command | prominent | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| save-project | glossary-command | prominent | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| screenshot-window | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| render-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| render-dotplot-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| inspect-feature-expert | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| render-feature-expert-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| panels import-isoform | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| panels inspect-isoform | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| panels render-isoform-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| panels validate-isoform | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| render-rna-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-info | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| render-lineage-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| protocol-cartoon list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| protocol-cartoon render-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| protocol-cartoon render-template-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| protocol-cartoon template-validate | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| protocol-cartoon render-with-bindings | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| protocol-cartoon template-export | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| gibson preview | glossary-command | prominent | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| gibson apply | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cache inspect | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| cache clear | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| render-pool-gel-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| arrange-serial | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| arrange-set-ladders | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks create-from-arrangement | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks place-arrangement | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks move | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks move-samples | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks move-blocks | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| racks labels-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks fabrication-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks isometric-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks openscad | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks carrier-labels-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks simulation-json | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks set-profile | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks apply-template | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks set-fill-direction | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks set-custom-profile | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| racks set-blocked | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| ladders list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ladders export | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| proteases list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| proteases show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| proteases digest | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| proteases digest-gel-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| export-pool | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| export-run-bundle | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import-pool | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| agents list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| agents preflight | glossary-command | shell-only | prominent | prominent | n/a | n/a | n/a | JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| agents discover-models | glossary-command | shell-only | prominent | prominent | n/a | n/a | n/a | JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| agents ask | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| agents plan | glossary-command | shell-only | prominent | prominent | n/a | n/a | n/a | JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| agents execute-plan | glossary-command | shell-only | prominent | prominent | n/a | n/a | n/a | JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ui intents | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ui open | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ui focus | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ui prepared-genomes | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ui latest-prepared | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genbank fetch | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| dbsnp fetch | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| variant annotate-promoters | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| variant promoter-context | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| variant reporter-fragments | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| variant materialize-allele | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| uniprot fetch | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| uniprot import-swissprot | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| uniprot list | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| uniprot show | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot map | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| uniprot projection-list | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot projection-show | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot feature-coding-dna | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot resolve-ensembl-links | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot transcript-accounting | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot compare-ensembl-exons | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot compare-ensembl-peptide | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot audit-projection | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| uniprot audit-parity | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| uniprot audit-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot audit-show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot audit-export | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot audit-parity-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot audit-parity-show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| uniprot audit-parity-export | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ensembl-gene fetch | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| ensembl-gene list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ensembl-gene show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ensembl-gene import-sequence | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| ensembl-region fetch | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| ensembl-protein fetch | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| ensembl-protein list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ensembl-protein show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ensembl-protein import-sequence | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| resources sync-rebase | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources sync-jaspar | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources sync-ucsc-rmsk | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources install-ucsc-rmsk | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources prepare-ucsc-rmsk-index | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources suggest-ucsc-rmsk-index | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources sync-jaspar-remote-metadata | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| resources summarize-jaspar | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| resources benchmark-jaspar | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| resources list-jaspar | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| resources resolve-tf-query | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| gene-groups list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| gene-groups show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| gene-groups resolve | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| gene-groups doctor | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| gene-groups draft | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources list-publication-datasets | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources status-publication-dataset | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources prepare-publication-dataset | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| resources inspect-jaspar | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| services status | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| services providers list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| services providers doctor | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| services project-preflight | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| services project-quote | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| services handoff | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| services guide | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes ensembl-available | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes install-ensembl | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes validate-catalog | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes preview-ensembl-specs | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes update-ensembl-specs | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes status | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes genes | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes prepare | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| genomes remove-prepared | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes remove-catalog-entry | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes blast | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes blast-start | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes blast-status | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes blast-cancel | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes blast-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| genomes blast-track | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| genomes extract-region | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| genomes extract-gene | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| genomes extract-promoter | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| genomes extend-anchor | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| genomes verify-anchor | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| helpers list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers vocabulary list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers vocabulary doctor | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers ensembl-available | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers install-ensembl | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| hosts list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers validate-catalog | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers preview-ensembl-specs | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers update-ensembl-specs | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers status | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers genes | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers prepare | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| helpers remove-prepared | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers remove-catalog-entry | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers blast | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers blast-start | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers blast-status | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers blast-cancel | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers blast-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| helpers blast-track | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| helpers extract-region | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| helpers extract-gene | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| helpers extract-promoter | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| helpers extend-anchor | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| helpers verify-anchor | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun list | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun status | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun prepare | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun project | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun interpret | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun list-read-reports | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun show-read-report | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun export-coverage | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| cutrun inspect-regulatory-support | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| tracks import-bed | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| tracks import-bigwig | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| tracks import-vcf | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| arrays inspect-microarray-track | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| arrays project-microarray-track | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| tracks tracked list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| tracks tracked add | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| tracks tracked remove | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| tracks tracked clear | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| tracks tracked apply | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| guides list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| guides show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| guides put | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| guides delete | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| guides filter | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| guides filter-show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| guides oligos-generate | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| guides oligos-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| guides oligos-show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| guides oligos-export | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| guides protocol-export | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| primers design | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| primers design-qpcr | glossary-command | prominent | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| primers test-cdna-pcr | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| primers test-cdna-qpcr | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| primers test-cdna-qpcr-fasta | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| primers prepare-restriction-cloning | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| primers seed-restriction-cloning-handoff | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers restriction-cloning-vector-suggestions | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers list-restriction-cloning-handoffs | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers show-restriction-cloning-handoff | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers export-restriction-cloning-handoff | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers preflight | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers seed-from-feature | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers seed-from-splicing | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers list-reports | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers show-report | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers export-report | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers list-qpcr-reports | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers show-qpcr-report | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| primers export-qpcr-report | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| features tfbs-summary | glossary-command | shell-only | shell-only | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| features repeat-query | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| features repeat-overlaps | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| features materialize-repeats | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| features repeat-cohort | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| features window-cohort-tfbs | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| features tfbs-score-tracks-svg | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| dotplot compute | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| dotplot overlay-compute | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| dotplot list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| dotplot show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| transcripts derive | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| transcripts exon-skip-plan | glossary-command | shell-only | prominent | prominent | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| transcripts exon-skip-materialize | glossary-command | shell-only | prominent | prominent | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| transcripts residue-genomic-coordinates | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| flex compute | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| flex list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| flex show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| splicing-refs derive | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| align compute | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-trace import | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-trace list | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-trace show | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-confirm run | glossary-command | prominent | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-confirm list-reports | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-confirm show-report | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-confirm export-report | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-confirm export-support-tsv | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| seq-primer suggest | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reverse-translate run | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reverse-translate list-reports | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reverse-translate show-report | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reverse-translate export-report | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| construct-reasoning build-protein-dna-handoff | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| construct-reasoning list-graphs | glossary-command | shell-only | prominent | prominent | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| construct-reasoning show-graph | glossary-command | shell-only | prominent | prominent | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| construct-reasoning set-annotation-status | glossary-command | shell-only | prominent | prominent | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| construct-reasoning write-annotation | glossary-command | shell-only | prominent | prominent | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| construct-reasoning export-graph | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reads acquire status | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reads acquire prepare | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reads acquire inspect | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| reads acquire cancel | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads interpret | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads batch-map | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads align-report | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads preflight-isoforms | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads list-reports | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| rna-reads show-report | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| rna-reads show-alignment | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| rna-reads summarize-gene-support | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads inspect-gene-support | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads inspect-alignments | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| rna-reads inspect-concatemers | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| rna-reads build-transcript-index | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| rna-reads materialize-hits | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-report | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-hits-fasta | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-sample-sheet | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-target-quality | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-paths-tsv | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-abundance-tsv | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-score-density-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-alignments-tsv | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-isoform-triage-tsv | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| rna-reads export-alignment-dotplot-svg | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| candidates delete | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates generate | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates generate-between-anchors | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| candidates metrics | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| candidates score | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates score-distance | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates score-weighted | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates top-k | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates pareto | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates filter | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates set-op | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| macros run | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| macros instance-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| macros instance-show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| macros template-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| macros template-show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| macros template-put | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| macros template-delete | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| macros template-import | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| macros template-run | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| routines list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| routines explain | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| routines compare | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning profile show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning profile set | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning profile clear | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning objective show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning objective set | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning objective clear | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning suggestions list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning suggestions accept | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning suggestions reject | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning sync status | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning sync pull | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| planning sync push | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| candidates macro | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates template-list | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| candidates template-show | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| candidates template-put | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates template-delete | glossary-command | shell-only | prominent | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| candidates template-run | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| set-param | glossary-command | shell-only | shell-only | gap | gap | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| op | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| workflow | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| batch plan | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| batch run | glossary-command | shell-only | prominent | n/a | n/a | n/a | n/a | MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| apply_operation | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| apply_workflow | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ask_agent_system | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| blast_helper_genome | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| blast_reference_genome | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| capabilities | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| digest | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| export_dna_ladders | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| export_rna_ladders | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| extend_genome_anchor | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| extract_genome_gene | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| extract_genome_region | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_genome_bed_track | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_genome_bigwig_track | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_genome_vcf_track | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_pool | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| inspect_dna_ladders | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| inspect_rna_ladders | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| is_reference_genome_prepared | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_agent_systems | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_dna_ladders | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_reference_genome_genes | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_reference_genomes | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_reference_catalog_entries | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_helper_catalog_entries | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_host_profile_catalog_entries | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_ensembl_installable_genomes | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_rna_ladders | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| load_dna | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| load_project | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| prepare_genome | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| render_pool_gel_svg | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| save_project | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| set_parameter | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| set_vcf_display_filter | glossary-command | gap | gap | gap | prominent | gap | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| state_summary | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| sync_jaspar | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| sync_rebase | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| write_gb | glossary-command | n/a | n/a | n/a | prominent | n/a | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>Lua: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| apply_operation | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| apply_workflow | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| ask_agent_system | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| blast_helper_genome | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| blast_reference_genome | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| capabilities | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| export_dna_ladders | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| export_rna_ladders | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| extend_genome_anchor | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| extract_genome_gene | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| extract_genome_region | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_genome_bed_track | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_genome_bigwig_track | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_genome_vcf_track | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| import_pool | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| inspect_dna_ladders | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| inspect_rna_ladders | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| is_reference_genome_prepared | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_agent_systems | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_dna_ladders | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_reference_genome_genes | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_reference_genomes | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_reference_catalog_entries | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_helper_catalog_entries | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_host_profile_catalog_entries | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_ensembl_installable_genomes | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| list_rna_ladders | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| load_dna | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| load_project | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| prepare_genome | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| render_pool_gel_svg | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| save_project | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| set_parameter | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| set_vcf_display_filter | glossary-command | gap | gap | gap | gap | prominent | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw glossary command rows. |
| state_summary | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| sync_jaspar | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| sync_rebase | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |
| write_gb | glossary-command | n/a | n/a | n/a | n/a | prominent | n/a | GUI: Local CLI/GUI command, not an engine route.<br>gentle_cli: Local CLI/GUI command, not an engine route.<br>MCP: Local CLI/GUI command, not an engine route.<br>JS: Local CLI/GUI command, not an engine route.<br>ClawBio: Local CLI/GUI command, not an engine route. |

## Engine Operations

| Capability | Source | GUI | gentle_cli | MCP | JS | Lua | ClawBio | Notes |
|---|---|---|---|---|---|---|---|---|
| LoadFile | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SaveFile | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderSequenceSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderDotplotSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderTfbsScoreTracksSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderTfbsScoreTrackCorrelationSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderFeatureExpertSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderIsoformArchitectureSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderRnaStructureSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderLineageSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderPoolGelSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderProteinGelSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderProteinGelReportsSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderProteaseDigestGelSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderProtein2dGelSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderProtocolCartoonSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderProtocolCartoonTemplateSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ValidateProtocolCartoonTemplate | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderProtocolCartoonTemplateWithBindingsSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportProtocolCartoonTemplateJson | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ApplyGibsonAssemblyPlan | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| CreateArrangementSerial | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetArrangementLadders | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetContainerDeclaredContentsExclusive | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| CreateRackFromArrangement | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PlaceArrangementOnRack | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MoveRackPlacement | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MoveRackSamples | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MoveRackArrangementBlocks | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetRackProfile | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ApplyRackTemplate | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetRackFillDirection | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetRackProfileCustom | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetRackBlockedCoordinates | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRackLabelsSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRackFabricationSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRackIsometricSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRackOpenScad | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRackCarrierLabelsSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRackSimulationJson | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportDnaLadders | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaLadders | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportPool | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportProcessRunBundle | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportLabAssistantInstructions | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PrepareGenome | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExtractGenomeRegion | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExtractGenomeGene | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExtractGenomePromoterSlice | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExtendGenomeAnchor | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportGenomeBedTrack | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportGenomeBigWigTrack | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportGenomeVcfTrack | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ProjectMicroarrayTrack | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ProjectGenomeInterval | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ListCutRunDatasets | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ShowCutRunDatasetStatus | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PrepareCutRunDataset | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ProjectCutRunDataset | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| InterpretCutRunReads | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ListCutRunReadReports | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ShowCutRunReadReport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportCutRunReadCoverage | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| InspectCutRunRegulatorySupport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportIsoformPanel | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportUniprotSwissProt | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FetchUniprotSwissProt | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FetchEnsemblGene | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FetchEnsemblRegion | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FetchEnsemblProtein | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FetchGenBankAccession | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FetchDbSnpRegion | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FetchUniprotLinkedGenBank | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportUniprotEntrySequence | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportEnsemblGeneSequence | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportEnsemblProteinSequence | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ProjectUniprotToGenome | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| QueryProteinResidueGenomicCoordinates | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| AuditUniprotProjectionConsistency | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| AuditUniprotProjectionParity | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ImportBlastHitsTrack | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DigestContainer | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MergeContainersById | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| LigationContainer | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FilterContainerByMolecularWeight | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| Digest | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FindRestrictionSites | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| QueryRepeatAnnotations | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| QueryRepeatOverlaps | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MaterializeRepeatFeatures | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| BuildRepeatEnvironmentCohort | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MergeContainers | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| Ligation | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| Pcr | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PcrAdvanced | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PcrMutagenesis | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DesignPrimerPairs | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DesignInsertionPrimerPairs | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportPrimerDesignReport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| AssessPrimerPairSpecificity | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PrepareRestrictionCloningPcrHandoff | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PcrOverlapExtensionMutagenesis | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DesignQpcrAssays | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| TestCdnaPcr | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| TestCdnaQpcr | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| TestCdnaQpcrFasta | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DeriveTranscriptSequences | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PlanExonSkippedIsoform | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MaterializeExonSkippedIsoform | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DeriveProteinSequences | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ReverseTranslateProteinSequence | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ProteaseDigestProteinSequence | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| BuildProteinToDnaHandoffReasoning | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ComputeDotplot | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ComputeDotplotOverlay | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ComputeFlexibilityTrack | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DeriveSplicingReferences | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| AlignSequences | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ConfirmConstructReads | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SuggestSequencingPrimers | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ListSequencingConfirmationReports | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ShowSequencingConfirmationReport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportSequencingConfirmationReport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportSequencingConfirmationSupportTsv | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ReadAcquireStatus | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ReadAcquirePrepare | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ReadAcquireInspect | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ReadAcquireCancel | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| InterpretRnaReads | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| AlignRnaReadReport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| PreflightRnaReadIsoforms | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ListRnaReadReports | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ShowRnaReadReport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeRnaReadGeneSupport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| InspectRnaReadGeneSupport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadIsoformTriageTsv | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RunRnaReadBatchMap | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeTfbsRegion | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeTfbsScoreTracks | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeTfbsTrackSimilarity | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeAlternativePromoterComparison | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizePromoterEvidenceMatrix | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeIsoformPromoterComparison | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizePromoterExpressionEvidence | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportPromoterArtifactManifest | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeMultiGenePromoterTfbs | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RenderMultiGenePromoterTfbsSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ScanTfbsHits | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| InspectJasparEntry | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeJasparEntries | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ResolveTfQueries | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| BenchmarkJasparRegistry | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ListJasparCatalog | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SyncJasparRemoteMetadata | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| AnnotatePromoterWindows | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SummarizeVariantPromoterContext | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SuggestPromoterReporterFragments | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MaterializeVariantAllele | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadReport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadHitsFasta | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadSampleSheet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadTargetQuality | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadExonPathsTsv | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadExonAbundanceTsv | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadScoreDensitySvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadAlignmentsTsv | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportRnaReadAlignmentDotplotSvg | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| MaterializeRnaReadHitSequences | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExtractRegion | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExtractAnchoredRegion | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SelectCandidate | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FilterByMolecularWeight | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FilterByDesignConstraints | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| GenerateCandidateSet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| GenerateCandidateSetBetweenAnchors | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DeleteCandidateSet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| UpsertGuideSet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DeleteGuideSet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FilterGuidesPractical | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| GenerateGuideOligos | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportGuideOligos | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportGuideProtocolText | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportFeaturesBed | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| InspectSequenceContextView | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ExportSequenceContextBundle | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ScoreCandidateSetExpression | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ScoreCandidateSetDistance | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| FilterCandidateSet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| CandidateSetOp | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ScoreCandidateSetWeightedObjective | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| TopKCandidateSet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ParetoFrontierCandidateSet | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| UpsertWorkflowMacroTemplate | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DeleteWorkflowMacroTemplate | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| UpsertCandidateMacroTemplate | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| DeleteCandidateMacroTemplate | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| Reverse | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| Complement | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| ReverseComplement | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| Branch | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetDisplayVisibility | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetLinearViewport | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetTopology | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| RecomputeFeatures | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| SetParameter | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |
| AnnotateTfbs | engine-operation | shell-only | shell-only | shell-only | shell-only | shell-only | n/a | ClawBio: ClawBio exposes curated skill intents rather than raw engine operation rows. |

## MCP Tools

| Capability | Source | GUI | gentle_cli | MCP | JS | Lua | ClawBio | Notes |
|---|---|---|---|---|---|---|---|---|
| capabilities | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| state_summary | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| exon_skip_plan | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| exon_skip_materialize | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| restriction_site_detail | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| agent_systems | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| agent_preflight | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| agent_models | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| agent_plan | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| agent_execute_plan | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| op | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| workflow | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| help | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| reference_catalog_entries | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| helper_catalog_entries | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| helper_semantics_vocabulary | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| host_profile_catalog_entries | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| ensembl_installable_genomes | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| construct_reasoning_graphs | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| construct_reasoning_graph | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| construct_reasoning_set_annotation_status | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| construct_reasoning_write_annotation | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| helper_interpretation | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| ui_intents | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| ui_intent | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| ui_prepared_genomes | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| ui_latest_prepared | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| blast_async_start | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| blast_async_status | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| blast_async_cancel | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |
| blast_async_list | mcp-tool | n/a | n/a | prominent | n/a | n/a | n/a | GUI: MCP tools/list descriptors are transport-specific metadata, not GUI actions.<br>gentle_cli: MCP tools/list descriptors are transport-specific metadata, not CLI subcommands.<br>JS: MCP tools/list descriptors are transport-specific metadata, not JavaScript wrapper functions.<br>Lua: MCP tools/list descriptors are transport-specific metadata, not Lua wrapper functions.<br>ClawBio: ClawBio consumes its skill intent descriptor rather than MCP tools/list metadata. |

## Open Gaps

| Capability | Source | Adapter | Engine operations |
|---|---|---|---|
| render-svg | glossary-command | MCP | RenderSequenceSvg |
| render-svg | glossary-command | JS | RenderSequenceSvg |
| render-svg | glossary-command | Lua | RenderSequenceSvg |
| render-dotplot-svg | glossary-command | MCP | RenderDotplotSvg |
| render-dotplot-svg | glossary-command | JS | RenderDotplotSvg |
| render-dotplot-svg | glossary-command | Lua | RenderDotplotSvg |
| render-feature-expert-svg | glossary-command | MCP | RenderFeatureExpertSvg |
| render-feature-expert-svg | glossary-command | JS | RenderFeatureExpertSvg |
| render-feature-expert-svg | glossary-command | Lua | RenderFeatureExpertSvg |
| panels import-isoform | glossary-command | MCP | ImportIsoformPanel |
| panels import-isoform | glossary-command | JS | ImportIsoformPanel |
| panels import-isoform | glossary-command | Lua | ImportIsoformPanel |
| panels render-isoform-svg | glossary-command | MCP | RenderIsoformArchitectureSvg |
| panels render-isoform-svg | glossary-command | JS | RenderIsoformArchitectureSvg |
| panels render-isoform-svg | glossary-command | Lua | RenderIsoformArchitectureSvg |
| render-rna-svg | glossary-command | MCP | RenderRnaStructureSvg |
| render-rna-svg | glossary-command | JS | RenderRnaStructureSvg |
| render-rna-svg | glossary-command | Lua | RenderRnaStructureSvg |
| render-lineage-svg | glossary-command | MCP | RenderLineageSvg |
| render-lineage-svg | glossary-command | JS | RenderLineageSvg |
| render-lineage-svg | glossary-command | Lua | RenderLineageSvg |
| protocol-cartoon render-svg | glossary-command | MCP | RenderProtocolCartoonSvg |
| protocol-cartoon render-svg | glossary-command | JS | RenderProtocolCartoonSvg |
| protocol-cartoon render-svg | glossary-command | Lua | RenderProtocolCartoonSvg |
| protocol-cartoon render-template-svg | glossary-command | MCP | RenderProtocolCartoonTemplateSvg |
| protocol-cartoon render-template-svg | glossary-command | JS | RenderProtocolCartoonTemplateSvg |
| protocol-cartoon render-template-svg | glossary-command | Lua | RenderProtocolCartoonTemplateSvg |
| protocol-cartoon template-validate | glossary-command | MCP | ValidateProtocolCartoonTemplate |
| protocol-cartoon template-validate | glossary-command | JS | ValidateProtocolCartoonTemplate |
| protocol-cartoon template-validate | glossary-command | Lua | ValidateProtocolCartoonTemplate |
| protocol-cartoon render-with-bindings | glossary-command | MCP | RenderProtocolCartoonTemplateWithBindingsSvg |
| protocol-cartoon render-with-bindings | glossary-command | JS | RenderProtocolCartoonTemplateWithBindingsSvg |
| protocol-cartoon render-with-bindings | glossary-command | Lua | RenderProtocolCartoonTemplateWithBindingsSvg |
| protocol-cartoon template-export | glossary-command | MCP | ExportProtocolCartoonTemplateJson |
| protocol-cartoon template-export | glossary-command | JS | ExportProtocolCartoonTemplateJson |
| protocol-cartoon template-export | glossary-command | Lua | ExportProtocolCartoonTemplateJson |
| gibson apply | glossary-command | MCP | ApplyGibsonAssemblyPlan |
| gibson apply | glossary-command | JS | ApplyGibsonAssemblyPlan |
| gibson apply | glossary-command | Lua | ApplyGibsonAssemblyPlan |
| render-pool-gel-svg | glossary-command | MCP | RenderPoolGelSvg |
| render-pool-gel-svg | glossary-command | JS | RenderPoolGelSvg |
| render-pool-gel-svg | glossary-command | Lua | RenderPoolGelSvg |
| arrange-serial | glossary-command | MCP | CreateArrangementSerial |
| arrange-serial | glossary-command | JS | CreateArrangementSerial |
| arrange-serial | glossary-command | Lua | CreateArrangementSerial |
| arrange-set-ladders | glossary-command | MCP | SetArrangementLadders |
| arrange-set-ladders | glossary-command | JS | SetArrangementLadders |
| arrange-set-ladders | glossary-command | Lua | SetArrangementLadders |
| racks create-from-arrangement | glossary-command | MCP | CreateRackFromArrangement |
| racks create-from-arrangement | glossary-command | JS | CreateRackFromArrangement |
| racks create-from-arrangement | glossary-command | Lua | CreateRackFromArrangement |
| racks place-arrangement | glossary-command | MCP | PlaceArrangementOnRack |
| racks place-arrangement | glossary-command | JS | PlaceArrangementOnRack |
| racks place-arrangement | glossary-command | Lua | PlaceArrangementOnRack |
| racks move | glossary-command | MCP | MoveRackPlacement |
| racks move | glossary-command | JS | MoveRackPlacement |
| racks move | glossary-command | Lua | MoveRackPlacement |
| racks move-samples | glossary-command | MCP | MoveRackSamples |
| racks move-samples | glossary-command | JS | MoveRackSamples |
| racks move-samples | glossary-command | Lua | MoveRackSamples |
| racks move-blocks | glossary-command | MCP | MoveRackArrangementBlocks |
| racks move-blocks | glossary-command | JS | MoveRackArrangementBlocks |
| racks move-blocks | glossary-command | Lua | MoveRackArrangementBlocks |
| racks labels-svg | glossary-command | MCP | ExportRackLabelsSvg |
| racks labels-svg | glossary-command | JS | ExportRackLabelsSvg |
| racks labels-svg | glossary-command | Lua | ExportRackLabelsSvg |
| racks fabrication-svg | glossary-command | MCP | ExportRackFabricationSvg |
| racks fabrication-svg | glossary-command | JS | ExportRackFabricationSvg |
| racks fabrication-svg | glossary-command | Lua | ExportRackFabricationSvg |
| racks isometric-svg | glossary-command | MCP | ExportRackIsometricSvg |
| racks isometric-svg | glossary-command | JS | ExportRackIsometricSvg |
| racks isometric-svg | glossary-command | Lua | ExportRackIsometricSvg |
| racks openscad | glossary-command | MCP | ExportRackOpenScad |
| racks openscad | glossary-command | JS | ExportRackOpenScad |
| racks openscad | glossary-command | Lua | ExportRackOpenScad |
| racks carrier-labels-svg | glossary-command | MCP | ExportRackCarrierLabelsSvg |
| racks carrier-labels-svg | glossary-command | JS | ExportRackCarrierLabelsSvg |
| racks carrier-labels-svg | glossary-command | Lua | ExportRackCarrierLabelsSvg |
| racks simulation-json | glossary-command | MCP | ExportRackSimulationJson |
| racks simulation-json | glossary-command | JS | ExportRackSimulationJson |
| racks simulation-json | glossary-command | Lua | ExportRackSimulationJson |
| racks set-profile | glossary-command | MCP | SetRackProfile |
| racks set-profile | glossary-command | JS | SetRackProfile |
| racks set-profile | glossary-command | Lua | SetRackProfile |
| racks apply-template | glossary-command | MCP | ApplyRackTemplate |
| racks apply-template | glossary-command | JS | ApplyRackTemplate |
| racks apply-template | glossary-command | Lua | ApplyRackTemplate |
| racks set-fill-direction | glossary-command | MCP | SetRackFillDirection |
| racks set-fill-direction | glossary-command | JS | SetRackFillDirection |
| racks set-fill-direction | glossary-command | Lua | SetRackFillDirection |
| racks set-custom-profile | glossary-command | MCP | SetRackProfileCustom |
| racks set-custom-profile | glossary-command | JS | SetRackProfileCustom |
| racks set-custom-profile | glossary-command | Lua | SetRackProfileCustom |
| racks set-blocked | glossary-command | MCP | SetRackBlockedCoordinates |
| racks set-blocked | glossary-command | JS | SetRackBlockedCoordinates |
| racks set-blocked | glossary-command | Lua | SetRackBlockedCoordinates |
| ladders export | glossary-command | MCP | ExportDnaLadders, ExportRnaLadders |
| ladders export | glossary-command | JS | ExportDnaLadders, ExportRnaLadders |
| ladders export | glossary-command | Lua | ExportDnaLadders, ExportRnaLadders |
| proteases digest | glossary-command | MCP | ProteaseDigestProteinSequence |
| proteases digest | glossary-command | JS | ProteaseDigestProteinSequence |
| proteases digest | glossary-command | Lua | ProteaseDigestProteinSequence |
| proteases digest-gel-svg | glossary-command | MCP | RenderProteaseDigestGelSvg |
| proteases digest-gel-svg | glossary-command | JS | RenderProteaseDigestGelSvg |
| proteases digest-gel-svg | glossary-command | Lua | RenderProteaseDigestGelSvg |
| export-pool | glossary-command | MCP | ExportPool |
| export-pool | glossary-command | JS | ExportPool |
| export-pool | glossary-command | Lua | ExportPool |
| export-run-bundle | glossary-command | MCP | ExportProcessRunBundle |
| export-run-bundle | glossary-command | JS | ExportProcessRunBundle |
| export-run-bundle | glossary-command | Lua | ExportProcessRunBundle |
| genbank fetch | glossary-command | MCP | FetchGenBankAccession |
| genbank fetch | glossary-command | JS | FetchGenBankAccession |
| genbank fetch | glossary-command | Lua | FetchGenBankAccession |
| dbsnp fetch | glossary-command | MCP | FetchDbSnpRegion |
| dbsnp fetch | glossary-command | JS | FetchDbSnpRegion |
| dbsnp fetch | glossary-command | Lua | FetchDbSnpRegion |
| variant annotate-promoters | glossary-command | MCP | AnnotatePromoterWindows |
| variant annotate-promoters | glossary-command | JS | AnnotatePromoterWindows |
| variant annotate-promoters | glossary-command | Lua | AnnotatePromoterWindows |
| variant promoter-context | glossary-command | MCP | SummarizeVariantPromoterContext |
| variant promoter-context | glossary-command | JS | SummarizeVariantPromoterContext |
| variant promoter-context | glossary-command | Lua | SummarizeVariantPromoterContext |
| variant reporter-fragments | glossary-command | MCP | SuggestPromoterReporterFragments |
| variant reporter-fragments | glossary-command | JS | SuggestPromoterReporterFragments |
| variant reporter-fragments | glossary-command | Lua | SuggestPromoterReporterFragments |
| variant materialize-allele | glossary-command | MCP | MaterializeVariantAllele |
| variant materialize-allele | glossary-command | JS | MaterializeVariantAllele |
| variant materialize-allele | glossary-command | Lua | MaterializeVariantAllele |
| uniprot fetch | glossary-command | MCP | FetchUniprotSwissProt |
| uniprot fetch | glossary-command | JS | FetchUniprotSwissProt |
| uniprot fetch | glossary-command | Lua | FetchUniprotSwissProt |
| uniprot import-swissprot | glossary-command | MCP | ImportUniprotSwissProt |
| uniprot import-swissprot | glossary-command | JS | ImportUniprotSwissProt |
| uniprot import-swissprot | glossary-command | Lua | ImportUniprotSwissProt |
| uniprot list | glossary-command | MCP | SummarizeTfbsRegion |
| uniprot list | glossary-command | JS | SummarizeTfbsRegion |
| uniprot list | glossary-command | Lua | SummarizeTfbsRegion |
| uniprot map | glossary-command | MCP | ProjectUniprotToGenome |
| uniprot map | glossary-command | JS | ProjectUniprotToGenome |
| uniprot map | glossary-command | Lua | ProjectUniprotToGenome |
| uniprot audit-projection | glossary-command | MCP | AuditUniprotProjectionConsistency |
| uniprot audit-projection | glossary-command | JS | AuditUniprotProjectionConsistency |
| uniprot audit-projection | glossary-command | Lua | AuditUniprotProjectionConsistency |
| uniprot audit-parity | glossary-command | MCP | AuditUniprotProjectionParity |
| uniprot audit-parity | glossary-command | JS | AuditUniprotProjectionParity |
| uniprot audit-parity | glossary-command | Lua | AuditUniprotProjectionParity |
| ensembl-gene fetch | glossary-command | MCP | FetchEnsemblGene |
| ensembl-gene fetch | glossary-command | JS | FetchEnsemblGene |
| ensembl-gene fetch | glossary-command | Lua | FetchEnsemblGene |
| ensembl-gene import-sequence | glossary-command | MCP | ImportEnsemblGeneSequence |
| ensembl-gene import-sequence | glossary-command | JS | ImportEnsemblGeneSequence |
| ensembl-gene import-sequence | glossary-command | Lua | ImportEnsemblGeneSequence |
| ensembl-region fetch | glossary-command | MCP | FetchEnsemblRegion |
| ensembl-region fetch | glossary-command | JS | FetchEnsemblRegion |
| ensembl-region fetch | glossary-command | Lua | FetchEnsemblRegion |
| ensembl-protein fetch | glossary-command | MCP | FetchEnsemblProtein |
| ensembl-protein fetch | glossary-command | JS | FetchEnsemblProtein |
| ensembl-protein fetch | glossary-command | Lua | FetchEnsemblProtein |
| ensembl-protein import-sequence | glossary-command | MCP | ImportEnsemblProteinSequence |
| ensembl-protein import-sequence | glossary-command | JS | ImportEnsemblProteinSequence |
| ensembl-protein import-sequence | glossary-command | Lua | ImportEnsemblProteinSequence |
| resources sync-jaspar-remote-metadata | glossary-command | MCP | SyncJasparRemoteMetadata |
| resources sync-jaspar-remote-metadata | glossary-command | JS | SyncJasparRemoteMetadata |
| resources sync-jaspar-remote-metadata | glossary-command | Lua | SyncJasparRemoteMetadata |
| resources summarize-jaspar | glossary-command | MCP | SummarizeJasparEntries |
| resources summarize-jaspar | glossary-command | JS | SummarizeJasparEntries |
| resources summarize-jaspar | glossary-command | Lua | SummarizeJasparEntries |
| resources benchmark-jaspar | glossary-command | MCP | BenchmarkJasparRegistry |
| resources benchmark-jaspar | glossary-command | JS | BenchmarkJasparRegistry |
| resources benchmark-jaspar | glossary-command | Lua | BenchmarkJasparRegistry |
| resources list-jaspar | glossary-command | MCP | ListJasparCatalog |
| resources list-jaspar | glossary-command | JS | ListJasparCatalog |
| resources list-jaspar | glossary-command | Lua | ListJasparCatalog |
| resources resolve-tf-query | glossary-command | MCP | ResolveTfQueries |
| resources resolve-tf-query | glossary-command | JS | ResolveTfQueries |
| resources resolve-tf-query | glossary-command | Lua | ResolveTfQueries |
| resources inspect-jaspar | glossary-command | MCP | InspectJasparEntry |
| resources inspect-jaspar | glossary-command | JS | InspectJasparEntry |
| resources inspect-jaspar | glossary-command | Lua | InspectJasparEntry |
| genomes prepare | glossary-command | MCP | PrepareReferenceGenome |
| genomes prepare | glossary-command | JS | PrepareReferenceGenome |
| genomes prepare | glossary-command | Lua | PrepareReferenceGenome |
| genomes blast-track | glossary-command | MCP | ImportBlastHitsTrack |
| genomes blast-track | glossary-command | JS | ImportBlastHitsTrack |
| genomes blast-track | glossary-command | Lua | ImportBlastHitsTrack |
| genomes extract-region | glossary-command | MCP | ExtractGenomeRegion |
| genomes extract-region | glossary-command | JS | ExtractGenomeRegion |
| genomes extract-region | glossary-command | Lua | ExtractGenomeRegion |
| genomes extract-gene | glossary-command | MCP | ExtractGenomeGene |
| genomes extract-gene | glossary-command | JS | ExtractGenomeGene |
| genomes extract-gene | glossary-command | Lua | ExtractGenomeGene |
| genomes extract-promoter | glossary-command | MCP | ExtractGenomePromoterSlice |
| genomes extract-promoter | glossary-command | JS | ExtractGenomePromoterSlice |
| genomes extract-promoter | glossary-command | Lua | ExtractGenomePromoterSlice |
| genomes extend-anchor | glossary-command | MCP | ExtendGenomeAnchor |
| genomes extend-anchor | glossary-command | JS | ExtendGenomeAnchor |
| genomes extend-anchor | glossary-command | Lua | ExtendGenomeAnchor |
| genomes verify-anchor | glossary-command | MCP | VerifyGenomeAnchor |
| genomes verify-anchor | glossary-command | JS | VerifyGenomeAnchor |
| genomes verify-anchor | glossary-command | Lua | VerifyGenomeAnchor |
| helpers prepare | glossary-command | MCP | PrepareReferenceGenome |
| helpers prepare | glossary-command | JS | PrepareReferenceGenome |
| helpers prepare | glossary-command | Lua | PrepareReferenceGenome |
| helpers blast-track | glossary-command | MCP | ImportBlastHitsTrack |
| helpers blast-track | glossary-command | JS | ImportBlastHitsTrack |
| helpers blast-track | glossary-command | Lua | ImportBlastHitsTrack |
| helpers extract-region | glossary-command | MCP | ExtractGenomeRegion |
| helpers extract-region | glossary-command | JS | ExtractGenomeRegion |
| helpers extract-region | glossary-command | Lua | ExtractGenomeRegion |
| helpers extract-gene | glossary-command | MCP | ExtractGenomeGene |
| helpers extract-gene | glossary-command | JS | ExtractGenomeGene |
| helpers extract-gene | glossary-command | Lua | ExtractGenomeGene |
| helpers extract-promoter | glossary-command | MCP | ExtractGenomePromoterSlice |
| helpers extract-promoter | glossary-command | JS | ExtractGenomePromoterSlice |
| helpers extract-promoter | glossary-command | Lua | ExtractGenomePromoterSlice |
| helpers extend-anchor | glossary-command | MCP | ExtendGenomeAnchor |
| helpers extend-anchor | glossary-command | JS | ExtendGenomeAnchor |
| helpers extend-anchor | glossary-command | Lua | ExtendGenomeAnchor |
| helpers verify-anchor | glossary-command | MCP | VerifyGenomeAnchor |
| helpers verify-anchor | glossary-command | JS | VerifyGenomeAnchor |
| helpers verify-anchor | glossary-command | Lua | VerifyGenomeAnchor |
| cutrun list | glossary-command | MCP | ListCutRunDatasets |
| cutrun list | glossary-command | JS | ListCutRunDatasets |
| cutrun list | glossary-command | Lua | ListCutRunDatasets |
| cutrun status | glossary-command | MCP | ShowCutRunDatasetStatus |
| cutrun status | glossary-command | JS | ShowCutRunDatasetStatus |
| cutrun status | glossary-command | Lua | ShowCutRunDatasetStatus |
| cutrun prepare | glossary-command | MCP | PrepareCutRunDataset |
| cutrun prepare | glossary-command | JS | PrepareCutRunDataset |
| cutrun prepare | glossary-command | Lua | PrepareCutRunDataset |
| cutrun project | glossary-command | MCP | ProjectCutRunDataset |
| cutrun project | glossary-command | JS | ProjectCutRunDataset |
| cutrun project | glossary-command | Lua | ProjectCutRunDataset |
| cutrun interpret | glossary-command | MCP | InterpretCutRunReads |
| cutrun interpret | glossary-command | JS | InterpretCutRunReads |
| cutrun interpret | glossary-command | Lua | InterpretCutRunReads |
| cutrun list-read-reports | glossary-command | MCP | ListCutRunReadReports |
| cutrun list-read-reports | glossary-command | JS | ListCutRunReadReports |
| cutrun list-read-reports | glossary-command | Lua | ListCutRunReadReports |
| cutrun show-read-report | glossary-command | MCP | ShowCutRunReadReport |
| cutrun show-read-report | glossary-command | JS | ShowCutRunReadReport |
| cutrun show-read-report | glossary-command | Lua | ShowCutRunReadReport |
| cutrun export-coverage | glossary-command | MCP | ExportCutRunReadCoverage |
| cutrun export-coverage | glossary-command | JS | ExportCutRunReadCoverage |
| cutrun export-coverage | glossary-command | Lua | ExportCutRunReadCoverage |
| cutrun inspect-regulatory-support | glossary-command | MCP | InspectCutRunRegulatorySupport |
| cutrun inspect-regulatory-support | glossary-command | JS | InspectCutRunRegulatorySupport |
| cutrun inspect-regulatory-support | glossary-command | Lua | InspectCutRunRegulatorySupport |
| tracks import-bed | glossary-command | MCP | ImportGenomeBedTrack |
| tracks import-bed | glossary-command | JS | ImportGenomeBedTrack |
| tracks import-bed | glossary-command | Lua | ImportGenomeBedTrack |
| tracks import-bigwig | glossary-command | MCP | ImportGenomeBigWigTrack |
| tracks import-bigwig | glossary-command | JS | ImportGenomeBigWigTrack |
| tracks import-bigwig | glossary-command | Lua | ImportGenomeBigWigTrack |
| tracks import-vcf | glossary-command | MCP | ImportGenomeVcfTrack |
| tracks import-vcf | glossary-command | JS | ImportGenomeVcfTrack |
| tracks import-vcf | glossary-command | Lua | ImportGenomeVcfTrack |
| arrays project-microarray-track | glossary-command | MCP | ProjectMicroarrayTrack |
| arrays project-microarray-track | glossary-command | JS | ProjectMicroarrayTrack |
| arrays project-microarray-track | glossary-command | Lua | ProjectMicroarrayTrack |
| guides put | glossary-command | MCP | UpsertGuideSet |
| guides put | glossary-command | JS | UpsertGuideSet |
| guides put | glossary-command | Lua | UpsertGuideSet |
| guides delete | glossary-command | MCP | DeleteGuideSet |
| guides delete | glossary-command | JS | DeleteGuideSet |
| guides delete | glossary-command | Lua | DeleteGuideSet |
| guides filter | glossary-command | MCP | FilterGuidesPractical |
| guides filter | glossary-command | JS | FilterGuidesPractical |
| guides filter | glossary-command | Lua | FilterGuidesPractical |
| guides oligos-generate | glossary-command | MCP | GenerateGuideOligos |
| guides oligos-generate | glossary-command | JS | GenerateGuideOligos |
| guides oligos-generate | glossary-command | Lua | GenerateGuideOligos |
| guides oligos-export | glossary-command | MCP | ExportGuideOligos |
| guides oligos-export | glossary-command | JS | ExportGuideOligos |
| guides oligos-export | glossary-command | Lua | ExportGuideOligos |
| guides protocol-export | glossary-command | MCP | ExportGuideProtocolText |
| guides protocol-export | glossary-command | JS | ExportGuideProtocolText |
| guides protocol-export | glossary-command | Lua | ExportGuideProtocolText |
| primers design | glossary-command | MCP | DesignPrimerPairs |
| primers design | glossary-command | JS | DesignPrimerPairs |
| primers design | glossary-command | Lua | DesignPrimerPairs |
| primers design-qpcr | glossary-command | MCP | DesignQpcrAssays |
| primers design-qpcr | glossary-command | JS | DesignQpcrAssays |
| primers design-qpcr | glossary-command | Lua | DesignQpcrAssays |
| primers test-cdna-pcr | glossary-command | MCP | TestCdnaPcr |
| primers test-cdna-pcr | glossary-command | JS | TestCdnaPcr |
| primers test-cdna-pcr | glossary-command | Lua | TestCdnaPcr |
| primers test-cdna-qpcr | glossary-command | MCP | TestCdnaQpcr |
| primers test-cdna-qpcr | glossary-command | JS | TestCdnaQpcr |
| primers test-cdna-qpcr | glossary-command | Lua | TestCdnaQpcr |
| primers test-cdna-qpcr-fasta | glossary-command | MCP | TestCdnaQpcrFasta |
| primers test-cdna-qpcr-fasta | glossary-command | JS | TestCdnaQpcrFasta |
| primers test-cdna-qpcr-fasta | glossary-command | Lua | TestCdnaQpcrFasta |
| primers prepare-restriction-cloning | glossary-command | MCP | PrepareRestrictionCloningPcrHandoff |
| primers prepare-restriction-cloning | glossary-command | JS | PrepareRestrictionCloningPcrHandoff |
| primers prepare-restriction-cloning | glossary-command | Lua | PrepareRestrictionCloningPcrHandoff |
| features repeat-query | glossary-command | MCP | QueryRepeatAnnotations |
| features repeat-query | glossary-command | JS | QueryRepeatAnnotations |
| features repeat-query | glossary-command | Lua | QueryRepeatAnnotations |
| features repeat-overlaps | glossary-command | MCP | QueryRepeatOverlaps |
| features repeat-overlaps | glossary-command | JS | QueryRepeatOverlaps |
| features repeat-overlaps | glossary-command | Lua | QueryRepeatOverlaps |
| features materialize-repeats | glossary-command | MCP | MaterializeRepeatFeatures |
| features materialize-repeats | glossary-command | JS | MaterializeRepeatFeatures |
| features materialize-repeats | glossary-command | Lua | MaterializeRepeatFeatures |
| features repeat-cohort | glossary-command | MCP | BuildRepeatEnvironmentCohort |
| features repeat-cohort | glossary-command | JS | BuildRepeatEnvironmentCohort |
| features repeat-cohort | glossary-command | Lua | BuildRepeatEnvironmentCohort |
| features window-cohort-tfbs | glossary-command | MCP | SummarizeWindowCohortTfbs |
| features window-cohort-tfbs | glossary-command | JS | SummarizeWindowCohortTfbs |
| features window-cohort-tfbs | glossary-command | Lua | SummarizeWindowCohortTfbs |
| features tfbs-score-tracks-svg | glossary-command | MCP | RenderTfbsScoreTracksSvg |
| features tfbs-score-tracks-svg | glossary-command | JS | RenderTfbsScoreTracksSvg |
| features tfbs-score-tracks-svg | glossary-command | Lua | RenderTfbsScoreTracksSvg |
| dotplot compute | glossary-command | MCP | ComputeDotplot |
| dotplot compute | glossary-command | JS | ComputeDotplot |
| dotplot compute | glossary-command | Lua | ComputeDotplot |
| dotplot overlay-compute | glossary-command | MCP | ComputeDotplotOverlay |
| dotplot overlay-compute | glossary-command | JS | ComputeDotplotOverlay |
| dotplot overlay-compute | glossary-command | Lua | ComputeDotplotOverlay |
| transcripts derive | glossary-command | MCP | DeriveTranscriptSequences |
| transcripts derive | glossary-command | JS | DeriveTranscriptSequences |
| transcripts derive | glossary-command | Lua | DeriveTranscriptSequences |
| transcripts exon-skip-plan | glossary-command | JS | PlanExonSkippedIsoform |
| transcripts exon-skip-plan | glossary-command | Lua | PlanExonSkippedIsoform |
| transcripts exon-skip-materialize | glossary-command | JS | MaterializeExonSkippedIsoform |
| transcripts exon-skip-materialize | glossary-command | Lua | MaterializeExonSkippedIsoform |
| flex compute | glossary-command | MCP | ComputeFlexibilityTrack |
| flex compute | glossary-command | JS | ComputeFlexibilityTrack |
| flex compute | glossary-command | Lua | ComputeFlexibilityTrack |
| splicing-refs derive | glossary-command | MCP | DeriveSplicingReferences |
| splicing-refs derive | glossary-command | JS | DeriveSplicingReferences |
| splicing-refs derive | glossary-command | Lua | DeriveSplicingReferences |
| align compute | glossary-command | MCP | AlignSequences |
| align compute | glossary-command | JS | AlignSequences |
| align compute | glossary-command | Lua | AlignSequences |
| seq-trace import | glossary-command | MCP | ImportSequencingTrace |
| seq-trace import | glossary-command | JS | ImportSequencingTrace |
| seq-trace import | glossary-command | Lua | ImportSequencingTrace |
| seq-trace list | glossary-command | MCP | ListSequencingTraces |
| seq-trace list | glossary-command | JS | ListSequencingTraces |
| seq-trace list | glossary-command | Lua | ListSequencingTraces |
| seq-trace show | glossary-command | MCP | ShowSequencingTrace |
| seq-trace show | glossary-command | JS | ShowSequencingTrace |
| seq-trace show | glossary-command | Lua | ShowSequencingTrace |
| seq-confirm run | glossary-command | MCP | ConfirmConstructReads |
| seq-confirm run | glossary-command | JS | ConfirmConstructReads |
| seq-confirm run | glossary-command | Lua | ConfirmConstructReads |
| seq-confirm list-reports | glossary-command | MCP | ListSequencingConfirmationReports |
| seq-confirm list-reports | glossary-command | JS | ListSequencingConfirmationReports |
| seq-confirm list-reports | glossary-command | Lua | ListSequencingConfirmationReports |
| seq-confirm show-report | glossary-command | MCP | ShowSequencingConfirmationReport |
| seq-confirm show-report | glossary-command | JS | ShowSequencingConfirmationReport |
| seq-confirm show-report | glossary-command | Lua | ShowSequencingConfirmationReport |
| seq-confirm export-report | glossary-command | MCP | ExportSequencingConfirmationReport |
| seq-confirm export-report | glossary-command | JS | ExportSequencingConfirmationReport |
| seq-confirm export-report | glossary-command | Lua | ExportSequencingConfirmationReport |
| seq-confirm export-support-tsv | glossary-command | MCP | ExportSequencingConfirmationSupportTsv |
| seq-confirm export-support-tsv | glossary-command | JS | ExportSequencingConfirmationSupportTsv |
| seq-confirm export-support-tsv | glossary-command | Lua | ExportSequencingConfirmationSupportTsv |
| seq-primer suggest | glossary-command | MCP | SuggestSequencingPrimers |
| seq-primer suggest | glossary-command | JS | SuggestSequencingPrimers |
| seq-primer suggest | glossary-command | Lua | SuggestSequencingPrimers |
| reverse-translate run | glossary-command | MCP | ReverseTranslateProteinSequence |
| reverse-translate run | glossary-command | JS | ReverseTranslateProteinSequence |
| reverse-translate run | glossary-command | Lua | ReverseTranslateProteinSequence |
| reverse-translate list-reports | glossary-command | MCP | ListReverseTranslationReports |
| reverse-translate list-reports | glossary-command | JS | ListReverseTranslationReports |
| reverse-translate list-reports | glossary-command | Lua | ListReverseTranslationReports |
| reverse-translate show-report | glossary-command | MCP | ShowReverseTranslationReport |
| reverse-translate show-report | glossary-command | JS | ShowReverseTranslationReport |
| reverse-translate show-report | glossary-command | Lua | ShowReverseTranslationReport |
| reverse-translate export-report | glossary-command | MCP | ExportReverseTranslationReport |
| reverse-translate export-report | glossary-command | JS | ExportReverseTranslationReport |
| reverse-translate export-report | glossary-command | Lua | ExportReverseTranslationReport |
| construct-reasoning build-protein-dna-handoff | glossary-command | MCP | BuildProteinToDnaHandoffReasoning |
| construct-reasoning build-protein-dna-handoff | glossary-command | JS | BuildProteinToDnaHandoffReasoning |
| construct-reasoning build-protein-dna-handoff | glossary-command | Lua | BuildProteinToDnaHandoffReasoning |
| construct-reasoning list-graphs | glossary-command | JS | ListConstructReasoningGraphs |
| construct-reasoning list-graphs | glossary-command | Lua | ListConstructReasoningGraphs |
| construct-reasoning show-graph | glossary-command | JS | ShowConstructReasoningGraph |
| construct-reasoning show-graph | glossary-command | Lua | ShowConstructReasoningGraph |
| construct-reasoning set-annotation-status | glossary-command | JS | SetConstructReasoningAnnotationStatus |
| construct-reasoning set-annotation-status | glossary-command | Lua | SetConstructReasoningAnnotationStatus |
| construct-reasoning write-annotation | glossary-command | JS | ConstructReasoningWriteAnnotation |
| construct-reasoning write-annotation | glossary-command | Lua | ConstructReasoningWriteAnnotation |
| construct-reasoning export-graph | glossary-command | MCP | ExportConstructReasoningGraph |
| construct-reasoning export-graph | glossary-command | JS | ExportConstructReasoningGraph |
| construct-reasoning export-graph | glossary-command | Lua | ExportConstructReasoningGraph |
| reads acquire status | glossary-command | MCP | ReadAcquireStatus |
| reads acquire status | glossary-command | JS | ReadAcquireStatus |
| reads acquire status | glossary-command | Lua | ReadAcquireStatus |
| reads acquire prepare | glossary-command | MCP | ReadAcquirePrepare |
| reads acquire prepare | glossary-command | JS | ReadAcquirePrepare |
| reads acquire prepare | glossary-command | Lua | ReadAcquirePrepare |
| reads acquire inspect | glossary-command | MCP | ReadAcquireInspect |
| reads acquire inspect | glossary-command | JS | ReadAcquireInspect |
| reads acquire inspect | glossary-command | Lua | ReadAcquireInspect |
| reads acquire cancel | glossary-command | MCP | ReadAcquireCancel |
| reads acquire cancel | glossary-command | JS | ReadAcquireCancel |
| reads acquire cancel | glossary-command | Lua | ReadAcquireCancel |
| rna-reads interpret | glossary-command | MCP | InterpretRnaReads |
| rna-reads interpret | glossary-command | JS | InterpretRnaReads |
| rna-reads interpret | glossary-command | Lua | InterpretRnaReads |
| rna-reads batch-map | glossary-command | MCP | RunRnaReadBatchMap |
| rna-reads batch-map | glossary-command | JS | RunRnaReadBatchMap |
| rna-reads batch-map | glossary-command | Lua | RunRnaReadBatchMap |
| rna-reads align-report | glossary-command | MCP | AlignRnaReadReport |
| rna-reads align-report | glossary-command | JS | AlignRnaReadReport |
| rna-reads align-report | glossary-command | Lua | AlignRnaReadReport |
| rna-reads preflight-isoforms | glossary-command | MCP | PreflightRnaReadIsoforms |
| rna-reads preflight-isoforms | glossary-command | JS | PreflightRnaReadIsoforms |
| rna-reads preflight-isoforms | glossary-command | Lua | PreflightRnaReadIsoforms |
| rna-reads summarize-gene-support | glossary-command | MCP | SummarizeRnaReadGeneSupport |
| rna-reads summarize-gene-support | glossary-command | JS | SummarizeRnaReadGeneSupport |
| rna-reads summarize-gene-support | glossary-command | Lua | SummarizeRnaReadGeneSupport |
| rna-reads inspect-gene-support | glossary-command | MCP | InspectRnaReadGeneSupport |
| rna-reads inspect-gene-support | glossary-command | JS | InspectRnaReadGeneSupport |
| rna-reads inspect-gene-support | glossary-command | Lua | InspectRnaReadGeneSupport |
| rna-reads materialize-hits | glossary-command | MCP | MaterializeRnaReadHitSequences |
| rna-reads materialize-hits | glossary-command | JS | MaterializeRnaReadHitSequences |
| rna-reads materialize-hits | glossary-command | Lua | MaterializeRnaReadHitSequences |
| rna-reads export-report | glossary-command | MCP | ExportRnaReadReport |
| rna-reads export-report | glossary-command | JS | ExportRnaReadReport |
| rna-reads export-report | glossary-command | Lua | ExportRnaReadReport |
| rna-reads export-hits-fasta | glossary-command | MCP | ExportRnaReadHitsFasta |
| rna-reads export-hits-fasta | glossary-command | JS | ExportRnaReadHitsFasta |
| rna-reads export-hits-fasta | glossary-command | Lua | ExportRnaReadHitsFasta |
| rna-reads export-sample-sheet | glossary-command | MCP | ExportRnaReadSampleSheet |
| rna-reads export-sample-sheet | glossary-command | JS | ExportRnaReadSampleSheet |
| rna-reads export-sample-sheet | glossary-command | Lua | ExportRnaReadSampleSheet |
| rna-reads export-target-quality | glossary-command | MCP | ExportRnaReadTargetQuality |
| rna-reads export-target-quality | glossary-command | JS | ExportRnaReadTargetQuality |
| rna-reads export-target-quality | glossary-command | Lua | ExportRnaReadTargetQuality |
| rna-reads export-paths-tsv | glossary-command | MCP | ExportRnaReadExonPathsTsv |
| rna-reads export-paths-tsv | glossary-command | JS | ExportRnaReadExonPathsTsv |
| rna-reads export-paths-tsv | glossary-command | Lua | ExportRnaReadExonPathsTsv |
| rna-reads export-abundance-tsv | glossary-command | MCP | ExportRnaReadExonAbundanceTsv |
| rna-reads export-abundance-tsv | glossary-command | JS | ExportRnaReadExonAbundanceTsv |
| rna-reads export-abundance-tsv | glossary-command | Lua | ExportRnaReadExonAbundanceTsv |
| rna-reads export-score-density-svg | glossary-command | MCP | ExportRnaReadScoreDensitySvg |
| rna-reads export-score-density-svg | glossary-command | JS | ExportRnaReadScoreDensitySvg |
| rna-reads export-score-density-svg | glossary-command | Lua | ExportRnaReadScoreDensitySvg |
| rna-reads export-alignments-tsv | glossary-command | MCP | ExportRnaReadAlignmentsTsv |
| rna-reads export-alignments-tsv | glossary-command | JS | ExportRnaReadAlignmentsTsv |
| rna-reads export-alignments-tsv | glossary-command | Lua | ExportRnaReadAlignmentsTsv |
| rna-reads export-isoform-triage-tsv | glossary-command | MCP | ExportRnaReadIsoformTriageTsv |
| rna-reads export-isoform-triage-tsv | glossary-command | JS | ExportRnaReadIsoformTriageTsv |
| rna-reads export-isoform-triage-tsv | glossary-command | Lua | ExportRnaReadIsoformTriageTsv |
| rna-reads export-alignment-dotplot-svg | glossary-command | MCP | ExportRnaReadAlignmentDotplotSvg |
| rna-reads export-alignment-dotplot-svg | glossary-command | JS | ExportRnaReadAlignmentDotplotSvg |
| rna-reads export-alignment-dotplot-svg | glossary-command | Lua | ExportRnaReadAlignmentDotplotSvg |
| candidates delete | glossary-command | MCP | DeleteCandidateSet |
| candidates delete | glossary-command | JS | DeleteCandidateSet |
| candidates delete | glossary-command | Lua | DeleteCandidateSet |
| candidates generate | glossary-command | MCP | GenerateCandidateSet |
| candidates generate | glossary-command | JS | GenerateCandidateSet |
| candidates generate | glossary-command | Lua | GenerateCandidateSet |
| candidates generate-between-anchors | glossary-command | MCP | GenerateCandidateSetBetweenAnchors |
| candidates generate-between-anchors | glossary-command | JS | GenerateCandidateSetBetweenAnchors |
| candidates generate-between-anchors | glossary-command | Lua | GenerateCandidateSetBetweenAnchors |
| candidates score | glossary-command | MCP | ScoreCandidateSetExpression |
| candidates score | glossary-command | JS | ScoreCandidateSetExpression |
| candidates score | glossary-command | Lua | ScoreCandidateSetExpression |
| candidates score-distance | glossary-command | MCP | ScoreCandidateSetDistance |
| candidates score-distance | glossary-command | JS | ScoreCandidateSetDistance |
| candidates score-distance | glossary-command | Lua | ScoreCandidateSetDistance |
| candidates score-weighted | glossary-command | MCP | ScoreCandidateSetWeightedObjective |
| candidates score-weighted | glossary-command | JS | ScoreCandidateSetWeightedObjective |
| candidates score-weighted | glossary-command | Lua | ScoreCandidateSetWeightedObjective |
| candidates top-k | glossary-command | MCP | TopKCandidateSet |
| candidates top-k | glossary-command | JS | TopKCandidateSet |
| candidates top-k | glossary-command | Lua | TopKCandidateSet |
| candidates pareto | glossary-command | MCP | ParetoFrontierCandidateSet |
| candidates pareto | glossary-command | JS | ParetoFrontierCandidateSet |
| candidates pareto | glossary-command | Lua | ParetoFrontierCandidateSet |
| candidates filter | glossary-command | MCP | FilterCandidateSet |
| candidates filter | glossary-command | JS | FilterCandidateSet |
| candidates filter | glossary-command | Lua | FilterCandidateSet |
| candidates set-op | glossary-command | MCP | CandidateSetOp |
| candidates set-op | glossary-command | JS | CandidateSetOp |
| candidates set-op | glossary-command | Lua | CandidateSetOp |
| macros template-put | glossary-command | MCP | UpsertWorkflowMacroTemplate |
| macros template-put | glossary-command | JS | UpsertWorkflowMacroTemplate |
| macros template-put | glossary-command | Lua | UpsertWorkflowMacroTemplate |
| macros template-delete | glossary-command | MCP | DeleteWorkflowMacroTemplate |
| macros template-delete | glossary-command | JS | DeleteWorkflowMacroTemplate |
| macros template-delete | glossary-command | Lua | DeleteWorkflowMacroTemplate |
| candidates macro | glossary-command | MCP | GenerateCandidateSet, GenerateCandidateSetBetweenAnchors, ScoreCandidateSetExpression, ScoreCandidateSetDistance, ScoreCandidateSetWeightedObjective, TopKCandidateSet, ParetoFrontierCandidateSet, FilterCandidateSet, CandidateSetOp, DeleteCandidateSet, UpsertCandidateMacroTemplate, DeleteCandidateMacroTemplate |
| candidates macro | glossary-command | JS | GenerateCandidateSet, GenerateCandidateSetBetweenAnchors, ScoreCandidateSetExpression, ScoreCandidateSetDistance, ScoreCandidateSetWeightedObjective, TopKCandidateSet, ParetoFrontierCandidateSet, FilterCandidateSet, CandidateSetOp, DeleteCandidateSet, UpsertCandidateMacroTemplate, DeleteCandidateMacroTemplate |
| candidates macro | glossary-command | Lua | GenerateCandidateSet, GenerateCandidateSetBetweenAnchors, ScoreCandidateSetExpression, ScoreCandidateSetDistance, ScoreCandidateSetWeightedObjective, TopKCandidateSet, ParetoFrontierCandidateSet, FilterCandidateSet, CandidateSetOp, DeleteCandidateSet, UpsertCandidateMacroTemplate, DeleteCandidateMacroTemplate |
| candidates template-put | glossary-command | MCP | UpsertCandidateMacroTemplate |
| candidates template-put | glossary-command | JS | UpsertCandidateMacroTemplate |
| candidates template-put | glossary-command | Lua | UpsertCandidateMacroTemplate |
| candidates template-delete | glossary-command | MCP | DeleteCandidateMacroTemplate |
| candidates template-delete | glossary-command | JS | DeleteCandidateMacroTemplate |
| candidates template-delete | glossary-command | Lua | DeleteCandidateMacroTemplate |
| set-param | glossary-command | MCP | SetParameter |
| set-param | glossary-command | JS | SetParameter |
| set-param | glossary-command | Lua | SetParameter |
| digest | glossary-command | GUI | Digest |
| digest | glossary-command | gentle_cli | Digest |
| digest | glossary-command | MCP | Digest |
| digest | glossary-command | Lua | Digest |
| export_dna_ladders | glossary-command | GUI | ExportDnaLadders |
| export_dna_ladders | glossary-command | gentle_cli | ExportDnaLadders |
| export_dna_ladders | glossary-command | MCP | ExportDnaLadders |
| export_dna_ladders | glossary-command | Lua | ExportDnaLadders |
| export_rna_ladders | glossary-command | GUI | ExportRnaLadders |
| export_rna_ladders | glossary-command | gentle_cli | ExportRnaLadders |
| export_rna_ladders | glossary-command | MCP | ExportRnaLadders |
| export_rna_ladders | glossary-command | Lua | ExportRnaLadders |
| extend_genome_anchor | glossary-command | GUI | ExtendGenomeAnchor |
| extend_genome_anchor | glossary-command | gentle_cli | ExtendGenomeAnchor |
| extend_genome_anchor | glossary-command | MCP | ExtendGenomeAnchor |
| extend_genome_anchor | glossary-command | Lua | ExtendGenomeAnchor |
| extract_genome_gene | glossary-command | GUI | ExtractGenomeGene |
| extract_genome_gene | glossary-command | gentle_cli | ExtractGenomeGene |
| extract_genome_gene | glossary-command | MCP | ExtractGenomeGene |
| extract_genome_gene | glossary-command | Lua | ExtractGenomeGene |
| extract_genome_region | glossary-command | GUI | ExtractGenomeRegion |
| extract_genome_region | glossary-command | gentle_cli | ExtractGenomeRegion |
| extract_genome_region | glossary-command | MCP | ExtractGenomeRegion |
| extract_genome_region | glossary-command | Lua | ExtractGenomeRegion |
| import_genome_bed_track | glossary-command | GUI | ImportGenomeBedTrack |
| import_genome_bed_track | glossary-command | gentle_cli | ImportGenomeBedTrack |
| import_genome_bed_track | glossary-command | MCP | ImportGenomeBedTrack |
| import_genome_bed_track | glossary-command | Lua | ImportGenomeBedTrack |
| import_genome_bigwig_track | glossary-command | GUI | ImportGenomeBigWigTrack |
| import_genome_bigwig_track | glossary-command | gentle_cli | ImportGenomeBigWigTrack |
| import_genome_bigwig_track | glossary-command | MCP | ImportGenomeBigWigTrack |
| import_genome_bigwig_track | glossary-command | Lua | ImportGenomeBigWigTrack |
| import_genome_vcf_track | glossary-command | GUI | ImportGenomeVcfTrack |
| import_genome_vcf_track | glossary-command | gentle_cli | ImportGenomeVcfTrack |
| import_genome_vcf_track | glossary-command | MCP | ImportGenomeVcfTrack |
| import_genome_vcf_track | glossary-command | Lua | ImportGenomeVcfTrack |
| prepare_genome | glossary-command | GUI | PrepareGenome |
| prepare_genome | glossary-command | gentle_cli | PrepareGenome |
| prepare_genome | glossary-command | MCP | PrepareGenome |
| prepare_genome | glossary-command | Lua | PrepareGenome |
| render_pool_gel_svg | glossary-command | GUI | RenderPoolGelSvg |
| render_pool_gel_svg | glossary-command | gentle_cli | RenderPoolGelSvg |
| render_pool_gel_svg | glossary-command | MCP | RenderPoolGelSvg |
| render_pool_gel_svg | glossary-command | Lua | RenderPoolGelSvg |
| set_parameter | glossary-command | GUI | SetParameter |
| set_parameter | glossary-command | gentle_cli | SetParameter |
| set_parameter | glossary-command | MCP | SetParameter |
| set_parameter | glossary-command | Lua | SetParameter |
| set_vcf_display_filter | glossary-command | GUI | SetParameter |
| set_vcf_display_filter | glossary-command | gentle_cli | SetParameter |
| set_vcf_display_filter | glossary-command | MCP | SetParameter |
| set_vcf_display_filter | glossary-command | Lua | SetParameter |
| export_dna_ladders | glossary-command | GUI | ExportDnaLadders |
| export_dna_ladders | glossary-command | gentle_cli | ExportDnaLadders |
| export_dna_ladders | glossary-command | MCP | ExportDnaLadders |
| export_dna_ladders | glossary-command | JS | ExportDnaLadders |
| export_rna_ladders | glossary-command | GUI | ExportRnaLadders |
| export_rna_ladders | glossary-command | gentle_cli | ExportRnaLadders |
| export_rna_ladders | glossary-command | MCP | ExportRnaLadders |
| export_rna_ladders | glossary-command | JS | ExportRnaLadders |
| extend_genome_anchor | glossary-command | GUI | ExtendGenomeAnchor |
| extend_genome_anchor | glossary-command | gentle_cli | ExtendGenomeAnchor |
| extend_genome_anchor | glossary-command | MCP | ExtendGenomeAnchor |
| extend_genome_anchor | glossary-command | JS | ExtendGenomeAnchor |
| extract_genome_gene | glossary-command | GUI | ExtractGenomeGene |
| extract_genome_gene | glossary-command | gentle_cli | ExtractGenomeGene |
| extract_genome_gene | glossary-command | MCP | ExtractGenomeGene |
| extract_genome_gene | glossary-command | JS | ExtractGenomeGene |
| extract_genome_region | glossary-command | GUI | ExtractGenomeRegion |
| extract_genome_region | glossary-command | gentle_cli | ExtractGenomeRegion |
| extract_genome_region | glossary-command | MCP | ExtractGenomeRegion |
| extract_genome_region | glossary-command | JS | ExtractGenomeRegion |
| import_genome_bed_track | glossary-command | GUI | ImportGenomeBedTrack |
| import_genome_bed_track | glossary-command | gentle_cli | ImportGenomeBedTrack |
| import_genome_bed_track | glossary-command | MCP | ImportGenomeBedTrack |
| import_genome_bed_track | glossary-command | JS | ImportGenomeBedTrack |
| import_genome_bigwig_track | glossary-command | GUI | ImportGenomeBigWigTrack |
| import_genome_bigwig_track | glossary-command | gentle_cli | ImportGenomeBigWigTrack |
| import_genome_bigwig_track | glossary-command | MCP | ImportGenomeBigWigTrack |
| import_genome_bigwig_track | glossary-command | JS | ImportGenomeBigWigTrack |
| import_genome_vcf_track | glossary-command | GUI | ImportGenomeVcfTrack |
| import_genome_vcf_track | glossary-command | gentle_cli | ImportGenomeVcfTrack |
| import_genome_vcf_track | glossary-command | MCP | ImportGenomeVcfTrack |
| import_genome_vcf_track | glossary-command | JS | ImportGenomeVcfTrack |
| prepare_genome | glossary-command | GUI | PrepareGenome |
| prepare_genome | glossary-command | gentle_cli | PrepareGenome |
| prepare_genome | glossary-command | MCP | PrepareGenome |
| prepare_genome | glossary-command | JS | PrepareGenome |
| render_pool_gel_svg | glossary-command | GUI | RenderPoolGelSvg |
| render_pool_gel_svg | glossary-command | gentle_cli | RenderPoolGelSvg |
| render_pool_gel_svg | glossary-command | MCP | RenderPoolGelSvg |
| render_pool_gel_svg | glossary-command | JS | RenderPoolGelSvg |
| set_parameter | glossary-command | GUI | SetParameter |
| set_parameter | glossary-command | gentle_cli | SetParameter |
| set_parameter | glossary-command | MCP | SetParameter |
| set_parameter | glossary-command | JS | SetParameter |
| set_vcf_display_filter | glossary-command | GUI | SetParameter |
| set_vcf_display_filter | glossary-command | gentle_cli | SetParameter |
| set_vcf_display_filter | glossary-command | MCP | SetParameter |
| set_vcf_display_filter | glossary-command | JS | SetParameter |
