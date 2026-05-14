# GUI / CLI / MCP Parity Audit

Snapshot: current local `gentle_rs_3`/`main` at `b421022b`.

This is a read-only action-surface audit. It does not close any parity gaps.

## Method

Commands and sources inspected:

- GUI action callsites:
  `rg -n "(\\.button\\(|\\.small_button\\(|\\.menu_button\\()" src --glob '*.rs'`
- Shared shell and CLI command catalog:
  `docs/glossary.json`, `src/engine_shell.rs`, `crates/gentle-shell/`, and
  `src/bin/gentle_cli.rs`
- MCP tool registry:
  `src/mcp_server.rs`, especially `tool_list()`
- JavaScript adapter:
  `src/js_interface.rs`
- Lua adapter:
  `src/lua_interface.rs`

Raw counts from this snapshot:

- GUI button/menu callsites: 617
- glossary commands: 409
- MCP tools from `tools/list`: 31
- JavaScript wrappers: 49
- Lua wrappers: 48

The GUI inventory is larger than 200 rows, so this matrix groups repeated
buttons and local controls into action families. A `yes` in the GUI column
means at least one visible GUI button, small button, or menu button exists for
that action family. A `yes` in MCP/JS/Lua may mean a direct named tool/wrapper,
or the documented generic `op`/`workflow`/`apply_operation` path when the action
is an engine-owned operation. The Notes column names that distinction.

`n/a` means the action is local GUI chrome, local adapter bootstrap, or another
surface-specific concern rather than a shared engine/shell contract.

## Findings

- No action family in this audit appears to be `GUI=yes` and `gentle_cli=no`
  for a non-inspection engine or workflow action.
- Most newer biology/workflow routes have CLI/shared-shell coverage, but MCP is
  still intentionally narrower: it provides a small typed tool set plus generic
  `op` and `workflow` mutation routes.
- JavaScript and Lua expose a thin curated wrapper set plus generic
  `apply_operation` / `apply_workflow`; they do not mirror every shell command.
- `docs/glossary.json` now includes the implemented agent setup/planning routes
  and construct-reasoning graph routes surfaced by this audit. The
  protein-to-DNA handoff build route remains CLI/shared-shell first for MCP
  callers, which can still use the generic `op` tool with
  `BuildProteinToDnaHandoffReasoning`.

## App And Project Shell

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| Project load/save | yes | yes | no | yes | yes | CLI has `load-project`/`save-project`; MCP works from state paths but has no explicit load/save tool. |
| Sequence import/export files | yes | yes | yes | yes | yes | CLI covers project/sequence/pool import paths; MCP/JS/Lua can use generic operation paths where the action is engine-owned. |
| Sequence SVG export | yes | yes | yes | yes | yes | CLI `render-svg`; MCP via `op`; JS/Lua via `apply_operation`. |
| DALG/lineage SVG export | yes | yes | yes | yes | yes | CLI `render-lineage-svg`; generic operation adapters can apply `RenderLineageSvg`. |
| Help, command glossary, capabilities | yes | yes | yes | yes | yes | MCP has `help` and `capabilities`; JS/Lua expose `capabilities`. |
| State summary | no | yes | yes | yes | yes | Not a visible GUI action; GUI renders equivalent state internally. |
| Operation history undo/redo/status | yes | yes | no | no | no | CLI/shared shell has `history status/undo/redo`; MCP intentionally does not expose session-local undo. |
| Command palette | yes | n/a | n/a | n/a | n/a | GUI navigation surface only. |
| Window open/focus/menu routing | yes | yes | yes | no | no | CLI has `ui intents`, `ui open`, `ui focus`; MCP has `ui_intents`/`ui_intent`. |
| GUI configuration/theme/window styling | yes | n/a | n/a | n/a | n/a | Local application preference surface, not shared engine state. |
| Disabled screenshot command | no | yes | no | no | no | CLI glossary reserves `screenshot-window`; GUI screenshot capture remains policy-disabled. |

## Agents, Planning, Routines, And Services

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| Agent-system list/config overview | yes | yes | yes | yes | yes | MCP `agent_systems`; JS/Lua `list_agent_systems`. |
| Agent setup preflight / Test Setup | yes | yes | yes | no | no | CLI route exists as `agents preflight`; MCP `agent_preflight`. |
| Agent model discovery | yes | yes | yes | no | no | CLI route exists as `agents discover-models`; MCP `agent_models`. |
| Agent ask | yes | yes | no | yes | yes | CLI `agents ask`; JS/Lua expose `ask_agent_system`. MCP exposes plan/execute, not `ask`. |
| Agent plan and execute-plan | yes | yes | yes | yes | yes | CLI routes exist; MCP `agent_plan`/`agent_execute_plan`; JS/Lua wrappers exist. |
| External-agent MCP handoff instructions | yes | n/a | n/a | n/a | n/a | GUI documentation/snippet around running `gentle_mcp`, not a biology operation. |
| Routine catalog list/explain/compare | yes | yes | no | no | no | CLI `routines ...`; GUI Routine Assistant wraps the same catalog concepts. |
| Routine macro preflight/run | yes | yes | no | no | no | CLI `macros template-run`; GUI Routine Assistant shells out for validation/run. |
| Planning profile/objective/suggestions/sync | yes | yes | no | no | no | CLI `planning ...`; no typed MCP/JS/Lua route beyond generic state mutation. |
| External service status/providers/doctor | yes | yes | no | no | no | CLI `services status`, `services providers ...`; GUI External Services panel. |
| External service quote/handoff/guide | yes | yes | no | no | no | CLI `services project-preflight`, `project-quote`, `handoff`, `guide`. |

## Reference Data, Resources, And Catalogs

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| REBASE/JASPAR resource sync | yes | yes | no | yes | yes | JS/Lua expose `sync_rebase` and `sync_jaspar`; MCP has no direct resource-sync tool. |
| UCSC RepeatMasker sync/index | yes | yes | no | no | no | CLI `resources sync-ucsc-rmsk`, `install-ucsc-rmsk`, index helpers. |
| Resource status/benchmark/list/inspect | yes | yes | no | no | no | GUI resource dialogs and JASPAR Expert; CLI `resources ...`. |
| JASPAR Expert list/inspect visible motifs | yes | yes | no | no | no | CLI `resources list-jaspar` and `resources inspect-jaspar`. |
| Gene-group list/show/resolve/doctor/draft | no | yes | no | no | no | Headless/catalog route only in this snapshot. |
| Reference genome catalog list/validate/status | yes | yes | yes | yes | yes | MCP has catalog-entry and prepared-genome tools; JS/Lua expose reference list/status helpers. |
| Helper catalog/vocabulary/status | yes | yes | yes | yes | yes | MCP has helper catalog, vocabulary, and interpretation tools. |
| Host-profile catalog | no | yes | yes | yes | yes | Primarily agent/construct-reasoning support data. |
| Ensembl installable genome discovery | yes | yes | yes | yes | yes | MCP `ensembl_installable_genomes`; JS/Lua `list_ensembl_installable_genomes`. |
| Ensembl prepare prepared refs | yes | yes | yes | yes | yes | `PrepareGenome` is operation-backed and available through generic adapters. |
| Ensembl install/update/remove catalog entries | yes | yes | no | no | no | CLI shell routes manage catalog/spec files; no typed MCP/JS/Lua route. |
| Genome/helper retrieve region/gene/promoter | yes | yes | yes | yes | yes | CLI `genomes/helpers extract-*`; JS/Lua have reference extraction wrappers. |
| Genome/helper extend/verify anchor | yes | yes | yes | yes | yes | JS/Lua expose `extend_genome_anchor`; MCP via generic `op`. |
| Genome/helper BLAST start/status/cancel/list | yes | yes | yes | no | no | MCP has direct async BLAST tools; JS/Lua expose direct synchronous BLAST wrappers only. |
| Genome/helper direct BLAST | yes | yes | no | yes | yes | JS/Lua `blast_reference_genome`/`blast_helper_genome`; MCP uses the async tool family instead. |
| BED/BigWig/VCF track import | yes | yes | yes | yes | yes | JS/Lua expose track import wrappers; MCP via generic `op`. |
| Cache inspect/clear | yes | yes | no | no | no | CLI `cache inspect/clear`; GUI clear-cache panels. |

## DNA Viewer, Features, And Visualization

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| Reverse/complement/revcomp/branch | yes | yes | yes | yes | yes | Engine operations through CLI `op`, MCP `op`, JS/Lua `apply_operation`. |
| Selection extraction and sequence export | yes | yes | yes | yes | yes | Engine-owned operation/export paths. |
| Restriction digest/merge/MW filtering | yes | yes | yes | yes | yes | CLI generic operations plus JS `digest`; other mutations via generic adapters. |
| Restriction-site detail inspection | yes | yes | yes | yes | yes | CLI `inspect-feature-expert`; MCP direct `restriction_site_detail`; JS/Lua feature expert wrappers. |
| Feature expert inspection/SVG | yes | yes | yes | yes | yes | CLI `inspect-feature-expert` and `render-feature-expert-svg`; MCP via specific tools/generic `op`. |
| RE scan / TFBS scan / TFBS score tracks | yes | yes | yes | yes | yes | Shared engine reports, CLI `features ...`, generic adapter operations. |
| TFBS similarity and repeat/materialize views | yes | yes | yes | yes | yes | CLI `features ...`; generic operation path for mutations. |
| Array/microarray track inspection/projection | yes | yes | yes | no | no | CLI `arrays ...`; MCP can use generic `op` where applicable. |
| CUT&RUN report/coverage/support inspection | yes | yes | yes | no | no | CLI `cutrun ...`; MCP generic `op` where operation-backed. |
| Dotplot compute/show/render | yes | yes | yes | yes | yes | CLI `dotplot ...` plus `render-dotplot-svg`; JS/Lua `render_dotplot_svg`; MCP via `op`. |
| Flexibility compute/list/show | yes | yes | yes | yes | yes | Engine operation path through generic adapters. |
| Pairwise alignment compute | yes | yes | yes | yes | yes | CLI `align compute`; generic operation path. |
| Transcript derivation and exon-skip planning | yes | yes | yes | yes | yes | MCP has direct `exon_skip_plan` and `exon_skip_materialize`. |
| Splicing reference derivation | yes | yes | yes | yes | yes | CLI `splicing-refs derive`; generic operation path. |
| RNA rendering/info | yes | yes | yes | yes | yes | CLI `render-rna-svg`/`rna-info`; generic operation path where applicable. |
| Protocol cartoon list/render/template actions | yes | yes | yes | yes | yes | Render/validate/export are operation-backed through generic adapters; list has no typed MCP/JS/Lua wrapper. |

## Cloning, Racks, Pools, And Ladders

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| Gibson preview | yes | yes | no | no | no | CLI `gibson preview`; GUI specialist has a richer preview surface. |
| Gibson apply | yes | yes | yes | yes | yes | CLI `gibson apply`; generic operation routes can apply persisted operation payloads. |
| DNA/RNA ladder list/export | no | yes | no | yes | yes | Adapter helper surface; not a visible GUI action family in this audit. |
| Pool import/export/run bundle | yes | yes | yes | yes | yes | JS/Lua `import_pool`; MCP via generic `op` for operation-backed mutations. |
| Pool gel SVG rendering | yes | yes | yes | yes | yes | CLI `render-pool-gel-svg`; JS/Lua `render_pool_gel_svg`; MCP via generic `op`. |
| Serial arrangements | yes | yes | yes | no | no | CLI `arrange-serial` and `arrange-set-ladders`; MCP via generic operation where available. |
| Rack create/place/move/profile | yes | yes | yes | no | no | CLI `racks ...`; MCP via generic operation where available. |
| Rack labels/fabrication/isometric/OpenSCAD exports | yes | yes | yes | yes | yes | CLI `racks labels-svg`, `fabrication-svg`, `isometric-svg`, `openscad`; generic operation adapters can apply the export operations. |

## Primers, qPCR, RNA Reads, And Sequencing

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| Primer-pair design | yes | yes | yes | yes | yes | CLI `primers design`; generic operation adapters. |
| qPCR assay design | yes | yes | yes | yes | yes | CLI `primers design-qpcr`; generic operation adapters. |
| cDNA PCR/qPCR testing on annotated transcripts | yes | yes | yes | yes | yes | CLI `primers test-cdna-pcr` and `test-cdna-qpcr`; underlying operations are generic-adapter reachable. |
| cDNA/ncRNA FASTA qPCR screen | yes | yes | yes | yes | yes | CLI `primers test-cdna-qpcr-fasta`; underlying operation is generic-adapter reachable. |
| Primer report list/show/export | yes | yes | no | no | no | CLI `primers list/show/export-*`; report-store shell routes. |
| Primer3 backend preflight | yes | yes | no | no | no | CLI `primers preflight`; GUI Probe Primer3. |
| Restriction-cloning handoff | yes | yes | yes | yes | yes | CLI `primers prepare/seed/list/show/export-restriction-cloning-*`; generic operation path for engine-owned pieces. |
| RNA-read interpretation | yes | yes | yes | yes | yes | CLI `rna-reads interpret`; generic operation path. |
| RNA-read phase-2 alignment/report refresh | yes | yes | yes | yes | yes | CLI `rna-reads align-report`; generic operation path. |
| RNA-read report/inspection views | yes | yes | no | no | no | CLI `rna-reads list/show/inspect-*`; mostly shell/report inspection. |
| RNA-read concatemer/gene-support review | yes | yes | no | no | no | CLI `inspect-concatemers`, `inspect-gene-support`, `summarize-gene-support`. |
| RNA-read materialize hits | yes | yes | yes | yes | yes | CLI `rna-reads materialize-hits`; operation-backed. |
| RNA-read TSV/FASTA/SVG/sample-sheet exports | yes | yes | yes | yes | yes | CLI `rna-reads export-*`; export operations are generic-adapter reachable. |
| Read acquisition status/prepare/inspect/cancel | yes | yes | no | no | no | CLI `reads acquire ...`; GUI task/cache controls. |
| Sequencing trace import/list/show | yes | yes | yes | yes | yes | CLI `seq-trace ...`; import is operation-backed. |
| Sequencing confirmation run/reports/exports | yes | yes | yes | yes | yes | CLI `seq-confirm ...`; export subroutes are shell-only but core run is operation-backed. |
| Sequencing primer suggestion | yes | yes | yes | yes | yes | CLI `seq-primer suggest`; generic operation path. |

## Protein, Projection, And Construct Reasoning

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| GenBank/dbSNP fetch | yes | yes | no | no | no | CLI `genbank fetch`, `dbsnp fetch`; GUI ingress dialogs. |
| UniProt fetch/import/list/show/map | yes | yes | no | no | no | CLI `uniprot ...`; no typed MCP/JS/Lua route. |
| UniProt projection/accounting/audits | yes | yes | no | no | no | CLI `uniprot projection-*`, `audit-*`, comparison commands. |
| Ensembl gene/region/protein fetch/import | yes | yes | no | no | no | CLI `ensembl-*`; no typed MCP/JS/Lua route beyond genome catalog helpers. |
| Protein Evidence expert open/inspect | yes | yes | yes | yes | yes | CLI feature/protein expert routes; generic operation/feature expert adapters. |
| Derived/protein mapping SVG exports | yes | yes | yes | yes | yes | CLI render/projection routes; JS/Lua `render_feature_expert_svg` where feature expert based. |
| Reverse translation | yes | yes | yes | yes | yes | CLI `reverse-translate ...`; generic operation path. |
| Protease list/show/digest/gel SVG | yes | yes | yes | yes | yes | CLI `proteases ...`; digest is operation-backed. |
| Construct-reasoning graph list/show/status/write | yes | yes | yes | yes | yes | MCP has direct list/show/status/write tools; JS/Lua wrappers exist. |
| Protein-to-DNA handoff reasoning | yes | yes | yes | yes | yes | CLI route exists and GUI button exists; MCP uses generic `op`, with no direct typed build tool in this snapshot. |
| Candidate-set and guide-set workflows | no | yes | no | no | no | Headless shell surface in this snapshot. GUI Routine Assistant may suggest macros but does not expose the full command family. |

## Adapter-Only Primitives

| Action family | GUI | gentle_cli | MCP | JS | Lua | Notes |
|---|---|---|---|---|---|---|
| Generic operation application | no | yes | yes | yes | yes | CLI `op`; MCP `op`; JS/Lua `apply_operation`. |
| Generic workflow application | no | yes | yes | yes | yes | CLI `workflow`; MCP `workflow`; JS/Lua `apply_workflow`. |
| JavaScript/Lua DNA load and GenBank write helpers | no | no | no | yes | yes | Adapter convenience functions, not GUI or CLI commands. |
| JavaScript/Lua VCF display filter helper | no | no | no | yes | yes | Adapter convenience wrapper around parameter operations. |

## Roadmap Classification

The strict parity-bug rule for this audit was:

> `GUI=yes` and `gentle_cli=no` for a non-inspection action.

No such row was found at the action-family level. The glossary parity gaps that
were found during this audit have since been closed, so this snapshot does not
leave a roadmap parity TODO.
