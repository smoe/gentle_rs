# MCP Server

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

## `gentle_mcp` (MCP stdio server)

`gentle_mcp` starts a Model Context Protocol server over stdio.
It provides both tool execution (`tools/call`) and capability
discovery/negotiation (`tools/list`, `capabilities`, `help`).

Current tools:

- `capabilities`
- `state_summary`
- `runtime_status`
- `restriction_site_detail` (shared restriction-site expert detail record,
  including `tooltip_lines[]`)
- `exon_skip_plan` (shared `transcripts exon-skip-plan` contract)
- `exon_skip_materialize` (shared `transcripts exon-skip-materialize`
  contract; requires `confirm=true` and accepts requested return payloads)
- `promoter_reporter_panel_plan` (read-only shared `promoters panel-plan`
  contract; accepts the exact request object)
- `promoter_reporter_panel_readiness` (read-only shared
  `promoters panel-readiness` contract; accepts a request or proposal and an
  optional exact approval digest)
- `promoter_reporter_panel_materialize` (shared `promoters panel-materialize`
  contract; requires `confirm=true` and the exact approved proposal digest)
- `op` (apply one operation; requires `confirm=true`)
- `workflow` (apply one workflow; requires `confirm=true`)
- `help`
- `agent_systems` (shared `agents list` contract)
- `agent_preflight` (shared `agents preflight` contract; optional
  `live=true` maps to `--live`)
- `agent_models` (shared `agents discover-models` contract)
- `reference_catalog_entries` (structured reference catalog rows via the shared `genomes list` contract)
- `helper_catalog_entries` (structured helper catalog rows with normalized helper `interpretation` when available)
- `host_profile_catalog_entries` (structured host-profile catalog rows via the shared `hosts list` contract)
- `ensembl_installable_genomes` (shared Ensembl discovery report for currently installable candidates)
- `helper_interpretation` (direct helper-construct interpretation lookup by id or alias)
- `ui_intents` (discover deterministic UI-intent contracts)
- `ui_intent` (run deterministic `ui open|focus|close` intent resolution path)
- `ui_prepared_genomes` (run deterministic prepared-genome query path)
- `ui_latest_prepared` (resolve latest prepared genome for one species)
- `blast_async_start` (start async BLAST job through shared shell route)
- `blast_async_status` (poll async BLAST job status and optional report)
- `blast_async_cancel` (request cancellation for async BLAST job)
- `blast_async_list` (list known async BLAST jobs in scope)

Tool-parity rule:

- MCP UI-intent tools execute through the same shared parser/executor route as
  CLI shell commands (`ui ...`), and return the same structured payload
  contracts.
- No MCP-only biology/UI logic branches are allowed.

Data-boundary and trust model:

- `gentle_mcp` reads/writes project state via the resolved `state_path`
  (default `.gentle_state.json`).
- It does not attach to another running GUI process's in-memory engine state.
- Unsaved GUI-only edits are therefore not visible to MCP until they are
  persisted to disk and read from that same state path.
- If MCP is pointed at the same project state file used by active GUI/CLI
  workflows, that shared visibility is intentional.
- Treat connected MCP clients as trusted local operators; `confirm=true`
  expresses explicit mutating intent but is not an authentication boundary.

Run:

```bash
cargo run --bin gentle_mcp
cargo run --bin gentle_mcp -- --state path/to/project.gentle.json
```

External-agent handoff:

- Use `gentle_mcp` when an outside agent already supports MCP and already has
  its own subscription/runtime.
- A ChatGPT/Codex subscription is not an OpenAI API key for GENtle's native
  OpenAI API transport.
- The MCP route exposes local GENtle tools over stdio and does not make GENtle
  spend OpenAI API usage by itself.

Minimum MCP JSON-RPC flow:

1. Initialize:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}
```

2. List tools:

```json
{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}
```

3. Call one tool:

```json
{
  "jsonrpc":"2.0",
  "id":3,
  "method":"tools/call",
  "params":{
    "name":"ui_intent",
    "arguments":{
      "action":"open",
      "target":"prepared-references",
      "catalog_path":"assets/genomes.json",
      "species":"human",
      "latest":true
    }
  }
}
```

`ui_intent` arguments:

- required:
  - `action`: `open|focus|close`
  - `target`: one of:
    - `sequence-window` (requires `seq_id`; controls a loaded DNA sequence
      viewer without mutating the sequence record)
    - `prepared-references`
    - `prepare-reference-genome`
    - `retrieve-genome-sequence`
    - `blast-genome-sequence`
    - `import-genome-track`
    - `pcr-design`
    - `sequencing-confirmation`
    - `protein-expression-handoff` (no arguments; read-only report inspector)
    - `agent-assistant`
    - `prepare-helper-genome`
    - `retrieve-helper-sequence`
    - `blast-helper-sequence`
- optional:
  - `state_path`
  - `seq_id`
  - `genome_id`
  - `helpers`
  - `catalog_path`
  - `cache_dir`
  - `filter`
  - `species`
  - `latest`

`ui_prepared_genomes` arguments:

- optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`, `filter`,
  `species`, `latest`

`ui_latest_prepared` arguments:

- required: `species`
- optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`

`reference_catalog_entries` arguments:

- optional: `catalog_path`, `filter`

`helper_catalog_entries` arguments:

- optional: `catalog_path`, `filter`

`helper_semantics_vocabulary` arguments:

- optional: `vocabulary_path`, `filter`

`host_profile_catalog_entries` arguments:

- optional: `catalog_path`, `filter`

`ensembl_installable_genomes` arguments:

- optional: `collection`, `filter`

`helper_interpretation` arguments:

- required: `helper_id` (catalog id or alias)
- optional: `catalog_path`

`restriction_site_detail` arguments:

- required: `seq_id`, `cut_pos_1based`
- optional: `state_path`, `enzyme`, `recognition_start_1based`,
  `recognition_end_1based`
- returns the same `kind = "restriction_site"` payload as
  `inspect-feature-expert`, including `data.tooltip_lines[]`

`blast_async_start` arguments:

- required: `genome_id`, `query_sequence`
- optional: `state_path`, `helpers`, `max_hits`, `task`, `options_json`,
  `catalog_path`, `cache_dir`

`blast_async_status` arguments:

- required: `job_id`
- optional: `helpers`, `with_report`

`blast_async_cancel` arguments:

- required: `job_id`
- optional: `helpers`

`blast_async_list` arguments:

- optional: `helpers`

`ui_intents` arguments:

- optional: `state_path` (accepted for API symmetry)

`runtime_status` arguments:

- optional `state_path` (used to inspect project-backed async registries such
  as BLAST jobs)
- returns the same `gentle.runtime_status.v1` shape as `introspect runtime`,
  tagged for the MCP process, with live `frames[]` plus observed
  `activities[]`

MCP tool result envelope behavior:

- `result.isError = false`: success; inspect `result.structuredContent`.
- `result.isError = true`: deterministic tool-level error; inspect
  `result.content[0].text`.
- UI-intent structured payloads match shared shell schemas:
  - `gentle.ui_intents.v1`
  - `gentle.ui_intent.v1`
  - `gentle.ui_recent_project_intent.v1`
  - `gentle.ui_tutorial_project_intent.v1`
  - `gentle.ui_configuration_intent.v1`
  - `gentle.ui_prepared_genomes.v1`
  - `gentle.ui_latest_prepared.v1`
- `gentle.ui_intents.v1` now includes per-target `target_details[]` rows with
  stable `title`, `detail`, `keywords`, `menu_path`, supported `actions`, and
  `optional_arguments` metadata so GUI command-palette and MCP discovery stay
  aligned with the shared shell target catalog.
- catalog-entry MCP tools match the shared shell list payload shape:
  - `catalog_path`
  - `filter`
  - `genome_count`
  - `genomes[]`
  - `entries[]`
- `helper_interpretation` returns:
  - `query`
  - `catalog_path`
  - `interpretation` (`null` when the helper exists but has no structured semantics)

Mutating tool safety:

- `op` and `workflow` reject execution unless `confirm=true` is explicitly set.
- UI-intent tools are currently non-mutating; if a routed command ever reports
  state mutation unexpectedly, MCP returns an explicit tool error.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Python Wrapper](python.md)
<!-- manual-page-nav:end -->
