# GENtle Agent Interface

Last updated: 2026-06-21

This guide explains how agents can control GENtle and how the available
interfaces differ.

For a self-standing narrative tutorial with role-based usage, step-by-step
flows, and explicit comparisons (CLI vs MCP vs Agent Assistant vs external
coding agents like Codex), start with:

- `docs/tutorial/01-01_agent_interfaces.md`

For the shortest working-tree loop an agent can run while editing GENtle, use:

- `docs/agent_dev_loop.md`

## Why this exists

GENtle exposes multiple agent-facing routes:

- deterministic command execution through CLI/shared shell
- machine-facing prose compilation through the typed planner
- deterministic tool execution through MCP
- guided assistant interaction through the Agent Assistant bridge

These routes are related, but they are not identical in transport, ergonomics,
and safety flow.

Plain-language note used in this file:

- "fixed format" means a command/tool with defined inputs and predictable
  outputs (usually JSON).

## First run from an empty project

When GENtle opens with an empty project, the inner Agent Assistant should be
treated as a command-suggestion layer, not as a database client that already
knows the current project. Start with a small, observable loop:

1. Open `File -> Agent Assistant...`.
2. Choose the provider profile (`Local Model`, `Codex Local`, OpenAI, Claude,
   Mistral, or another catalog entry).
3. For local OpenAI-compatible services such as Ollama, Jan, or Msty, set the
   base URL and click `Discover Models`; then pick one concrete discovered
   model.
4. Click `Test Setup`. This checks endpoint/model reachability without sending
   a generation request.
5. Leave `Project summary` / `Include state summary` unchecked for the first
   prompt. In an empty project there is no useful project context to send, and
   small local models often behave better with the shortest possible request.
6. Ask a response-format probe before asking biology:

```text
Introduce yourself briefly as GENtle's internal Agent Assistant.

Return strict gentle.agent_response.v1 JSON only.
Suggest 2-3 valid GENtle shared-shell commands only.
Do not invent slash commands.
```

The first pass is successful only when GENtle parses the reply and every shown
suggestion is a valid GENtle command. Good first suggestions include
`state-summary`, `capabilities`, or `/help`. If the status reports
`AGENT_RESPONSE_PARSE`, the provider answered in a form GENtle could not parse.
Native HTTP transports tolerate a single top-level Markdown `json` code fence,
but they still reject prose wrapped around the JSON.

The prompt box is still an Agent Assistant prompt, not the full GENtle shell.
For convenience it intercepts a small set of local control commands before
contacting the model: `/...` Agent Assistant slash commands plus bare `help`,
`state-summary`, and `capabilities`. Mistyped slash commands fail locally
instead of being reinterpreted as prose. Ollama REPL habits such as bare
`/path/to/file` attachments are not used here; use `/open file PATH`, `/import
file PATH`, or the GUI import path for local sequence files.

After the format probe passes, ask for the real task while keeping execution
reviewed:

```text
I want to retrieve the human FUS gene with isoform annotations from a public
database and present it in GENtle's DNA sequence viewer.

Suggest only valid GENtle shared-shell commands.
Ask before any network/database retrieval.
```

For an empty project, expect the model to suggest a discovery/import path, not
to refer to an existing `seq_id`. Network/database actions should remain
explicitly confirmation-gated. If a model suggests an invalid command, treat
that as a model output problem; GENtle should mark the row as invalid and avoid
running it.

## Agent-facing routes

### 1) CLI and shared shell

Primary binaries and routes:

- `gentle_cli` direct commands (`op`, `workflow`, `genomes`, `tracks`, etc.)
- `gentle_cli shell '...'` shared shell command parser/executor
- GUI Shell panel (same shared shell parser/executor as `gentle_cli shell`)

What it is good for:

- deterministic automation scripts
- CI jobs and reproducible local workflows
- explicit operation sequencing with low overhead

Key properties:

- same engine behavior as GUI/JS/Lua
- machine-readable JSON outputs
- strong reproducibility and auditability

### 2) Machine-facing prose compiler (`agents plan` / `agents execute-plan`)

The planner is the new non-chat route that still accepts prose, but compiles it
into structured executable candidates instead of conversational suggestions.

Entry points:

- CLI/shared shell:
  - `agents plan SYSTEM_ID --prompt TEXT`
  - `agents execute-plan PLAN_JSON_OR_@FILE --candidate-id ID [--confirm]`
- JavaScript/Lua:
  - `plan_agent_system(...)`
  - `execute_agent_plan(...)`
- MCP:
  - `agent_plan`
  - `agent_execute_plan`
- ClawBio wrapper:
  - `mode=agent-plan`
  - `mode=agent-execute-plan`

What it is good for:

- prose-to-typed-action compilation
- auditable compile-then-execute loops
- external orchestrators that want GENtle-native action candidates without
  depending on GUI chat UX

Key properties:

- accepts free prose just like the local assistant
- returns `gentle.agent_plan_result.v1`
- executes stored plans through `gentle.agent_execution_result.v1`
- shell candidates may execute directly through the shared shell executor
- `op`/`workflow` candidates require explicit confirmation when marked mutating
- execution never silently re-plans

### Documentation context for inner helpers

GENtle's inner Agent Assistant should be grounded in GENtle's own
documentation before it proposes commands. The minimum bundle is:

- `docs/glossary.json` for command paths and syntax skeletons
- `docs/cli.md` for operand conventions and shell examples
- `docs/protocol.md` for request/response schemas and execution rules
- `docs/ai_prompt_contract.md` for agent behavior expectations
- `docs/ai_cloning_primer.md`, `docs/ai_task_playbooks.md`, and
  `docs/examples/ai_cloning_examples.md` for biology-first context
- optional compact terminology: `docs/ai_glossary_extensions.json`

This matters most for small local models: glossary placeholders such as
`QUERY`, `ID`, `SEQ_ID`, `ENTRY_ID`, and `PATH` are not enough on their own.
If the model has not been given the relevant documentation, or if the operand
semantics remain unclear, it should ask for the missing identifier or exact path
instead of inventing a GENtle command.

### 3) MCP (`gentle_mcp`)

MCP is the tool-based route for external AI clients that speak JSON-RPC over
stdio.
It is not only transport: MCP is also the standardized discovery/negotiation
surface (`tools/list`, `capabilities`, `help`) used to understand available
GENtle functionality before execution.

Current MCP tool families include:

- engine/state/help tools:
  - `capabilities`
  - `state_summary`
  - `op` (requires `confirm=true`)
  - `workflow` (requires `confirm=true`)
  - `help`
- catalog/introspection tools:
  - `agent_systems`
  - `agent_preflight`
  - `agent_models`
  - `agent_plan`
  - `agent_execute_plan`
  - `reference_catalog_entries`
  - `helper_catalog_entries`
  - `helper_semantics_vocabulary`
  - `host_profile_catalog_entries`
  - `ensembl_installable_genomes`
  - `exon_skip_plan`
  - `exon_skip_materialize`
  - `construct_reasoning_graphs`
  - `construct_reasoning_graph`
  - `construct_reasoning_set_annotation_status`
  - `construct_reasoning_write_annotation`
  - `helper_interpretation`
  - `restriction_site_detail`
  - `blast_async_start`
  - `blast_async_status`
  - `blast_async_cancel`
  - `blast_async_list`
- UI-intent tools:
  - `ui_intents`
  - `ui_intent`
  - `ui_prepared_genomes`
  - `ui_latest_prepared`

What it is good for:

- integrating GENtle into MCP-compatible agent runners
- explicit tool-call loops (`tools/list`, `tools/call`)
- explicit capability negotiation before/alongside execution
- deterministic tool result handling with standard MCP envelopes

Key properties:

- MCP tools are thin wrappers over existing shared engine/shell paths
- no MCP-only biology logic branch
- `tools/list` descriptors expose `mutating: false`, `mutating: true`, or
  `mutating: "external"` so agents can apply their own safety boundary before
  calling a route
- `op`, `workflow`, and materialization-style routes require explicit
  confirmation (`confirm=true`) at the tool boundary where supported
- UI-intent tools are currently non-mutating query/intent routes

#### Intentionally MCP-excluded shell commands

`tools/list` exposes the curated typed MCP surface, not every shared shell
command as a one-command/one-tool mirror. The following glossary-backed shell
commands are intentionally excluded from dedicated MCP tool descriptors for now:
they either require richer route-specific schemas, remain CLI/GUI workflow
surfaces, or are expected to travel through typed `op`/`workflow` payloads once
an agent has selected a deterministic operation.

- `agents ask`
- `align compute`
- `arrange-serial`
- `arrange-set-ladders`
- `arrays inspect-microarray-track`
- `arrays inspect-probe-region-output`
- `arrays import-apt-probe-region-output`
- `arrays probe-regions`
- `arrays project-probe-region-output`
- `arrays project-microarray-track`
- `arrays render-probe-region-output-svg`
- `batch plan`
- `batch run`
- `cache clear`
- `cache inspect`
- `candidates delete`
- `candidates filter`
- `candidates generate`
- `candidates generate-between-anchors`
- `candidates list`
- `candidates macro`
- `candidates metrics`
- `candidates pareto`
- `candidates score`
- `candidates score-distance`
- `candidates score-weighted`
- `candidates set-op`
- `candidates show`
- `candidates template-delete`
- `candidates template-list`
- `candidates template-put`
- `candidates template-run`
- `candidates template-show`
- `candidates top-k`
- `construct-reasoning build-protein-dna-handoff`
- `construct-reasoning export-graph`
- `cutrun export-coverage`
- `cutrun gene-set-regulatory-support`
- `cutrun inspect-regulatory-support`
- `cutrun interpret`
- `cutrun list`
- `cutrun list-read-reports`
- `cutrun prepare`
- `cutrun project`
- `cutrun show-read-report`
- `cutrun status`
- `dbsnp fetch`
- `dotplot compute`
- `dotplot list`
- `dotplot overlay-compute`
- `dotplot show`
- `ensembl-gene fetch`
- `ensembl-gene import-sequence`
- `ensembl-gene list`
- `ensembl-gene show`
- `ensembl-protein fetch`
- `ensembl-protein import-sequence`
- `ensembl-protein list`
- `ensembl-protein show`
- `ensembl-region fetch`
- `export-pool`
- `export-run-bundle`
- `features materialize-repeats`
- `features repeat-cohort`
- `features repeat-overlaps`
- `features repeat-query`
- `features tfbs-score-tracks-svg`
- `features tfbs-summary`
- `features window-cohort-tfbs`
- `flex compute`
- `flex list`
- `flex show`
- `genbank fetch`
- `gene-groups doctor`
- `gene-groups draft`
- `gene-groups list`
- `gene-groups resolve`
- `gene-groups show`
- `gene-sets promoter-cohort`
- `gene-sets resolve`
- `genomes blast`
- `genomes blast-track`
- `genomes extend-anchor`
- `genomes extract-gene`
- `genomes extract-promoter`
- `genomes extract-region`
- `genomes genes`
- `genomes install-ensembl`
- `genomes prepare`
- `genomes preview-ensembl-specs`
- `genomes remove-catalog-entry`
- `genomes remove-prepared`
- `genomes status`
- `genomes update-ensembl-specs`
- `genomes validate-catalog`
- `genomes verify-anchor`
- `gibson apply`
- `gibson preview`
- `guides delete`
- `guides filter`
- `guides filter-show`
- `guides list`
- `guides oligos-export`
- `guides oligos-generate`
- `guides oligos-list`
- `guides oligos-show`
- `guides protocol-export`
- `guides put`
- `guides show`
- `helpers blast`
- `helpers blast-track`
- `helpers ensembl-available`
- `helpers extend-anchor`
- `helpers extract-gene`
- `helpers extract-promoter`
- `helpers extract-region`
- `helpers genes`
- `helpers install-ensembl`
- `helpers prepare`
- `helpers preview-ensembl-specs`
- `helpers remove-catalog-entry`
- `helpers remove-prepared`
- `helpers status`
- `helpers update-ensembl-specs`
- `helpers validate-catalog`
- `helpers verify-anchor`
- `helpers vocabulary doctor`
- `history redo`
- `history status`
- `history undo`
- `import-pool`
- `inspect-feature-expert`
- `ladders export`
- `ladders list`
- `load-project`
- `macros instance-list`
- `macros instance-show`
- `macros run`
- `macros template-delete`
- `macros template-import`
- `macros template-list`
- `macros template-put`
- `macros template-run`
- `macros template-show`
- `mirna catalog-show`
- `mirna explain-seed`
- `mirna scan-target`
- `orthologs promoter-comparison`
- `orthologs resolve-promoter-cohort`
- `panels import-isoform`
- `panels inspect-isoform`
- `panels render-isoform-svg`
- `panels validate-isoform`
- `planning consult cloning`
- `planning protein-expression-handoff`
- `planning objective clear`
- `planning objective set`
- `planning objective show`
- `planning profile clear`
- `planning profile set`
- `planning profile show`
- `planning suggestions accept`
- `planning suggestions list`
- `planning suggestions reject`
- `planning sync pull`
- `planning sync push`
- `planning sync status`
- `primers design`
- `primers design-qpcr`
- `primers export-qpcr-report`
- `primers export-report`
- `primers export-restriction-cloning-handoff`
- `primers list-qpcr-reports`
- `primers list-reports`
- `primers list-restriction-cloning-handoffs`
- `primers preflight`
- `primers prepare-restriction-cloning`
- `primers restriction-cloning-vector-suggestions`
- `primers seed-from-feature`
- `primers seed-from-splicing`
- `primers seed-restriction-cloning-handoff`
- `primers show-qpcr-report`
- `primers show-report`
- `primers show-restriction-cloning-handoff`
- `primers test-cdna-pcr`
- `primers test-cdna-qpcr`
- `primers test-cdna-qpcr-fasta`
- `proteases digest`
- `proteases digest-gel-svg`
- `proteases list`
- `proteases show`
- `protocol-cartoon list`
- `protocol-cartoon render-svg`
- `protocol-cartoon render-template-svg`
- `protocol-cartoon render-with-bindings`
- `protocol-cartoon template-export`
- `protocol-cartoon template-validate`
- `racks apply-template`
- `racks carrier-labels-svg`
- `racks create-from-arrangement`
- `racks fabrication-svg`
- `racks hero-svg`
- `racks isometric-svg`
- `racks labels-svg`
- `racks move`
- `racks move-blocks`
- `racks move-samples`
- `racks openscad`
- `racks place-arrangement`
- `racks set-blocked`
- `racks set-custom-profile`
- `racks set-fill-direction`
- `racks set-profile`
- `racks show`
- `racks simulation-json`
- `reads acquire cancel`
- `reads acquire inspect`
- `reads acquire prepare`
- `reads acquire status`
- `render-dotplot-svg`
- `render-feature-expert-svg`
- `render-lineage-svg`
- `render-pool-gel-svg`
- `render-rna-svg`
- `render-svg`
- `reporters export-corpus`
- `reporters list`
- `reporters recommend`
- `resources benchmark-jaspar`
- `resources inspect-jaspar`
- `resources install-ucsc-rmsk`
- `resources list-jaspar`
- `resources list-publication-datasets`
- `resources prepare-publication-dataset`
- `resources prepare-ucsc-rmsk-index`
- `resources resolve-tf-query`
- `resources status-publication-dataset`
- `resources suggest-ucsc-rmsk-index`
- `resources summarize-jaspar`
- `resources sync-jaspar`
- `resources sync-jaspar-remote-metadata`
- `resources sync-rebase`
- `resources sync-ucsc-rmsk`
- `reverse-translate export-report`
- `reverse-translate list-reports`
- `reverse-translate run`
- `reverse-translate show-report`
- `rna-info`
- `rna-reads align-report`
- `rna-reads batch-map`
- `rna-reads build-transcript-index`
- `rna-reads export-abundance-tsv`
- `rna-reads export-alignment-dotplot-svg`
- `rna-reads export-alignments-tsv`
- `rna-reads export-hits-fasta`
- `rna-reads export-isoform-triage-tsv`
- `rna-reads export-paths-tsv`
- `rna-reads export-report`
- `rna-reads export-sample-sheet`
- `rna-reads export-score-density-svg`
- `rna-reads export-target-quality`
- `rna-reads inspect-alignments`
- `rna-reads inspect-concatemers`
- `rna-reads inspect-gene-support`
- `rna-reads interpret`
- `rna-reads list-reports`
- `rna-reads materialize-hits`
- `rna-reads preflight-isoforms`
- `rna-reads show-alignment`
- `rna-reads show-report`
- `rna-reads summarize-gene-support`
- `routines compare`
- `routines explain`
- `routines list`
- `save-project`
- `screenshot-window`
- `seq-confirm export-report`
- `seq-confirm export-support-tsv`
- `seq-confirm list-reports`
- `seq-confirm run`
- `seq-confirm show-report`
- `seq-primer suggest`
- `seq-trace import`
- `seq-trace list`
- `seq-trace show`
- `sequence create`
- `services delivery-route`
- `services guide`
- `services handoff`
- `services project-preflight`
- `services project-quote`
- `services providers doctor`
- `services providers list`
- `services route-project-source`
- `services status`
- `set-param`
- `splicing-refs derive`
- `tracks import-bed`
- `tracks import-bigwig`
- `tracks import-vcf`
- `tracks tracked add`
- `tracks tracked apply`
- `tracks tracked clear`
- `tracks tracked list`
- `tracks tracked remove`
- `transcripts derive`
- `transcripts residue-genomic-coordinates`
- `uniprot audit-export`
- `uniprot audit-list`
- `uniprot audit-parity`
- `uniprot audit-parity-export`
- `uniprot audit-parity-list`
- `uniprot audit-parity-show`
- `uniprot audit-projection`
- `uniprot audit-show`
- `uniprot compare-ensembl-exons`
- `uniprot compare-ensembl-peptide`
- `uniprot feature-coding-dna`
- `uniprot fetch`
- `uniprot import-swissprot`
- `uniprot list`
- `uniprot map`
- `uniprot projection-list`
- `uniprot projection-show`
- `uniprot resolve-ensembl-links`
- `uniprot show`
- `uniprot transcript-accounting`
- `variant annotate-promoters`
- `variant materialize-allele`
- `variant promoter-context`
- `variant reporter-fragments`

### 4) Agent Assistant bridge (`agents ...` and GUI Agent Assistant)

Agent Assistant runs configured external/internal AI systems and can return:

- assistant text
- follow-up questions
- suggested shell commands with a visible intent title, optional
  `preconditions[]`, optional postcondition-like `expected_outcomes[]`, and
  execution mode `chat|ask|auto`

Entry points:

- CLI/shared shell: `agents list`, `agents ask`, `agents preflight`,
  `agents preflight --live`, `agents discover-models`
- GUI: `File -> Agent Assistant...`

What it is good for:

- interactive planning and translation from human request to deterministic
  commands
- human-in-the-loop execution control
- quick exploration with optional state summary context

Key properties:

- not a direct replacement for deterministic interfaces
- produces suggestions that can be executed through shared shell commands
- optimized for human-in-the-loop chat and suggestion review, not as ClawBio's
  primary machine-facing planning API
- recursion guardrail blocks nested `agents ask`, `agents plan`, and
  `agents execute-plan` execution from suggested commands
- suggested commands can execute shared BLAST routes (`genomes/helpers blast`,
  `genomes/helpers blast-track`) because they use the same parser/executor as
  CLI shell
- current limitation: long-running suggested commands execute synchronously;
  dedicated async job-handle/progress/cancel flow for agent-driven BLAST (and
  future primer-pair multi-BLAST selection) is planned
- ChatGPT/Codex subscriptions do not authenticate the OpenAI API. Native OpenAI
  API mode needs `OPENAI_API_KEY`. The separate `Codex Local` catalog entry
  delegates to the logged-in local Codex CLI through `scripts/codex-agent-bridge`
  and uses Codex/ChatGPT plan limits instead of API billing.
- Native Claude mode uses the Anthropic API directly and needs
  `ANTHROPIC_API_KEY`; Claude Code or Claude.ai subscription/login tokens do
  not authenticate direct Anthropic API calls.
- Native Mistral mode uses the Mistral API directly and needs
  `MISTRAL_API_KEY`; Le Chat or Mistral account login tokens do not
  authenticate direct Mistral API calls.
- `agents preflight` remains config-only by default. `--live` adds the optional
  `gentle.agent_preflight.v1.live_probe` model-list probe without sending a
  generation request.

## CLI vs planner vs MCP vs Agent Assistant

| Topic | CLI/shared shell | Planner | MCP | Agent Assistant |
|---|---|---|---|---|
| Transport | process args/stdin/stdout | shared shell + stored JSON plan | JSON-RPC over stdio | catalog-driven agent transport + shared shell execution |
| Best use | scripts and deterministic automation | prose-in typed compile/execute loops | tool-based external agent integration | interactive assistant with suggestion flow |
| Mutating safety gate | command-level intent | candidate-level confirm for mutating `op`/`workflow` | explicit `confirm=true` on mutating tools | per-suggestion execution policy (`ask`/`auto`) |
| Output model | direct JSON/text/markdown | typed plan/result JSON | MCP envelope + structuredContent | assistant payload + optional execution reports |
| Deterministic parity target | canonical | canonical shell/op/workflow execution | canonical via wrappers | uses canonical shell routes for execution |

## The "agent prompt" offered by GENtle

GENtle provides prompt templates in the GUI Agent Assistant. These templates
help users write better requests, but they are not an execution interface by
themselves.

Current template set includes:

- Structured (recommended)
- Compact intro (no state)
- Candidate between anchors
- BLAST specificity check
- Track import + prioritization
- Macro/template authoring

`Compact intro (no state)` is optimized for live demonstrations and "what can
you do?" prompts. Selecting it in the GUI disables project-state summary
injection for that draft request so the model can answer quickly without
spending tokens on an empty or irrelevant project snapshot.

Important distinction:

- Interface (CLI/MCP/shared shell): executable command format and deterministic
  result
- Prompt template: guidance text for the model request

In short:

- prompt quality affects assistant output quality
- execution determinism comes from underlying shell/engine command formats

## Typical usage patterns

### Pattern A: deterministic script-first

1. Use `gentle_cli` or `gentle_cli shell` directly.
2. Store commands/workflows in version control.
3. Re-run unchanged for reproducible outputs.

### Pattern B: machine-facing planner compile/execute

1. Call `agents plan` or `agent_plan` with prose.
2. Inspect `gentle.agent_plan_result.v1`.
3. Choose one candidate.
4. Execute it later with `agents execute-plan` / `agent_execute_plan`.

### Pattern C: external MCP agent orchestrator

1. Connect to `gentle_mcp`.
2. Discover tools with `tools/list`.
3. Call tools deterministically (`ui_*`, `state_summary`, `op`, `workflow`).
4. Require explicit `confirm=true` for mutating calls.

Use this when the external agent already has its own MCP-capable runtime or
subscription. Point `gentle_mcp --state PATH` at the saved GENtle state file the
agent should inspect; without `--state`, the default is `.gentle_state.json`.

### Pattern D: interactive assistant then execute

1. Use GUI Agent Assistant or `agents ask`.
2. Ask for explicit command suggestions.
3. Execute selected suggestions through shared shell.
4. Keep `ask-before-run` by default for safety.

Example command suggestions (valid through Agent Assistant execution path):

- `genomes blast "Human GRCh38 Ensembl 116" ACGTACGT --task blastn-short --max-hits 20`
- `helpers blast-track "Plasmid pUC19 (online)" ACGTACGT seq1 --track-name primer_offtarget --clear-existing`

## Safety and governance notes

- Keep business/biology logic in shared engine paths only.
- Keep adapters thin (CLI/MCP/GUI assistant should not fork biology behavior).
- Use explicit confirmation for mutating routes.
- ClawBio/OpenClaw should prefer the typed planner boundary or deterministic
  shell/op/workflow routes rather than driving GENtle through `agents ask`.
- Prefer deterministic machine-readable outputs for automation.

## Related manuals

- `docs/gui.md` (GUI usage and Agent Assistant UI details)
- `docs/cli.md` (CLI and MCP operational commands)
- `docs/quickstart_claude.md` (Claude-specific internal and external driving
  scenarios)
- `docs/protocol.md` (protocol details and schema-level definitions)
- `docs/architecture.md` (architecture invariants and parity rules)
