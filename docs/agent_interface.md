# GENtle Agent Interface

Last updated: 2026-08-26

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

1. Open `Configuration -> Agent Systems`.
2. Choose the provider profile (`Local Model`, `Codex Local`, OpenAI, Claude,
   Mistral, or another catalog entry).
3. For local OpenAI-compatible services such as Ollama, Jan, or Msty, set the
   base URL and click `Discover Models`; then pick one concrete discovered
   model.
4. Click `Test Setup`. This checks endpoint/model reachability without sending
   a generation request.
5. Close Configuration and open `File -> Agent Assistant...`.
6. Leave `Project summary` / `Include state summary` unchecked for the first
   prompt. In an empty project there is no useful project context to send, and
   small local models often behave better with the shortest possible request.
7. Ask a response-format probe before asking biology:

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
All model-facing transports tolerate a single top-level Markdown `json` code
fence, but they still reject prose wrapped around the JSON. External stdio
adapters remain responsible for returning the versioned response schema and
should emit bare JSON whenever possible.

Successful turns appear under `Conversation`. GENtle stores that transcript
with the current project and supplies a bounded recent portion to later agent
requests, so a follow-up answer such as `homo_sapiens` remains available even
when Codex Local starts a fresh process. `Clear Conversation` removes the
stored transcript. API keys and session-only provider settings are never part
of project conversation metadata. Prompt text is stored verbatim, so do not
paste secrets into the prompt.

The prompt box is still an Agent Assistant prompt, not the full GENtle shell.
For convenience it intercepts a small set of local control commands before
contacting the model: `/...` Agent Assistant slash commands plus bare `help`,
`state-summary`, and `capabilities`. Mistyped slash commands fail locally
instead of being reinterpreted as prose. Ollama REPL paths are not a biological
sequence-import command here; use `/open file PATH`, `/import file PATH`, or
the GUI import path for local sequence files.

Requests made from the GUI also receive `gentle.agent_gui_context.v1`. This is
the authoritative bounded mirror of the same Open Recent Project and Open
Tutorial Project catalogs visible in the File menu, plus the five global
Configuration sections. Consequently, questions such as "Which previous
projects are known?" or "Which tutorial projects are known?" should return the
actual rows rather than infer an empty list from an empty current project.
Recent-project rows contain path-free display metadata, file availability,
size/time, and an opaque id; only the live GUI can resolve that id:

```text
ui open recent-project RECENT_ITEM_ID
ui open tutorial-project CHAPTER_ID
ui open configuration agent-systems
ui open configuration microarrays
```

These remain reviewed UI intents. A missing recent file is listed as
unavailable, a stale id fails closed, and opening Configuration only navigates
to the selected tab. The user still applies credentials, executable paths, and
other global settings through the existing Configuration confirmation model.
Headless `agents ask` calls do not invent a GUI recent-project catalog; they can
still parse and report the same intent commands for a future GUI host.

For textual guidance documents, use the separate bounded document-context
route. Click `Attach text document...`, or include an exact absolute path such
as `/Users/name/checkout/docs/roadmap.md` in a natural-language prompt. GENtle
reads that supported UTF-8 text file itself and adds `x_local_documents` to the
request with byte limits, truncation state, and SHA-256 provenance. One level of
prioritized local Markdown guide/runbook links can be included under the same
directory tree. Pi/Codex/local HTTP models still receive no filesystem tool and
cannot browse adjacent files. Review the visible paths and selected provider:
the copied contents are sent to that provider.

A useful walkthrough request is:

```text
Please guide me through the GUI smoke test described in
/absolute/path/to/gentle_rs/docs/roadmap.md.

Use validated GENtle ui/display commands where they exist. Otherwise give me
one exact manual step at a time and ask what I see before continuing.
```

The resulting command suggestions remain review rows in GENtle. Running one
uses the same parser and confirmation boundary as any other Agent Assistant
suggestion. If a checklist action has no registered GUI intent, the agent can
coach the human but must not claim to have pressed or observed the control.

For visual GUI problems, use `Agent help` in the window where the problem is
visible. A normal click captures that exact GENtle viewport before opening Agent
Assistant; on macOS, the right-click menu can optionally capture the complete
native GENtle window after Screen Recording permission is granted. Agent
Assistant previews the image locally and asks for explanatory prompt text. The
image is sent only when the user clicks `Ask agent`, and only to catalog systems
that explicitly declare image support. Bare prompt paths are deliberately not
treated as attachments.

An image-capable agent may return one optional `screenshot_request` containing
only a bounded stable `id` and a human-readable `reason`. This creates a consent
card, not a capture. The user chooses one currently registered GENtle content
window and may allow exactly one cross-platform egui capture or decline. The
approval is bound to the response turn, provider, project generation, and
selected viewport; it expires on change and is never persisted. The resulting
image enters the same local preview and still requires `Ask Agent` before
transport. Agents cannot select arbitrary windows, paths, coordinates, native
window ids, desktop capture, or capture commands. Auto-run cannot approve it.

The separate controls do not imply that pixels are inherently more private
than structured state. A complete project summary can expose more scientific
context, while a screenshot is often necessary to interpret layout, colour,
direction, and disabled controls. Both are disclosures to the selected provider
and should be reviewed for the current professional context.

### Reproducible screenshot-consent smoke

Use an Agent Assistant catalog system with
`supports_image_attachments: true`. On Linux, the same desktop flow can be run
under Xvfb after building the GUI:

```bash
xvfb-run -a -s "-screen 0 1440x900x24" target/debug/gentle
```

Open a main/sequence content window and Agent Assistant, then ask the configured
agent: `For this consent-flow test, request one screenshot with id
consent-smoke-1 and explain that you need to inspect the visible controls.` The
expected sequence is:

1. The response shows a consent card with the reason, provider, and registered
   window selector; there is no preview or pending attachment yet.
2. `Decline` removes the card without a capture. Repeat the request, choose
   `Main Window`, and click `Allow one screenshot` once.
3. One local preview appears with source title, `egui.viewport` backend, pixel
   dimensions, and digest-backed attachment metadata. A second click/replay
   cannot create another capture.
4. The image remains local and no agent task starts until `Ask Agent` is clicked.
   A provider lacking image support instead keeps approval disabled.

Changing provider/project, clearing the conversation, closing the selected
window, or starting another request before approval should expire the card.
This smoke does not use or enable `screenshot-window` and does not exercise the
user-only native macOS context-menu capture.

Local command results are presented directly beneath the prompt status. In the
GUI, `/help` also opens the built-in Help window at `Shell Commands`; an
optional topic is copied into Help search. `/list` renders a compact current-
project overview with sequence ids, names, lengths, topology, container counts,
and arrangement counts. Its object rows have context menus with concrete shared-
shell actions: sequence rows can open the sequence, query annotations, or scan
restriction sites, while container rows expose the same exhaustive/subset toggle
as the main project overview. A `Show complete project overview` action focuses
the main table/graph for analyses and the other project-object context menus.
`Copy JSON` preserves the complete engine-owned result when a machine-readable
handoff is needed.

Session history is also available locally: `/history` reports undo/redo
availability, while `/undo` and `/redo` apply the corresponding session-local
transition. A command typed and run by the user is explicit. A command proposed
by an agent must be reviewed and run by the user; GENtle refuses to
automatically execute undo or redo even if the provider marks it `auto`. In the
GUI, history transitions remain disabled while background jobs are active and
refresh dependent sequence windows and caches after success.

Viewer-window requests are GUI intents, not project-data mutations. To control a
catalogued GENtle tool/dialog, use `ui open TARGET`, `ui focus TARGET`, or
`ui close TARGET` such as `ui close pcr-design`. To control a DNA sequence
viewer while keeping the sequence loaded, use
`ui open sequence-window SEQ_ID`, `ui focus sequence-window SEQ_ID`, or
`ui close sequence-window SEQ_ID`; `/open sequence-window SEQ_ID` and
`/close sequence-window SEQ_ID` are convenience aliases. Do not use
delete/remove commands for window-close requests. For Ensembl gene retrievals
from the Agent Assistant, add `--no-open` when the fetched sequence should be
loaded into project state without automatically opening a DNA sequence viewer.
To work with the active DNA-window selection, use
`ui selection sequence-window SEQ_ID --range START..END` (0-based,
end-exclusive) or `ui selection sequence-window SEQ_ID` to inspect it. To
toggle display layers, use `display show TARGET` or `display hide TARGET`, for
example `display show tfbs` or `display hide restriction-enzymes`.

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
- `docs/biological_extension_guide.md` for coding agents/contributors adapting
  a biological workflow rather than merely invoking it
- `docs/ai_cloning_primer.md`, `docs/ai_task_playbooks.md`, and
  `docs/examples/ai_cloning_examples.md` for biology-first context
- optional compact terminology: `docs/ai_glossary_extensions.json`

This matters most for small local models: glossary placeholders such as
`QUERY`, `ID`, `SEQ_ID`, `ENTRY_ID`, and `PATH` are not enough on their own.
If the model has not been given the relevant documentation, or if the operand
semantics remain unclear, it should ask for the missing identifier or exact path
instead of inventing a GENtle command.

When a request is biologically clear but no registered route implements it,
the inner assistant must distinguish that capability gap from missing project
input. It should classify the request as parameterization, composition, or a
missing primitive and return a concise implementation brief rather than a
fictional command. A coding-capable outer agent can then follow
`docs/biological_extension_guide.md`; the inner assistant itself remains an
execution/planning client and does not edit GENtle source.

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
  - `runtime_status`
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
- `runtime_status` returns the same `gentle.runtime_status.v1` shape as
  `introspect runtime` for the MCP process: process-local live `frames[]` plus
  observed `activities[]` from existing persisted/project ledgers. It does not
  create status files or see another process's unsaved GUI-local frames.
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
- `collections run digest`
- `collections run export-pool`
- `collections run primer-specificity`
- `collections run restriction-scan`
- `collections run tfbs-scan`
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
- `digest`
- `display`
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
- `features restriction-scan`
- `features tfbs-scan`
- `features tfbs-score-tracks-svg`
- `features tfbs-summary`
- `features window-cohort-tfbs`
- `flex compute`
- `flex list`
- `flex show`
- `facts eval`
- `facts graph`
- `genbank fetch`
- `gene-groups doctor`
- `gene-groups draft`
- `gene-groups list`
- `gene-groups resolve`
- `gene-groups show`
- `gene-sets produce co-regulated`
- `gene-sets produce direct-list`
- `gene-sets produce ontology-assignment`
- `gene-sets create-pool`
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
- `genomes adopt-blast-resource`
- `genomes prepare-blast-resource`
- `genomes inspect-blast-resource`
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
- `introspect all`
- `introspect capabilities`
- `introspect facts`
- `introspect readiness`
- `introspect verify-effects`
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
- `primers design-terminal-exon-rt-pool` (the typed
  `DesignTerminalExonRtPrimerPool` operation remains available through the
  generic MCP `op` tool)
- `primers design-transcript-assay-panel`
- `primers compose-gene-isoform-study-workflow-batch`
- `primers inspect-gene-isoform-study-reuse`
- `primers execute-gene-isoform-study-workflow`
- `primers execute-gene-isoform-study-workflow-batch`
- `primers design-qpcr`
- `primers experimental-handoff` (the typed `BuildExperimentalAssayHandoff` operation remains available through the generic MCP `op` tool)
- `primers export-qpcr-report`
- `primers export-report`
- `primers export-restriction-cloning-handoff`
- `primers export-transcript-assay-fallback`
- `primers export-transcript-assay-panel`
- `primers import-external-pairs`
- `primers list-qpcr-reports`
- `primers list-reports`
- `primers list-restriction-cloning-handoffs`
- `primers list-transcript-assay-fallbacks`
- `primers list-transcript-assay-panels`
- `primers primerbank search`
- `primers primerbank show`
- `primers primerbank test-cdna`
- `primers screen-variants` (the typed `ScreenPrimerVariants` operation remains available through the generic MCP `op` tool)
- `primers preflight`
- `primers prepare-restriction-cloning`
- `primers restriction-cloning-vector-suggestions`
- `primers seed-from-feature`
- `primers seed-from-splicing`
- `primers seed-restriction-cloning-handoff`
- `primers show-qpcr-report`
- `primers show-report`
- `primers show-restriction-cloning-handoff`
- `primers show-transcript-assay-fallback`
- `primers show-transcript-assay-panel`
- `primers transcript-assay-specificity-plan`
- `primers transcript-assay-specificity-finalize`
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
- `resources import-co-regulated-cache`
- `resources import-gene-list-cache`
- `resources import-ontology-assignment-cache`
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
- `rna-reads allele-hash-screen`
- `rna-reads align-report`
- `rna-reads batch-map`
- `rna-reads build-transcript-index`
- `rna-reads export-abundance-tsv`
- `rna-reads export-alignment-dotplot-svg`
- `rna-reads export-alignments-tsv`
- `rna-reads export-dexseq-annotation-gff`
- `rna-reads export-dexseq-counts-tsv`
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
- `rna-reads verify-dexseq`
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
- `ui selection`
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

For `rna-reads allele-hash-screen`, agents may supply
`--from-rna-report REPORT_ID` to let the engine resolve target-gene-aligned
retained rows and their stored sequences through the same gene-support cohort
logic used by the audit command. Salmon `unmapped_names` and mapping SAM files
are ID selectors, not complete sequence stores: when either Salmon flag is
used, the command also requires explicit `--read-file` or `--read-pair`
inputs. The v3 report exposes overlapping per-read `source_origins[]` and
aggregate `source_provenance[]` counts so downstream imbalance interpretation
does not lose its cohort basis. Its representation verdict is threshold-based
on phased allele-informative depth (`--min-informative-reads` and the balanced
band flags); report-derived retained-read support is contextual expression
weight, and the optional binomial value is advisory rather than a significance
gate. Agents must describe the result as sequence representation, not
biological significance or causation.

### 4) Agent Assistant bridge (`agents ...` and GUI Agent Assistant)

Agent Assistant runs configured external/internal AI systems and can return:

- assistant text
- follow-up questions
- suggested shell commands with a visible intent title, optional
  `preconditions[]`, optional machine-readable `precondition_expr`, optional
  postcondition-like `expected_outcomes[]`, optional machine-readable
  `expected_effects[]`, and execution mode `chat|ask|auto`

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
- `/undo` and `/redo` suggestions require explicit confirmation and cannot be
  selected by automatic-suggestion execution
- suggested commands can execute shared BLAST routes (`genomes/helpers blast`,
  `genomes/helpers blast-track`) because they use the same parser/executor as
  CLI shell
- suggested-command `precondition_expr` fields are evaluated advisory-only
  against the current `gentle.project_fact_graph.v1`; the GUI shows readiness
  as `ready`, `blocked`, or `unknown` while the shared shell still owns hard
  validity checks
- the Agent Assistant system prompt receives its known project-fact vocabulary
  from the same registry used by the evaluator, instead of a separately
  maintained prose-only list
- when project-context injection is enabled, Agent Assistant and `agents ask`
  also receive `x_introspection` from the same engine snapshot as
  `state_summary`; it contains at most 128 concrete facts, complete fact-type
  counts, explicit truncation metadata, and the read-only routes used to obtain
  more detail
- the prompt's compact introspection card explains the distinction between a
  registered fact name and a currently projected fact. This guidance is
  embedded at runtime because native HTTP models and the isolated Codex Local
  bridge cannot be assumed to read files from the GENtle checkout
- use `facts graph` and `facts eval ...` in shared shell/CLI contexts when an
  external agent needs the same deterministic project-fact projection
- open-world facts such as restriction-site absence require explicit evidence
  bundles today; until a persisted evidence ledger lands, GUI readiness without
  attached evidence remains `unknown` for those facts
- current limitation: long-running suggested commands execute synchronously;
  dedicated async job-handle/progress/cancel flow for agent-driven BLAST (and
  future primer-pair multi-BLAST selection) is planned
- ChatGPT/Codex subscriptions do not authenticate the OpenAI API. Native OpenAI
  API mode needs `OPENAI_API_KEY`. The separate `Codex Local` catalog entry
  delegates to the logged-in local Codex CLI through `scripts/codex-agent-bridge`
  and uses Codex/ChatGPT plan limits instead of API billing. The bridge honors
  `CODEX_BIN`, then `codex` on `PATH`, and on macOS checks both the bundled
  ChatGPT `Contents/Resources/codex` CLI and `Contents/MacOS/ChatGPT`
  executable before the legacy standalone Codex app location. Selecting Codex
  Local immediately presents `Codex default`; GENtle reads the CLI's visible
  model ids from its local metadata cache and forwards a selected id through
  `GENTLE_AGENT_MODEL` to `codex --model`.
- Pi can be used as a second local Agent Assistant harness through the
  `Pi Local (uses Pi provider login)` catalog entry. GENtle invokes
  `scripts/pi-agent-bridge`, which runs one ephemeral `pi --print` request with
  tools, sessions, extensions, skills, prompt templates, and project-context
  discovery disabled. GENtle keeps the conversation and project-fact context;
  Pi supplies the selected provider/model. Model discovery uses
  `pi --list-models`, and authentication remains in Pi's own `/login` flow.
  `agents preflight pi_local_stdio --live` and GUI `Test Setup` also run a
  non-generating command-shape probe: Pi is invoked with the same ephemeral
  no-tools/no-session flags used for real Agent Assistant requests plus
  `--help`, so unsupported local Pi CLI flags are reported before any prompt is
  sent to a model. GENtle does not read or copy Pi's credential store.
- This Pi entry is not a source-editing developer agent. A future coding mode
  would need a separate workspace, permission, diff-review, and test contract;
  it must not inherit the Agent Assistant's reviewed-command permissions.
- Native Claude mode uses the Anthropic API directly and needs
  `ANTHROPIC_API_KEY`; Claude Code or Claude.ai subscription/login tokens do
  not authenticate direct Anthropic API calls.
- Native Mistral mode uses the Mistral API directly and needs
  `MISTRAL_API_KEY`; Le Chat or Mistral account login tokens do not
  authenticate direct Mistral API calls.
- `agents preflight` remains config-only by default. `--live` adds the optional
  `gentle.agent_preflight.v1.live_probe`: native HTTP systems use model-list
  probes, while Pi Local uses the no-generation command-shape probe described
  above.

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
- `ui selection sequence-window fus_live --range 100..250`
- `display show tfbs`

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
