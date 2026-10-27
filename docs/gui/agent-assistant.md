## Agent Assistant

<!-- manual-counterparts:start -->
Related CLI reference: [Agent Bridge](../cli/agents.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

GENtle provides a standalone `Agent Assistant` window for structured support from
configured external/internal agent systems.

### Background Commands

Prepared Run buttons and typed commands use the same execution policy. Genome
preparation, `op`, `workflow`, workflow/candidate macros and template runs, and a
verified single gene-study workflow run on an engine snapshot in the background.
Run does not replace your draft or ask the model to reinterpret the command.
You can type canonical commands such as `genomes prepare ID` or `workflow @file`
directly into the prompt for these routes; no invented slash alias is needed.
You can navigate, inspect the live project, or send a new model request while
these commands run. Direct commands remain available during a model request.

For work lasting more than 200 ms, the assistant shows a command preview, phase,
available completed-step count and Cancel. Counts are not a time estimate.
The local execution log retains the exact command; the running row limits its
preview so an inline JSON workflow cannot consume the whole window.
Cancel requests cooperative termination; it does not claim the worker has
already stopped. The terminal receipt distinguishes computed operations from
committed operations. Failure/cancellation discards unpublished engine changes,
but downloaded resources and exported files can remain. Concurrent structural
edits or project replacement reject stale results; viewport-only edits are kept.

The host accepts at most four active commands and retains 32 results/receipts.
If the project lock or capacity is unavailable, submission reports that it was
not admitted, rather than waiting. Results above 16 MiB require a file-export
workflow. These are initial audited routes, not a claim that every shell command
is nonblocking. Broader command/adapter migration remains
in the [responsive-command plan](../asynchronous_command_execution_plan.md).

Agent Assistant additionally submits `genomes blast` / `helpers blast` through
the shared BLAST job scheduler, for both Run buttons and typed commands. The
query/options are unchanged; the response is an admitted job, not a completed
search. Explicit `blast-start`, `blast-status`, `blast-list` and `blast-cancel`
commands use the same lifecycle and remain available during a model request.
Tool preflight happens on the worker, with diagnostics exposed by status.
The BLAST FIFO allows 64 waiting jobs beyond the configured running slots;
the general four-worker command service above still has no waiting queue.
Use `blast-status JOB_ID --with-report` for the final result. Compact project
receipts do not retain that full report, so export it before closing the host.

Conceptual/tutorial companion:

- `docs/tutorial/01-01_agent_interfaces.md` (who runs what where, and how Agent
  Assistant differs from CLI, MCP, and external coding agents).
- external MCP clients discover available GENtle tools/capabilities with
  `tools/list` and `capabilities` before calling operations.

Access:

- main menu: `File -> Agent Assistant...`
- command palette action: `Agent Assistant`
- shortcut: `Cmd+Shift+A`
- `Agent help` in every GENtle-owned window, including the main project and
  sequence windows, specialist-window navigation, Help/History/Jobs/JASPAR
  diagnostics, and sequence-owned Dotplot, Splicing Expert, RNA-read Mapping,
  Isoform Expert, and Promoter design workspaces

Window screenshot help:

- Click `Agent help` in the window where the problem is visible. The capture
  request is bound to that viewport before Agent Assistant is opened, so the
  assistant window cannot accidentally replace the intended target.
- GENtle shows the captured image, source-window title, dimensions, size, and
  capture backend in Agent Assistant. The user can remove it or add explanatory
  text before sending.
- Only catalog systems declaring image-attachment support can receive the
  screenshot. GENtle disables submission for other systems rather than silently
  dropping the image.
- A successful conversation turn retains image digest, dimensions, source
  window, and backend but omits the temporary local file path and image bytes.
- The normal viewport capture does not require macOS Screen Recording access.
  The optional right-click `Capture complete macOS window` path includes native
  window chrome, is restricted to the current GENtle process, and opens a clear
  permission/settings recovery path when macOS blocks it.
- An agent can ask for visual context by returning one validated request id and
  reason, but it cannot name a path, window id, coordinate, or capture command.
  The consent card lists only GENtle's current open-window registry and excludes
  Agent Assistant itself. Approval expires when the response/provider/project
  changes, the selected window closes, another agent request starts, or the
  conversation is cleared; replay and double-click do not capture twice.
- Project-state context and screenshots are separate provider disclosures, not
  a privacy ranking. The state summary may contain richer scientific content;
  the screenshot contributes visible layout, colour, direction, and control
  state. Review the provider and both controls according to the project.

Behavior:

- Prompt template `Evidence-guided regulatory reporter study` helps the agent
  collect perturbation context, responsive genes, assembly, TSS policy,
  evidence sources, and exact vector identity before suggesting commands. It
  keeps response, motif, occupancy, and reporter-result evidence separate and
  stops before mutation or materialization.
- before the first conversation, a small `Start here` card offers two direct
  GUI actions: open one or more existing sequence documents, or open GENtle's
  general Configuration window; neither action is sent to an agent
- provider setup is separated from conversation:
  - `Configuration -> Agent Systems` contains catalog, provider, key, endpoint,
    model, timeout/retry, model-discovery, and `Test Setup` controls
  - `File -> Agent Assistant...` contains the project conversation, prompt,
    response details, and reviewed command execution
- loads systems from catalog JSON (default `assets/agent_systems.json`)
- `Agent Systems` is the first heading in its Configuration tab; system
  selection is the primary dropdown, while the optional `Quick start` shortcuts
  are collapsed below it
- the selected system id is persisted in `~/.gentle_gui_settings.json` and
  restored on the next launch; API keys and other session overrides are not
  persisted
- unavailable systems remain selectable and show the reason
- `Quick start` buttons expose common first-run routes:
  - `Use OpenAI API`
  - `Use Claude API`
  - `Use Mistral API`
  - `Use Local Model (no OpenAI API billing)`
  - `Use Demo Echo`
- `External Agent / MCP` shows the `gentle_mcp` command for users whose
  external agent already supports MCP and already has its own account/runtime
  subscription
  - the command includes the active saved project path when available
  - otherwise it shows the default MCP state path, `.gentle_state.json`
- the OpenAI quick-start path uses `OPENAI_API_KEY` and talks to the API
  directly; it does not reuse a ChatGPT/Codex chat session as authentication
- the `Codex Local (uses Codex CLI login)` catalog entry is different: it calls
  the local Codex CLI through `scripts/codex-agent-bridge`, reusing the user's
  existing Codex CLI/App login instead of `OPENAI_API_KEY`
- the `Pi Local (uses Pi provider login)` entry calls
  `scripts/pi-agent-bridge`; it uses provider authentication configured by Pi,
  not the API-key field in GENtle, and runs Pi without built-in tools or direct
  project-file access. The Agent Assistant offers a session-only `Allow public
  internet research` checkbox for this capable system. When checked, only
  GENtle's co-shipped public search/page tools are enabled for that request
- `Attach text document...` beside the prompt inserts an exact absolute path
  and shows every detected document path before submission. A supported text
  path written directly in the prompt is detected in the same way.
  - GENtle copies bounded UTF-8 content into the request; it does not grant the
    provider a filesystem tool
  - the warning above the prompt is intentional: copied document content is
    sent to the currently selected provider, including cloud providers
  - limits are four documents, 128 KiB each, and 256 KiB total; one level of
    prioritized local Markdown guide/runbook links can be included under the
    selected document's directory tree
  - sequence files remain project inputs and should be opened/imported, not
    attached as prose context
- the Claude quick-start path uses `ANTHROPIC_API_KEY` and talks to the
  Anthropic API directly
- the Mistral quick-start path uses `MISTRAL_API_KEY` and talks to the Mistral
  API directly
- this window is intentionally the local chat-oriented assistant surface; the
  typed prose compiler/executor (`agents plan` / `agents execute-plan`) is the
  headless CLI/MCP/ClawBio route for machine-facing compile/execute loops
- `OpenAI API key` / `Anthropic API key` / `Mistral API key` in
  `Configuration -> Agent Systems` is a session-only override for the selected
  native provider
  - OpenAI keys usually begin with `sk-...`; Anthropic keys usually begin with
    `sk-ant-...`; Mistral keys should be generated in La Plateforme
  - click `Clear Key` to remove it from current session
  - the key is not persisted to disk by GENtle settings
- native-provider credential precedence is: session key field, provider
  environment variable, then the provider's conventional token file
  (`~/.codex_token`, `~/.claude_token`, or `~/.mistral_token`)
  - GENtle reads token files into memory for the current process only and never
    copies their contents into project or GUI settings
  - the Configuration tab reports the active source, a missing/invalid token
    file, and token files readable beyond their owner; `Reload credentials`
    refreshes the file snapshot
  - Codex Local is intentionally excluded from this fallback: it uses the
    installed Codex CLI/App login rather than an OpenAI API token
- `Base URL override` field is a session-only endpoint override for
  `native_openai`, `native_anthropic`, `native_mistral`, and
  `native_openai_compat`
  - use this for local endpoints such as `http://localhost:11964` or
    `http://localhost:11964/v1`
  - click `Clear URL` to remove it from current session
- `Model override` field is a session-only model-name override for
  `native_openai`, `native_anthropic`, `native_mistral`, and
  `native_openai_compat`
  - use this to force a concrete model id
  - `unspecified` means no override
  - click `Clear Model` to remove it from current session
- `timeout_sec` is a session-only request timeout override
  - applies to stdio and native agent transports
  - maps to `GENTLE_AGENT_TIMEOUT_SECS`
  - empty/`0` means default timeout
- `connect_timeout_sec` / `read_timeout_sec` are session-only transport controls
  - map to `GENTLE_AGENT_CONNECT_TIMEOUT_SECS` and `GENTLE_AGENT_READ_TIMEOUT_SECS`
  - empty/`0` keeps defaults
- `max_retries` / `max_response_bytes` are session-only guardrails
  - map to `GENTLE_AGENT_MAX_RETRIES` and `GENTLE_AGENT_MAX_RESPONSE_BYTES`
  - `max_retries=0` disables retries for the request
  - empty keeps defaults
- if model remains `unspecified`, GENtle blocks requests until you pick a
  discovered model or set a concrete override
- `Test Setup` runs an inline preflight against the currently selected system
  and session overrides before you send a prompt
  - shows availability, resolved base URL/model, endpoint candidates, runtime
    limits, and warnings such as missing `OPENAI_API_KEY` or `ANTHROPIC_API_KEY`
    or `MISTRAL_API_KEY`
  - for `native_openai`, `native_anthropic`, `native_mistral`, and
    `native_openai_compat`, it also runs a live non-generating model-list probe
    (`GET /models`, with `/v1/models` fallback for OpenAI-compatible or root
    provider URLs)
  - for `Pi Local`, it runs a non-generating command-shape probe that asks Pi
    for local help after the baseline ephemeral no-tools/no-session flags;
    this catches unsupported Pi CLI flags before an assistant prompt is sent
  - the live probe reports `status_class`, attempted/selected endpoints or
    local command shape, reachability, authentication, model-list parsing, and
    whether the selected model was present when the transport exposes that
    information
  - common next actions are shown for missing key, authentication failure,
    endpoint unreachable, unspecified model, and selected model not found
- `Discover Models` queries the current endpoint and populates a model dropdown
  for explicit selection
- prompt templates are available via one-click `Insert` / `Append` buttons
  before the prompt editor
- the `Compact intro (no state)` prompt template is meant for live demos and
  quick "what can you do?" questions; inserting or appending it disables
  project-state summary injection for that draft request
- optional `Include state summary` injects current project summary context and
  a bounded `x_introspection` fact projection from the same engine snapshot
  - at most 128 concrete facts are included, together with complete fact-type
    counts and explicit truncation metadata
  - the agent receives read-only routes for `introspect facts`, `facts graph`,
    `facts eval`, `introspect readiness`, and `introspect verify-effects` when
    it needs information beyond the bounded projection
  - missing open-world facts are unknown rather than false; this context does
    not bypass shared parser/engine validation or confirmation requirements
- every request includes a bounded `x_local_references` inventory of prepared,
  manifest-backed reference genomes even when `Include state summary` is off
  - only catalog identity and component readiness are exposed; local paths are
    omitted and inventory construction performs no downloads
  - the agent is instructed to prefer a compatible local gene extraction plus
    strand-aware anchor-extension chain over external retrieval
  - catalog-only rows are not presented as installed references
- every request includes a prompt-matched view of GENtle's helper/vector
  catalog. This lets the assistant recognize existing records such as Promega
  pGL4.10[luc2] (E6651, GenBank AY738222.1) before guessing or searching the
  web; a catalog record is not presented as a loaded project sequence
- successful public-web research appears under a collapsible source section
  with the actual searches and pages recorded by the Pi bridge
- public-web research is off by default because prompt-derived search terms may
  leave the machine; keep it off for sensitive projects, and do not ask the
  agent to transmit sequences, local paths, credentials, personal data, or
  confidential project details
- GUI requests also include `gentle.agent_gui_context.v1`, independently of
  `Include state summary`
  - `recent_projects[]` mirrors the current `File -> Open Recent Project...`
    list with an opaque item id, path-free filename/parent label, existence,
    size, modification time, current-project flag, and exact `open_command`
  - `tutorial_projects[]` comes from the same generated tutorial manifest used
    by `File -> Open Tutorial Project...`, including chapter/example identity,
    title/summary, group/tier, use cases, learning objectives, concepts,
    prerequisites, expected outcomes, GUI-acceptance profile, online status,
    review status, and explicit truncation counts
  - `tutorial_guides[]` contains the remaining catalogued walkthroughs and
    references with exact `ui open tutorial-guide TUTORIAL_ID` commands; these
    open teaching text in Help rather than constructing a worked project
  - `tutorial_recommendations[]` is a deterministic shortlist of at most five
    rows matched against the current request, with recent user intent used for
    context-light follow-ups; every row reports its retrieval score, matched
    terms/fields, prerequisites, and exact command. The score is not biological
    confidence
  - `configuration_sections[]` covers External Applications, Agent Systems,
    Microarrays, Graphics, and Language with exact `ui open configuration
    SECTION` commands
  - `ui open recent-project ITEM_ID` resolves only against the live host list;
    missing files and stale ids fail closed, and no absolute saved-project path
    is sent to the provider
  - opening a Configuration section is navigation, not a silent settings
    mutation; existing Apply/Cancel and credential handling remain in force
- optional `Allow auto execute` only applies to suggestions marked with `auto`
  that the shared parser classifies as non-mutating; mutating suggestions stay
  visible and require the user to click `Run`
- suggestion cards with `precondition_expr` are evaluated against the live
  project fact graph plus the current GUI-host availability on every frame; an
  unmet or unknown precondition dims the complete card and disables `Run`
  until it becomes satisfied
  - `ui.host_available` is therefore true in the live Agent Assistant window,
    so exact tutorial/configuration commands from `x_gui_context` remain
    runnable after confirmation even when the current project is empty
  - the same fact gate is checked again immediately before manual or automatic
    execution, so stale UI state and `execution=auto` cannot bypass it
  - a downstream card can therefore remain visible while an earlier command
    creates its required sequence/report, then become runnable without asking
    the agent for a replacement response
- after the user states an interest, GENtle ranks the tutorial catalog before
  contacting the model; the agent should offer only one to three rows from
  that shortlist with short reasons and ask one focused question when no row
  matches. Reviewed, non-stale, offline material is preferred when relevance
  is comparable
- introductory replies describe Agent Assistant as a goal-directed workflow
  guide before listing commands: it knows the documented GUI controls and can
  invoke their underlying functionality through shared parity interfaces,
  bind tutorials to user data and lab constraints, and compose a new
  reviewable path when no established workflow fits
  - parity concerns the underlying action, not literal mouse-click emulation;
    missing GUI/CLI/script/agent projections are parity gaps to report and fix
  - confirmations and other safety boundaries continue to require the user
  - biological knowledge supplied by the selected model remains advisory and
    distinct from GENtle results, supplied evidence, and confirmed constraints
- successful conversation turns are shown in chronological order and stored as
  `gentle.agent_conversation.v1` metadata with the current project
  - the next request receives the 12 most recent turns, so ephemeral transports
    such as Codex Local remember explicit follow-up answers such as species
  - up to 50 successful turns are retained; oldest turns are discarded first
  - API keys, endpoint overrides, and other provider configuration are never
    stored in the conversation record; prompt text is stored verbatim, so do
    not paste secrets into the prompt
  - a conversation-only unsaved project counts as user content, so GENtle asks
    before closing it without saving
  - `Clear Conversation` removes the stored turns and latest response
- `Ask Agent` runs in background and reports status in `Background Jobs`
- The bounded fact context prioritizes subjects named in the current prompt
  and, when known, the active DNA window's sequence. It preserves a fair share
  of every available fact type; included rows carry deterministic selection
  explanations. A named subject outranks the active view. No sequence is
  guessed when the active window is not a DNA viewer, and selection priority
  is not evidence of biological importance or readiness.
- Helper/vector suggestions prioritize the current question, then the most
  recent user message, then its assistant reply. Older topics do not influence
  catalog matching; brief follow-up questions can still refer to the last reply.
- The next project-context request includes bounded, session-only execution
  receipts. Blocked, failed, partially successful, dispatched, running and
  completed actions remain distinct. Ensembl metadata retrieval followed by a
  failed sequence import is `partial`; starting a BLAST job is not completion.
  Deferred UI actions remain `dispatched` until separately inspected. Feedback
  does not yet poll all GUI/background activities automatically.
  - Receipts send hashes and typed outcomes, not raw command outputs or errors.
    Inspect current results with GENtle's read-only commands.
  - At most 100 receipts are retained, including rejections; suggestions are
    correlated with the latest 12 turns. Project switches and Clear Conversation
    discard them. Saving a project does not persist them.
  - Manually typed commands still populate the result panel and local execution
    log. Their receipts have no suggestion binding, so they are not sent to
    the model as anonymous completions; missing feedback is not proof of failure.
  - Recording a receipt cannot mark a project dirty or clear redo. Later edits
    can require rechecking a result without erasing its historical completion.
- when an accepted UI intent is opening a DNA Sequence Viewer, Agent Assistant
  shows that opening state until deferred sequence loading and the first
  content paint complete; the sequence loading placeholder polls at a bounded
  rate instead of forcing an unbounded repaint loop
- `Ctrl+Return` while the prompt editor is focused is equivalent to clicking
  `Ask Agent`; plain `Return` still inserts a new line
- GENtle-local slash aliases are deliberately small and parser-validated:
  - `/help [TOPIC]` opens the built-in Help window at `Shell Commands`; an
    optional topic is placed in the Help search field
  - `/list` shows a readable `Command result` panel with current sequences,
    lengths/topology, containers, and arrangements; right-clicking a sequence or
    container shows parser-validated commands for that object, and `Show complete
    project overview` focuses the main table/graph with its broader object context
    menus
  - successful local read-only commands retain their latest structured output
    in that panel and offer `Copy JSON`; they no longer disappear behind the
    generic execution-log status
  - `/history` shows session-local undo/redo availability
  - `/undo` and `/redo` use the same guarded transitions as `Edit -> Undo` and
    `Edit -> Redo`; they are disabled while background jobs are active and
    refresh open sequence windows after a successful transition
  - agent-suggested `/undo` or `/redo` commands always require the user to
    click `Run`; the automatic-suggestion option cannot execute them
  - `/open` and `/import` open the same sequence-file picker as
    `File -> Open Sequence...`
  - `/open file PATH [--id ID]` and `/import file PATH [--id ID]` import an
    exact user-provided sequence file and open the resulting sequence window
  - `ui open TARGET`, `ui focus TARGET`, and `ui close TARGET` control
    catalogued GENtle tool/dialog windows without changing project data
  - `ui open sequence-window SEQ_ID`, `ui focus sequence-window SEQ_ID`, and
    `ui close sequence-window SEQ_ID` control the DNA sequence viewer for a
    loaded sequence without deleting that sequence from project state;
    `/open sequence-window SEQ_ID` and `/close sequence-window SEQ_ID` are
    convenience aliases
  - `ui selection sequence-window SEQ_ID --range START..END` sets the
    DNA-window selection for a loaded sequence without changing sequence data;
    `ui selection sequence-window SEQ_ID` reports the current selection
  - `display show TARGET`, `display hide TARGET`, and
    `display visibility TARGET on|off` update shared project display state for
    feature/display targets such as `tfbs`, `restriction-enzymes`,
    `gene-features`, `cds-features`, `repeat-features`, and
    `methylation-sites`
  - `/paste sequence --sequence-text DNA [--id ID]` creates a sequence from
    explicit pasted IUPAC DNA text
  - `/fetch genbank`, `/fetch ncbi`, `/fetch uniprot`, `/fetch ensembl*`, and
    `/fetch dbsnp` normalize to existing external fetch routes and should be
    confirmation-gated
  - `/fetch ensembl ... --assembly NAME --flank-bp N` verifies the resolved
    assembly and requests strand-aware genomic context around the gene;
    `--no-open` imports the record without automatically opening a DNA sequence
    viewer
  - `/grep`, `/find`, `/ls`, `/new`, `/example`, and vague file-discovery
    requests are rejected rather than treated as ClawBio/OpenClaw commands
- response panel can include:
  - assistant message text
  - follow-up questions
  - width-bounded suggestion cards whose message, command, preconditions,
    expected outcomes, and rationale wrap instead of requiring horizontal
    scrolling
  - a `Run` action at the start of every suggestion card; explanatory
    (`execution=chat`), empty, or invalid commands show a disabled action and
    the reason directly in the card
  - invalid or invented suggestions are shown as `Invalid GENtle command` and
    cannot be run
  - `Copy Response JSON` for copying the latest strict agent-response JSON
    payload to the clipboard
- execution is always per suggestion (row-run, explicit all, or guarded
  non-mutating auto); there is no global always-execute mode
- each executed suggestion is logged with status/output in the same window

OpenAI setup (explicit):

1. Open `Configuration -> Agent Systems`.
2. Choose system `OpenAI GPT-5 (native HTTP)` from the dropdown.
3. Paste your API key into `OpenAI API key` (format `sk-...`), export
   `OPENAI_API_KEY`, or place the token in `~/.codex_token`. The Configuration
   tab shows which source is active.
4. Click `Test Setup` to confirm the key/base URL/model resolve correctly.
   This checks model discovery only; it does not intentionally send a
   token-generating request.
5. Open `File -> Agent Assistant...`, enter prompt text, and click `Ask Agent`.
6. If you prefer environment setup instead of GUI key field, launch GENtle with:

```bash
export OPENAI_API_KEY=sk-...
cargo run --bin gentle
```

Codex Local setup (uses Codex CLI login, no API key):

1. Install/sign in to Codex outside GENtle, for example with the Codex app/CLI
   login flow. This may use a browser or device-code style confirmation.
2. Ensure the bridge in this checkout is reachable. The default source-tree
   catalog entry runs `scripts/codex-agent-bridge`, so source-checkout use works
   when GENtle is launched from the repository root. If GENtle is launched as a
   packaged `.app` or from another working directory, either add the repo
   `scripts/` directory to `PATH` and set the catalog command to
   `codex-agent-bridge`, or use an overlay catalog with the absolute bridge
   path.
3. If the Codex CLI itself is not on the GUI process `PATH`, set `CODEX_BIN` to
   the executable path. Current macOS ChatGPT app installs normally expose the
   bundled CLI at `/Applications/ChatGPT.app/Contents/Resources/codex`; GENtle
   also checks `/Applications/ChatGPT.app/Contents/MacOS/ChatGPT`. Older
   standalone Codex app installs may expose
   `/Applications/Codex.app/Contents/Resources/codex`.
4. Open `Configuration -> Agent Systems` and choose
   `Codex Local (uses Codex CLI login)` from the system dropdown.
5. The model selector appears immediately with `Codex default`. GENtle then
   reads the visible model ids from the logged-in Codex CLI's local
   `models_cache.json`; choose one of those ids, keep `Codex default`, or enter
   an explicit custom model override. `Discover Models` repeats the local-cache
   read without sending a prompt.
6. Click `Test Setup`, then open `File -> Agent Assistant...` and send a small
   prompt. `Test Setup` confirms the
   bridge executable/runtime/model settings. Reading the local model metadata
   does not contact the provider, so Codex login/quota is confirmed by the first
   actual `Ask Agent` request.

Notes:

- Codex Local uses the user's Codex/ChatGPT plan limits, not OpenAI API billing.
- The bridge runs Codex in an empty temporary directory by default so project
  files are not implicitly sent as working-tree context; GENtle passes only the
  explicit Agent Assistant prompt, optional state summary, and bounded stored
  conversation.
- The bridge asks Codex for strict `gentle.agent_response.v1` JSON and then
  validates it before GENtle sees the response.
- A selected Codex model is passed to the CLI as `codex --model MODEL`.
  `Codex default` omits that flag and lets the installed CLI choose its
  configured default.
- If a Finder/Spotlight-launched macOS app cannot find Codex or cannot access
  the same login/keychain context, launch GENtle from a terminal as a
  workaround or set explicit `CODEX_BIN`/catalog paths.

Pi Local setup (multiple providers, no GENtle-managed API key):

1. Install Pi and run it once outside GENtle. In Pi, use `/login` to configure
   the provider you want. GENtle does not inspect `~/.pi/agent/auth.json`.
2. Confirm `pi --list-models` shows at least one provider/model row. If Pi is
   not on the GUI process `PATH`, set `PI_BIN` to its executable path; a common
   Apple Silicon Homebrew path is `/opt/homebrew/bin/pi`.
3. Open `Configuration -> Agent Systems` and choose
   `Pi Local (uses Pi provider login)`.
4. Click `Discover Models` and choose a provider-qualified id such as
   `mistral/codestral-latest`, or leave `Pi default` selected to use Pi's
   configured default.
5. Click `Test Setup`, then open `File -> Agent Assistant...` and ask a small
   orientation question.
6. For documentation-guided help, click `Attach text document...` or name the
   exact absolute text-document path in the prompt. For example, ask Pi to guide
   the manual GUI smoke described by `docs/roadmap.md`; GENtle supplies the
   bounded roadmap/runbook text while Pi remains unable to browse other files.
7. When current external information is needed, enable `Allow public internet
   research` before sending that prompt. Leave it off for turns that need only
   GENtle and project context.

The Pi bridge uses `--no-session` because GENtle stores the project conversation
and sends its bounded recent turns with each request. It also disables built-in
tools, context files, extension discovery, skills, and prompt templates; the
Agent Assistant can therefore suggest validated GENtle commands but cannot edit
files or run operating-system commands through Pi. Public-web opt-in explicitly
loads only `gentle_web_search` and `gentle_web_fetch`; those tools reject local
and private networks and do not expose files, credentials, or shell execution.
GENtle supplies Pi with the same shared command contract used by native
providers, including exact local genome-extraction/anchor-extension and Ensembl
assembly/flank syntax, so Pi does not need to guess parser routes. It can return
reviewable `ui ...` or `display ...` suggestions for exposed
GUI intents. Where no intent exists, it must describe one manual step and ask
what the user sees rather than claiming to have clicked or observed it.
Execution receipts intentionally contain status and hashes, not scientific
result rows. Each local structured command result therefore offers **Use
reviewed result in next prompt**: after checking the JSON for sequences, local
paths or other sensitive data, the user can add the complete bounded result to
the retained draft. The pretty-printed JSON is limited to 128 KiB, and the entire
resulting draft to 192 KiB of UTF-8, including the existing draft, command and
all framing text. Rejection leaves the draft unchanged. Nothing is sent until
**Ask Agent** is clicked, oversized results are rejected rather than truncated,
and the handoff grants no approval for a subsequent command. This supports
iterative feature discovery,
transcript-panel feasibility/design and result review without letting the model
infer identifiers or coverage from receipt hashes. TSS previews retain their
narrower typed handoff. A source-editing inner coding agent would be a separate
future mode with a different approval boundary.

Claude setup (explicit):

1. Open `Configuration -> Agent Systems`.
2. Click `Use Claude API` or choose `Claude Sonnet (native Anthropic HTTP)`.
3. Paste your Anthropic Console API key into `Anthropic API key` (format
   `sk-ant-...`), export `ANTHROPIC_API_KEY`, or place that API key in
   `~/.claude_token`. Claude Code or Claude.ai subscription/login tokens are
   not Anthropic API keys; GENtle warns when a file has the shape of such a
   login token.
4. Click `Test Setup` to confirm the key/base URL/model resolve correctly.
   This checks model discovery only; it does not intentionally send a
   token-generating request.
   - GENtle flags obvious Claude Code/Claude.ai OAuth tokens before contacting
     Anthropic; final validation still comes from the live model-list probe.
5. Open `File -> Agent Assistant...`, enter prompt text, and click `Ask Agent`.
6. If you prefer environment setup instead of GUI key field, launch GENtle with:

```bash
export ANTHROPIC_API_KEY=sk-ant-...
cargo run --bin gentle
```

Mistral setup (explicit):

1. Open `Configuration -> Agent Systems`.
2. Click `Use Mistral API` or choose `Mistral Large (native Mistral HTTP)`.
3. Paste your Mistral La Plateforme API key into `Mistral API key`, export
   `MISTRAL_API_KEY`, or place that API key in `~/.mistral_token`. Le Chat or
   Mistral account login tokens are not Mistral API keys.
4. Click `Test Setup` to confirm the key/base URL/model resolve correctly.
   This checks model discovery only; it does not intentionally send a
   token-generating request.
5. Open `File -> Agent Assistant...`, enter prompt text, and click `Ask Agent`.
6. If you prefer environment setup instead of GUI key field, launch GENtle with:

```bash
export MISTRAL_API_KEY=...
cargo run --bin gentle
```

External MCP handoff:

1. Use this path when your outside assistant already supports MCP.
2. Save the current GENtle project if the external assistant should inspect the
   same state.
3. Copy the `External Agent / MCP` command from
   `Configuration -> Agent Systems`, for example:

```bash
gentle_mcp --state .gentle_state.json
gentle_mcp --state path/to/project.gentle.json
```

Notes:

- ChatGPT/Codex subscriptions are not OpenAI API keys.
- MCP lets the external agent call GENtle tools over stdio; GENtle is not
  making OpenAI API calls for that route.
- Unsaved GUI-only edits are not visible to MCP until saved to the state path.

Local LLM setup (Jan/Msty/OpenAI-compatible endpoint):

1. Open `Configuration -> Agent Systems`.
2. Select one of:
   - `Local Llama (OpenAI-compatible)`
   - `Jan Local (template)`
   - `Msty MLX Local (template)`
   - `Msty Local (template)`
3. Set `Base URL override` to your local endpoint, e.g.
   `http://localhost:11973/v1` for Msty MLX or `http://localhost:11964`
   for the Msty gateway profile.
4. Optionally set `timeout_sec` for slow local models (for example `600`).
5. Click `Discover Models` and select one discovered model from the dropdown
   (or set `Model override` directly, for example `deepseek-r1:8b`).
6. Click `Test Setup` to confirm the endpoint/model resolve correctly.
7. Optionally set persistent defaults in `assets/agent_systems.json`.
8. If your local service expects no key, keep `OpenAI API key` empty.
9. Open `File -> Agent Assistant...` and ask the agent as usual.
10. For local root URLs (such as `http://localhost:11964`), GENtle will try both:
   - `/chat/completions`
   - `/v1/chat/completions`
11. If `http://localhost:11964/v1/models` returns `data: null` while an MLX
    server exposes models on `http://localhost:11973/v1/models`, choose
    `Msty MLX Local (template)` or set that URL as the base override.

Common failure interpretation:

- `AGENT_ADAPTER_UNAVAILABLE ... status=429 ... code=insufficient_quota`
  - connection/auth path works, but OpenAI API project quota is exhausted
  - the Agent Assistant, Background Jobs, and `Test Setup` live-probe panels
    show direct links for:
    - OpenAI usage
    - OpenAI billing
  - `Test Setup` treats non-OK live probes as the visible overall status even
    when the static catalog/key configuration is otherwise available
  - fix billing/quota at:
    - `https://platform.openai.com/usage`
    - `https://platform.openai.com/settings/organization/billing/overview`
  - this is separate from ChatGPT plan credits

What to send the agent (recommended):

- Domain bootstrap docs for local models:
  - `docs/ai_cloning_primer.md`
  - `docs/ai_task_playbooks.md`
  - `docs/ai_prompt_contract.md`
  - `docs/examples/ai_cloning_examples.md`
- Optional compact glossary extension:
  - `docs/ai_glossary_extensions.json`
- For a related biological workflow that may not yet be a command:
  - choose prompt template `Adapt a biological workflow`
  - the assistant first classifies the request as an existing parameter,
    composition of existing capabilities, or missing engine primitive
  - unsupported requests produce an implementation brief, not a fabricated
    command; coding agents can continue with
    `docs/biological_extension_guide.md`

Copy/paste prompt template for users:

```text
Objective:
<what you want, in one sentence>

Context:
<project sequence/genome/helper IDs and why this region matters>

Inputs:
- seq_id / genome_id / helper_id: ...
- anchors or coordinates: ...
- feature labels/kinds: ...

Constraints:
- length: ...
- GC range: ...
- motifs/restriction sites to require/avoid: ...
- strand assumptions: ...

Output wanted:
- plan
- exact gentle_cli shell commands
- verification checks

Execution policy:
chat-only | ask-before-run | allow-auto-exec
```

Prompt examples:

- \"Generate 20 bp candidates between TP53 start/end on `grch38_tp53`, keep GC 40-80%, maximize distance from CDS boundaries, return top 25 with exact commands, ask-before-run.\"
- \"Check BLAST specificity for sequence `ACGT...` against `homo_sapiens_grch38_ensembl_116`, task `blastn-short`, max hits 20, chat-only.\"
- \"Import `/path/peaks.bed.gz` onto anchored `tp53_region`, then keep candidates within 200 bp of imported track features, provide commands and stop before execution.\"

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: GENtle Shell (GUI)](gentle-shell-gui.md) | [Next: About GENtle](about-gentle.md)
<!-- manual-page-nav:end -->
