# Agent Bridge

<!-- manual-counterparts:start -->
Related GUI reference: [Agent Assistant](../gui/agent-assistant.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Agent-assistant capability status:

- `gentle_cli`: supported via shared-shell command family (`agents list`,
  `agents ask`, `agents plan`, `agents execute-plan`) and direct forwarding
  (`gentle_cli agents ...`)
- GENtle-local slash aliases are available through the shared shell parser for
  agent convenience, not as OpenClaw/OS commands. They include `/help`, `/list`,
  session history (`/history`, `/undo`, `/redo`),
  GUI file-picker intents (`/open`, `/import`), exact file import
  (`/open file PATH [--id ID]`, `/import file PATH [--id ID]`), non-mutating
  viewer control (`/open sequence-window SEQ_ID`,
  `/close sequence-window SEQ_ID`), pasted sequence creation
  (`/paste sequence --sequence-text DNA [--id ID]`), and
  explicit external fetch aliases (`/fetch genbank|ncbi|uniprot|ensembl*|dbsnp
  ...`). Unknown slash commands return a typed
  `gentle.shell_alias_rejection.v1` diagnostic with supported alternatives.
- `gentle_js`: supported via helper wrappers (`list_agent_systems`,
  `ask_agent_system`, `plan_agent_system`, `execute_agent_plan`) over shared
  shell execution
- `gentle_lua`: supported via helper wrappers (`list_agent_systems`,
  `ask_agent_system`, `plan_agent_system`, `execute_agent_plan`) over shared
  shell execution
- GUI: supported via standalone `Agent Assistant` viewport using the same bridge and command executor


- Shared shell:
  - Commands:
    - `agents list [--catalog PATH]`
    - `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--allow-web-research] [--execute-all] [--execute-index N ...] [--no-state-summary]`

Agent bridge commands:

Conceptual/tutorial companion:

- `docs/tutorial/01-01_agent_interfaces.md` (role map + interface comparison:
  CLI/shared shell vs MCP vs Agent Assistant vs external coding agents).

Recommended Pi Local workflow:

1. Start `pi` outside GENtle once and use Pi's `/login` command to authenticate
   the provider you intend to use. Confirm that `pi --list-models` lists at
   least one provider/model pair. GENtle neither reads nor stores Pi's
   credentials.
2. Check the adapter without starting a model request:
   `cargo run --bin gentle_cli -- agents preflight pi_local_stdio`.
3. Read the exact provider-qualified model ids that GENtle can select:
   `cargo run --bin gentle_cli -- agents discover-models pi_local_stdio`.
4. Ask for help, selecting one of those ids when desired:
   `cargo run --bin gentle_cli -- agents ask pi_local_stdio --model mistral/codestral-latest --prompt "Summarize this project and suggest one safe next GENtle command"`.
   To ground the request in one local text document, include its exact absolute
   path in the prompt, for example:
   `cargo run --bin gentle_cli -- agents ask pi_local_stdio --prompt "Guide me through the GUI smoke test described in $PWD/docs/roadmap.md"`.
   Add `--allow-web-research` when Pi should search and read public pages for
   current information. This is a per-request permission, not a persistent
   provider setting.
5. Inspect the returned command before running it. Keep the default confirmation
   boundary for mutating or network-backed work; use `--execute-index` or
   `--execute-all` only when the corresponding suggestions have been reviewed.

Pi supplies the selected language model, while GENtle supplies the bounded
project summary/facts, validates the structured response, and controls command
execution. Each Pi invocation is ephemeral and runs without sessions, built-in
tools, extension discovery, skills, prompt templates, or project-file discovery.
With `--allow-web-research`, GENtle explicitly loads only its public search and
page-reading extension and allowlists those two tools. Therefore this adapter
helps operate GENtle; it is deliberately not an inner source-code editor.
Every Agent Assistant request also carries bounded `x_local_references`
readiness for prepared, manifest-backed reference genomes. Pi and other
providers receive catalog ids and reusable-component status, but no local
paths. They should prefer a compatible `gene_extraction_ready` reference and a
reviewable `genomes extract-gene` plus `genomes extend-anchor` command chain
before proposing confirmation-gated web retrieval. Merely catalogued genomes
are not represented as installed.

Every request also carries a prompt-matched `x_helper_catalog` projection from
GENtle's bundled helper/vector catalog. It includes up to six matching records
with exact known product/catalogue/accession identity, constraints, source URLs,
and valid inspection/preparation/fetch routes. Matching includes the bounded
recent conversation, so a short follow-up such as “verify that vector” retains
the earlier vector context. This catalog is consulted before web research and
does not imply that a sequence is loaded in the current project.

`agents ask` execution rows include typed, hash-bound `feedback` and an
`execution_status_counts` summary. `executed=false` includes commands rejected
before dispatch. Starting an asynchronous job does not mean the analysis
completed; inspect its status/result separately. CLI invocations remain
stateless and do not automatically recall previous execution receipts. See the
[agent request contract](../protocol.md) for GUI session feedback and privacy limits.

Agent Assistant's initial background-command routes now use the host-owned
service described under [Background Commands](../gui.md#background-commands).
Standalone CLI commands remain blocking by caller policy; this does not add a
cross-process job server. Shared workflow/macro progress includes completed/total
step counts and checks cooperative cancellation between short steps. Cancelled
execution returns an error, not a successful partial result. Existing synchronous
nontransactional prefixes and transactional rollback semantics remain unchanged.

Its `x_introspection` projection uses the current prompt to prioritize matching
subject/report/enzyme identities within the existing per-type fact budget.
`selection_explanations` records each included fact's priority and matches.
The shared CLI/MCP path has no active-view assumption; relevance ordering
neither changes canonical facts nor supplies execution approval.

For Pi Local, `--allow-web-research` adds `x_web_access.enabled=true`. Pi may
search the public web and read arbitrary public HTTP/HTTPS text pages, not just
an allowlist of scientific domains. Localhost, private/reserved networks, URL
credentials, binary responses, unsafe redirects, oversized responses, shell,
and filesystem access remain blocked. GENtle records the actual queries and
pages supplied to Pi as `x_web_research`; these sources are shown separately
from model prose and do not grant execution authority. Search terms derived
from the prompt may leave the machine; do not enable this option for sensitive
projects. The agent contract forbids transmitting sequences, local paths,
credentials, personal data, or confidential project details.

GUI-originated Agent Assistant requests additionally carry a bounded
`x_gui_context` with the live recent-project list, executable tutorial
projects, tutorial guides/references, a deterministic prompt-matched tutorial
shortlist, and Configuration sections. The shortlist uses the current request
and, for a context-light follow-up, up to three recent user messages. It reports
matched terms/fields and is a retrieval aid rather than a biological-confidence
score. The CLI parser and MCP expose the
corresponding UI-intent commands, including `ui open tutorial-guide
TUTORIAL_ID`, but a headless `agents ask` process has no private GUI
recent-project list and therefore cannot resolve an opaque recent item id by
itself. Tutorial ids and Configuration section names remain stable; the actual
open/focus action is applied only when a GUI host receives the intent.

When an `agents ask` prompt explicitly names a supported absolute text-document
path, the shared request builder includes bounded UTF-8 content as
`x_local_documents` (four files maximum, 128 KiB each, 256 KiB total, with
digest/truncation provenance and one prioritized Markdown-link level). This is
the same behavior as the GUI attachment control. The content is sent to the
selected provider; it is not a general filesystem permission. FASTA, GenBank,
and other sequence inputs still use GENtle import commands.
Use Pi separately, outside GENtle, for repository-editing work. If the GUI does
not inherit Pi's executable path, set `PI_BIN` to the exact executable before
starting GENtle.

- `agents list [--catalog PATH]`
  - Lists configured agent systems from catalog JSON.
  - Default catalog path: `assets/agent_systems.json`.
  - Includes availability status (`available`, `availability_reason`) so callers
    can skip systems that are not currently runnable.
- `agents preflight SYSTEM_ID [--live] [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N]`
  - Returns read-only transport/runtime metadata for one configured system.
  - This is the preferred CLI parity route for external orchestrators that want
    availability, endpoint, and runtime-limit diagnostics without starting a
    chat/planning request.
  - By default this is config-only and does not contact the provider.
  - `--live` adds `live_probe` to `gentle.agent_preflight.v1` for
    `native_openai`, `native_anthropic`, `native_mistral`,
    `native_openai_compat`, and `pi_local_stdio`.
  - Native HTTP live probes only use model-discovery endpoints (`GET /models`,
    with `/v1/models` fallback for OpenAI-compatible roots). They do not make a
    chat/completion/responses request and do not intentionally generate tokens.
  - `pi_local_stdio` uses a local command-shape probe: GENtle invokes Pi with
    the baseline ephemeral no-tools/no-session flags used by `scripts/pi-agent-bridge`
    and appends `--help`, so Pi parses the intended flags but no prompt is sent
    to a model.
  - `live_probe.status_class` is one of:
    `ok`, `missing_key`, `auth_failed`, `quota_or_billing`, `model_missing`,
    `endpoint_unreachable`, `unsupported_transport`, `provider_error`.
  - Quota/billing is reported only when the provider returns that error during
    the model-list probe.
- `agents discover-models SYSTEM_ID [--catalog PATH] [--base-url URL]`
  - Discovers model ids for one `native_openai`, `native_anthropic`,
    `native_mistral`, or `native_openai_compat` system using its effective
    endpoint configuration.
  - For `codex_local_stdio`, reads visible model ids from the logged-in Codex
    CLI's local `models_cache.json`; it does not read credentials or send a
    provider request.
  - For `pi_local_stdio`, runs `pi --list-models` and returns the
    provider-qualified ids visible to the installed Pi CLI. If none are
    available, start Pi outside GENtle and use `/login` before trying again.
- `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--allow-web-research] [--execute-all] [--execute-index N ...] [--no-state-summary]`
  - Invokes one configured agent system and returns message/questions/suggested shell commands.
  - `--base-url` sets a per-request runtime endpoint override (maps to
    `GENTLE_AGENT_BASE_URL`) for native transports.
  - `--model` sets a per-request model override (maps to `GENTLE_AGENT_MODEL`)
    for native transports, `codex_local_stdio`, and `pi_local_stdio`; each
    local bridge forwards it as its CLI's `--model MODEL` argument.
  - `--allow-web-research` opts one request into the selected system's declared
    public-web tools. It is supported by `pi_local_stdio`; other systems reject
    the flag rather than silently ignoring it.
  - `--timeout-secs` sets per-request timeout override (maps to
    `GENTLE_AGENT_TIMEOUT_SECS`) for stdio/native transports.
  - `--connect-timeout-secs` sets per-request HTTP connect timeout (maps to
    `GENTLE_AGENT_CONNECT_TIMEOUT_SECS`) for native transports.
  - `--read-timeout-secs` sets per-request read timeout (maps to
    `GENTLE_AGENT_READ_TIMEOUT_SECS`) for stdio/native transports.
  - `--max-retries` sets per-request transient retry budget (maps to
    `GENTLE_AGENT_MAX_RETRIES`; `0` means no retries).
  - `--max-response-bytes` caps adapter response size per attempt (maps to
    `GENTLE_AGENT_MAX_RESPONSE_BYTES`).
  - for `native_openai_compat`, model must resolve to a concrete value (catalog
    model or `--model`); `unspecified` is rejected.
  - by default, project context includes both `state_summary` and a bounded
    `x_introspection` fact projection from the same engine state; the latter
    carries complete fact-type counts, explicit truncation metadata, and
    read-only retrieval routes
  - `--no-state-summary` disables both forms of project-context injection in
    the request
  - External adapter responses must be strict `gentle.agent_response.v1` JSON.
  - Unknown canonical response fields are rejected (extensions must use `x_` or `x-` prefix).
  - Adapter calls retry transient failures with exponential backoff before returning an error.
  - Suggested commands are executed only when explicitly selected
    (`--execute-all`, `--execute-index`) or when `--allow-auto-exec` is enabled,
    the suggestion intent is `auto`, and the shared parser classifies it as
    non-mutating. Mutating `auto` suggestions remain unexecuted until an
    explicit selection.
  - Recursive `agents ask`, `agents plan`, and `agents execute-plan` execution
    from suggested commands is blocked by design.
  - Failures use deterministic error prefixes for scripting, e.g.
    `AGENT_ADAPTER_UNAVAILABLE`, `AGENT_ADAPTER_TRANSIENT`,
    `AGENT_RESPONSE_PARSE`, `AGENT_RESPONSE_VALIDATION`.
  - Supported transports:
    - `builtin_echo` (offline/demo)
    - `external_json_stdio` (local bridge executable)
    - `native_openai` (built-in OpenAI HTTP adapter, requires `OPENAI_API_KEY`)
    - `native_anthropic` (built-in Anthropic Claude HTTP adapter, requires
      `ANTHROPIC_API_KEY`)
    - `native_mistral` (built-in Mistral HTTP adapter, requires
      `MISTRAL_API_KEY`)
    - `native_openai_compat` (built-in local OpenAI-compatible adapter for
      Jan/Msty/Ollama-style `/chat/completions` endpoints; key optional)
      - endpoint host/port come from catalog `base_url` (or `--base-url` if set);
        GENtle does not silently switch to a different host/port
      - `msty_mlx_local_compat_template` targets Msty/MLX Knife servers at
        `http://localhost:11973/v1`; use it when the Msty gateway on
        `http://localhost:11964/v1/models` is reachable but returns no model
        ids
- `agents plan SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--max-candidates N] [--no-state-summary] [--no-mutating-candidates]`
  - Accepts prose like `agents ask`, but returns typed `gentle.agent_plan_result.v1`
    candidates instead of chat-oriented suggestion rows.
  - Current candidate kinds are `shell`, `op`, `workflow`, and `ui_intent`.
  - This is the recommended machine-facing prose boundary for MCP/ClawBio-style
    orchestration.
- `agents execute-plan PLAN_JSON_OR_@FILE --candidate-id ID [--confirm]`
  - Executes one stored plan candidate without re-planning.
  - `shell` candidates run through the shared shell executor.
  - `op` / `workflow` candidates require `--confirm` when marked mutating.
  - Nested `agents ask` / `agents plan` / `agents execute-plan` recursion from
    candidate shell payloads is blocked.

Prompting users should send to agents:

- Prefer compact structured prompts so local models stay reliable.
- Use this template:

```text
Objective:
<one clear goal>

Context:
<sequence/genome/helper ids>

Inputs:
- anchors/coordinates:
- feature labels/kinds:

Constraints:
- length:
- GC range:
- motifs/sites:
- strand assumptions:

Output wanted:
- exact commands
- validation checklist

Execution policy:
chat-only | ask-before-run | allow-auto-exec
```

- For local models without cloning background, preload:
  - `docs/glossary.json`
  - `docs/cli.md` (especially Operand metavariable conventions)
  - `docs/protocol.md`
  - `docs/ai_cloning_primer.md`
  - `docs/ai_task_playbooks.md`
  - `docs/ai_prompt_contract.md`
  - `docs/examples/ai_cloning_examples.md`
  - optional: `docs/ai_glossary_extensions.json`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Screenshot Capture Policy](screenshot-capture.md) | [Next: ClawBio and OpenClaw](clawbio.md)
<!-- manual-page-nav:end -->
