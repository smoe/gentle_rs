# GENtle Agent Assistant and Agent Interfaces Tutorial

> Type: `operational reference tutorial`
> Status: `manual/reference`
> Audience: users operating GENtle through the in-app Agent Assistant, CLI/shared shell, MCP, or external coding agents.

Last updated: 2026-05-18

This tutorial explains how to let an AI assistant help with GENtle without
giving up reproducibility. The most important idea is simple:

- the model proposes or selects GENtle commands
- GENtle executes reviewed commands through the shared shell or engine
- replayability comes from those concrete commands, not from the free-form chat

You can open this page inside GENtle from:

1. `Help -> Tutorials -> Agent Assistant and Agent Interfaces Tutorial`
2. `Help -> Agent Interface` for the shorter protocol overview

Protocol-level details live in `docs/agent_interface.md` and
`docs/protocol.md`. This page is the practical walk-through.

## 1) The mental model

GENtle has one deterministic execution core.

The same command should mean the same thing when it is reached from:

- GUI Shell
- `gentle_cli shell "..."`
- MCP tools that delegate to shared shell or engine operations
- Agent Assistant suggested commands
- JavaScript/Lua shell wrappers

An agent can help you find, combine, or explain commands. It should not become
a second biology implementation.

## 2) The internal path: Agent Assistant

Use this when you are already in GENtle and want the assistant to propose next
steps.

Open it from:

```text
File -> Agent Assistant...
```

The in-app Agent Assistant sends a request to one configured agent system and
expects a `gentle.agent_response.v1` reply:

```json
{
  "schema": "gentle.agent_response.v1",
  "assistant_message": "Short explanation for the user.",
  "questions": [],
  "suggested_commands": [
    {
      "title": "Inspect restriction sites",
      "rationale": "This is a read-only check.",
      "command": "features restriction-scan --sequence-text GAATTC... --enzyme EcoRI",
      "execution": "ask"
    }
  ]
}
```

Important: `suggested_commands[].command` contains GENtle shared-shell
commands, not operating-system shell commands. GENtle runs them internally after
you review them.

## 3) First internal test: offline demo

This does not contact any provider.

1. Open `File -> Agent Assistant...`.
2. Click `Use Demo Echo`.
3. Click `Test Setup`.
4. Enter a prompt such as:

```text
ask: features restriction-scan --sequence-text GAATTCGGGCCCGGGCCCGAGCTCGAATTC --enzyme EcoRI --enzyme SmaI
```

5. Click `Ask Agent`.
6. Review the suggested command.
7. Run it only if the command is the GENtle command you intended.

The demo echo is not intelligent. It is useful because it exercises the same
response and command-review loop without requiring an API key.

## 4) Provider quick starts

The quick-start buttons configure the built-in native HTTP transports.

### OpenAI

1. Click `Use OpenAI API`.
2. Paste an OpenAI Platform API key into `OpenAI API key`, or set
   `OPENAI_API_KEY` before launching GENtle.
3. Click `Test Setup`.

ChatGPT or Codex subscriptions are not OpenAI API keys.

### Claude

1. Click `Use Claude API`.
2. Paste an Anthropic Console API key into `Anthropic API key`, or set
   `ANTHROPIC_API_KEY`.
3. Click `Test Setup`.

Claude Code or Claude.ai subscription/login tokens are not Anthropic API keys.
GENtle tries to catch obvious wrong-token shapes before contacting Anthropic,
but the live model-list probe is the final check.

### Mistral

1. Click `Use Mistral API`.
2. Paste a Mistral La Plateforme API key into `Mistral API key`, or set
   `MISTRAL_API_KEY`.
3. Click `Test Setup`.

Le Chat or Mistral account login tokens are not Mistral API keys.

### Local model

1. Click `Use Local Model (no OpenAI API billing)`.
2. Set `Base URL override` to the local OpenAI-compatible endpoint, for example
   `http://localhost:11964`.
3. Click `Discover Models`.
4. Pick a concrete discovered model.
5. Click `Test Setup`.

Local roots may expose `/chat/completions` or `/v1/chat/completions`; GENtle
tries both for OpenAI-compatible local services.

## 5) What Test Setup actually tests

`Test Setup` is intentionally non-generating for native HTTP transports.

It checks:

- catalog system selection
- provider/key availability
- resolved base URL
- resolved model
- runtime limits
- model-list endpoint reachability
- whether the selected model appears in the returned model list

It does not intentionally send a chat/completion/responses request. Quota or
billing is reported only if the provider returns that error during model-list
probing.

Common outcomes:

| Status | Meaning | Next action |
|---|---|---|
| `ok` | Endpoint, auth, model list, and selected model are consistent. | Ask the assistant. |
| `missing_key` | GENtle has no provider API key. | Paste a session key or set the provider env var. |
| `auth_failed` | Provider rejected the key. | Check key type and provider account. |
| `model_missing` | Model list worked but the selected model was absent or unspecified. | Discover/pick a model or set model override. |
| `endpoint_unreachable` | Base URL could not be reached. | Start the local server or correct Base URL. |
| `provider_error` | Model-list response was malformed or unexpected. | Inspect the provider message. |

## 6) A safe first real prompt

Use a read-only request first:

```text
Please inspect this short sequence for EcoRI and SmaI sites. Return one GENtle
shared-shell command only. Use execution "ask".

Sequence:
GAATTCGGGCCCGGGCCCGAGCTCGAATTC
```

A good command suggestion looks like:

```text
features restriction-scan --sequence-text GAATTCGGGCCCGGGCCCGAGCTCGAATTC --enzyme EcoRI --enzyme SmaI
```

The command is read-only and uses inline sequence text. You can run it from the
Agent Assistant suggestion row, the GUI Shell, or the terminal:

```bash
cargo run --quiet --bin gentle_cli -- shell 'features restriction-scan --sequence-text GAATTCGGGCCCGGGCCCGAGCTCGAATTC --enzyme EcoRI --enzyme SmaI'
```

That replay command is the audit trail. The model prompt is only how you got
there.

## 7) Execution policy

Prefer this policy while learning:

- keep `Auto-run suggestions marked as 'auto'` off
- ask the model for `execution: "ask"`
- run one suggestion at a time
- read the command before running it
- copy important commands into lab notes or issue comments

GENtle blocks recursive agent execution. Suggested commands may not silently run
`agents ask`, `agents plan`, or `agents execute-plan` again.

## 8) The external path: Claude, Codex, or another coding agent

Use this when the assistant is outside GENtle and can run terminal commands or
edit the checkout.

Recommended discovery commands:

```bash
scripts/dev-gentle-cli doctor --agent
scripts/dev-gentle-cli capabilities
scripts/dev-gentle-cli shell "help"
```

For persistent project state, give the external agent an explicit state file:

```bash
scripts/dev-gentle-cli --state path/to/project.gentle.json shell "state-summary"
```

External coding agents are useful for:

- testing documented workflows
- patching docs or source
- adding missing CLI/shared-shell routes
- comparing GUI, CLI, and MCP behavior

They are not the same as the in-app Agent Assistant. They operate your
development environment; the in-app assistant operates through GENtle's agent
catalog and reviewed shared-shell suggestions.

## 9) The MCP path

Use MCP when your external assistant supports typed tools.

Start the server with an explicit state path:

```bash
gentle_mcp --state .gentle_state.json
```

The MCP client should first call:

1. `tools/list`
2. `capabilities`
3. the specific read-only or mutating tool

Mutating operation/workflow routes require explicit confirmation. MCP is best
when the outside agent should see a typed tool palette instead of inventing
plain text commands.

## 10) Choosing the right interface

| Need | Use |
|---|---|
| Human asks for help inside GENtle | Agent Assistant |
| Human knows exact command | GUI Shell or `gentle_cli shell` |
| Script or CI | `gentle_cli` / `scripts/dev-gentle-cli` |
| Tool-calling external agent | `gentle_mcp` |
| Repository edits or missing route implementation | Codex/Claude Code-style external agent |
| Private/local inference | Local OpenAI-compatible Agent Assistant profile |

## 11) Troubleshooting

`auth_failed` for OpenAI:

- Use an OpenAI Platform API key.
- A ChatGPT/Codex subscription is not enough.

`auth_failed` for Claude:

- Use an Anthropic Console API key.
- A Claude Code/Claude.ai login token is not enough.

`auth_failed` for Mistral:

- Use a Mistral La Plateforme API key.
- A Le Chat login token is not enough.

Local model says `model_missing`:

- Run `Discover Models`.
- Pick a returned model.
- Or set `Model override` to an exact model id.

The assistant suggests an OS shell command:

- Do not run it through the Agent Assistant suggestion row.
- Ask again for a GENtle shared-shell command only.

MCP and CLI disagree:

- Treat that as a parity bug.
- Record the exact command/tool call and state path.

## 12) Related docs

- `docs/agent_interface.md` - concise interface overview and schema notes
- `docs/quickstart_claude.md` - Claude-specific internal/external loop notes
- `docs/gui.md` - Agent Assistant UI details
- `docs/cli.md` - CLI and shared-shell command syntax
- `docs/protocol.md` - schemas and adapter contracts
- `docs/gui_cli_mcp_parity.md` - current parity matrix
