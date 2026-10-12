# Agent feedback prompt regressions

The execution-feedback requests in `test_codex_agent_bridge.py` and
`test_pi_agent_bridge.py` are hand-authored synthetic JSON, not real user,
provider, or biological data. Each test constructs its request inline and
calls only the bridge's prompt renderer. No provider credentials or network
access are needed for these regressions.

Recreate and run with:

```bash
python3 -m pytest -q tests/test_codex_agent_bridge.py tests/test_pi_agent_bridge.py
```

The tests verify context visibility and the command-completion/QA/approval
boundary. They do not measure a model's compliance or certify agent reasoning.
