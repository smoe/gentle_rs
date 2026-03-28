# Contributing to GENtle

Thank you for your interest in contributing to GENtle.

For installation and first local runs, start with [INSTALL.md](INSTALL.md).

## Start Here

Depending on the area you want to change, these documents are the most useful
entry points:

- `docs/architecture.md`
- `docs/roadmap.md`
- `docs/protocol.md`
- `docs/gui.md`
- `docs/cli.md`
- `docs/testing.md`

## Minimal Validation Loop

Before opening a pull request, a useful baseline is:

```sh
cargo check -q
cargo test -q workflow_examples -- --test-threads=1
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
```

Disk-usage note:

- This repository configures a shared Cargo target directory in
  `.cargo/config.toml` (`../../.gentle_target_shared`) so sibling worktrees
  reuse compiled dependencies instead of each creating a large private
  `target/`.
- You can override this locally with `CARGO_TARGET_DIR` if your environment
  needs a different path.

## Contribution Notes

- Keep behavior aligned across GUI, CLI, and other interfaces when possible.
- Update documentation when user-visible behavior changes.
- Prefer deterministic checks and reproducible examples when validating a
  change.

## Agentic Engineering

It is perfectly fine to use agentic engineering to contribute to GENtle.

Examples of acceptable tools include:

- OpenAI Codex
- Anthropic's Claude
- OpenCode.ai
- kilo.ai

Expectations:

- You remain responsible for reviewing the generated changes.
- Please run deterministic validation before opening a pull request.
- Keep prompts and agent instructions focused on reproducible engineering work.

## Pull Request Expectations

Every pull request should include:

- a short summary of the change
- the validation commands you ran
- the prompts used to describe or direct the development

If no agentic tooling was used, say `Not used` in the prompt section.
