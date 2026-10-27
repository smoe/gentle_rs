# CLI Options and Progress

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

## `gentle_cli` (automation/agent CLI)

`gentle_cli` exposes operation-based control with JSON input and output.
It is intended for scripted workflows and AI agents.

State is persisted to `.gentle_state.json` by default (override with `--state PATH`).


Global CLI options:

- `--state PATH`: use a non-default project state file
- `--progress` or `--progress-stderr`: print live progress events to `stderr`
- `--progress-stdout`: print live progress events to `stdout`
- `--allow-screenshots`: rejected; the CLI screenshot bridge remains disabled,
  independently of [user-approved GUI screenshot help](../gui.md#screenshot-capture-policy)

Place global progress options before the command name, for example
`gentle_cli --state "$STATE" --progress-stderr rna-reads interpret ...`.
Current progress events include TFBS annotation updates, genome-prepare
updates (download/index phases), genome-track import updates, RNA-read
interpretation/alignment updates, and primer/qPCR design updates.
Primer/qPCR progress lines expose backend choice/fallback plus candidate and
search counts, for example `progress primers ... stage=pair_search ...` and
`progress primers ... stage=assay_search_complete ...`.
Progress-capable Primer3 builds additionally emit
`stage=primer3_search primer3_record=N primer3_work=COMPLETED/BOUND
primer3_percent=P`. Transcript-assay runs also attach the GENtle bounded
search-record id and ordinal, deterministic record/operation pair-work bounds,
elapsed time, and a separately labelled uncertain linear runtime projection.
The native work values count bounded candidate evaluations; they are not
seconds and must not themselves be presented as an ETA. GENtle runs older
Primer3 binaries in legacy mode when `--help` omits `--progress`, so those
installations continue to work but expose only the coarse Primer3 stages.
When `--progress-stdout` is used, progress lines are emitted before the final JSON output.

When the configured Primer3 is a bare name such as the default `primer3_core`,
GENtle considers every match on `PATH` rather than taking the first. It selects
a build advertising `--progress` over one that does not, because bounded search
reporting and SIGUSR1 status depend on it; then the highest version, compared
component by component so `2.10.0` outranks `2.6.1`; then the most recently
created file, falling back to the modification time where the filesystem
records no birth time. A configured value naming a path is authoritative and is
used unchanged, so an operator's pin is never silently overridden. `primers
preflight` reports the candidates considered, each one's `--progress` support,
version and creation time, and a `selection_reason` naming the winner.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: GUI Launcher](gui-launcher.md) | [Next: File Formats](file-formats.md)
<!-- manual-page-nav:end -->
