# Project and Sequence State

<!-- manual-counterparts:start -->
Related GUI reference: [Engine Ops State Persistence](../gui/engine-ops-state-persistence.md) | [Loading and creating sequences](../gui/loading-and-creating-sequences.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

`state-summary` output includes:

- sequences
- containers
  - each container row now includes `declared_contents_exclusive`
- display visibility flags

Project aliases:

- `save-project PATH` aliases `export-state PATH`
- `load-project PATH` aliases `import-state PATH`

### Unavailable Candidate Sidecars

A permissive load can open unaffected sequences when a candidate sidecar is
missing or malformed. Its original reference and structured
`candidate_sets_load_warning` remain available in project metadata; this is not
an empty candidate collection. Candidate commands and project Save/Save As,
including implicit CLI/MCP saves, refuse with `candidate_store_unavailable`.
The saver creates no staging/output directories and does not replace the project
or remove its sidecars on this refusal. A command may already have changed its
transient engine or written a separate approved artifact; this is not rollback.
One-shot CLI/MCP errors explain that transient changes are discarded when the
invocation ends, while separately written files remain.

Restore the original dependency or intact project/sidecar backup and reopen.
Old warning-only states cannot recover a lost reference by guessing a filename.
There is no automatic repair/discard command. `state-summary` and unaffected
inspection remain available. Strict-load opt-in retains its existing hydration
errors; unknown-schema and old warning-only projects remain loadable, but
candidate use and saving stay blocked. No new count validation is implied.


- Shared shell:
  - Commands:
    - `containers set-exclusive CONTAINER_ID true|false`
    - `load-project PATH`
    - `save-project PATH`
    - `sequence create --sequence-text DNA [--output-id ID] [--name TEXT] [--topology linear|circular]`
      - creates a persistent project sequence from inline bases so an
        Agent Assistant or external agent can promote a paper/vendor/database
        sequence into the same state used by GUI, GUI Shell, and later commands

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Error Behavior](errors.md) | [Next: Collections](collections.md)
<!-- manual-page-nav:end -->
