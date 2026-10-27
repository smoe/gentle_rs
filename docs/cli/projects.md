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
