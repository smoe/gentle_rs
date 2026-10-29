## Engine Ops State Persistence

<!-- manual-counterparts:start -->
Related CLI reference: [Project and Sequence State](../cli/projects.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Engine Ops panel input state is persisted in project metadata per active
sequence id. This includes panel visibility and operation-form inputs.

The same persistence payload now also stores:

- shell panel visibility
- last shell command text

Metadata key format:

- `gui.engine_ops.<seq_id>`

### Saving After A Candidate-Sidecar Warning

If a permissive project load reports an unavailable candidate store, you can
still inspect unaffected DNA. Candidates are **unavailable, not an empty list**;
the candidate panel displays the refusal instead of offering replacement data.
Save, Save As and shared implicit saves return `candidate_store_unavailable`
before save-owned files/directories change. A failed save must not be treated
as a saved project. This does not roll back earlier edits or separate exports.

Keep unsaved work separate. Restore the original sidecar or an intact
project/sidecar backup, then reopen the original project; restoring files does
not unlock the already degraded session. Historical warning-only files need a
backup containing the lost reference. No automatic repair or discard is offered.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Command Palette and History/Jobs Panels](command-palette-and-historyjobs-panels.md) | [Next: Loading and creating sequences](loading-and-creating-sequences.md)
<!-- manual-page-nav:end -->
