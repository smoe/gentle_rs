## Protein Expression Handoff

<!-- manual-counterparts:start -->
Related CLI reference: [Host Profiles](../cli/hosts.md) | [Experiment Planning](../cli/planning.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

`Services > Protein Expression Handoff...`, the command palette, and
`ui open|focus|close protein-expression-handoff` address one read-only inspector.
The intent takes no arguments. Product ID, optional pasted reviewed-requirements
JSON and optional provider ID are session-local inputs; none is selected from
the active viewer or a provider template automatically. A blank product permits
inspection of unanswered product-definition questions.

**Inspect** runs the existing `planning protein-expression-handoff` command on
a history-free background snapshot with the effective project objective/profile.
The window displays the shared v2 readiness/review gate, unanswered requirements,
sequence and annotation evidence, CDS/translation and tags, unchanged chassis/
vector candidates, neutral draft, selected-provider projection and typed
withholding reasons. **Copy report JSON/text** retains the shared report;
suggested actions are copy-only. There is no submit, quote-export or order action
in this inspector. A requested yield or local projection is not predicted yield,
expression readiness or provider acceptance.

Changing inputs clears the report. Project edits or reopening/replacing a
project discard pending/displayed results; viewport-only changes may also require
reinspection. Closing abandons the result receiver without committing worker
state. Reopening/focusing preserves local inputs and does not choose a provider.
No worker commits, project edits, output files or provider calls occur. External
catalog changes require an explicit fresh inspection; the window does not poll
providers. Tutorial [09.03](../tutorial/09-03_protein_expression_handoff_cli.md#gui-inspection-of-the-same-handoff)
provides synthetic neutral/withheld/projected checkpoints. Native GUI acceptance
remains separate from source and unit checks.

The named MCP `protein_expression_handoff` tool and inner-agent suggestions
reach the same builder. MCP compares saved `state_path` only, so save GUI edits
before comparing reports. Native tutorial screenshot capture remains pending;
the inspector exposes window/control/status semantic IDs under
`gui-test-support` for locating the actual rendered controls without injecting
input or disclosing pasted requirements in the semantic snapshot.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Reference Genomes (Main Menu)](reference-genomes-main-menu.md) | [Next: Patterns And gRNA Setup](patterns-and-grna-setup.md)
<!-- manual-page-nav:end -->
