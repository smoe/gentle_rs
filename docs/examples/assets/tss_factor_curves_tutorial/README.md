# TSS factor-curve parser templates

Origin: hand-authored copies of the request JSON blocks in sections 2 and 3
of `docs/tutorial/08-17_tp73_dnp73_factor_curves.md`. No new biological data,
production output, or approval has been added.

Recreation: copy those two JSON blocks into `inventory.template.json` and
`materialize.template.json`, retaining both `PASTE_FROM_YOUR_PREVIEW`
placeholders exactly. Object key order, whitespace and LF/CRLF do not matter.

Use: `agent_parity.parser_payload` in source 08.17 binds each learner-created
`@file` to its template. The shared tutorial validator compares semantic JSON
with the guide and parses it without executing a command. These are not runnable
materialization fixtures: real execution still requires the learner's current
preview digest and selected TSS IDs. A successful parser check is not approval.
