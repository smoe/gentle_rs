## Superseded Tutorial

The former TP73 promoter/luciferase walkthrough is no longer the canonical
GUI-first reporter tutorial.

Use this page instead:

- [`docs/tutorial/vkorc1_warfarin_promoter_luciferase_gui.md`](./vkorc1_warfarin_promoter_luciferase_gui.md)

Why:

- the new tutorial starts from a concrete ClawBio-style pharmacogenomic alert
  (`rs9923231`)
- it uses the reverse-strand `VKORC1` promoter geometry explicitly
- it replaces the older "promoter to luciferase" placeholder with one
  biologically better-defined TSS-centered reporter fragment

The TP73 material can still remain as a secondary example in workflow history,
but the maintained GUI-first path is now the VKORC1 / warfarin reporter
tutorial above.

1. The same promoter candidate can be regenerated from the same anchor/extract
   parameters.
2. The vector import ID is stable and traceable.
3. Primer/qPCR reports refer to the same construct ID.
4. GUI actions and CLI commands do the same thing at each step.

## Reference Workflow JSON (Appendix)

Reference example file:

- `docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json`

Associated metadata in that file:

- format tag: `gentle.workflow_example.v1`
- `test_mode: "skip"`

Meaning:

- The JSON is a machine-readable workflow example used for reproducible
  GUI/CLI/agent consistency checks.
- It is maintained in-repository as reference documentation for this
  tutorial scenario.
- `test_mode: "skip"` indicates external dependencies are intentionally involved
  (for example online GenBank fetch and prepared-genome assumptions), so default automated tests only
  check the JSON structure for this example instead of always executing it
  end-to-end.
