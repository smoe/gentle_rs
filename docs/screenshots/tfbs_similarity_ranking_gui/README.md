# TFBS similarity-ranking screenshots

These native GENtle captures were recorded on 2026-10-01 with source revision
`3c6226337db450313f3b669007dc90e133a96122`, after loading the committed 46 bp
FASTA into an explicit temporary state. They show the reviewed motif selection,
ranking settings, and whole-sequence action immediately before execution.

The action did **not** produce a result-table capture: selecting it terminated
the exact reviewed GUI build with `thread '<unknown>' has overflowed its stack`.
The offline CLI workflow and ClawBio wrapper completed independently, so these
images are evidence of GUI setup and the reproducible boundary, not of GUI
ranking acceptance.

| File | SHA-256 |
| --- | --- |
| `01-motif-selection.png` | `4e73ca6d1a8ce8810963b2a099dfa84fa4d81af1d504f1ccb54ca510b8794b62` |
| `02-ranking-settings.png` | `ee85822ccbcc0bc92e1b677ac73e3180c3603eaee36a1588d301b108f6b66fa9` |
| `03-whole-sequence-action.png` | `a2fba61cbd11ba50f6a43d346b4563f43558adf1a73105766b03128891b4230d` |
