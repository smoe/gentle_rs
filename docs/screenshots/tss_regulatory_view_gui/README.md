# Native TSS / Regulatory-view tutorial evidence

This directory retains publication-safe manual/hybrid evidence for tutorial
08.16. The run used only the committed synthetic-minus GenBank record and its
hash-bound profile report.

The two tutorial images are untouched captures of the native DNA viewer:

- `01-annotated-context.raw.png` shows the annotation-only TSS view before a
  profile is attached.
- `02-attached-profile.raw.png` shows the lower part of the same viewer after
  the report was attached. The status line identifies the report producer and
  the bottom lane is the report-provided `MA0004.1` score trace. This deliberately
  tiny teaching sequence has zero display values after negative scores are
  clipped; the lane is therefore flat, not missing.

`agent-assistant-report-command.raw.png` records the intervening hosted action.
GENtle recognized the literal `ui open tss-view --report ...` prompt as a
reviewed shared command, kept the explicitly launched DNA viewer as its subject,
and queued the report for that viewer. No external language-model transport was
invoked for this direct command.

The adjacent semantic snapshots retain the exact native-window geometry used
for the tutorial captures. [`evidence.json`](evidence.json) binds the display
revision, producer revision, binary, inputs, command checkpoint, screenshots,
and snapshots by SHA-256. This is live manual/hybrid evidence, not an automated
GUI-acceptance verdict.
