# Feature Editor GUI evidence

This directory retains publication-safe manual/hybrid evidence for tutorial
02.05. The capture used only the committed synthetic 120 bp GenBank fixture and
an exact-revision native Linux GUI build.

- `01-location-preview.raw.png` shows the exact 21..40 to 21..45 Location
  preview, its feature fingerprint and review-only related annotations.
- `02-create-overlap-preview.raw.png` shows the `review_patch` Create preview,
  ordered qualifiers and two overlap/shared-identifier review candidates.

No mutation was applied while making these captures. They support the two
preview claims only; the CLI regression remains the oracle for the complete
five-operation state transition, and complete live GUI Apply/Undo/Redo/Split/
Merge/Delete acceptance remains unclaimed.

See [`evidence.json`](evidence.json) for hashes and scope.
