# Provenance-bound GUI tutorial screenshots

These screenshots come from the exact live-green three-tutorial Linux run at
GENtle revision `3de4803aa7cd4f960b01f71af8a7b4b272b7d576`.

These are intentionally revision-pinned historical teaching captures, not
acceptance evidence for today's HEAD. `--check` verifies archive integrity only.
For new runs, the [candidate-bound tutorial gate](../../testing.md#62-candidate-bound-tutorial-gate-for-external-auditors)
derives a verified selection automatically and supports explicit staging into a
new review directory before any repository update. The legacy manual selection
procedure below remains available for historical archives.

For every selected checkpoint, the directory retains:

- an untouched X11-root `*.raw.png` capture;
- a whole-screen `*.orientation.svg` teaching view;
- a padded `*.context.svg` interaction view;
- the original `gentle.tutorial_gui_screenshot_evidence.v1` sidecar; and
- the exact typed semantic snapshot used to resolve the control.

The SVGs are deterministic annotations of the adjacent raw PNG. They are not
independent screenshots and are not scientific verification. The GUI
acceptance ledger's typed comparison with the independently generated oracle
remains the scientific check.

`publication-manifest.json` binds the repository copies by path, size and
SHA-256. The copied sidecars deliberately remain byte-identical to the capture
records, so their absolute paths describe where the original evidence was
created rather than where a later checkout stores the published copy.

To update the images after a GUI change, run the documented three-tutorial
acceptance profile first, review that it is green, update `selection.json` to
the new exact revision, then publish and verify the chosen checkpoints:

```bash
python3 scripts/publish_tutorial_gui_screenshots.py \
  --evidence-root /path/to/green/acceptance-run
python3 scripts/publish_tutorial_gui_screenshots.py --check
```

Do not edit the raw PNG or either SVG by hand. Change the semantic target or
crop policy in the GUI acceptance runner, rerun the tutorial, and republish.
