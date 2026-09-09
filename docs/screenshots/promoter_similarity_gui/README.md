# Promoter-similarity GUI tutorial screenshots

These screenshots come from the live-green, isolated Linux/X11 GUI tutorial at
GENtle revision `b087cb22ca54f0c92b8e47796b6772a7cb7a04ae`.

The four retained checkpoints show how to open the saved upstream-region
manager, launch `Promoter similarity...`, inspect the target/index universe,
and review the alignment and promoter-window constraints before starting the
search.

For every checkpoint, this directory retains the untouched X11-root PNG, the
whole-screen orientation SVG, the focused context SVG, the original screenshot
evidence sidecar, and the exact semantic snapshot. The SVGs are deterministic
teaching projections of the adjacent raw PNG; they are not independent
screenshots or scientific evidence.

`publication-manifest.json` binds every repository copy by path, size, and
SHA-256. To refresh the images after a GUI change, rerun the typed chapter from
the exact clean revision, review the raw captures, update `selection.json`, and
run:

```bash
python3 scripts/publish_tutorial_gui_screenshots.py \
  --evidence-root /path/to/green/acceptance-run \
  --selection docs/screenshots/promoter_similarity_gui/selection.json \
  --output-root docs/screenshots/promoter_similarity_gui

python3 scripts/publish_tutorial_gui_screenshots.py \
  --output-root docs/screenshots/promoter_similarity_gui \
  --check
```

Do not edit the PNG or SVG files by hand. Change the semantic target or crop
policy, rerun the GUI tutorial, and republish.
