# macOS egui Window Investigation Pack

This folder keeps the macOS `egui`/`eframe` child-viewport investigation
separate from GENtle's main biology and GUI feature work.

## Goals

- isolate native child-viewport behavior from GENtle-specific hosted-window
  workarounds
- keep the canonical repro command, manual test matrix, and upstream-report
  notes in one place
- make it easier to decide whether a window bug belongs in GENtle or upstream
  `egui`/`eframe`

## Entry Points

- Repro binary: `cargo run --bin gentle_egui_window_repro`
- Repro source: `src/bin/gentle_egui_window_repro/main.rs`
- Hosted-window wrapper under investigation contrast:
  `src/egui_compat.rs`
- Main-app comparison path:
  `GENTLE_MACOS_NATIVE_CHILD_VIEWPORTS=1 cargo run --bin gentle`

## Files In This Pack

- `README.md`: overview and command entry points
- `manual_matrix.md`: the canonical manual scenario list
- `upstream_report_template.md`: concise capture template for upstream issue
  filing or local triage notes

## Suggested Workflow

1. Run `cargo run --bin gentle_egui_window_repro`.
2. Exercise the scenarios in `manual_matrix.md`.
3. Compare the repro outcome with GENtle's hosted-window behavior and, when
   relevant, with `GENTLE_MACOS_NATIVE_CHILD_VIEWPORTS=1 cargo run --bin gentle`.
4. If the failure reproduces in the minimal harness, capture it with
   `upstream_report_template.md` before filing or updating an upstream issue.

## Guardrails

- Keep this pack focused on viewport lifecycle, focus ordering, resize,
  maximize, and stale-window behavior.
- Do not add GENtle biology or engine logic here.
- Do not rely on screenshot capture unless that is explicitly approved by
  project policy for the investigation session.
