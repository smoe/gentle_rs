# macOS egui Upstream Report Template

Use this when the bug reproduces in `cargo run --bin gentle_egui_window_repro`
and looks upstream rather than GENtle-specific.

## Environment

- macOS version:
- machine / CPU:
- Rust toolchain:
- `eframe` version:
- `egui` version:
- renderer backend:

## Repro command

```bash
cargo run --bin gentle_egui_window_repro
```

## Repro mode

- `immediate native viewport`
- `deferred native viewport`
- `embedded egui::Window`

## Steps

1.
2.
3.

## Expected behavior

-

## Actual behavior

-

## Scope check

- Reproduces in minimal repro: yes / no
- Reproduces in embedded mode too: yes / no
- Reproduces in full GENtle with hosted fallback: yes / no
- Reproduces in full GENtle with `GENTLE_MACOS_NATIVE_CHILD_VIEWPORTS=1`: yes / no

## Notes

- stale twin windows:
- wrong focus order:
- maximize / restore specific:
- resize specific:
- close / reopen specific:
