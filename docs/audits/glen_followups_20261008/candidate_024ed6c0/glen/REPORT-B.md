# Help/startup comparison B — build profile

## Boundary

- Exact source: `024ed6c0520d38c4221517219d44c792084de4e5`
  (tree `3982cc6059530c9caa21ba61e2fd63f4405bf0cc`).
- Same Linux host, Xvfb/Openbox path, source tree, fonts, documentation assets,
  public Gibson project and fixed interaction sequence for both profiles.
- `dev` binary SHA-256:
  `635620123bb64a492869eab82070031bf96a53fb77197a22cd9d253f96242a7a`.
- `package-opt1` binary SHA-256:
  `ab71e0c1c96e56f81276de1e37ec8fee6f9928d4d71da9265fad64240b0bd541`.
- Both binaries self-report source revision `024ed6c0520d38c4221517219d44c792084de4e5`.
- The package profile's captured final `rustc` command has `opt-level=1`,
  256 codegen units, debug assertions enabled and no LTO. The ordinary dev
  build used default features and no profile override.
- Per profile: three cold-cache and three warm-cache sessions. One explicit
  warm-cache filling session per profile is retained but excluded by the
  analyzer through `metadata.excluded=true`.
- All 12 included sessions completed the same interaction sequence with 52
  retained events, zero dropped events, zero dropped image observations and
  zero failed spans. Inclusive trace spans are reported independently and are
  never summed.

## Results

Values are median `[min–max]`, milliseconds, `n=3` per cell.

| Measurement | dev cold | package-opt1 cold | dev warm | package-opt1 warm |
| --- | ---: | ---: | ---: | ---: |
| launch → usable main window | 2646.717 `[2642.747–2701.962]` | 2389.852 `[2384.833–2396.339]` | 2684.394 `[2671.679–2685.118]` | 2414.108 `[2363.981–2440.165]` |
| launch → subject window | 3730.752 `[3688.287–3861.372]` | 2798.369 `[2796.648–2800.357]` | 3772.771 `[3759.116–3826.339]` | 2827.483 `[2767.771–2949.935]` |
| `help_preparation` | 523.850 `[522.606–525.740]` | 83.479 `[82.577–83.490]` | 53.805 `[53.592–54.720]` | 12.198 `[12.193–12.321]` |
| `help_manuals` | 43.349 `[42.605–43.389]` | 9.396 `[9.288–9.588]` | 42.971 `[42.885–43.885]` | 9.347 `[9.347–9.446]` |
| `help_shell_reference` | 5.912 `[5.833–6.016]` | 1.433 `[1.431–1.453]` | 5.919 `[5.847–5.926]` | 1.455 `[1.455–1.464]` |
| `help_tutorial_discovery` | 3.586 `[3.521–3.798]` | 1.123 `[1.109–1.127]` | 3.576 `[3.553–3.576]` | 1.116 `[1.113–1.129]` |
| `help_tutorial_selected_load` | 471.883 `[469.747–472.518]` | 71.512 `[70.431–71.651]` | 1.324 `[1.303–1.331]` | 0.279 `[0.278–0.280]` |
| first Help click → window | 189.765 `[187.811–303.486]` | 119.603 `[118.237–119.795]` | 120.468 `[117.093–178.091]` | 120.499 `[117.924–120.559]` |
| Tutorial click → window | 5317.236 `[5247.781–5385.015]` | 250.956 `[190.421–256.593]` | 5177.171 `[5126.299–5182.604]` | 253.110 `[193.993–253.117]` |
| Topic switch click → window | 454.099 `[394.514–456.452]` | 127.557 `[127.393–128.927]` | 451.941 `[381.271–464.279]` | 122.741 `[120.907–130.395]` |
| tutorial menu hover → observation | 1002.491 `[1002.410–1002.563]` | 1003.461 `[1002.446–1003.794]` | 1002.444 `[1002.387–1002.503]` | 1002.420 `[1002.407–1002.496]` |

Cold sessions generated one PNG cache entry in both profiles. Warm sessions
started and ended with that entry. Both cold profiles recorded one completed
rasterization and no failure; warm profiles recorded none. The dev cold
rasterization median was 470.456 ms versus 71.154 ms for `package-opt1`.

## Interpretation

`package-opt1` improves the median usable-window time by 256.865 ms cold
(9.7%) and 270.286 ms warm (10.1%). The first subject-correct window improves
by about 0.93–0.95 s (25%). Help preparation and topic switching are also
materially faster. Most conspicuously, the externally observed Tutorial-window
latency falls from about 5.2 s to about 0.25 s in both cache conditions.

The profile benefit is therefore not merely a shift of startup work into the
first Help use: it remains visible at the subject and Help interactions. The
fixed one-second menu-hover observation is harness dwell time and appropriately
does not improve. The measurements do not justify new shared-font, failed-SVG
or menu-cache work: no font-specific failure was observed, SVG failures were
zero, and the menu discovery span itself was sub-millisecond to low-millisecond.

This virtual-X11 comparison is performance evidence only. It does not establish
physical desktop responsiveness, package acceptance or release acceptance.

## Retained evidence

- Raw runs: `help-startup-B/runs/`
- Per-run commands: `help-startup-B/*.command`
- Included-run table: `help-startup-B/analysis-B/runs.csv`
- Machine-readable statistics: `help-startup-B/analysis-B/statistics.json`
- Exclusion-aware analyzer: `analyze_help_runs.py`
- Exact binaries/build logs: `builds/dev/`, `builds/package-opt1/`
- Host/toolchain/profile records: `host-*.txt`, `profile-definitions.txt`,
  `cargo-metadata-no-deps.json`
