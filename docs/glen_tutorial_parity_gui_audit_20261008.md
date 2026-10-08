# GENtle tutorial parity and GUI audit

## Candidate and verdict

- Exact source: `30f23aa4cb084694870a139714d0f0fef1726378` (`smoe/main`).
- Review date: 2026-10-08, Linux/X11, software Mesa, Xvfb + Openbox.
- Scope newly added since the prior accepted splicing candidate `a02ef777`: tutorial 08.04, **VKORC1 rs9923231 promoter/luciferase planning**.
- Verdict: the documented deterministic workflow, direct CLI route, shared-shell route and native GUI project loading are functional on the exact candidate. The tutorial does **not** yet have typed `agent_parity` or `gui_acceptance` metadata, so live Agent Assistant quality and semantic click-path parity are not closed by repository automation.

## Source and binary identity

The exact candidate was built with locked offline dependencies. Preserved stripped binaries all self-report `0.1.0-internal.12+git.30f23aa4...`:

| Binary | SHA-256 |
| --- | --- |
| `gentle` | `68a4632eb34ec03cd23c40363798ce3f68d32f7e1d1dff7e3d429705c937b958` |
| `gentle_cli` | `09965f06cdfb012300328c00c3a824aab0e9ff9a6dc5ff09c7b86a9ea2f3155c` |
| `gentle_examples_docs` | `3f0a812cd2bde2420806e312c3de5ca768ab735eaeb91e168cdb934ff3ae44ab` |
| `gentle_mcp` | `5677dea826245f119ea1ab86dd30a1fb61b0ffd7dfe36fb8556206f9673bb7a5` |

The initial four-binary build took 8m49.64s and peaked at 9,412,832 KiB RSS; the matching MCP binary then took 21.09s and peaked at 4,106,400 KiB RSS. These are cold development-build costs, not application latency.

## Documentation and executable parity

Passed on the exact source/binaries:

- `git diff --check` and `cargo fmt --check`.
- tutorial parity-matrix check.
- tutorial catalog: 61 entries valid.
- generated tutorial manifest: 29 chapters valid.
- generated workflow examples: 51 valid.
- tutorial checkout test: 18/18 passed; two pre-existing review-age warnings remain.
- matching-revision CLI/MCP tutorial walkthrough suite: 9/9 passed.
- final splicing regression verification: protocol 10 passed; renderer 9 passed plus one intentionally ignored manual visual smoke; root splicing suite 108 passed serially, including the post-`a02ef777` regressions and ATtRACT windowed-submatrix case.

For the VKORC1 tutorial itself:

- The complete workflow replay succeeded in 32.18s, peak RSS 546,988 KiB, after the pinned human GRCh38/Ensembl 116 cache had been prepared.
- Resulting project: 11 sequences, 9 containers, no arrangements.
- rs9923231 is represented at local 0-based position 3000 with reference `C` and alternates `A,G,T`.
- The chosen transcript is `ENST00000498155` on the negative strand; promoter overlap is true and signed TSS distance is -388 bp.
- Direct CLI and shared shell produced byte-identical normalized promoter-context and reporter-candidate reports and identical explicitly materialized `T` sequence payloads.
- Both routes reject an omitted alternate identically: `MaterializeVariantAllele found multiple alternate alleles 'A,G,T'; choose one explicitly`.
- Complete saved states differ only in ordering of restriction-enzyme groups/sites; DNA and sorted set content agree. This is a serialization-order issue, not a biological difference.

### Agent boundary

The built-in echo adapter exercised request/response transport only. The default Codex adapter failed closed because its search list does not include OpenClaw's packaged Codex binary. Pointing `CODEX_BIN` at that binary launched it, but the isolated adapter lacked a bearer credential and received HTTP 401. Consequently the shared command contract and explicit-alternate guard are verified, but live model-planning quality is **not** accepted. Tutorial 08.04 also lacks an executable `agent_parity` block.

### GUI boundary

The exact `gentle` binary opened the replayed project under a separate Xvfb/Openbox session with software Mesa. The project graph, VKORC1 sequence window, promoter-window features and rs9923231 feature rendered without stderr output. The tutorial metadata is currently `manual_uncovered`: it has neither semantic GUI acceptance targets nor a recorded GUI acceptance SHA. Coordinate-driven screenshots therefore demonstrate rendering, not a durable automated click-path contract.

## White-background captures

The virtual X root and ordinary application chrome were set to white/light. GENtle's lineage graph and sequence-map canvases remain deliberately dark in the product itself. Changing those to white would require a GENtle theme/rendering feature; screenshot tooling should not falsify the application appearance.

Useful evidence:

- `gui-vkorc1/03-sequence-window-maximized.png` — native sequence window with promoter and variant annotations, SHA-256 `1bd02da89760473674fe755b67bfb30786521f084b934e5140babe18810a4a8e`.
- `tutorial-images/vkorc1_rs9923231_promoter_context.png` — deterministic 2400-pixel exported figure on white, SHA-256 `8d7bf10cfb02af7cc4ba07354c5ea74ffe0330f468440ce79522f35c4fd0e61f`.
- reference and alternate reporter-map PNGs are intentionally visually identical at this map scale (same SHA-256 `98965522a6736eb3f23fa5fbcada538396a6f6fbebda6dd0f1c9cf5a5b9fda6d`); the one-base allele difference remains in the sequence/project data.

## Scrolling and resizing

Native smoke:

- six requested X11 geometries from 820x520 through 1920x1080 were acknowledged exactly in 63-66 ms each; this measures window-manager acknowledgement, not repaint completion.
- 100 wheel events, 10 ms apart, dispatched in 1.06s; the application remained responsive and the before/after root screenshots differed by 450,668 pixels, demonstrating an actual viewport change. This is a bounded input smoke, not a frame-latency distribution.

The maintained Criterion target measures the real CPU-side DNA-window constructor, deferred hydration, first and steady embedded egui frames at 820x520, 1200x800, 1600x1000 and 1920x1080, plus the first frame after resizing from 1200x800. It explicitly does not measure native compositor, GPU, event delivery or perceived interaction stalls.

The exact `bench-audit` build took 30m09s (30m18s including execution), peaked at 21,135,284 KiB RSS and made no swaps. Reusing that exact executable, the complete `--quick --noplot` series including the opt-in PATZ1 public stress fixture ran in 15.23s with 157,584 KiB peak RSS. The generated PATZ1 state and report hashes are `42c216e44f855dd9bc66923836145666465b1edf8429131624d9b526f34799fd` and `2ea1c00df161d8174945cc16d88971909b367a875259807819c8d3af11df7f5a`. Criterion midpoint estimates:

| Fixture | First frame: 820 / 1200 / 1600 / 1920 | Steady frame: 820 / 1200 / 1600 / 1920 | First frame after resize: 1200->820 / 1600 / 1920 |
| --- | --- | --- | --- |
| synthetic 120 kb, 0 features | 15.881 / 15.514 / 15.462 / 16.652 ms | 0.423 / 0.425 / 0.425 / 0.425 ms | 1.811 / 1.274 / 1.303 ms |
| TP73 83,686 bp, 140 features | 50.001 / 50.370 / 59.002 / 56.755 ms | 1.894 / 1.893 / 1.892 / 1.900 ms | 7.021 / 6.487 / 6.576 ms |
| PATZ1 20,802 bp, 75 features, 13 transcripts, 17 joined source records | 26.909 / 27.412 / 27.244 / 27.360 ms | 0.859 / 0.838 / 0.805 / 0.786 ms | 4.064 / 4.005 / 4.054 ms |

These are one short audit series, not a historical regression verdict. The TP73 first frame crosses 50 ms at all sizes, so first-open latency remains materially heavier than steady painting; resizing itself stays below 7.1 ms in this CPU-side proxy. Window size is not the dominant factor for these fixtures.

## Actionable gaps

1. Add typed `gui_acceptance` metadata and stable semantic widget IDs for the VKORC1 promoter/variant controls; coordinate scripts are not a maintainable acceptance surface.
2. Add an `agent_parity` case that requires explicit allele `T`, validates the proposed shared-shell command, and proves that ambiguous A/G/T input is never silently resolved.
3. Extend the maintained GUI benchmark with a deterministic scroll-event + next-frame series. The current target covers sizes/resizes but no scroll input.
4. Decide whether lineage/sequence canvases should gain a printable light theme. Until then, retain honest dark canvases and use white exported figures for tutorials intended for print.
5. Consider canonical ordering for restriction-enzyme groups/sites if byte-stable whole-state equivalence is desired across direct CLI and shared shell.

## Side effects and limits

- No GENtle source, branch, tag or remote was changed.
- The dirty primary checkout was not cleaned or switched. Its pinned reference cache was populated/rebuilt by the documented workflow (large downloaded/indexed GRCh38/Ensembl 116 assets); source changes were untouched.
- Rebuildable Cargo output was cleaned only from the audit-owned target after preserving source-identified binaries; the exact optimized benchmark output and Criterion records were then recreated and retained.
- Native GUI rendering is covered; human scientific approval, release approval and live Agent Assistant quality remain separate decisions.
