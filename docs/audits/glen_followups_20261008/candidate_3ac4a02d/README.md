# Repaired Candidate Evidence

Runtime/workflow source: `3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da`.
Lockfile: `3dfe44c7a08bf779e32289cf1312103ea1fb5883e44bba5f6e1912ce3ed77745`.
These are new receipts, not replacements for historical evidence. Local status
commits do not change the frozen remote candidate or inherit runtime acceptance.

## Scroll CPU Cases

`scroll_cpu_smoke.json` is the unchanged receipt from artifact `11581332454` of
[CI run 37848423833](https://github.com/smoe/gentle_rs/actions/runs/37848423833).
Its raw SHA-256 is
`87a274adaf3c8046f4458b8d0374757ea7db3cc63fab2a5eed7e9dbc5094c00d`.
Origin: actual Ubuntu 24.04/x86-64 dev-profile CI execution of the scoped
`gui_operations --test` benchmark, not fabricated timing or model output.
The fixture origins/recreation remain documented in the benchmark source and
PATZ1 preparation script; synthetic control, TP73 and public PATZ1 are separate.

On 2026-10-09 (Europe/Berlin), independently verify all 25 retained file hashes,
including executable
`daacf5896aa076df63746fafd90bd03f48efb5a83568678f83d8b9d468e2ff39`.
Verify all 24 observations against raw `smoke.log`: three fixtures, four screen
sizes and both point-wheel directions; exact source/version, bounded viewport,
unchanged span/length and one unchanged content hash per fixture. The tested
source also checks the no-drift next frame and full content equality.
`scroll_verification.json` retains the independent verification result; digest
`58a16df1a020904e7b7e7ea4573e1eb314f9aec48c253094cd6bf70281fbb3ab`.
Its verifier at `/private/tmp/verify_gentle_scroll_3ac4a02d.py` has SHA-256
`04b1f56b53c58adf0c53a66d784d278f21aa3dc2ebe155175afb5bbbf5c9acbf`.
No timed regression, native input-to-paint latency or package-profile claim.

Retrieve only the named artifact, then verify this receipt and every listed
hash before inspection:

```bash
gh run download 37848423833 -R smoe/gentle_rs \
  -n gui-scroll-cpu-smoke-3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da \
  -D DIR
```
The retained local bundle is `/private/tmp/gentle-scroll-3ac4a02d`.

## CLI Build Boundary

`agent_cli_build.json` is the unchanged macOS/arm64 dev build receipt from
artifact `11580808183` of the same run. Its raw SHA-256 is
`8f28d2f10887acd283234d40a5440a27d55fb17e88d174bd57c3b3e39986c341`.
Independently verify both listed artifacts and execute the downloaded CLI's
`--version`; its complete output exactly matches `binary_version`, including
the full frozen source SHA. Binary digest:
`a6ea8462b52dcb45463d419463460dab1c8fbbdd5d27f3ff102b2f50f5d253fb`.
Actual observed-version digest:
`10d1c9dd5aa154ee81b1243006b7f7f6e127d5f60d8b7f1189bab1e31e6a0550`.
Local bundle: `/private/tmp/gentle-agent-cli-3ac4a02d`.

This build receipt verifies a source-bound dev CLI, not a model invocation.
Its `live_agent_accepted:false` and `package_accepted:false` remain unchanged.
Actual fresh public/model acceptance is recorded separately below. No local
Rust build or credential inspection was performed. No public-locus GUI,
scientific or package approval.

## Synthetic Native Contract

`native_gui.json` is the unchanged receipt from artifact `11582285369` of
the same CI run; SHA-256
`bef0707566d9f68f19f51d57cf760d47c8fbef30818b2f05dcc22bf4c75a121d`.
Actual Linux/X11 dev execution passed 19 focused Rust tests, scoped generation
checks and eight ordinary native input steps in an offline network namespace.
The independent check verifies all 115 hashes and safe exact file membership,
the unsatisfied-to-satisfied completion transition, explicit T, all three full
oracle `Seq` records and sole position-6 C/T change in the synthetic 20-bp pair.
It does not sort qualifier arrays or compare DNA alone.

`native_verification.json` retains this independent result; SHA-256
`cf6b2500740d97071f7e010f53b3fdf97d6de68f1880d9dc3b1e0b4559c4e6f0`.
Verifier `/private/tmp/verify_gentle_native_3ac4a02d.py` has SHA-256
`f614b2c87099243cbab6074e84cb9b2f5361c11b6a8cd0753126549cd5c4a9ec`.
Local raw bundle: `/private/tmp/gentle-vkorc1-gui-3ac4a02d`.
Retrieve the named `vkorc1-gui-audit-3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da`
artifact from run `37848423833`; verify every listed hash before inspection.

On 2026-10-09, visually inspect unchanged `open_fragment`, `review_t`,
`scroll_to_pair` and `materialize_pair` raw PNGs: light application chrome,
dark DNA canvas, visible T and reachable pair action. Layered windows after
materialization are retained, not cosmetically repaired. These are synthetic
native control captures, not public-locus GUI, base-pair image proof, clinical,
wet-lab, timed responsiveness, package or human scientific approval. The
separate public projection is verified and inspected below.

## Public Reference And Base Projection

`public.json` is the unchanged receipt from artifact `11583615420` of the
same run; SHA-256
`4ae522db4eac214ba7877c775eb1698750424f6932c4710dbb351c5eaf9bb9b9`.
It retains all 44 hashes, full GRCh38/Ensembl 116 reference-input hashes and
public refSNP response origin. Independently verify all 44 files and all four
full ordered `Seq` records against both CLI/shared-shell routes. The large
reference cache is not uploaded or downloaded locally.

`base_comparison.json` SHA-256:
`5ced393927f191ccd78bb3f81a26597688a1286eaaefb2b5ee7ce30d92dc4313`.
`base_comparison.svg` SHA-256:
`9aee9045ce6e107fd43c0c3f68946907f4e8ac7aaed275367507ca2b528b2c9d`.
Recompute every proof field from `direct.project.gentle.json` with
`scripts.prepare_public_vkorc1_audit.pair_content`; reproduce the exact SVG
with `public_pair_svg`, full revision and receipt binary hash. Both checks
pass. Render the unchanged SVG with `rsvg-convert` and visually inspect the
33-bp window [572,605): highlighted position 588 is C/C/T. JSON retains the
exact full 1,089-bp inserts and case; uppercase is display-only. This labelled
data projection is not a native screenshot or reporter-construct map.

Local bundle: `/private/tmp/gentle-vkorc1-public-3ac4a02d`. Retrieve the named
`vkorc1-public-audit-3ac4a02dd81a7b2dabe3c2bd82b28d76f6ce21da` artifact from
run `37848423833`; verify every listed hash before use. GENtle uses these
receipts to audit public retrieval/extraction, allele choice and projection,
not as evidence of biological function, clinical effect or wet-lab suitability.

## Genuine Live Agent

`live_agent.json` is a new executed receipt, SHA-256
`89f4fdba74d027d813fd5daecb9a15527cb6b9a9d3fb04b23110567c441c79d7`.
On macOS/arm64, use the verified CI-built dev CLI and exact fresh public
starter. Call `agents ask codex_local_stdio --catalog assets/agent_systems.json
--timeout-secs 180 --max-retries 0` twice with the prompts retained in the
receipt, using the existing authorized Codex CLI login. No credentials were
read or copied, and no particular model ID is attested. The bridge uses an
empty read-only sandbox; no Claude invocation or agent suggestion auto-runs.

The actual first response asks which of A/G/T to choose and suggests nothing.
The second returns exactly two explicit genomic-forward T/reference `ask`
commands, both `not_run`. Raw `cmp` after each request passes. The bare-
alternate negative control exits one with `InvalidInput` and preserves raw
project bytes. Review and manually execute the two explicit suggestions;
both exit zero. All four full records equal both fresh CI routes, including
ordered feature/qualifier provenance, and both source records are unchanged.
The 1,089-bp pair differs only at position 588, C/T; nonvariant bytes/case stay
unchanged. No DNA-only or sorted-array parity shortcut is used.

Raw directory: `/private/tmp/gentle-agent-public-fresh-3ac4a02d`; its retained
`verify_live_acceptance.py` SHA-256 is
`f1337f8d249b172aec45ca03173265dd8c71b3fbd3eca466fa3ea908877f1953`.
Recreate with the receipt's prompts from a fresh copied public starter, review
suggestions before manual execution and rerun full-record verification. Model
response bytes are not deterministic; protocol/effect acceptance is checked.
The earlier starter's failed provenance comparison remains historical evidence,
not a hidden or inherited pass. This is not public native GUI, package, clinical,
functional, wet-lab, performance or human scientific acceptance.

## Checkout Contract

GENtle uses these receipts as the six-item acceptance ledger's evidence, not
runtime scientific inputs. Nine scoped LF attributes and the fast checkout
regression preserve their raw digests in LF/CRLF clones, with separate failing
negative controls when each attribute is removed. Do not normalize or regenerate
the receipts to conceal a checkout difference.
