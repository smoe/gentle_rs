# FLNA Splicing Readability Handoff

## Scope And Review

Working base: `b1c47ed6d9f6bd9c29e49096632c20528ef2a697`, branch
`gentle_rs_2_main`. Fetched `origin/main`:
`04debbe0eb229e9758090d1964663849621470b9`. The local transcript-CDS repair
is preserved, not reverted. The checks below cover the readability patch on
that base; record the final committed/rebased SHA and rebuild before independent
acceptance rather than treating these pre-commit checks as a new candidate run.

The original Codex direction was UniProt grouping plus dense SVG readability.
Claude's read-only review identified the important constraint: payload order
drives derivation IDs, fingerprints and other consumers. The revised minimal
plan uses one pure protocol presentation layout, exact feature-ID joins,
shared GUI row bands and emitted-text geometry tests. The implementation
follows that plan; it does not change the existing biological classifier.

Changed source: `crates/gentle-protocol/src/{lib,splicing_presentation}.rs`,
`crates/gentle-render/src/feature_expert.rs`, `src/main_area_dna.rs`,
`src/main_area_dna/{auxiliary_workspaces,tests}.rs`. Scoped documentation:
decisions, GUI, CLI, protocol, roadmap, changelog and this acceptance handoff.
No routes, dependencies, schemas, tutorial catalogs or scientific fixtures
are changed. Synthetic fixtures are explicitly marked in-memory test data.

## Checks And Limits

Local host: macOS 27.0, arm64; Rust 1.100.0-beta.1
(`e3feeb59c`, 2026-09-27), Cargo 1.100.0-beta.1 (`3d7cf6e9`, 2026-09-25).
`Cargo.lock` SHA-256:
`3dfe44c7a08bf779e32289cf1312103ea1fb5883e44bba5f6e1912ce3ed77745`.
Rust runs use the unoptimized test/dev profile, default features, locked
offline dependencies and one build job. These are correctness checks, not
performance or native GUI acceptance.

```bash
cargo test -p gentle-protocol splicing_ --locked --offline -j 1
cargo test -p gentle-render splicing_ --locked --offline -j 1
cargo test --lib splicing --locked --offline -j 1 -- --test-threads=1
cargo check -q --locked --offline -j 1
cargo fmt --check
git diff --check
python3 scripts/maintenance_chore.py session-close --plan FLNA-splicing-readability-display-only
```

- Protocol: eight focused tests pass, including legacy evidence, grouping,
  shuffled/missing/duplicate ID joins, payload/fingerprint preservation,
  badges, oriented boundary summaries and inert header geometry.
- Renderer: five tests pass; the manual saved-view smoke is intentionally
  ignored in ordinary CI, and was separately run successfully below. Dense
  geometry checks inspect every emitted text element at minimum 8 px, all
  pairwise text overlaps and exclusion from lane exon boxes.
- The new interaction test uses the existing texture-cleaning egui helper and
  simulates click, delayed hover after pointer movement, and secondary click.
  It exposed a real cross-frame context-menu defect: transient response
  positions disappear after opening. The bounded repair uses egui's retained
  menu anchor with the same row hit-testing, not a new cache or biological key.
  The final serial root run passes all 107 focused splicing tests, with zero
  failures or ignored tests (87.26 s test runtime), including the interaction,
  engine and shared-shell paths.
- An intermediate parallel run also failed the unchanged ATtRACT windowed
  submatrix test (`exact_match_bp` instead of `llr_bits_windowed`). The same
  test passes in the serial run. This suggests shared snapshot interference;
  no ATtRACT code, fixture or synchronization policy is changed here, and a
  serial pass is not evidence that the parallel-test issue is fixed.
- Final locked offline Cargo check, formatting and whitespace checks pass.
  Session-close reports four OK, two warnings (intentional dirty files and
  manual scope confirmation), zero failures. Root debug linking emits the existing
  large `__eh_frame` warning; do not infer performance from these checks.
- Full-workspace and native Windows/Linux GUI tests are unrun. Existing
  TP73 CDS and other release acceptance remain separate pending gates.

## Public FLNA Renderer Smoke

Network retrieval was optional and explicit; no renderer fetches evidence.
Temporary directory: `/private/tmp/gentle-flna-splicing-20261007`.
The pre-existing input-producing `target/debug/gentle_cli` self-identifies as
`cee1aef45212d868a9de01e6732319618dfcff2f` (internal.12), binary SHA-256
`cf928dadf4d47b01a6ad2bab23c7add274958df00f6614257bc62ab18f999941`.
This is not a rebuilt CLI for the current patch. It produced the pinned view
and the before-image; a freshly compiled `gentle-render` test rendered the
same saved view after the change, asserting that its JSON remains unchanged.

```bash
target/debug/gentle_cli --state /private/tmp/gentle-flna-splicing-20261007/state.json ensembl-gene fetch FLNA --species homo_sapiens --assembly GRCh38 --flank-bp 2000 --entry-id flna
target/debug/gentle_cli --state /private/tmp/gentle-flna-splicing-20261007/state.json ensembl-gene import-sequence flna --output-id flna_locus
target/debug/gentle_cli --state /private/tmp/gentle-flna-splicing-20261007/state.json shell 'uniprot fetch FLNA_HUMAN --entry-id flna_uniprot'
target/debug/gentle_cli --state /private/tmp/gentle-flna-splicing-20261007/state.json inspect-feature-expert flna_locus splicing 300 > /private/tmp/gentle-flna-splicing-20261007/view.json
target/debug/gentle_cli --state /private/tmp/gentle-flna-splicing-20261007/state.json render-feature-expert-svg flna_locus splicing 300 /private/tmp/gentle-flna-splicing-20261007/flna_before.svg
GENTLE_SPLICING_VIEW_JSON=/private/tmp/gentle-flna-splicing-20261007/view.json GENTLE_SPLICING_VISUAL_DIR=/private/tmp/gentle-flna-splicing-20261007 cargo test -p gentle-render splicing_manual_saved_view_visual_smoke --locked --offline -j 1 -- --ignored --nocapture
/opt/homebrew/bin/rsvg-convert -w 2400 -o /private/tmp/gentle-flna-splicing-20261007/splicing_after.png /private/tmp/gentle-flna-splicing-20261007/splicing_after.svg
```

Returned Ensembl identity: `ENSG00000196924`, gene version 21, GRCh38,
canonical `ENST00000369850.10` / FLNA-202, feature 300. The genomic gene is
negative-strand at X:154348524..154374643; the imported local sequence and
splicing axis are gene-oriented `+`. Do not interpret the local `+` as a
change to genomic strand. The REST response did not pin an annotation-release
number; none is invented here.

The fetched Swiss-Prot record verifies accession P21333, reviewed, sequence
version 4 (2007-01-23), entry version 280 (2026-09-02). The saved view has
36 transcripts, 103 unique exons, two referenced and 34 other evaluated rows.
All 1,858 boundary markers remain; oriented text is summarized in 146 rows.
This is exact loaded-xref evidence, not a validity ranking of the transcripts.

Retained temporary artifacts and SHA-256:

| Artifact | SHA-256 |
| --- | --- |
| `state.json` | `7e1017233a19e4a16e4bc03441e10c1b779a9284bd82cb8b34f458cc28b1678c` |
| `view.json` | `747faaa7dc0e9d879fae339d76f2fb1e7b0df876c338b357499165c671fa09e1` |
| `flna_before.svg` | `6b83f8f5f192077c0d75db7f785f02da65dd423c7bc22624d3ea3e38eb1fc025` |
| `splicing_after.svg` | `8ecb8e6a5124ee2e689c22b4024591c596707d0ade9e5be7aff73727373ce5bc` |
| `splicing_after.png` | `3e97636aff6a79a4a4ad34b23d539dc81c1e48576388967e2e287a7b64f2fc41` |

The full 2400-px PNG, top-of-chart and matrix crops were visually inspected.
This exposed matrix backgrounds covering the ends of accepted headers; all
header backgrounds now precede labels, and cell numbers must fit inside their
own cells. Geometry and paint-order/width regressions protect both repairs.
The main chart uses fixed lane/header heights, not collision-driven expansion.
The complete SVG is still a long report because matrices and the complete
unique-boundary appendix are retained; it is not a one-page publication layout.
Downloaded bytes are not committed and may disappear with temporary-directory
cleanup. Preserve these inputs or retrieve fresh inputs with new provenance
before reproducing this particular comparison.

## Independent Request To Glen

Use the final committed candidate, record source/binary/toolchain/profile
identity, and rebuild the CLI and GUI rather than treating the old input
producer as current acceptance. Repeat the focused commands above and the
public FLNA import/export with returned transcript/feature IDs, not assumed
ones. Retain inputs, view JSON, SVG, untouched 2400-px PNG and hashes.

In the native Splicing Expert, verify shared group/matrix order, target
highlighting, accession/review badge hovers and all boundary ticks. Click a
group header (no selection), then a reordered lane, exon and intron; inspect
hover and right-click intron actions against exact feature IDs. Check that
missing/duplicate matrix joins show unknown rather than absence, and that a
project without relevant UniProt records remains unevaluated without headers.

Report correctness and readability separately from timing, with exact
candidate/platform, individual failures, artifacts and hashes. Preserve
ambiguity and local-evidence wording; do not claim global UniProt absence or
biological validity. Native acceptance, human scientific approval and release
status remain yours/the owner's separate decisions. No push was performed and
no pre-existing user changes were overwritten.
