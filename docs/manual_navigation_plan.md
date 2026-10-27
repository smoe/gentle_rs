# Manual Navigation Implementation

Owner-authorized 2026-10-27 after a manually relayed read-only Claude review.
Scope: retain the complete GUI manual, introduce a hyperlinked document tree
and contextual window Help, then give the CLI manual the corresponding topics.
No engine/protocol behavior changes. No local Rust/Cargo execution.

## Baseline And Evidence

The working branch was advanced to local main at `1804605f` before editing;
its GUI manual is identical to `a6ba0da6`, the explicit text-comparison baseline.
The latest Dotplot additions are included. The apparent loss is primarily a
navigation problem: RNA text displaced the introduction, and the sectioned
reader showed that preamble first. Historical heading/list checks found no
deletion, but do not prove that every historical paragraph is retained.

Normal Help buttons already exist in the DNA and specialist navigation bars.
They open the general manual, separately from Agent Help. Reuse them; never
route manual navigation through capture or a model request. `window.specialist`
and translated/window-instance titles are not suitable documentation identities.

## Delivery Steps

- [x] S0: retain an explicit Git baseline and a read-only content-migration check.
  The check preserves nonblank text lines including repeated occurrences and
  code/list indentation; it deliberately does not claim Markdown semantics.
- [x] S1: move the RNA/catalog preamble into named sections and introduce topic
  headings within Main window layout, without deleting original content.
- [x] S2: author GUI same-page and cross-page anchors, bounded section rendering,
  visible tree navigation, breadcrumbs and section history, reusing existing
  search/render caches. UTF-8 fragments, duplicate headings, LF/CRLF and real
  link-click/no-capture regressions are authored, not Rust-executed. Cross-book
  history and non-GUI fragments follow the book migration, not this first step.
- [x] S3 (authored; Rust/native checks pending): make `docs/gui/` topic files authoritative, with `docs/gui.md` as an
  index and legacy-link compatibility surface. Retain offline package and
  embedded loading, independent of current working directory; no index-only
  fallback or mixed-source page collection.
  There are 73 topic files and 83 legacy destinations. `docs/manuals.json`
  owns their ordering, groups and aliases; `scripts/manual_docs.py` generates
  navigation and the embedded source inventory. Missing chapters or a different
  external catalog select the complete embedded book, never mixed chapters.
  Linux/Windows sibling-resource lookup and runtime tests for all three package
  layouts are authored. Package staging/smoke now reject missing/untracked
  catalogued chapters, including a checksum-consistent index-only package.
- [x] S4 (authored; Rust/native checks pending): bind existing Help buttons to stable window topics; cover all window
  families, actual dispatch, target existence and no-capture behavior.
  `HelpTopic` supplies 40 explicit, typed destinations. All existing ordinary
  Help buttons use it; shared specialist bars require it as an argument, so
  embedded windows cannot accidentally inherit the root viewport's topic.
  Analysis windows with only Agent Help now also expose ordinary Help. Registry
  entries require a topic, including dynamic sequence auxiliary windows. The
  existing public three-field auxiliary listing remains compatible. New native
  menu handling uses the active registered window; missing context uses Main.
  Python checks all 40 destination paths/fragments. Authored Rust tests click
  every topic button, consume its real request and check unchanged project/no
  capture, plus registry and active-window routing; they have not run.
- [x] S5 (authored; Rust/native checks pending): reorganize CLI content in the
  shared catalog with explicit bidirectional counterpart links; retain
  standalone CLI topics and commands. There are 69 CLI topic files and 31
  relocated legacy headings. Both books share pagination, file-scoped anchors,
  Find/Copy, lazy image preparation and cross-book Back/Forward. Every topic
  has reciprocal counterparts or an explicit absence reason, with no assertion
  that related interfaces have identical capability coverage.
- [ ] External verification: focused Rust tests, Cargo check and formatting;
  native window navigation, packaged Linux/Windows/macOS lookup from unrelated
  working directories, links, images, Find and wheel/trackpad interaction.

## Verification

Compare against the integrated main source (read-only). Earlier S1/S5 counts
below record their original baselines before the owner-requested rebase:

```sh
python3 scripts/check_manual_migration.py --base-ref 1adad235 --relocated-links --repair-map docs/manuals.json docs/gui.md docs/gui/*.md
python3 scripts/check_manual_migration.py --base-ref 1adad235 --source docs/cli.md --relocated-links docs/cli.md docs/cli/*.md
python3 scripts/manual_docs.py
python3 -m unittest scripts.test_tutorial_checkouts.DocumentationStatusTests
python3 -m unittest scripts.test_package_desktop
```

S1 must retain all 7,164 nonblank source lines including multiplicity and
indentation. Later path rewrites require separately reviewed link-equivalence
checks; do not normalize all whitespace or update expected hashes to hide loss.
The migration aid is not a permanent freeze on future documentation edits.
Permanent guards protect orientation, topic/window coverage and navigation.

All static results must remain distinct from unrun Rust and native acceptance.
Existing LF/CRLF checkout tests remain in the fast no-Rust CI path. A generated
navigation/embedding freshness guard, if needed, must include its scoped LF
policy and real checkout regression in that path.

Initial verification: the S1 migration check retained all 7,164 nonblank lines,
and all 13 `DocumentationStatusTests` passed. Added GUI behavior and pagination
tests have not been executed because local Rust/Cargo is prohibited. This is
an in-flight implementation record, not a final runtime acceptance report.

S3 static verification retained all 7,164 lines after relocating links. Six
pre-existing, unnumbered tutorial URLs were repaired against the actual catalog;
their exact mapping is recorded in `legacy_resources`, not silently normalized.
All 17 targeted documentation/checkout tests and all 35 synthetic package tests
passed. Package executable calls in these tests are mocked: no Rust binary ran.
The new real-click, embedded-book, file-scoped-anchor and fallback Rust tests
remain unexecuted. Images/tutorials still need installed resources; embedding
all manual text does not embed screenshots or the tutorial artifact archive.

S5 source verification compares against `1804605f`, including the newer
observed-transcript reference: all 8,862 nonblank lines survive with exact
whitespace and multiplicity after link relocation. The CLI RNA preamble is in
named chapters; the index retains discovery, saved-project instructions, operand
conventions and the task index. Nested command lists retain their indentation
and explicit parents when separated by new headings. The reciprocal links
are generated from a single pair list; missing/unknown/duplicate bindings fail
validation. Existing command-coverage and example tests now read their
canonical topics or the complete catalogued CLI book, not just the index.

Final scoped Python run: 100 tests passed across documentation status, manual
freshness/counterparts/fences, the real LF/CRLF checkout test, synthetic desktop
packages, primer payload/trace, observed-transcript and specialist tutorial
source checks. This includes negative controls for a missing scoped LF rule,
missing counterpart and checksum-consistent packages missing a GUI/CLI chapter.
All package executable invocations are mocked. No Rust/GENtle binaries ran.
Both content-migration checks, `manual_docs.py`, `check_documentation_status.py`
and `git diff --check` pass. Rust cross-book click/history, search/cache and
whole-book fallback tests are authored, not executed.

## External Acceptance Handoff

At the eventual candidate SHA, Glen/CI should run the focused Rust tests and
normal formatting/type checks. This list is a handoff, not local execution:

```sh
cargo test -q --lib help_manual
cargo test -q --lib help_docs
cargo test -q --lib cli_docs_mention_all_cli_glossary_paths
cargo test -q --lib sequence_design_window_gc_documented_walkthrough_is_replayable
cargo test -q --lib genome_edit_intent_documented_operation_is_the_tested_synthetic_request
cargo check -q
cargo fmt --check
```

Review narrow/wide Help windows, tree selection, both directions of a GUI/CLI
counterpart, old bookmarks, Find across chapters, Copy Manual, Reload and wheel
scrolling. Check ordinary Help in sequence, specialist and detached windows,
with no capture or model request. Test installed macOS/Linux/Windows packages
from an unrelated working directory, with intact and missing chapter resources;
confirm a complete embedded fallback plus a visible notice, not mixed chapters.
Images and tutorial artifacts require installed resources. Native Windows and
native render/performance acceptance remain explicitly unverified. Tutorial
and other non-book guides still disclose their existing section-scroll limit.

The earlier `scripts.test_tutorial_checkouts` run executed 77 tests and had two
failures, both LF/CRLF subcases of capability-maturity freshness. The committed
`1804605f` already records glossary SHA-256 `2f1a98cb...`, while its committed
`docs/glossary.json` hashes to `c324fe17...`. Neither that input nor the maturity
catalog was changed here. Main subsequently corrected that hash independently;
retain the earlier failures as historical evidence, not current manual failures.

## Rebase Integration, 2026-10-27

Integrate local main `1adad235` without invoking Rust/Cargo. Keep both changelog
histories. Transfer main's new fresh-only sequence-design writer guarantees and
limits into `docs/gui/detailed-reference.md` and `docs/cli/sequence-design.md`,
retaining the navigation-only indexes. The integrated source comparison now
retains all 7,171 GUI and 8,869 CLI nonblank lines, with exact indentation and
multiplicity after the documented link moves. The extra seven lines per book
come from main's intentional policy update, not recovered or invented content.

The expanded integration check passes 164 Python tests, including the full
`scripts.test_tutorial_checkouts` suite and the earlier package/source-contract
checks. The historical maturity-hash failure is absent with main's correction.
Manual freshness, documentation status and whitespace checks also pass. Rust,
Cargo and native GUI execution remain unrun.
