# Codex prompt: native TSS view export and scoring follow-ups

Self-contained brief. Six scoped items on the native annotated-TSS view. Items
1 and 2 change what Glen would be benchmarking, so they come first; items 3–5 may
be deferred without blocking them. Item 6 is a small, independent docs fix that
should land regardless of the others.

## Verified starting state

- Repository: `/Users/u005069/GitHub/gentle_rs`, branch `main`.
- HEAD at briefing time: `f523a757` ("docs(release): record Glen's tutorial
  checkout validation"). This brief was first written against `f70be167`; the 45
  commits since are tutorial reviews, release preparation and a ClawBio skill,
  and **none of them touch the code targeted below** — every anchor was
  re-verified at `f523a757`. **Re-check `git rev-parse HEAD` before starting**;
  this repository moves under long sessions.
- The routes under discussion landed in `1a107452` ("feat(tss): expose local
  scoring and SVG export to agents") and were corrected in follow-ups.
- Relevant files: `src/tss_sequence_view/svg.rs`,
  `src/engine/ops/operation_handlers.rs` (`Operation::ExportTssViewSvg`),
  `src/tss_sequence_view/local_scoring.rs`, `src/engine_shell.rs`
  (`promoters tss-view-svg`), tutorial
  `docs/tutorial/08-17_tp73_dnp73_factor_curves.md`.

Already done, do not redo: Glen's three Linux screenshots are captured and
registered (`review_status: unreviewed`, 3 graphics). Benchmarks, Windows/macOS
acceptance and merged-candidate GUI/headless equivalence remain **Glen's**, not
yours — see `docs/tp73_dnp73_factor_curves_glen_handoff.md`.

## Release context: this is not a `.12` gate

`v0.1.0-internal.11` is published (confirmed 2026-09-30 at `a51bbc06`, all three
native packages plus the container). The current candidate is
`v0.1.0-internal.12`, and its active aim (`docs/roadmap.md` line 25) is
**responsive startup and DNA-feature presentation with measured opt-level=1
packaging**; TSS scientific and GUI checks are explicitly carried forward, not
gating.

So: none of the items below may block, delay or compete with the `.12` latency
work. Do not add any of them to the release gate, and do not change the
DNA-feature presentation or startup paths while doing them. If an item turns out
to need such a change, stop and report it.

## Interaction with the outer-agent (ClawBio) delegate

`integrations/clawbio/skills/gentle-tss-collection/` (commit `3bd49de6`) routes
08.15's six TSS collection intents — preview, materialize, list, validate, open
windows, forget — with mutating steps behind confirmation-gated routes and
approval digests that must never be invented. It does **not** yet route
`promoters tss-view-svg` or local scoring.

`docs/roadmap.md` line 132 plans exactly that next: extend the 13 cases from
08.16/08.17 through the same delegate with artifact-write and selection
approvals. `promoters tss-view-svg` writes an artifact, so it is one of those
cases, and items 1, 2 and 5 below change its surface or semantics.

- Do **not** implement the ClawBio extension in this task; it is a separate
  roadmap item with its own acceptance.
- Keep the shell surface of `promoters tss-view-svg` backward compatible:
  existing flags keep their meaning, new behaviour arrives through new flags or
  new defaults that are disclosed in the output.
- Item 5 adds a collection target, which changes what the delegate must route
  and approve. Say so explicitly in your changelog entry so the delegate work
  starts from the settled shape.
- Settling items 1, 2 and 5 first lets the delegate be written once, against a
  stable route, instead of twice.

## Non-negotiable invariants

This subsystem's value is its refusal to overclaim. Every item below must
preserve all of these; if an item cannot be done without breaking one, stop and
report that instead of relaxing the invariant.

- **Unavailable is not zero.** `TfbsScoreTrackValidity` masks, `score_at()`
  returning `Option<f64>`, `fully_evaluated()` false when the mask is absent,
  the amber `unavailable-score` bands (upper half local `+`, lower half local
  `-`) and the grey `terminal-unavailable` band all stay exactly as they are.
  Zero remains a real scored value.
- **Separate, non-cross-calibrated scales** for report-backed, locally computed
  and imported DuckDB lanes. Nothing may rank factors across lanes.
- **Exact accessions only.** No factor-name aliases, no `ALL`, no consensus
  fallback on these routes.
- **No new scoring hidden in a presentation path**, no sequence retrieval, no
  motif-database query, and no mutation of the loaded sequence.
- Scores never become occupancy, binding or reporter-activity claims, in code
  comments, messages, docs or tutorial prose.

## 1. Bound the exported SVG's hover-title density

**Problem.** `src/tss_sequence_view/svg.rs` writes one
`<circle data-role="score-window">` with a `<title>` per scored window start
(around the `write!(hovers, ...)` call near line 296). Measured: a 701 bp window
with three matrices produces **1.2 MB** of SVG containing 4,134 such titles.
The documented admission ceiling is 50,000 bp and 32 matrices, which extrapolates
to a file that would "succeed" and then be unopenable. Today SVG size, not
scoring cost, is the practical ceiling — and it is an accidental one.

**Do.** Introduce an explicit, disclosed bound on emitted per-window titles.
Pick one policy and state it in the figure itself:

- prefer local maxima per lane and strand, or a score threshold, or a hard cap;
- when titles are omitted, say so in the lane's own summary text, with the
  emitted and total counts, so a reader never mistakes a missing title for a
  missing window;
- keep the curve geometry and every band unchanged — this reduces *hover
  annotation*, never displayed evidence.

**Design decision for the owner, do not pick silently.** Whether the default
policy is maxima-only or a cap, and what the number is. Propose one, implement
it behind a named constant, and say plainly in the changelog what hover
information a reader loses.

**Also.** Add a regression asserting the title count stays bounded as window
length grows, and that the omission is disclosed. Do not assert an exact byte
size; assert the bound and the disclosure.

## 2. Make `--span` reduce scoring work, not just clip the view

**Problem.** In the `Operation::ExportTssViewSvg` handler, `compute_local_scores`
is called on the whole window sequence, and `start_0based` /
`end_0based_exclusive` only populate `TssViewSvgOptions`. The tutorial's narrow
proximal export (`--span 401..701`) therefore costs exactly as much as the full
export, which is surprising and makes Glen's benchmark ladder measure the wrong
thing.

**Do.** Pass the requested span into the scan target so a partial export scores
only the windows it needs.

**Watch these two traps.**

- `TssLocalScoreRequest::validate_budget` currently validates against the whole
  window length. Decide deliberately whether the budget applies to the scored
  span or the window, and make the refusal message say which. Do not silently
  let a 50,000 bp window through because only 300 bp are exported, unless that is
  the stated intent.
- Motif windows straddling the span edge: a matrix starting before
  `start_0based` can still overlap the exported region. Decide and document
  whether those are scored, and keep the reported evaluated-window counts and
  the `terminal-unavailable` band honest for the span actually scored. A count
  that silently changes meaning between full and partial exports is worse than
  the current redundant work.

If the edge semantics cannot be made unambiguous, implement nothing and report
the ambiguity — a wrong count here undermines the whole "counts differ only
because matrix lengths differ" teaching point in the tutorial.

## 3. Make the TP73 P1 / TAp73 promoter reachable

**Context.** `test_files/tp73.ncbi.gb` is `NC_000001.11 REGION:
3652516..3736201`. The P1 / TAp73 start sits at local base 1, so the inventory
correctly reports `missing_flanks` and no upstream window can be built. That
refusal is good behaviour and the tutorial already teaches it — keep it.

**Do.** Add a route for the interesting comparison: P1 versus the internal
ΔNp73 P2 promoter, which is where the TAp73 → ΔNp73 feedback actually lives.
Either
(a) document an `ExtendGenomeAnchor` step from a prepared genome, clearly marked
as requiring a reference GENtle cannot supply offline, or
(b) retain a second public excerpt that includes upstream flank, with full
provenance in `test_files/README.md` per the repository's fixture rule.

Prefer (b) if a suitable public record exists, because the tutorial's offline
guarantee is worth keeping. Either way, do not weaken the `missing_flanks`
refusal to make P1 appear usable.

This also gives a real-data minus-strand or cross-promoter case, which the
synthetic 08-16 fixture currently carries alone.

## 4. Per-lane strongest-window summary

**Problem.** The view deliberately refuses to rank factors, which is right. But
"where does each factor score highest relative to the TSS" is the next question
any biologist asks, and it is currently answerable only by hovering.

**Do.** Add a per-lane top-N list of window starts with TSS-relative offsets,
in the side summary and in the SVG. Strictly within one lane and one matrix:
no cross-lane table, no shared ordering, no implied ranking between factors.
Say in the text that these are the highest-scoring positions *for that matrix*,
not the most likely binding sites.

Note that `SummarizeTfbsTrackSimilarity` already exists for cross-matrix
questions; do not reimplement any part of it here.

## 5. Collection-wide headless export

**Problem.** `ui open tss-view --collection` opens up to 32 member windows, but
`promoters tss-view-svg` takes a single `seq_id`. A real screen therefore needs
32 invocations and has no coherent scale across them.

**Do.** Accept a collection as an export target, reusing the existing validated
membership path (`GetTssCollection`) rather than a new membership notion.

**Design decision for the owner, do not pick silently.** Whether member figures
share one scale or keep per-window scales. Shared scales make windows visually
comparable but require the cross-window calibration this subsystem otherwise
refuses; per-window scales stay safe but are not comparable. Propose one, make
the choice explicit in each figure's own text, and keep the existing bounded
member count and refusal behaviour.

## 6. Name the outer-agent route in tutorial 08.17

Independent of items 1–5, docs-only, and small. Do it even if the others slip.

**Problem.** The owner's original requirement was that 08.17 be executable by
the inner **or** the outer agent. Since 08.17 was written, `4f6621bd` formalised
the distinction in `docs/tutorial/01-01_agent_interfaces.md`: the **inner agent**
(Agent Assistant) works on the live project and may arrange the GUI through
registered `ui ...` commands; the **outer agent** (MCP or ClawBio/OpenClaw) has
no live project access, works against explicit state, produces artifacts and
receipts, and must not claim to have seen the live GUI. 08.17 says "inner agent"
six times and never names the outer agent, MCP or ClawBio, so the outer half of
the requirement is met in practice but not stated.

**Do.** Using 01-01's vocabulary and linking to it, state per step which agent
can do what:

- Steps 1–3 (`/open file`, `promoters tss-inventory`, `promoters tss-materialize`)
  are headless shared operations: both agents can run them.
- Step 4 and the `ui open tss-view --local-score` form of step 5 are GUI-hosted
  intents: inner agent only. An outer agent receives `applied=false` and must
  report that, not claim the view opened.
- Steps 5–6 via `promoters tss-view-svg` are the outer agent's route, producing
  the artifact and its SVG SHA-256 as the receipt.

Keep the existing "Ask the inner agent" blocks. Add outer-agent guidance only
where it differs; do not duplicate every step. If 01-01 already supplies a reusable
pattern for this, follow it rather than inventing new section headings.

Generated outer-agent replay guidance (`002de212`) is written into
`docs/tutorial/generated/chapters/` only. 08.17 is `hand_written_markdown` and
has no generated section, so edit it directly — but never hand-edit anything
under `docs/tutorial/generated/`.

## Boundaries and verification

- **No local builds or tests** — consistent with recent Codex sessions in this
  repository. State plainly in the changelog which checks did not run.
- **No push, no tag, no workflow dispatch.** Leave commits for the owner.
- Do not touch `Cargo.lock`. It drifted twice during the original work and was
  reverted both times; lockfile digests are bound into TSS report provenance.
- `docs/tutorial/catalog.json` is **generated** from
  `docs/tutorial/sources/*.json` by `gentle_examples_docs`. Edit the source and
  regenerate; never hand-edit the catalog. Keep `review_manifest.json` in step,
  and use only existing `review_status` values.
- Update `docs/cli.md`, `docs/gui.md`, the inner-agent guidance in
  `src/agent_bridge.rs`, the capability descriptor in `src/engine_shell.rs`, and
  the tutorial wherever behaviour changes. A changed default that is not in the
  tutorial's prose is a defect.
- Every item needs at least one deterministic test. For items 1 and 2 the test
  must pin the new semantics (bounded titles with disclosure; span-scoped counts),
  not merely that the command exits zero.
- Record in `docs/CHANGELOG.md` and adjust the native-TSS follow-up sentence in
  `docs/roadmap.md` line 93, which currently states the SVG-size question as
  open and pending Glen's measurements. If item 1 settles it, say so there. If
  the route's shape changes, also note it beside the delegate plan on line 132.
- Line numbers drift; locate those two roadmap sentences by content, not number.

## What not to conclude

None of this work produces biological evidence. The tutorial's three factors
(E2F1 `MA0024.3`, PATZ1 `MA1961.2`, TP73 `MA0861.2`) were chosen because the
TAp73 → ΔNp73 feedback and E2F1's documented role at this locus make them a
coherent teaching set. That is motivation for a figure, not a result from one.
If any change makes the view easier to read as "these factors bind the ΔNp73
promoter", the change is wrong even if the code is correct.
