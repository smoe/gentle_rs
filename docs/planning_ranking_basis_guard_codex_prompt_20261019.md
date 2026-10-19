# Codex prompt: pin the cloning consultation ranking order against silent drift

## Status: Implemented; Glen/CI Verification Pending

Updated 2026-10-19. Implemented in commit
`5a8bd8cec472f68bc73b88448ed135b1131eb792`
(`test(planning): pin cloning consultation ranking policy`). The brief below
is historical design/review provenance, not a pending implementation request.
Its original line anchors and test counts refer to the stated older revisions.

Seven new test-only guards cover default/waiting absolute snapshots, tie
structure, model/weight invariance, policy-time linearity, a genuine reorder
negative control and each row's policy basis/consumers. The existing 13 tests
remain unchanged. The default snapshot has six distinct scores: two tied
blocks plus four singletons, correcting the historical count below. Production
scoring, contracts, assets, adapters, roadmap and release scope are unchanged.

**Verification remains pending with Glen/CI.** No Rust build, test or check was
run locally for these new guards. Earlier 13-test passes cited below do not
verify the expanded 20-test suite or establish current-source acceptance.
Offline source/fixture checks, 17 specialist source tests, scoped rustfmt and
whitespace checks passed; these are not substitutes for Rust execution.

Run on the explicitly selected integrated revision, without re-implementing
this prompt or automatically re-pinning unexpected results:

```bash
cargo test --locked --lib planning_wetlab_
cargo test --locked --lib execute_planning_consult_cloning_returns_ranked_report
cargo check -q --locked
```

Native, cross-platform, scientific and release acceptance remain separate. The
text report still shows only the top five candidates; the guards inspect all
eleven JSON rows. The previous paused goal is not resumed by this status note.

---

Prepared 2026-10-19 against `main` at
`8af42b9792da05e53422a53e04e12509bff1f429`
("fix(planning): separate bench estimates from ranking policy").

Re-verified 2026-10-19 against `main` at
`84ebede4ba144eb065d43ba19aa00b9c2f1ce47a`
("test(expression): align CLI preflight with missing-adapter refusal"). That
commit touches only `tests/protein_expression_handoff_cli.rs`, its fixture
README, the changelog and one roadmap line, so every anchor below is unchanged:
`estimate_routine_planning_with_bench_catalog` at `src/engine_shell.rs:4706`,
`policy_time_hours` at `:4813`, the guardrail at `:4887`, the score at `:4904`,
the shared display line at `src/engine_shell/planning_wetlab.rs:362`, the
comparison axis at `:404`, the CLI report at `src/engine_shell.rs:34560`, the
GUI at `src/app/routine_and_agent_assistant_ui.rs:4838` and `:4924`, and the
two landed tests at `src/engine_shell/planning_wetlab_tests.rs:103` and `:194`.
The 13 `planning_wetlab_` tests pass at this SHA on macOS.

`main` moves under long-running sessions in this repository. Re-verify
`git rev-parse HEAD` before starting, and re-locate every anchor below by
symbol name rather than trusting the line numbers if HEAD has moved again.

**Tests and documentation only, in one commit.** No scoring change, no
contract change, no asset change, no new bench model. If you find yourself
editing `estimate_routine_planning`, you have misread this prompt.

## Read this first: the cross-basis ranking defect is already fixed

This prompt was commissioned to answer a question that `8af42b97` had already
answered, and the investigation below is recorded so the work is not done
twice.

The commissioned question was: *should a routine whose time estimate comes
from a curated bench-step model be ranked against routines whose estimate is
still the legacy family scalar?* At `293819fd` that question was live —
`composite_meta_score` subtracted `weight_time * estimated_time_hours`, and
`restriction.digest_ligate_extract_sticky`'s displayed estimate had just
fallen from 6.0 h to 3.2 h because it gained the only bench-step model.

`8af42b97` answered it by construction. `estimate_routine_planning_with_bench_catalog`
(`src/engine_shell.rs:4706`) now captures

```rust
let policy_time_hours = estimated_time_hours;   // src/engine_shell.rs:4813
let model = bench_catalog.model(&routine.routine_id);   // :4814
```

one line *before* the bench catalog is consulted, and both policy consumers
read that captured value rather than the model-derived one:

- `composite_meta_score` at `src/engine_shell.rs:4902`, which subtracts
  `objective.weight_time * policy_time_hours` (`:4904`),
- the `max_time_hours` guardrail at `src/engine_shell.rs:4887`.

The explanation records the separation at `src/engine_shell.rs:4911`-`4916`
as `time_policy` with `hours`, a `legacy_family_scalar` basis and
`used_by: ["composite_meta_score", "max_time_hours"]`. The invariant is
written down in [`docs/decisions.md`](decisions.md) ("Bench-model presence and
duration must not change policy outcomes") and
[`docs/architecture.md`](architecture.md).

So of the four options the commission offered:

- **(a) rank only within a basis, legacy candidates as a separate group** —
  unnecessary. There is no longer more than one ranking basis. Every candidate
  is ranked on the legacy family/base scalar plus summed delays, whether or not
  it has a bench model. There are no incomparable groups to separate.
- **(b) surface the basis on every row** — already done. Every row's estimate
  carries `estimated_time_hours_basis`; the shared display lines emit
  `Ranking/time-guardrail proxy: …` (`src/engine_shell/planning_wetlab.rs:362`-`365`),
  consumed by the CLI text report (`src/engine_shell.rs:34560`) and the GUI
  (`src/app/routine_and_agent_assistant_ui.rs:4838`, `:4924`); `routines compare`
  carries a `time_policy_proxy_hours` axis
  (`src/engine_shell/planning_wetlab.rs:404`-`411`).
- **(c) neutralize the basis in scoring** — effectively what landed, but in the
  honest direction. Rather than deriving a comparable figure for the legacy
  routines, `8af42b97` keeps the legacy figure for *every* routine and labels it
  a policy proxy rather than a laboratory duration. Do not revisit this.
- **(d) a guard test pinning the current ordering** — **the only option with
  work left in it, and the subject of this prompt.**

## What the investigation measured

Commands run on macOS at `8af42b97`, against a fresh state file outside the
repository:

```bash
cargo test --locked --lib planning_wetlab_
cargo build --locked --bin gentle_cli
./target/debug/gentle_cli --state "$WORK/s1.json" \
  shell 'planning consult cloning --format json'
```

`cargo test --locked --lib planning_wetlab_` passes 13/13, including
`planning_wetlab_model_presence_and_duration_do_not_change_ranking_or_guardrails`
(`src/engine_shell/planning_wetlab_tests.rs:103`) and
`planning_wetlab_consultation_order_matches_the_no_model_counterfactual`
(`:194`).

The measured ranking under the default objective, no profile override:

| rank | family | score | policy h | displayed h | basis |
| --- | --- | --- | --- | --- | --- |
| 1 | sequence | -0.1500 | 1.0 | 1.0 | legacy_family_scalar |
| 2 | crispr | -3.0000 | 4.0 | 4.0 | legacy_family_scalar |
| 3 | golden_gate | -3.0000 | 4.0 | 4.0 | legacy_family_scalar |
| 4 | pcr | -3.0000 | 4.0 | 4.0 | legacy_family_scalar |
| 5 | ta_gc | -3.0000 | 4.0 | 4.0 | legacy_family_scalar |
| 6 | topo | -3.0000 | 4.0 | 4.0 | legacy_family_scalar |
| 7 | infusion | -3.1500 | 4.0 | 4.0 | legacy_family_scalar |
| 8 | nebuilder_hifi | -3.1500 | 4.0 | 4.0 | legacy_family_scalar |
| 9 | gibson | -4.1500 | 5.0 | 5.0 | legacy_family_scalar |
| 10 | **restriction** | **-5.3000** | **6.0** | **3.2** | **step_model_critical_path** |
| 11 | gateway | -7.1500 | 8.0 | 8.0 | legacy_family_scalar |

Three findings matter for the test design.

**1. The reorder the commission feared is real in magnitude, and is prevented.**
Sweeping `weight_time` over `{0.0, 0.5, 1.0, 2.0}` confirms the score is exactly
linear, `score = K - weight_time * policy_hours`, with `K` constant per family
(`K`: crispr/golden_gate/pcr/ta_gc/topo `1.00`; gateway/gibson/infusion/
nebuilder_hifi/sequence `0.85`; restriction `0.70`). Scoring the restriction
routine on its *displayed* 3.2 h instead of its policy 6.0 h would give
`0.70 - 3.2 = -2.5`, placing it **rank 2 of 11** — an eight-place jump, earned
by nothing but being curated first. Solving the same linear law pairwise, that
leak would change an adjacency at any `weight_time > 0.03125` (restriction
overtaking gateway), pass infusion/nebuilder_hifi above `0.1875` and the
five-family `-3.0` cluster above `0.375`. It is not a knife-edge that the
default weight happens to sit on; it bites almost immediately.

At `8af42b97` no `weight_time` can produce that reorder, because the quantity
`weight_time` multiplies is itself basis-independent and bit-identical with and
without a model (asserted at
`src/engine_shell/planning_wetlab_tests.rs:157`-`160`). This is a stronger
statement than "the order is unaffected at default weights": it is unaffected
at *every* weight.

**2. Seven of the eleven rows are not ranked by score at all.** Only four
distinct scores exist across eleven families: `-3.0` is a bit-exact tie across
five families (`0xc008000000000000`) and `-3.15` across two
(`0xc009333333333333`). Those seven rows are ordered by the alphabetical
tiebreak in `strategy_candidates.sort_by` (`src/engine_shell.rs:35367`-`35372`:
score descending via `total_cmp`, then `family`, then `routine_id`). Any guard
that pins the order is therefore mostly pinning an alphabetical accident, and
must say so in its own comments — otherwise a future reader will mistake the
pinned sequence for a ranking judgement. A change to one legacy family scalar
would break a five-way tie and move a block of rows at once.

**3. The existing order test is differential, not absolute — this is the gap.**
`planning_wetlab_consultation_order_matches_the_no_model_counterfactual`
builds its expectation by calling the *same production scoring function* with
an empty bench catalog, then compares rank-for-rank. That correctly pins "a
bench model cannot change the order" and nothing else. A change to
`default_routine_base_time_hours`, to `routine_family_alignment_bonus`, to
local-fit scoring, or to the `PlanningObjective::default` weights
(`src/engine.rs:3016`) moves both sides of that comparison together: the
visible ranking reorders and the test still passes.

The commission already identified that
`execute_planning_consult_cloning_returns_ranked_report`
(`src/engine_shell/tests.rs:14742`) cannot catch a reorder — it collapses the
families into a `BTreeSet` and asserts only that there are 11 of them. That
remains true. Across both files, no test asserts an absolute rank sequence or
any pinned score value for the bundled catalog.

A fourth, smaller observation: the CLI text report prints only
`strategy_candidates.iter().take(5)` (`src/engine_shell.rs:34551`), so a reorder
below rank 5 is invisible in `--format text` even to a reader who is looking.
The guard must work on the JSON report.

## Constraints

- **Tests and documentation only.** No change to `estimate_routine_planning`,
  `estimate_routine_planning_with_bench_catalog`, `composite_meta_score`,
  the guardrail arithmetic, the selection or sort in
  `src/engine_shell.rs:35349`-`35372`, or any adapter.
- No change to the landed `gentle.planning_estimate.v2` contract, to
  `gentle.cloning_wetlab_steps.v1`, or to **any** field name. In particular
  `time_policy`, `estimated_time_hours_basis`, `hours`, `basis` and `used_by`
  stay exactly as they are.
- No change to `assets/cloning_wetlab_steps.json`,
  `assets/cloning_routines.json`, or any other asset. **Do not add a bench-step
  model for any further routine.** The tripwire at
  `src/engine_shell/planning_wetlab_tests.rs:553`
  (`planning_wetlab_all_bundled_routines_reconcile_and_only_one_has_a_model`)
  is deliberate; leave it able to fire.
- Deeper wet-lab process modeling is in the [roadmap](roadmap.md) parking lot
  and `.13` is an active release candidate. This work must **not** enter `.13`
  acceptance, release notes, or any candidate ledger. **Do not edit
  `docs/roadmap.md`** — it must stay exactly 400 lines and the release-candidate
  checks assert that count.
- Per `AGENTS.md`, no `Done` entry in the roadmap. The changelog entry goes in
  [`docs/CHANGELOG.md`](CHANGELOG.md) in this same commit.
- **Always pass `--locked`.** A bare `--offline` rewrites `Cargo.lock` in this
  workspace and produces fabricated compile errors.
- `main` carries pre-existing unrelated test failures and a broken
  `cargo fmt --all`. Establish your baseline at `8af42b97` before attributing
  any failure to this change.
- Nothing here may imply that a curated bench model makes a protocol
  biologically better, more likely to succeed, or better ranked. The pinned
  numbers are a synthetic policy snapshot, not laboratory evidence and not a
  method recommendation.

## The guard to add

Place the new tests beside the existing bench-model tests in
`src/engine_shell/planning_wetlab_tests.rs` (wired as a `#[path]` module from
`src/engine_shell/tests.rs:14722`).

Pin the ranking as a **policy snapshot** with these properties:

- Assert the `(rank, family, routine_id)` sequence exactly, as an ordered
  literal, against the bundled `assets/cloning_routines.json` and
  `assets/cloning_wetlab_steps.json` through the real
  `ShellCommand::PlanningConsultCloning`, `--format json`, with
  `objective_json: None` and `PlanningProfileScope::Effective`.
- Pin each row's `composite_meta_score` to a literal, compared with the `1e-9`
  tolerance the existing comparison axes already use. Do not compare raw `f64`
  with `==` across platforms, and do not pin `to_bits()` for the snapshot —
  bit-equality is the right tool for the *same-run counterfactual* at `:157`,
  not for a committed cross-platform literal.
- Record the tie structure in the test itself: the five-family `-3.0` tie and
  the two-family `-3.15` tie are resolved by the documented alphabetical
  tiebreak, not by score. State in a comment that re-sorting those rows is a
  tiebreak artefact rather than a ranking change.
- Fail with a readable diff — print the expected and actual sequences
  side by side, not a bare `assert_eq!` on two eleven-element vectors.
- Say in a comment that updating these literals is an owner decision requiring
  a `docs/CHANGELOG.md` line, so a future session cannot quietly re-bless a
  reorder by editing the expectation.

Assert the state is unchanged either side of the call, as the neighbouring
tests do (`report.state_changed` false and a before/after `engine.state()`
comparison) — `planning consult cloning` is read-only.

## Deliberately out of scope

- Any scoring, weighting or ranking-policy change, including a calibrated
  cross-method ranking. [`docs/decisions.md`](decisions.md) makes that a
  separately approved policy decision.
- Grouping, hiding or re-presenting candidates by basis — options (a) and (b)
  above, both settled.
- Bench-step models for any of the remaining 21 routines.
- Raising the `take(5)` limit in the text report, or any other adapter change.
  Note it as a known limitation in the handoff; do not fix it here.
- Package sizes, vendor offers, attempt/repeat multipliers, success
  probability or yield prediction.
- Re-pinning or relabelling the historical reports that predate `time_policy`.
  Absent metadata stays unknown and must not be reconstructed.

## Tests

1. **Absolute order snapshot, default objective.** Exactly 11 rows; the
   `(rank, family, routine_id)` sequence matches the pinned literal; each
   `composite_meta_score` is within `1e-9` of its pinned literal. Expected
   values are those measured above: `sequence -0.15`, then the `-3.0` block
   `crispr, golden_gate, pcr, ta_gc, topo`, then the `-3.15` block
   `infusion, nebuilder_hifi`, then `gibson -4.15`,
   `restriction.digest_ligate_extract_sticky -5.3`, `gateway -7.15`.
   Re-measure before committing; do not trust this table if HEAD has moved.
2. **Absolute order snapshot under the waiting profile.** The same pin with
   [`docs/examples/planning/cloning_waiting.json`](examples/planning/cloning_waiting.json)
   applied as a `ProjectOverride`, since queue and purchasing delays enter
   `policy_time_hours` and so enter the score. Pin that order and those scores
   separately; do not assume they match test 1.
3. **The tie structure is what the pin says it is.** Assert that exactly five
   families share the top tied score and exactly two share the next, that each
   tied block is in ascending `family` order, and that the four remaining rows
   have pairwise-distinct scores. If a future change makes all eleven scores
   distinct, this test must fail and force the comment in test 1 to be revised
   rather than silently becoming wrong.
4. **`weight_time` cannot reintroduce a basis difference.** For
   `weight_time` in at least `{0.0, 0.03125, 0.375, 1.0, 2.0, 1e6}`, the
   restriction routine's `composite_meta_score` is **bit-identical** with the
   bundled bench catalog and with an empty one, reusing
   `estimate_routine_planning_with_bench_catalog` as `:103` already does. This
   turns "no reorder at any weight" from an argument into an executed
   invariant, and the chosen weights straddle the crossover thresholds derived
   above.
5. **The score is linear in `weight_time` with a basis-independent slope.**
   For each ranked family, assert
   `score(w) == K - w * explanation.time_policy.hours` within `1e-9` at two
   distinct weights, with `K` recovered at `w = 0.0`. This pins the mechanism
   that makes test 4 generalize beyond the sampled weights.
6. **The guard is load-bearing.** Prove the pin can fail: recompute the ranking
   with one deliberately perturbed input — for example the same consultation at
   a different `weight_time`, or a locally cloned routine catalog with one
   family's base time altered — and assert the resulting sequence differs from
   the pinned one. A snapshot test that cannot fail is worse than none, and
   this is the assertion that distinguishes the two.
7. **Every ranked row still declares its policy basis.** For all 11 rows,
   `explanation.time_policy.basis == "legacy_family_scalar"` and
   `used_by` contains both `composite_meta_score` and `max_time_hours`,
   including the one row whose `estimated_time_hours_basis` is
   `step_model_critical_path`.
8. **`execute_planning_consult_cloning_returns_ranked_report`
   (`src/engine_shell/tests.rs:14742`) still passes unchanged.** Do not edit
   it; its `BTreeSet` family check is order-blind by design and the new tests
   cover ordering.
9. **The existing bench-model suite still passes unchanged** — all 13
   `planning_wetlab_` tests, in particular `:103`, `:194` and `:553`. The new
   absolute pin complements the differential counterfactual at `:194`; it does
   not replace it, and `:194` must not be deleted or weakened.

Run `cargo check -q --locked` and `cargo test --locked --lib planning_wetlab_`
before handoff. Full-suite, native-GUI, Windows and scientific acceptance stay
with Glen and CI; a local targeted pass is not release acceptance and must not
be reported as one.

## Documentation to update in the same commit

- [`docs/CHANGELOG.md`](CHANGELOG.md) — that the consultation ranking order is
  now pinned as a policy snapshot, that the pin is mostly tiebreak-ordered, and
  that no scoring behavior changed.
- [`docs/testing.md`](testing.md) — a short note on the pinned-ranking policy
  snapshot: what it asserts, that updating the literals is an owner decision,
  and that it is a drift guard rather than a correctness or biological claim.

In the handoff note list files changed, tests run, and explicitly: that no
scoring, contract or asset change was made; that the pinned scores are a
synthetic policy snapshot rather than laboratory timings; that seven of the
eleven rows are decided by an alphabetical tiebreak rather than by score; that
the text report still shows only the top five candidates; and that no wet-lab
or platform acceptance was obtained.
