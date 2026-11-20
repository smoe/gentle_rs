# Owner-supplied Claude read-only review

Received 2026-11-16. The following is the supplied review, retained verbatim.
Its execution claims are reviewer-reported; Codex reconciliation is in
[the revised plan](transcription_termination_annotation_plan.md).

Read-only review complete at `b62d21c9` — no files edited, no code run, working tree untouched.

## 1. Major changes needed

**a. Phase 1 duplicates work the roadmap already parks.** [docs/roadmap.md:398](docs/roadmap.md:398) parks "Report-first curated feature-library auto-annotation with a separately confirmed apply step" and already links your plan; [docs/legacy_gentle_m_intake.md:48](docs/legacy_gentle_m_intake.md:48) says it needs "source/version binding, overlap/ambiguity rules and a separate confirmed apply step" — that is Phase 1's contract verbatim, minus the catalog. "Audit the existing feature-library proposal" ([plan line 120](docs/transcription_termination_annotation_plan.md:120)) defers the only decision that matters. Choose: **(i)** Phase 1 *is* the generic exact-recognition primitive, with termination elements as its first catalog, or **(ii)** termination is a *consumer* and you ship only the catalog plus evidence typing now. Otherwise you build a second private matcher that the library work must later absorb.

**b. No curated element catalog exists anywhere in the tree.** I found no terminator or cleavage/polyA element sequences under `src`, `crates`, `test_files` or `assets`. Phase 1 is therefore "author a redistribution-cleared reference collection, *then* recognize" — the licensing clearance is the real cost and it is owner work, not agent work. Reuse the existing provenance vocabulary rather than inventing one: `redistribution_status` at [src/genomes.rs:516](src/genomes.rs:516), and `HelperConstructSemantics`/`HelperConstructComponent` at [src/genomes.rs:743](src/genomes.rs:743), which already types `kind: "terminator"` → `transcription_termination` ([src/genomes.rs:4031](src/genomes.rs:4031)). Note explicitly that `HelperConstructComponent` carries **no sequence** — it is a naming layer, not a matchable catalog, so the plan's citation of genomes.rs as "classification, not a detector" understates how far it is from one.

**c. Phase 2's folding wording hides a module that cannot be reused as written.** [src/rna_structure.rs](src/rna_structure.rs) already shells out to ViennaRNA (`DEFAULT_RNAFOLD_BIN`, env override through `tool_overrides`). But `ensure_rna_input` ([src/rna_structure.rs:195](src/rna_structure.rs:195)) refuses anything whose `molecule_type` is not RNA/ssRNA, so it cannot be pointed at a plasmid DNA window; it folds the *whole* sequence with no windowed or stem-geometry API; and `fold_rna_sequence` records executable, args and stdout but **no ViennaRNA version**. The consequence your plan misses: the precedent's materialize recomputes the screen and refuses unless `content_sha256` matches ([src/engine/analysis/tata_boxes.rs:622](src/engine/analysis/tata_boxes.rs:622)). Any external MFE inside the hashed report body makes that digest install-dependent, so materialize would spuriously refuse across machines and CI. Write the rule into the contract: external folding is optional annotation **outside** the hashed body, or the report pins a tool version and refuses mismatches — and the deterministic detector stays pure engine code either way.

**d. Circularity in Phase 1 contradicts the nearest precedent.** [src/engine/analysis/cryptic_splicing.rs:511](src/engine/analysis/cryptic_splicing.rs:511) refuses an origin-wrapping span outright ("materialize a linear span first"), and the TATA screen ignores topology entirely. The only real wrap semantics live in [src/restriction_enzyme.rs:349](src/restriction_enzyme.rs:349) (`can_cleave`, `recognition_bounds_for_topology`, test at line 576). Make Phase 1 **refuse** origin-crossing search with a structured error and report the refusal; wrap becomes its own slice reusing those bounds rather than a third wrap implementation. This deletes a whole acceptance row from the first slice.

**e. You point at the wrong writeback seam first.** `construct_reasoning_annotation_candidate_is_writeback_eligible` ([src/engine.rs:22305](src/engine.rs:22305)) requires `source_kind == "generated_annotation"` plus `EvidenceScope::SequenceSpan`; the role→kind map (22313-22338) and the `regulatory_class` push (22455-22475) are construct-reasoning-graph specific. Your plan warns against forcing records onto it, then still lists `engine.rs` first. Invert: copy [crates/gentle-protocol/src/tata_boxes.rs](crates/gentle-protocol/src/tata_boxes.rs) + [src/engine/analysis/tata_boxes.rs](src/engine/analysis/tata_boxes.rs) end to end — schema const, non-claim const, `content_sha256`, `expected_report_sha256` + `row_ids`, per-row idempotence key as a qualifier, `prepare_sequence` after the feature push ([tata_boxes.rs:734](src/engine/analysis/tata_boxes.rs:734)). `engine.rs`'s only role is two `Operation` variants (cf. [4752](src/engine.rs:4752)) and the read-only classification list (cf. [10504](src/engine.rs:10504)).

## 2. Unsupported biological claims still in the revised doc

- **"Use RNA pairing semantics, including any permitted G-U pairing"** ([line 154](docs/transcription_termination_annotation_plan.md:154)). A wobble allowance inside a DNA inverted-repeat screen is not RNA folding; calling it a pairing *semantic* implies structural validity the heuristic lacks. Demote it to a declared, versioned mismatch/wobble policy knob.
- **Phase 3 cleavage windows.** The doc never states that a hexamer plus downstream U/GU-rich element predicts *site usage* only with expression evidence. PMC7615648 supports your cleavage-vs-termination separation, not a positional predictor. Add: thresholds are authored defaults with no measured sensitivity or specificity until benchmarked, and hexamer absence does not imply absent 3′ processing.
- **"Prokaryotic" in the title versus bacterial-only delivery.** Archaea run a Pol II-like machinery with its own termination; [line 158](docs/transcription_termination_annotation_plan.md:158) correctly excludes them, so the title promises what Phase 2 refuses. Say "bacterial" wherever the deliverable is bacterial.
- **Rho.** "Imported evidence only" is right, but silence is the wrong output. Bind [line 70](docs/transcription_termination_annotation_plan.md:70)'s "report unsupported mechanisms" to a concrete unsupported-mechanism report row, so a user who asks why no Rho call appeared gets an answer.
- **Unannotated-UTR scope is still open-ended.** "Do not assign a hit to the nearest CDS" says what not to do, not what to do. Minimal fix: the request takes one `seq_id` + interval **or** one resolved transcript/CDS feature id, refusing ambiguity the way `resolve_single_transcript_design_template` already does ([docs/biological_extension_guide.md:120](docs/biological_extension_guide.md:120)); whole-plasmid both-strand scans return **unassigned** candidates with no transcript association at all. Association is reported only when the caller named the transcript.

## 3. Hidden scope and acceptance gaps

**Parity is not Phase-4 polish.** The TATA precedent needed glossary paths ([docs/glossary.json:2650](docs/glossary.json:2650), 2693), capability descriptors with the mutating flag ([src/engine_shell.rs:22514](src/engine_shell.rs:22514), 22516), the registry name list (crates/gentle-protocol/src/lib.rs ≈5762), MCP wiring ([src/mcp_server.rs:6712](src/mcp_server.rs:6712), 6724), a GUI module (`src/main_area_dna/tata_box_ui.rs`), a persistence field ([src/engine/protocol.rs:5838](src/engine/protocol.rs:5838)) and protocol docs. Per DEC-003 that surface is mandatory and roughly as large as the biology. Cost it as a phase, and fix one command noun now — `terminators screen` / `terminators materialize` — so the vocabulary is not negotiated twice.

Four acceptance gaps: **(i)** "benchmark against independently curated examples" names no source; with no catalog in the tree and no web access in the implementing session that criterion is currently unachievable — mark it owner-blocked, not acceptance. **(ii)** Fixture provenance for the catalog itself is unplanned (AGENTS.md mandatory rule, DEC-017); the TATA tests solve it with hand-crafted non-biological sequences ([tata_boxes.rs:745](src/engine/analysis/tata_boxes.rs:745)) — do the same and say so. **(iii)** If the catalog is hashed it *is* a byte-bound input, so the `scripts.test_tutorial_checkouts` LF/CRLF gate applies from Phase 1, not "before adding file I/O" later. **(iv)** Nothing tests that an existing `terminator` annotation is recognized and not duplicated — that test must use a real feature with `regulatory_class=terminator` as produced at [src/engine.rs:22463](src/engine.rs:22463) and parsed at 22886-22909.

## 4. Revised minimal phases

| Phase | Deliverable | Acceptance |
|---|---|---|
| **A** Contract + catalog | Library-ownership decision; schema; catalog format with `redistribution_status`-style provenance; hand-crafted fixtures. No matcher. | Schema round-trip, unknown-field rejection, enumerated refusals (wrap, ambiguous transcript, unsupported mechanism, insufficient flank). Zero biological claim. |
| **B** Exact recognition + confirmed apply | `ScreenTranscriptionTermination` + `MaterializeTranscriptionTerminationFeatures`; exact matching, both strands, **linear only**. | Determinism, strand geometry, stale-digest refusal, idempotence, existing-annotation dedupe, no-hit row, read-only classification, glossary/capability/MCP parity. |
| **C** Bacterial intrinsic candidates | Versioned pure hairpin + U-tract heuristic; components named as heuristic, never ΔG; ViennaRNA optional and outside the hashed body. | Synthetic positive/negative, ambiguity, work caps, stable ordering, **digest stable with and without RNAfold installed**. |
| **D** One mammalian Pol II profile | Explicit subject binding; signal, context, window, cleavage and termination kept as separate geometry; unassigned candidates on plasmid scans. | Motif-only negative, multiple candidates, reverse strand, annotated-UTR vs candidate-region distinction, missing-flank refusal. |
| **E** Presentation + parity completion | One evidence model across report, map overlay, exports. Deferred list unchanged (approximate matching, Pol III, other hosts, read-end import). | Adapter parity totals, export/undo, no new biology in any adapter. |

## 5. Remaining decisions (yours, not mine)

1. Feature-library ownership — generic primitive, or termination as consumer.
2. Command noun and whether it groups under an existing one.
3. Whether external folding is permitted in v1 at all. My advice: no.
4. Catalog home — protocol-side catalog, or helper-construct semantics extended with sequences.
5. Who supplies the curated elements and their redistribution clearance.
6. Pol III oligo-T in or out of the first biological release. My advice: out.
7. The doc's "Claude feedback: unavailable" statement at [lines 19-31](docs/transcription_termination_annotation_plan.md:19) is now stale if you adopt any of this, and by the doc's own standard it must record which items were accepted or rejected.

Want me to write this into `docs/transcription_termination_annotation_plan.md` (that needs edit authorization, which this task did not grant), or leave the review here in the transcript?


---

## Owner-supplied continuation (sections 6–11)

Received 2026-11-16; retained verbatim.

Continuing — the parts Codex actually needs to act on: the exact parity surface and the ownership question the plan leaves open. Let me verify those seams.

Continuing the review — sections 6-11 below extend, not replace, the earlier critique.

## 6. The protocol-ownership question is already settled — stop hedging

The plan says "choose the existing protocol ownership pattern before placing new records; prefer a small typed protocol module rather than growing a central file" ([plan lines 103-105](docs/transcription_termination_annotation_plan.md:103)). That inspection is done, and the answer is unambiguous: `crates/gentle-protocol/src/lib.rs` declares **34 dedicated `pub mod` modules**, with `pub mod tata_boxes;` at line 44. Records are referenced as `gentle_protocol::tata_boxes::TataBoxScreenRequest` from `src/engine.rs:4753`. `src/engine/protocol.rs` holds exactly **one** line for the precedent — the persistence field at [src/engine/protocol.rs:5838](src/engine/protocol.rs:5838).

So for Codex: a new `crates/gentle-protocol/src/transcription_termination.rs`, one `pub mod` line, one persistence field. `src/engine/protocol.rs` is **not** a seam for the new records, despite [docs/biological_extension_guide.md:103](docs/biological_extension_guide.md:103) listing it as "stable request/result records" — the guide's table predates the crate extraction and is stale on this point. Worth a one-line correction to the guide separately; don't let Codex rediscover this.

Required conventions, copied verbatim from the precedent: a `*_SCHEMA` const and a `*_NON_CLAIM` const ([tata_boxes.rs:7-8](crates/gentle-protocol/src/tata_boxes.rs:7)); `#[serde(default, deny_unknown_fields)]` on the request (line 11) — this is what satisfies the guide's "unknown-field rejection" test; an explicit `impl Default` with the policy defaults visible (lines 30-48); `content_sha256` documented as "hash of this report with this field empty" (line 114); and the materialize request carrying only `screen` + `expected_report_sha256` + `row_ids` (lines 142-146).

## 7. The parity surface Codex must budget for

Phase 4 treats this as presentation polish. Here is what the precedent actually required, so the estimate is honest:

| Surface | File and anchor |
|---|---|
| Protocol records | `crates/gentle-protocol/src/tata_boxes.rs` (146 lines) + `pub mod` at `lib.rs:44` |
| Operation variants | [src/engine.rs:4752](src/engine.rs:4752), 4757 |
| Read-only classification | [src/engine.rs:10504](src/engine.rs:10504) — screen in, materialize **out** |
| Operation-name registry | `crates/gentle-protocol/src/lib.rs` ≈5762 |
| Handler dispatch | [src/engine/ops/operation_handlers.rs:41558](src/engine/ops/operation_handlers.rs:41558), 41572 |
| Analysis implementation | `src/engine/analysis/tata_boxes.rs` (1111 lines, ~370 of them tests) |
| Report persistence | [src/engine/protocol.rs:5838](src/engine/protocol.rs:5838) |
| Capability descriptors ×4 | [src/engine_shell.rs:22513-22516](src/engine_shell.rs:22513) — registered **twice per operation**, once under the command path and once under the operation name, each with `mutating` + `requires_confirmation` |
| Readiness preconditions | [src/engine_shell.rs:32624](src/engine_shell.rs:32624) |
| Shell parser + subcommand error text | [src/engine_shell.rs:44232](src/engine_shell.rs:44232), 44390-44402 |
| Glossary | [docs/glossary.json:2650](docs/glossary.json:2650), 2693 |
| MCP | generic typed `op` route only — **no dedicated tool needed** ([src/mcp_server.rs:6712](src/mcp_server.rs:6712)) |
| GUI | `src/main_area_dna/tata_box_ui.rs` (543 lines) + command palette [src/app.rs:5788](src/app.rs:5788) + [src/app.rs:9088](src/app.rs:9088) + `src/app/subject_selection.rs:190`, 308 |

**The gate the plan never mentions:** `docs/gui_cli_mcp_parity.md` is a **generated, byte-exact-checked** Markdown projection. `run_parity_matrix_check_mode` ([src/bin/gentle_examples_docs.rs:730](src/bin/gentle_examples_docs.rs:730)) fails with "is stale; run `scripts/regenerate_parity_matrix.sh`", and the precedent occupies four rows (245, 251, 662, 663). Adding two operations adds four rows. Per AGENTS.md, byte-exact generated-file freshness checks — "including Markdown projections" — need scoped checkout policy and LF/CRLF coverage in the fast `scripts.test_tutorial_checkouts` gate. Codex must regenerate this file in the same change. It needs **no** `docs/parity_matrix_overrides.json` entry: the precedent has none, so the generic routes suffice.

One useful relief: the generated matrix *records* prominence rather than demanding it. The precedent's engine-operation rows read `shell-only` across every adapter while its glossary-command rows read `prominent`. Combined with [docs/biological_extension_guide.md:168](docs/biological_extension_guide.md:168) ("Prominent buttons or dedicated MCP tools are optional. Reachability through the shared operation is not optional"), the 543-line GUI module is genuinely deferrable to a later slice — the matrix will simply say `shell-only` and be accurate. Say that explicitly in the plan so Codex doesn't treat the GUI as a gate on Phase B.

## 8. Acceptance rows missing from the plan

The precedent's own tests supply four rows the plan's verification section lacks:

1. **MCP refuses the mutating op without `confirm`, and the state file is byte-identical afterwards.** [src/mcp_server.rs:6730-6736](src/mcp_server.rs:6730) asserts `isError == true` *and* `std::fs::read(&path) == before`. This is the concrete DEC-013 test; "selected undoable writeback" in the plan does not cover it.
2. **Glossary entries are parse-verified generically.** `glossary_cli_usage_smoke_commands_parse` and `glossary_cli_usage_flags_parse_one_by_one` ([src/engine_shell/tests.rs:3018](src/engine_shell/tests.rs:3018), 3043) iterate every glossary command, so a glossary entry without a working parser breaks the suite automatically. Good news; just know it will fire.
3. **Parity-matrix freshness**, per section 7.
4. **Digest stability with and without RNAfold installed** — the test that proves finding (c) from my earlier critique was handled rather than merely acknowledged.

Two more gaps specific to this biology. The dedupe test must use a feature that the existing code actually produces — kind `terminator` with `regulatory_class=terminator` as written at [src/engine.rs:22463](src/engine.rs:22463) and parsed back at 22886-22909 — otherwise it proves nothing about real GenBank input. And the readiness precondition needs a decision: the precedent registers `Some(vec![])` (no required project facts) because EPD is optional. If the curated catalog ships bundled, empty is right; if it is user-supplied, it needs a readiness fact and a structured "catalog unavailable" refusal. The plan is silent, and this is the difference between a screen that works out of the box and one that always refuses.

## 9. What to cut from the first slice

- Origin-crossing search → structured refusal, per section 1d.
- Approximate/bounded matching → already deferred; keep it deferred and delete the identity/coverage/indel fields from the v1 schema so they cannot be half-implemented.
- The GUI module → Phase E, with the matrix honestly reading `shell-only`.
- External folding → out of v1 entirely (my recommendation on decision 3 below).
- "Separate scientific benchmark against independently curated examples" → relabel owner-blocked, not acceptance; with no catalog in the tree and no web access in the implementing session it cannot be met.

That leaves a first slice of roughly: one protocol module, two operation variants, one analysis module with exact matching, parser + glossary + descriptors, regenerated parity matrix, and the deterministic test set. Comparable to the precedent, which is the right size.

## 10. Codex-ready first-slice brief

> Add exact curated-element recognition for transcription-termination and 3′-processing elements, plus a separately confirmed apply step, as one shared engine decision. Classify as **class C then B**: the missing primitive is a provenance-bound curated element catalog; the composition is screen-then-materialize, which already exists as a pattern.
>
> Copy `crates/gentle-protocol/src/tata_boxes.rs` and `src/engine/analysis/tata_boxes.rs` end to end: schema const, non-claim const, `deny_unknown_fields` request with explicit `Default`, report with `content_sha256`, materialize taking `expected_report_sha256` + `row_ids`, per-row idempotence key stored as a feature qualifier, `prepare_sequence` after the feature push. Do **not** route through `construct_reasoning_annotation_candidate_*` in `src/engine.rs` — that path requires `source_kind == "generated_annotation"` plus `EvidenceScope::SequenceSpan` and is graph-specific.
>
> Scope: exact matching only; both strands; **linear spans only** — refuse origin-crossing search with a structured error, as `src/engine/analysis/cryptic_splicing.rs:511` does. Subject is one `seq_id` + interval, or one explicitly resolved transcript/CDS id; refuse ambiguity rather than choosing an isoform. Whole-sequence scans return unassigned candidates with no transcript association. Keep element match, predicted motif, cleavage window, annotated RNA end and polymerase termination region as five distinct geometries; none implies another, and none implies efficiency. Catalog entries carry redistribution and host/mechanism provenance in the vocabulary of `GenomeCatalogEntry.redistribution_status`. Fixtures hand-crafted and non-biological, documented per AGENTS.md.
>
> Deliver with: parser + glossary + four capability descriptors + readiness preconditions, regenerated `docs/gui_cli_mcp_parity.md`, protocol docs, and the acceptance rows in section 8. No GUI, no approximate matching, no folding backend, no prediction of any kind in this slice. `cargo check -q` before handoff.

## 11. Remaining decisions, with my recommendation

| # | Decision | Recommendation |
|---|---|---|
| 1 | Feature-library ownership: generic primitive, or termination as consumer | **Generic primitive**, termination as its first catalog — it is the parked entry's own contract and avoids a second matcher |
| 2 | Command noun | `terminators screen` / `terminators materialize`; don't nest under `promoters` |
| 3 | External folding in v1 | **No.** It breaks digest reproducibility and `src/rna_structure.rs` cannot take a DNA window anyway |
| 4 | Catalog home | New protocol-side catalog. `HelperConstructComponent` carries no sequence field, so extending it means changing a shipped record for an unrelated consumer |
| 5 | Who supplies elements and clearance | Owner, explicitly. Not an agent task |
| 6 | Pol III oligo-T in first biological release | Out |
| 7 | Catalog bundled or user-supplied | Decide before Phase A ends; it determines the readiness fact and the refusal path |
| 8 | The doc's "Claude feedback: unavailable" note at [lines 19-31](docs/transcription_termination_annotation_plan.md:19) | Now stale. By the plan's own standard it must record which of these items were accepted and which rejected, attributed to this review |

Nothing was edited and no code ran; the working tree is as I found it, with your unrelated modified and untracked files intact. The one correction worth making beyond the plan itself is [docs/biological_extension_guide.md:103](docs/biological_extension_guide.md:103), which still points new records at `src/engine/protocol.rs` after the protocol crate extraction moved them — stale guidance that will mislead the next agent too.
