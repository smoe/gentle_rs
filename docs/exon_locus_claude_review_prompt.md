# Claude Review: General Exon-Locus Pileups

Please review read-only in `/Users/u005069/GitHub/gentle_rs`. Do not edit files,
run Rust/builds/tests/shell tools, contact providers, process private data,
push, dispatch CI, or interpret this as scientific/release acceptance. Leave
`paper`, private sequencing data, `output/`, unrelated prompt files and stashes
alone. Read current working-tree changes, not only HEAD; implementation may not
yet be committed. Report the source boundary you actually inspected.

## Task And Original Codex Plan

Owner authorized implementing Glen's general exon-locus extension after the
transcript-bound oriented-pileup pass at `7fa6b15a`. The canonical prompt is:
<https://gentle.functional.domains/openclaw/canvas/codex-exon-locus-prompt-20261109-d72152e0/report.md>
SHA-256 `d72152e0197392e2f820312e5a338e860e56cc7a5737f94ef057fa610d79df6d`.
Do not fetch inaccessible URLs or bypass an allowlist; the local plan and tests
record the contract. TP73/TP73-AS are examples, never fixed product identifiers.

Original four-part Codex plan, retained for critique:

1. Pure versioned locus/read/assignment records before gene filtering; stable
   sample/unit IDs; all hypotheses and seven classes; deterministic splice-first
   ranking, explicit conflicts and unresolved ties. No BAM-only RNA-origin call.
2. One authenticated service and two projections, preserving old gene-only
   coverage and v1/v2 wires. Explicit bindings, bounded work/memory, cancellation,
   confirmation and publication; missing annotation keeps genome inspection,
   changed annotation recomputes and bad raw hashes refuse.
3. Fixed genomic tracks for both annotation strands and explicit mature mRNA
   projection of the same unique reads. Shared renderer/native virtualization,
   clips/off-window links, original literal A/T, complete labels and pixel guards.
4. General synthetic Rust/service/real-VM/actual-CLI cases, five formats,
   deterministic exports, LF/CRLF and quoted paths, no-clobber, unchanged project
   and inputs. Record pending native Windows/GUI and scientific boundaries.

## Actual Implementation To Inspect

- `src/observed_transcript_structures/read_ends/locus.rs`,
  `locus/presentation.rs`, `locus/tests.rs`: membership, geometry ranking,
  assignments, unknown/conditional molecule-strand status, pages/counters/scenes.
- `src/observed_transcript_structures/read_ends.rs`, `oriented.rs`, `pileup.rs`:
  factored exon map/read key and preserved older gene-filtered/original-base path.
- `src/rna_read_ends.rs`, `src/rna_read_ends/locus.rs`, `locus_tests.rs`:
  bound loaders, missing annotation, completion verification, same exporter,
  request-bound cancellation, managed result retention and no-clobber publisher.
- `src/command_execution/rna_artifacts.rs`,
  `src/app/rna_panel_ui/read_ends.rs`, `src/engine_shell.rs`:
  observation admission, thin native controls, view staleness and discovery.
- `crates/gentle-render/src/read_end_pileup.rs`: annotation tracks, hypothesis
  sublanes, actual checked dimensions, shared native/full-export drawing.
- Updated `read_end_` tests in JS/Lua/MCP; synthetic
  `scripts/prepare_exon_locus_rehearsal.py` and
  `tests/test_rna_read_end_locus.py`; provenance in
  `tests/README_rna_read_end_metrics.md`.
- `docs/decisions.md` (Exon-Locus Assignments Before Gene Filtering),
  `docs/rna_read_end_metrics.md`, `docs/rna_read_end_integration_plan.md`,
  `docs/glen_next_verification_handoff.md`, CLI/GUI/protocol references.

Evidence: previous `7fa6b15a` pass is scoped to older transcript-bound code,
not this change. The local Python fixture/source checks pass; actual CLI tests
skip without an explicitly source-matched binary. Scoped rustfmt parses Rust,
but no local Rust compilation, VM execution or native acceptance is claimed.
Glen/CI owns those commands in the active handoff. Lockfile/dependencies remain
unchanged. No Claude review has run before this prompt was prepared.

## Requested Critique

Give concise risk-ordered findings with file/line references and concrete
counterexamples. Prioritize defects over feature expansion. Check specifically:

1. Can a shared-exon read, reverse-complement or mirrored locus acquire a false
   unique sense/antisense assignment from BAM direction, endpoint, A/T or ties?
   Can lower evidence outweigh a structural conflict? Are alternative legal
   hypotheses collapsed honestly? Supplementary alone must not imply chimera.
2. Are all candidates, original A/T coordinates/letters, partial/unmapped runs,
   genomic/selected-transcript geometries, IDs, counters, model/rule bindings and
   missing/changed annotation facts retained without double counting or caching?
3. Are work, candidate/geometry/run bounds checked before costly allocation;
   cancellation, confirmation, raw-path/hash checks and no-clobber intact?
   Does the actual-layout raster guard include tracks/sublanes/wrapped labels?
4. Are GUI controls truly thin, stale results honest, and virtualization bounded?
   Do synthetic tests cover prompt 11-21 meaningfully rather than restating code?
   Identify plausible compile/API errors too; do not claim you compiled them.
5. Suggest the smallest justified improvements and a revised minimal repair/
   verification plan. Separate must-fix bugs from optional future work; avoid
   new biology, dependencies, scoring calibration, queues or release expansion.

Present: original-plan critique, findings/evidence, then revised minimal plan.
Explicitly state missing runtime evidence. A no-finding review is scoped,
not a pass for native Windows, real data, physical tails, APA or causality.
