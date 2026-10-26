# Primer-Anchored Read Trace, Revised Phase 0/1

Implementation baseline: `adacb2bc`, 2026-10-26. This revises the owner-forwarded
Claude prompt `primer_trace_phase01_codex_prompt_20261026.md` in the primary
checkout. That untracked original is not modified or treated as acceptance.
The revisions below are Codex's source-grounded changes, not a separate
Claude-reviewed or scientifically accepted contract.

## Revisions Before Implementation

- Use the explicitly supplied oligo and components, not an assumed 44-mer.
  The retained September review supplies newer annealing cores, not proof that
  those complete tailed oligos were ordered. Require a terminal target core;
  allow a core-only oligo. Sequences are canonical DNA, at most 64 bases.
- Recover observed intervals and core edits from the approximate alignment,
  including insertions/deletions. Report primer/read orientation independently
  of transcript orientation. The payload starts at the observed primer end.
- A full hit suppresses only its own redundant component hits. Keep other
  primers, disjoint occurrences and equivalent template matches. Collapse
  overlapping approximate representations deterministically, not exact repeats.
- Validate core plus adjacent payload on one contiguous template interval.
  A seed only nominates a site; it cannot establish template concordance.
  A short local probe is not whole-read, isoform or mechanism evidence.
- A short/ambiguous payload, missing core site or insufficient template flank
  is unresolved, not boundary divergence. Define edit tolerance and exact-prefix
  concordance explicitly. Intentional antisense priming is not an artifact.
- Reuse packed-word postings once per bound reference. Seeds only nominate
  bounded candidate intervals, including indel drift; verify the whole probe.
  Bind seed length/posting limits as policy/budgets and require a seed length
  compatible with the probe's unit-edit budget. Do not rescan full references
  with an approximate aligner on every unresolved read.
- No-anchor fields are absent/null, not invented terminal/single classifications.
  Handle observations have no attributed primer/target. Common handles are
  searched once. Per-read classes and all-hit histograms have separate units.
- Expose explicit read-only inline-draft sealing via panel validation. Ordinary
  validation and execution remain strict; a nonempty mismatched digest is never
  silently replaced. Agents can obtain a real digest without OS tooling.
- Preserve FASTQ quality in the existing streaming record; existing nomination
  consumers continue to ignore it. No legacy serialization or policy changes.
- Add a named artifact-only managed execution policy for this operation, with
  its own typed result binding. Do not infer persistence from registry ReadOnly
  or generalize the still-proposed lifecycle framework.

## Implementation Order

1. Versioned primer panel, explicit policy/budgets, request, hit/row, summary,
   report and result binding; synthetic orientation/alignment tests.
2. Shared streaming intake, EOF checksums, cancellation/marker errors, staged
   publication and bounded inspection. Existing RNA-screen bytes stay unchanged.
3. Shared operations, shell routes, generic MCP/agent operation parity and
   external-write confirmation; no GUI or private-data execution.
4. Protocol/CLI/tutorial explanation and one consolidated Glen handoff.

The axis model and biological defaults remain proposed pending owner review.
The poly(T) observation is configurable and does not authenticate kit identity.
Payload FASTA defaults off. No segmentation graph, junction classifier, aligner
adapter, mechanism call or automatic input discovery is included.

## Verification Boundary

Use synthetic sequences only. No local Rust build/test/check is authorized by
the forwarded task. Author Rust regressions, run scoped rustfmt and available
Python source/document/checkout checks, and disclose execution gaps. Glen/CI
must execute the focused tests and compiler checks at the final integrated SHA.
Real oligos, private reads, barcode controls, throughput, native platforms and
scientific acceptance remain separate and unclaimed. No push is authorized.

## Implementation And Local Evidence

Phase 0/1 source is implemented on `gentle_rs_2_main` based on `adacb2bc`.
The policy/biological decision remains proposed. Local `main` advanced to
`4acee734` during this session (RNA-screen profiling/summary work); this diff
was not rebased or compiled against that newer source. Integrate before the
single consolidated Glen/CI run; do not inherit either branch's old verdicts.

Changed areas: new protocol/matching/publication modules; the existing FASTQ
record's quality field; operation/result wiring and required constructor
defaults; explicit managed artifact/observation/control classification; shared
shell, generic MCP command-path metadata and approval guards; glossary/parity
projection, protocol/CLI/decision/roadmap, tutorial 10.04 and the consolidated
Glen checklist. No new GUI control, dependency, private fixture or index change
to the RNA screen. The pre-existing untracked `output/` remains untouched.

Checks executed on macOS, 2026-10-26:

```sh
python3 -m unittest scripts.test_primer_trace_contract scripts.test_rna_panel_tutorial scripts.test_capability_maturity scripts.test_specialist_tutorials scripts.test_tutorial_gap_lessons scripts.test_tutorial_checkouts -q
python3 scripts/check_documentation_status.py
python3 scripts/generate_capability_maturity.py --check
rustfmt --edition 2024 --check --config skip_children=true src/engine/protocol/primer_trace.rs src/engine/analysis/rna_panel/primer_trace.rs src/engine/analysis/rna_panel/primer_trace/matching.rs src/engine/analysis/rna_panel/primer_trace/tests.rs
git diff --check
```

The combined Python run discovered 167 tests: 162 passed, five optional
binary-dependent checks skipped. The new tutorial-preparation test really
executes its synthetic Python snippet in a disposable directory and checks
the panel/input hashes and orientation arithmetic; it does not execute GENtle.
Later focused source checks also pass. Documentation status and maturity
freshness pass. `generate_capability_maturity.py --write` refreshed only the
four required glossary-input digest projections, not maturity attestations.
The matrix's added rows/counts were reconciled from the canonical registry
rules and checked for count consistency; the Rust generator's byte-exact check
remains external work. Existing LF policy/checkout tests already cover that
matrix artifact; no new byte-exact freshness gate or committed data fixture
was introduced.

Rustfmt parsed all changed tracked Rust files without rewriting them; the four
new Rust modules pass scoped formatting. This is **not** Rust type-checking.
Twenty-one new `primer_trace` Rust regressions are authored but unrun. No
Cargo build/test/check, native Windows, full-suite, scientific, performance or
release acceptance is claimed. `Cargo.lock` is unchanged. Session-close reports
four OK, two warnings (intentional dirty/generated files and manual plan-fidelity
reminder), zero failures. The dedicated named artifact policy and explicit
read-only draft sealing are intentional, documented parts of this slice; the
general lifecycle proposal remains unimplemented.
