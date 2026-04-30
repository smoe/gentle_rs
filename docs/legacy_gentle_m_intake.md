# Legacy `gentle-m` Intake Priorities

Purpose: capture which legacy `GENtle-persons/gentle-m` capabilities are still
worth carrying into the Rust rewrite, in what order, and which upstream source
files are the best seed material.

Inspected upstream snapshot:

- Repository: `https://github.com/GENtle-persons/gentle-m`
- Commit: `ad80adfd837b2eebc1bc0e678d49fb88e39f7ba6`
- Main references reviewed:
  - `docs/manual.adoc`
  - `src/TSequencingAssistantDialog.cpp`
  - `src/ABItype.cpp`
  - `src/SCFtype.cpp`
  - `src/MiscDialogs.cpp`
  - `src/TSilmutDialog.cpp`
  - `src/AutoAnnotate.cpp`
  - `src/TLigationDialog.cpp`

Core rule for intake:

- Port biology/algorithmic logic only.
- Do not port wxWidgets dialog structure, direct database plumbing, or
  adapter-local behavior that bypasses shared engine contracts.
- Every adopted behavior must become deterministic shared-engine state,
  reports, or explanation artifacts that GUI/CLI/JS/Lua can all reuse.

## Priority Order

### 1. Sequencing confirmation and raw-trace intake

Why this is first:

- It is the clearest legacy capability that improves user trust in whether a
  clone is actually correct.
- It directly matches the current Rust-roadmap gap for construct validation
  from read evidence.
- The legacy code contains both parser seeds and task framing:
  - `src/ABItype.cpp`
  - `src/SCFtype.cpp`
  - `src/TSequencingAssistantDialog.cpp`
  - `docs/manual.adoc` sequencing and alignment sections

Recommended Rust delivery slices:

1. Sanger/amplicon called-sequence confirmation first.
   - Expected construct + one or more read sequences.
   - Deterministic confirmation report with target-level status:
     `confirmed`, `contradicted`, `insufficient_evidence`.
   - Reuse current `AlignSequences` / `SequenceAlignmentReport` machinery.
2. Raw ABI/AB1/SCF import next.
   - Parse traces and called bases into an engine-owned trace/read record.
   - Non-mutating inspection only at first.
3. Trace-aware confirmation after the report contract is stable.
   - Connect called-base evidence and trace support into the same confirmation
     report rather than building a standalone viewer-first subsystem.

Detailed implementation plan:

- `docs/sequencing_confirmation_plan.md`

What to mine from legacy code:

- ABI/AB1 record parsing approach from `src/ABItype.cpp`
- SCF v3 parsing approach from `src/SCFtype.cpp`
- Reverse-complement/orientation sanity check idea from
  `src/TSequencingAssistantDialog.cpp`

What not to port:

- wx-based ABI viewer UI
- alignment-window-centric workflow where confirmation exists only as an
  informal human interpretation step

### 2. Sequencing-primer suggestion overlay

Why this is second:

- It is small, bench-friendly, and synergistic with sequencing confirmation.
- It helps users plan and explain read coverage before or after confirmation.
- The legacy implementation is simple enough to port cleanly without dragging
  UI baggage along.

Legacy reference:

- `src/MiscDialogs.cpp`
  - `findBestMatch`
  - `matchToVector`
  - `addSequencingPrimer`
- `docs/manual.adoc` "Sequencing Primers"

Rust intake shape:

- Engine operation or helper that scans a construct against a primer set using
  exact 3' anneal length thresholds.
- Emit candidate sequencing-primer features/annotations with direction,
  annealing span, and original primer sequence retained as structured metadata.
- Surface in GUI as optional overlays and in CLI as machine-readable results.

### 3. Silent-mutation verification helper

Why this is third:

- It is a genuinely useful bench-biology verification trick.
- It complements PCR and confirmation workflows instead of competing with the
  newer cloning-specialist architecture.
- The old logic is algorithmic and portable.

Legacy reference:

- `src/TSilmutDialog.cpp`
- `docs/manual.adoc` "Silent Mutation"

Rust intake shape:

- Engine contract that enumerates synonymous edits within a selected coding
  window that introduce a chosen or candidate restriction site.
- Return a deterministic report with:
  - nucleotide edits
  - amino-acid preservation proof
  - before/after restriction-site counts
  - predicted fragment lengths for verification digest

### 4. Auto-annotation library scan

Why this is fourth:

- Still valuable, but less urgent than trust-building confirmation features.
- The Rust roadmap already wants this, and the legacy code confirms user value.
- The old implementation is useful as heuristic inspiration, but it needs a
  more explicit evidence model than the legacy "annotate if matched" flow.

Legacy reference:

- `src/AutoAnnotate.cpp`
- `docs/manual.adoc` "Automatic annotation"

Rust intake shape:

- Shared engine scan against curated vector/feature libraries.
- Output confidence/overlap diagnostics as a report first.
- Separate "apply these annotations" from "scan and suggest them".

### 5. Legacy ligation assistant behavior

Why this is lower:

- Useful historically, but the Rust rewrite already has a stronger typed
  cloning-routine direction than the old generic ligation dialog.
- We should not regress into a loose combinatorial UI-first model.

Legacy reference:

- `src/TLigationDialog.cpp`
- `docs/manual.adoc` "Ligation"

Rust intake shape:

- Only lift sticky-end/blunt-end product-enumeration ideas if they fit shared
  engine routine contracts.
- Do not recreate the legacy dialog as-is.

## Recommended Execution Sequence In This Repo

1. Implement sequencing-confirmation report contracts and CLI inspection/export.
2. Add GUI confirmation specialist using the same report, not adapter-local
   logic.
3. Add ABI/AB1/SCF import and trace inspection as engine-owned evidence intake.
4. Add sequencing-primer suggestion overlays tied to confirmation coverage.
5. Add silent-mutation verification helper for PCR confirmation workflows.
6. Add auto-annotation library scan with report-first/apply-second semantics.

## Explicit Non-Goals For Intake

- Porting old wx dialogs one-to-one
- Recreating legacy database UX
- Recreating legacy BLAST web wrappers
- Keeping functionality that only exists as display state rather than
  engine-owned records
