# gRNA Design Spec (Draft v0.1)

This document defines a first implementation-oriented spec for gRNA
design in GENtle, focused on:

1. practical design-constraint filters
2. oligo/export for lab execution
3. a macro/template pattern that fits GENtle's operation/workflow model

Status: implementation draft. Guide-candidate operations are now implemented in
the shared engine/shell layer (`UpsertGuideSet`, `DeleteGuideSet`,
`FilterGuidesPractical`, `GenerateGuideOligos`, `ExportGuideOligos`,
`ExportGuideProtocolText`). Remaining scope is macro/template packaging and
advanced scoring/off-target models.

## 1. Goals and non-goals

Goals (this phase):

- Standardize practical guide-quality filters so users get predictable pass/fail
  behavior and explainable rejections.
- Standardize oligo generation + export for cloning and ordering.
- Keep everything compatible with existing `Operation` + `Workflow` execution,
  so GUI/CLI/JS/Lua can share the same behavior.

Non-goals (later phases):

- Full on-target efficacy model integration.
- Full off-target genome search and ranking.
- Base/prime editor-specific edit window simulation.

## 2. Data model (implemented baseline)

The following structures are proposed as protocol-level JSON payloads.

### 2.1 Guide candidate

```json
{
  "guide_id": "g_0001",
  "seq_id": "target_seq",
  "start_0based": 1234,
  "end_0based_exclusive": 1254,
  "strand": "+",
  "protospacer": "GACCTGTTGACGATGTTCCA",
  "pam": "AGG",
  "nuclease": "SpCas9",
  "cut_offset_from_protospacer_start": 17
}
```

### 2.2 Practical filter config

```json
{
  "gc_min": 0.30,
  "gc_max": 0.70,
  "max_homopolymer_run": 4,
  "max_homopolymer_run_per_base": {
    "A": 5,
    "C": 4,
    "G": 4,
    "T": 4
  },
  "reject_ambiguous_bases": true,
  "avoid_u6_terminator_tttt": true,
  "u6_terminator_window": "spacer_plus_tail",
  "max_dinucleotide_repeat_units": 4,
  "forbidden_motifs": ["GAATTC", "AAGCTT"],
  "required_5prime_base": "G",
  "allow_5prime_g_extension": true
}
```

### 2.3 Filter result (per guide)

```json
{
  "guide_id": "g_0001",
  "passed": false,
  "reasons": [
    {
      "code": "u6_terminator_t4",
      "message": "Contains TTTT in spacer/tail window"
    }
  ],
  "metrics": {
    "gc_fraction": 0.25,
    "max_homopolymer_run": 4,
    "max_dinucleotide_repeat_units": 2
  }
}
```

## 3. Practical design-constraint filters (spec)

Filters are deterministic and evaluated in a fixed order. Each filter produces:

- `pass/fail`
- machine-readable `code`
- human-readable `message`
- optional measurement in `metrics`

Hard filters reject a guide. Soft filters add warnings only (future option).

### 3.1 Required core filters

1. `gc_bounds`
- Compute GC over protospacer only.
- Reject if outside `[gc_min, gc_max]`.

2. `homopolymer_run`
- Compute longest run of identical base in protospacer.
- Reject if above global `max_homopolymer_run`.
- If `max_homopolymer_run_per_base` is set, per-base cap overrides global cap.

3. `ambiguous_base_check`
- Reject if protospacer contains non-ACGT when `reject_ambiguous_bases=true`.

4. `u6_terminator_tttt`
- For U6-driven sgRNA use, reject guides with `TTTT` in configured window.
- `u6_terminator_window`:
  - `spacer_only`
  - `spacer_plus_tail` (default)
- `tail` here means sequence added immediately downstream of spacer in the
  sgRNA expression transcript context.

5. `dinucleotide_repeat`
- Reject when max repeat units for any dinucleotide exceed
  `max_dinucleotide_repeat_units`.
- Example: `ATATATAT` has 4 repeat units of `AT`.

6. `forbidden_motif`
- Reject if protospacer contains any configured motif (exact match, ACGT only).

7. `required_5prime_base`
- If set, spacer 5' base must equal requested base.
- If `allow_5prime_g_extension=true` and required base is `G`, operation may
  mark candidate as "rescuable by 5' G extension" instead of hard reject.

### 3.2 Reason codes (initial set)

- `gc_too_low`
- `gc_too_high`
- `homopolymer_run_exceeded`
- `contains_ambiguous_base`
- `u6_terminator_t4`
- `dinucleotide_repeat_exceeded`
- `forbidden_motif_present`
- `required_5prime_base_missing`

### 3.3 Engine operation (implemented)

`FilterGuidesPractical`

```json
{
  "FilterGuidesPractical": {
    "guide_set_id": "set_001",
    "config": { "gc_min": 0.30, "gc_max": 0.70 }
  }
}
```

Behavior:

- Reads guide candidates from persisted guide-set metadata.
- Writes per-guide filter report into metadata.
- Optionally writes a reduced passed-only guide set.

## 4. Oligo/export layer (spec)

## 4.1 Oligo generation model

Guide-to-oligo generation should be template-driven.

Template object:

```json
{
  "template_id": "lenti_bsmbi_u6_default",
  "description": "U6 sgRNA cloning oligos with BsmBI overhangs",
  "forward_prefix": "CACC",
  "forward_suffix": "",
  "reverse_prefix": "AAAC",
  "reverse_suffix": "C",
  "reverse_uses_reverse_complement_of_spacer": true,
  "uppercase_output": true
}
```

Per guide output:

```json
{
  "guide_id": "g_0001",
  "forward_oligo": "CACCGACCTGTTGACGATGTTCCA",
  "reverse_oligo": "AAACTGGAACATCGTCAACAGGTC",
  "notes": ["5' G extension applied"]
}
```

### 4.2 Export targets

Minimum export targets:

1. `csv_table`
- One row per guide with filter metrics and oligos.

2. `plate_csv`
- 96/384 layout with `well`, `guide_id`, `forward_oligo`, `reverse_oligo`.
- Deterministic fill order: sorted by guide rank, then guide id.

3. `fasta`
- FASTA entries for oligos (`>guide_id|forward`, `>guide_id|reverse`).

4. `protocol_txt`
- Human-readable step list for technical assistant execution.
- Includes template, annealing setup, cloning overhang assumptions, and QC notes.

### 4.3 Engine operations (implemented)

1. `GenerateGuideOligos`

```json
{
  "GenerateGuideOligos": {
    "guide_set_id": "set_001",
    "template_id": "lenti_bsmbi_u6_default",
    "apply_5prime_g_extension": true
  }
}
```

2. `ExportGuideOligos`

```json
{
  "ExportGuideOligos": {
    "guide_set_id": "set_001",
    "format": "csv_table",
    "path": "exports/guides.set_001.csv"
  }
}
```

3. `ExportGuideProtocolText`

```json
{
  "ExportGuideProtocolText": {
    "guide_set_id": "set_001",
    "path": "exports/guides.set_001.protocol.txt",
    "include_qc_checklist": true
  }
}
```

## 5. Macro/function template for GENtle

The design should treat guide design as a reusable operation pipeline, not a
special-case UI feature.

## 5.1 Pipeline template concept

Define reusable workflow templates that compile to standard `Workflow` JSON:

1. candidate generation
2. practical filtering
3. ranking/selection
4. oligo generation
5. export

This aligns with current shared execution:

- `Operation` for atomic logic
- `Workflow` for composition
- shell/CLI/GUI all call the same engine core

## 5.2 Macro file shape (proposed)

```json
{
  "schema": "gentle.macro.v1",
  "macro_id": "crispr_u6_basic",
  "parameters": {
    "target_seq_id": "tp73",
    "guide_set_id": "tp73_guides",
    "gc_min": 0.35,
    "gc_max": 0.70,
    "max_guides": 50
  },
  "workflow": {
    "run_id": "macro:${macro_id}:${target_seq_id}",
    "ops": [
      { "FindGuides": { "seq_id": "${target_seq_id}", "nuclease": "SpCas9" } },
      { "FilterGuidesPractical": { "guide_set_id": "${guide_set_id}", "config": { "gc_min": "${gc_min}", "gc_max": "${gc_max}" } } },
      { "GenerateGuideOligos": { "guide_set_id": "${guide_set_id}", "template_id": "lenti_bsmbi_u6_default" } }
    ]
  }
}
```

## 5.3 Shell UX (proposed)

Add shell wrappers that still route through `Workflow`:

- `macro validate <macro-json-or-@file>`
- `macro run <macro-json-or-@file> [key=value ...]`
- `macro list` (future, if stored in project metadata)

The shell should print:

- resolved parameters
- expanded workflow JSON preview
- operation results with deterministic op ids

## 6. Determinism and provenance requirements

All guide/oligo operations should log:

- filter config snapshot
- template id + concrete oligo prefix/suffix values
- selected guide ids and ranks
- exported artifact paths and checksums (if practical)

This keeps guide design reproducible across GUI/CLI/shell.

## 7. Remaining implementation order (recommended)

1. guide-design macro/template packaging over shared shell/workflow layer
2. advanced ranking hooks (on-target efficacy model integration points)
3. off-target model/search integration
4. GUI panel focused on guide-design operation composition

## 8. Acceptance criteria for phase 1

- Same guide/filter result from GUI, CLI, JS, Lua for identical inputs.
- `csv_table` export round-trips in tests (row count and oligo sequences).
- `protocol_txt` export contains deterministic step numbering and guide ids.
- Failures return structured `ErrorCode` + actionable message.
