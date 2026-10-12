# 04.09 - Check One Primer Pair Against Explicit References

**Purpose:** decide what one primer pair's specificity evidence actually supports,
without confusing a completed search, a recent report or a missing database with
a biological pass. This is a shared-shell walkthrough, not a new GUI specialist
window or a native GUI acceptance record.

## Why This Matters At The Bench

Imagine that a pair amplifies the intended transcript in your loaded locus. That
does not yet answer two different questions:

1. Could it amplify another mature transcript in the prepared cDNA collection?
2. Could residual genomic DNA produce a compatible product?

A junction-spanning primer may reduce genomic carryover while still matching
another transcript. Conversely, a pair that separates two isoforms in the local
matrix may have a compatible product elsewhere in the transcriptome. Keep these
questions separate. A favorable result is conditional on the reference sequences,
annotation and thresholds used; it does not predict laboratory performance.

| Feature | Why it is important | What you inspect |
| --- | --- | --- |
| Exact assembly and release binding | `chr1` in two assemblies is not the same coordinate system. | Reference ID, assembly/release, index kind, content and annotation fingerprints. |
| Caller-provided mapping per reference | A different assembly or cDNA collection needs its own target geometry; copying genomic coordinates can misclassify an off-target as intended. | Subject IDs and one-based inclusive product/binding ranges; mappings are labelled caller-provided, not independently established liftover. |
| Full oligo and annealing/tail boundary | An adapter is not the sequence expected to anneal in the initial PCR. | Full oligo, annealing segment, tail sequence and tail length; dimer/oligo QC remains a separate gate. |
| Required versus optional references | An unavailable essential reference is missing evidence, not a negative search. | Each reference's requirement, availability and verdict; at least one must be required. |
| Complete, retained search evidence | A process exit or truncated hit list does not prove specificity. | Both searches, exit codes, completeness, raw TSV sizes/hashes, policy and intended-target binding. |
| Separate genomic/transcriptome verdicts | Genomic carryover and transcript cross-amplification have different experimental consequences. | `genomic` and `transcriptome` independently. |
| Exact panel selection receipts | A newer passing report for another reference must not erase a required failure. | The selected reference/policy/pair/report hashes for each existing panel dimension, and replaced selections retained as history. |
| Immutable standalone summaries | Results must remain inspectable after a reference is updated. | Stable summary ID and evidence snapshots; historical display is not a current-resource revalidation. |
| Shared engine routes and approval | GUI, scripts and agents must not invent different biological interpretations or silently write evidence. | The same typed operations; approve file-writing preparation/import explicitly. |

## 1. Establish Your Starting State

Use an existing saved primer-design report and already prepared, validated local
indexes. Inspect the report with `primers show-report REPORT_ID`. Decide which
pair to assess: rank is one-based; index is zero-based. Supply **one or the other**,
not both. This workflow assesses one pair, not every assay in a panel.

Choose explicit catalog IDs for each reference, its actual `genomic_dna` or
`transcriptome_cdna` kind and whether it is required. Preparation does not download,
index or search anything. A missing required index stops before bundle files are
written; a missing optional index remains an explicit unassessed row. Wrong kinds,
wildcards, duplicate resolved references and more than eight references are refused.

The new routes work in the GUI Shell on its current project. Terminal commands
below use a separately saved state file; a CLI does not inherit unsaved GUI state.
Use a disposable copy for learning. The existing single-reference panel workflow
is still available; this standalone workflow does not replace it.

## 2. Make The Reference Request Explicit

Save this as `specificity-request.json`, substituting your actual saved report and
catalog IDs. The placeholders are deliberate: they are not real references or
validated biological examples.

```json
{
  "schema": "gentle.primer_pair_multi_reference_request.v1",
  "pair": {
    "kind": "saved_pair",
    "primer_report_id": "YOUR_SAVED_PRIMER_REPORT",
    "pair_rank": 1
  },
  "policy": {},
  "references": [
    {
      "target_genome_id": "YOUR_PREPARED_CDNA_REFERENCE",
      "expected_index_kind": "transcriptome_cdna",
      "required": true
    },
    {
      "target_genome_id": "YOUR_PREPARED_GENOMIC_REFERENCE",
      "expected_index_kind": "genomic_dna",
      "required": true
    },
    {
      "target_genome_id": "YOUR_OPTIONAL_COMPARISON_REFERENCE",
      "expected_index_kind": "genomic_dna",
      "required": false
    }
  ]
}
```

`policy: {}` requests the existing shared defaults; review the effective policy
in the returned handoff before executing. It covers product ceilings, mismatch
criteria, hit retention and full alignment. Do not put a single
`specificity_target_genome_id` in this common policy: the reference list selects
the targets. Alternate catalog/cache roots are explicit top-level `catalog_path`
and `cache_dir` fields, not an invitation to silently choose another database.

### Mapping An Intended Target

GENtle admits saved-source geometry only when its original reference, assembly
and release agree with the searched reference. Otherwise it remains unknown.
No mapping is inferred from the first BLAST hit or a chromosome alias.

For a cDNA reference or a different assembly, provide an `intended_target` within
that reference row **only after independently establishing the subject and
coordinates**. Here is the shape of a wholly invented cDNA mapping; do not copy
its identifiers or coordinates into a real experiment:

```json
{
  "model": "transcript_set",
  "expected_products": [
    {
      "target_space": "transcriptome_cdna",
      "subject_id": "SYNTHETIC_TRANSCRIPT.1",
      "expected_product_range": {"start_1based": 100, "end_1based": 219}
    }
  ],
  "source": "caller_checked_synthetic_example"
}
```

Ranges here are one-based inclusive **subject** coordinates, not positions in
the displayed genomic locus. Genomic mappings use `model: "genomic_interval"`,
the genomic subject, forward/reverse binding ranges and product geometry. A
junction-spanning mapping must explicitly describe its genomic expectations;
do not pretend a contiguous genomic product exists merely because cDNA does.
Preparation stamps caller mappings with the inspected exact reference identity
and labels the evidence caller-provided. That binding is not proof the mapping
itself is biologically correct. With no admissible intended mapping, hits remain
inspectable but cannot earn a specificity pass.

### Assessing Explicit Oligos Instead

Replace `pair` with `kind: "explicit_pair"` and `forward`/`reverse` objects.
Each object gives `role`, `full_sequence`, `annealing_sequence`,
`annealing_length_bp`, `non_annealing_5prime_tail` and
`non_annealing_5prime_tail_bp`. Use normalized uppercase DNA and consistent
boundaries. With no tail, its string is empty and its length is zero. These
details are hash-bound; changing an adapter boundary requires a new handoff.
An explicit pair has no saved-template design provenance, so its per-reference
intended mapping must be supplied separately for a passing interpretation.

## 3. Prepare, Inspect, Then Execute Externally

From the repository root, using an existing state and a **new** bundle directory:

```sh
gentle_cli --state study.gentle.json primers specificity-multi-handoff @specificity-request.json evidence-one-pair
```

Inspect `evidence-one-pair/handoff.json` and `execution_manifest.json`. The parent
binds the effective, sorted reference list, exact pair, current saved-template
snapshot, policy, child handoffs and structured commands. Each available reference
has forward and reverse annealing-segment searches. Unavailable optional references
have no commands. The manifest initially says `pending`, never `pass`.

Use your trusted external runner to execute each child command's `program` and
`args` as an argument vector, not an interpolated shell string. GENtle starts no
BLAST here and adds no scheduler. Review commands before authorizing the runner;
retain its logs as independent process evidence. Do not edit the bound arguments,
move the bundle or copy another run's TSVs into it.

The runner fills each existing manifest command with:

- `state: "completed"` only for a finished successful process, otherwise
  `failed` or `cancelled`; pending work stays `pending`;
- the real `exit_code`;
- `output_size_bytes` for the retained output bytes;
- `output_sha256` as `sha256:` followed by their lowercase hexadecimal SHA-256.

Keep the supplied command ID/digest, path and parent/pair bindings unchanged.
Hash the bytes actually retained, not newline-normalized text. Empty TSVs are
allowed only with genuine successful complete searches and their own size/hash.
A manifest is content-bound evidence, **not an authenticated signature of a
process**. You remain responsible for trustworthy external execution.

## 4. Import And Interpret The Result

Import is a project mutation and requires explicit approval when suggested by
an agent or called through MCP. In your authorized terminal session:

```sh
gentle_cli --state study.gentle.json primers specificity-multi-import evidence-one-pair/handoff.json evidence-one-pair/execution_manifest.json
gentle_cli --state study.gentle.json primers specificity-multi-list
gentle_cli --state study.gentle.json primers specificity-multi-show YOUR_RETURNED_SUMMARY_ID
```

Import checks command identities, exits, raw-byte hashes/sizes, query sequences,
pair/tail interpretation, policy, complete-search limits and source/reference
currency, including annotation fingerprints. It parses those same retained TSV
bytes. Tampering is refused before any new summary is persisted. Stale references
and unfinished searches remain incomplete; valid child evidence is retained for
inspection. Do not interpret a successful import as a biological pass.

| Situation | Dimension verdict | Meaning |
| --- | --- | --- |
| A validated required reference fails | `fail` | A newer passing reference does not cancel this failure. |
| All required references pass and execution is complete | `pass` | Conditional evidence for this declared dimension and policy. |
| Required evidence is stale, missing, unknown or incomplete | `incomplete` | Not biological absence and not a negative experimental result. |
| No references in a dimension | `not_requested` | No claim was sought. |
| Only optional references in a dimension | `not_required` | Not a vacuous pass; individual results remain visible. |

An optional **biological failure** remains visible without overriding passing
required references. A partial/cancelled recorded execution cannot create an
aggregate pass, even if another child finished successfully. Genomic and
transcriptome dimensions are reported separately, not combined into a tier score.

## 5. Keep History And Panel Readiness Distinct

Standalone summaries are immutable and content-derived. Reimporting the same
scientific evidence returns the same ID; changed evidence produces another
historical result. Show/list do not probe databases. Their `current_at_import`
label describes the time of import, **not current applicability today**. Reimport
against changed resources checks currency again and cannot turn stale output into
a pass. A changed saved pair/template requires fresh preparation.

Existing panel finalization uses a separate exact selection receipt for each
genomic/transcriptome dimension. Complete, validated re-finalization can deliberately
replace that dimension's active selection; the old receipt and reports remain
inspectable in persisted evidence. Incomplete finalization cannot replace it.
Legacy reports without binding remain readable but cannot earn a passing readiness
gate. This is why timestamp ordering is insufficient.

**The multi-reference import does not attach its summary to panel readiness.**
Keep using the existing approved panel workflow for those gates. Neither a standalone
pass nor a panel's generated primers establishes isoform discrimination, oligo QC,
laboratory validation or order readiness. Those need their own evidence and review.

## 6. Use The Same Contract Through Agents And Scripts

MCP exposes `primer_specificity_multi_handoff`, `primer_specificity_multi_import`,
`primer_specificity_multi_show` and `primer_specificity_multi_list`. Preparation
and import require `confirm: true`. Raw `op` calls, JS/Lua/Python and workflows use
the same four engine operations documented in [the protocol](../protocol.md).
There is no adapter-local biological rule or automatic external execution.

Suggested inner-agent message:

> Inspect my saved pair and propose a multi-reference request for these exact
> catalog IDs. Explain which references are required, any missing intended-target
> mappings, both specificity dimensions and the effective policy. Do not download,
> prepare resources, start BLAST, write a handoff or import evidence until approved.

Have Glen test the shared deterministic routes independently of language quality.
This page is a source/readability review, not an inner-agent or native-GUI test.

## Review Questions

Can you identify the exact reference behind every passing row? Can you explain
why an unavailable optional reference is not "zero off-targets"? Does each cDNA
subject mapping really belong to that reference? Which separate experiments would
still be needed to claim isoform discrimination and laboratory performance?

For the preceding design/coverage workflow, return to
[04.08: primer-pair markers for PATZ1 isoforms](04-08_gene_assay_study_gui.md).
For interface details, see [CLI](../cli.md) and [protocol](../protocol.md).

## Verification Scope

Authored synthetic Rust regressions cover admission, retained LF/CRLF bytes, mappings,
receipt tampering, independent verdicts, incomplete execution, immutable history,
save/reload and stale references. Execution remains with Glen/CI under the owner's
no-local-build instruction. No real-data specificity, native GUI acceptance or
order approval is claimed. The examples above are learning templates, not fixtures
with an assumed passing result.

Glen/CI should run these on one frozen revision, keeping results distinct from
real-reference and native-GUI acceptance:

```sh
cargo test --locked --lib specificity_reference_
cargo test --locked --lib specificity_multi_
cargo test --locked --lib writing_routes_require_confirmation_before_state_or_file_access
cargo test --locked --lib primer_specificity_handoff_plans_without_running_and_imports_completed_outputs
cargo test --locked --lib transcript_assay_panel_specificity_finalization_is_atomic_and_distinguishes_outcomes
cargo test --locked --test capability_registry_parity
cargo test --locked --test parity_matrix_freshness
cargo run --locked --bin gentle_examples_docs -- tutorial-catalog-check
cargo run --locked --bin gentle_examples_docs -- tutorial-check
cargo check -q --locked
python3 -m unittest scripts.test_primer_specificity_multi_tutorial scripts.test_tutorial_checkouts
```

The Unix fake-tool import fixture is deterministic process-contract evidence,
not a real BLAST search. Retained-byte/path tests are portable; native Windows
filesystem behavior still needs Windows CI. This chapter is manual/hybrid, so
the offline executable tutorial replay does not certify its external-runner step.
