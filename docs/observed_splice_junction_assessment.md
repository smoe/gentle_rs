# Observed Splice-Junction Assessment (J1)

Assess **observed intron intervals** against the exact references used by a
saved [S1 report](observed_transcript_structures.md). This experimental, offline
service adds offsets, optional motifs and declared sample grouping. It never
realigns reads, infers expression, assigns proteins or establishes complete RNA.

## Commands And Bindings

```text
rna-reads splice-junctions-assess REQUEST_JSON_OR_@FILE
rna-reads splice-junctions-show REPORT [--table TABLE.tsv] [--rows N] [--binding JSON_OR_@FILE] [--verify-artifacts]
rna-reads splice-junctions-cancel REQUEST_JSON_OR_@FILE
```

Shared-shell/CLI routes map to `AssessObservedSpliceJunctions`,
`InspectObservedSpliceJunctions` and `CancelObservedSpliceJunctions`.
Generic MCP `op` still requires `confirm:true`, including inspection.
Existing JS/Lua generic operations and Python's CLI wrapper share the engine.
Assessment/cancel write artifacts, not project content. Managed execution
discards snapshot audit entries; raw engine application retains its normal audit
entry without adding undo/redo content checkpoints. Cancel does not confirm
worker termination.

Schemas: `gentle.observed_splice_junction_assessment_request.v1`,
`gentle.observed_splice_junction_assessment_report.v1` and
`gentle.observed_splice_junction_assessment_result_binding.v1`.
Requests reject unknown fields and require explicit bindings and budgets.

| Field | Contract |
| --- | --- |
| `source` | `{"kind":"s1_report","result":S1_RESULT_BINDING,"report":{"path":PATH,"sha256":RAW_SHA256,"bytes":RAW_BYTES}}`. Verify raw/semantic hashes and every S1 artifact to EOF. Non-completed S1 gives unavailable, counts absent. |
| `reference_sequence` | Optional `{file,format:"fasta_plain_v1",authority}`. Authority is `{"kind":"full_assembly"}` requiring S1's exact FASTA SHA-256, or `{"kind":"declared_excerpt","contig":NAME,"offset_0based":N,"length":N}` with only a declared assembly relationship. |
| `sample_manifest` | Optional raw binding of `gentle.observed_sample_manifest.v1`. Exact context mapping described below. |
| `near_window_bp` | Required 1..50, no default. |
| `recurrence_policy` | Optional `{unit,min_units,certainty,require_qualified_flanks}`. Unit: `library`, `cdna_prep`, `biological_sample`; min_units: 2..1,000,000; certainty: `definite`, `definite_or_possible`; qualified flanks must be required. No policy means no label. |
| `budgets` | Eight positive limits below; exhaustion fails, never truncates. |
| `output_dir`, `cancel_path` | Fresh output and optional exact-request marker outside it. No overwrite. |

Original hash-bound normalized annotation paths supply the full model geometry
that S1's summary does not duplicate. Missing catalogs yield source-specific
`not_assessable`; changed bytes or contradictory bindings are refused. There is
no annotation replacement argument or download.

### Resource Ceilings

| Required budget | Admitted maximum |
| --- | ---: |
| `max_input_bytes` | 1,000,000,000,000 |
| `max_reference_bytes` | 64,000,000,000 |
| `max_manifest_bytes` | 8,000,000 |
| `max_distinct_junctions` | 100,000 |
| `max_annotated_junctions_compared` | 1,000,000 |
| `max_near_hits_per_junction` | 100,000 |
| `max_total_work` | 1,000,000,000,000 |
| `max_output_bytes` | 128,000,000 |

Input charges include repeated verification reads. FASTA streams to EOF;
structured source JSON is capped at 128 MB/file. FASTA lines: 1 MiB; sequence
names: 4096 bytes, 100,000 headers and 8 MB total names. Manifests: 100,000 rows.
Selected inspection: 1000 rows and 8 MB; each emitted table row: 8 MB. Output
budget includes the final report and all artifacts. These engineering ceilings
are not evidence that every assembly fits a caller's chosen byte/work budget.

## Interpretation

Physical identity binds assembly/FASTA hash, contig and exact zero-based
half-open intron `[start,end)`. Full SHA-256 IDs are authoritative; 8-hex display
IDs can collide. Offsets are **annotated minus observed** at each end. Near
never becomes exact; coordinates are never snapped. Sources retain all nearby
pairs and nonexclusive geometry events. The `novel_combination` event is a
literal structural label, **not a biological novelty claim**.

Model strands, never query orientation or motifs, determine transcription
orientation. Disagreement/unresolved sources leave it unassigned. Raw
reference-forward first/last dinucleotides and forward/reverse motif classes
are retained. Classification handles soft masking and reverse complements;
ambiguity/excerpt-edge gaps give `not_assessable`. Verified full-assembly and
declared excerpt authority stay distinct.

Each `(source_sha256,record_index,read_id)` contributes at most once per junction
at highest supported certainty. `definite` and `possible_only` are disjoint.
Each partitions into qualified_only, insufficient_only, unassessed_only and
mixed; mixed evidence is not recurrence-eligible. Incomplete search never gives
definite support. Separate reads never manufacture a junction. Saved model
relations, end constraints and clipping accompany `rna_end_status:"not_verified"`.

The manifest has `schema:"gentle.observed_sample_manifest.v1"` and `rows[]`
containing `context:{sample_id,barcode,priming,file_id}`, `library_id`,
`cdna_prep_id`, `biological_sample_id` and optional `cell_line`. Match contexts
exactly. One library maps to one prep and one prep to one biological sample;
contradictions are refused. Unmatched contexts and unused rows stay visible.
Repeated libraries do not establish independent biological samples. One record
assigned to multiple claimed units cannot establish recurrence. Grouping axes
overlap and must not be summed.

Library denominators distinguish total cohort, gene-coherent, excluded and
unassessed records and context coverage. Zero gene-coherent records is valid
non-detection, not absence. Policy results are `recurrent_candidate`,
`not_recurrent_under_policy` or `not_assessable`; none establishes novelty,
translation, secretion, molecules or expression.

## Saved Output

Publish `report.json` last. Hashed `artifacts/` contain `assessment.json`
(`gentle.observed_splice_junction_assessment_data.v1`), `junction_assessment.tsv`,
`chain_assessment.tsv`, `sample_support.tsv` and `unavailable.tsv`.
Tables use one quoted `row_json` column. Sample rows distinguish denominator,
junction support and unused manifest rows. All derive from one typed result.
Fault cleanup removes only owned identity-matched files; unowned/orphan output
is disclosed rather than recursively deleted.

Inspection defaults to metadata/zero rows. Bounded rows and `--verify-artifacts`
are explicit; inspection never reopens S1, annotation, FASTA or manifest inputs.
See the [plan and execution ledger](observed_junction_consequences_plan.md).
The synthetic CD274-like fixture is not real RNA evidence. Historical legacy
`5466827 -> 5468193` used mixed display; canonical storage is
`[5466826,5468192)`, exon-boundary label `5466826 -> 5468193`.
CDS/ORF, UniProt and soluble-form hypotheses remain owner-gated J2/J3 work.

The Phase 2/S1 compatibility fix admits explicitly complete searches with
`declared_complete=true` and an exactly empty `competitor_search_limitations`.
S1's unchanged certainty rules can then produce definite support. Nonempty
limitations, including `none`, remain conservative; incomplete declarations do
not become definite. Whitespace-only descriptions are rejected rather than
normalized into an empty declaration. Existing reports are not rewritten or
promoted. The synthetic fixture exercises the saved Phase 2 -> S1 -> J1 path;
it still does not establish real search exhaustiveness or RNA biology.
