# Bound Coding Transcript Models (C1)

C1 derives inspectable coding models from explicitly supplied, checksum-bound
annotation/reference files. It does not reinterpret cDNA allele calls or produce
HGVS/consequences. Existing comparison and normalized annotation v1 bytes remain
unchanged.

```sh
gentle_cli rna-reads coding-models-derive @coding.request.json
gentle_cli shell 'rna-reads coding-models-derive @coding.request.json'
```

Both routes return `result.coding_transcript_models`. GUI Shell, generic MCP
`op`, JS/Lua shell bindings and Python CLI delegation share the engine service.
No project, output path or provider contact is required. These expensive
read-only observations consume execution capacity, not a project baseline; they
do not save, publish files or create undo checkpoints. No new GUI window is added.
Generic MCP `op` still requires its existing `confirm: true` before input I/O;
C1 does not bypass that adapter gate.

## Request And Binding

`gentle.coding_transcript_models_request.v1` rejects unknown fields:

| Field | Required content |
| --- | --- |
| `gtf`, `normalized_annotation` | Local raw-byte `{path,sha256,bytes}` bindings |
| `reference` | Comparison-compatible `assembly`, `contig`, `start_0based`, declared `genomic_reference_sha256`, bound single-record `locus_fasta` |
| `selector` | Exact `source_id` and `gene_id`, no symbol guessing |
| `provider` | Declared `name`, `release`, explicit `cds_row_convention` |
| `reference_kind`, `translation_table` | Explicit `nuclear` and `1`; other contexts refused in v1 |
| `protein_fasta` | Optional bound protein FASTA with exact versioned protein IDs |
| `budgets` | Optional positive limits with hard ceilings below |

The convention is `stop_codon_excluded_from_cds_rows` or
`stop_codon_included_in_cds_rows`. Provider names never select it. No contig-name
heuristic licenses a code: the producer must declare reference context honestly.
Organelle/non-nuclear interpretation is outside v1. Discovery templates are
typed placeholders, not executable bindings.

Supplied bytes are reverified under existing local-file/symlink/path rules.
Annotation assembly, release, reference identity and source must agree. Selected
GTF transcript IDs, strands and exact ascending genomic exons must equal the
normalized annotation. Unselected GTF rows consume parsing/work limits but are
not checked against the selected locus's reference length. Hashes bind the
supplied snapshots, not proof that the annotation was generated from that GTF.
The whole-genome digest is declared/cross-bound, not reauthenticated by rereading
a whole genome.

## Coordinates And Biological Boundaries

Raw CDS rows retain GTF frame. Projection reuses the existing read-end exon map
and supported-block geometry. Minus strands become mature-transcript 5'-to-3'
sequence. Phase is checked across segments; bases are **not** trimmed anew at
every segment. Start/stop rows may span exon junctions. Known stop geometry must
agree with the explicit convention. Missing stops never gain three guessed bases.

`coding` retains raw CDS/start/stop/UTR rows, normalized genomic segments and
the mature-cDNA half-open coding interval including an explicitly known stop,
plus 1-based inclusive CDS coordinates, codon strings and protein binding.
Missing stop evidence leaves the supplied boundary unchanged and translation
incomplete. Transcript `mRNA_start_NF`/`mRNA_end_NF` tags stay separate and alone
do not invalidate a complete CDS. Absence of such tags does not verify a complete
RNA. CDS-specific tags, missing/partial start/stop or initial nonzero phase
withhold verified translation. Absent CDS remains `noncoding` or
`coding_not_defined`.

Binding/structural errors refuse the request: exon disagreement, overlapping or
out-of-exon segments, noncontiguous mature CDS or contradictory stop convention.
Biological irregularities are per-model `translation_reasons`: incomplete
start/end, inconsistent phase, non-ATG start, non-triplet length, ambiguous bases,
internal stop or non-stop terminal sequence. Nothing is repaired.
`translation_status` is the first reason, or `complete_verified`; inspect the
full reasons array, not only the first label.

Translation uses `AMINO_ACIDS.codon2aa` table 1. `genetic_code_sha256` binds the
actual `assets/codon_tables.json`, not the separate codon-usage catalog.
`translated_protein_sha256` hashes uppercase amino-acid letters without a
terminal stop or newline; absent when sequence/completeness checks fail. Optional
protein FASTA permits one terminal `*` and reports `not_supplied`, `identical`,
`differs`, `missing` or `not_assessable` (translation withheld). A supplied
mismatch is flagged, not reconciled; its sequence-derived digest is retained.
This is a sequence prediction, not expression evidence.

Version states are `absent`, `declared` or `inconsistent`. Accessions are emitted
only for consistent numeric evidence (including an explicit numeric ID suffix).
Conflicting suffix/attributes or mixed missing/present version attributes
withhold the accession. Unversioned geometry can translate but does not license
future HGVS notation.

## Determinism And Limits

Report schema: `gentle.coding_transcript_models.v1`. `request_sha256` hashes the
typed request. `content_sha256` hashes compact `serde_json::Value` text with that
field omitted. `sequence_sha256` hashes uppercase concatenated locus sequence,
excluding headers/newlines, as in cDNA comparison. Raw LF/CRLF hashes remain
different while model content agrees.

| Budget | Default | Hard ceiling |
| --- | ---: | ---: |
| All input bytes | 64,000,000 | 128,000,000 |
| GTF/protein line bytes | 65,536 | 1,048,576 |
| Rows per text input | 500,000 | 1,000,000 |
| Selected transcripts | 256 | 4,096 |
| Selected GTF rows / exons per model | 100,000 | 1,000,000 |
| Mature transcript/protein bases | 2,000,000 | 16,000,000 |
| Work units | 100,000,000 | 100,000,000 |
| Serialized output bytes | 8,000,000 | 64,000,000 |

The shared locus reader also caps FASTA bytes at 16,000,000. Work covers parsing,
projection and translation; it is not a measured RSS/throughput guarantee.
Cooperative cancellation refuses rather than publishing a partial report.

Inline Rust fixtures and `scripts/test_coding_transcript_models.py` invent
two-/single-exon loci, MKP proteins and stop rows; no biological sample is copied.
The fast checkout gate recreates raw bindings in LF/CRLF checkouts. Actual CLI
tests require `GENTLE_CODING_MODELS_CLI`; skips are not runtime acceptance.
See the [C1 Glen checklist](glen_next_verification_handoff.md#c1-bound-coding-models-2026-11-11).

HGVS/consequences, codon combinations, same-unit isoform linkage and frozen known
labels belong to the separate downstream [C2-C4 contract](cdna_allele_interpretation.md),
not to this C1 report. C1 establishes no genotype, zygosity, expression, CNV,
somatic/germline, pathogenicity or NMD verdict. No private reads or completed
cDNA analyses are rerun.
