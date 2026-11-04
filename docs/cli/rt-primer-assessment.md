# Prospective RT-Primer Assessment

<!-- manual-counterparts:start -->
No direct GUI counterpart: Shared-shell RT background assessment has no dedicated native GUI form in this slice.
<!-- manual-counterparts:end -->

`primers assess-rt-candidates REQUEST_JSON_OR_@FILE` screens explicitly supplied
RT-primer tails against explicitly supplied FASTA backgrounds. Use it through
the GUI Shell, `gentle_cli shell`, or generic MCP/JS/Lua/Python operation access:
`{"AssessRtPrimerCandidates":{"request":...}}`.

This is an unweighted compatibility report, not a probability of RT initiation,
full-tail specificity verdict, primer design or ordering recommendation. A
unique full-tail genomic match does not imply a rare terminal seed. No private
reads, aligner, reference download or automatic discovery is involved.

## Request

This synthetic request uses invented sequences, not a TP73 assay. Replace the
file path, byte size and raw-byte SHA-256 with a real completed uncompressed
FASTA binding. The placeholder deliberately fails validation.

```json
{
  "schema": "gentle.rt_primer_assessment_request.v1",
  "primers": [{
    "primer_id": "invented",
    "assay_label": "RT compatibility example",
    "sequence": "TTTTTTGGGGACGTACGA",
    "components": [
      {"component_id":"handle","role":"common_handle","start":0,"end":6},
      {"component_id":"tail","role":"target_core","start":6,"end":18}
    ],
    "core_orientation": "antisense_3prime",
    "expected_targets": [],
    "variant": "synthetic-v1",
    "source": "Hand-crafted example; not a physical oligo",
    "source_sha256": null
  }],
  "references": [{
    "reference_id": "synthetic_rna",
    "scope": "custom_rna",
    "source": "Explicit synthetic reference revision",
    "fasta": {"path":"/REPLACE/background.fa","sha256":"REPLACE_RAW_SHA256","bytes":0},
    "unavailable_reason": null,
    "circular_record_ids": []
  }],
  "intended_sites": [],
  "budgets": {
    "max_input_bytes": 16000000,
    "max_reference_bases": 1000000,
    "max_record_bases": 1000000,
    "max_records": 10000,
    "max_base_comparisons": 100000000,
    "max_examples_per_cell": 8
  }
}
```

Primers reuse the trace-panel oligo contract, without changing its serialized
bytes or digest semantics. Components must tile the complete canonical DNA
oligo, with exactly one terminal `target_core`, length at least 10, and at most
one preceding `common_handle`. A single core component explicitly supplies a
core-only primer. Complete oligos remain limited to 64 bases. Primer IDs must
be unique; different handles may share the same tail. Orientation must be
`antisense_3prime`, independent of FASTA genomic strand. Provenance, assay
label and variant are retained; existing `expected_targets` are citations and
do not automatically exclude any search sites.

Scopes are `genome`, `transcriptome`, `mitochondrial_rna`, `ribosomal_rna` and
`custom_rna`. Genome scans both template strands. All RNA scopes expect
biological 5′-to-3′ sequences using the DNA alphabet (`T`, not `U`) and scan
only the complementary primer orientation. A mitochondrial DNA reference
therefore uses `genome`, with its origin topology explicitly declared.

Every reference has a unique ID and nonempty source/revision description.
Supply exactly one of a `fasta` binding or an `unavailable_reason`. A named
missing scope produces `not_assessed`, with no clean-negative implication.
Unlisted scopes have not been searched. Labels are not reference validation:
the binding hashes the exact bytes consumed to EOF and checks file stability.
LF and CRLF FASTA are accepted with their own raw-byte hashes. IDs use the
first FASTA header token and must be unique; sequence lines accept canonical
or ambiguous IUPAC DNA, without spaces. Empty records, blank sequence lines,
invalid characters, missing circular IDs and checksum changes fail explicitly.

## Read The Matrix

The report includes its complete request, request/content SHA-256, verified
reference summaries and one cell per primer/reference/terminal length 6/8/10/12.
Counts are exhaustive within those FASTA bytes when a report returns; filling
the example list never truncates counts. It does not certify that the supplied
FASTA contains an entire genome or transcriptome.

For each physical sequence window, `mismatch_histogram[k]` counts the minimum
number of positional mismatches across searched template orientations.
Thus ≥6/8 positional matches are bins 0, 1 and 2 of the eight-base cell. The
`terminal_run_histogram[k]` independently counts the longest exact contiguous
run starting at the primer's 3′ end. A terminal mismatch can give a zero run
despite ≥6/8 positional matches. These are marginal histograms, not a joint
distribution: different genomic orientations may supply the two best values.
There is no indel or thermodynamic binding model.

Both genomic orientations at one position count as one physical window;
palindromic exact examples retain both `template_strands`. Repeated sites in
one record count as multiple positions but one `exact_matching_records`.
Scopes may duplicate records or overlap biologically. Results are never pooled
across references or interpreted as molecule counts.

`assessed_windows` is the histogram denominator. Ambiguous windows are counted
separately in `ambiguous_windows`; short records in
`records_shorter_than_seed`. A short tail or no unambiguous windows yields
`not_assessed` and a reason. `completed` means the declared scan finished;
it does not remove these explicit unknown windows. Zero exact matches is
meaningful only alongside these denominators and the declared input scope.

Linear records have `length - seed_length + 1` windows when long enough.
Records named in `circular_record_ids` have one window per start position,
including origin crossings. Seeds longer than the record are unassessed;
multi-lap circular binding is not modelled. Example end coordinates may exceed
the circular record length. All coordinates are 0-based, half-open.

Optional `intended_sites` supply `primer_id`, `reference_id`, `record_id`,
`tail_start_0based`, `tail_end_0based_exclusive`, and `template_strand` (`+`/`-`).
They must match the entire tail exactly in the declared reference and template
orientation; unknown, mismatching and out-of-range sites fail. In this slice
intended full-tail intervals must be linear, non-wrapping intervals even on a
circular record. RNA sites permit only `+`. A `+` template maps the primer's
terminal seed to the beginning of this interval; a `-` template maps it to
the end. `excluding_intended_orientations` removes that orientation only;
the opposite genomic orientation remains eligible. Always compare its own
assessed denominator with `all_sites`. No gene-name or approximate exclusion
is inferred. Exact examples describe `all_sites`, not the excluded subset.

## Work Limits And Effects

Requests admit 1–16 primers and 1–16 references, with up to 1,024 intended
sites. Defaults are shown above. Hard ceilings are 256 MB total declared input,
16 million total reference bases, 16 million bases per record, 100,000 total
records, one billion base comparisons, and 32 retained exact examples per cell.
FASTA lines are capped at 1 MiB including their newline. Byte/base/record/work
budgets are cumulative across all references, including overlapping inputs.

A single candidate over a million-base both-strand background requires at
most `2 × (6+8+10+12) × 1,000,000 = 72,000,000` base comparisons. This motivates
the default 100-million work cap; multiplying candidates/references multiplies
work, not the cap. These are admission bounds, not measured latency promises
or whole-human-genome support. Exceeded limits return an error with no partial
report. Ambiguous windows still consume the comparison budget. Input parsing
retains at most one bounded record; the FASTA is read/hash-verified once per
reference for all supplied candidates.

The operation cooperatively polls the existing progress callback, including
before work and before returning the final report. Cancellation or a source
failure returns an error, with no partial report, file write or readiness
change. Shared-shell assessment leaves project content/history unchanged.
Generic engine apply records its ordinary operation journal, without an undo
checkpoint or content revision. Managed command execution uses the existing
expensive observation lane and needs no project snapshot. MCP's generic-op
confirmation policy is unchanged.

The returned portable JSON is the assessment artifact; this slice does not
add a saved-report store or a file exporter. `abundance_status` and
`retrospective_status` are `not_assessed`. The separate comparison operation
below adds per-isoform context and descriptive oligo metrics. Hairpin prediction,
automatic redesign and native GUI forms remain later work. No assessment changes
existing primer-order readiness.

## Compare Candidates Across Annotated Isoforms

`primers compare-rt-candidates REQUEST_JSON_OR_@FILE` combines the prospective
assessment with exact RT-tail binding searches against explicitly selected
loaded transcript sources. Its generic operation is
`{"CompareRtPrimerCandidates":{"request":...}}`; the portable
`gentle.rt_primer_comparison.v1` report is returned in
`OpResult.rt_primer_comparison`. This is a separate contract: the assessment
request/report v1 shapes and their counting semantics remain unchanged.

The comparison request has schema `gentle.rt_primer_comparison_request.v1`,
embeds an ordinary assessment **request** in `assessment`, and adds `sources`
and comparison budgets. GENtle runs that assessment once through the existing
shared operation. It does not trust supplied saved counts or silently reuse
an assessment from another candidate/reference revision. A previously saved
assessment request is suitable input; a saved assessment report is not.

For example, adapt a fully prepared `assessment-request.json` from the earlier
example after loading an annotated DNA record as `locus`:

```python
import json
from pathlib import Path

request = {
    "schema": "gentle.rt_primer_comparison_request.v1",
    "assessment": json.loads(Path("assessment-request.json").read_text()),
    "sources": [{
        "source": {
            "seq_id": "locus",
            "source_feature_id": 0,
            "coverage_universe": {"kind": "all_annotated_cdna_classes"},
            "annotation_release": "synthetic-v1"
        },
        "expected_source_record_sha256": None
    }],
    "budgets": {
        "max_transcript_base_comparisons": 100000000,
        "max_bindings": 10000
    }
}
Path("comparison-request.json").write_text(json.dumps(request, indent=2))
```

```text
primers compare-rt-candidates @comparison-request.json
```

The Python only prepares request metadata. Feature `0`, record `locus` and
`synthetic-v1` are explicit example bindings, not claims about a real gene or
loaded annotation. GENtle interprets the transcript geometry and biology.
The embedded assessment still needs valid explicit FASTA bindings, or named
missing-reference scopes whose results remain `not_assessed`.

### Isoform Presence And Region Context

Each source uses the existing `TranscriptCaptureSource` and shared coverage
resolver. The default universe retains all selected annotated transcript
members. An `explicit_transcripts` universe supplies
`required_transcript_ids`; unresolved or ambiguous required IDs fail rather
than disappearing from the denominator. The retained source bindings include
coverage resolution, available genome anchors and the full loaded DNA-record
hash, including annotation. `expected_source_record_sha256`, when supplied,
refuses a changed source record. It is not merely a sequence-base checksum.
Each member also retains its biological strand, mature cDNA hash and length,
transcript identity and transcript-local annotated CDS ranges. UniProt inventory
joins are not added here; resolve that inventory to an explicit supported
transcript selection first.

RT tails bind complementary RNA. The search therefore uses the reverse
complement of the supplied terminal target tail against biological 5′-to-3′
mature transcript sequence. The handle is not part of this binding search.
Source genomic strand does not reverse the RT assay meaning. Every exact
occurrence is retained, including exon-junction-spanning sites; source ranges
remain 0-based, half-open, ascending within each exon and in biological exon
order on either genomic strand. No single representative isoform is chosen.

The report has an `isoforms` row for every candidate/member pair. `presence`
and each binding's `region` answer different questions:

| Field | Value | Meaning within this supplied transcript scope |
| --- | --- | --- |
| `presence` | `present` | At least one exact, canonical full-tail binding site exists. |
| `presence` | `absent` | No exact site or ambiguity-compatible possible site exists. This does not establish absence outside the selected transcripts. |
| `presence` | `unresolved` | No proven exact site, but one or more ambiguous windows could match the tail. |
| binding `region` | `cds` | The complete site lies in the mapped annotated CDS union. |
| binding `region` | `five_prime_utr` / `three_prime_utr` | The complete site lies before / after the annotated CDS in biological transcript coordinates. |
| binding `region` | `mixed` | The site spans more than one known region class. |
| binding `region` | `unknown` | Any site bases lack a resolvable CDS/UTR classification. |

`possible_ambiguous_sites` counts ambiguity-compatible windows separately
from exact sites. A `present` row may still have additional possible sites.
An ambiguous symbol is not exact evidence, and a canonical incompatibility
elsewhere in that window rules it out. A transcript shorter than the full
tail has no such full-tail site; its member is still retained.

An empty CDS result means unknown coding annotation, including transcripts
without a CDS model. It never makes the whole transcript UTR. CDS intervals
are treated as a union: gaps inside their envelope remain unknown. A binding
with any unknown bases has region `unknown`, even if other bases are known
CDS or UTR; its individual base counts remain visible. `mixed` is reserved
for crossings among verified known classes. Exact CDS/UTR boundary endpoints
use half-open coordinates. Unresolvable source geometry and conflicting CDS associations fail explicitly;
an absent CDS model remains unknown. The operation does not borrow a gene-wide
coding model.

Comparison-specific admission also refuses an empty mature transcript,
multiple `cds_ranges_1based` qualifiers on one selected transcript, more than
4,096 nonempty comma-separated explicit CDS tokens, or a qualifier mixing
valid and malformed tokens. This prevents silently dropping part of a CDS
model before classification. An entirely malformed or empty explicit CDS
qualifier remains unknown under the existing resolver. These guards do not
change the shared parser's legacy policy or silently repair source annotation.

### Candidate Comparison Rows

`candidates` retain input order without an automatic rank or winner. For each
primer, `total_members`, `members_with_exact_site`,
`members_with_possible_only_site` and `members_without_site` preserve the
selected universe denominator. Multiple exact occurrences in one transcript
count as one member with a site. `region_member_counts` counts a member once
per binding-region class; a member with separate CDS and UTR sites contributes
to both classes. These overlapping counts must not be summed as distinct
isoforms. The detailed binding rows explain every classification.

The embedded `assessment` supplies the terminal-risk matrix, joined by
`primer_id`. Its reference IDs, candidate definitions, hashes and explicit
unknowns remain visible alongside isoform coverage. Different background
scopes are not pooled into one specificity score. A low burden in the
old-primer-selected cohort is not substituted for this prospective evidence.

Both `target_tail` and `full_oligo` contain descriptive length, GC fraction,
Tm, longest homopolymer and shared self/general 3′-complementary-run heuristics.
The report's `tm_model` states the shared estimator's assumptions; these are
not measured reaction temperatures. Full-oligo Tm is a sequence diagnostic,
not the RNA-binding Tm of a noncomplementary handle plus tail. Complementary
runs are not hairpin energies or calibrated structure probabilities.
`secondary_structure_status` is `not_assessed`: no new Primer3 invocation or
hairpin prediction is implied. Concentration, abundance weighting, isoform
expression, new candidate generation and order readiness remain unassessed.

### Comparison Limits And Effects

The comparison accepts 1–16 loaded sources, at most 256 materialized
transcripts and one million total mature cDNA bases. Source admission is
bounded to 16 million loaded DNA bases and 100,000 features, cumulatively
across source bindings, including repeated underlying records. Before source
hashing or Splicing Expert view construction, annotation preflight admits at
most 16,384 location nodes plus nonempty exon/CDS qualifier tokens, 16 MB of
feature-qualifier and external-location text, and location depth at most 32.
It also bounds cumulative conservative geometry units to 25 million: each
source contributes `U × U × max(S, 1)`, where `U` is the number of unique
source/exon/CDS interval pairs from splicing-seed features and `S` is their
feature count. These guards cover the shared view's work before cDNA
materialization, not just the later binding scan. They conservatively inspect
the loaded source annotation; a narrow selected transcript group does not
remove other source annotation from admission.

The transcript search defaults to 100 million base comparisons and 10,000
exact binding records, with hard ceilings of one billion comparisons and
50,000 bindings.
These budgets are separate from the embedded prospective assessment's FASTA
budgets; neither expands the other. Excess work or binding output fails
explicitly without a partial comparison report.

The existing callback supports cancellation before work, during scanning and
before the final return. The comparison needs a snapshot of the explicitly
loaded transcript sources for managed execution. It does not inherit the
file-only assessment's project-free observation classification or introduce a
new artifact-publication exception. Shared-shell execution leaves project
content/history unchanged. Generic engine apply keeps its ordinary journal
without an undo checkpoint or content revision. No file is written and no
primer-order readiness is upgraded. Native GUI forms, real-reference
performance and biological acceptance require separate evidence.

## Inspect Saved Primer–Payload Junction Evidence

`rna-reads primer-junctions REQUEST_JSON_OR_@FILE` joins a saved primer trace,
its prepared payload queries and saved Phase 2 alignment evidence. Use it in
the GUI Shell or CLI, or through the shared
`{"InspectPrimerPayloadJunctions":{"request":...}}` operation. The portable
`gentle.primer_junction_inspection.v1` result is returned in
`OpResult.primer_payload_junctions`.

This answers which observed primer interval led to each saved payload query,
which partner alternatives were retained, and how many distinct source reads
support each saved profile. It does not rerun trimming, an aligner or an RT
experiment. It never reopens original FASTQ, original PAF, reference FASTA or
annotation files. Previously verified producer declarations remain historical
evidence, not a new exhaustive-search or biological validation.

### Inspection Request And Saved Inputs

Use the exact alignment result binding returned for the saved Phase 2 report,
plus the report file's raw-byte SHA-256 and size. The example below contains
intentional placeholders and cannot run until they are replaced.

```json
{
  "schema": "gentle.primer_junction_inspection_request.v1",
  "source": {
    "schema": "gentle.primer_payload_result_binding.v1",
    "kind": "alignments",
    "request_sha256": "REPLACE_SAVED_REQUEST_SHA256",
    "report_content_sha256": "REPLACE_SAVED_REPORT_CONTENT_SHA256",
    "report_path": "/REPLACE/payload-run/report.json",
    "outcome": "completed"
  },
  "source_report": {
    "path": "/REPLACE/payload-run/report.json",
    "sha256": "REPLACE_RAW_REPORT_FILE_SHA256",
    "bytes": 0
  },
  "file_id": null,
  "primer_id": null,
  "offset": 0,
  "rows": 100
}
```

The semantic report-content digest and raw file digest serve different checks;
do not substitute one for the other. The inspector verifies the saved report
binding and the saved trace metadata, query manifest/FASTA and trace-row
artifacts needed to establish anchor identity, query geometry and masks. A
historical `retained_primary_fasta` query preparation also requires its saved
retained payload FASTA. A missing retained dependency fails explicitly; the
inspector does not fall back to original-read replay. A historical
`replay_verified_fastq` source mode does not authorize a new FASTQ read here.

All three Phase 2 sidecars are hash-verified to EOF before returning even a
small preview. The inspector binds evidence to its saved query, batch, read
record, primer alternatives and sample context, requires the complete saved
population for completed batches, and recomputes the existing partner-profile
and breakpoint-cluster aggregates. Inconsistent saved arithmetic, missing or
duplicated evidence, and hash failures refuse the inspection. These checks
preserve the saved evidence contract; hashes are not an authenticated producer
signature or proof that an alignment search was scientifically complete.

### Read The Joined Rows

`primer_panel` retains the complete ordered oligos and component definitions.
Each returned row has `file_id`, source-file SHA-256 and sample/barcode/priming
context, plus its original `trace`, prepared `queries`, saved payload
`evidence` and `unassessed` dispositions. `payload_previews` maps each joined
query ID to its first at most 64 saved query bases. These bases are already
masked and oriented into the observed primer direction; they are not newly
trimmed read sequence or reconstructed reference flanks.

Trace hit coordinates remain 0-based, half-open coordinates of the original
read. For a forward anchor, the observed payload is after `anchor.end`; for a
reverse-complement anchor it is before `anchor.start`. Prepared query evidence
is already oriented into the observed primer direction. Do not reverse it a
second time or equate this orientation with a transcript's biological strand.
The saved masks use query coordinates; trace intervals use original-read
coordinates. Both coordinate systems remain visible.

An approximate full-primer hit uses its observed aligned interval, including
indel effects. A core-only hit does not establish that the complete handle was
observed or removed. Query anchors are representatives of a shared geometry;
other primers at that geometry may have different full/core classes or edit
counts. Inspect all retained `trace.hits`, rather than assigning the
representative hit's class to every `primer_id`. Handle-only observations
remain unattributed and do not become gene-specific primer evidence.

Ambiguous primer sets, partner-gene sets, equivalent segments and primer-edge
alternatives remain sets. The inspector does not choose a new winning primer,
partner or breakpoint. Terminal seed mismatch counts and unavailable reasons
are displayed as saved. Phase 2 does not retain literal reference-flank bases,
so no local RNA sequence or new base-by-base flank alignment is invented.

### Denominators, Selection And Outcomes

Distinct read identity is `(file_id, record_index)`, not the display `read_id`.
Two records with the same read name remain two source records. Multiple anchors
on one record can contribute multiple `anchor_observations` while counting once
in that group's distinct-read support. Duplicate source-file SHA-256 values in
the selected verified cohort are refused so file aliases cannot inflate counts.
These are record counts, not counts of independent original RNA molecules.

`cohort_distinct_reads` includes every trace record in the selected, verified
query-source units, including handle-only and excluded/unassessed records.
Each profile keeps its own reference, sample, barcode, priming, complete primer
ambiguity set and complete partner-equivalence set. The eligible, assessed,
ambiguous, unmapped, unassessed and excluded denominators retain the existing
Phase 2 meanings. A read can contribute to overlapping groups or buckets; do
not sum these as mutually exclusive reads, genes or molecules. Reference scopes
are not pooled. Internal-boundary and primer-adjacent edge clusters remain
separate, with the saved ambiguity and cluster-policy restrictions.

Optional `file_id`/`primer_id`, `offset` and `rows` select joined trace previews
only. They never filter the verified cohort used for profiles, clusters or
denominators. `matched_trace_rows` reports the matching row population before
pagination, and `rows_truncated` discloses a capped preview. Row selection is
deterministic; requesting no rows still verifies the saved evidence population.

Read `source_status`, `trace_status`, `query_status`, their control errors and
named missing inputs/files/batches independently. A successful inspection does
not turn a cancelled or control-error analysis into a completed analysis, nor
turn an unobserved input into zero support. `saved_artifacts_verified` describes
saved-artifact validation only. `original_reads_reverified` and
`reference_bases_reverified` remain false.

### Inspection Limits And Effects

The optional `budgets` field reuses `PrimerPayloadBudgets`. Defaults are exposed
by `primer_junction_inspection_default_budgets()`: 64 MB cumulative saved input,
8 MB response, 4,096 source records and queries, 1 MB lines and 128 alignments
per query. Inspection hard ceilings are 128 MB saved input, 8 MB response,
4,096 records/queries and 16 source batches and references. Each of the three
sidecars—evidence, partner profiles and breakpoint clusters—is independently
limited to `budgets.max_records` rows. Preview limits are at most 1,000 rows
and an offset at most 4,096. Retained dependencies count against the cumulative
input budget even when the display is filtered. Known report/sidecar/query
manifest bytes are reserved before trace validation, which receives the
remaining metadata byte allowance; the prepared-query validator then checks
the complete retained dependency total before consuming its data artifacts. Oversized output or evidence fails
explicitly; only the requested preview can be truncated.

The existing callback is checked before work, during streamed verification and
before returning. Saved-query validation is a bounded synchronous phase;
this is not a native responsiveness or large-cohort throughput guarantee.
Cancellation and failures return no partial inspector report. The operation
writes no files and changes no primer-order readiness or project content.
Shared-shell inspection leaves project history unchanged; generic engine
apply keeps its ordinary journal without a content revision or undo checkpoint.
Managed execution uses the existing expensive observation lane without a project
snapshot. Generic MCP confirmation policy is unchanged.

Native inspector forms, candidate-specific flank rescoring, abundance,
constrained redesign and conclusions about which molecule initiated RT remain
separate work.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Primer-Free Payload Evidence](primer-payload.md) | [Next: Observed Transcript Structures And Junctions](observed-transcripts.md)
<!-- manual-page-nav:end -->
