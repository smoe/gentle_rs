# cDNA SNP And Local Haplotype Comparison

The `--compare-cdna` mode of the existing `allele-hash-screen` family compares
explicit cDNA base/quality evidence with a bounded genomic reference. It reports
each sample separately, projects differences into every selected transcript,
and compares candidate bases between samples. It does not assume that a sample
expresses exactly two distinguishable alleles.

Current worktree: `codex/cdna-allele-comparison`, rebased onto main
`249173ee3b35bfcf6102f3b93b39feb357df553e`. Glen's supplied real fixture is pinned at
[`42e819c9`](https://github.com/smoe-bot/gentle_rs/tree/42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649/test_files/fixtures/cdna_haplotype_tp53_tp73). It contains TP53/TP73
barcode15/SK-MEL-28 and barcode16/WM793 primary-only cDNA selections against
GRCh38/Ensembl 116. The initial implementation and its pre-rebase verification
at `8df29b7d` remain historical receipts below, not acceptance of the current
partition and transitive-phasing changes. See the current analysis and source
receipt in [TP53/TP73 analysis receipt](tp53_tp73_cdna_fixture_analysis_20261110.md).
The scope is TP53/TP73; MDM4 is not included.

## Commands And Execution

```sh
gentle_cli allele-hash-screen --compare-cdna @tp53-comparison.request.json
gentle_cli shell 'rna-reads allele-hash-screen --compare-cdna @tp53-comparison.request.json'
```

The shared command accepts one inline JSON value or `@FILE`, lowering to
`CompareCdnaAlleles { request }`. Its optional result is
`OpResult.cdna_allele_comparison`; shell/direct CLI output exposes the report at
`result.cdna_allele_comparison`. Generic MCP `op`, GUI Shell and the JS/Lua/Python
shell wrappers use this same contract. There is no dedicated native GUI form.

The service verifies files and uses shared saved-transcript geometry before
calling the pure `gentle_engine::allele_hash_screen::cdna_comparison` kernel.
The normalized kernel request/report use the distinct
`gentle.cdna_allele_comparison_evidence_request.v1` and
`gentle.cdna_allele_comparison_evidence.v1` schemas; the latter is nested at
`comparison` in the file-bound report.
Existing v3 allele-hash screening keeps its supplied-variant/hash semantics.
The new mode reads explicit files and returns portable JSON on stdout. It writes
no files, requires no loaded project, and changes no project contents or history.
The managed command policy is an expensive read-only observation, with no
artifact-publication exception or implicit MCP state save.

## File Request

The strict schema is `gentle.cdna_allele_comparison_request.v1`. Unknown fields
are refused. Required fields are:

| Field | Meaning |
|---|---|
| `gene` | Human-readable gene label; not an implicit selector. |
| `source_id`, `gene_id` | Exact saved reference-gene identities; select the stable gene in the normalized annotation. |
| `reference.assembly`, `.contig`, `.start_0based` | Exact reference context and genomic offset of the supplied locus sequence. No contig aliases are guessed. |
| `reference.genomic_reference_sha256` | Declared parent reference hash, cross-bound to saved analysis and annotation. |
| `reference.locus_fasta` | Bound uncompressed single-record FASTA. Its header ID must equal `contig`. |
| `normalized_annotation` | Bound `gentle.primer_payload_annotation.v1` JSON with ascending genomic exons and explicit strands. |
| `samples` | One to sixteen biological-sample records, each with `{sample_id,label,analysis,basecalls}` and optional disjoint partitions below. Sample IDs and saved analyses must be distinct. |

Every bound file record is `{path,sha256,bytes}`, using exact raw file bytes and
a `sha256:`-prefixed digest. Each `samples[].analysis` is the unchanged
`gentle.observed_sam_structures_result_binding.v1` record:
`{schema,request_sha256,report_content_sha256,report_path,outcome}`. Preserve the
existing report and all artifacts it declares, including their relative layout.
The existing saved-report inspector verifies those artifacts before consumption.

A sample may add `additional_partitions: [{analysis,basecalls}, ...]`, with at
most 64 partitions including its required first pair. Multiple partitions require
`partition_source_sha256`, the canonical original BAM digest. Every saved
conversion receipt must declare that same BAM source. Stable producer unit IDs
name original BAM records, rather than local chunk ordinals; repeated original
source/unit keys and repeated raw saved units are rejected. The preparer verifies
unique original read IDs, full BAM/FASTQ populations and unchanged bases/qualities.
No QNAME-only merge, read multiplication or silent deduplication occurs.

All units from the sample's verified partitions enter **one** SNP discovery and
local-phasing pass. Input, base, unit, projection, work and output budgets remain
cumulative; partitioning does not remove a model or downsample reads. Outer
`inputs` retain each `partition_index`, individual bindings and completeness
counts. Original BAM identity and original-record ownership remain producer
declarations cross-bound to receipts, not independent proof of molecule identity.

Optional `supplied_sites` is a list of reference-forward **genomic zero-based
positions**, not transcript coordinates or a catalog of prescribed ALT alleles.
The routine also searches actual supported, quality-qualified observations for
candidate substitutions; `discovery_assessed_positions` describes that observed
scope. An empty list does not mean all genomic positions were callable. Discovery
and checking use the same reads and do not constitute independent confirmation.

Optional `phase_marker_positions` scopes only the phase graph, templates,
assignments and switch question. It must be a nonempty, unique list of genomic
zero-based positions within the reference and `supplied_sites`, with at most
`max_phase_sites_per_block` markers (hard ceiling 64). The normalized basis is
sorted and echoed in the core report. Omitting it admits all reported eligible
sites. It never filters discovery, per-site calls, transcript annotation, evidence
or reads. A panel selected from these same RNA reads is data-driven exploratory,
not prospective validation; using that same panel at Q30 is a sensitivity check.

Optional `policy` and `budgets` have explicit pilot defaults below.
`max_input_bytes` defaults to 64,000,000 (hard ceiling 128,000,000), cumulatively
including saved report artifacts. `max_basecall_bases` defaults to 2,000,000
(hard ceiling 16,000,000). The outer response is capped at 8,000,000 bytes.
The locus FASTA is capped at 16,000,000 bytes; all selected transcripts must fit
inside its interval. There must be 1–256 selected transcripts. Resource failure
or cancellation returns an error, not a shortened negative result.

Raw locus FASTA bytes and its normalized uppercase sequence digest are verified
separately. The full parent genome is not reopened: its hash is a declaration
cross-bound to the saved analysis/annotation, not a new full-genome verification.

## Original Bases And Qualities

Each sample supplies a separate, immutable sidecar:

```json
{
  "schema": "gentle.cdna_allele_basecalls.v1",
  "analysis_context_sha256": "sha256:REPLACE_WITH_SAVED_ANALYSIS_CONTEXT",
  "orientation": "original_sequencing_5_to_3",
  "provenance": "REPLACE_WITH_PRODUCER_VERSION_INPUT_HASHES_AND_EXTRACTION_METHOD",
  "reads": [
    {
      "unit": {"source_sha256": "sha256:REPLACE_SOURCE", "unit_id": "REPLACE_UNIT"},
      "read_id": "REPLACE_EXACT_SAVED_QNAME",
      "sequence": "ACGT",
      "qualities_phred33": "IIII"
    }
  ]
}
```

This is a structural example, not a usable sample or biological default. Copy
the unit key from the saved analysis. Each sequence must cover that unit's
entire original unclipped query axis; its `read_id` must match the saved QNAME
witnesses. Qualities, when present, are exact ASCII Phred+33 bytes in the same
original sequencing order and must have the same length. `null` preserves
unavailable qualities and cannot provide a quality-qualified allele call.
Retain genuinely different supplied basecall alternatives for a unit; duplicates
are refused. Missing unit evidence stays unassessed.

No original SAM/BAM/FASTQ is automatically reopened or joined by QNAME. Saved
SAM reports retain geometry and length witnesses, not the original SEQ/QUAL
strings. The existing read-end original-bases v1 sidecar remains unchanged and
does not gain qualities. A producer must explicitly establish the new sidecar's
original-record binding; valid hashes alone do not prove that provenance.

## Calls, Local Haplotypes And Sample Differences

Only supported aligned blocks are projected, using the existing shared query
geometry. Original and reference-forward bases, original query positions,
quality and alignment orientation remain in evidence. Reverse-aligned bases are
complemented once for reference comparison. Transcript annotations use each
transcript's strand and report mature cDNA positions **one-based**; a site absent
from that transcript remains `not_present`. Protein consequences, HGVS annotation,
indels and clinical classification are not implemented.

A definite unit/site call requires quality-qualified agreement across every
eligible retained mapping/basecall alternative and assessment of every supplied
query view. This is conditional on the supplied hypotheses, not proof that all
mapping alternatives or supplementary views were supplied.
Possible bases, conflicts, missing qualities and reasons remain explicit.
`N` skips, deletions, clips, masks, missing alternatives and uncovered positions
never become reference bases. The saved exclusion/mask declaration is required
for a base to be admitted. Counts are distinct producer-bound units per sample
and site; technical reads are not declared independent molecules. Search/view
completeness stays explicit and is not proven by a producer declaration.
Each input reports `inference_scope: supplied_eligible_mapping_hypotheses` and
separate counts for incomplete search declarations, incomplete view declarations
and unassessed supplied views. A primary-only selection may have every supplied
view assessed while both completeness declarations remain false. Omitted
alternatives, including recorded missing SA companions, remain limitations.

The local phaser forms blocks from supported co-observed marker links. It can
join repeatedly witnessed pair and larger partial-pattern families transitively;
a single read need not cover the entire marker set. Each local link needs the
recorded minimum distinct-unit support. Weak bridges do not silently join blocks.
Shared-allele couplings must be identifiable: unique or one-to-many star joins
can extend a supported local pattern; ambiguous many-to-many couplings remain
unresolved. Connectivity alone never proves a unique phase or K. Disconnected
components are not concatenated, and no missing second allele is synthesized.

Each candidate records `phase_method` (`direct_complete_vectors` or
`transitive_supported_local_patterns`) and any `support_links`, carrying
`sites_0based`, `alleles` and the distinct witnessing `unit_ids` for each necessary
local pattern. `exact_seed_units` still lists actual whole-vector reads; it can
be empty for a transitively supported candidate. Link witnesses are nonadditive:
several links may share a read, and their union is not an assembled molecule.
Partial units may receive bounded error-tolerant assignments; ties, unsupported
patterns and unresolved K remain visible. H1/H2/H3 names are local to a
sample/block; compare exact reference/marker/base fingerprints, never local labels.
The earlier complete-vector-only pilot is superseded by this transitive method.

Same-hypothesis observed skips accompany haplotype support. A switch-like marker
pattern can have several explanations; retain their marker brackets, competing
whole-template explanations and junction evidence. A bracket is not an exact
breakpoint. These are uncalibrated candidate observations, not a demonstration
of cross-allelic splicing, in-vivo recombination or trans-splicing. RT/PCR
switching, editing, mismapping and extra expressed sequence classes remain
competing explanations.

More than two supported templates sets an excess-expressed-haplotype candidate
flag in that connected block, with `copy_number_status: not_assessed`. RNA data
alone do not establish DNA copy number or copies per cell. One observed class
does not exclude an unexpressed allele; two classes do not exclude duplication.
Separate sample candidate-base differences retain callable denominators; lack
of support in the other sample is not genomic absence.

### Uncalibrated Pilot Policy

| Policy field | Default |
|---|---:|
| `minimum_base_quality` / `minimum_mapping_quality` | 20 / 20 |
| `minimum_alternate_units` / `minimum_alternate_fraction` | 3 / 0.10 |
| `minimum_haplotype_units` / `minimum_assignment_sites` | 3 / 2 |
| `maximum_assignment_mismatches` / `maximum_assignment_phred_penalty` | 1 / 40 |
| `maximum_haplotypes` | 8 |

These are inspectable software pilot settings, not calibrated variant-calling
or biological thresholds. SNP, phase-K and switch error models need independent
review and appropriate controls. No allele-frequency detection limit or
binomial confidence claim is made from correlated PCR reads.

Assignments rank threshold-compatible templates by summed mismatch penalty,
then mismatch count; equal best assignments remain ambiguous. A mismatch uses
the unit call's minimum qualifying Phred, clamped to 1–40. Missing sites supply
no match or mismatch. Switch scoring applies the same callable-site costs to
whole and mosaic templates, including disagreements at shared markers, and
requires support after leaving the inspected unit out: whole-vector seed support
for a direct candidate, or minimum support for **every necessary proof link** of
a transitive candidate. Allele-switch markers and splice junctions must occur on
that same original unit and retained mapping hypothesis. Transitive phase joins
never manufacture a read carrying another read's splice event. These scores are
inspectable heuristics, not genotype or switch probabilities.

Default kernel budgets are 10,000 units, 512 reported sites, 5,000,000 base
observations, 128 hypotheses per unit, 64 phase sites per block, 100,000,000
charged comparisons and 10,000 switch explanations. The file service additionally
caps projected base observations at 2,000,000 and mapping/basecall combinations
at 64 per unit. The kernel's default output ceiling is 64,000,000 bytes, with the
stricter 8,000,000-byte outer file-service response ceiling still applying.
Limits bound work; reaching one cannot prove the absence of
another allele, haplotype or candidate switch.

## Glen's Asset Bundle: TP53, Then TP73

Glen's real assets are now available in the pinned [fixture](https://github.com/smoe-bot/gentle_rs/tree/42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649/test_files/fixtures/cdna_haplotype_tp53_tp73).
All 1,129 selected BAM records were checked against the supplied FASTQ bases and
qualities in original orientation; all 18 listed raw-file checksums match. The
cohorts have TP53 482/481 and TP73 156/10 reads (barcode15/barcode16), with 40/20
selected transcript models. These are regular-cDNA libraries without TP53/TP73
target primers. The original unpartitioned TP53 producer hit its existing
projection-memory guard; deterministic 32-record preparation kept every record
and model in 38 bounded partitions, without weakening the guard. This does not
close the separate general exon-locus/native-capacity acceptance work.

The reproducible adapter is `scripts/prepare_cdna_fixture_comparison.py` with an
explicit source-matched CLI and external pysam installation. Its source/quality,
ordinal and population checks precede engine annotation preparation and saved-SAM
production. It creates only a fresh requested output directory; original fixture
files remain unchanged. The inventory/audit and interpretation belong with the
current source receipt, not the historical test counts below.

Prepare two independent gene requests using the same declared reference and
analysis policies. Start with TP53 as the characterization/control case; only
then interpret TP73. Each request contains these samples separately:

| Sample ID | Label |
|---|---|
| `barcode15` | `SK-MEL-28` |
| `barcode16` | `WM793` |

For **each gene and sample**, provide the exact saved observed-SAM result binding,
completed report and every declared artifact, plus the new base/quality sidecar
bound to its saved analysis context. Retain all query views, optimal mapping
hypotheses, masks, completeness/selection declarations and alternatives. Do not
preselect a favorable mapping or merge the cell lines. If the saved analysis
lacks complete view or exclusion declarations, retain that limitation rather
than repairing evidence by inventing a declaration.

For **each gene**, provide one shared normalized annotation with the exact
`source_id`/stable `gene_id` used by both analyses, plus one bounded single-record
reference FASTA whose contig interval covers every selected transcript. Include
assembly, contig, genomic start, raw byte counts/digests, the parent genome SHA,
annotation release/provenance and all request JSON. Original files must preserve
their raw hash-bound bytes; changing newline conventions requires new bindings.
For a new dataset, provide local accessible paths or a documented authorized transfer method.

Include a small provenance manifest naming extraction software/versions,
basecaller, sequencing platform, library/RT/PCR preparation, read selection and
any barcode assignment, deduplication/UMI or matched-DNA information available.
Keep missing controls explicit. Same-sample reviewed SNP/DNA truth can guide
`supplied_sites` and validation, but this slice does not ingest a phased VCF or
Sniffles result. Installed Sniffles remains optional for a separate matched-DNA
SV/CNV evidence branch; it is not called by this routine.

At the tested exact SHA, retain commands/platform/features, input and report
hashes, pass/fail/skip counts and limits/policy. Review per-sample callability,
transcript projections, local candidate patterns and cross-sample differences;
report every available switch candidate and unresolved block. Require independent
scientific review before interpreting an excess-haplotype or switch finding.
The supplied fixture contains only selected primary alignments. Search/view
completeness declarations remain false; omitted SA companions and low TP73 depth
in WM793 remain explicit. Candidate counts describe RNA/cDNA evidence against
GRCh38, not germline/somatic calls or RNA-editing classification. In particular,
14 pooled TP73 candidate sites are not 14 supported cell-line differences:
report each sample's callability, depth and base support separately.

## Verification Ledger

### Historical Initial-Implementation Receipt At 8df29b7d

The following receipt predates the rebase, real-data partition intake and
transitive-phasing correction. It does not accept those later changes or their
scientific results.

Tested locally on 2026-11-10, macOS/Darwin arm64, Rust
`1.100.0-beta.13 (3b0202dcd 2026-11-07)` and Cargo
`1.100.0-beta.13 (3acc79c0c 2026-11-05)`. The initial `8df29b7d` commit plus the
uncommitted source/contract changes were tested, not a new committed SHA.
The canonical changed-source manifest (baseline plus 28 path/byte-count/SHA-256
records, sorted paths and compact sorted-key JSON) hashes to
`e76673ef7f3e056fdb02a246ca8bf562ed07caf3b669c77a72e200145fffa050`.
The local run record is `/private/tmp/gentle-cdna-verification-20261110/source.json`;
the digest excludes its final newline. It covers changed Rust sources/tests,
the glossary, compiled CLI manual and new CLI replay script; unchanged sources/dependencies are
bound by the starting commit. Cargo manifests and lockfile are unchanged.

The owner authorized local execution after `cargo clean`. Cleaning the fresh
worktree's isolated `target` removed zero files; no other worktree target was
cleaned. Every Cargo command below used `--locked --offline`.

| Check | Result |
|---|---|
| `cargo check -q` with default features | Pass. |
| `cargo test -p gentle-engine --lib` | 82 passed, including 17 new cDNA tests and legacy allele-screen tests. |
| `cargo test --lib cdna_comparison` with default features | 14 passed: 12 ingress plus discovery and execution-policy checks. |
| Same default-feature libtest, filter `sam_intake_` | 12 passed, including byte binding, LF/CRLF, unit identity and cancellation. |
| Same libtest, filter `observed_structures_s2a_` | 4 passed: frozen bytes/IDs/work/order, context and errors. |
| Same libtest, filter `read_end_service_` | 12 passed, including original bases, saved-only reads and PNG/PDF preservation. |
| `cargo build --no-default-features --bin gentle_cli` | Pass; explicit headless CLI. |
| `cargo test --no-default-features --test capability_registry_parity cdna_comparison` | 1 passed, checking actual MCP tool discovery and read-only metadata. |
| `GENTLE_CDNA_CLI=.../target/debug/gentle_cli python3 -m unittest scripts.test_cdna_allele_comparison -v` | 2 passed: genuine producer artifacts, direct/shared parity, synthetic TP53/TP73 display labels, LF/CRLF, quoted Unicode paths, qualities/reference refusal and unchanged files/state. |
| `python3 -m unittest scripts.test_capability_maturity` | 17 passed. |
| `python3 -m unittest scripts.test_primer_payload_contract` | 9 passed, 1 expected optional legacy CLI replay skip (`GENTLE_PAYLOAD_CLI` unset). |
| `python3 -m unittest scripts.test_primer_trace_contract` | 7 passed. |
| `python3 -m unittest scripts.test_tutorial_checkouts` | 87 passed, including LF/CRLF projections and negative controls. |

Totals are **125 unique Rust tests and 122 Python tests passed, one Python test
skipped**; reruns are not added again. The default-feature libtest was
`target/debug/build/GENtle/4094ae3ac2f87558/out/gentle-4094ae3ac2f87558`, SHA-256
`c88e3a5672efe9265b860bff86392199f4eb33646dc23e170edb85c092915a4c`.
The headless CLI SHA-256 was
`f3de1fe39028cd3aa23b89535b0d0e510a8fb9951edf55bdfaa926a0a123f5f0`.

Initial failures were corrected before the final green run: a missing test
trait import, a macOS temporary-path alias, and a fixture that changed producer
prose instead of semantic qualification settings. Production intake guards and
context identities were preserved. Initial expected generated-projection
freshness failures were resolved through the canonical protocol renderer and
`scripts/generate_capability_maturity.py --write`; the four maturity changes
only refresh the glossary digest, without elevating evidence verdicts.
Builds emitted headless unused-code warnings in unchanged sections and the
macOS linker's large unwind-table warning; these were not treated as runtime
or scientific failures. Scoped formatting, whitespace and documentation-status
checks passed. Session-close warnings concern the intentional uncommitted diff,
required generated projections and manual plan-fidelity review.
Manual diff review confirms the bounded scope described in the implementation
ledger; broader statistical inference and biological acceptance remain open.

At that historical receipt, native Windows, full-suite CI, real-sample execution,
statistical calibration and independent scientific acceptance remained open.
Subsequent real-data work is separately source-bound; it does not inherit this
receipt or an earlier Glen receipt as scientific acceptance. The
existing saved-artifact inspector has no caller cancellation callback;
cancellation is rechecked after that bounded verification finishes.

See the [original proposal and implementation ledger](rna_cdna_haplotype_plan.md)
and the [optional independent review prompt](rna_cdna_haplotype_claude_review_prompt.md).

### Current Post-Rebase Verification Follow-Up

The validated analysis baseline was main `249173ee3b35bfcf6102f3b93b39feb357df553e` plus the
recorded working-diff source/contract manifest
`b55a282dfff02de1ebd8aa9535db27a5fa4fb25a37c3d16d46957b04130ea5db`
(30 sorted path/byte-count/SHA records, compact sorted-key JSON, excluding its
final LF). The complete [current analysis receipt](tp53_tp73_cdna_fixture_analysis_20261110.md#verification-and-artifact-receipt)
binds this source, the final CLI/libtest hashes and every real-data output.

Current-source checks passed: **140 unique Rust tests**, **45 Python tests**,
with one expected optional legacy payload-CLI skip. The Rust scope is 91 engine
library tests, 20 headless comparison/ingress tests, 12 SAM-intake, 4 frozen-S2A,
12 read-end service tests and one actual MCP discovery test. Python covers 17
maturity, 9 payload-contract (+1 skip), 7 trace-contract, 2 actual cDNA CLI and
10 fixture-preparer tests. The default-feature `cargo check -q --locked --offline`
passed. Cargo builds/tests used locked/offline dependencies; no new core
package dependency was added. These totals replace no earlier source receipt.

The six-record whole-versus-three-partitions control passed: stable original
IDs, forward/reverse bases and original qualities, six distinct units, site
observations, phase evidence and two separate same-read skip sets are semantically
identical after excluding only provenance digests. Every real original record
was independently checked against its BAM/FASTQ and global-ordinal proof. All
1,129 records and 40/20 transcript models remain in 38 partitions.

All eight final full/panel Q20/Q30 runs completed. Their per-site calls,
annotations, original-base/quality evidence and sample counts are exactly
unchanged by the transitive-phasing correction. Independent witness checks
verified 520 complete-seed and 50 link-witness memberships across repeated
analyses; those are not counts of distinct molecules. All actual witness units
support the recorded qualified allele patterns. See the numerical report for
local classes, weak links, undercoverage and the absence of an informative
cross-splicing test.

Native Windows, full-suite CI, native GUI, calibrated error/K inference,
matched-DNA corroboration and independent scientific acceptance remain open.
Headless unused-code and macOS unwind-table linker warnings were retained;
no warning was promoted to a scientific failure. The preserved unpartitioned
TP53 memory refusal was resolved through bounded composition, not a weakened
cap or changed read/model universe. An intermediate check during active source
editing failed on a then-missing initializer; final stable-source checks pass.
Session-close reports no failures and two expected reminders for intentional
uncommitted/generated files and manual plan fidelity. Scope review confirms
TP53/TP73 only, with MDM4 excluded. Claude was not invoked.
