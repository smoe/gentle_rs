# cDNA SNP Evidence, Local Haplotypes And Allele-Linked Splicing

Current status, 2026-11-10: the owner subsequently authorized a bounded
implementation and scoped local tests after `cargo clean`. See the
[implementation ledger](#implementation-ledger-2026-11-10) and
[implemented file contract](cdna_allele_comparison.md). The worktree is now rebased
onto main `249173ee3b35bfcf6102f3b93b39feb357df553e` and Glen's real TP53/TP73 fixture is pinned
at `42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649`. See the post-rebase follow-up at the end of
this document. Claude review remains optional, not a prerequisite. The original
proposal and initial implementation receipt below are dated history; the initial
complete-vector-only phasing requirement has been superseded by supported
transitive local phasing, without joining disconnected or ambiguous components.

Original proposal status, 2026-11-10: owner-requested proposal. No implementation, read analysis,
package installation or new execution campaign is authorized by this plan.
Inspected HEAD: `7fa6b15adc28d0adfe3e1edce0190a2621b54705`. The pending MCP
roadmap/changelog correction and unrelated `outputs/` are preserved. No local
Rust builds or tests run. An optional [Claude review prompt](rna_cdna_haplotype_claude_review_prompt.md)
is prepared at the owner's request; no Claude invocation or waiting is required.

## Biological Questions And Interpretation

1. Which supplied SNP alleles are observed in each cDNA read, and which sites
   are actually callable? Separate a supported reference base from an unobserved
   site, skipped exon, failed quality filter or ambiguous alignment.
2. Which SNP combinations occur together on reads within a sample and connected
   phase block? Infer expressed sequence haplotypes without forcing two clusters.
3. Which splice chains/events are associated with those haplotypes? Separately
   identify individual reads with a change between established haplotype patterns.
4. Are more than two reproducible, distinguishable sequence haplotypes supported
   in the same marker block? Report an extra-haplotype finding and the evidence
   needed to investigate copy number, subclones, editing or artifacts.

For this proposal, a **haplotype-switch candidate** is one read with ordered SNP
evidence consistent with one established haplotype on one side and another on
the other side. A **junction-associated/cross-allelic splicing candidate** also
requires an actually supported junction in the marker bracket on the same
mapping hypothesis; a switch without such a junction remains reportable under
the generic label. This differs from two alleles expressing different splice patterns. Even a
recurrent switch is not proof that two RNA molecules were joined in vivo.
Reverse transcriptase can create apparently trans-spliced cDNA; PCR can also
produce additional apparent haplotypes. Report the observation before assigning
a mechanism. [RT template-switching study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0012271),
[long-read PCR phasing study](https://www.nature.com/articles/srep21746).

More than two error-tolerant RNA clusters is **a flag for further investigation,
not a copy-number call**. RNA editing and aneuploidy are both relevant to
multi-haplotype transcript analysis, and errors can create false variant sets.
[FLAIR2](https://link.springer.com/article/10.1186/s13059-024-03301-y).
Two simple counterexamples set the inference boundary: an A/A/B duplication
still has only two distinguishable haplotypes; a mixture of diploid subclones
can have three sequence haplotypes without every cell carrying three copies.
RNA read fractions are not DNA copy fractions. A one-haplotype result does not
exclude another unexpressed or insufficiently sampled allele.

## Existing Owners And Missing Primitives

This combines existing evidence services with new allele-observation and phasing
primitives; it is not an adapter-only addition. Preserve [DEC-027](decisions.md#dec-027-rna-allele-evidence-does-not-invent-phase):
unphased variants and disconnected blocks cannot become an invented gene-wide
haplotype. A physical splice-chain identity remains independent of its SNPs.

| Reuse | Boundary or required extension |
|---|---|
| `src/observed_transcript_structures/{sam_report,core}.rs` | Saved unit hypotheses, supported blocks, skips, MAPQ, orientation and exact chain IDs; keep frozen contracts unchanged. |
| `geometry.rs::project` | Project a reference-forward site through a supported block to the original query axis. Complement a reverse-aligned observed base exactly once for reference-forward comparison. |
| `read_ends/oriented.rs::{Document,OriginalRead,prepare}` | Explicit context/unit/read-name/full-axis binding and sequence alternatives. Its v1 document has **no base qualities**. |
| `sam.rs::Row`, `sam_text.rs`, saved SAM row assessments | Retain SEQ/QUAL lengths, not their strings. Old saved reports alone cannot supply SNP bases or quality. Never reconstruct either from a reference or `=/X` operations. |
| `indexed_locus.rs` and `src/rna_read_ends.rs` | Configured bounded alignment selection, saved-input validation, cancellation, input protection and separate publication. Do not add another BAM parser or aligner. |
| `crates/gentle-engine/src/allele_hash_screen.rs` | Reviewed variant-table/VCF-plus-transcript-map intake, reference checks and supplied block-local allele evidence. Existing h1/h2 k-mer screening is not an arbitrary-K phaser or proof of a site's genotype. |

Keep new biology in a portable engine/helper module, with a separate versioned
report. GUI/CLI/MCP/JS/Lua/Python consume that same report. Do not reinterpret
existing S1/SAM/Phase 2 schemas, change their digest recipes, add SNPs to
`CoreFamily`, or alter current UTR/coverage/pileup denominators.

## Inputs And Counting Units

- Explicit sample, biological source, barcode, library/RT/PCR preparation,
  sequencing platform/basecaller and available UMI/deduplication metadata.
  Analyze samples separately; pooling cell lines or donors can manufacture an
  apparent excess of haplotypes. Technical reads are not independent molecules.
- Bound saved analysis/data plus reference assembly, contig dictionary and
  reference FASTA/locus bytes. Validate every variant's reference allele against
  that exact authority; retain explicit transcript projections and isoform IDs.
  No guessed contig aliases, allele strands or reference-version substitution.
- Initially a supplied SNP catalog: preferably same-sample DNA genotypes with
  phase-set IDs, or a reviewed SNV table. A population catalog is candidate
  annotation, not this sample's genotype. Distinguish DNA-confirmed, externally
  phased, RNA-only candidate and editing-suspect sites.
- Explicit base-and-quality evidence associated with each saved unit and
  sequence alternative: original sequencing orientation, full axis, sequence
  SHA, quality encoding and provenance. Prototype a separate immutable sidecar
  or explicit extractor receipt; preserve the existing original-bases v1 format.
  Missing qualities remain absent, not invented high scores. Such observations
  may be displayed but cannot silently enter a quality-qualified inference tier.
- Optional matched DNA/phased VCF, RNA-editing annotation, independently prepared
  RNA libraries and external DNA copy-number/SV evidence, each with its own
  hashes, sample identity and unavailable status.

Identity is initially `(analysis/source SHA, unit_id)` plus observation and
basecall-alternative IDs. QNAME is a display/check field. Do not join independent
analyses, loci or trimmed queries merely by QNAME. Cross-analysis deduplication
requires an explicit original-record linkage manifest. UMI-supported molecule
counts are optional, separately reported; identical cDNA sequence/endpoints
alone do not justify collapsing independent molecules.

## Staged Delivery

### H0 — Input Contract And Known-SNP Callability

Implement the smallest first slice: a bounded single-locus, single-sample
consumer of saved geometry, explicit bases/qualities and a supplied SNP list.
No de novo variant caller, clustering, switch mechanism or CNV inference yet.

For every unit/site retain constituent observation rows and all mapping/basecall
alternatives. Allele inference uses only the saved unit's retained optimal
hypotheses; nonoptimal/rejected observations remain audit evidence and cannot
create a haplotype or inflate support. A definite unit-level reference/alternate
call requires callable agreement across every eligible mapping and basecall
alternative. An unresolved eligible branch prevents a clean negative or definite
assignment; preserve positive possible support separately. Observation states are:

| Observation | Meaning |
|---|---|
| `reference_supported` | Quality-qualified observed reference base at a supported genomic position. |
| `alternate_supported` | Quality-qualified observed catalog alternate, with its allele ID. |
| `other_base` | Quality-qualified observed canonical base outside the supplied allele set; retain it for review. |
| `ambiguous` | Conflicting admissible observations or insufficient base identity/assignment. |
| `unassessed` | Uncovered, exon skipped under that hypothesis, `D`, clipped, missing bases/qualities, masked/low-quality or unsupported mapping. Preserve the exact reason. |

An `N` skip, deletion, absent exon or missing base is never a reference call.
Keep a separate callable consensus/status across alternatives, without choosing
the best-looking allele. Count a unit at most once per site; possible-allele
sets and their overlapping support are not additive molecule counts.

Reference-forward genomic coordinates are the canonical SNP basis; retain
original-query coordinates, observed literal base, quality and orientation.
Use existing supported `M/=/X` geometry; do not convert `N` into `D`, strip clips
or fill sequence to accommodate a DNA caller. Record local indel/homopolymer,
junction-proximity, primer/adapter and ambiguous-paralog context as filters or
warnings with explicit policy. A primer-derived base must not become a genomic
allele witness.

Emit reference/alternate/other/unassessed counts and denominators per sample,
site and declared cohort. A site with zero alternate support is reported as
**not observed among N callable units**, with sampling uncertainty, not absent
from the genome. DNA homozygosity and variant absence are not inferred from
RNA-only lack of observations. Always show raw callable-unit denominators, but
calculate a detection limit only under an explicit sampling/dependence model.
Do not treat PCR-correlated reads as independent binomial trials; if molecule
independence or mapping-bias control is unavailable, quantitative detection
power remains unassessed rather than receiving a deceptively tight interval.

### H0D — Optional RNA-Only SNP Discovery

If no adequate same-sample SNP catalog exists, benchmark an RNA-aware SNV
discovery route, such as isoLASER, against declared truth before admitting it
as an optional external producer. Restrict the initial scope to substitutions;
retain the exact supported reference bases searched, callable coverage and
uncallable/masked intervals. A catalog-only run must not be labelled de novo
discovery or a complete negative screen.

Import candidate sites and per-read evidence through the same H0 contract,
with tool/model/version, reference/input hashes and parameters. Review splice-
alignment, indel/homopolymer and reference-mapping artifacts, and annotate
editing-suspect sites separately. Prefer independent DNA or library validation;
without it, keep the label RNA sequence-variant candidate, not proven genomic
SNP. Use held-out validation where feasible and disclose discovery/validation
overlap. No stock DNA caller is presumed valid on spliced cDNA without an
RNA-specific benchmark. H0D does not change the original alignment or its CIGAR
semantics. Its catalog can then feed H1, with RNA-only and DNA-validated results
kept distinguishable.

### H1 — Connected, Error-Aware Local Haplotypes

Construct a sparse unit-by-site allele/quality matrix. Initially phase eligible
supplied heterozygous SNVs; keep RNA-only or editing-suspect variation in a
separate exploratory layer rather than treating every mismatch as a DNA allele.
Build phase connectivity from co-observed sites. Split weak/disconnected blocks
and preserve separate phase-set identities; optional trusted external DNA phase
can bridge them only through an explicit provenance binding.

Use a bounded likelihood/consensus model with an explicit error profile and
candidate K values including K > 2. Assess missing observations only where
callable: missing positions are not reference alleles and shared isoform labels
are not phase evidence. Compare clusters on enough shared informative sites;
avoid transitive single-linkage merging through partly observed bridge reads.
Retain posterior/compatibility alternatives, unresolved units, a noise/residual
class, and phase-block-specific stable H1/H2/H3 labels. They are not parental
labels without a separate parental/DNA anchor.

Require distinguishable marker patterns, calibrated support and stability under
site/error-profile sensitivity checks. Connectivity alone does not identify K
when sites are missing. Require jointly informative shared markers or explicit
identifiable constraints across the component; otherwise leave K unresolved.
Test whether a one/two-haplotype error model still explains the data rather than
counting local cluster labels as additional haplotypes. Record underpowered versus supported K,
minor-cluster support and compatible K ranges; never turn an exhausted search
cap into evidence that only two haplotypes exist. Where data allow, separate
discovery from validation by independent library or declared read family. Do
not select variants/clusters and then claim independent confirmation from the
same evidence without disclosing that dependence.

RNA-aware read-backed phasing and haplotype/splicing association are established
approaches; benchmark suitable tools against a shared truth corpus rather than
assuming one fits GENtle's contracts. [phASER](https://pmc.ncbi.nlm.nih.gov/articles/PMC5025529/),
[isoLASER](https://www.nature.com/articles/s41467-025-64605-6).
FLAIR2 is also a useful multi-haplotype benchmark. Do not add a backend or new
dependency until its data/error/ploidy assumptions and deterministic receipts
are reviewed. The existing allele-hash screen remains an independent screen.

### H2 — Haplotype/Splice Associations And Switch Candidates

Join local haplotype evidence to the **same unit and same compatible mapping
hypothesis** as the saved splice geometry. Never combine an allele from one
alternative with a junction from another into an artificial transcript. Produce
haplotype-by-chain and haplotype-by-junction/exonic-event tables with callable,
unassigned, ambiguous and excluded units retained. Report descriptive counts
first; biological association tests need declared replication, dependence and
multiple-testing policy, not an assumed 50:50 RNA ratio.

Keep ordinary allele-associated alternative splicing separate from a physical
cross-allelic switch. Search within reads for a change between independently
established local marker patterns, initially labelled `haplotype_switch_candidate`.
A junction-associated label additionally requires supported splice geometry
inside the marker bracket on that same hypothesis, with its flank qualification
retained. Never invent a junction between separate alignment segments. A switch
without a qualifying junction stays in the generic candidate table. Compare a no-switch explanation with an
error-aware switch explanation, accounting for the number of positions searched.
Use leave-one-read/family-out templates or an external phased anchor so that the
candidate does not create the haplotypes used to validate itself. Retain RNA
editing, a third consistent haplotype and DNA rearrangement as competing models.

Report **every screened switch-like observation** with an evidence tier and
reason: low-information, ambiguous, artifact-suspect, supported candidate or
independently corroborated candidate. A lone discordant base cannot justify a
high-confidence switch; the planned stronger tier needs multiple informative
markers on each side (at least two per side as a proposed structural minimum),
plus calibrated quality/support criteria. We do not silently discard weaker
observations; they remain in the audit table and cannot count as confirmed events.

For each candidate retain:

- Original unit/observation/basecall IDs, library/sample and distinct read/family
  support; independent-library recurrence is a separate count.
- Ordered marker calls, qualities, competing haplotypes, score/uncertainty,
  callable sites and reasons for exclusions.
- The interval between the last informative marker on one side and the first
  on the other, in original-query and genomic coordinates. Do not invent an
  exact breakpoint across an unobserved interval or disconnected phase block.
- Any splice junctions bracketed by that interval, exact chain IDs, supported
  blocks, alignment alternatives and overlap with repeats, primers or clips.
  Genomic position ordering and mRNA-direction presentation stay distinct.
- Artifact flags: RT/PCR template switching, concatemer/adapters, repeated
  endpoints, mismapping/paralog, low-complexity/indel errors, editing, phase
  errors, sample mixing and possible genomic recombination/rearrangement.

Validation should use separately prepared libraries/RT controls where available,
UMI-aware evidence if present, an orthogonal RNA-level method not sharing the
same RT/PCR dependence, and matched DNA to test genomic alternatives. Replicate
PCRs from the same cDNA are not independent RT confirmation. A junction-adjacent
switch is reportable but not itself proof of trans-splicing. In isoLASER,
"trans-directed" splicing refers to regulation; it is not a claim of joining
two allelic RNA molecules.

### H3 — Extra-Haplotype Findings And Optional DNA Evidence

Count supported sequence classes only within one sample and overlapping,
connected, jointly informative marker block. Never sum the two labels from
several disconnected blocks or count splice isoforms as extra haplotypes.
Emit `excess_expressed_haplotype_diversity` only when a calibrated, stable model
supports more than two identifiable classes on a shared marker basis or a
sufficiently constrained component; otherwise K remains unresolved. Include K
uncertainty, marker provenance and read/library support. Separate findings based on DNA-validated
markers from RNA-only/editing-sensitive patterns. Preserve below-threshold or
unresolved extra clusters as such; do not relabel them as absent.

The RNA-only report keeps `copy_number_status: not_assessed`. Extra classes can
motivate CNV/subclonality investigation but cannot determine copies per cell.
Likewise a strong allelic RNA imbalance is not evidence of a particular DNA
dosage. Correlate with separately supplied DNA read-depth/BAF, genotypes,
SV/CNV calls or a targeted DNA copy assay, with sample/purity/normalization
provenance. Preserve the external method's own uncertainty and scope.

**Sniffles branch:** a reproducible Debian installation could help if matched
genomic long-read BAM/CRAM is available. Sniffles is a structural-variant caller,
not the cDNA SNP caller or arbitrary-K RNA phaser. It can contribute external
duplication/rearrangement evidence; it does not convert RNA cluster count into
copy number, and spliced RNA alignments are not its default DNA assay input.
A negative SV call is not a complete negative CNV test; retain DNA coverage and
assay callability separately.
[Official Sniffles documentation](https://github.com/fritzsedlazeck/Sniffles),
[Sniffles2 protocol](https://www.nature.com/articles/s41596-026-01362-w).
Record exact package/tool/dependency versions, parameters, reference, inputs and
VCF hashes. A specific Debian package/version has not been verified or installed
here. Keep this optional and outside the core dependency graph.

### H4 — Discovery, Inspection And An Auditable Pilot

After the private kernel is tested, add one strict engine-owned request/report
and thin shared-shell/generic adapter routes. Names are proposals, not currently
executable commands. Prefer saved-report inspection and approved fresh-directory
publication over hidden session state; do not mutate source/project evidence.
Cancellation/control failures must remain distinct from empty or absent findings.

Proposed output artifacts, all versioned and request/source/policy-bound:

| Artifact | Purpose |
|---|---|
| `site_observations.tsv` | Unit/site/base/quality/callability and all alternatives/exclusions. |
| `site_summary.tsv` | Explicit per-site cohort/callable denominators and allele counts. |
| `phase_blocks.tsv`, `haplotypes.tsv`, `read_assignments.tsv` | Connected scope, consensus alleles, K uncertainty, support and alternatives. |
| `haplotype_splicing.tsv` | Same-hypothesis allele-to-chain/event linkage and eligibility denominators. |
| `switch_candidates.tsv` | Every screened switch-like observation, including insufficient/artifact-suspect cases, marker brackets and provenance. |
| `extra_haplotype_findings.tsv` | More-than-two findings, competing explanations and separate DNA-evidence status. |
| `report.json` and optional locus SVG | One portable summary plus coverage, marker, haplotype and junction inspection; no GUI-owned calculation. |

Raw sequences/qualities and identifiable read names should not be re-exported
by default. Evidence IDs and bounded requested previews are sufficient for
ordinary inspection; explicit data export follows the existing approval rules.

## Error Calibration And Deterministic Acceptance

No universal biological quality, support, minor-cluster fraction or switch-score
cutoff is selected by this plan. A versioned policy must state base/MAP quality,
context masks, minimum co-observed markers, read/family support, error model,
K/model selection, stability and tested-site multiplicity. Calibrate on synthetic
truth plus an appropriate declared control, then freeze the pilot policy before
looking for target-specific events. Show sensitivity results without selecting
the most interesting outcome. Mapping bias is an explicit control, including
allele-swap/personalized-reference sensitivity on an independently recorded
alignment branch when needed. [WASP mapping-bias study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4626402/).

Author deterministic synthetic cases for:

1. Reference/alternate/other versus unobserved, spliced-out, deleted, clipped,
   masked, low-quality and unknown-quality sites; no false reference filling.
2. Both genomic strands/read orientations, exact query coordinates, junctions,
   nearby indels/homopolymers, primer-derived bases and ambiguous paralogs.
3. Repeated QNAMEs, multiple mappings/basecall alternatives, duplicate record
   bindings, UMIs and separate analyses/samples; unchanged original denominators.
4. One/two/three true sequence haplotypes plus calibrated errors and editing;
   minor clusters, imbalance, partial coverage, weak bridges, incompatible
   clusters, deterministic tie ordering and disconnected phase blocks.
5. Allele-specific splice usage with no switch; same-read A-to-B and B-to-A
   marker switches; no jump assembled across mapping alternatives; candidates
   with zero/one/multiple junctions in their marker bracket.
6. RT/PCR-like chimeras, adapter/concatemer and mismapping controls; high-scoring
   artifacts must remain reportable suspects rather than confirmed biology.
7. A/A/B duplication with two RNA haplotypes, diploid-clone mixture with three,
   editing-generated extra classes and an unexpressed allele. None yields an
   RNA-derived absolute copy-number claim.
8. Raw hash/REF/assembly/quality-axis refusal, cancellation/control errors,
   LF/CRLF, portable paths, output protection, no-clobber and resource boundaries.

Proposed first-slice resource admission is single-locus and explicit: at most
10,000 units, 512 supplied SNV sites, 1,000,000 admitted unit/site observations,
64 sites per inference component, eight candidate haplotypes, 100 million
charged comparisons, 64 MB input and an 8 MB inline response. These are reviewable
engineering ceilings, not biological detection thresholds. Lower caller limits
apply cumulatively; stream evidence sidecars, poll cancellation and refuse or
label an unresolved bound instead of silently truncating counts, candidates or K.
No exponential full-gene enumeration or unbounded whole-genome materialization.

Execution belongs in the next authorized consolidated Glen/CI scope after each
implemented slice. Keep exact source, package/model versions, hashes, commands,
features/platform and nonzero pass/fail/skip counts. Add actual adapter/CLI and
native/platform checks only when those boundaries exist. Do not reissue completed
MCP/parity, containment, BLAST/reader or private 15/16 campaigns as unrun.

## Recommended Next Decision

Approve H0 alone first: explicit known-SNP read evidence and callability, with
quality-aware binding and synthetic tests. Before an actual-data pilot, obtain
the exact sample/locus/reference/read-unit manifest, original bases/qualities,
library preparation and a decision on matched DNA availability. No data transfer
is needed just to review this plan. H0D is the separate RNA-aware discovery
route when a supplied catalog is insufficient; without it, H0/H1 assess supplied
sites only and do not claim all variants have been found.

H1 follows only after the observation contract is reliable; H2 adds allele/splice
and switch reporting; H3 adds extra-haplotype triage with optional DNA corroboration;
H4 exposes the proven records. All candidate switches and excess-haplotype
findings remain inspectable even when mechanistic or copy-number confirmation
is unavailable. Experimental trans-splicing, DNA CNV and scientific acceptance
remain independent outcomes, not statuses conferred by software completion.

## Preparation Receipt

The plan and optional prompt are prepared without Claude execution, package
installation, read analysis or local Rust builds/tests. Documentation-status
and whitespace checks pass. Session-close reports four OK, two expected
worktree/manual-plan reminders, zero failures. The roadmap remains 400 lines;
one existing RT paragraph was reflowed solely to admit the new proposal line.
The earlier pending MCP-roadmap correction and unrelated `outputs/` are preserved.
These checks concern document consistency, not algorithmic or scientific validity.

## Implementation Ledger (2026-11-10)

Subsequent owner instruction: implement a reusable comparison routine now,
starting with TP53 and then TP73, with SK-MEL-28/barcode15 and WM793/barcode16
reported separately and compared. Glen's advice was to extend the existing
`allele-hash-screen` owner at a clean exact pushed SHA, preserve bases/qualities
and alternatives, and treat independent Claude review as a later review gate.
At this initial implementation stage, real read/reference assets were not yet
accessible, and the owner asked to be asked again once Glen prepared them.
The later delivery and analysis are recorded separately below. Sniffles is already installed; no installation or invocation is needed
for this RNA comparison slice.

Implementation starts from `8df29b7dcc1da8bc14f4b5fdc9a7c945661af77a` in a clean
`codex/cdna-allele-comparison` worktree. The owner lifted the previous local Rust
restriction for this task and required `cargo clean` first; the fresh worktree
clean removed zero bytes. This authorization supersedes the original proposal's
planning-only/no-local-Rust status, without broadening unrelated test campaigns.

The bounded implementation composes one new pure kernel under the existing
`gentle_engine::allele_hash_screen` module, a saved-evidence/file-verification
service, and the existing family's `--compare-cdna` mode. It includes observed
substitution nomination, per-unit quality-aware callability, strand-aware
transcript positions, conservative connected-block candidate templates,
error-tolerant assignments, same-hypothesis observed-junction support, switch
candidate explanations and per-sample candidate-base comparisons. It returns
portable read-only JSON, without new file publication or source/project changes.
The original v3 hash-screen, observed-transcript/SAM and original-bases v1
contracts are preserved.

This initial implementation was narrower than the proposed H0–H4 research
workflow. Its repeated-complete-vector seed requirement was too conservative for
transitively linked partial reads and is superseded by the post-rebase method
below. The initial receipt applies to that earlier method. The initial Q20/MAPQ20,
three-unit/10% alternate and three-unit template defaults, with at most one
assignment mismatch, are explicit **uncalibrated pilot policy**. The routine
does not implement calibrated likelihood/model selection, independent discovery
validation, UMI-aware molecule inference, paired DNA/Sniffles import, protein
consequences, indels or clinical calls. Candidate switches and excess haplotypes
are not trans-splicing mechanisms or CNV diagnoses; no second allele is forced.
The original proposal's broader controls and scientific acceptance remain open.

The [implementation contract](cdna_allele_comparison.md) records strict request,
base/quality sidecar, reference/source/hash identities, budgets, denominator and
orientation rules, plus Glen's exact asset checklist. Raw files and supplied
locus sequence are reverified; parent-genome identity is declared/cross-bound
to saved analysis and annotation, without rereading the full genome. No implicit
join of original reads by QNAME is added.

Historical pre-rebase verification at `8df29b7d`: **125 unique Rust and 122 Python tests passed locally**, with one
expected optional legacy CLI skip. This includes all 82 engine tests (17 new),
14 comparison ingress/routing/policy tests, 28 saved-intake/frozen/read-end
preservation tests, one actual MCP discovery test and two new process-level CLI
tests. The default-feature locked/offline check passed; actual CLI replay used
the source-matched headless binary. Canonical parity and the four necessary
maturity digest projections are current. The [full ledger](cdna_allele_comparison.md#verification-ledger)
records the baseline plus source binding, binary hashes, commands/features,
corrected fixture failures, warnings and limits. No synthetic or earlier
source receipt establishes actual TP53/TP73 results, native Windows behavior,
full-suite acceptance or scientific validity. Claude has not been invoked;
the [original review prompt and implementation supplement](rna_cdna_haplotype_claude_review_prompt.md)
are ready for the owner to pass along without waiting.


## Post-Rebase Real-Data And Transitive-Phase Follow-Up (2026-11-10)

The implementation was rebased onto main `249173ee3b35bfcf6102f3b93b39feb357df553e`, preserving
main's CF01 and exon-locus review history. Glen supplied the pinned
[real fixture](https://github.com/smoe-bot/gentle_rs/tree/42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649/test_files/fixtures/cdna_haplotype_tp53_tp73) at `42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649`.
The scope remains TP53/TP73; MDM4 is excluded. The 1,129 original BAM records
match supplied FASTQ bases/qualities in original orientation, with TP53 482/481
and TP73 156/10 selected reads for barcode15/barcode16. All 18 listed raw file
checksums match; these are regular-cDNA primary-only selections, not full mapping
universes or matched DNA. All 40 TP53 and 20 TP73 transcript models remain.

The unpartitioned real TP53 producer hit its existing projection-memory guard.
Deterministic bounded preparation now uses 38 disjoint partitions, preserving
all original records, qualities and original-BAM/global-ordinal proofs. Stable
source/unit identities reject reuse; the shared service collects each sample's
parts before one comparison/phasing pass. It does not merge separately inferred
haplotypes, silently deduplicate QNAMEs, weaken budgets or remove models. A small
whole-versus-partitioned control checks semantic equivalence beyond provenance,
including reverse qualities and same-read splice sets. The preparer suite has
10 passing tests. The final source-bound receipt records 140 unique Rust and
45 Python passes, one expected optional legacy CLI skip, the actual semantic
partition control and eight completed real-data replays.

Glen's review correctly required transitive local phase rather than a read
spanning all markers. The revised engine uses supported pair and larger partial-
pattern families, requiring the minimum distinct-unit support per necessary
link. Identifiable shared-allele joins can propagate phase; ambiguous many-to-
many coupling remains unresolved, and weak/disconnected components stay separate.
Candidates expose their phase method and nonadditive read-witness support links.
Real complete-vector seed lists may be empty. Switch leave-one-unit-out support
is checked for every necessary link, and allele/splice evidence always remains
on one original read and retained hypothesis.

Optional explicit `phase_marker_positions` scopes only phase, leaving discovery,
all per-site calls, annotations, read counts and original evidence intact. The
TP53 and TP73 exploratory panels were selected post hoc from the same Q20 RNA
observations, using the recorded exonic/support rule, and reused unchanged at
Q30. They are not prospectively fixed or independently validated. The broad
baseline and data-derived sensitivity runs must both remain inspectable.
Candidate counts describe RNA/cDNA evidence against GRCh38; they do not classify
germline/somatic variation or RNA editing. TP73 site counts are not a count of
robust cell-line differences: each sample's callability, depth and actual base
support must be stated, especially for ten-read WM793/TP73. The source-bound
numerical report is [TP53/TP73 analysis receipt](tp53_tp73_cdna_fixture_analysis_20261110.md).

Current verification: default-feature locked/offline check, 140 unique Rust
tests and 45 Python tests passed, with one expected optional legacy CLI skip.
The [current receipt](tp53_tp73_cdna_fixture_analysis_20261110.md#verification-and-artifact-receipt)
records exact source/CLI/libtest hashes, scoped commands, the full/partitioned
control, original-unit preservation, all eight real-data runs and witness
validation. Earlier 125/122 totals and earlier complete-vector outputs remain
historical. Native Windows, full CI, statistical/error calibration, orthogonal
DNA and independent scientific acceptance remain separate. Claude has not been
invoked; the prepared read-only prompt is available without blocking work.
