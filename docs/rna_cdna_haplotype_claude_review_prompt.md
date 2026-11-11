# Optional Read-Only Claude Review: cDNA SNPs And Local Haplotypes

Prepared 2026-11-10 at the owner's request. This prompt has not been run;
the plan is not Claude-reviewed. Do not wait for a review to deliver the plan.
The initial prompt below is preserved as proposal history. For the subsequently
authorized implementation, pass the [latest real-data/transitive-phase supplement](#real-data-and-transitive-phase-review-supplement-2026-11-10); the earlier implementation supplement is dated history.

If separately invoked, follow the repository's noninteractive read-only form:

```text
claude -p --permission-mode plan --allowed-tools "Read Grep Glob" --output-format text -- '<paste the prompt below>'
```

No shell/write tools, Rust builds/tests, package installation, model-provider
setup, private-data analysis or auth/token inspection. Do not use bypass modes.

## Prompt To Pass

```text
Please give a concise read-only scientific and architecture critique of this
Codex plan, followed by a revised minimal plan. Treat it as a proposal, not an
implemented or experimentally validated workflow.

Repository: /Users/u005069/.codex/worktrees/faae/gentle_rs
Inspected HEAD: 7fa6b15adc28d0adfe3e1edce0190a2621b54705.
Read AGENTS.md and the current plan before reviewing. If files have changed,
identify the actual differences instead of assuming this snapshot is current.
Do not modify files or run builds/tests/tools outside Read, Grep and Glob.

Owner task:
Test cDNA reads for SNP presence/absence, derive error-tolerant haplotypes,
report potential cross-splicing between alleles, and investigate more than two
haplotypes as a possible copy-number signal. The owner also asked whether a
Debian Sniffles package would help. This is planning only; no private read input
or matched DNA data has been supplied for execution.

Original Codex plan:
- H0: quality-aware per-unit/site observations and callability at supplied SNVs.
  Reference support differs from uncovered/spliced-out/deleted/clipped/low-quality
  or ambiguous evidence. Inference uses retained optimal hypotheses, with callable
  agreement across all eligible mapping/basecall alternatives for definite calls.
  Keep nonoptimal observations audit-only and preserve original denominators.
- H0D: if no adequate catalog exists, benchmark an optional RNA-aware de novo
  SNV producer and import its sites through H0. RNA-only variants remain distinct
  from DNA-confirmed SNPs and editing-suspect sites.
- H1: read-backed, connected local phase blocks and error-aware arbitrary-K
  clustering with uncertainty/noise, without forcing diploidy or bridging gaps.
  Connectedness is not sufficient to identify K with missing sites; shared-marker
  evidence or identifiable component constraints are required, else K is unresolved.
- H2: join haplotypes to the same read/mapping hypothesis's splice chain/events.
  Report ordinary allele-associated splicing separately from within-read A-to-B
  marker switches. Call generic switches haplotype_switch_candidate; a splicing
  association requires a supported same-hypothesis junction in the marker bracket.
  Preserve all screened switch-like candidates, including weak
  and artifact-suspect observations. Stronger tiers need multiple markers on
  each side and calibrated controls; do not claim an exact breakpoint in gaps.
- H3: report stable excess expressed haplotype diversity, not RNA-derived copy
  number. Evaluate RNA editing, sample/subclone mixtures, paralogs and RT/PCR
  chimeras. Use independent RNA-level and matched DNA evidence when available.
  Sniffles is optional for genomic-DNA structural-variant corroboration, not a
  cDNA SNP caller, primary RNA phaser or copy-number result from RNA clusters.
- H4: one portable engine-owned report with thin adapters, bounded sidecars and
  saved-result inspection; preserve existing schema bytes, chain IDs and counts.
- Implement only H0 first after owner approval; no production cutoffs, package
  dependency, runtime command, release gate or scientific acceptance is chosen.

Full proposal:
docs/rna_cdna_haplotype_plan.md

Existing evidence and relevant seams:
- docs/decisions.md, DEC-027: unphased genotypes/disconnected phase sets must not
  become a synthetic gene-wide haplotype; read/fragment/evidence units differ.
- src/observed_transcript_structures/{core,sam_report,geometry}.rs:
  saved unit hypotheses, supported blocks/skips and reference/query projection.
  geometry::project already handles reverse query axes; allele comparison is
  reference-forward and must complement a reverse-aligned base exactly once.
- sam.rs::Row, sam_text.rs and saved SAM assessments retain SEQ/QUAL lengths,
  not base/quality strings. Existing reports are not sufficient SNP inputs.
- read_ends/oriented.rs::OriginalRead/Document/prepare binds explicit original
  sequences to analysis/unit/QNAME/full axis but has no base qualities. Preserve
  sequence alternatives. Plan needs explicit quality-aware input/receipt.
- src/rna_read_ends.rs and indexed_locus.rs provide saved-input and configured
  bounded alignment-selection/publication seams, not automatic FASTQ joins.
- crates/gentle-engine/src/allele_hash_screen.rs supports reviewed VCF/transcript
  maps and supplied phase-block k-mer evidence; it does not infer arbitrary-K
  read-backed haplotypes. Do not silently expand its h1/h2 semantics.
- docs/roadmap.md and docs/rna_read_end_integration_plan.md: oriented pileup is
  implemented but its new runtime scope is external. No local Rust execution.
- docs/glen_next_verification_handoff.md: preserve completed MCP/parity,
  containment and reader receipts; add only changed-slice checks later.

Primary research linked in the plan: phASER, isoLASER, FLAIR2, WASP, the
Houseley/Tollervey RT template-switching study, amplicon phasing pitfalls and
official Sniffles documentation. Do not confuse isoLASER's trans-regulation
terminology with physical cross-allelic RNA trans-splicing.

Please focus on:
1. Missing evidence or circular inference in SNP calls, local phasing, cluster-K
   selection and switch detection, especially partial/nonoverlapping transcripts.
2. Controls that distinguish true alleles from RNA editing, systematic error,
   mapping bias, mixed samples/clones and RT/PCR artifacts.
3. Whether a bulk RNA finding can support the proposed interpretation without
   overstating CNV, molecule independence or a trans-splicing mechanism.
4. Smallest reusable engine boundary, required quality/base identity contract,
   feasible resource limits and meaningful synthetic/orthogonal validation.
5. A revised minimal first slice and the essential decisions still needed.

Return concrete blockers first, then concise improvements and a revised minimal
plan. Do not redesign unrelated modules, choose biological defaults by fiat,
request a broad new data campaign or imply that your review ran any tests.
```

## Implementation Review Supplement (2026-11-10)

Historical initial-implementation context; the later supplement supersedes its
complete-vector-only requirement and unavailable-assets status.

```text
The owner subsequently authorized implementation and scoped local tests after
cargo clean. Review the current diff read-only; do not apply the original
planning-only recommendation as a claim that implementation is unauthorized.

Repository/worktree: /private/tmp/gentle-cdna-allele-comparison-20261110
Exact clean starting SHA: 8df29b7dcc1da8bc14f4b5fdc9a7c945661af77a.

Original plan and dated implementation ledger:
docs/rna_cdna_haplotype_plan.md
Current contract and asset checklist: docs/cdna_allele_comparison.md
Code: crates/gentle-engine/src/allele_hash_screen/cdna_comparison.rs,
src/allele_hash_screen/comparison.rs and their tests; existing family parser in
src/engine_shell/command_parsers.rs, read-only execution/policy and discovery.

Current minimal implementation extends the existing allele-screen owner with a
pure kernel and one bound-file ingress service. It verifies saved observed-SAM
artifacts, explicit original bases/Phred+33 qualities, exact annotation/source
gene bindings and a bounded locus FASTA; parent-genome SHA is cross-bound but
the full genome is not reread. Shared supported-block geometry preserves
orientation, masks, alternatives, query axes and distinct read-unit counts.

The kernel nominates observed substitutions, retains callable/ambiguous/missing
evidence, projects each site into individual transcripts and compares samples.
Connected-block templates require repeated complete marker vectors; partial
units receive bounded error-tolerant assignments. K remains unresolved without
adequate complete vectors/support. Policy is explicitly an uncalibrated pilot
(Q20/MAPQ20; alternate >=3 units and >=10%; template >=3; assignment <=1 mismatch).
It does not force exactly two alleles or join disconnected blocks. Same-
hypothesis observed junction support and switch explanations remain candidates.
Extra templates do not establish DNA copy number or a trans-splicing mechanism.
No protein/indel/clinical calls, second pipeline, implicit original-read join or
file/project mutation. Existing v3/SAM/original-bases schemas remain unchanged.

Real TP53/TP73 barcode15 SK-MEL-28 and barcode16 WM793 inputs are unavailable;
Glen will prepare exact bound assets. No actual sample findings are claimed.
Sniffles is already installed, but not invoked for this RNA-only slice.

Please identify concrete statistical/logic/provenance blockers, especially
variant selection using the same evidence, incomplete marker vectors,
uncalibrated errors/K, seed support, switch-template circularity, phase block
boundaries, cross-hypothesis mixing and denominator preservation. Suggest the
smallest corrections and essential independent controls. Separate source-only
findings from runtime/scientific acceptance. Return a concise critique plus a
revised minimal plan; do not write, build, test or install anything.
```


## Real-Data And Transitive-Phase Review Supplement (2026-11-10)

```text
Read-only critique; do not modify files, execute shell/build/tests, install tools,
or infer that this review establishes runtime or biological acceptance.
Repository: /private/tmp/gentle-cdna-allele-comparison-20261110
Rebased main baseline: 249173ee3b35bfcf6102f3b93b39feb357df553e
Real fixture: https://github.com/smoe-bot/gentle_rs/tree/42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649/test_files/fixtures/cdna_haplotype_tp53_tp73
Pinned fixture SHA: 42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649

Task: compare real TP53 and TP73 RNA/cDNA differences against GRCh38 for
barcode15/SK-MEL-28 and barcode16/WM793, annotate per-transcript alleles,
derive supported local haplotypes, and report switch/excess-pattern candidates.
MDM4 is out of scope. RNA candidates are not germline/somatic, editing, CNV or
trans-splicing diagnoses.

Read docs/cdna_allele_comparison.md, docs/decisions.md DEC-027,
docs/rna_cdna_haplotype_plan.md and the current source-bound analysis receipt
in docs/tp53_tp73_cdna_fixture_analysis_20261110.md. Relevant code and tests:
crates/gentle-engine/src/allele_hash_screen/cdna_comparison.rs,
src/allele_hash_screen/comparison.rs, scripts/prepare_cdna_fixture_comparison.py
and their adjacent tests. Original Codex plan is preserved above and in the plan
file; do not apply its original planning-only restriction to authorized work.

Original Codex plan: H0 preserves quality-qualified per-unit/site callability;
optional H0D imports or nominates RNA-only SNVs; H1 infers error-aware arbitrary-K
local phase without joining disconnected components; H2 links the same read's
alleles and splice events; H3 reports excess RNA-pattern candidates with matched
DNA needed for copy-number inference; H4 exposes one portable shared report.
The first pilot narrowed H1 to complete-vector seeds; Glen requested the current
transitive correction. The original proposal/history remains readable above.

Current minimal revised plan/implementation:
- Verify explicit BAM/FASTQ source populations, original orientations and Phred
  bytes. All 1,129 supplied records match; the primary-only selection keeps both
  completeness declarations false, separately from assessment of supplied views.
- Prepare the full 40/20 annotation models through existing engine GTF intake.
  Preserve the unpartitioned memory refusal. Use 38 bounded disjoint saved-SAM
  partitions, stable original-BAM/global-record IDs and one unified sample/core
  pass; do not downsample, weaken guards or merge separately phased outputs.
- Extend the existing phaser with supported pair/larger partial-pattern families.
  Repeated shared-marker links can propagate phase without any whole-vector read.
  Unique/star couplings are allowed; many-to-many ambiguity, weak bridges and
  disconnected components remain unresolved/separate. Proof links retain their
  distinct read witnesses; exact_seed_units remains actual full-vector evidence.
- Preserve arbitrary K, raw alternatives and denominators. Check switch-template
  support leave-one-read-out on every required proof link. Allele/splice evidence
  must stay on one original read and one retained mapping hypothesis.
- Explicit phase_marker_positions changes only the phase question. Preserve all
  discovery, calls, annotations and reads. Panels are post hoc/data-driven from
  these same Q20 RNA reads; using the identical panel at Q30 is sensitivity, not
  prospective fixation or independent confirmation. Preserve broad baselines.
- Report per-sample callability/depth/base support before comparing cell lines;
  TP73 has 156 versus 10 selected reads, so pooled 14-candidate wording must never
  imply 14 robust differences. No missing second allele or negative CNV/switch
  conclusion may be manufactured by low power or an unresolved block.

Critique concrete identifiability/statistical/provenance bugs first, especially:
transitive joins and weak-bridge boundaries; missing/error-correlated sites;
K and minor-template sensitivity; same-read splice linkage; partition equivalence
and identity; support-link leave-one-out circularity; and data-driven selection.
Distinguish software conformance, candidate observations and scientific claims.
Return a concise critique and a revised minimal plan with the smallest necessary
corrections. Do not widen scope.
```
