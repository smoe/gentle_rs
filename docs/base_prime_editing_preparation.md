# Base And Prime Editing: Preparation Brief

Status: preparation only, 2026-10-20. No editing capability, executable protocol,
approved editor profile, active tutorial or new release gate is introduced.
The implementation sequence below is proposed, not authorized by its existence.

## Scope And Source Identity

Inspected local `main` and `gentle_rs_2_main` at
`050c1030ab133ecaebc1da9fdcfd1527d04e8c49`, after fast-forwarding main from
`09639034f83de6f985c761c5871b8b886639da84`. No fetch was performed; this names
the inspected local source, not an assertion about the current remote tip.
Work continues in the linked `gentle_rs_2_main` worktree. Existing `paper/`,
`output/` and unrelated untracked primary-checkout documents are untouched.

Initial scope is conventional, guide-directed nuclear-DNA research in cultured
mammalian cells. Other organisms, organelles, chemistries, multiplex/large edits
and complex editor systems need separate applicability review. Sequence design
does not establish delivery, editing efficiency, specificity, function or
therapeutic suitability. Laboratory and human scientific approval stay separate.

This source/literature review made no Rust build, model/provider invocation,
sequence upload, order, tag or publication. Browsing public references is not
contacting a service with a design. No external design/analysis tool was run.
An optional read-only Claude consultation was offered; none was invoked and
this brief must not be described as Claude-reviewed.

## Verified Gap Map

These are source-inspection findings, not newly executed runtime tests.

| Need | Existing capability and evidence | Reusable composition | Missing primitive / deferred work |
| --- | --- | --- | --- |
| Reference identity | `GenomeExtractionProvenance`, `BiologicalContext`, `GenomicRegionReference` and genomic half-open intervals in [protocol](protocol.md); [context](../crates/gentle-protocol/src/biological_context.rs) and [region types](../crates/gentle-protocol/src/genomic_regions.rs) | Bind one loaded sequence and its exact extraction mapping, reference/annotation identity and hashes | Strict original-to-desired edit intent, versioned contig, haplotype and stale-reference checks; existing optional context fields alone are insufficient |
| Guide intake | `GuideCandidate` / `UpsertGuideSet` in [engine protocol](../src/engine/protocol.rs), [normalizer](../src/engine.rs) and [guide handlers](../src/engine/ops/candidate_guides.rs) | Guide IDs, explicit coordinates and practical-filter diagnostics | Intake validates strings/interval order/cut-offset range, but does not prove the spacer/PAM matches loaded DNA or an editor; no editor-aware enumeration |
| Practical constraints / oligos | `FilterGuidesPractical`, `GenerateGuideOligos`, exports and [guide specification](rna_guides_spec.md) | Reuse diagnostic helpers and explicit oligo roles after biological validation | No CBE/ABE windows, bystanders, PBS/RTT or prime nick model. Built-in prefix/suffix templates and optional added 5' G are not pegRNA design |
| Guide workflows | [Filter workflow](examples/workflows/guides_filter_and_generate_oligos.json), [export workflow](examples/workflows/guides_export_csv_and_protocol.json), tutorial 04.04/04.05 projections | Existing parser/engine/adapter path and contract-test pattern | The `tp73` example labels are not evidence-bound TP73 targets; do not relabel these guide strings as a public editing example |
| Alleles / annotations | `MaterializeVariantAllele`, shared validator in [variant_promoter.rs](../src/engine/analysis/variant_promoter.rs), [explicit allele decision](decisions.md#explicit-variant-allele-selection); splicing/CDS interpretation | Fresh reference refusal, source-forward alleles and annotation projections | Current materialization is SNV-only and mutating, not a generic read-only indel-intent validator; no silent widening of its v1 contract |
| Sequence redesign | [Sequence-design proposal/apply](../src/engine/analysis/sequence_design.rs), [portable records](../crates/gentle-protocol/src/sequence_design.rs) | Content-bound proposal, separate approval, cancellation/budget and derived-lineage patterns | Synonymous constraint solving does not model editor chemistry, nicking, repair or measured editing outcomes |
| Cloning / samples | Existing oligo, construct, arrangement, container and lineage contracts in [architecture](architecture.md) | Explicit component records and later reviewed materialization/storage handoff | Editor-specific component roles, RNA chemistry and reviewed cloning compatibility; no automatic ordering/execution |
| Verification | `ConfirmConstructReads`, ABI/SCF traces, expected/baseline alignments and saved reports in [sequencing contract](sequencing_confirmation_plan.md) | Amplicon/primer context, confidence and local expected-edit inspection | No read-count editing-outcome importer, denominator ledger or clone-genotype inference from bulk frequencies |
| Off-target evidence | Prepared-reference BLAST and primer-specificity paths | Reference/tool provenance and bounded-search patterns only | PAM/bulge-aware editor assessment and guide-independent DNA/RNA evidence; generic BLAST/primer scores are not editing safety scores |
| Interfaces | Shared `Operation`, shell parser/executor, canonical reachability/readiness and GUI intents | One biological calculation and saved report across CLI/MCP/JS/Lua/Python/GUI/agents | No base/prime editing request or route found in the inspected `src/`/`crates/` paths; proposed names below are not discoverable capabilities |

## Short Literature And Correction Ledger

Public source pages checked 2026-10-20. The following compact records are
bibliographic/applicability notes, not retained copies of protocols. No PDF,
HTML payload, figure or experimental parameter set is committed here. Hence no
raw-document hash is claimed. Before a document/profile becomes executable
evidence, retain its permitted original bytes with URL/DOI, retrieval date,
version, license/access level, SHA-256, correction links and reviewed section/
figure locators. Hash a corrected document separately; never overwrite its
predecessor. Parsed sequence hashes and document-byte hashes are distinct.

### L1: Conventional Base Editing

[Huang, Newby and Liu](https://www.nature.com/articles/s41596-020-00450-9),
DOI `10.1038/s41596-020-00450-9`, online 2021-01-18, Nature Protocols 16,
1089-1128. Publisher access is a subscription preview, not full-protocol access.
Scope: CBE/ABE architectures, target/editor/delivery selection in mammalian
cells; the described experimental example is HEK293T. Do not infer a universal
delivery method or transfer parameters to another cell type from the abstract.
No experimental quantities or editor-window values are extracted here.
Accessible full protocol and exact construct identity remain prerequisites
before activating conventional CBE/ABE parameter profiles.

[Publisher correction](https://www.nature.com/articles/s41596-021-00525-1),
DOI `10.1038/s41596-021-00525-1`, 2021-11-17: figure-label legibility,
nucleotide alignment in Figure 2 and received-date corrections. The notice is
readable; it does not grant full-protocol access. Treat strand diagrams as
version-sensitive evidence, not authority to guess a window.

### L2: Prime-Editing Architecture And Workflow

[Doman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9799714/), DOI
`10.1038/s41596-022-00724-4`, online 2022-08-08, issue November 2022; full author
manuscript accessible in PMC. Covers Cas9(H840A)-RT/pegRNA systems, PE2/PE3
and later variants; RNA stabilization and protein architecture are separate
choices. Describes plasmid transfection in HEK293T and RNA electroporation in
primary human fibroblasts. Use component orientation, controls and validation
as design-contract evidence, not universally transferable efficacy/delivery
settings. No linked correction was identified on the consulted page; this is
not proof that no later correction exists. Recheck at profile review.

### L3: Cell-Specific Prime-Editing Applicability

[Wu et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10519846/), DOI
`10.1016/j.xpro.2023.102583`, online 2023-09-21, issue 2023-12-15; accessible
full text, CC BY-NC-ND. Human pluripotent stem cells; PE3 workflow with
cell-specific interventions and plasmid electroporation, followed by amplicon
analysis. Its editor variants and delivery context are not GENtle defaults.
Do not silently adopt p53/MMR interventions or infer safety from this example.

[2026 correction](https://pmc.ncbi.nlm.nih.gov/articles/PMC13158587/), DOI
`10.1016/j.xpro.2026.104582`, online 2026-05-04, issue 2026-06-19; open notice,
CC BY-NC-ND. Corrects Scheme 2 scaffold/extension cloning overhangs, including
their written directions. The PMC landing page intermittently presented a
browser check; its indexed primary notice and [PubMed erratum record](https://pubmed.ncbi.nlm.nih.gov/42085185/)
were available. No schematic bytes were retained or approved as a template.
The corrected notice must accompany any later use of that cloning scheme.

### L4: External Design Comparison

[Hsu et al., PrimeDesign](https://www.nature.com/articles/s41467-021-21337-7),
DOI `10.1038/s41467-021-21337-7`, 2021-02-15; open access. The 2021-03-03
update corrects author-name indexing, not a biological parameter. Describes
pegRNA/ngRNA design and demonstrations in human cells; not a delivery SOP.
The [official repository](https://github.com/pinellolab/PrimeDesign) provides
a local CLI with bounded PBS/RTT options and optional additional substitutions.
Compare a pinned commit/container digest with matched settings, rather than
reimplement every feature. Disable unrequested PAM-disrupting changes and
retain all external changes as reviewable data. No version is selected or run
here; tool ranking is not an efficacy or safety oracle.

### L5: Observed Outcomes, Not Design Feasibility

[Clement et al., CRISPResso2](https://www.nature.com/articles/s41587-019-0032-3),
DOI `10.1038/s41587-019-0032-3`, 2019-02-26, Nature Biotechnology 37, 224-226.
Publisher text was accessible through search; direct loading was intermittent.
Amplicon analysis includes base-editing and allele-aware quantification; this
is an analysis architecture, not a prescribed editor or delivery method. No
correction was established by this limited check; recheck before tool pinning.
The [official service](https://crispresso2.pinellolab.org/) is a reference only:
GENtle must not upload sequences implicitly. Prefer local, version-pinned
result import from the [official tool](https://github.com/pinellolab/CRISPResso2).

[Official parameter documentation](https://docs.crispresso.com/latest/parameters.html)
currently identifies itself as v2.3.4 and documents prime-editing component
orientation, reference overrides and quantification windows. It is a moving
page, not a pinned runtime. Do not attribute these later prime-editing options
to the 2019 paper, import its default thresholds silently, or use its
orientation-check override to make a GENtle validation pass.

## Protocol-Applicability Matrix

These are proposed applicability states, not current runtime enum variants.
`reference_only`, `needs_review`, `unsupported` and `reviewed_for_context` must
stay distinct. A citation alone cannot produce the last state.

| Evidence / choice | Permitted preparation use | Missing applicability binding | Initial disposition |
| --- | --- | --- | --- |
| L1 plus correction | Conventional CBE/ABE feasibility questions | Accessible full text, exact construct, PAM/spacer/window/context and cell/delivery review | `reference_only`; no activated window |
| Explicit conventional CBE profile | C-to-T on the declared exposed strand | Versioned construct/source and position-specific activity/context evidence | `needs_review`; no transversion/editor substitution |
| Explicit conventional ABE profile | A-to-G on the declared exposed strand | Same bindings, separately reviewed from CBE | `needs_review`; not a CBE preset renamed |
| L2 / conventional prime editor | PBS/RTT/nick reconstruction and strategy separation | Exact editor/scaffold/geometry and bounded design policy | `needs_review`; PE2 first, PE3/PE3b later |
| L3 plus corrected Scheme 2 | Cell-context/correction teaching and later cloning review | hPSC, delivery and intervention-specific human approval | `reference_only`; not a generic workflow |
| PrimeDesign / CRISPResso2 | Local interchange/comparison and saved outcome interpretation | Pinned versions, formats, parameters, inputs and complete artifact hashes | `needs_review`; no tool execution in preparation |
| Other organisms, organelles, nonconventional chemistries or therapeutic use | Record the request and missing evidence | Separate applicability brief and approval | `unsupported` in initial implementation |

## Proposed Typed Contract

All names in this section are **proposals only**, not Rust definitions, registered
schemas or callable commands. Prefer a small `gentle-protocol` module and one
shared engine analysis path, wired through existing `Operation`/shell adapters.
Reuse context/region primitives above without treating their optional metadata
as a completed editing-reference check. Do not extend `GuideCandidate` into an
unvalidated bag of optional editor fields or introduce another adapter registry.

| Proposed type / operation | Minimum typed content and invariant |
| --- | --- |
| `EditingReferenceBinding` | Biological context, assembly accession/version, contig accession/version, taxon, exact loaded/inline sequence and genomic-forward reference window with separate semantic SHA-256 values, local-to-genomic mapping and orientation, optional raw-source hash, explicit selected haplotype; synthetic references carry a distinct synthetic namespace |
| `GenomeEditIntent` | Reference binding, one genomic-forward `start_0based..end_0based_exclusive`, exact `deleted`/`inserted` A/C/G/T strings, original and reconstructed desired sequence digests; reject unchanged intent, reference mismatch, unresolved phase or mixed/ambiguous targets |
| `EvaluateGenomeEditIntent` / `GenomeEditIntentReport` | Read-only validated context, independent reconstruction check, structured refusals and unresolved annotation evidence; no persisted derived DNA or implicit fetch |
| `BaseEditorProfile` | Required profile ID/version/digest, exact editor construct, conversion, PAM pattern and side, spacer length, exposed-strand relation, position-origin and editable positions, source-specific sequence-context rules/evidence, applicability state and correction bindings; no universal window or default editor |
| `EvaluateBaseEditing` / `BaseEditingFeasibilityReport` | Intent/profile hashes; candidate spacer/PAM intervals and both strands, intended-position check, potential bystanders, transcript/CDS/splice annotation references with their release; practical diagnostics separate from feasibility, probabilities absent without a validated model |
| `PrimeEditorProfile` | Exact protein/scaffold identity, nicked strand and interbase geometry, RNA component order/alphabet, strategy-independent architecture and evidence; never the generic nuclease cut offset |
| `EvaluatePrimeEditing` / `PrimeEditingFeasibilityReport` | Intent/profile hashes; explicit spacer, scaffold, PBS, RTT, optional stabilizer and nicking guide; per-component sequence/direction/digest, nick coordinates, independently reconstructed edited strand and allele, bounds and incompatibilities |
| `EditingSearchSummary` | Explicit input/search limits, deterministic ordering/ties, examined and retained counts, `complete`/`budget_exhausted`/`cancelled`, truncation reason and omitted-count availability; incomplete search cannot claim no feasible design or best design |
| `EditingOutcomeImportRequest` / `EditingOutcomeReport` | Local tool/version/parameter and raw-file bindings, sample/replicate/control IDs, amplicon and reference/desired digests, allele sequences plus integer counts, classification policy and named denominators; preserve unknowns and source labels |
| Later `EditingComponentHandoff` / explicit apply | Frozen report/request digest and exact human-approved component set; refuse stale source, revalidate before materialization, preserve original and lineage/undo; file export and project mutation remain separately declared effects |

Keep separate report dimensions for design feasibility, off-target evidence,
laboratory execution, observed allele and functional validation. Initial
off-target status is `not_assessed`, not zero risk. Read-only evaluation must
leave project bytes, mutation revision and undo/redo unchanged; an operation
journal receipt does not by itself mean the project was mutated.

### Coordinates And Orientations

Canonical edits use genomic-forward half-open intervals. Substitution/deletion
requires `end-start == deleted.len()` and an exact reference match; an insertion
uses `start == end`, empty deleted string and an explicitly anchored interbase
position. Resolve the loaded extraction to a genomic-forward reference window
`R` and map absolute coordinates to window-local `s,e` before reconstructing
`R[..s] + inserted + R[e..]`. Never index a reverse-oriented local sequence with
absolute genomic coordinates. Bind both the original extraction and `R` hashes;
the report retains the transforms back to the loaded display. Preserve submitted
and canonical forms when normalizing repeat indels; never shift coordinates
without a recorded reference-bound transform. Synthetic examples below use
an explicitly zero-origin forward reference. Initially accept one SNV or one
explicitly bounded short indel; overlapping/multi-edit/haplotype combinations
require a later contract. Selecting a reference haplotype is a design basis,
not a claim about the genotype of the actual cell line.

Transcript/HGVS input must resolve through one exact versioned transcript and
exon map; no symbol-only resolution, guessed complement, assembly-family
fallback or silent liftover. Distinguish gene strand, local-sequence orientation,
spacer sequence direction, hybridized DNA strand and chemically edited strand.
For a local length `n`, reverse projection maps `[s,e)` to `[n-e,n-s)` and an
interbase position `k` to `n-k`; apply the actual extraction anchor separately.

On genomic-forward alleles, conventional CBE on exposed plus yields C-to-T;
on exposed minus it yields G-to-A. ABE yields A-to-G on exposed plus and
T-to-C on exposed minus. Do not conflate the guide-hybridized strand with
the exposed/editable strand. Profile positions must declare their origin
explicitly (for example, spacer 5' base = position 1), not reuse map labels.

First base-editing output lists all potentially editable bystander positions,
not invented combination probabilities. Annotation consequences name each
versioned transcript and available CDS phase/splice evidence. Unavailable
annotations stay unknown; a synonymous label is not a safety verdict.

Prime RNA components are represented in RNA 5'-to-3' order: spacer, scaffold,
RTT, PBS, then an optional explicitly chosen stabilizing element. Store DNA
encoding separately. Check PBS complementarity to the nick-exposed DNA 3' end
and that reverse transcription of the RTT plus flanking homology reconstructs
the declared allele on the correct strand. Reverse complement is not a display
toggle. Check both the original and length-changed coordinate frames for indels.
PE2 has no second nicking guide; PE3 requires a separately validated opposite-
strand nick. PE3b needs an explicit edited-versus-reference recognition rule,
not a promise of allele specificity. Protein architecture and RNA stabilization
remain orthogonal to that strategy. No extra PAM-disrupting change is implicit.

### Outcome Denominators And Non-Claims

Do not turn Sanger discrepancy rows or a pooled percentage into clone genotypes.
Import original allele counts before derived percentages, retaining input,
quality-passed, aligned, fully callable and excluded/unclassified counts with
reasons. Denominators and quantification windows are required, versioned data.
Missing counts and zero denominators yield unavailable fractions, not zero.

Define disjoint categories: intended-only, intended-plus-declared-bystander,
other substitutions, unexpected indels, unedited and unclassified. An intended
insertion/deletion itself is not an unexpected indel. Preserve orthogonal flags
such as scaffold incorporation and presence of all intended changes; define
precedence for combined outcomes instead of double-counting reads. Partial
coverage or ambiguous haplotype/alignment cannot establish an intended allele.
Require counts to reconcile within each declared denominator. Report untreated
controls separately; do not subtract background silently or treat reads as
independent biological replicates. Local amplicons cannot establish genome-wide
specificity, guide-independent DNA/RNA editing absence or biological function.

## Independently Reviewable Slices

Stop after each slice's own bounded evidence; do not freeze routine development
or add per-branch release/Glen campaigns. The first implementable request is S1,
not a combined editing workbench.

| Slice | Smallest deliverable and dependencies | Exit evidence / limitations |
| --- | --- | --- |
| S0: evidence/profile preparation | Obtain permitted full L1 text and correction; select one exact conventional CBE and one ABE construct for later review. Bind complete PAM/window/context provenance without guessing numbers; select one conventional prime architecture separately | Human applicability review and raw-byte source ledger. Until then use explicitly synthetic test profiles only, never biological presets |
| S1: reference-bound intent | Shared types, validator and read-only report for explicit sequence/anchor plus desired substitution/short indel; one generic-op route with canonical shell/parity/readiness tests | Independent literal reconstruction, both orientations, stale/ref/haplotype refusal, no mutation/files/network. No editor recommendation or GUI designer |
| S2: conventional base feasibility | S1 plus explicitly selected reviewed CBE/ABE profile (or clearly synthetic test profile); strand-correct candidate/PAM/window checks and bystander/annotation report | Exhaustive tiny cases, absent PAM/impossible conversion, annotation-unknown behavior and bounded-search status. No efficacy/off-target ranking |
| S3a: prime component validation | S1 plus one explicit profile and caller-supplied spacer/scaffold/PBS/RTT; read-only PE2 reconstruction before any enumeration | Independent strand/PBS/RTT/nick/indel reconstruction; fail closed on bad orientation, added edits or insufficient flanks |
| S3b: bounded prime design | S3a plus explicit finite search ranges, then separately PE3/PE3b nick-guide rules; local PrimeDesign interchange/comparison | Retained normalized candidates and mismatch explanations at pinned settings. External output never bypasses S3a; tool installation/execution requires later authorization |
| S4: inspectable component handoff | Reuse oligo/cloning/sample contracts; maps/GenBank/EMBL/JSON from saved reports with thin DNA-viewer inspection. No scoring in rendering | Targeting bases, promoter-added bases, overhangs, scaffold/RTT/PBS and modifications visibly distinct; read-only handoff first. Separately approved materialization rechecks exact digests, lineage, undo and stale-job rejection |
| S5: local observed-result import | Import version-pinned CRISPResso2/explicit allele-count exports against S1/reference/amplicon bindings; do not build a second general read mapper | Conserved integer counts, disjoint categories, denominators, mixed/partial/null cases and control separation. No clone or specificity inference |
| S6: learner/native acceptance | Activate the challenge below only for capabilities actually present; same request/report in GUI/CLI/MCP/scripting/agent examples | Offline replay/parity, separate Glen native checks and human scientific/laboratory review. No current-release blocker |

Schema/shell names are selected only when a slice is implemented. Proposed
verbs such as `editing intent`, `editing base-feasibility`,
`editing prime-feasibility` and `editing import-outcomes` are **not executable**.
Do not register them in the glossary/capability catalog or teach an agent to
emit them now. Initial generic-op reachability precedes dedicated UI polish.

## Deterministic Test Matrix (Planned, Not Run)

New fixture files belong with a nearby origin/recreation/use README: explicit
synthetic sequences/profiles, no misleading gene names, hand-written expected
results independent of production transforms. Public examples need exact source
accessions/releases/licenses/raw hashes. Byte-bound fixtures and any generated
freshness contract also need scoped checkout policy plus LF/CRLF coverage.

| Boundary | Independent expectation / negative control |
| --- | --- |
| Literal edit reconstruction | On synthetic `ACGTAC`: `[1,2)` C-to-T gives `ATGTAC`; insert `GG` at `[3,3)` gives `ACGGGTAC`; delete `[2,4)` `GT` gives `ACAC`. Do not derive expected strings with the engine helper |
| Strand round-trip | Hand-written reverse-complement mates and a nonzero genomic window origin; map first/last base, both interval ends and insertion boundaries 0/n. Same allele after projecting back, no second reversal or use of absolute coordinates as local indices |
| Reference/haplotype | Wrong deleted base, stale semantic/raw digest, wrong accession version/assembly, unresolved phase and ambiguous target all refuse; matching a base alone cannot fix a bad anchor |
| Conversion / window | Both exposed strands for CBE and ABE, every edge of a synthetic declared window, outside-window base, incompatible transversion and absent/wrong-side PAM. No candidate is not a failed request when a complete valid search finds none |
| Bystanders / annotations | Zero/one/several editable sites; intended base plus editable neighbor retained, not assigned probability. Two transcript consequences stay separate; missing CDS phase/splice evidence stays unknown |
| Prime orientation | Both nicked strands with asymmetric components and independent PBS pairing/RTT reconstruction. Incorrect swapped PBS/RTT, reversed-only extension, wrong nick strand, inadequate homology and generic nuclease offset refuse |
| Indels / extra changes | Intended insertion/deletion across nick/homology boundaries; wrong post-edit mapping, repeated-sequence placement and extra PAM-disrupting substitution cannot silently change intent |
| Strategy / components | PE2 with a second nick guide, PE3 without one, and an unjustified PE3b claim refuse; stabilizer/protein architecture cannot silently change strategy |
| Bounds / cancellation | Explicit sequence/profile/search/candidate/byte caps, cancel before/during work, empty complete result versus truncated/budget/cancelled result. No partial best/absence proof; deterministic ties and counters |
| Local tools | Pin binary/container/version and effective settings; malformed/unknown versions, unsupported chemistry and incompatible input hashes refuse. Cross-check normalized tiny outputs, not ranking agreement alone |
| Counts / denominator | Synthetic counts 20/10/5/5/45/5 sum to 90 retained reads for the six categories; record 100 input/10 QC-excluded independently. Negative/fractional/overflow counts, duplicate alleles and inconsistent totals refuse; zero denominator is unavailable |
| Bulk / controls | Same bulk fraction can have different cell genotypes; never infer zygosity. Missing/unequal control, low-coverage and ambiguous reads do not become untreated-zero or independent replicate evidence |
| State / adapters | Same request through shared op/shell/GUI-host/MCP/scripting; read-only leaves project/revisions/undo unchanged. Approval denial, tampered/stale report, cancellation and project replacement never materialize/adopt results |
| Persistence / exports | Save/reopen retains profiles, sequence and source hashes, typed unknowns and component order; rendering reads saved results without recomputation. LF/CRLF inputs, quoted/Windows paths and atomic no-clobber exports use existing portability rules |

## Inactive Learner Challenge

Working title: **Can this exact allele be proposed, built and observed?**
The [standalone learner draft](tutorial_drafts/base_prime_editing.md) expands
these checkpoints with hand-written answers and tested synthetic teaching data.
It stays outside the active tutorial tree/source/catalog/manifest registries;
Help's missing-catalog fallback recursively scans that tree. No editing command
sequence, screenshot or published protocol is promised yet.

1. Compare the synthetic reference/desired strings above. Explain half-open
   bases versus insertion boundaries; intentionally reject a stale reference.
2. With an explicitly synthetic CBE/ABE profile, inspect both strands, PAM,
   window and all bystanders. Explain why feasible is not effective or safe.
3. Inspect a synthetic PE2 component map, reject a reversed extension, then
   inspect the reconstructed substitution and short-indel alleles. Explain
   what PE3/PE3b would add without inventing an extra guide or sequence change.
4. Review a component handoff without ordering or mutating the source. Only
   when separately implemented, demonstrate exact approval and derived lineage.
5. Interpret the synthetic 90-read result and its exclusions/untreated control.
   Compare intended-only with intended-plus-bystander; discuss unknowns and
   why neither fraction determines a clone's genotype or off-target safety.
6. Later substitute a separately selected public nuclear-DNA example with
   pinned reference, editor, component and outcome provenance. A paper's gene
   label is not a GENtle-validated target; require human selection/review before
   adding a real locus. Keep wet-lab conditions outside this software lesson.

Activate only after relevant slices, independent expected results, canonical
interface/checkout checks and explicit tutorial review metadata exist. Codex
readability review, Glen's native evidence and human scientific approval are
different records; absent review dates must remain absent.

## Ownership And Next Decision

- Automated tests: exact sequence/coordinate/chemistry geometry, bounded search,
  counts, read-only/approval/stale state, persistence and adapter/checkout parity.
- Local external comparison: version-pinned PrimeDesign and CRISPResso2 formats,
  settings and original outputs. Differences need explanation, not forced hashes
  or changed intent; no automatic upload or model call.
- Glen: native both-strand maps, component/warning readability, selection,
  cancellation/project-switch rejection, save/reopen and raw screenshot/report
  evidence at a consolidated candidate. No request to repeat it now.
- Human scientific/laboratory reviewers: editor/cell/delivery applicability,
  public locus/haplotype choice, controls, intervention approval, experimental
  execution, clone genotype, specificity/function and publication decisions.

Recommended next authorization: S1 only, with synthetic references and no live
tools. In parallel, obtain lawful full L1 access and review exact conventional
profiles for S2. No need to choose an experimental delivery recipe or efficacy
model to implement the reference validator. This preparation does not authorize
either slice, change the `.13` aim or waive any existing acceptance gate.

Preparation checks: tracked and new-file whitespace checks pass; session-close
reports four OK, two expected warnings (intentional documentation changes and
manual plan fidelity), no failures. Roadmap remains 400 lines. No Rust/runtime,
external-tool, native GUI or laboratory test result is claimed by these checks.
