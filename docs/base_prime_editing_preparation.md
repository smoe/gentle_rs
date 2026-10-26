# Base And Prime Editing: Preparation Brief

Status: S1 read-only SNV-validation implementation, 2026-10-23; later slices
remain proposals. No editor/guide design, executable laboratory protocol,
approved profile, active editing tutorial or new release gate is introduced.

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

The initial source/literature review made no Rust build, model/provider invocation,
sequence upload, order, tag or publication. Browsing public references is not
contacting a service with a design. No external design/analysis tool was run.
No direct Claude invocation was made. The initial 2026-10-20 preparation was
not Claude-reviewed; the owner-forwarded 2026-10-22 software-contract review
and Codex reconciliation below have a narrower scope than scientific approval.

### Base-Protocol Follow-Up Identity

The 2026-10-22 follow-up inspected `gentle_rs_2_main` at
`2ea67b1bc35c9a3a662546fa871aad675b8cd169`, matching local `main`.
No fetch was performed. The [base-protocol workstream below](#base-editing-protocol-workstream)
extends this brief; it does not replace the shared S1/S2 reference/feasibility
contracts, activate the learner draft, or make prime editing a dependency of
base editing. The owner requested a plan, not implementation, then forwarded
Claude's read-only software-contract reply to the
[review prompt](base_editing_protocol_claude_review_prompt_20261022.md).
The review and revised minimal contract are reconciled below. Claude did not
independently access the literature or run tests; no human applicability or
scientific approval is claimed. Existing untracked `output/` bundles remain
untouched. The owner subsequently authorized revised S1 implementation and
focused local Rust tests/checks; this does not authorize S2 or laboratory work.

The implementation was rebased onto local main
`2351fdcfe6534dacd3980228e4ac918c74dadbaf` on 2026-10-23, retaining its
Nanopore operation/result additions. No fetch or push was performed. Its new
RNA-summary tests called nonexistent `journal()` twice; the bounded prerequisite
repair uses existing `journal_len()` without changing assertions or RNA behavior.

### Implemented S1 Boundary

`EvaluateGenomeEditIntent` now composes the existing reference/context records
through `gentle_protocol::genome_editing`, with one shared validator in
`src/engine/analysis/genome_editing.rs`. The [protocol contract](protocol.md#read-only-snv-intent)
and [tested CLI example](cli.md#read-only-snv-intent) describe the implemented
request/report rather than a proposed command. Reachability is generic `op`
across Shell/CLI/MCP/scripting; no dedicated editing verb or GUI designer.

The bound window is 1..100000 linear A/C/G/T bases, with explicit case-only
normalization and independent source/genomic-forward/desired digests. One SNV
is checked on either orientation, including transversions as sequence intents,
not editor feasibility. Required identities/digest/anchor/haplotype, optional
paired annotation labels, ambiguity, bounds and stale state fail closed.
Inline input invents no stored sequence/projection identity. Loaded genomic
input reuses checked ROI projection, but additionally requires typed retained
provenance identifying the same versioned assembly/contig. The newest typed
record must itself be complete and matching: no older-record fallback or
mixing records through caller-supplied aliases. Unversioned catalog names or
aliases alone are deliberately insufficient; no filesystem reference
resolution or network fallback is introduced. Caller-declared inline reference
association remains labelled as such.
The reused local projection retains its raw stored-byte digest; it is checked
separately from the request/report's uppercase semantic hash, including for
lowercase imports. No legacy projection or source bytes are rewritten.

The returned `OpResult.genome_edit_intent` is transient; no project sequence,
report store, file, undo entry or mutation/structural revision is created.
The operation also skips generic lineage/container repair, preserving even
unreconciled project state on success and refusal rather than weakening the
state-preservation tests.
Generic MCP retains its `confirm=true` tool gate but does not create/rewrite
state for this operation; the CLI's legacy direct `op` path also withholds its
usual state-file save. Existing `OpResult` JSON can omit the new optional
field; other operations' serialization and legacy ROI deserialization remain
unchanged. Scientific/editor suitability and all later workflow slices remain
unimplemented, not implied by arithmetic or parser checks.

### Landed S1 Review Reconciliation (2026-10-25)

Source check: local main `09e62b3f`, no fetch, build or scientific review.
The owner forwarded Claude's source-only review of `472c84d6` at `b1bae75e`.
Original S1 plan: validate one reference-bound SNV without project mutation;
keep S2 and broader lifecycle changes outside that authorization.

| Claude feedback | Source check and bounded follow-up |
| --- | --- |
| A: no test pins the no-reconcile guard | Existing `loaded()` fixtures insert DNA directly into a new engine without lineage/containers, and `evaluate_unchanged()` compares state bytes on success/refusal. They already cover the guard indirectly. Add a named regression with independent controls proving each reconcile routine changes the fixture, for both strands and inline/loaded targets. |
| B: derive one inspection helper from registry `ReadOnly` | Do not apply this literally: the registry also labels `SaveFile` read-only, and currently labels `PlanDnaSequenceDesign` mutating despite its history/confirmation exemptions. Its optional path can write a report. One boolean cannot replace these distinct lifecycle policies safely. A shared, explicitly reviewed policy remains a prerequisite before S2, not part of this narrow patch. |
| C: clarify whether S1's guarantee applies to every preview | Scope no-reconcile/no-state-file-save explicitly to S1 in the decision record; do not silently change legacy preview or file-writer behavior. |
| D: mixed coordinate names look like an off-by-one | Name the numerically equivalent half-open window end at the comparison. Preserve arithmetic, checked bounds and existing edge tests. |
| E: submitted versus canonical allele wording | State that S1 uppercases returned intent alleles without separate raw-allele fields. The original inline target spelling remains in `request.target`. Preserve submitted/canonical placement only as a requirement for future repeat-indel normalization, not an implemented S1 promise. |

The named regression and case-normalization controls now have
[headless macOS execution at `98c9331a`](operation_lifecycle_policy_proposal.md#executed-local-verification-2026-10-26):
14 S1 tests pass, including MCP state-file preservation. This is not native GUI,
Windows/full-suite or scientific acceptance, or authority to begin S2. No
serialization or biology contract is broadened.
The owner requested a [request-aware lifecycle proposal](operation_lifecycle_policy_proposal.md)
and [forwardable Claude review prompt](operation_lifecycle_claude_review_prompt_20261025.md)
on 2026-10-25. Claude's source-only review is now reconciled in that proposal.
The owner authorized only optional report-write confirmation and conservative
static metadata; their focused approval/parity execution is recorded there.
The S1 decision is narrowed to existing apply-time and CLI/MCP disk behavior:
generic constructors can still repair JS/Lua state-out. Construction/workflow/
autosave coverage and destination-alias/writer rules remain separately gated,
as does S2. No full lifecycle or cross-adapter preservation guarantee is adopted.

### Software Review Reconciliation (2026-10-22)

Original proposal: S1 reference-bound intent, S2 editor feasibility, S4 saved
protocol/inspection, S5 measured results, then S6 learner activation; S0 profile
review is parallel and prime-editing S3 is not a base-editing dependency. The
original proposal remains in the review prompt rather than a second pending
work list. Codex rechecked the following source at unchanged `2ea67b1b`.

| Claude feedback | Source check and revised Codex decision |
| --- | --- |
| A defaulted `GenomicRegionLocalProjection` accepts `{}` with `Current` | The serde defaults are permissive, but this is not proof of a live validation bypass. `src/engine/state/genomic_regions.rs::validate_local_projection` already validates digest, explicit anchor strand and geometry; regulatory-fragment consumers also refresh status. S1 must validate required request inputs, construct projections internally and never accept a caller's status as proof. Preserve legacy deserialization. |
| Drift statuses need explicit refusal | Accepted. The S1 status-to-error table below forbids warning-only success on stale/missing projection evidence. Use existing error classes and explicit refusal reasons. |
| 1-based anchors need an explicit half-open conversion | Accepted. Add exact anchor/window formulas, checked arithmetic and local-index conversion tests; no zero anchor or saturation-to-zero rescue. |
| Reuse reference/context/projection rather than a parallel flat binding | Accepted as composition through a strict editing-specific wrapper. Existing optional fields stay backward-compatible. `GenomicRegionReference` supplies contig identity; `BiologicalContext` supplies annotation identity. Validate and reconcile both, including numeric/string taxon IDs, before `semantically_matches`; that method compares context fields, not the whole approval basis. |
| Reject guides whose nuclease was defaulted to SpCas9 | Accept the no-inference requirement, not an impossible provenance test. Normalization discards whether SpCas9 was explicit. S2 requires a separate exact editor profile and validates actual geometry; a legacy guide alone is insufficient, but a matching explicitly selected profile need not be rejected. |
| Optional practical filtering lets ambiguous bases through | Accepted. Exact S1/S2 inputs need their own A/C/G/T validation, independently of practical filters. Distinguish observed PAM letters from an intentionally degenerate profile PAM pattern. |
| Reuse DNA-design approval verbatim | Reuse its content-bound-preview/fresh-state-check pattern, not the private DNA-design report hash as a general approval API. It does not currently bind an editing profile, cell/SOP or laboratory reviewer. Context equality alone cannot authorize a changed candidate or SOP. |
| Require annotation release and real accessions for every S1 request | Narrowed. Sequence-only validation does not require annotation evidence. If bound, annotation source/release must be complete and consistent; otherwise consequences remain unavailable. Synthetic references use a disjoint namespace and do not invent biological accessions or taxon IDs. |

Claude confirmed the literal edit/strand arithmetic, count example and absence
of an editing route by source inspection. These are not newly executed Rust
checks or independent chemistry verification. The existing inactive-tutorial
guard already rejects teaching nonexistent commands; do not add a duplicate.
The owner accepted S1 only, narrowed to one SNV; implementation follows that
boundary, not the entire proposed workstream.

## Verified Gap Map

These are preparation-baseline source findings at the identities above, not a
current absence claim for the implemented S1 boundary. Remaining editor/guide,
protocol and measured-outcome gaps still apply.

| Need | Existing capability and evidence | Reusable composition | Missing primitive / deferred work |
| --- | --- | --- | --- |
| Reference identity | `GenomeExtractionProvenance`, `BiologicalContext`, `GenomicRegionReference` and engine ROI validators/projection helpers in [genomic_regions.rs](../src/engine/state/genomic_regions.rs) | Compose existing metadata with a strict editing-request wrapper; reuse identity/digest/geometry validation without trusting a serialized status | Original-to-desired edit validation, stronger required-field policy, explicit synthetic/genomic distinction and stale-reference refusal; current ROI types alone are not an edit validator |
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

### L6: Owner-Requested Base-Editing Review

[Porto, Komor, Slaymaker and Yeo](https://www.nature.com/articles/s41573-020-0084-6),
DOI `10.1038/s41573-020-0084-6`, published 2020-10-19, Nature Reviews Drug
Discovery 19, 839-859. Read the [full author manuscript in PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7721651/)
on 2026-10-22; the publisher landing page supplies only a preview. Consulted
DNA/RNA classes, specificity, product purity, targeting/bystanders, Box 1,
delivery, applications and limitations, including Figures 1-4 and Tables 1-2.
No original bytes or figures were retained or licensed for redistribution.

The review motivates explicit editor/guide pairing, context-dependent windows,
bystander inspection and distinct guide-dependent DNA, guide-independent DNA
and RNA evidence. Its delivery discussion motivates a separately reviewed
experimental context. It is a historical review, not a current construct
catalog, validated predictor or laboratory SOP. Do not copy its efficacy
percentages, clinical expectations or broad technology-limit statements into
GENtle defaults. Recheck primary construct evidence and corrections before
profile activation. Diagram inspection does not approve coordinates.

### L7: Accessible Practical Methods Alternative

[Osgood et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11948366/),
*Genome editing with programmable base editors in human cells*, DOI
`10.1016/bs.mie.2025.01.001`, Methods in Enzymology 712, 351-404 (2025);
PMC lists manuscript availability as 2026-02-05. Full-text methods were
accessible on 2026-10-22. Sections 3-4 distinguish editor choice from guide
design, with subsequent component, delivery and analysis methods. Its window
discussion also cautions that activity can occur outside a typical window.
Therefore a displayed activity window must not imply zero activity elsewhere.

This is an alternative source for preparing a context-specific protocol,
not an approved GENtle preset. If chosen instead of L1, review its exact
constructs, supplements, corrections and cell/delivery applicability on their
own merits; do not borrow L1's review status. A permitted full-text source
exists, so obtaining L1 specifically need not block the software workstream.
No full-text redistribution, method parameters or editor sequences are added
here, and no complete contemporary literature/correction search is claimed.

## Protocol-Applicability Matrix

These are proposed applicability states, not current runtime enum variants.
`reference_only`, `needs_review`, `unsupported` and `reviewed_for_context` must
stay distinct. A citation alone cannot produce the last state.

| Evidence / choice | Permitted preparation use | Missing applicability binding | Initial disposition |
| --- | --- | --- | --- |
| L1 plus correction, or independently reviewed L7 | Conventional CBE/ABE feasibility and protocol questions | Exact source/construct, PAM/spacer/window/context and cell/delivery review; full L1 access if that method is selected | `reference_only`; no activated window or SOP |
| Explicit conventional CBE profile | C-to-T on the declared exposed strand | Versioned construct/source and position-specific activity/context evidence | `needs_review`; no transversion/editor substitution |
| Explicit conventional ABE profile | A-to-G on the declared exposed strand | Same bindings, separately reviewed from CBE | `needs_review`; not a CBE preset renamed |
| L2 / conventional prime editor | PBS/RTT/nick reconstruction and strategy separation | Exact editor/scaffold/geometry and bounded design policy | `needs_review`; PE2 first, PE3/PE3b later |
| L3 plus corrected Scheme 2 | Cell-context/correction teaching and later cloning review | hPSC, delivery and intervention-specific human approval | `reference_only`; not a generic workflow |
| PrimeDesign / CRISPResso2 | Local interchange/comparison and saved outcome interpretation | Pinned versions, formats, parameters, inputs and complete artifact hashes | `needs_review`; no tool execution in preparation |
| Other organisms, organelles, nonconventional chemistries or therapeutic use | Record the request and missing evidence | Separate applicability brief and approval | `unsupported` in initial implementation |

## Proposed Typed Contract

S1 request/report types and `EvaluateGenomeEditIntent` are implemented in the
module named above; later editor, search, protocol and outcome names remain
**proposals only**, not callable commands. Use one shared engine analysis path
and existing `Operation`/shell adapters.
Reuse context/region primitives above without treating their optional metadata
as a completed editing-reference check. Do not extend `GuideCandidate` into an
unvalidated bag of optional editor fields or introduce another adapter registry.

| Proposed type / operation | Minimum typed content and invariant |
| --- | --- |
| Validated editing binding | A small strict request/engine wrapper over `GenomicRegionReference`, optional explicitly bound `BiologicalContext`, exact loaded/inline operand, required expected sequence digest, anchor and selected design haplotype. Genomic identity fields are validated as required, not made nonoptional globally. Construct a validated projection internally where the stored-sequence contract applies; never take one as caller-certified input. Synthetic identity is disjoint and does not invent accessions. No duplicate flat reference registry |
| `GenomeEditIntent` | Validated binding, one genomic-forward `start_0based..end_0based_exclusive`, exact `deleted`/`inserted` A/C/G/T strings and reconstructed desired digest. First S1 accepts one changed base only; later explicit short-indel support preserves the same coordinate convention. Reject unchanged/mismatched/ambiguous targets |
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

### Minimal S1 Software Contract After Review

The implemented `EvaluateGenomeEditIntent` operation takes a strict, versioned
request with exactly one loaded-sequence or inline operand, an expected
`sha256:` digest, explicit reference kind/identity and selected design haplotype,
an explicit plus/minus anchor, and one genomic-forward SNV interval plus
single-base `deleted`/`inserted` strings. No input projection/status is accepted.
No serde whole-record defaults may turn absent required fields into a binding;
validate nested legacy records explicitly as well. Reject empty, malformed or
nonmatching digests rather than merely checking that a string exists.

For genomic mode, require explicit versioned assembly/contig accessions and
taxon identity in the existing reference record. Reconcile supplied context
assembly/taxon against that record, checking numeric taxon conversion. A loaded
genomic source must match its retained anchor/provenance, not a caller override.
For inline genomic data, report the external reference association as caller-
supplied unless independently retained evidence verifies it; internal sequence
consistency is not authentication against a reference provider. Synthetic mode
requires an explicit synthetic namespace, origin and sequence identity; it must
not fabricate real accessions or silently reclassify a rejected genomic input.
Choosing a reference haplotype is not a cell-line genotype claim.

Annotation metadata is not required for arithmetic S1. If annotation evidence
is explicitly bound, require source/release together and reject conflicting or
half-specified bindings. With neither, return `annotation_not_supplied` and no
CDS/splicing conclusion. Preserve the existing context/reference types; do not
make their optional fields mandatory for unrelated consumers.

Use canonical A/C/G/T DNA with an explicitly documented case policy for the
bound window and alleles; reject ambiguity rather than cleaning it away. Bind
the loaded-or-inline sequence and its genomic-forward representation with
separate digest roles. Stored-input projections may reuse
`GenomicRegionLocalProjection` after engine checks; inline operands must not
invent a stored `seq_id` merely to populate that type. Reuse checked geometry
helpers where their preconditions fit, not the nonempty-ROI restriction for a
future insertion boundary. Require `end - start == 1`, one base in each allele,
an exact reference match and `deleted != inserted` in the first delivered S1.

| Existing projection status / failure | Required S1 handling |
| --- | --- |
| `Current` | Necessary but insufficient: recompute/validate digest, anchor/reference, bounds and geometry before constructing a successful report |
| `SequenceDigestMismatch` | `InvalidInput` with an explicit digest-mismatch reason; no successful intent report |
| `AnchorMismatch` | `InvalidInput` with an explicit anchor/reference-mismatch reason |
| `OutsideSequence` | `InvalidInput` with an explicit outside-window reason |
| `SequenceUnavailable` | `InvalidInput` if discovered while validating a bound projection; a requested loaded sequence absent at initial lookup remains the ordinary `NotFound` case |
| Malformed/omitted input, context conflict, wrong alphabet, invalid anchor or arithmetic/index overflow | `InvalidInput`; do not repair by defaults, saturation, context fallback or warnings |
| Well-formed insertion/deletion, multi-edit, transcript/HGVS or circular target | `Unsupported` in initial SNV-only S1, not a claim of biological infeasibility |

Return a `GenomeEditIntentReport` containing the checked binding, canonical
intent and reconstructed sequence/digests, with explicit unknown annotation
and reference-verification limits. No editor recommendation, persistence,
sequence creation, file writes, network, undo entry or mutation revision.
Follow `plan_dna_sequence_design`'s read-only function shape, not its output-ID,
synonymous-design policy, private approval routine or solver. Wire the operation
through the existing generic-op/shell discovery, readiness and effect-audit
paths; dedicated CLI verbs and GUI design work are not S1 requirements.

Planned first-slice tests are required-field/unknown-field parsing (including
`{}` refusal and no caller status), wrong/empty digest, genomic/synthetic and
loaded/inline separation, conflicting/missing identities, annotation-absent
success versus half-bound annotation refusal, ambiguous bases, nonzero-anchor
first/last-base transforms, explicit plus/minus reference reconstruction and
all error classes above. Hand-write expected strings/coordinates independently
of the production helpers. Assert unchanged project bytes, mutation/structural
revision and undo/redo on success and refusal; a generic operation receipt is
not automatically a project mutation. Indel arithmetic remains future teaching/
contract evidence, not an implemented S1 promise. This was the reviewed test
specification; executed implementation checks are recorded separately below.

### Coordinates And Orientations

Canonical edits use genomic-forward half-open intervals. Substitution/deletion
requires `end-start == deleted.len()` and an exact reference match; an insertion
uses `start == end`, empty deleted string and an explicitly anchored interbase
position. Resolve the loaded extraction to a genomic-forward reference window
`R` and map absolute coordinates to window-local `s,e` before reconstructing
`R[..s] + inserted + R[e..]`. Never index a reverse-oriented local sequence with
absolute genomic coordinates. Bind both the original extraction and `R` hashes;
the report retains the transforms back to the loaded display. S1 returns uppercase
intent alleles without separate raw-allele fields; inline target spelling stays
as submitted. A future repeat-indel normalization must preserve submitted and
canonical placements; never shift coordinates without a recorded reference-bound
transform. Synthetic examples below use
an explicitly zero-origin forward reference. Initial S1 accepts one SNV only;
the insertion/deletion rules describe later extensions, not first-slice scope.
Short indels and overlapping/multi-edit/haplotype combinations require later
authorization. Selecting a reference haplotype is a design basis,
not a claim about the genotype of the actual cell line.

Existing projections use **1-based inclusive** `u64` anchors. For
`A1 = anchor_start_1based` and `E1 = anchor_end_1based`, require `1 <= A1 <= E1`,
then compute `A0 = A1 - 1`, `E0 = E1` (exclusive) and
`n = E1 - A1 + 1` with checked arithmetic. Require `n` to equal the bound
sequence length. For absolute genomic `[g0,g1)`, local plus coordinates are
`[g0-A0,g1-A0)`; local minus coordinates are `[E0-g1,E0-g0)`.
Keep genomic arithmetic in `u64`; convert only validated local indices with
`usize::try_from` before indexing. Refuse zero/reversed anchors, unstranded
anchors, length mismatch, out-of-window intervals and overflow. Do not use
saturating subtraction to legitimize a missing anchor.

Independent example: anchor `1001..1006` inclusive has `A0=1000`, `E0=1006`
and length 6. Genomic `[1000,1001)` maps to local `[0,1)` on plus and `[5,6)`
on minus; `[1005,1006)` maps to `[5,6)` and `[0,1)`. Future insertion boundaries
1000 and 1006 map to 0/6 on plus and 6/0 on minus, without fabricating a base.

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
S2 requires an explicitly selected exact `BaseEditorProfile`; a guide's
`nuclease` string (including the legacy SpCas9 default) is never an editor or
chemistry selection. Test that a guide alone refuses, and separately that the
same guide with a matching explicit profile still needs sequence/PAM/window
validation. Non-ACGT observed reference/spacer/PAM bases refuse independently
of `GuidePracticalFilterConfig`; an IUPAC profile PAM pattern such as NGG is
deliberate pattern syntax, not ambiguous observed DNA.

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
or add per-branch release/Glen campaigns. S1 is the implemented boundary,
not authorization for a combined editing workbench or the remaining slices.

| Slice | Smallest deliverable and dependencies | Exit evidence / limitations |
| --- | --- | --- |
| S0: evidence/profile preparation | Review permitted full L1 text plus correction, or an independently selected full-text method such as L7; select one exact conventional CBE and one ABE construct. Bind PAM/window/context provenance without guessing numbers; select prime architecture separately | Human applicability review and raw-byte source ledger. Until then use explicitly synthetic test profiles only, never biological presets |
| S1: reference-bound intent | Implemented strict wrapper over existing reference/context primitives; one explicit SNV, checked mapping and transient report via generic op | Authored both-orientation/anchor/omitted-field/refusal/unknown-annotation, history and MCP file-preservation checks; see executed verification separately. No editor, GUI or indel promise |
| S2: conventional base feasibility | S1 plus explicitly selected reviewed CBE/ABE profile (or clearly synthetic test profile); strand-correct candidate/PAM/window checks and bystander/annotation report | Exhaustive tiny cases, absent PAM/impossible conversion, annotation-unknown behavior and bounded-search status. No efficacy/off-target ranking |
| S3a: prime component validation | S1 plus one explicit profile and caller-supplied spacer/scaffold/PBS/RTT; read-only PE2 reconstruction before any enumeration | Independent strand/PBS/RTT/nick/indel reconstruction; fail closed on bad orientation, added edits or insufficient flanks |
| S3b: bounded prime design | S3a plus explicit finite search ranges, then separately PE3/PE3b nick-guide rules; local PrimeDesign interchange/comparison | Retained normalized candidates and mismatch explanations at pinned settings. External output never bypasses S3a; tool installation/execution requires later authorization |
| S4: inspectable component/protocol handoff | Reuse oligo/cloning/sample contracts and document backends; base-specific steps follow the workstream below without waiting for S3. Maps/GenBank/EMBL/JSON read saved reports; no scoring in rendering | Targeting bases, added bases, overhangs and RNA roles visibly distinct; protocol draft and human SOP review remain separate. Approved materialization rechecks exact digests, lineage, undo and stale-job rejection |
| S5: local observed-result import | Import version-pinned CRISPResso2/explicit allele-count exports against S1/reference/amplicon bindings; do not build a second general read mapper | Conserved integer counts, disjoint categories, denominators, mixed/partial/null cases and control separation. No clone or specificity inference |
| S6: learner/native acceptance | Activate the challenge below only for capabilities actually present; same request/report in GUI/CLI/MCP/scripting/agent examples | Offline replay/parity, separate Glen native checks and human scientific/laboratory review. No current-release blocker |

Schema/shell names are selected only when a slice is implemented. Proposed
verbs such as `editing intent`, `editing base-feasibility`,
`editing prime-feasibility` and `editing import-outcomes` are **not executable**.
Do not register them in the glossary/capability catalog or teach an agent to
emit them now. `op` with `EvaluateGenomeEditIntent` is the implemented S1 route;
initial generic-op reachability precedes dedicated UI polish.

## Base Editing Protocol Workstream

This is the base-specific implementation order within S0/S1/S2/S4/S5/S6,
not a second editing roadmap. Start with conventional nuclear-DNA CBE/ABE
research in cultured mammalian cells, one explicit SNV and one selected
reference haplotype. Keep RNA/mitochondrial editing, transversion/dual editors,
prime editing, multiplex screens, clinical and germline use outside this
initial applicability contract. Unsupported here does not mean biologically
impossible. The 2020 review is a design rationale, not a present-day limit.

### What The User Should Be Able To Do

From a selected DNA base or explicit allele, state the desired change; compare
reviewed editor/guide options; inspect windows and potential bystanders; save a
protocol draft with controls and unresolved questions; approve exact components
separately; and later attach measured outcomes. The same saved records must
explain each step in the DNA viewer, CLI, inner agent and exported report.

GENtle should distinguish four questions: **Can this sequence be targeted?**,
**Is this protocol reviewed for this experimental context?**, **What was
actually observed?**, and **What functional claim is supported?** One successful
check must never answer the other three by implication.

### Existing Seams And Missing Boundary

`GuideCandidate` and its normalizer in `src/engine.rs` validate intake fields,
not a reference-matching editor/guide pair. Existing practical filters and
oligo operations in `src/engine/ops/candidate_guides.rs` are reusable after that
validation. `ExportGuideProtocolText` currently emits an oligo-cloning checklist,
not a base-editing experiment. Do not relabel it or reuse its nuclease cut offset
as an editing window. Its storage/export effects also make it unsuitable as
the read-only feasibility evaluator.

`ExportLabAssistantInstructions` currently derives its v2 handoff from the
operation journal and sequence/container state. Reuse the material/checkpoint
presentation and document writer in `src/engine/state/sequence_ops.rs` and
`src/engine/state/lab_assistant_export.rs` where appropriate, but do not silently make
that journal-based operation accept or certify a new experimental protocol.
Only the editing-specific composition is missing; a general laboratory
automation framework or another routine registry is not required.

### Ordered Deliverables

1. **S1: reference-bound intent first.** Implement the shared validator already
   specified above, admitting one SNV only with synthetic reference tests.
   Compose existing metadata through the strict wrapper, not a flat duplicate.
   Bind original/desired alleles, source and available annotation versions, orientation and
   sequence digests. Refuse reference mismatch, stale state and unresolved
   locus/haplotype selection. No editor or guide is selected automatically.
2. **S0/S2: editor profiles and feasibility.** Independently review one exact
   CBE and ABE profile; use explicitly synthetic profiles until then. Enumerate
   and verify spacer/PAM/window relationships on both strands. Distinguish
   hard incompatibility, context preference and unavailable evidence. A finite
   search reports completeness and rejection reasons, not an efficacy ranking.
   Potential bystanders retain per-transcript consequences, including
   annotation-backed CDS/splice effects; unavailable annotations remain unknown.
   Show supported window evidence and its limitations, not an assertion that
   edits outside it cannot occur. Do not multiply independent per-base
   probabilities to invent an allele-outcome distribution.
3. **S4a: read-only protocol draft.** Compose a proposed versioned
   `BaseEditingProtocolDraft` from one selected feasibility candidate and an
   explicitly selected context/SOP binding. Minimum fields are reference,
   intent, candidate, profile and source digests; objective; cell/sample and
   delivery context; exact component roles; ordered stages; control/sample
   roles; assay/readout plan; unresolved fields; and scoped review records.
   Proposed stages are reference verification, component preparation,
   delivery under the selected SOP, sample collection, sequence analysis and
   interpretation. Each stage cites its inputs, outputs, dependencies and
   checkpoint. These are software workflow stages, not invented bench settings.
   Values with units, provenance and reviewer identity may come from an
   approved local SOP; missing quantities, timing, kit/vector identities or
   acceptance thresholds stay absent. A partial draft remains inspectable
   with explicit blocking reasons. Software consistency, editor applicability,
   laboratory approval and execution observations are separate states; changing
   a bound component, SOP or input invalidates the corresponding approval.
   Reuse existing content-bound preview/stale-check mechanisms, with an editing-
   specific digest basis; do not interpret the private DNA-design hash or
   `BiologicalContext::semantically_matches` as laboratory approval. Require
   annotation source/release only when that evidence is actually bound.
4. **S4b: inspect and hand off.** Add a DNA-viewer action with shared contextual
   availability, not a detached biological implementation. Draw reference and
   proposed allele, spacer/PAM, exposed strand, evidence-backed window, intended
   change, potential bystanders and per-transcript annotations. Use distinct
   glyphs/labels as well as colors, exact coordinates and readable side summaries
   plus hovers. Candidate selection only selects a report; it never edits the
   source DNA. Reuse background jobs, cancellation and stale-result rejection.
   Export saved draft/results to JSON and human-readable Markdown/PDF with the
   same warnings and provenance, without rerunning design. Reuse document/
   vector backends, keeping format parity separate from laboratory approval.
   Explicit oligo/construct creation is a later confirmed action, using exact
   reviewed vectors/scaffolds and separate targeting/added bases; no default
   lenti template, automatic ordering, delivery or execution.
5. **S5: results close the loop.** Bind local allele-count/CRISPResso2 imports
   to the saved intent, amplicon, sample and control identities. Preserve raw
   counts, callable/excluded denominators, intended-only versus bystander and
   other outcomes. Keep untreated, delivery/mock and any reviewer-required
   controls separate; report missing controls instead of manufacturing them.
   Bulk sequence data cannot prove clone genotype, genome-wide specificity or
   phenotype. Guide-dependent DNA, guide-independent DNA and RNA evidence get
   separate assessment slots, initially `not_assessed`; a genome search or
   clean amplicon must not convert the latter two into passed checks.
6. **S6: activate one honest lesson.** Extend the existing inactive draft with
   a synthetic CBE/ABE feasibility-to-protocol walkthrough and identical saved
   requests/reports across interfaces. Then add one separately approved public
   reference/construct example. Keep expected sequence checks, native GUI
   screenshots, observed experimental data and human interpretation distinct.
   Publish the lesson only once the taught capabilities actually exist.

The first useful milestone is **S1 plus synthetic S2 and a saved draft
preview**: the user can inspect why an edit is or is not targetable and see
exactly what laboratory information is still missing. It is not an executable
SOP. Reviewed biological profiles and experimental context are additional gates,
not prerequisites for testing the software with clearly synthetic data.

### Shared Contract And Acceptance

Use the proposed editing records above in a small `gentle-protocol` module,
with one shared engine analysis/orchestration path. Add each implemented
operation to the existing shell, effect classification, capability/readiness
and parity contracts. CLI, GUI, scripting, generic MCP and the inner agent all
invoke that path; no adapter calculates chemistry or fills absent SOP fields.
The agent may explain and propose exact available actions, not grant scientific
approval. Proposed names in this brief remain non-callable until registered.
Dedicated GUI affordances follow the first headless contract, not vice versa.

The deterministic matrix below remains authoritative. Add protocol-specific
guards for missing/SOP-mismatched parameters, wrong units and source versions,
changed candidate invalidating approval, and identical saved content across
JSON/Markdown/PDF. Test that a visually reversed window still describes the same
canonical change, that preference evidence is not treated as impossibility,
and that an empty complete search differs from a truncated search. Save/reopen,
project replacement and replay must preserve unknowns and reject stale binding.
Evaluation/preview writes no files and changes no project or undo history;
explicit report retention, exports and component materialization each declare
their own effects and use existing approval/rollback semantics.

Optional local BE-Designer/beditor/BE-Hive comparisons can be a later evidence
adapter, not a runtime dependency. Pin source/model/version, parameters and
applicability before comparing exact candidate geometry. Model predictions must
remain labelled predictions and separate from measured outcomes; no implicit
sequence upload or network dependency. No predictor is selected or run here.

Open owner decisions are which exact conventional constructs and cell/SOP
context to review, which public example to teach, and who can approve that
context. Their absence does not block S1/synthetic tests. Glen's native and
candidate-bound replay checks belong in one consolidated future handoff, not a
new per-branch audit or current-release gate. Claude's software-contract review
is reconciled above; it did not independently verify the literature or grant
scientific, native or implementation acceptance.

## Deterministic Test Matrix (Planned, Not Run)

New fixture files belong with a nearby origin/recreation/use README: explicit
synthetic sequences/profiles, no misleading gene names, hand-written expected
results independent of production transforms. Public examples need exact source
accessions/releases/licenses/raw hashes. Byte-bound fixtures and any generated
freshness contract also need scoped checkout policy plus LF/CRLF coverage.

| Boundary | Independent expectation / negative control |
| --- | --- |
| Literal edit reconstruction | On synthetic `ACGTAC`: S1 `[1,2)` C-to-T gives `ATGTAC`; future insertion `GG` at `[3,3)` gives `ACGGGTAC`; future deletion `[2,4)` `GT` gives `ACAC`. First S1 refuses the latter two as unsupported. Do not derive expected strings with the engine helper |
| Required fields / permissive legacy records | `{}`, empty/malformed digest, missing anchor, and a supplied projection/status cannot produce a successful S1 report. Leave legacy ROI deserialization unchanged; test the new operation's validation boundary, not only serde behavior |
| Inclusive anchor conversion | `1001..1006` inclusive becomes `[1000,1006)` and length 6. Test genomic first/last bases in plus/minus display; refuse zero/reversed/unstranded anchors, unequal sequence length, checked-arithmetic and local-index overflow |
| Reference kind / context completeness | Synthetic namespace does not require invented accessions/taxon or annotation release. Genomic mode requires its declared reference identity; reconcile numeric/string taxon and assembly values. No annotation binding succeeds with unavailable consequences; conflicting or half-bound claimed annotation refuses |
| Strand round-trip | Hand-written reverse-complement mates and a nonzero genomic window origin; map first/last base, both interval ends and insertion boundaries 0/n. Same allele after projecting back, no second reversal or use of absolute coordinates as local indices |
| Reference/haplotype | Wrong deleted base, stale semantic/raw digest, wrong accession version/assembly, unresolved phase and ambiguous target all refuse; matching a base alone cannot fix a bad anchor |
| Conversion / window | S2 requires an exact separate profile even for a normalized SpCas9 guide. Test both exposed strands, window edges, unsupported outside-window position, incompatible conversion and absent/wrong-side PAM. Unsupported-by-profile is not proof of zero biological activity. No candidate is not a failed request when a complete valid search finds none |
| Bystanders / annotations | Zero/one/several editable sites; intended base plus editable neighbor retained, not assigned probability. Reject N/R/other ambiguity in observed reference/spacer/PAM independently of practical-filter flags; retain IUPAC syntax in profile PAM patterns. Two transcript consequences stay separate; missing CDS phase/splice evidence stays unknown |
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

S1 was separately authorized and implemented. The next possible development
slice is S2 with explicitly synthetic profiles, only after new authorization;
review exact conventional profiles in parallel using lawful full L1 access plus
correction or independently reviewed L7. No experimental recipe or efficacy
model is needed for arithmetic validation. This does not select a new `.13`
aim, authorize laboratory use or waive an existing acceptance gate.

Initial 2026-10-20 preparation checks: tracked and new-file whitespace checks pass; session-close
reports four OK, two expected warnings (intentional documentation changes and
manual plan fidelity), no failures. Roadmap remains 400 lines. No Rust/runtime,
external-tool, native GUI or laboratory test result is claimed by these checks.

2026-10-22 plan follow-up checks: all 20 tests in
`python3 -m unittest scripts.test_tutorial_gap_lessons -q` pass, as does
`git diff --check`. These exercise existing teaching/source contracts, not an
editing implementation. Session-close reports four OK, two expected pending-
change/manual-plan warnings, zero failures; roadmap stays at 400 lines.
No local Rust, new CLI/native run or biological acceptance is claimed.

### S1 Executed Verification (2026-10-23)

Platform: macOS 27.0.1 (26A434), aarch64-apple-darwin, Rust
1.100.0-beta.5 (8ea38001b). Tested local main `2351fdcf` plus this S1 working
change and the two-call RNA test prerequisite repair described above. All Cargo
commands used the repository's default features/profile and configured 16-MiB
test-thread stack, with `--locked --offline -j1`; no dependency changes.

- `cargo test -q --locked --offline -p gentle-protocol -j1`: 113 passed,
  zero failed/ignored; includes both new strict-schema/capability tests.
- `cargo test -q --locked --offline --lib genome_edit_intent -j1`: 13 passed,
  zero failed/ignored, including exact state/history preservation, both strands,
  newest-record/alias refusal, lowercase stored bytes, explicit unknown
  annotations, actual documented request and MCP file-preservation boundaries.
- `cargo test -q --locked --offline --bin gentle_cli genome_edit_intent -j1`:
  one passed, zero failed/ignored. Exercises the legacy CLI's actual request-file
  loader/executor/output wrapper with LF/CRLF and quoted-space paths, both
  progress modes, byte-preserving success/refusal and ordinary-save control.
- `cargo test -q --locked --offline --lib rna_panel_summary -j1`: seven passed,
  zero failed/ignored after the accessor-only upstream test repair.
- Existing reverse-anchor projection and sequence-design MCP confirmation/apply
  regressions each passed separately (filters
  `genomic_region_reverse_anchor_projects_once_and_keeps_local_strand_explicit`
  and `sequence_design_mcp_confirmation_preview_and_exact_apply_parity`).
- `cargo check -q --locked --offline -j1`, `cargo fmt --all -- --check`,
  `git diff HEAD --check` and
  `python3 scripts/generate_capability_maturity.py --check`: passed.
  The current `gentle-protocol` canonical parity renderer reproduced
  `docs/gui_cli_mcp_parity.md` byte-for-byte. The macOS linker warned about its
  oversized unwind table; compilation and tests completed successfully.
- `python3 -m unittest scripts.test_tutorial_checkouts scripts.test_tutorial_gap_lessons -q`:
  95 run, 94 passed, one failed. The existing-main failure is
  `test_canonical_tutorial_sources_have_review_manifest_entries` because
  `targeted_nanopore_fastq_panel` has no review-manifest entry. It was already
  reproduced after rebasing and does not touch S1 sources; leave unrelated
  tutorial metadata for its owner rather than asserting an unperformed review.
- Session-close: four OK, two warnings (intentional changes/untracked owner
  output and manual plan fidelity), zero failures; roadmap remains 400 lines.

An earlier run exposed generic reconciliation mutating an unreconciled project
before inspection; the production exception fixes that boundary, retaining the
strict state-byte assertions. The final test counts above include that fix.
These are focused software checks, not a full workspace/test-suite, native
Windows/GUI, provider/reference authentication, editor or biological acceptance.
No new external evidence request, release claim, tag, push or S2 authorization.
