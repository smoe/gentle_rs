Read-only second opinion. Repository: /Users/u005069/GitHub/gentle_rs
Task: prepare an implementation plan for prokaryotic/eukaryotic transcription termination annotation of plasmids or an unannotated 3' UTR. Only planning is authorized; do not edit files or run code. Use Read/Grep/Glob only. Give concise critique and a revised minimal plan, with concrete source paths. Identify unsupported biological claims, hidden scope and acceptance gaps.

Original Codex plan:
1. Accept whole plasmid, interval or explicit transcript/CDS; require intended host/polymerase, support both strands and circular origin crossings. Missing transcript boundaries mean downstream-of-CDS search region, not definite UTR.
2. Small provenance-documented reference collection of bacterial terminators and eukaryotic cleavage/polyadenylation elements; exact/bounded approximate matches with orientation, coverage, identity and host context; identity does not prove efficiency.
3. Separate predictors: bacterial intrinsic hairpin plus U tract, documented/imported Rho-dependent evidence only, host-specific Pol II polyadenylation signals/context/cleavage windows, separate Pol III oligo-T rules. Archaeal/Pol I/special non-polyadenylated mechanisms unsupported initially.
4. Reuse transcript grouping/strand-aware mapping; distinguish mature/genomic spans, alternative cleavage and transcript-specific UTRs; motifs alone cannot create exact endpoints; downstream flank outside UTR permitted.
5. Versioned portable report: input hash, topology, host/polymerase, transcript identity, coordinates, model/reference versions, score components, limits and warnings. Separate reference match/prediction/imported annotation/experimental support, and no candidates/unsupported/insufficient context. Literal A/T read runs do not establish poly(A) sites.
6. Shared engine schemas/analysis/operations; reuse construct annotation and feature materialization where appropriate; same behavior across GUI/shell/CLI/MCP/JS/Lua/Python; selected undoable writeback preserving annotations/deduplication.
7. First milestone reference recognition+bacterial intrinsic+narrow mammalian Pol II profile. Deterministic strand/circular/ambiguous/truncated/motif-negative/multiple-cleavage/spliced-UTR/persistence/export/parity tests, provenance fixtures and separate scientific calibration. Targeted tests and cargo check for implementation.

Evidence already gathered:
- docs/biological_extension_guide.md requires biological contract first, smallest primitive/composition, shared engine, portable versioned report and thin adapters. It identifies src/engine/protocol.rs, src/engine.rs, src/engine/ops/operation_handlers.rs, src/engine/analysis/feature_expert_ops.rs and transcript projection helpers as seams.
- src/engine.rs around 22310-22475 maps ConstructRole::Utr3Prime to 3'UTR and Terminator to terminator, adds regulatory_class=terminator; existing annotation writeback eligibility has its own generated_annotation/SequenceSpan semantics. Do not assume arbitrary new prediction records fit unchanged.
- src/engine.rs has MaterializeTataBoxFeatures and other feature materializers; inspect for minimal precedent.
- src/genomes.rs around 4031 recognizes terminator as transcription_termination.
- docs/protocol.md around 2767/5846 mentions terminator overlap in promoter evidence, not a standalone predictor.
- docs/decisions.md 'Original-Read Orientation And Literal Terminal Runs' explicitly says literal A/T are not physical tails, polyadenylation sites or transcript completeness.
- docs/roadmap.md .13 development/release scope is already broad; this plan must not alter release gates or imply implementation authorization. Parking-lot link only if needed.
- Relevant docs: docs/architecture.md, docs/decisions.md, docs/roadmap.md, docs/biological_extension_guide.md, docs/protocol.md, docs/cli.md, docs/gui.md, docs/testing.md; local AGENTS.md applies.
- Existing working tree has unrelated modified docs/roadmap.md, CHANGELOG.md, primer_design_specialist_plan.md, glen_next_verification_handoff.md and unrelated untracked files/submodule changes. Preserve them.

Biological grounding from primary work: bacterial intrinsic structural study https://pmc.ncbi.nlm.nih.gov/articles/PMC10091898/ ; Pol II cleavage versus downstream polymerase termination https://pmc.ncbi.nlm.nih.gov/articles/PMC7615648/ . You have no web tool, so avoid pretending independent external verification.

Please especially critique: scope of unannotated UTR recovery, bacterial versus archaeal distinction, transcript-end versus actual termination claims, host profile calibration, whether reference matching should precede new prediction, and exact minimal reusable architecture seams. Return (1) major changes needed, (2) revised phases with acceptance criteria, (3) remaining decisions. Treat recommendations as advice, not authority. Do not implement anything.
