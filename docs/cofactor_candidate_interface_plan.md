# Simple Cofactor Candidate Interface

## Objective

Answer one question first: **Which predicted factor sites near an experimentally
supported p73 region should we consider for follow-up experiments?**

The first screen should be a small candidate list, not a regulatory atlas or
decision-tree exercise. Regulatory annotation, alternative transcripts, TATA
analysis and promoter similarity remain optional context, not prerequisites.
Motif proximity generates a hypothesis; it does not identify a bound protein
or establish cooperation, enrichment or reporter activity.

## Prepared Now

`Genome > Gene Set Inspector... > Regulatory partners` now presents the
existing shared `gene-sets regulatory-partner-screen` result as a candidate
table. It shows the predicted anchor and partner sites, exact matrix IDs,
scores, signed distances, the engine's proximity outcome and promoter-level
CUT&RUN state. Expression reads `Not evaluated`, never zero or unexpressed.

The table preserves engine order and every displayed tuple. It does not merge
motifs into protein identities or count multiple tuples as independent peaks.
Nearby-only filtering is reversible; excluded/unresolved genes remain visible
through a status section. Clicking a gene opens its existing promoter DNA
detail; the decision tree and exact-hit audit are expandable. TSV copy includes
join IDs and source-report digests; full JSON retains the analysis policy.

This is a presentation slice, **not yet the full experimental-peak workflow**:

- Current regions are defined around a resolved transcript TSS.
- Anchor-to-partner distances are between predicted motif centres.
- CUT&RUN support is associated with the promoter, not an individual motif.
- `Use 150 bp` changes only the declared motif proximity parameter. It must not
  masquerade as selecting a 500 bp experimental neighbourhood.

The helper lives in `crates/gentle-render/src/regulatory_partners.rs`; the GUI
stays thin in `src/app/gene_set_ui.rs`. Engine scanning, protocol payloads and
shared-shell defaults are unchanged. CLI/MCP/scripts retain the same full
ledger. Native acceptance and the added Rust tests still need Glen/CI.

## Next Slice: Experimental Anchor

Keep the simple interaction to three steps:

1. **Choose evidence.** Select a gene/region and a content-bound CUT&RUN source,
   with assembly, condition, replicate/control and peak-calling policy visible.
2. **Review the anchor and window.** Default to a 500 bp neighbourhood and
   a 150 bp proximity rule, but require an explicit anchor type and position.
3. **Review candidates.** Show the compact table, open sequence context on
   selection, and export a reviewed shortlist for experimental planning.

Before wiring these controls, extend shared engine request/report contracts:

- An experimental peak interval, a reported summit, a user-reviewed point and a
  predicted p73 motif are different types of anchor. A BED interval without a
  summit must not silently acquire one. A user-selected interval midpoint is
  labelled as such, not described as an observed binding position.
- The 500 bp extraction is centred on the declared point, not the TSS. Record
  requested and actual bounds, clipping and reference/sequence digests. If no
  point is available, ask for review rather than choosing another anchor type.
- Internal coordinates remain zero-based half-open, with one-based genomic
  intervals for human display. Declare the direction of signed distance;
  preserve half-base motif centres, reverse-strand projections and overlaps.
- If a p73 motif is selected within a supported region, retain both the
  experimental interval and motif record. Display which one defines distance;
  occupancy must not be reassigned to the motif by proximity alone.
- Reuse the exact shared TFBS scanner. Bind the motif catalog/version, resolved
  matrix IDs, thresholds and background/scoring policy. Unknown or truncated
  coverage is not a negative finding. No qualifying motif is a valid outcome.
- Gene and transcript associations are context, not distance anchors. Allow
  distal regions. Preserve multiple gene/transcript associations without
  duplicating a genomic peak or silently choosing the longest transcript.
- Surface the same operation and portable report through shared shell, MCP,
  scripts and the inner agent before adding a prominent GUI action. Do not
  introduce a GUI-only calculation or guess command names in agent prompts.

No new operation/schema names are advertised until that contract is reviewed.

## Expression and Prioritisation

Add optional expression evidence through exact factor identifiers and a
declared biological context: organism, cell model, condition, contrast, unit,
source report and digest. Shared-family motifs can map to multiple possible
factors; do not guess one expression row from a display name. Distinguish
missing, ambiguous, measured zero and below-detection observations.

Initially show independent reasons: proximity under the stated rule,
experimental region support, and expression status. Do not invent a composite
confidence score. Aggregate views must distinguish unique motif sites,
anchor-partner pairs, experimental peaks and genes, with their denominators.
Claims of overrepresentation require a separately declared background and
statistical method; raw counts alone are descriptive.

Reporter follow-up starts only after the user selects a candidate region and
reviews its boundaries. Existing region/report and reporter-planning contracts
should carry that choice; no automatic sequence materialisation or order is
part of this interface preparation.

## Acceptance

- Run the added headless render tests: `cargo test -p gentle-render regulatory_partners`.
- Run the existing/new GUI tests: `cargo test --lib regulatory_partner_`.
- Before handoff, Glen/CI should run `cargo check -q` and native GUI acceptance:
  candidate list first, reversible filtering, no-pair/unresolved states,
  TSV copy, gene selection and advanced evidence. Check a narrow window and
  a large ledger; off-screen rows must not all be painted.
- For the next engine slice, cover absent/ambiguous anchors, assembly/digest
  mismatch, clipping, forward/reverse orientation, half-base centres, duplicate
  transcript associations, no-hit and incomplete scans, expression ambiguity,
  deterministic JSON/TSV and GUI/CLI/MCP parity on the same fixture.

This preparation was not Claude-reviewed; an optional read-only consultation
was offered. It makes no real-data candidate ranking or biological sign-off.
