# FLNA Splicing Interpretation: Draft, Activation Pending

This is a teaching outline, **not an active tutorial or acceptance record**.
Do not put it in the discovery catalog until the activation gate below is met.
The owner has now supplied the retained genuine FLNA bundle, independently
checked in the [archive inspection](../../glen_tutorial_inputs_review_20261012.md).
It supplies source bytes and a stable transcript selector, but no Ensembl
release record, approved teaching exon or native screenshots. Numeric feature
`300` from that snapshot is not a reusable selector.

## Biological Question

Which FLNA exon boundaries and annotated isoforms support a proposed splice
interpretation, and what do exact loaded UniProt links actually establish?

## Required Retained Inputs

The received archive supplies public FLNA ENSG00000196924 annotation, canonical
transcript ENST00000369850.10, loaded reviewed UniProt P21333, complete project
state, view/report and checksum manifests. Before activation, complete the
release and teaching/native evidence rather than fetching an unlabeled substitute.
Record retrieval URLs/date, Ensembl annotation release, exact response/file
hashes and producer SHA/binary/profile. Prefer a release-pinned archive/GTF.
If using live REST, capture `/info/data` before and after and retain the raw
responses; the release query alone does not pin a later response. Do not label
live data as release-pinned without supporting evidence.

The historical run had 36 transcripts and a genomic-minus/local-plus locus.
Neither count nor exact transcript version is a timeless expected result.
Discover the selected transcript by stable annotation identity in the retained
bundle, then resolve its current feature ID. Never copy `300` into learner steps.

## Proposed Learner Checkpoints

1. Open the retained locus; read genomic strand separately from the local
   sequence axis. Compare Ensembl/RefSeq only when both sources are supplied.
2. Resolve the intended transcript and open the existing Splicing expert.
   Preserve the loaded view/report fingerprint and record the feature selector.
3. Inspect donor/acceptor boundary motifs and both canonical flags. Explain
   exceptional markers, annotations and loaded-evidence diagnostics; a `?`
   is not automatically transcript ambiguity or a noncanonical splice site.
4. Compare exact loaded UniProt xrefs and grouped transcripts. Missing loaded
   evidence is unknown, not proof of no protein or an invented equivalence.
5. Select an exon; explain the engine's projected sequence/frame consequence.
   Do not infer biological exon skipping merely from a selectable feature.
6. Save/reopen and export the same view. Check long IDs in hover, coincident
   junctions, mixed strands, dense matrix retention and summary bounds.

## Activation Gate and Glen Return

Return a branch/full SHA, public bundle paths and SHA-256 manifest, selector
resolution receipt, actual shared command lines and a per-check verdict.
Capture native Splicing expert at 2400 px plus a normal usable window size;
retain raw images and semantic snapshots separately from annotated teaching
copies. State Linux/macOS/Windows coverage explicitly. Real font readability
and human scientific interpretation require independent acceptance.

After those inputs land, fill in exact learner commands and expected observations,
link 06.02/07.03/04.08, add a catalog source and provenance, and only then mark
this as a usable manual guide. No numeric feature IDs, fake retrieval dates or
screenshots are supplied by this draft.
