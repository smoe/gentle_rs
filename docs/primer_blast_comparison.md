# GENtle and NCBI Primer-BLAST

## Purpose

Primer-BLAST and GENtle overlap, but they are not interchangeable. Primer-BLAST
is an excellent hosted workflow for designing primers and checking them against
NCBI databases. GENtle is a project-state and evidence engine that can use
Primer3 and local BLAST resources while retaining transcript, annotation,
experimental, approval, and artifact provenance.

This comparison is grounded in the original Primer-BLAST paper, the current
NCBI interface guide, and NCBI's related-sequence group-target documentation:

- Ye et al. (2012), [Primer-BLAST: A tool to design target-specific primers for
  polymerase chain reaction](https://doi.org/10.1186/1471-2105-13-134)
- NLM, [Primer-BLAST interface and result
  guide](https://www.nlm.nih.gov/ncbi/workshops/2023-09_Primer-BLAST/interface.html)
- NCBI, [Designing primers for a group of related
  sequences](https://ncbiinsights.ncbi.nlm.nih.gov/2020/09/25/primer-blast-related-sequences/)
- NCBI, [current Group Target
  interface](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?GROUP_TARGET=on)

## Capability Overview

| Capability | Primer-BLAST | GENtle |
| --- | --- | --- |
| Primer3 primer generation | Yes, through the hosted interface | Yes, through internal or Primer3 backends with exact request provenance |
| Existing primer-pair specificity | Yes | Yes, against prepared local genomic-DNA or whole-cDNA databases |
| Related-sequence common primers | Yes; longest sequence is aligned to the others and common regions constrain design | Yes; `DesignPrimerGroupTarget` records every alignment, exact common interval, budget, backend request, and pair-by-member product |
| Organism/database selection | NCBI database and organism controls | Content-fingerprinted local resource catalogs with assembly/release/index-kind provenance |
| Exon/junction-aware assays | Transcript options in Primer-BLAST | Transcript equivalence, endpoint/SYBR/TaqMan modes, required junctions, product matrices, and oligo-dT interpretation |
| Variant/repeat/low-complexity planning | Interface/database dependent | Explicit interval evidence with hard exclusion or advisory priority and independent provenance |
| Fixed-primer completion | Yes | Yes, through side constraints and fixed-primer completion paths |
| Off-target explanation | Primer-BLAST specificity result | Complete-primer alignments, genomic/cDNA target-space separation, reviewed exact allowances, HTML projection, and redesign dispositions |
| Multi-objective alternatives | Primarily one ranked result list | Existing ranking plus a bounded non-dominated Pareto alternative surface |
| Reproducible offline execution | No; hosted NCBI service | Yes, with local Primer3/BLAST, content hashes, request digests, and structured reports |
| Multi-isoform experimental planning | Limited | Gene/transcript coverage universes, endpoint matrices, virtual gels, assay handoff, oligo orders, and publication dossiers |
| Approval and resumable batches | No project workflow contract | Content-bound approval digests, exact operation batches, progress, checkpoints, and reviewed reuse |

## What GENtle Does Not Replace

GENtle does not mirror NCBI's hosted database coverage, Entrez convenience, or
server-maintained sequence corpus. A local GENtle result is only as current and
complete as the prepared resources whose fingerprints it records. Primer-BLAST
also remains a useful independent comparison for ordinary primer pairs.

GENtle does not infer that supplied related sequences are orthologs, paralogs,
or biologically equivalent. The group-target operation consumes the caller's
declared group and reports exact sequence geometry. It does not turn common
sequence into evidence of shared regulation or function.

GENtle deliberately separates design guidance from final specificity. A
whole-cDNA similarity map can steer primer placement away from paralogous
regions, but every selected pair must still enter the complete-cDNA and genomic
specificity gates required by its assay policy.

## Flexibility Adopted From Primer-BLAST

The useful mechanics adopted in GENtle are:

1. design from existing primers or complete a fixed side;
2. expose Primer3 chemistry and positional controls rather than hiding them;
3. derive common binding regions for groups of related sequences;
4. use short-query BLAST as a broad locus finder, followed by complete-primer
   alignment and explicit product construction;
5. distinguish ordinary off-target products from narrowly reviewed exceptions;
6. retain credible alternatives instead of presenting one scalar rank as the
   only answer.

GENtle adds stricter reproducibility around those mechanics: search budgets are
checked before Primer3, member and request content are hashed, intermediate
evidence remains machine-readable, and incomplete searches cannot become a
specificity pass.

## Choosing A Route

- Use `primers design` for one template and one bounded target.
- Use `primers design-group-target` when one pair should amplify a caller-
  supplied group of related sequences.
- Use `primers design-transcript-assay-panel` for annotated isoform coverage,
  discrimination, endpoint matrices, SYBR, or TaqMan work.
- Use `primers specificity-plan` and the transcript-panel specificity routes
  for complete local genomic/cDNA QA and wrapper-owned BLAST execution.
- Use Primer-BLAST as a convenient hosted design/check or independent review,
  while recording its input and output as external evidence rather than
  silently treating it as GENtle-owned provenance.
