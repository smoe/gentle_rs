# Example Asset Provenance

## `cdna_assay_demo.gb`

- Origin: synthetic, hand-crafted GENtle demo locus.
- Recreation: create a 32 bp linear GenBank record with one two-exon `TEST1`
  transcript (`join(1..12,21..32)`) and the sequence
  `ATGAAACCCGGGTTTTTTTTCCCAAATTTGGG`.
- Used by: cDNA PCR/qPCR assay-test workflow and ClawBio direct request
  examples.

## `cdna_assay_nonspecific_demo.gb`

- Origin: synthetic, hand-crafted GENtle demo locus.
- Recreation: create a 72 bp linear GenBank record with two `NONSPEC`
  transcripts sharing exon 1 (`1..12`) and using alternative second exons
  (`21..32` and `41..72`) so the same `AAACCC` / `CCCAAA` primer pair detects
  two transcript-derived cDNA products of different sizes.
- Used by: cDNA PCR/qPCR product-materialization and product-gel workflow and
  ClawBio examples for nonspecific assay product visualization.

## `tp73_promoter_artifact_demo.gb`

- Origin: synthetic, hand-crafted GENtle promoter-signoff demo locus.
- Recreation: create a 249 bp linear GenBank record with one `TP73`-labeled
  gene, two transcript features sharing the same 5' boundary
  (`ENSTTP73DEMO1`, `ENSTTP73DEMO2`), one alternative-start transcript
  (`ENSTTP73DEMO3`), one exact shared promoter annotation (`61..116`), one
  synthetic `SP1` TFBS (`81..89`), one synthetic promoter variant (`95`), one
  synthetic LINE/L1 repeat interval (`105..116`), one TP73-like TFBS
  (`150..158`), and one CUT&RUN-labeled track interval (`130..152`).
- Exact sequence:
  `AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGCGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGATGTGTAACTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA`
- Used by: offline promoter-design artifact slice workflow for
  alternative-promoter grouping, promoter evidence matrix review, TFBS
  score-track SVG export, and TFBS similarity JSON export.

## Promoter cohort tutorial fixtures

Files:

- `promoter_cohort_tutorial_human.fa`
- `promoter_cohort_tutorial_human.gtf`
- `promoter_cohort_tutorial_mouse.fa`
- `promoter_cohort_tutorial_mouse.gtf`
- `promoter_cohort_tutorial_genomes.json`
- `promoter_cohort_tutorial_gene_groups.json`
- `promoter_cohort_tutorial_orthologs.json`

Origin: synthetic, hand-crafted offline tutorial fixtures.

Recreation: create one 300 bp human `chr1` sequence carrying a single
GC-rich `GGGGCGGGGCGGGG` patch inside the 40 bp upstream / 10 bp downstream
TP73 promoter window, one 300 bp human `chr2` sequence with a deliberately
T-rich divergent promoter for `TP73D`, and one 300 bp mouse `chr1` sequence
with a T-rich `Trp73` promoter. Add one-exon GTF gene/transcript rows at
`101..260` with TSS `151`, then point `promoter_cohort_tutorial_genomes.json`
at those local FASTA/GTF files. The gene-group catalog defines the reviewed
`tutorial_p73_promoter_cohort` set (`TP73`, `TP73D`). The ortholog resource
maps synthetic human `TP73` to mouse `Trp73` as a local one-to-one row.

Used by: offline gene-set promoter cohort and ortholog promoter cohort tutorial
generation. The fixtures are intentionally artificial and demonstrate
relationship-flag mechanics only; they do not make claims about real TP73
regulation or orthology evidence.

## Promoter-reporter panel tutorial fixtures

Files:

- `promoter_reporter_panel_demo_source.fasta`
- `promoter_reporter_panel_demo_candidates.json`
- `promoter_reporter_panel_demo_helper_vectors.json`
- `promoter_reporter_panel_demo_request.json`

Origin: synthetic, hand-crafted offline tutorial fixtures. The vector sequence
itself is the existing repository-owned
`test_files/fixtures/reporter_vectors/synthetic_mcs_backbone.gb` fixture and is
explicitly not Promega pGL4.10.

Recreation: concatenate the 52 bp prefix
`ATCGTACGATGCTAGCTACGTTAGCGATCGTACCTGACTGATCGTAGCTAGC`, the 36 bp spacer
`TCGATGACCTAGTCGATCGTAGCATGCTACGATCGT`, the canonical synthetic p53-family
response-element sequence `GGGCATGCCCGGGCATGCCC`, and the 43 bp suffix
`CTAGTCGATGCTAGCATCGATCGTACGATGCTAGCTACGATCG`. The resulting 151 bp FASTA is
selected as one whole-fragment candidate with its motif at zero-based interval
`88..108`. The helper catalog repeats the synthetic vector's exact 240 bp,
circular, MCS `1..70`, and `luc2_demo` `100..180` fixture assertions. The
insert deliberately omits those MCS restriction sites so the panel exercises
one shared directional strategy. The request file selects those records and
uses a disposable `/tmp` output path so the CLI tutorial can plan, inspect, and
then explicitly approve the exact proposal digest.

Used by: `promoter_reporter_panel_planning_offline`, which exercises read-only
proposal generation and exact-vector validation. These records demonstrate
workflow and approval mechanics only; they are not biological TP73 promoter
evidence, a functional reporter design, or a substitute for pGL4.10.

## Regulatory-fragment panel tutorial regions

File: `regulatory_fragment_panel_demo_regions.json`

Origin: deterministic tutorial selections projected onto the repository's
pinned Ensembl-116 SERPINE1 locus fixture. The underlying locus sequence and
transcript geometry come from
`test_files/fixtures/genomic_regions/serpine1_offline/ensembl_116_serpine1_entry.json`;
the three selected spans and their candidate/partner/control roles are
hand-crafted for instruction and are not public regulatory annotations.

Recreation: run `PrepareGeneLocusEvidence` offline for `ENSG00000106366` with
1,000 bp 5' flank, no 3' flank, and output sequence ID
`regulatory_panel_locus`. Capture plus-strand local intervals `500..650`,
`2742..2923`, and `2923..3043` at fixed timestamp `1706000000000` as,
respectively, `reporter_candidate`, `candidate_cis_regulatory_region`, and
`promoter_region` in set `regulatory_panel_regions`, then export the set as
lossless JSON. The content digests in the file are the engine-generated values.
Use region ids `candidate_a_roi`, `partner_b_roi`, and
`minimal_promoter_roi`, retaining the labels, descriptions, notes, and set
label recorded in the asset; those presentation fields also affect content
digests. The executable workflow pins the complete imported records and the
same digests in its planner request, so drift fails validation before planning.

Used by: `regulatory_fragment_panel_planning_offline`, which declares the
first span as candidate A, the second as partner B, and the third as an
explicit promoter-context control. Those roles define testable comparisons;
they do not assert endogenous sufficiency, enhancer/silencer activity, or a
biological partnership.

The workflow regression test
`workflow_examples_regulatory_fragment_panel_is_bounded_read_only_and_replayable`
uses these same records to check eight-member coverage, exact variant geometry,
independent evidence states, JSON/SVG replay, and unchanged project state.
