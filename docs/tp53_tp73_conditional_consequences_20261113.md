# TP53 / TP73 Conditional cDNA Consequences — 2026-11-13

The owner-authorized run used the existing bound TP53 artifacts first, followed
by TP73. All four saved views per gene were retained: full/panel at Q20/Q30.
There was no new alignment, SNP discovery, phasing, read selection or model
selection. These are RNA/cDNA observations against the supplied GRCh38 loci,
with consequences conditional on a named transcript; they do not establish DNA
alleles, protein expression, copy number or cross-allelic splicing.

**C1/C2 completed for all eight views. Full-cohort C3 isoform linkage is
unavailable:** both linked Q20 attempts reached the existing 100,000,000-work
ceiling. Completed reports explicitly use `linkage: none`, null linkage fields
and `normalized_evidence_reverified: false`. Those fields do not mean zero
compatible reads. No cap, policy, read set or transcript universe was changed.

## Inputs And Model Coverage

| Gene | SK-MEL-28 / barcode15 selected units | WM793 / barcode16 selected units | Models retained | Complete CDS | Noncoding | CDS undefined |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| TP53 | 482 | 481 | 40 | 37 | 1 | 2 |
| TP73 | 156 | 10 | 20 | 16 | 3 | 1 |

Counts refer to distinct producer-bound units, not independently verified
molecules or expression measurements. The original 38 partitions remain one
sample per gene/barcode. Primary-only mapping and incomplete search/view
declarations remain unchanged. Locus overlap alone does not identify RNA origin.

C1 used the existing checksum-bound Ensembl-116 GTF, normalized annotations and
locus FASTAs, with nuclear table 1 and explicit stop-excluded CDS rows. All
models retained versioned accessions; unavailable coding geometry stayed
unavailable. Complete CDS status is separate from incomplete mRNA-end tags.
Protein FASTA and a frozen known-variant manifest were not supplied: predictions
have no independent protein-sequence corroboration or external variant labels.
The parent-genome digest remains declared/cross-bound; supplied locus bytes,
rather than the full genome, were verified.

## TP53: Reference And Alternate Patterns Separately

Representative reference: `ENST00000269305.9`, predicted protein
`ENSP00000269305.4`. Every site is also projected into all other supplied
models; the example does not identify an expressed source isoform.

The adjacent genomic changes chr17:7675178 A>C and 7675177 C>A correspond to
`c.434T>G` and `c.435G>T`. Their actual same-unit combination is
`c.434_435delinsGT`, predicting `p.(Leu145Arg)`. The latter single substitution
alone predicts `p.(Leu145=)`; the joint pattern is therefore essential.

| Quality | Sample | Joint Leu145Arg pattern | Joint reference pattern | Other complete pattern |
| --- | --- | ---: | ---: | --- |
| Q20 | SK-MEL-28 | 292/297 | 5/297 | none |
| Q20 | WM793 | 0/211 | 210/211 | 1/211 synonymous, below support threshold |
| Q30 | SK-MEL-28 | 246/251 | 5/251 | none |
| Q30 | WM793 | 0/187 | 187/187 | none |

These denominators are units callable at both reported codon-marker positions,
not all selected reads. Incomplete units remain present in the reports. Five
reference-pattern units in SK-MEL-28 are a minor RNA pattern, not confirmation
of a second DNA allele. Genomic vectors in the JSON follow transcript order
(7675178, 7675177 for this minus-strand model).

At chr17:7676154 G>C, the same model projects `c.215C>G`, predicting
`p.(Pro72Arg)`:

| Quality | Sample | REF units | ALT units | Callable | ALT status |
| --- | --- | ---: | ---: | ---: | --- |
| Q20 | SK-MEL-28 | 254 | 3 | 258 | below inherited threshold; one other base retained |
| Q20 | WM793 | 29 | 145 | 174 | inherited candidate |
| Q30 | SK-MEL-28 | 121 | 3 | 125 | below inherited threshold; one other base retained |
| Q30 | WM793 | 15 | 67 | 82 | inherited candidate |

Thus SK-MEL-28 has the dominant Leu145Arg RNA pattern, while WM793 has strong
Pro72Arg support alongside reference support at that position. The local phase
patterns remain those of the earlier comparison. No full-length protein allele,
chromosomal cis/trans or zygosity is assembled from separate codon groups.
Protein numbering and availability remain model-specific in the full JSON.

## TP73 And WM793 Undercoverage

Representative reference: `ENST00000378295.9`, predicted protein
`ENSP00000367545.4`. The principal nominated observations on this model are:

| Genomic change | Conditional consequence | SK-MEL-28 Q20 ALT/callable | WM793 Q20 ALT/callable |
| --- | --- | ---: | ---: |
| 1:3722110 C>T | `c.519C>T`, `p.(Thr173=)` | 73/73 | 2/2 |
| 1:3733753 C>A | `c.*674C>A`, 3′ UTR | 41/67 | 0/1 |
| 1:3734534 C>T | `c.*1455C>T`, 3′ UTR | 20/20 | 0/0 |
| 1:3734562 G>T | `c.*1483G>T`, 3′ UTR | 15/15 | 0/0 |
| 1:3734754 A>G | `c.*1675A>G`, 3′ UTR | 18/18 | 0/0 |

The synonymous T observation is shared: WM793 has 2/2 at both Q20 and Q30,
below the three-unit nomination minimum, not absent ALT. SK-MEL-28 has 62/62 at
Q30. At `c.*674C>A`, WM793 has one reference observation at Q20 and no callable
unit at Q30, while SK-MEL-28 has A41/C26 at Q20 and A27/C11 at Q30. This cannot
support an absence or allele-frequency comparison between the lines.

WM793 has only **10 selected TP73 units**: **9/14 nominated sites are uncallable
at Q20; 11/15 at Q30**. Zero callability is unassessed, never a reference call.
Its three nominated positions 3672337, 3672437 and 3672844 are absent from all
20 selected mature transcript models; they are not three mature-TP73 allele
differences. Fourteen Q20 pooled candidates are not fourteen cell-line
differences. No inherited-candidate nonsynonymous TP73 protein prediction was
emitted at Q20 or Q30; this is limited to these models and callable observations,
not evidence that TP73 coding variants are absent. No multi-site TP73 codon
group was emitted; local RNA phase and UTR patterns remain in inherited evidence.

## Execution, Preservation And Remaining Scope

Verification documentation was committed first as `114126e0`. The first linked
attempt at source `5d310958` exposed order-only `linkage_input_drift`: comparison
hashes sorted evidence, while C3 hashed raw lowering order. Repair `984ce08d`
reuses exactly the existing kernel normalization and retains strict input
equality. New regressions preserve comparison bytes and still reject changed
calls. The subsequent linked attempts passed that boundary and reached the
unchanged work ceiling; their requests, stderr and execution receipts remain.

Final executable: source `984ce08df9eec5d6bf0ef956b3ace9b1e3866536`, SHA-256
`f812626ea7620e6c522c27c9b007bc973d75258d5025bb1e38c74b6b0ee0785f`.
TP53 C1 was produced by the original source-matched `5d310958` executable and
rederived/authenticated by the repaired executable; TP73 C1 used the repair.
Original comparison receipt/fixture identities remain in the
[2026-11-10 report](tp53_tp73_cdna_fixture_analysis_20261110.md).

After cleaning only the task-owned target, **109 engine-library and 63 selected
headless cDNA tests passed**. The six Python companion checks passed, including
both actual CLI process cases with zero skips. Default-feature locked/offline
`cargo check -q` passed after allowing the Swift compiler-cache access blocked
by the first sandboxed attempt. Compiler: Rust/Cargo 1.100.0-beta.13, macOS
arm64. Headless unused-code and large unwind-table linker warnings remain;
no native GUI, Windows/new-source CI, package or scientific acceptance is claimed.

Independent audit verified **6,640 site/sample/model/ALT rows**, exact inherited
samples/phase/policy/reference/warnings and actual same-unit codon witnesses.
All **440 bound input files** are byte-identical before/after execution.
The audit's 71,896 complete-witness and 74,480 incomplete-group memberships
repeat across models/views and are not molecule counts. Full/panel Q20/Q30
views reuse the same reads; panels remain post hoc and Q30 is sensitivity.

The run-owned report, eight full engine JSONs, model ledgers, requests,
verification/source/audit receipts and reconstruction scripts are retained at
`outputs/tp53_tp73_consequences_20261113` in the original owner workspace.
Captured paths retain their original bindings; inspection copies do not change
those paths or authorize silently rehashing old evidence.

Next software work is bounded C3 capacity improvement over the same artifacts,
without increasing caps or dropping models/reads. Matched DNA, more informative
WM793/TP73 coverage, calibrated error/K inference and scientific review remain
separate. Absence of prior switch/excess-K flags establishes neither absence of
cross-allelic splicing nor normal copy number. MDM4 remains outside scope.
