# TP53 / TP73 cDNA evidence against GRCh38 — 2026-11-10

The supplied regular-cDNA fixture was analyzed separately for SK-MEL-28 (barcode15) and WM793 (barcode16). These are RNA/cDNA sequence observations against the supplied GRCh38 locus sequences. They are not automatically germline, somatic, RNA-editing, copy-number or trans-splicing calls.

Fixture source: `smoe-bot/gentle_rs`, commit `42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649`, directory `test_files/fixtures/cdna_haplotype_tp53_tp73`. At the 2026-11-10 verification, implementation branch `codex/cdna-allele-comparison` was an uncommitted working diff on main `249173ee3b35bfcf6102f3b93b39feb357df553e`.

## Inputs and counting

| Gene | SK-MEL-28 selected reads | WM793 selected reads | Ensembl-116 transcript models | Producer partitions |
|---|---:|---:|---:|---:|
| TP53 | 482 | 481 | 40 | 16 + 16 |
| TP73 | 156 | 10 | 20 | 5 + 1 |

All 18 inventory checksums and all 1,129 BAM/FASTQ original-orientation sequences and Phred qualities matched. No hard clips or duplicate read names were present. Every original record appears exactly once in the 38 producer partitions; stable original-BAM ordinals and digest bindings preserve its identity. Producer partitioning bounds the existing SAM/model projection memory contract; all units are reunited into one sample before SNP assessment and phasing. No cohort was downsampled and no transcript model was removed.

The BAMs contain primary mapped locus-overlapping alignments only. Search/view completeness declarations remain false. All supplied views were assessed; that is a separate property. Missing SA companions affect 17 / 10 TP53 reads and 6 / 2 TP73 reads (barcode15 / barcode16); the latter WM793 reads include references to chr11 and MT. Locus-overlap selection does not establish transcript origin: 9 / 6 TP53 reads and 5 / 5 TP73 reads have no supported block overlapping an annotated exon of the selected gene.

Counts below use distinct producer-bound read units, not verified independent molecules. Selected read counts are not differential-expression measurements. Full parent reference identity is declared and cross-bound; the supplied locus/annotation bytes are reverified, not the complete multi-gigabyte genome.

## Substitution candidates and sample differences

Baseline pilot policy: base Q≥20, MAPQ≥20, at least 3 qualified ALT units and ALT fraction≥0.10 among callable units. A reference-supported base, an ALT observation below nomination threshold and an uncallable position are different states. The 10 TP53 / 14 TP73 positions are union RNA/cDNA substitution candidates, not 24 validated differences between the cell lines.

Genomic positions are one-based below; input/API site positions are zero-based. A/C/G/T counts show all qualified bases at the selected positions. The denominator is the number of callable units at that site, not all selected reads. A zero callability count is not a reference allele.

### TP53 — Q20

| Genomic site, REF→nominated ALT | Models containing site | SK-MEL-28 A/C/G/T | Callable | WM793 A/C/G/T | Callable |
|---|---:|---|---:|---|---:|
| 17:7670141 T→C | 1/40 | 0 / 0 / 0 / 6 | 6 | 0 / 3 / 0 / 0 | 3 |
| 17:7675177 C→A | 38/40 | 293 / 17 / 0 / 0 | 310 | 0 / 215 / 0 / 1 | 216 |
| 17:7675178 A→C | 38/40 | 5 / 293 / 0 / 0 | 298 | 214 / 0 / 0 / 0 | 214 |
| 17:7676154 G→C | 30/40 | 0 / 3 / 254 / 1 | 258 | 0 / 145 / 29 / 0 | 174 |
| 17:7676279 A→T | 0/40 | 2 / 0 / 0 / 3 | 5 | 1 / 0 / 0 / 0 | 1 |
| 17:7676282 G→T | 0/40 | 0 / 0 / 1 / 3 | 4 | 0 / 0 / 1 / 0 | 1 |
| 17:7676483 G→C | 3/40 | 0 / 0 / 4 / 0 | 4 | 0 / 3 / 0 / 0 | 3 |
| 17:7676520 C→T | 3/40 | 0 / 4 / 2 / 6 | 12 | 0 / 3 / 1 / 4 | 8 |
| 17:7687272 C→T | 1/40 | 0 / 6 / 0 / 0 | 6 | 0 / 4 / 0 / 3 | 7 |
| 17:7687509 G→C | 13/40 | 0 / 14 / 74 / 0 | 88 | 0 / 5 / 21 / 0 | 26 |

### TP73 — Q20

| Genomic site, REF→nominated ALT | Models containing site | SK-MEL-28 A/C/G/T | Callable | WM793 A/C/G/T | Callable |
|---|---:|---|---:|---|---:|
| 1:3652352 A→G | 3/20 | 17 / 0 / 3 / 0 | 20 | 0 / 0 / 0 / 0 | 0 |
| 1:3656824 C→G | 0/20 | 0 / 1 / 3 / 0 | 4 | 0 / 0 / 0 / 0 | 0 |
| 1:3656826 T→G | 0/20 | 1 / 0 / 3 / 0 | 4 | 0 / 0 / 0 / 0 | 0 |
| 1:3656829 A→T | 0/20 | 0 / 1 / 0 / 3 | 4 | 0 / 0 / 0 / 0 | 0 |
| 1:3656847 T→A | 0/20 | 4 / 0 / 0 / 0 | 4 | 0 / 0 / 0 / 0 | 0 |
| 1:3672337 G→C | 0/20 | 0 / 0 / 0 / 0 | 0 | 0 / 5 / 0 / 0 | 5 |
| 1:3672437 T→A | 0/20 | 0 / 0 / 0 / 0 | 0 | 5 / 0 / 0 / 0 | 5 |
| 1:3672844 T→C | 0/20 | 0 / 0 / 0 / 0 | 0 | 0 / 4 / 0 / 0 | 4 |
| 1:3722110 C→T | 17/20 | 0 / 0 / 0 / 73 | 73 | 0 / 0 / 0 / 2 | 2 |
| 1:3732315 G→A | 0/20 | 3 / 0 / 0 / 0 | 3 | 0 / 0 / 0 / 0 | 0 |
| 1:3733753 C→A | 13/20 | 41 / 26 / 0 / 0 | 67 | 0 / 1 / 0 / 0 | 1 |
| 1:3734534 C→T | 13/20 | 0 / 0 / 0 / 20 | 20 | 0 / 0 / 0 / 0 | 0 |
| 1:3734562 G→T | 13/20 | 0 / 0 / 0 / 15 | 15 | 0 / 0 / 0 / 0 | 0 |
| 1:3734754 A→G | 13/20 | 0 / 0 / 18 / 0 | 18 | 0 / 0 / 0 / 0 | 0 |

The clearest TP53 sequence differences are the adjacent chr17:7675177 C→A and 7675178 A→C observations in SK-MEL-28, and chr17:7676154 G→C in WM793. At Q30 the corresponding ALT support is 247/255, 250/255 and 67/82, respectively. Adjacent substitutions are correlated observations and require joint/multibase interpretation; this routine does not normalize a complex variant or assign a protein consequence.

For TP73, chr1:3722110 T is observed in **both** samples: 73/73 and 2/2 at Q20. WM793 misses the three-ALT nomination threshold but does not lack the ALT observation. At chr1:3733753 SK-MEL-28 has A41/C26 of 67 callable units; WM793 has only C1/1 at Q20 and no callable unit at Q30. This does not support an allele-frequency or absence comparison. All three WM793-nominated TP73 positions are absent from the 20 selected mature transcript models; they must not be presented as three mature-TP73 allele differences.

## Per-transcript nucleotide annotations

All 40 TP53 and 20 TP73 models are retained. The exported per-sample transcript-allele table projects reference and nominated alternate bases into each model using engine-provided strand-aware coordinates; model absence remains explicit. It does not assign every genomic observation to that transcript as its RNA of origin.

For TP53 ENST00000269305, genomic 7675177 C→A maps to mature-cDNA nt577 G→T, 7675178 A→C to nt576 T→G, and 7676154 G→C to nt357 C→G. These are mature-cDNA positions, not CDS-relative HGVS c. coordinates. For TP73 ENST00000378295, genomic 3722110 C→T, 3733753 C→A and 3734534 C→T map to mature-cDNA nt678, nt2744 and nt3525. Protein consequences, pathogenicity and clinical interpretation are outside this contract.

## Local phase and sensitivity

The final engine forms components from distinct-read-supported local links and phases through identifiable overlapping partial patterns. A whole-block read is not required. Many-to-many overlap ambiguity, contradictory supported patterns and insufficient links remain unresolved; disconnected components are never concatenated. Weak single-read bridges do not force a larger phase block.

Full discovery still leaves the broad exonic/mixed-marker components non-identifiable for both genes; this is now overlap ambiguity, not a requirement for a complete-marker read. The full reports also contain sparse single-site blocks and intronic TP73 sequence classes; these must not be counted as mature-transcript alleles. The scoped local panel results are:

| Gene | Q | Sample | Genomic marker vector | Assigned units | Actual complete-vector units | Minimum transitive-link units | Method |
|---|---:|---|---|---:|---:|---:|---|
| TP53 | 20 | SK-MEL-28 | ACG | 293 | 198 | — | direct_complete_vectors |
| TP53 | 20 | SK-MEL-28 | CAG | 14 | 3 | — | direct_complete_vectors |
| TP53 | 20 | WM793 | CAC | 126 | 121 | — | direct_complete_vectors |
| TP53 | 20 | WM793 | CAG | 27 | 27 | — | direct_complete_vectors |
| TP53 | 30 | SK-MEL-28 | ACG | 248 | 84 | — | direct_complete_vectors |
| TP53 | 30 | SK-MEL-28 | CAG | 7 | 1 | 3 | transitive_supported_local_patterns |
| TP53 | 30 | WM793 | CAC | 53 | 52 | — | direct_complete_vectors |
| TP53 | 30 | WM793 | CAG | 14 | 13 | — | direct_complete_vectors |
| TP73 | 20 | SK-MEL-28 | TAT | 30 | 3 | — | direct_complete_vectors |
| TP73 | 20 | SK-MEL-28 | TCT | 20 | 3 | — | direct_complete_vectors |
| TP73 | 20 | WM793 | unassessed | — | — | — | insufficient supported markers |
| TP73 | 30 | SK-MEL-28 | TAT | 21 | 1 | 7 | transitive_supported_local_patterns |
| TP73 | 30 | SK-MEL-28 | TCT | 8 | 2 | 7 | transitive_supported_local_patterns |
| TP73 | 30 | WM793 | unassessed | — | — | — | insufficient supported markers |

Vectors follow ascending genomic marker order, not contiguous sequence or transcript order. TP53 uses 7675177 / 7675178 / 7676154; TP73 uses 3722110 / 3733753 / 3734534. Assigned counts include qualified partial reads and limited error-tolerant assignments, so they can exceed complete-vector counts. Link witness sets are nonadditive; their union is not a count of whole molecules. H1/H2 display names are local and are not compared across samples.

WM793/TP53 has two well-supported local RNA classes, CAC and CAG, at both Q20 and Q30. SK-MEL-28 is dominated by ACG; its CAG candidate is weaker. At Q30 CAG has one actual complete-vector read plus two transitive supporting links with 5 and 3 units. This meets the pilot threshold but remains fragile and does not confirm a second DNA allele.

SK-MEL-28/TP73 retains TAT and TCT at Q30 through overlapping reads: actual complete-vector counts are 1 and 2, while the necessary link supports are 15/7 and 7/7. This is the material correction to the earlier complete-vector-only pilot. WM793/TP73 lacks a supported exonic panel phase, rather than having zero alleles.

For TP53 ENST00000269305, the three RNA patterns correspond to mature-cDNA bases at nt357 / nt576 / nt577: ACG → C/G/T; CAG → C/T/G; CAC → G/T/G. For TP73 ENST00000378295, TAT and TCT correspond to nt678 / nt2744 / nt3525. Complete per-haplotype, per-model annotations retain absence where a model lacks a marker.

Uniqueness is within the observed supported local-pattern relations. It is not a calibrated estimate of true K, ploidy, or absence of unseen classes.

The optional three-marker panels are **post hoc, data-driven exploratory selections from the same Q20 reads**, not prospectively specified or independently validated. The rule selects exonic reported candidates with at least 20 qualified ALT units for one base in at least one sample. The TP53 panel also has at least 20 callable units in both samples. It selects genomic positions 7675177 / 7675178 / 7676154 (shared by 30 models) and 3722110 / 3733753 / 3734534 for TP73 (shared by 13 models). The identical marker sets are reused for Q30; no K target is imposed. All SNP discovery, raw evidence, counts and models outside each phase panel remain present.

Q30 additionally nominates chr17:7687528 and chr1:3652475 under the same fraction/count thresholds; a stricter Q filter changes the callable denominator and can change a thresholded nomination. This is not a new independent experiment.

## Allele–splice links, switches and excess diversity

Allele/splice associations must use each original read and retained mapping hypothesis. Transitive phase-support witnesses do not become an assembled physical RNA or a union of splice chains. The exports preserve nonadditive same-read skip sets and all emitted switch explanations. CIGAR skips are not automatically qualified junctions or resolved isoforms.

No switch candidate and no excess-K>2 flag was emitted in any of the eight final full/panel Q20/Q30 reports. No inferred K or zero switch count establishes genomic copy number or absence of cross-allelic splicing. The WM793/TP53 classes differ at one site, the TP73 classes at one site, and SK-MEL-28/TP53 only at the two adjacent dinucleotide positions, without a junction between them. The exploratory panels contain only three sites, with even fewer allele-distinguishing markers; they cannot supply two informative markers on both sides of a strong switch. Sparse WM793/TP73 evidence remains particularly uninformative. Any excess RNA-pattern finding remains a candidate for sequencing/alignment/RT/PCR, editing, mixture or copy-number follow-up, not a diagnosis.

## Reproduction and outputs

Use [the fixture preparer](../scripts/prepare_cdna_fixture_comparison.py) with an explicit source-matched headless CLI, a fresh output directory and externally supplied pysam. It invokes the existing GTF normalizer and saved-SAM producer; it contains no SNP caller, GTF/CIGAR geometry engine or phaser. The [contract](cdna_allele_comparison.md) defines partition/source identity and phase-marker scope.

```sh
python scripts/prepare_cdna_fixture_comparison.py \
  --fixture /ABS/PINNED/cdna_haplotype_tp53_tp73 \
  --fixture-commit 42e819c9668e5f5d31fd58e50d3bdd6d1d7f6649 \
  --gentle-cli /ABS/gentle_cli --output /ABS/FRESH_OUTPUT --compare
```

Raw-source audit, the preserved unpartitioned memory-limit refusal, producer receipts, all requests, Q20/Q30 outputs, exploratory selection metadata and table exports live under the run-owned `/private/tmp` paths recorded in the final receipt. Scientific review remains open. MDM4 is outside this analysis.

## Verification and artifact receipt

Tested source is main `249173ee3b35bfcf6102f3b93b39feb357df553e` plus the
uncommitted source manifest
`b55a282dfff02de1ebd8aa9535db27a5fa4fb25a37c3d16d46957b04130ea5db`.
The manifest binds 30 changed source/compiled-contract files; unchanged code and
dependencies are bound by that baseline. It is compact sorted-key JSON, with
sorted path/byte-count/SHA records, hashed without its final newline.

Platform: macOS arm64; Rust `1.100.0-beta.13 (3b0202dcd 2026-11-07)`, Cargo
`1.100.0-beta.13 (3acc79c0c 2026-11-05)`. The isolated target was cleaned before
the initial authorized builds; the compiler remained unchanged during this
follow-up. No other worktree target was cleaned. Cargo used locked/offline
inputs. Evidence preparation used external Python 3.15.2 / pysam 0.24.1,
without adding a core dependency or invoking Sniffles.

| Current-source check | Result |
|---|---:|
| `cargo test -p gentle-engine --lib` | 91 passed, including 26 cDNA tests |
| `cargo test --no-default-features --lib cdna_comparison` | 20 passed |
| Same headless libtest, `sam_intake_` | 12 passed |
| Same libtest, `observed_structures_s2a_` | 4 passed |
| Same libtest, `read_end_service_` | 12 passed |
| `cargo test --no-default-features --test capability_registry_parity cdna_comparison` | 1 passed |
| Default-feature `cargo check -q --locked --offline` | passed |
| Python capability maturity / payload / trace contracts | 17 / 9 / 7 passed; 1 optional legacy CLI skip |
| Python actual cDNA CLI / fixture preparer | 2 / 10 passed |

Totals are **140 unique Rust tests and 45 Python tests passed**, with one
expected optional legacy payload-CLI test skipped. Reruns are not added again.
No full-suite CI or native Windows/GUI acceptance is claimed. Statistical error
calibration, independent scientific review and orthogonal DNA remain open.

The explicit six-record full-SAM versus three-part control passed. Original
read IDs/sequences/qualities (including reverse orientation), six distinct
units, SNP evidence, phase classes and two separate same-read skip sets are
semantically identical after excluding container provenance digests. The real
preparer independently preserved all 1,129 original records and all 23 pinned
fixture files; all 18 data-inventory checksums still match. Duplicate original
source/unit IDs are refused, not counted twice. QNAME alone never joins records.

All eight final real-data runs completed. In each, SNP counts, transcript
annotations, total units and original base/quality evidence are identical to
the corresponding earlier outputs. The witness audit checked 520 actual
complete-seed memberships and 50 transitive-link memberships across repeated
analyses, each against qualified unit calls. Those totals are not molecule
counts. The audit SHA-256 is
`3881627c70db0d7a196900b01075ae67a6a0ea5b94bb544907ec090370020a41`.

Final CLI SHA-256:
`c645bdeaf13af3925bbb59d26b3d11c1789af821164170ba6916002029e2ef3f`.
Headless libtest SHA-256:
`2754876980ccbc71899f11f73be8b98a95dc1326961b3305344ada6ded1322b8`.
The final run is `/private/tmp/gentle-cdna-real-analysis-transitive-20261110`;
its `verification.json`, `source.json`, per-run execution records and output
content hashes bind the receipt. Earlier complete-vector outputs and the first
unpartitioned projection-memory refusal remain separately preserved; neither
is presented as current-method acceptance.

The readable artifact copy is
`outputs/tp53_tp73_cdna_20261110` in the original owner workspace. It contains
all final JSONs, 13 TSV projections, per-model/sample allele evidence, phase and
splice tables, empty switch-explanation table, raw-source and witness audits,
selection metadata, this report and a source snapshot. Serialized file locators
retain their captured run paths; rerun the preparer with fresh paths to replay
the workflow. It is an inspection/reproduction bundle, not an assertion that
moving hash-bound locators preserves their provenance hashes.

Builds retained headless unused-function and macOS unwind-table linker warnings.
One intermediate check during active editing saw an incomplete initializer;
all final stable-source checks pass. Session-close has zero failures and two
expected reminders for intentional uncommitted/generated changes and manual
plan fidelity. No unrelated owner files, main review receipts, model universe
or data thresholds were changed. MDM4 remains excluded; Claude was not invoked.
