# Base And Prime Editing: From An Exact Allele To Interpretable Evidence

**Inactive teaching draft, 2026-10-20.** This is a pen-and-paper tutorial, not
an active GENtle walkthrough or laboratory procedure. GENtle does not yet
execute the editing-design or editing-outcome-import steps below. No editor
profile, guide, delivery method or experimental result is approved here.

## What You Will Learn

Before choosing an editor, write down what must change and what must remain
unchanged. Then ask three different questions:

- Can an explicitly selected editing system represent the intended change?
- What other changes or unresolved evidence must the design disclose?
- What did the measured alleles actually show, and which reads were counted?

By the end, you should be able to describe a substitution or short indel
without reversing its meaning on a minus-strand display, distinguish an
intended-only allele from one with bystanders, and explain why an amplicon
percentage is neither a clone genotype nor a safety verdict.

Allow about 20-30 minutes. You need only this page and a calculator; no project,
download, account, model, ordering service or laboratory sample is needed.
All DNA and counts here are hand-crafted synthetic teaching data, not PATZ1,
TP73, a public genome, a real editor's targeting window or a measured experiment.
The short strings deliberately cannot serve as complete editing reagents.

## Where GENtle Stands Today

Existing [guide filtering and oligo formatting](../tutorial/generated/chapters/04-04_guides_filter_and_generate_oligos.md)
teach practical constraints, not editor-specific feasibility. The separate
[synthetic sequence-design tutorial](../tutorial/06-07_synthetic_sequence_design.md)
teaches in-silico redesign, not a mechanism that makes the cell acquire it.
[Sequencing confirmation](../tutorial/10-01_sequencing_confirmation_trace_cli.md)
can inspect local expected/read evidence, but does not supply the editing
allele-count analysis proposed here.

These existing chapters do not form a completed editing pipeline. Do not run
proposed names from the [preparation brief](../base_prime_editing_preparation.md)
as commands, or ask an agent to invent them. In a future implementation, GUI,
CLI, scripting and agent interactions must inspect the same engine report.

## 1. State The Exact Change Before Choosing A Tool

Write the reference sequence in genomic-forward orientation:

```text
Synthetic reference:  A C G T A C
Local base index:    0 1 2 3 4 5
```

Coordinates below are zero-based and half-open: `[1,2)` contains base 1,
not two bases. Each row is a separate exercise starting from the original
`ACGTAC`; do not apply them successively.

| Request | Reference bases replaced | Inserted bases | Expected complete sequence |
| --- | --- | --- | --- |
| Replace `[1,2)` | `C` | `T` | `ATGTAC` |
| Insert at boundary `[3,3)` | Empty | `GG` | `ACGGGTAC` |
| Delete `[2,4)` | `GT` | Empty | `ACAC` |

**Checkpoint:** Why is the insertion interval empty? Because it identifies
a gap between reference bases; it does not consume a reference base.
The reconstructed sequence is the unchanged prefix, inserted bases and
unchanged suffix. A request claiming `G` at `[1,2)` must be rejected, not
"fixed" by choosing a different position or complement.

For real data, add organism, assembly and contig accession/version, sequence
digest and selected haplotype. A chromosome name or matching gene symbol
does not identify the exact allele. A reference haplotype is a design basis,
not evidence of the actual cell line's genotype.

## 2. Read A Minus-Strand Display Without Changing The Allele

Imagine the six-base teaching reference occupies synthetic genomic interval
`[1000,1006)`. This is a made-up coordinate system, not a genome accession.
The reference `C` is at genomic `[1001,1002)` (one-based position 1002).

| View | Sequence in that view's forward direction | Interval containing the target | Letters shown before/after |
| --- | --- | --- | --- |
| Genomic-forward reference window | `ACGTAC` | Local `[1,2)` | C-to-T |
| Reverse-complemented local display | `GTACGT` | Local `[4,5)` | G-to-A |

The edited reverse display is `GTACAT`, the reverse complement of `ATGTAC`.
Those are two descriptions of the same requested allele, not two different
edits. Do not index a six-base window with the absolute coordinate 1001.

**Checkpoint:** A length-`n` reverse view maps `[s,e)` to `[n-e,n-s)`.
An insertion boundary `k` maps to `n-k`. Use the original reference frame
to describe an edit; after an indel, the edited sequence has a different
length and needs its own mapping. Gene strand, view orientation and the strand
to which a guide hybridizes must not be treated as interchangeable labels.

## 3. Ask What A Base Editor Could Change, Not Just What You Want

Conventional cytosine base editing (CBE) and adenine base editing (ABE) act
on a declared exposed DNA strand. In genomic-forward allele notation:

| Selected editor class | Exposed plus strand | Exposed minus strand |
| --- | --- | --- |
| Conventional CBE | C-to-T | G-to-A |
| Conventional ABE | A-to-G | T-to-C |

This is a conversion table, not an editor profile or prediction of efficiency.
A C-to-G request is not provided by these conventional conversions. That does
not mean no other editing chemistry exists; it is outside this initial scope.

For a letter-count exercise only, use `ACCCGTAA` and pretend that a synthetic
rule marks zero-based positions 1, 2 and 3 as potentially editable. Choose
position 2 as intended. Positions 1 and 3 are potential **bystanders**: changing
the intended `C` does not authorize silently dropping the two neighboring `C`s
from the report. This eight-base example has no validated spacer, PAM or
biological editing window. It cannot establish even a feasible guide.

**Checkpoint:** A real feasibility report would need the exact editor construct,
PAM and spacer geometry, editable positions and sequence-context evidence,
with source/correction provenance. It would list bystanders and available
transcript/CDS/splice annotations separately. Missing annotations remain
unknown; "synonymous" does not mean safe. Possible edits are not measured
probabilities, and a practical-filter pass is not an off-target assessment.

## 4. Separate Prime-Editing Components From Cloning Decorations

Prime editing adds a template-bearing extension to the guide RNA. Read the
following as a component diagram, not an actual reagent sequence:

```text
RNA 5' -> spacer | scaffold | RTT | PBS | optional stabilizing element -> 3'
```

The primer-binding site (PBS) pairs with DNA exposed at the nicked 3' end.
The reverse-transcription template (RTT) encodes the intended edit and relevant
homology. Its interpretation must reconstruct the desired DNA on the correct
strand. Swapping the PBS and RTT, or merely reversing a displayed extension,
does not establish a valid design. Keep the DNA encoding of an RNA separate
from the RNA itself.

PE2 has no second nicking guide. PE3 adds a separately checked opposite-strand
nick; PE3b adds an edited-versus-reference recognition condition, not a
guarantee of allele specificity. Protein architecture and RNA stabilization
are separate choices. Do not substitute an ordinary nuclease's cut offset
for a prime editor's nick geometry.

**Checkpoint:** Before accepting any component map, require the declared PBS
pairing, nicked strand and independently reconstructed edited allele to agree.
An extra PAM-disrupting substitution is still an extra edit requiring review.
Promoter-added bases, cloning overhangs and chemical modifications are not
part of the genomic targeting sequence merely because they share an oligo.
There is deliberately no cloning recipe, PBS/RTT length preset or materialized
component in this lesson.

## 5. Interpret Counts With The Denominator Still Attached

Suppose a synthetic run begins with 100 reads. Ten are excluded by its declared
quality policy, leaving 90 retained reads. Five retained reads remain
unclassified because the complete allele cannot be called; 85 are callable.

| Disjoint category | Retained reads |
| --- | --- |
| Intended-only | 20 |
| Intended-plus-declared-bystander | 10 |
| Other substitutions | 5 |
| Unexpected indels | 5 |
| Unedited | 45 |
| Unclassified | 5 |
| Total retained | 90 |

For intended-only alleles among callable reads, the fraction is `20/85`,
approximately **23.53%**. Among retained reads it is `20/90`, approximately
**22.22%**. Among all input reads it is `20/100`, or **20%**. None should be
printed as an unqualified "editing efficiency".

The ten intended-plus-bystander reads are not intended-only. The two "other"
categories do not tell us whether they also contain the intended edit; without
per-allele flags, total intended-edit presence cannot be calculated from this
table. An intended insertion/deletion is not itself an unexpected indel.
Classification needs explicit precedence so combined outcomes are not counted
twice, and unclassified reads must remain visible.

**Checkpoint:** With zero callable reads, the intended-only fraction is
unavailable, not 0%. No untreated control is supplied in this exercise, so no
background correction can be inferred. Reads are not independent biological
replicates. The same bulk frequencies can arise from different cell genotypes;
they cannot establish a clone's zygosity, genome-wide specificity, absence of
guide-independent DNA/RNA editing or biological function.

## 6. Decide What Evidence Is Still Needed

Before turning this conceptual exercise into a real GENtle workflow:

- Bind an exact public reference and explicit allele/haplotype; have a human
  select and review the example rather than borrowing a paper's gene label.
- Select a reviewed editor profile and its cell/delivery applicability.
  Preserve protocol corrections; a citation is not an approved procedure.
- Validate the candidate and reconstructed allele with the shared engine.
  A bounded or cancelled search cannot prove that no design exists.
- Keep a read-only component handoff separate from an approved mutation that
  creates derived records and lineage. A design is not permission to order.
- Import version-pinned local sequencing results with raw-file, parameter,
  reference/amplicon and sample/control bindings. Do not upload sequences
  automatically to a web tool.
- Ask Glen for consolidated native presentation/reopen/cancellation checks
  when those capabilities exist. Human scientific and laboratory validation
  remain separate from software tests and screenshots.

The [preparation brief](../base_prime_editing_preparation.md) records the staged
implementation and primary literature, including the corrected base-editing
figures and prime-editing cloning schematic. For the component concepts, see
[Doman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9799714/); for local
analysis conventions, see the [official CRISPResso2 documentation](https://docs.crispresso.com/latest/parameters.html).
No external tool is installed, run or accepted by this lesson.

## Teaching Data And Review Boundary

The following JSON is **teaching data, not a GENtle request**. It repeats the
small exercises so fast independent tests can check arithmetic and coordinate
consistency. To recreate it, transcribe the reference, separate edit rows,
reverse-view example, synthetic bystander rule and count table above. No
retrieval, random seed, scientific fixture or measured read file is involved.
It is consumed only by `scripts.test_tutorial_gap_lessons`, not by the engine.
Whitespace and LF/CRLF have no semantic meaning for these teaching data.

<!-- editing-teaching-data:start -->
```json
{
  "reference": "ACGTAC",
  "genomic_start_0based": 1000,
  "edits": [
    {"start": 1, "end": 2, "deleted": "C", "inserted": "T", "expected": "ATGTAC"},
    {"start": 3, "end": 3, "deleted": "", "inserted": "GG", "expected": "ACGGGTAC"},
    {"start": 2, "end": 4, "deleted": "GT", "inserted": "", "expected": "ACAC"}
  ],
  "reverse_display": "GTACGT",
  "reverse_target_interval": [4, 5],
  "reverse_edited_display": "GTACAT",
  "toy_bystanders": {"sequence": "ACCCGTAA", "editable_positions": [1, 2, 3], "intended_position": 2, "expected_bystanders": [1, 3]},
  "counts": {
    "input": 100,
    "quality_excluded": 10,
    "retained": 90,
    "callable": 85,
    "categories": {"intended_only": 20, "intended_plus_bystander": 10, "other_substitutions": 5, "unexpected_indels": 5, "unedited": 45, "unclassified": 5}
  }
}
```
<!-- editing-teaching-data:end -->

Codex review, 2026-10-20: readability, preparation-brief consistency and synthetic
arithmetic only. Engine/source inspection baseline: `050c1030`; the separately
committed preparation brief records its own source scope. Human scientific
review, native GUI acceptance, actual editor
feasibility and laboratory validation are **pending**, not implied by this
review. The source/catalog/manifest registries remain unchanged.

This file lives outside `docs/tutorial` because Help's missing-catalog fallback
recursively scans that tree. Do not activate it, assign a numbered chapter or
add executable GUI/CLI/agent steps until their capabilities and corresponding
acceptance boundaries exist.
