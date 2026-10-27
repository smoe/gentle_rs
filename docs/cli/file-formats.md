# File Formats

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

## File format expectations

GenBank and EMBL are peer annotated nucleotide exports through shared `SaveFile`;
FASTA deliberately omits the feature table. See [format parity, exact commands
and limitations](../sequence_format_parity.md).

Current CLI workflows rely on sequence files supported by internal loaders:

- SnapGene DNA (`.dna`)
- GenBank
- EMBL
- FASTA
  - default interpretation: synthetic blunt `dsDNA`
  - optional FASTA-header metadata tokens for synthetic oligos:
    - `molecule=ssdna` for single-stranded DNA
    - `molecule=rna` for RNA (input `T` is normalized to `U`)
    - `molecule=dsdna` plus optional overhangs:
      - `forward_5=...` (alias `f5=...`)
      - `forward_3=...` (alias `f3=...`)
      - `reverse_5=...` (alias `r5=...`)
      - `reverse_3=...` (alias `r3=...`)
- NCBI sequence XML (`GBSet/GBSeq`, `INSDSet/INSDSeq`)
  - other XML dialects are rejected with explicit diagnostics

Example FASTA headers:

- `>oligo_ss molecule=ssdna`
- `>oligo_rna molecule=rna`
- `>oligo_ds molecule=dsdna f5=GATC r5=CTAG`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: CLI Options and Progress](cli-adapter.md) | [Next: Shared Shell](shell.md)
<!-- manual-page-nav:end -->
