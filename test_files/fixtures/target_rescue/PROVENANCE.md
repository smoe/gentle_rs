# target_rescue synthetic fixtures

These are **entirely synthetic** fixtures for the `gentle_cli rescue-screen`
command (the RNA-seq "target-region rescue screen"). They contain no real
genomic data and are safe to redistribute.

## How they were generated

All files were produced by [`generate_fixtures.py`](generate_fixtures.py) with
a fixed RNG seed (`20260704`). Re-running the script reproduces the files
byte-for-byte. The script asserts, at generation time, that the GENEA and GENEB
transcripts share **zero** 11-mers, so a read drawn purely from one gene cannot
accidentally hit the other.

## Files

| File | Role |
| --- | --- |
| `mini_transcripts.fasta` | Three Ensembl-style transcripts with `gene_symbol:` tags: `GENEA` (`ENSTA1.1`), `GENEB` (`ENSTB1.1`), and `GENEC` (`ENSTC1.1`, present but not a requested target). |
| `mini_reads.fastq` | Five reads: `readA_geneA` + `readA2_geneA` (GENEA only), `readB_geneB` (GENEB only), `readAB_ambiguous` (spans GENEA and GENEB), `readN_nohit` (matches nothing). Tests gzip this committed FASTQ into a temporary `.fastq.gz` file so gzip decoding remains covered without committing an ignored `.gz` artifact. |
| `read_allowlist.txt` | A generic read-ID allowlist selecting the two GENEA reads and the ambiguous read. |
| `salmon_unmapped_names.txt` | Emulates Salmon `unmapped_names.txt` (`<name> <code>`): the two GENEA reads, the ambiguous read, and the no-hit read. |
| `salmon_mappings.sam` | Emulates `salmon quant --writeMappings`: `readB_geneB` maps to the GENEB target transcript, `readA_geneA` maps to `ENSTC1.1` (GENEC, not a requested target, so excluded from the target-mapped universe). |

## Intended screen configuration

Requested genes: `GENEA,GENEB,MISSINGX` with `--kmer-len 11`. `MISSINGX` has no
transcript representation and must be reported as missing. `GENEC` is present in
the FASTA but is not requested, so it is ignored.
