# Repeat-cohort software exercise

`local_project.rmsk.txt` is hand-crafted UCSC rmsk-shaped text, not RepeatMasker
output. LINE/L1 and SINE/Alu labels are deliberately invented parser inputs on
the bacterial AB011549.2 contig. They are **not** human repeat annotations or
evidence that this contig contains LINEs. No transcript/TSS/UTR is invented.

The reference is the existing `LocalProject` in `assets/genomes.json`, using
`../genomes/AB011549.2.fa` and `../genomes/AB011549.2.gb`. Their original public
NCBI provenance is documented in `test_files/fixtures/README.md`: E. coli
O157:H7 Sakai plasmid pO157, accession AB011549.2, 92,721 bp. The two
LINE-labelled BED intervals `[1843,1849)` and `[3877,3883)` both contain `AGCCAG`
in the genomic forward sequence, but have different coordinates and strands.
The third row tests family filtering. Identical DNA must not collapse loci.

Recreation: copy these three literal rows with LF line endings; no downloader
or biological repeat caller is involved. Verify the forward substrings of the
retained FASTA before using a different reference revision. Recreate derived
reports using tutorial 08.18, in a fresh directory with an isolated genome cache.
Do not commit generated reports as invented execution receipts.

Raw rmsk SHA-256 (LF bytes):
`2fe38f29ab6ebdfb20e1ef3a0f5c28f4e01fffed5dfba32d0503c40af4722063`.
The checkout regression verifies this value in both LF/CRLF modes and includes
an unprotected-line-ending negative control; no raw-byte normalization occurs.

Consumers:

- `engine::tests::repeat_cohort_local_project_shell_pipeline_keeps_distinct_loci`
  exercises preparation, shared-shell filtering/cohort construction, saved JSON
  consumption and nonempty motif-window summaries, plus unavailable geometry.
- `scripts.test_tutorial_gap_lessons` checks raw rows/reference substrings.
- `docs/tutorial/08-18_repeat_environment_cohorts.md` teaches the software path
  before using independently prepared, assembly-matched genuine rmsk data.

The GenBank annotation is not a prepared transcript-index source. Its unavailable
projection is an explicit limitation of this exercise, not proof of no bacterial
CDS and not an inference about human promoters. RNA ranking is unimplemented.
Rust execution and native/resource acceptance remain Glen/CI checks.
