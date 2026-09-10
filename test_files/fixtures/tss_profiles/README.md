# Synthetic TSS Bundle

## Origin And Use

All identifiers, coordinates, transcript memberships and sequences here are
hand-crafted synthetic data, created on 2026-09-10. They were not downloaded
from a genome or annotation resource and are not representations of real genes.
There are no third-party sequence inputs in this fixture. The manually written
`panel.json` pins the real `MA0004.1` / `Arnt` matrix from GENtle's bundled JASPAR
registry (`assets/jaspar.motifs.json`, source URL recorded there). It contains
no copied matrix counts; strict resolution binds the actual counts at runtime.

`src/tss_fasta_bundle.rs` unit tests consume `manifest.json`, `plus.fa`,
`minus.fa`, `selection.json` and `SHA256SUMS`. They test the Stage A reader and
create temporary, checksum-rebound mutations for corruption, geometry,
membership, selection and path-safety cases. The reader returns data only and
does not alter this fixture. Pure validation tests also live in
`crates/gentle-engine/src/tss_profiles.rs`; their inline panel names and
accessions test syntax only, not registry resolution or scientific identity.

This fixture uses the explicit `gentle.tss_fasta_bundle.v1` development schema.
The existing `gentle.target_tss_fasta_export.v1` input format is checked
separately in `src/tss_fasta_bundle/target_export/tests.rs`, using inline
synthetic variants and an opt-in read-only check of Glen's pinned real bundle.
The fixture's one-matrix panel also exercises the typed compute/export path in
`src/engine/analysis/tss_profiles_tests.rs`, including SVG/PNG/PDF and report-only
replay. It is a software smoke input, not a promoter prediction. The panel is
bound separately by the computation receipt, not by the FASTA `SHA256SUMS`.
Recreate it verbatim from the committed JSON; it requires no network retrieval.

## Deterministic Recreation

The committed text files are the canonical source, encoded as UTF-8 with LF
line endings and a final LF. No generator or external retrieval is needed.
Recreate the records using these exact values and the identifiers and reference
strings in the manifest:

| Promoter | Strand | TSS | Inclusive interval | Window | Canonical sequence |
| --- | --- | --- | --- | --- | --- |
| `synthetic-plus` | `+` | 100 | 97..107 | -3..+7 | `ACGTNACGTAC` |
| `synthetic-minus` | `-` | 300 | 293..303 | -3..+7 | `TTACGACGTAA` |

Both records are 11 bases long, including the TSS base. The plus FASTA stores
the sequence as the two lines `acgtn ac` and `gtac` to test case and whitespace
normalization. Its header lists `synthetic.plus.tx1,synthetic.plus.tx2`, while
the manifest deliberately preserves the opposite membership order. Equality
is by transcript set; returned membership order is the manifest's order. The
minus FASTA stores `TTACGACGTAA` on one line, already oriented from transcript
5-prime to 3-prime. Do not reverse-complement it during reading. For local
indices 0, 3 and 10 its genomic coordinates are respectively 303, 300 and 293.

All FASTA headers have exactly these keys, separated by `|`, with nonempty
values and no duplicate or unknown keys:

```text
gene_symbol=...|gene_id=...|promoter_id=...|assembly=...|chromosome=...|strand=...|tss_1based=...|genomic_1based=START..END|window=-U..+D|orientation=transcript_5prime_to_3prime|transcripts=comma-list|sequence_sha256=...
```

To recreate normalized sequence digests, remove ASCII whitespace from the
sequence lines, uppercase A/C/G/T/N and hash those bytes without a trailing
newline. No other ambiguity codes, RNA U, substitutions or Unicode whitespace
are accepted. All SHA-256 strings are 64 lowercase hexadecimal digits with no
prefix. For example:

```sh
printf '%s' ACGTNACGTAC | shasum -a 256
printf '%s' TTACGACGTAA | shasum -a 256
```

Set `sequence_sha256` to those hashes in the manifest and matching FASTA
headers. Then hash the exact FASTA file bytes with `shasum -a 256` and put those
values in `manifest.json`'s `fasta_files`. Finally regenerate `SHA256SUMS`
using this command from this directory:

```sh
shasum -a 256 manifest.json minus.fa plus.fa selection.json
```

Those output lines, in that order, are the checksum file's exact contents.
Verify with `shasum -a 256 -c SHA256SUMS`. Do not include `SHA256SUMS` in its own
inventory. This README is not a runtime input and is not checksum-bound.

## Reader Contract

The manifest's parent directory is the bundle root. Member paths are exact
portable relative paths without `.`/`..`, absolute prefixes, backslashes or
symlink escapes. Caller FASTA/selection arguments use those names or absolute
paths inside that root. If a FASTA list is supplied, its set must exactly match
`fasta_files`; duplicate requested names fail. No subset or basename matching
is inferred. Returned records preserve manifest order, not filename or selected
order. Without an explicit selection argument every record is unselected;
with `selection.json` only `synthetic-minus` is selected.

`SHA256SUMS` must bind the manifest, every FASTA member and any requested
selection. Every listed auxiliary file is also contained, size-bounded and
hash-verified. Exact file hashes and normalized sequence hashes are separate
checks. Duplicate JSON object keys are rejected, including repeated filenames
in `fasta_files`.

Read limits are 8 MiB for the manifest and requested selection, 1 MiB for
checksums, 32 MiB per other member, and 64 MiB for all consumed bytes together.
The limits include 1024 FASTA files, 4096 checksum members, 10000 TSS records,
one million bases per window, 64 KiB per FASTA header and 4096 distinct
transcript memberships per record. Coordinates must fit checked arithmetic
and the protocol's signed relative axis; contig-end clipping is not inferred.

Repeated promoter IDs or physical TSSs (reference, chromosome, strand, TSS)
fail. Identical sequences at distinct valid loci survive with a diagnostic;
the tests construct this case without adding an extra biological fixture.
The v1 record has one gene ID, so multi-gene membership at one physical TSS is
an explicit unsupported case, not full multi-gene support. Such duplicate rows
fail rather than silently merging or discarding gene/transcript memberships.

Bundle verification proves internal consistency only, not reference genome or
annotation authenticity. No downloads or independent genome extraction occur.
Annotated transcript starts are TSS candidates, not experimentally established
initiation sites. Any later JASPAR profiles are sequence-model predictions,
not measured binding, affinity, occupancy, cofactor interaction, promoter
activity or functional regulation. Cross-matrix magnitudes are not comparable
without documented calibration.
