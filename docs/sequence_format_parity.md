# Annotated Sequence Format Parity

GenBank and EMBL are peer annotated nucleotide interchange formats in GENtle.
Both serializers consume the same sequence record and INSDC feature table;
choosing EMBL must not remove annotations or change the biology. FASTA is useful
for sequence-only exchange, but is not an equivalent annotated format.

| Content | GenBank | EMBL | FASTA |
| --- | --- | --- | --- |
| Nucleotide sequence | Yes | Yes | Yes |
| Linear/circular topology | Yes | Yes | Not a standard FASTA field |
| Compound, reverse-strand and partial feature locations | Yes | Yes | No |
| Repeated, quoted and valueless feature qualifiers | Yes | Yes | No |
| Stored feature scores, source identifiers and provenance notes | Yes | Yes | Not an annotation table |
| Project reports, review state and GUI layout | Use the project/JSON records | Use the project/JSON records | No |

## Ordinary Sequence Export

In the DNA viewer, **Export Seq** offers GenBank, EMBL and FASTA. The filename
extension selects the format: `.gb`/`.gbk`, `.embl`/`.emb`, or `.fa`/`.fasta`.
The shared operation is also available through GUI Shell, CLI, MCP `op`, and
the operation interfaces in JS/Lua/Python:

```text
op '{"SaveFile":{"seq_id":"my_sequence","path":"sequence.embl","format":"Embl"}}'
op '{"SaveFile":{"seq_id":"my_sequence","path":"sequence.gb","format":"GenBank"}}'
```

`supported_export_formats` advertises `GenBank`, `Embl`, and `Fasta`.
`EMBL` and `embl` are accepted aliases for `Embl` in `SaveFile` JSON. Reopen
either annotated file with the normal `LoadFile`/`load_dna` path. SaveFile is a
file-writing operation and retains the adapters' existing confirmation policy.

Legacy functions explicitly named `write_gb` and specialized report fields
named `genbank` still mean GenBank. They are not renamed or silently changed;
use `SaveFile` for format choice on a loaded sequence.

## Reporter/TSS Exports

Request either annotated format, or both, from the same bound TSS report:

```sh
gentle_cli features tss-tfbs-profiles-export \
  --report /path/to/context-enriched/report.json \
  --output-dir /fresh/tss-formats --formats svg,genbank,embl
```

Every window needs verified context bases. Older contexts without bases must
be enriched from their original score report and context manifest as described
in [Annotated TSS Sequences](tss_tfbs_profiles.md#annotated-tss-sequences).
The export index maps promoter IDs to `.gb` and/or `.embl` files. The receipt
binds their hashes; verification also checks their exact projection from the
source report, not just a self-consistent file hash.

The [integrated-locus compositor](integrated_locus_tss_profiles.md) accepts
`--output-genbank FILE.gb` and/or `--output-embl FILE.embl`. It concatenates the
engine's selected records, validates their DNA against the selected FASTA
digests, and binds all outputs in its receipt. It does not construct a second
annotation model. No existing scientific output is overwritten, and no
rescoring is needed merely to select a different serialization format.

## Fidelity and Limits

Parity means preserved biological content, not identical bytes or an archival
round trip of every database header. Line wrapping, sequence letter case and
insignificant metadata whitespace can normalize. The shared record retains one
primary accession and date; it is not a complete ENA submission model.

EMBL comments retain GenBank DBLINK, reference-description and original SOURCE
metadata where the internal record has no unambiguous native crosswalk. These
are labelled `GENtle original GenBank ...` with JSON payloads and restored by
GENtle on import. In particular, gb-io flattens ORGANISM and taxonomy; GENtle
does not guess a taxonomy boundary to populate EMBL's OC field. Native EMBL
PR/DR cross-reference text is retained as labelled comments, not discarded.

EMBL export requires nonempty IUPAC nucleotide data and a nucleotide molecule
type; protein and gap-containing records fail explicitly. Local identifiers
must be nonempty ASCII without whitespace or semicolons. This slice does not
add contig-only, protein or submission validation. An absent accession version
uses the local-record `SV 0` marker, not a claimed public sequence version.
The deterministic GenBank date `01-JAN-1970` means no date was supplied, not an
experiment date. Missing molecule type follows GENtle's default of DNA in both
writers; an empty GenBank type could otherwise hide circular topology from
token-based readers. Neither writer certifies submission compliance.

Feature predictions remain predictions. A stored PWM score, CUT&RUN signal or
TSS note has exactly the same interpretation in either file; changing formats
does not make it a measured binding probability or validated reporter construct.

Synthetic inline tests exercise both topologies and strands, compound/partial
locations, metadata and qualifiers, shared-shell export, GUI format selection,
TSS report immutability, and receipt tampering. Real private bundles still need
fresh export and independent inspection before publication.

```sh
cargo test --locked --lib annotated_format_parity
cargo test --locked --lib write_annotated_format_parity -- --ignored --nocapture
python3 -m unittest discover -s scripts -p test_compose_locus_tss_profile_pdf.py
```

The explicitly requested interoperability examples print temporary directories
with synthetic linear/circular `.gb`/`.embl` pairs and a receipt-verified
plus/minus TSS bundle for independent readers such as Biopython. They write no
repository fixtures or scientific outputs.
