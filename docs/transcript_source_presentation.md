# Source-Coherent Transcript Presentation

GENtle can join explicitly supplied Ensembl and RefSeq annotations. This is
source annotation, not a consensus, preferred TSS or transcript-selection rule.
CUT&RUN values, TFBS scoring, selected windows and primer design are unchanged.

## Workflow

Add `transcript_annotation_sources` to the existing `gene-locus prepare
@request.json` request, or `GeneLocusEvidenceDisplayRequest` in a shared engine
operation. Each source is a local, hash-bound GFF3 or saved
`gentle.ensembl_gene_entry.v1` JSON. No download occurs. One illustrative member:

```json
{
  "path": "annotations/source.gff3",
  "sha256": "<64 lowercase hex digits: exact annotation file bytes>",
  "provider": "ref_seq",
  "format": "gff3",
  "assembly": "<exact locus assembly>",
  "release": "<source release declaration>",
  "accession": "<accession-pinned source resource identifier>",
  "chromosome": "<sequence ID as written in the annotation>",
  "locus_sequence_sha256": "<SHA-256 of the loaded locus sequence>",
  "gene_ids": ["<exact GFF3 gene ID, including provider prefix>"]
}
```

Providers: `ensembl`, `ref_seq`. Formats: `gff3`, `ensembl_gene_entry`. Input
hashes accept `sha256:`; report hashes are unprefixed. Use the existing locus
`sequence_binding.sequence_sha256`, not the FASTA file hash: headers and wrapping
are not sequence. Supply gene-scoped annotation with parent links intact;
limits are 128 MiB/file, 16 sources and 256 transcript records per presentation.

The locus JSON carries `transcript_presentation`, schema
`gentle.transcript_structure_presentation.v1`. Use **that same JSON** as
`locus_report` in the existing TSS `--context-manifest`. Overview and detailed
pages then consume the same canonical presentation. FASTA and locus JSON keep
their independent existing hash bindings.

For detail-only inspection, the context manifest's gene entry can itself carry
`transcript_annotation_sources`, with paths relative to the manifest. This
joins in memory and does not modify its source locus JSON or old overview.
For a coherent complete dossier, prepare/export the enriched locus first.
An already-attached presentation cannot be silently replaced by another join.
The compositor rejects differing overview/detail presentations or SVG pages
without the matching presentation digest, even if their individual file receipts
have been refreshed.

These fields travel through existing shared operations in CLI/shell, workflows
and MCP/JS/Lua consumers. No new command or inner-agent behavior is introduced.
Original per-transcript GenBank/EMBL annotations and the native flat-file viewer
are not rewritten by this SVG presentation layer.

## Identity And Display

- Physical exon IDs bind assembly, chromosome, strand and exact inclusive
  interval. Source exon IDs stay attached to every member record.
- Exon-chain IDs bind the full ordered chain, not its visible/cropped part.
- CDS IDs separately bind intervals, known phases and availability. Unavailable
  CDS differs from explicitly noncoding. Known versus unavailable phase keeps
  separate rows conservatively; it is not silently treated as equal.
- Structure IDs combine chain and CDS IDs. Identical structures share a row;
  one shared exon cannot collapse different chains. Every source transcript,
  designation and content hash is retained.

Shared physical exon boxes are drawn once, stacked where intervals overlap.
Distinct CDS geometries also appear once. Structure rows reference the full
chain via figure-local `E1`, `E2`, etc.; these are not biological exon numbers
or stable IDs. SVG hover exposes all attached source transcript/exon IDs,
designations, release/accession, geometry and hashes. Off-window differences
still keep distinct structure rows.

Separate colour-labelled Ensembl/RefSeq annotated-start lanes use `=` for exact
agreement, `s` for a source-only coordinate in supplied annotations, and `?` when
the other provider was not supplied. All supplied Ensembl x RefSeq start pairs
have a typed delta: `(RefSeq - Ensembl) * strand`. Thus +10 genomic bp is -10 bp
on a minus-strand gene. No nearest/preferred start is chosen. Comparisons
touching the view are printed (first 16, explicitly labelled); all remain in
JSON/hover. Ticks use the same base-centre axis as CUT&RUN and TFBS.

## Coverage And Mapping Limits

`detail_context.transcript_payload_coverage` lists requested, included and
unassessed exact IDs with one sentence. An omitted payload record is not a
missing transcript or biological absence. A missing legacy coverage field is
not a measured zero. Unresolved gene/source IDs, invalid geometry/strand,
assembly mismatch and failed hashes fail the join; partial joins are not published.

- GFF3 resolves direct gene -> transcript -> exon/CDS parents and retains
  versions/source exon IDs. Explicit `tag` values supply Ensembl canonical,
  MANE Select/Plus Clinical and RefSeq Select, retaining field and value. No
  designation is inferred from an accession, transcript length or rank.
- RefSeq reference accessions map through a source GFF3 `region` record's
  `chromosome` attribute. Otherwise only exact IDs and the `chr` prefix convention
  match. Include that region record; there is no guessed accession map/liftover.
- Ensembl entry JSON supplies `is_canonical` and exon/translation geometry but
  lacks normalized release, MANE and CDS phase. Release/accession are explicit
  request declarations bound to the file; phases stay unavailable. Use GFF3
  tags for MANE/RefSeq Select rather than assuming them.
- A GFF3 `#!genome-build` header must match the assembly when present. Release
  declarations are provenance, not live-provider verification. Partial,
  start/end-range and exception attributes are retained as notes; annotated
  starts do not prove a complete experimentally observed TSS.
- IDs are exact, including version. Unversioned/differently versioned links are
  not silently equated. Structural identity does not prove cDNA/protein identity.

Serialized derived identities are rebuilt and compared at ingestion/rendering.
This detects stale/inconsistent content, not an attacker replacing every input
and receipt with a new self-consistent dataset.

## PDF And Interactive Pages

Keep the existing `scripts/compose_locus_tss_profile_pdf.py` command and PDF.
Optionally add:

```text
--output-svg-directory /new/output/gene-pages
--output-svg-zip /new/output/gene-pages.zip
```

The ZIP is optional. The directory contains ordered `page_0001.svg`, etc.,
`index.html`, `pages.json` and selected-sequence FASTA. HTML displays native SVG
objects with hover, without recomputing annotations. PDF receives the same
staged SVG bytes in the same order. The outer receipt hashes the files/order
and ZIP; it is not recursively included in its own archive. Existing outputs
are not overwritten; failed publication leaves no success receipt.

Real-data regeneration and scientific acceptance remain Glen's task. Tests use
hand-crafted plus/minus records and temporary annotation files, not private data.
