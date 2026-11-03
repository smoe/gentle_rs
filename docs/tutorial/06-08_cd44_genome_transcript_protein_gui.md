# Compare CD44 Genome, Transcripts, and Protein Products

> Type: `GUI walkthrough`
> Status: `manual/hybrid`
> Network: required for a fresh Ensembl/UniProt fetch

This tutorial uses CD44 to compare three linked coordinate systems in GENtle:

1. genomic exon and CDS coordinates
2. an annotated transcript product
3. the resulting protein axis and optional curated protein features

It contrasts canonical membrane CD44 with the short CD44RC product. It does
**not** treat every extracellular CD44 measurement as evidence for alternative
splicing: proteolytic shedding and extracellular-vesicle-associated CD44 are
separate mechanisms.

## Expected Reference Result

The bounded acceptance run on 2026-11-03 used Ensembl gene
`ENSG00000026508` and UniProt `P16070`:

- canonical `ENST00000428726.8` (`CD44-208`): 742 aa, including the curated
  transmembrane interval at residues 650–670
- `ENST00000278386.11` (`CD44-203`, CD44RC): 139 aa, signal peptide retained,
  canonical transmembrane interval absent

Public annotations can change. Record the returned source versions and explain
any changed transcript version or length rather than forcing these values.

## Step 1: Import the CD44 Locus

Open `File -> Protein Evidence...` and use the Ensembl gene section:

1. query `CD44`
2. species `homo_sapiens`
3. fetch the gene entry
4. import its genomic sequence as `cd44_ensembl_locus`

The resulting sequence should contain mRNA features and annotated CDS ranges.
Do not continue from an N-only or annotation-only placeholder if you intend to
inspect translated products.

Equivalent shared-shell commands:

```text
ensembl-gene fetch CD44 --species homo_sapiens --entry-id cd44_ensembl_gene
ensembl-gene import-sequence cd44_ensembl_gene --output-id cd44_ensembl_locus
```

## Step 2: Open the Canonical Comparison

Use the newly explicit route:

1. `Genome -> Genome / Transcript / Protein Comparison...`
2. choose `cd44_ensembl_locus`
3. enter transcript `ENST00000428726`
4. click `Open Genome / Transcript / Protein View`

Read the view from top to bottom:

- the upper lanes retain genomic exon/intron and CDS geometry
- matching exon-family colors connect genomic structure to the product
- the lower lane uses the transcript-local protein axis
- the details table reports the derived sequence length and derivation mode

The protein axis must end at 742 aa for the reviewed reference result. A
terminal stop codon belongs to CDS geometry but is not an amino acid and must
not appear as residue 743.

## Step 3: Add UniProt Context

Return to the same specialist window:

1. fetch reviewed UniProt accession `P16070`
2. project it to `cd44_ensembl_locus`
3. restrict the projection to `ENST00000428726`
4. open `Open Protein Expert`

The optional external overlay should show the signal peptide, extracellular
region, transmembrane helix, and cytoplasmic region without replacing GENtle's
transcript-derived product as the primary model.

Equivalent shared-shell commands:

```text
uniprot fetch P16070 --entry-id CD44_UNIPROT_P16070
uniprot map CD44_UNIPROT_P16070 cd44_ensembl_locus --projection-id cd44_uniprot_p16070 --transcript ENST00000428726
inspect-feature-expert cd44_ensembl_locus uniprot-projection cd44_uniprot_p16070
```

## Step 4: Compare the Short CD44RC Product

Open the transcript-native view again with transcript
`ENST00000278386`. For the reviewed reference result, verify:

- derived protein length is 139 aa, not 140
- the N-terminal signal peptide is retained
- the canonical transmembrane region is absent because the product ends far
  before residues 650–670
- the unique C terminus is structurally compatible with a secreted product

Optionally create a second UniProt projection restricted to this transcript:

```text
uniprot map CD44_UNIPROT_P16070 cd44_ensembl_locus --projection-id cd44_uniprot_p16070_cd44rc --transcript ENST00000278386
inspect-feature-expert cd44_ensembl_locus uniprot-projection cd44_uniprot_p16070_cd44rc
```

The comparison is evidence for a splice-defined product lacking the canonical
membrane anchor. It is not evidence that a particular sample expresses or
secretes that product.

## Interpretation Boundaries

Keep these mechanisms distinct:

- **splice-defined soluble product:** transcript structure changes the protein
  and removes the canonical transmembrane region
- **ectodomain shedding:** a membrane protein is cleaved by proteases such as
  ADAM or MMP family members
- **vesicle-associated CD44:** extracellular location does not imply a freely
  soluble molecule

GENtle can establish the annotated structure-to-sequence link. Sample-specific
claims still require RNA junction evidence and, where possible, discriminating
protein evidence.

## Checkpoints

- The comparison is reachable from the `Genome` menu and command palette.
- Genomic exon/CDS context and transcript-local protein lanes are visible in
  one native window.
- The canonical reviewed product ends at 742 aa; CD44RC ends at 139 aa.
- Terminal stops remain part of CDS accounting but never extend a protein axis.
- UniProt is visibly an optional external opinion, not the source of transcript
  identity.
- The conclusion distinguishes alternative splicing from shedding and
  vesicle association.

## Related Tutorials

- [06.02 Transcript-Native Protein Expert Sanity Check](./06-02_protein_transcript_native_expert_gui.md)
- [06.04 TP53 UniProt domain mapping](./generated/chapters/06-04_tp53_uniprot_projection_online.md)
