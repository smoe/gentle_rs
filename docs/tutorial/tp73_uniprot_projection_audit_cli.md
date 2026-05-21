# TP73 UniProt Projection Audit

See also: executable reference chapter
[`19 Audit a TP73 UniProt Projection Against Ensembl and Derived Coding Sequence`](./generated/chapters/19_tp73_uniprot_projection_audit_cli.md).

This walkthrough shows two ways to inspect a projected TP73 UniProt entry:

1. the integrated audit path, and
2. the same result rebuilt from the reusable public primitives.

The audit stores a local report plus an unsent maintainer-email draft. It does
not send email.

## Preconditions

- TP73 genomic context is already in the project as a sequence with transcript
  features.
- a TP73 UniProt entry is imported or fetched into the project
- a TP73 Ensembl protein entry is fetched into the project

Example fetch/import sequence:

```text
uniprot fetch Q9H3D4 --entry-id TP73_UNIPROT
ensembl-protein fetch ENSP00000264724 --entry-id TP73_ENS
uniprot map TP73_UNIPROT tp73_locus --projection-id tp73_uniprot_projection
```

## Integrated audit

Run the high-level audit:

```text
uniprot audit-projection tp73_uniprot_projection \
  --ensembl-entry TP73_ENS \
  --report-id tp73_projection_audit
```

Then inspect the stored report:

```text
uniprot audit-show tp73_projection_audit
```

The report includes per-transcript rows with:

- contributing exon spans in transcript order
- exon nucleotide sum
- untranslated 5' nt
- untranslated 3' nt
- translated nt
- divisibility-by-3 result
- expected amino-acid count
- UniProt amino-acid count
- direct-compare vs global-alignment branch
- mismatch reasons
- a local unsent maintainer-email draft for failing isoforms

## Reusable primitive composition

The same TP73 result can be inspected step by step:

```text
uniprot resolve-ensembl-links tp73_uniprot_projection
uniprot transcript-accounting tp73_uniprot_projection
uniprot compare-ensembl-exons tp73_uniprot_projection --ensembl-entry TP73_ENS
uniprot compare-ensembl-peptide tp73_uniprot_projection --ensembl-entry TP73_ENS
```

Those commands expose the same building blocks a calling AI can orchestrate
outside GENtle.

## Direct-vs-composed parity

To verify that the integrated Rust audit matches the public primitive
composition, run:

```text
uniprot audit-parity tp73_uniprot_projection \
  --ensembl-entry TP73_ENS \
  --report-id tp73_projection_audit_parity
```

Then inspect it:

```text
uniprot audit-parity-show tp73_projection_audit_parity
```

The parity report compares:

- transcript row set
- accounting numbers
- divisibility/expected-aa outputs
- mismatch classification
- direct-compare vs alignment path
- failing transcript set used in the email draft

## Export

Reports can be exported as JSON:

```text
uniprot audit-export tp73_projection_audit tp73_projection_audit.json
uniprot audit-parity-export tp73_projection_audit_parity tp73_projection_audit_parity.json
```

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context copied from GENtle Help -> Tutorial -> Copy Feedback Context.

- Tutorial title:
- Tutorial/chapter id:
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
