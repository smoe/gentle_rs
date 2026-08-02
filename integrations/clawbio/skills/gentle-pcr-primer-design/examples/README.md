# Primer-Design Skill Examples

These files are deterministic examples for the descriptor-only
`gentle-pcr-primer-design` ClawBio/OpenClaw skill.

## Provenance

- `synthetic_imported_primer_pairs.json` is hand-crafted test data. Its short
  sequences match the existing synthetic cDNA-assay tutorial fixture and do
  not represent commercial, literature, patient, or unpublished laboratory
  primers.
- The conventional-PCR and PATZ1-like requests delegate to existing committed
  GENtle workflows. Those workflows document their own synthetic fixture
  provenance.
- The missing-Primer3 request deliberately names a nonexistent executable to
  demonstrate an actionable dependency failure without installing software.

## Recreation

The request files are ordinary `gentle.clawbio_skill_request.v1` JSON. They can
be recreated by pairing the documented shared shell/workflow commands with the
generic `gentle-cloning` wrapper. No network access or generated timestamp is
required.

Every request enables `claim_attribution_mode: strict` and
`presentation_profile: pcr_primer_design`. The resulting report therefore
retains the raw GENtle payload, prefixes display prose by source, and writes a
machine-readable `claim_ledger.json`.

Every request also carries a static `delegation` identity matching the owning
intent and plan step. The generic wrapper verifies that identity against this
directory's catalog and descriptor and writes the hashes to
`reproducibility/execution_manifest.json`. The imported-primer example uses an
explicit `input_bindings[]` row so both its workflow and synthetic primer batch
are content-addressed before execution.

The PATZ1 transcript-panel request additionally records oligo-dT priming and
the absence of cap-dependent 5-prime-completeness evidence in `input_claims[]`.
Those statements remain caller/fixture assumptions and are rendered as
`[input]`, while GENtle's transcript geometry and oligo-dT reach results remain
`[gentle]`.

The imported-primer batch can be recreated exactly with:

```json
{
  "schema": "gentle.external_primer_pair_batch.v1",
  "batch_id": "synthetic_cdna_assay_primer_pair",
  "pairs": [
    {
      "source_kind": "laboratory",
      "provider": "Synthetic GENtle fixture",
      "catalogue_id": "",
      "source_url": "",
      "claimed_accession": "TX1",
      "aliases": ["synthetic_TX1_pair"],
      "forward_sequence_5_to_3": "AAACCC",
      "reverse_sequence_5_to_3": "CCCAAA",
      "claimed_target": "Synthetic TX1 cDNA tutorial target",
      "validation_claims": [],
      "annotations": {
        "fixture_status": "synthetic_only"
      }
    }
  ]
}
```

## Use In GENtle

- The specialized skill tests validate metadata, intent ownership, delegation,
  demo gating, and evidence-state wording.
- `imported_primer_review_offline.workflow.json` exercises the existing
  `ImportExternalPrimerPairs` operation against
  `docs/examples/assets/cdna_assay_demo.gb`.
- The examples never establish genomic specificity or order readiness.
- The default experimental-assay report asks GENtle for a shared-condition
  virtual gel with one primer pair per lane; it is predicted visualization, not
  a laboratory gel.
