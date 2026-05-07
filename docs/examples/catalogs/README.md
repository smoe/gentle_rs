# Catalog Extension Examples

This directory contains non-executable catalog fragments that demonstrate how
to extend GENtle's catalog-backed knowledge layer outside the source tree.

## External-Service Providers

`external_service_providers_project_overlay.json` demonstrates a project-local
provider override for the config schema
`gentle.external_service_provider_config.v1`. To use the same shape in a real
project, place it under:

```text
.gentle/catalogs/external_service_providers.d/metabion_local_policy.json
```

Provider discovery is overlay-based: built-in assets load first, then system,
user, and project catalogs. Later files override the same `provider` id, so a
lab can adjust WOP/email channels, template field names, validation rules,
delivery defaults, or required follow-up text without changing GENtle source
code.

Validate the active chain with:

```bash
gentle_cli services providers doctor
```

The example keeps the v1 safety posture: it prepares deterministic quote/
handoff artifacts only and does not submit orders, scrape portals, or store
PO/account/shipping/billing details.

## Helper Semantics Vocabulary

`helper_semantics_vocabulary_extension.json` is a project-local vocabulary
overlay example. To use the same shape in a real project, place it under:

```text
.gentle/catalogs/helper_semantics_vocabulary.d/solubility_tags.json
```

Then add helper catalog entries whose `semantics.components[].kind` or other
semantic fields use one of the overlay aliases, for example:

```json
{
  "Project MBP vector": {
    "description": "MBP fusion helper for soluble bacterial expression",
    "sequence_remote": "https://example.invalid/project_mbp.fa.gz",
    "annotations_remote": "https://example.invalid/project_mbp.gb.gz",
    "helper_kind": "plasmid_vector",
    "semantics": {
      "schema": "gentle.helper_semantics.v1",
      "components": [
        {
          "id": "mbp",
          "kind": "MBP tag",
          "label": "MBP"
        }
      ]
    }
  }
}
```

`helpers list --filter mbp` and downstream structured helper-catalog consumers
will then see a normalized `component_kind=solubility_tag` term enriched with
the vocabulary label/description and any vocabulary-provided routine hints.

## Gene-Group Drafts

`neoneurogenesis_metastatic_cancer_neuroimmune_axis_draft.json` is a
review-gated trial output from:

```bash
gentle_cli gene-groups draft \
  --description "Genes and receptor pathways induced or co-opted in metastatic cancers that support tumor innervation, perineural migration, neural niche formation, and nerve-associated immune evasion..." \
  --candidate "NGF=tumors can release nerve growth factor to stimulate nerve growth and invasion-associated neoneurogenesis" \
  --candidate "EFNB1=exosome-packaged EphrinB1 potentiates tumor innervation/axonogenesis in head-and-neck cancer models" \
  --agent-provider Codex \
  --agent-model GPT-5 \
  --output docs/examples/catalogs/neoneurogenesis_metastatic_cancer_neuroimmune_axis_draft.json
```

The file is intentionally still `curation_status=draft`. It demonstrates how an
AI/literature-assisted session can produce a local catalog fragment with
candidate symbols, evidence notes, agent provenance, and a typo/alias bridge
(`neuneurogenesis`) without making those candidates trusted GENtle facts.

The trial candidate set was informed by open review/primary sources on
tumor-induced neoneurogenesis, perineural invasion, tumor exosome-mediated
innervation, beta-adrenergic immune suppression, semaphorin/Plexin-B1 tumor
biology, and the TAp73alpha-HDAC2/REST-GABBR2 melanoma metastasis axis.
Expert review should decide which candidates are retained, split into
subgroups, or rejected before adding this fragment to a project/user catalog.
