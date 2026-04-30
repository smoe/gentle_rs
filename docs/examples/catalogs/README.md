# Catalog Extension Examples

This directory contains non-executable catalog fragments that demonstrate how
to extend GENtle's catalog-backed knowledge layer outside the source tree.

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
