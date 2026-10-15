# Gibson Planning Examples

The companion `../workflows/gene_set_storage_starter.json` is hand-crafted:
two literal eight-base source records and unannotated logical GENEA/GENEB
labels. Tutorial 03.10 consumes it; no accession, genome annotation or physical
stock is represented. Recreate it from the workflow literals. Its normalized
stable IDs are `symbol:genea`/`symbol:geneb`, not the display labels.

## Synthetic Protein-Expression Review Inputs

`protein_expression_requirements_in_house.json` and
`protein_expression_requirements_review_only.json` are hand-crafted tutorial
09.03 inputs, not measured requirements or provider approval. Recreate them
from their literal fields. The latter permits a **local** explicit GeneArt
review projection for the synthetic starter only; it grants no network,
commercial or real-project consent. All eight reviewed topics stay explicit.
`../workflows/protein_expression_review_starter.json` creates the literal
27-bp `ATG` + six His codons + Ala + stop product. No accession or natural gene
is represented. Consumed by tutorial 09.03 and its source/MCP parity tests;
the existing CLI integration matrix independently constructs the same shape.

This folder contains canonical example payloads for
`gentle.gibson_assembly_plan.v1`.

Current examples:

- `gibson_destination_first_single_insert.json`
- `gibson_destination_first_multi_insert.json`

Current status:

- `gibson_destination_first_single_insert.json`
  - matches the engine-consumed single-insert Gibson preview model
  - becomes executable once the referenced `seq_id` values exist in the active
    project state
- `gibson_destination_first_multi_insert.json`
  - remains a draft design-resource example for the future multi-fragment
    expansion

Purpose:

- document the destination-first Gibson planning model
- make the input-vs-derived boundary explicit
- document ordering and internal-junction adaptation policies
- document one-sided versus split overlap selection at each junction
- capture advisory setup guidance separately from hard assembly-junction logic
- model Gibson primers explicitly as 5' overlap plus 3' gene-specific segments
- provide canonical payloads for routine/cartoon/primer-design work
