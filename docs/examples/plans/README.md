# Gibson Planning Examples

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
