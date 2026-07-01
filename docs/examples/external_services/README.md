# External-Service Example Requests

These examples exercise the provider-neutral
`gentle.external_service_request.v1` contract used by `services
project-preflight` and `services project-quote`, plus the
`gentle.external_service_delivery_route_request.v1` contract used by
`services delivery-route`.

The sequences are synthetic tutorial inputs. They are intentionally named as
demo material and must not be treated as ready-to-order biology without human
review.

Current examples:

- `geneart_cloned_gene_request.json`: synthetic cloned-gene/plasmid
  construction quote handoff for GeneArt.
- `geneart_protein_expression_request.json`: synthetic protein-expression
  quote handoff for GeneArt, including a `return_spec` that asks for
  amino-acid sequence, codon-targeted DNA, quote metadata, and a handoff
  bundle.
- `metabion_oligo_single_tube_request.json`: multi-line-item DNA oligo handoff
  for the Metabion single-tube oligo route.
- `metabion_mblock_request.json`: one synthetic DNA fragment handoff for the
  Metabion m-block route.
- `sequence_delivery_short_oligo_route_request.json`: generic "deliver this
  sequence" routing input that resolves to a Metabion single-tube oligo
  candidate.
- `sequence_delivery_fragment_route_request.json`: generic delivery routing
  input that resolves to a Metabion m-block fragment candidate.
- `sequence_delivery_cloned_gene_route_request.json`: generic delivery routing
  input with vector/construct context that resolves to a GeneArt cloned-gene
  candidate.
- `sequence_delivery_protein_route_request.json`: generic delivery routing
  input that resolves to a GeneArt protein-expression candidate.

All examples are local, deterministic, and non-mutating. GENtle prepares
review artifacts only; it does not submit orders, scrape WOP, store
credentials, or persist PO/shipping/billing details.

To classify a generic sequence-delivery request before quote handoff, use:

```bash
cargo run --quiet --bin gentle_cli -- services delivery-route \
  @docs/examples/external_services/sequence_delivery_short_oligo_route_request.json
```

To write a local handoff bundle from either request, use:

```bash
cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/metabion_oligo_single_tube_request.json \
  --output-dir artifacts/external_services/metabion_oligo_demo
```

The output directory is created if needed and the returned quote report records
generated files in `service_ready_bundle.local_files[]`.
