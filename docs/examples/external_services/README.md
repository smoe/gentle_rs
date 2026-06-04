# External-Service Example Requests

These examples exercise the provider-neutral
`gentle.external_service_request.v1` contract used by `services
project-preflight` and `services project-quote`.

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

Both examples are local, deterministic, and non-mutating. GENtle prepares
review artifacts only; it does not submit orders, scrape WOP, store
credentials, or persist PO/shipping/billing details.

To write a local handoff bundle from either request, use:

```bash
cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/metabion_oligo_single_tube_request.json \
  --output-dir artifacts/external_services/metabion_oligo_demo
```

The output directory is created if needed and the returned quote report records
generated files in `service_ready_bundle.local_files[]`.
