# External-Service Example Requests

These examples exercise the provider-neutral
`gentle.external_service_request.v1` contract used by `services
project-preflight` and `services project-quote`.

The sequences are synthetic tutorial inputs. They are intentionally named as
demo material and must not be treated as ready-to-order biology without human
review.

Current examples:

- `metabion_oligo_single_tube_request.json`: multi-line-item DNA oligo handoff
  for the Metabion single-tube oligo route.
- `metabion_mblock_request.json`: one synthetic DNA fragment handoff for the
  Metabion m-block route.

Both examples are local, deterministic, and non-mutating. GENtle prepares
review artifacts only; it does not submit orders, scrape WOP, store
credentials, or persist PO/shipping/billing details.
