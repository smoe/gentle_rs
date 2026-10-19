# Protein Expression Handoff: Manual macOS Captures

These four untouched window-only PNGs illustrate tutorial 09.03 using the
owner-rebuilt GUI reporting revision
`9f01098b749c79fb5c9eccfcfb184bb085c63357` on 2026-10-18. The exact binary,
input and image hashes, capture times, window geometry and observed states
are in [evidence.json](evidence.json). Source cleanliness at build time was
not independently attested. No image was re-encoded, annotated or simulated;
the native window shadow remains part of the raw capture.

## Origin And Reproduction

The product is hand-crafted synthetic DNA, not a natural gene:
`ATGCACCACCACCACCACCACGCTTAA` (27 nt; M + six His + Ala + stop).
The retained [synthetic.gentle.json](synthetic.gentle.json) is the exact starter
project and remained byte-identical after inspection. It was produced by the
already-built CLI at `1c5be3175323eec41359b980b5a41e22c7edaec1` using
`docs/examples/workflows/protein_expression_review_starter.json`. That older
CLI prepared input only; it is not evidence of same-revision adapter parity.
Replaying the workflow reproduces its scientific content, not necessarily
the byte order of serialized computed-feature caches.

1. Start a supplied GUI binary with this project and fresh HOME/XDG/cache/tmp
   and catalog roots. The session used a temporary native app bundle executing
   the owner-built binary with an explicit `env -i` environment, without
   inherited provider credentials. No network namespace was enforced.
2. Open **Services > Protein Expression Handoff...**, leave every field blank,
   click **Inspect**, and retain `missing_product.raw.png` after completion.
3. Explicitly select `expression_review_product`. Paste the complete
   `protein_expression_requirements_in_house.json` fixture from
   `docs/examples/plans/`, leave provider blank, inspect, and expand CDS
   evidence for `neutral_product.raw.png`.
4. Enter `geneart`, inspect again, and expand binding reasons for
   `provider_withheld.raw.png`: `outsourcing_disallowed` is visible.
5. Replace requirements with `protein_expression_requirements_review_only.json`,
   inspect, and expand binding details for `provider_projected.raw.png`.
   The request draft remains local. Do not execute suggested commands.

The JSON was reformatted before pasting. The receipt retains observed source
file hashes and canonical JSON hashes; source verification uses the latter
so a CRLF checkout does not masquerade as changed requirements.

Capture the inspector window alone with macOS `screencapture -l WINDOW_ID`
(the screenshot skill's `take_screenshot.py` wrapper was used here). Inspect
the resulting PNG, not just the live preview: an initially undersized final
capture was discarded and retaken. Use new filenames/receipts for a new run.

## Evidence Limits

This is `captured_not_accepted`, not an automated or scientific acceptance
verdict. Semantic snapshots and copied report JSON were not retained. Native
accessibility exposed the inspector window but not its internal controls;
ordinary GUI coordinates were used after visual inspection. The new
`gui-test-support` locators and MCP/inner-agent follow-up are not in this
binary and still need Glen/CI execution at the integrated revision.

The images demonstrate visible product review, explicit provider binding and
withholding. They do not certify four-interface runtime parity, live-model
behavior, provider readiness, yield or permission to order. No model/provider
was contacted, no suggested service command was run, and only the isolated
synthetic test instance was closed. Existing screenshot-policy gates remain
unchanged.
