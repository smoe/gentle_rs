# Expression Draft CLI Fixture Provenance

`protein_expression_handoff_cli.rs` constructs temporary, hand-authored inputs:
a 27-bp synthetic coding sequence (`ATG` + six His codons + Ala + stop), or
the corresponding eight-amino-acid protein, with no biological accession.
It reuses the explicitly synthetic requirements in
`docs/examples/plans/protein_expression_requirements_in_house.json`, changing
permission/allowlist facts only inside the temporary case. An active-catalog
negative control copies the bundled provider record and removes its expression
capability, disables quote support or requires a synthetic absent source field;
it makes no claim about an actual provider's current services.

Recreate with:

```bash
cargo test -q --locked --test protein_expression_handoff_cli
```

Glen/CI should also run the `protein_expression` and `external_service` library
filters with `--no-default-features --locked --offline`, the protocol tests and
the existing tutorial/parity checks. Local Python source/checkout checks do not
execute these Rust contracts or accept an expression/provider workflow.

The real CLI and shared shell read LF/CRLF requirements in paths containing
spaces/apostrophes. Cases include neutral/explicit selection, DNA/protein,
in-house permission, exclusions, missing product/permission, partial answers
and unknown providers, unusable products, unavailable quotes and blocked local
preflights under project overlays. HOME/XDG/catalog/temp roots are
child-local; no process-global environment or working directory is changed.
Before/after directory and raw-file snapshots prove no output files or input
rewrites for both the handoff and explicitly invoked local preflight. Adapter
reports must match except their generation timestamp.
Eligible, blocked and quote-disabled cases compare text action lines with the
structured report, while invoking the same canonical preflight separately.
Core unit tests additionally assert `state_changed == false`, state equality,
local-preflight eligibility and no direct submission. No network, credential,
pricing, submission or order is involved. Authored tests are not executed
acceptance until Glen/CI records their exact source and platform.
