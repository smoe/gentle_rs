# Expression Draft CLI Fixture Provenance

`protein_expression_handoff_cli.rs` constructs temporary, hand-authored inputs:
a 27-bp synthetic coding sequence (`ATG` + six His codons + Ala + stop), or
the corresponding eight-amino-acid protein, with no biological accession.
It reuses the explicitly synthetic requirements in
`docs/examples/plans/protein_expression_requirements_in_house.json`, changing
permission/allowlist facts only inside the temporary case. An active-catalog
negative control copies the bundled provider record and removes its expression
capability; it makes no claim about an actual provider's current services.

Recreate with:

```bash
cargo test -q --locked --test protein_expression_handoff_cli
```

The real CLI and shared shell read LF/CRLF requirements in paths containing
spaces/apostrophes. Cases include neutral/explicit selection, DNA/protein,
in-house permission, exclusions, missing product/permission, partial answers
and a capability-removing project overlay. HOME/XDG/catalog/temp roots are
child-local; no process-global environment or working directory is changed.
Before/after directory and raw-file snapshots prove no output files or input
rewrites. Adapter reports must match except their generation timestamp.
Core unit tests additionally assert `state_changed == false`, state equality,
local-preflight eligibility and no direct submission. No network, credential,
pricing, submission or order is involved. Authored tests are not executed
acceptance until Glen/CI records their exact source and platform.
