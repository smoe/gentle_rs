# Synthetic PrimerBank HTML fixture

`synthetic_primerbank_search_result.html` is a hand-written, deterministic HTML
fixture modeled only on the field labels and table nesting exposed by the
public PrimerBank search form in July 2026. It contains invented identifiers,
sequences, descriptions, and measurements; no PrimerBank database record or
experimental result is copied.

The fixture is used by `src/primerbank.rs` unit tests to verify HTML parsing,
forward/reverse coding-sequence coordinates, multiple pair retention, source
links, and the rule that GENtle does not infer validation from catalog
presence. Recreate it by editing the HTML directly while retaining the two
table classes exercised by the parser: one gene-description table followed by
one or more primer-pair tables, nested inside a non-semantic layout table as on
the public site.
