# Synthetic Read-End Gene-Ranking Tutorial

Origin: Glen / smoe-bot's hand-crafted six-sample teaching generator,
`scripts/build_rna_read_end_tutorial.py`, introduced in the published UTR stack
ending at `3ff288027ba3bede9a3aeaadb951459f4d418832`.

Regenerate `rna_read_end_gene_ranking/` and its two teaching SVGs with:

```sh
python3 scripts/build_rna_read_end_tutorial.py
```

The generator intentionally replaces only its generated fixture directory.
Rows, transcript names and repeated digest characters are synthetic. Named
genes/cell lines illustrate the explicit contrast design; these numbers are not
measurements or biological evidence. The fixture exercises descriptive ranking
and separate library-preparation contrasts, not alignment or UTR projection.
The inline tests in `tests/test_rna_read_end_boundaries.py` independently
exercise those projection boundaries.

The generated report binds exact TSV input bytes. Scoped LF checkout rules and
`scripts.test_tutorial_checkouts` protect both those hashes and regeneration
freshness in real LF/CRLF-policy Git checkouts, with a missing-policy control.
