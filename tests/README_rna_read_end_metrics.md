# Synthetic RNA Read-End Fixtures

Origin: Glen's inline synthetic fixtures and the hand-crafted Codex boundary
controls in `test_rna_read_end_boundaries.py`. No public/private sequencing
sample, real GTF or reference FASTA is embedded. Stable names and repeated
digest characters are synthetic identifiers, not authenticated source evidence.

Recreate and exercise the temporary fixtures with:

```sh
python3 -m unittest discover -s tests -p 'test_*rna_read_end*.py' -v
```

These tests exercise the standalone GTF preparer, transcript-oriented projector,
descriptive ranking/contrasts and unchanged comparator. Boundary controls cover
both strands, split/included/separate stops, missing candidate geometry, legacy
models, true anchored zero, censoring and actual Python CLI subprocesses under
LF/CRLF and spaced/apostrophe/Unicode paths. They establish no Rust/GUI/native
Windows or biological acceptance on a macOS run.
