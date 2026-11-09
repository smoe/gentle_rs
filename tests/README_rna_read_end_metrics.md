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

`test_rank_rna_read_end_genes.py` additionally creates inline synthetic detailed
rows for status/value validation under LF/CRLF and quoted/Unicode paths. Both
public ranking writers must refuse malformed rows before creating outputs,
even when gene assignment is unavailable/ambiguous. Literal controls retain
observed zero, a zero censored bound and blank missing values; censored bounds
never enter the median. No external data or persistent fixture is needed.

`test_read_end_rehearsal.py` uses `scripts/prepare_read_end_rehearsal.py` to
create wholly synthetic original sequences (12 A, C spacer, 18 T), plus/minus
transcript annotations and forward/reverse SAM representations. The original
strings are generated independently of SAM reverse complementation. Its fast
test checks raw hashes, explicit unit bindings and LF/CRLF semantics; two real
CLI process tests require `GENTLE_READ_END_BIN_DIR` and otherwise skip explicitly.
The oriented process case covers both transcript strands, softclips, both
directions, JSON/TSV/SVG/PNG/PDF, no-clobber and absence of a project file.
Inline Rust fixtures in `read_ends/oriented_tests.rs`, `pileup_tests.rs`, the
service tests and renderer tests use the same explicitly synthetic boundaries,
adding internal runs, hard clips, alternatives, counters, pixel caps and dense
labels. Recreate externally with the `read_end_` Rust filters in Glen's handoff;
Python/source tests do not execute Rust or certify native GUI/biology.

`test_rna_read_end_core_contract.py` additionally checks the hand-crafted
31-case corpus in `test_files/rna_read_end_core/cases.json` against the unchanged
Python projector under LF and CRLF. Its nearby README documents origin and
recreation. The same literal inputs/outcomes are consumed by the new Rust
`read_ends::tests`; running the Python half does not execute or accept Rust.
