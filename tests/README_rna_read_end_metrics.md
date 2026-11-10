# Synthetic RNA Read-End Fixtures

The new `test_rna_read_end_locus.py` uses
`scripts/prepare_exon_locus_rehearsal.py`: three hand-crafted reads, shared-exon
sense/antisense models and independently created A12/Cn/T18 originals. No real
genes, barcodes or sequencing files. Recreate with
`python3 -m unittest tests.test_rna_read_end_locus -v`; its actual CLI test needs
`GENTLE_READ_END_BIN_DIR` and otherwise explicitly skips. It binds fresh SAM,
manifest, annotation, CDS and original-base files before using the shared
operation, both projections and five export formats.
The fast checks exercise fixture generation/binding behavior only, not GENtle
semantics. The actual process test replays both real LF and CRLF input graphs
through GENtle and also exercises a narrow-window JSON-only export; it stays
unrun/skipped without the explicit source-matched binary. No source-text grep
is offered as confirmation/service acceptance.
Inline Rust locus/service/renderer tests add junction preference, mirrored/RC geometry, novel/unresolved
cases, missing/changed annotation, chimeras, counters, request-bound cancellation
and no-clobber. Run externally with the handoff's `read_end_` filters. These are
synthetic contract checks, not model calibration or real-data biological evidence.

`read_ends/locus/review_tests.rs` adds hand-crafted source-review counterexamples:
one/three-base and intronic containment, mirrored chains plus an unrelated model,
deletion-split full junction chains, query-disjoint segment hypotheses, positive
novel compatibility, differing hypotheses/CDS metrics, unavailable units, named
budgets, two distinct samples, run-variant unions and membership. Its scale case
creates 30 models (12 exons, 10 minus strand) and 20 reads (40 deletion-split
blocks each) and prints consumed work; no original bases are supplied in that
geometry-scale fixture. It is not a measured real-read capacity benchmark.
Service tests create declared-absent paths without deleting source fixtures,
retain wrong-hash refusal and check absent originals and JSON-only export.
Recreate externally with `cargo test --locked --lib -j1 read_end_locus_ -- --nocapture --test-threads=1`.

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
