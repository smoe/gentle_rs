# Synonymous motif removal core

Private, dependency-free Rust prototype, independently versioned at 0.1.0.
Original GENtle implementation, not code copied from DNA Chisel. MSRV: Rust 1.85.
The library has no serialization, GUI, Python, filesystem or network dependency.

`solve` accepts one explicit synthetic DNA/CDS/protein problem plus a complete
codon-to-residue mapping. It returns a deterministic, budget-bound result;
`validate_output` independently checks a chosen output without rerunning search.
`solve_with_strategy` optionally selects `SearchStrategy::ConflictDirected`;
`solve` and `FullEnumeration` keep the original enumeration behavior.
Code 1 only, literal frozen ATG and terminal stop, uppercase A/C/G/T, one forward
contiguous CDS, frozen flanks/protected bases, finite IUPAC patterns. Coordinates
are zero-based half-open. Reverse matches use reverse-complement patterns;
palindromes requested on both strands produce one `Both` match. All overlapping
windows and flank/CDS boundary windows are evaluated.

Limits: 12,000 bases, 16 motifs, 32 bases/motif, 256 protected intervals,
100,000 requested evaluations, 50,000,000 conservative full-evaluation work
units, 4,096 reported matches. The work bound can reduce the effective candidate
budget. Search-space overflow is explicitly unknown, not infeasible. Default
search is full candidate enumeration, initially violating codons first.
Explicit conflict search instead branches over synonymous changes in the first
current violating window and backtracks using a bounded, nonrecursive stack.
Every candidate and final output receives full validation; no local-update
shortcuts or evaluation cache are used. A completed conflict search can prove
minimum edits without enumerating unrelated synonyms; budget interruption still
does not. No large-insert efficiency or runtime-speedup claim is made.
Cancellation never publishes a candidate, even when one was found.

The synthetic inline fixtures in `src/tests.rs` document their literal recreation
and use. The four-state coupled-edit oracle explicitly rules out strict greedy
repair. An independently constructed eight-state F/K/F oracle adds degenerate
motifs, strands and protection; a literal CTG/GAA/TTC fixture requires repairing
a newly introduced conflict outside the original match. The repeated-GCT fixture
compares candidate counts/completion only, not measured performance. No fixture
represents PATZ1 or another natural gene.

```sh
cargo test -p gentle-sequence-design
```

Preserved translation does not establish expression, splicing, folding,
regulatory function or experimental suitability. General DNA Chisel
specifications, GC/codon adaptation, nonstandard codes and publication are deferred.
