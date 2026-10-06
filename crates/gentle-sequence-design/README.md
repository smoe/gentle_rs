# Synonymous motif/GC design core

Private, dependency-free Rust prototype, independently versioned at 0.1.0.
Original GENtle implementation, not code copied from DNA Chisel. MSRV: Rust 1.85.
The library has no serialization, GUI, Python, filesystem or network dependency.

`solve` accepts one explicit synthetic DNA/CDS/protein problem plus a complete
codon-to-residue mapping. It returns a deterministic, budget-bound result;
`validate_output` independently checks a chosen output without rerunning search.
`solve_with_strategy` optionally selects `SearchStrategy::ConflictDirected`;
`solve` and `FullEnumeration` keep the original enumeration behavior.
Optional `gc_content` adds inclusive basis-point bounds on the complete declared
CDS, including frozen start/stop, not flanks. Bounds are converted once to integer
GC counts; GC-only requests require explicit bounds. Both solvers and the fresh
output validator use one full motif/GC predicate. `assess_gc` returns integer
input facts and output facts only for present, validated DNA. Absent GC has no
assessment. GC-enabled algorithm identities are distinct; legacy identities stay
unchanged. Reachability extrema can prove impossibility, never feasibility.

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
For a GC conflict it selects unassigned codons capable of changing GC in the
needed direction, but tries all their synonyms, including neutral/worsening
ones. Frames retain branch cursors, not copied whole-CDS branch lists, keeping
stack storage linear in assigned-codon depth. Every candidate and final output
receives full validation; no local-update
shortcuts or evaluation cache are used. A completed conflict search can prove
minimum edits without enumerating unrelated synonyms; budget interruption still
does not. No large-insert efficiency or runtime-speedup claim is made.
Cancellation never publishes a candidate, even when one was found.

The synthetic inline fixtures in `src/tests.rs` document their literal recreation
and use. The four-state coupled-edit oracle explicitly rules out strict greedy
repair. An independently constructed eight-state F/K/F oracle adds degenerate
motifs, strands and protection; a literal CTG/GAA/TTC fixture requires repairing
a newly introduced conflict outside the original match. The repeated-GCT fixture
compares candidate counts/completion only, not measured performance. GC fixtures
add exact count boundaries, empty windows, frozen contributions, a literal
unrestricted eight-state oracle, a 36-state L/R oracle with mixed-direction
synonyms and coupled motif/GC repairs. No fixture
represents PATZ1 or another natural gene.

```sh
cargo test -p gentle-sequence-design
```

Preserved translation does not establish expression, splicing, folding,
regulatory function or experimental suitability. General DNA Chisel
specifications, sliding-window GC/codon adaptation, nonstandard codes and
publication are deferred.
