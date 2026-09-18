# Synthetic TSS Collection Fixture

Origin: hand-crafted, entirely synthetic; no private study data or real genome.
The FASTA is `AACCGTGA` repeated 250 times (2,000 bases), wrapped at 80 bases.
The GTF has four single-exon transcripts, all labelled TOY:

| Transcript | Gene ID | Genomic interval (1-based inclusive) | Strand |
| --- | --- | --- | --- |
| plus_a | plus_gene | 601..1100 | + |
| plus_b | plus_gene | 601..1200 | + |
| plus_c | plus_gene | 901..1250 | + |
| minus_a | minus_gene | 1201..1500 | - |

Recreate by repeating the sequence and writing these four rows as GTF exon
records with source `synthetic`, gene_name `TOY`, and exon_number `1`.
The catalog uses local relative paths only. No network or BLAST evidence is
needed; optional index availability does not establish any scientific result.

Used by `tss_collection_gui_starter`, the TSS tutorial, and deterministic
engine/tutorial tests. Expected 500-upstream/200-downstream windows:
101..801 (+), 401..1101 (+), and 1300..2000 (-), each 701 bp.
The negative-strand window is reverse-complemented: local base 501 maps to
genomic 1500, then increasing local coordinates decrease genomic coordinates.
The two plus transcripts sharing base 601 form one window. The two gene IDs
are never merged; their common label only selects both for this artificial demo.

`inventory.request.json` fixes the preview request. The independent
`tss_collection_gui_oracle` workflow stores the resulting explicit approval and
selected IDs, not a general-purpose approval override. Reproduce/check it with
`cargo test --lib tss_tutorial_starter -- --nocapture`; compare the emitted
request and IDs before updating the oracle after a deliberate fixture or
grouping-contract change. The 2026-09-18 refresh uses resolved gene IDs rather
than display labels in TSS identity; DNA, starts, strands, memberships and
window lengths are unchanged. This is explicit reapproval of this public
synthetic oracle only, not a migration of saved user collections.
`tss_gui_acceptance_starter_and_oracle_are_independent_and_report_bound` checks
the fixed approval, report assertions, persistence, stale rejection, forgetting
and undo, including the sequence IDs named by the GUI state verifiers.
No screenshot or human sign-off is bundled.
