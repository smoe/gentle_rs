Synthetic promoter cohort and ortholog fixtures for GENtle tutorials.

Origin: hand-crafted, deterministic test/tutorial data created for the offline
promoter cohort tutorial. The sequences are artificial 700 bp chromosomes with
short motif-rich or motif-poor promoter windows around synthetic TP73/E2F1/PATZ1
labels. They are not derived from a real genome and must not be interpreted as
biological evidence.

Recreation:

1. Start each toy chromosome as 700 `A` bases.
2. Insert the motif-rich segment
   `GGGGCGGGGTTTGGGGCGGGGTTTGGGGCGGGGTTTGGGGCGGGG` at 1-based positions
   111 and 311 for `human_toy.fa`.
3. Insert the alternating motif-poor segment
   `ATATATATATATATATATATATATATATATATATATATATATATATATATA` at 1-based
   position 511 for `human_toy.fa` and at 111 for `mouse_toy.fa`.
4. Use the same motif-rich segment at 111 for `rat_toy.fa`; extra E2F1 rows
   keep the ortholog resource a small multi-gene fixture.
5. GTF coordinates place plus-strand transcripts at TSS positions 151, 351, and
   551 so `upstream_bp=40` / `downstream_bp=10` resolves compact promoter
   windows.

Used by:

- `docs/examples/workflows/promoter_gene_set_ortholog_cohort_offline.json`
- `docs/tutorial/sources/08-07_promoter_gene_set_ortholog_cohort_offline.json`

