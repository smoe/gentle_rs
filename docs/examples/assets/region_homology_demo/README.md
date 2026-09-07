# Synthetic Region-Homology Tutorial Genomes

These FASTA files are hand-crafted, deterministic tutorial data. They are not
derived from a real organism, commercial sequence, or experimental result.

- `query_toy.fa` contains one 80 bp query region and a second copy of one
  24 bp subregion elsewhere on the same synthetic chromosome.
- `ortholog_toy.fa` contains two query-like blocks separated by a one-base
  target insertion and substitutions. Its expected locus is declared by the
  tutorial request; GENtle does not infer orthology from the alignment.
- `other_toy.fa` contains an unassigned cross-species-like match.
- `genomes.json` is the local genome catalog consumed by `PrepareGenome` and
  `ScreenGenomicRegionHomology`.
- The workflow also saves a hand-crafted occupancy-purpose span inside the
  query so the GUI module-assessment trace can be exercised without claiming
  that an experiment measured occupancy.
- `*_toy.gtf` files provide one hand-crafted placeholder gene annotation per
  synthetic genome so the ordinary genome-preparation contract is exercised.

To recreate the files, use the literal FASTA records and catalog object tracked
in this directory. The sequences intentionally remain small enough to inspect
by eye. They are used by
`docs/examples/workflows/region_homology_promoter_modules_offline.json` and the
generated conservation tutorial. Building their BLAST indexes requires the
ordinary local BLAST+ tools; no network access is used.

The expected-locus declaration is synthetic evidence for contract testing. It
does not establish evolutionary orthology or regulatory function.
