# Synthetic Region-Homology Tutorial Genomes

These FASTA files are hand-crafted, deterministic tutorial data. They are not
derived from a real organism, commercial sequence, or experimental result.

- `query_toy.fa` contains one 80 bp query region and a second copy of one
  24 bp subregion elsewhere on the same synthetic chromosome.
- `ortholog_toy.fa` contains a 31-base prefix, then query bases 0..30 (with
  one extra A after query base 19), the replacement `TTCCGGAAACCC` for query
  bases 30..42, and query bases 42..80, followed by a synthetic suffix.
  This deliberately creates two exact blocks with a measurable shared-locus
  gap; the insertion is inside the left block. Its expected locus is declared by the
  tutorial request; GENtle does not infer orthology from the alignment.
- `other_toy.fa` contains an unassigned cross-species-like match.
- `genomes.json` is the local genome catalog consumed by `PrepareGenome` and
  `ScreenGenomicRegionHomology`.
- The workflow also saves a hand-crafted occupancy-purpose span inside the
  query so the GUI module-assessment trace can be exercised without claiming
  that an experiment measured occupancy.
- `*_toy.gtf` files provide one hand-crafted placeholder gene annotation per
  synthetic genome so the ordinary genome-preparation contract is exercised.
- `module_cases.json` supplies independent synthetic evidence spans and
  expected standalone, paired, insufficient, and strict-repetition decisions.
  The 1% repetition policy is an intentionally strict contrasting scenario,
  not a recommended biological cutoff. No expected verdict is written into a
  GENtle report by the runner.

To recreate the files, use the literal FASTA records and catalog object tracked
in this directory. The sequences intentionally remain small enough to inspect
by eye. They are used by
`docs/examples/workflows/region_homology_promoter_modules_offline.json` and the
generated conservation tutorial. Building their BLAST indexes requires the
ordinary local BLAST+ tools; no network access is used.

The regression test
`workflow_examples_region_homology_uses_real_isolated_blast_indexes` executes
the same committed workflow with scoped native BLAST tools and checks its
query-width alignment rows and explicitly declared orthology evidence.
It skips explicitly if any of `makeblastdb`, `blastdbcmd`, or `blastn` is
missing; installed-but-broken tools or failed workflows remain errors. Nothing
installs BLAST automatically. A separate, real-sequence smoke fixture lives in
[`test_files/fixtures/blast_tp73_isoforms`](../../../../test_files/fixtures/blast_tp73_isoforms/README.md).

Complete the walkthrough, including SVG and assessed module reports, with:

```sh
python3 docs/examples/run_region_homology_tutorial.py --gentle target/debug/gentle_cli --output /tmp/gentle-conservation-tutorial
```

Use an empty output directory. The runner invokes only the GENtle CLI, retains
the exact screen and assessment/render workflows, checks all four decisions,
query width/insertion provenance, and source-project preservation, and writes
`tutorial_receipt.json` with command/output/binary hashes. The same reports can
be inspected in the GUI. This is a headless executable acceptance path; it does
not claim a completed Xvfb GUI acceptance run. `tutorial-check` checks the
starter/screen workflow when BLAST+ is available, separately from this full
acceptance command. Documentation generation never runs this optional native
workflow, keeping committed pages identical with and without BLAST+.
The assessment/render operations consume their complete report inputs in an
empty scratch project. They do not reopen the source project; this also avoids
the CLI's ordinary save round-trip reordering cached restriction-site groups.

The expected-locus declaration is synthetic evidence for contract testing. It
does not establish evolutionary orthology or regulatory function.
