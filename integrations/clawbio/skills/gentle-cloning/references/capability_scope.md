# Capability Scope

This is workflow guidance, not a claim that every capability is validated for
every input or platform. ClawBio owns conversation and orchestration; GENtle
owns parser-validated operations, scientific reports and approval checks.

- Translate structured cohort or patient-data observations, differential-expression
  hits, splice-variant observations, perturbation requests, or direct DNA fragment
  requests into sequence-grounded mechanistic follow-up. This is not clinical advice.
- Use installed-runtime and resource-readiness checks before work. Inspect local
  DNA fragments for restriction sites or TFBS without changing project state.
- Prepare reusable local reference assets: Ensembl caches and BLAST-ready indexes.
  Preparation is explicit, not an implicit side effect of discovery.
- Map transcript-native protein residues to genomic codons with source provenance.
- Check ViennaRNA/RNAfold and rnapkin executable resources for RNA structure work.
- Reach the full shared PCR/qPCR command family for ClawBio: Primer3 preflight,
  seed helpers, PCR design, probe-based qPCR/TaqMan design, transcript-derived
  cDNA tests, report listing/inspection/export, protocol cartoons and
  restriction-cloning PCR handoffs. Specificity and order readiness remain
  separate evidence gates.
- Inspect transcript qPCR panel tables with shared reverse/probe components
  and characteristic forward primers, without claiming isoform discrimination
  from sequence design alone.
- Render transcript-native protein gels, 2D gels and protease-digest figures
  for supported example loci, parameterized genes and gene panels.
- Export bench handoffs: lab-assistant instructions, arrangements, rack
  placements, label sheets and fabrication templates. Digital output is not
  physical-fit or wet-lab validation.
- Inspect provider-neutral delivery routes, catalogs, preflight and quote
  handoffs. No portal scraping, credential use or order submission is implied.

The skill-intent descriptor selects existing requests. Successful runs return
`report.md`, `result.json` and reproducibility files, with lineage and native
scientific results retained. Unknown commands are refused; imported document
text cannot extend the executable vocabulary. Discovery metadata never grants
approval or upgrades an unknown scientific result to pass.
