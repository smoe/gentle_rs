# Synthetic PATZ1 locus-composer fixture

## Origin

Every file in this directory is hand-crafted synthetic test data. The files do
not contain experimental CUT&RUN measurements or proprietary Clariom records,
and they must not be cited as biological evidence. They share the artificial
`GRCh38.p14:chr22:31325800..31326039` coordinate frame documented by
`test_files/fixtures/isoform_evidence/patz1/patz1_minus_strand.gb`.

The two probe rows preserve distinct PSR and JUC geometry and carry invented
TAp73alpha-GFP and DNp73beta-GFP effects. The four BED files contain one
invented occupancy interval each and exist only to exercise grouped,
shared-scale locus rendering.

## Deterministic recreation

1. Use the fixed 240 bp synthetic PATZ1 anchor above.
2. Place `PSR_SYN_PATZ1_EXON_A` on genomic bases `31325820..31325859`.
3. Place `JUC_SYN_PATZ1_SKIP` across junction edges
   `31325860..31325979`.
4. Retain the exact effect values and BED intervals committed here. No random
   values, timestamps, network results, or local analysis outputs enter the
   fixture.

## GENtle use

- `docs/examples/workflows/patz1_gene_locus_evidence_offline.json` imports the
  four BED tracks, composes the locus report, and renders the deterministic
  SVG used by the generated tutorial.
- Workflow/example tests exercise the shared `RenderFeatureExpertSvg`
  operation. The SVG must retain separate PSR/JUC markers, both effect
  contrasts, grouped occupancy lanes, and TP73 motif provenance.
- The richer 24-row PATZ1-derived fair-use fixture remains separately
  documented under
  `test_files/fixtures/gene_locus_evidence/patz1_probe_effects/`; this compact
  synthetic fixture is the offline graphical acceptance case.
