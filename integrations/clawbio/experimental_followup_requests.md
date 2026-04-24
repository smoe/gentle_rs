# Requests for ClawBio Experimental Follow-Up and Perturbation Planning

Audience: ClawBio/OpenClaw developers and skill-orchestration maintainers.

Status: proposed integration behavior. GENtle already exposes useful
sequence-context, dbSNP, promoter, reporter, TFBS, splicing, isoform,
routine-catalog, and planning-estimate routes. This note describes how ClawBio
should orchestrate them for natural-language questions such as "determine the
effect of this SNP", "what should we do with this differentially expressed
gene", or "how do we characterize this splice variant".

## Goal

ClawBio should treat biological observations as prompts for a structured
experimental follow-up menu, not as isolated annotation tasks.

The observation may be:

- a SNP, small variant, VCF row, or rsID,
- a differentially expressed gene,
- a splice isoform or exon-usage change,
- a protein/domain change,
- a pathway or gene-set hit,
- a direct user request to overexpress, knock down, inhibit, or disrupt a
  gene.

For each case, ClawBio should:

1. normalize the biological target and evidence type,
2. optionally run external predictors or knowledge sources,
3. ask GENtle for deterministic sequence-grounded context and assay-planning
   artifacts,
4. compare plausible perturbation and readout families,
5. include practical planning information such as time, cost, local material
   availability, and missing-resource procurement delays.

GENtle should be used as the deterministic sequence, construct, routine, and
planning engine. ClawBio and its language-model layer should decide how to
weigh predictor evidence, literature, cell models, disease context, and
experimental strategy.

## Requested ClawBio Intent Routing

Route user prompts to the `gentle-cloning` skill when they contain an rsID,
genomic coordinate, gene symbol, transcript/isoform id, VCF-like variant, or
pasted sequence plus intent language such as:

- "determine the effect of this SNP"
- "what should we do with this differentially expressed gene"
- "characterize this splice variant"
- "overexpress this gene"
- "knock down this gene"
- "make a loss-of-function construct"
- "design an assay for this regulatory variant"
- "test this promoter/enhancer/UTR"
- "compare plasmid versus viral delivery"
- "which follow-up is cheapest or fastest"

The skill route should be used even if ClawBio also invokes external tools.
The intended split is:

- ClawBio/external predictors: biological likelihood, disease/literature
  context, omics evidence, predictor evidence, and strategic prioritization.
- GENtle: sequence context, local feature geometry, construct/routine
  planning, cloning artifacts, readout suggestions, and graphics.

## Requested Evidence Classes

ClawBio should map each prompt into one or more evidence classes before choosing
follow-up actions:

- `variant_context`: SNP, indel, VCF row, or rsID.
- `differential_expression`: upregulated/downregulated gene from RNA-seq,
  proteomics, or user-supplied table.
- `splice_isoform_context`: isoform switch, exon inclusion/skipping, junction
  usage, or transcript-specific protein change.
- `protein_context`: coding change, domain disruption, tag/fusion/isoform
  characterization, protein size/pI expectations.
- `regulatory_context`: promoter, enhancer, TFBS, UTR, eQTL, or chromatin
  evidence.
- `perturbation_intent`: explicit overexpression, knockdown, knockout,
  CRISPRi/CRISPRa, antisense, siRNA/shRNA, base/prime editing, reporter, or
  rescue experiment.

Evidence classes should remain additive. For example, a splice variant in a
differentially expressed gene may produce both isoform-characterization and
expression-perturbation follow-up options.

## Requested External Predictor Orchestration

ClawBio should be free to run or suggest predictors before or alongside
GENtle. The result should be summarized as evidence, not treated as final
truth.

Useful predictor/source classes:

- Variant consequence and transcript context: VEP, Ensembl lookup,
  SnpEff-like annotation.
- Coding effect: SIFT, PolyPhen, AlphaMissense-like scores when available.
- Splicing effect: SpliceAI or equivalent splice-impact predictors.
- Regulatory evidence: RegulomeDB-like annotations, GTEx/eQTL links, ENCODE or
  local epigenomic context when available.
- Expression evidence: differential-expression tables, fold change, adjusted
  p-value, expression baseline, tissue/cell-type relevance.
- Isoform evidence: transcript abundance, junction support, ORF/protein-domain
  consequences, isoform-specific peptide/protein detectability.
- Population/clinical context: gnomAD frequency, ClinVar, pharmacogenomic
  resources, and relevant literature.

ClawBio should pass the predictor summary into the final response separately
from GENtle's deterministic outputs. This keeps uncertainty visible.

## Requested GENtle Calls Available Now

Current GENtle/ClawBio scaffold examples already cover useful building blocks:

- `examples/request_skill_info.json`
- `examples/request_dbsnp_fetch_rs9923231.json`
- `examples/request_inspect_sequence_context_rs9923231_vkorc1.json`
- `examples/request_export_sequence_context_bundle_rs9923231_vkorc1.json`
- `examples/request_workflow_vkorc1_context_svg_auto_prepare.json`
- `examples/request_workflow_vkorc1_planning.json`
- `examples/request_render_svg_rs9923231_vkorc1_linear.json`
- `examples/request_export_bed_rs9923231_vkorc1_context_features.json`
- `examples/request_seed_qpcr_tp53_splicing.json`
- `examples/request_workflow_tp53_splicing_expert_svg.json`
- `examples/request_workflow_tp73_isoform_protein_gel.json`
- `examples/request_workflow_tp73_isoform_protein_2d_gel.json`
- `examples/request_genomes_extract_gene_tp53_auto_prepare.json`
- `examples/request_protocol_cartoon_gibson_svg.json`
- `examples/request_protocol_cartoon_qpcr_svg.json`

For a new rsID, the first generic path should adapt the dbSNP fetch request to
the user's variant, then follow with sequence-context or graphics requests.

For a gene or isoform named from differential expression or splicing evidence,
the first generic path should adapt a genome/gene extraction or Ensembl-gene
fetch request, then follow with routine planning, qPCR/splicing/isoform, or
protein-gel/2D-gel requests as appropriate.

## Requested Perturbation and Readout Menu

ClawBio should map evidence classes into intervention/readout families. GENtle
can provide sequence-grounded artifacts for many of these routes, but the final
ranking should include external predictor evidence, user goals, local lab
constraints, and planning estimates.

Expression increase or overexpression:

- transient plasmid overexpression,
- stable plasmid integration when appropriate,
- adenoviral delivery for efficient transient expression in difficult cells,
- lentiviral delivery for stable expression or hard-to-transfect models,
- isoform-specific ORF design,
- tagged or untagged expression constructs,
- rescue experiment after knockdown or knockout,
- protein gel / 2D-gel expectations for isoform or tag characterization.

Expression decrease or knockdown:

- antisense oligo or morpholino-style knockdown where applicable,
- siRNA for transient knockdown,
- shRNA or lentiviral shRNA for longer-term knockdown,
- CRISPRi for transcriptional repression without cutting,
- promoter or enhancer perturbation when regulatory control is the question,
- qPCR/RT-PCR primer planning to measure knockdown or isoform effects.

Genomic disruption or editing:

- CRISPR knockout by frameshift indel,
- deletion of promoter/enhancer/UTR elements,
- base or prime editing for precise nucleotide changes,
- CRISPRa for endogenous activation,
- endogenous editing plus rescue as a stronger causal design when feasible.

Reporter assays:

- promoter or enhancer luciferase reporter,
- allele-paired promoter fragment design,
- UTR or translation-efficiency reporter,
- minigene splicing reporter,
- TFBS delta scan and local sequence graphics.

Splicing and isoform characterization:

- RT-PCR or qPCR isoform assay,
- minigene splicing reporter,
- RNA-read mapping and Splicing Expert follow-up,
- isoform architecture and protein-domain inspection,
- protein molecular-weight or 2D-gel figure for isoform separation.

Coding or protein-level changes:

- allele-paired overexpression construct,
- protein-domain or isoform-context report,
- tag/fusion construct planning,
- expression vector and host/delivery comparison,
- endogenous editing when exogenous overexpression may distort biology.

Deep intergenic or weakly contextualized signals:

- nearest-gene and enhancer/context report,
- TFBS or epigenomic overlap scan when annotations are available,
- cautious "insufficient local evidence" wording unless external evidence
  supports a specific target,
- prioritize external predictor/literature/context retrieval before expensive
  wet-lab planning.

## Requested Economics and Practical Planning Layer

ClawBio should make GENtle's routine/planning layer visible whenever the user
asks what to do next, which method to choose, or how expensive/slow a plan may
be.

GENtle already exposes:

- `routines list [--query TEXT] [--seq-id SEQ_ID]`
- `routines explain ROUTINE_ID [--seq-id SEQ_ID]`
- `routines compare ROUTINE_A ROUTINE_B [--seq-id SEQ_ID]`
- `planning profile show`
- `planning objective show`
- `planning objective set JSON_OR_@FILE`

Routine comparison payloads can include:

- `estimated_time_hours`,
- `estimated_cost`,
- `local_fit_score`,
- `composite_meta_score`,
- guardrail pass/fail information,
- missing local material classes,
- procurement delays using the deterministic business-day heuristic.

ClawBio should surface this as practical planning evidence, not as exact
pricing. Good wording:

- "GENtle estimates this as faster with current local materials."
- "This route is likely delayed by missing material procurement."
- "This comparison uses the current planning profile and should be adjusted for
  your lab's actual prices and lead times."

ClawBio should ask for or infer practical constraints when the choice is not
obvious:

- fastest path,
- lowest-cost path,
- highest biological fidelity,
- easiest first-pass screen,
- stable versus transient perturbation,
- native-context versus exogenous expression,
- available cell model and delivery method,
- local vector/helper availability.

## Requested Presentation Behavior

The first response should not overclaim causality or feasibility. Prefer this
structure:

1. "Here is what the observation appears to represent."
2. "Here is what external predictors or omics evidence say, if available."
3. "Here are GENtle-supported follow-up families."
4. "Here are practical tradeoffs, including delivery, time, cost, and local
   materials where known."
5. "Here is the most direct next command/artifact to run."

Use cautious labels:

- "candidate follow-up",
- "assay family",
- "supports testing",
- "sequence-grounded context",
- "predictor evidence",
- "planning estimate".

Avoid phrasing such as:

- "this SNP causes...",
- "this gene is proven to...",
- "GENtle proves...",
- "clinically actionable..." unless a dedicated clinical source and workflow
  explicitly support that framing.

## Requested Structured Output Shape

When ClawBio asks GENtle or composes the final response, preserve a
machine-readable follow-up menu. A useful ClawBio-side normalized shape would
be:

```json
{
  "schema": "clawbio.experimental_followup.v1",
  "observation": {
    "kind": "variant_context",
    "id": "rs9923231",
    "target_gene": "VKORC1",
    "genome": "GRCh38"
  },
  "evidence_classes": [
    "variant_context",
    "regulatory_context"
  ],
  "external_evidence_summary": [],
  "gentle_artifacts": [],
  "followup_candidates": [
    {
      "family": "promoter_luciferase_reporter",
      "intervention_direction": "reporter_readout",
      "priority": "high",
      "rationale": "Variant is promoter-proximal to VKORC1.",
      "caveats": [
        "Reporter assays test regulatory potential outside native chromatin context."
      ],
      "suggested_gentle_request": "examples/request_workflow_vkorc1_planning.json",
      "planning_axes": [
        "estimated_time_hours",
        "estimated_cost",
        "local_fit_score"
      ]
    }
  ]
}
```

This does not need to become a GENtle engine schema immediately. It is a
ClawBio orchestration shape that can later align with a dedicated GENtle
experimental-follow-up command if that proves useful.

## Requested User-Confirmation Rules

ClawBio should ask before launching long-running, network-heavy, or materially
planning-heavy actions:

- preparing a full reference genome,
- running live external predictors,
- downloading large reference assets,
- executing a multi-step workflow that writes many artifacts,
- selecting a viral delivery or genome-editing strategy as the primary plan,
- changing the planning profile or objective.

Short local inspections can run directly when the user has already asked for a
follow-up and the skill runtime is configured:

- `skill-info`,
- `services status`,
- prepared-reference status,
- already-local sequence context inspection,
- small stateless sequence scans,
- routine catalog search and explanation.

## Minimal Acceptance Tests for ClawBio

ClawBio-side tests should cover:

1. A prompt like "determine the effect of rs9923231" routes to
   `gentle-cloning`.
2. A prompt like "what should we do with this differentially expressed gene"
   routes to `gentle-cloning` and produces expression-perturbation candidates.
3. A prompt like "characterize this splice variant" produces splicing and
   isoform-characterization candidates.
4. The planner can run or propose `request_skill_info.json` before deeper work.
5. The planner can adapt a dbSNP fetch request for a different rsID.
6. The final response separates external predictor evidence from GENtle
   sequence-grounded assay suggestions.
7. The final response includes at least one concrete follow-up command or
   request object.
8. Routine/planning outputs are surfaced when the user asks for fastest,
   cheapest, easiest, or best-fit follow-up.
9. Long-running reference preparation and genome-editing/viral-delivery
   strategy choices are confirmation-gated.

