# Plan a High-Yield Protein-Expression Handoff

> Type: `CLI / agent walkthrough`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but every operational step routes
> through the shared `planning ...` and `services ...` shell/CLI contracts.

This tutorial shows how to turn the deliberately vague request "give me the
maximal amount of protein" into a reviewable GENtle planning report.

The key lesson is that "maximal protein" is not one biological objective.
Maximum total expression, soluble expression, active protein, purified protein,
secreted protein, and membrane-localized protein can point to different hosts,
vectors, tags, induction policies, purification endpoints, and providers.

GENtle therefore does not answer the request by choosing the strongest promoter
or one default expression host. It emits a read-only handoff report that makes
the missing decisions explicit.

## Safety Boundary

The `planning protein-expression-handoff` route does not:

- create constructs,
- optimize codons,
- select a final vector,
- submit a provider request,
- upload data to GeneArt or any other vendor,
- claim that a sequence is expression-ready.

It is a planning bridge for a human or agent who needs to decide what should
happen next.

## Step 1: Run The Minimal Handoff

From the repository root:

```bash
cargo run --quiet --bin gentle_cli -- planning protein-expression-handoff \
  --objective '{"schema":"gentle.planning_objective.v1","biological_intent":"protein_expression_max_yield"}' \
  --format json
```

Expected outcome:

- `schema` is `gentle.protein_expression_handoff.v1`
- `biological_intent` is `protein_expression_max_yield`
- `status` is `needs_product_definition`
- `product_definition.sequence_present` is `false`
- `product_readiness.status` is `needs_product_definition`
- the report contains `host_chassis_candidates[]`,
  `vector_route_candidates[]`, `missing_questions[]`,
  `service_handoff_candidates[]`, `warnings[]`,
  `suggested_next_actions[]`, and defaulted sequence-readiness records

This is the safest first response when no product sequence has been loaded.
GENtle can still explain the choices, but it refuses to pretend that the
product is already defined.

## Step 2: Read The Missing Questions First

Start with `missing_questions[]`. With no usable product context, the first
questions should ask which coding sequence, ORF, CDS annotation, or target
protein boundaries define the product.

Once a CDS or protein context is inferable, a high-yield protein plan should
answer:

- which yield metric matters: total, soluble, active, purified, secreted, or
  membrane-localized protein
- which purification endpoint matters: crude lysate, affinity capture, polished
  purity, activity assay, buffer, or delivery format
- which tag strategy is preferred, including N- vs C-terminal position,
  cleavage, retention/removal, and whether annotated tags or signal peptides
  should be preserved
- which chassis is acceptable: bacterial, yeast, insect, mammalian, cell-free,
  or provider-managed expression
- whether the protein needs disulfides, glycosylation, cofactors, secretion,
  membrane insertion, low-temperature expression, chaperones, or solubility
  tags
- whether toxicity is expected and what induction/repression policy is needed
- what scale, tag, purity, buffer, and delivery endpoint define success

Do not proceed to cloning or vendor handoff until those answers are explicit
enough for the biology.

## Step 3: Interpret Chassis And Route Candidates

The first V1 report ranks the following chassis candidates conservatively:

- `e_coli`
- `yeast`
- `hek293`
- `insect_baculovirus`
- `cell_free`

These are review candidates, not final recommendations. `E. coli` appears
first because it is often the fastest and cheapest high-biomass screening
route, but the warning row should be read just as seriously: eukaryotic
processing, solubility, toxicity, secretion, membrane insertion, and activity
requirements can make a slower route more biologically correct.

The route candidates follow the same idea. They name plausible expression
routes and their missing inputs, while leaving the execution decision to the
reviewer.

## Step 4: Add A Product Sequence When One Exists

If the protein-coding product is already present in the GENtle project state,
pass its sequence id:

```bash
cargo run --quiet --bin gentle_cli -- --state path/to/project.gentle.json \
  planning protein-expression-handoff \
  --seq-id target_cds_or_protein_product \
  --objective '{"schema":"gentle.planning_objective.v1","biological_intent":"protein_expression_max_yield"}' \
  --format text
```

Expected difference:

- `product_definition.seq_id` records the supplied id
- `product_definition.sequence_present` tells whether the id is loaded
- if loaded, `sequence_context` records sequence name/type, nucleotide or
  protein length, feature count, GC percent/range, ambiguous bases, and
  relevant CDS/protein/tag annotations
- `cds_assessment` reports whether annotated CDS context or whole-sequence
  fallback context was used, plus nucleotide length, inferred protein length,
  start/stop sanity, internal stops, ambiguous codons, translation table, and
  warnings
- `tag_assessment` summarizes annotated affinity, solubility, epitope, or
  signal-tag context while keeping tag policy explicit
- if no usable CDS/protein context is found, the report asks for product
  boundaries instead of selecting an expression route
- if a CDS/protein context is inferable, the report shifts toward the
  expression-specific review questions
- `product_definition.readiness.status` determines the next actions: CDS
  candidates can point to GeneArt preflight and cloning consultation,
  protein-only targets point to reverse-translation or provider-target review,
  and unresolved sequences point back to CDS/ORF boundary inspection

This analysis is read-only. GENtle does not infer that the product is
expression-ready, codon-optimize it, mutate the sequence, create a construct,
choose a final host/vector/tag, or submit anything to a provider.

## Step 5: Inspect The GeneArt Preflight Scaffold

The handoff report includes a service scaffold pointing at the existing
GeneArt protein-expression example:

```bash
cargo run --quiet --bin gentle_cli -- services project-preflight \
  @docs/examples/external_services/geneart_protein_expression_request.json
```

Expected outcome:

- schema is `gentle.external_service_preflight.v1`
- provider is `geneart`
- service kind is `protein_expression`
- direct submission remains unavailable
- warnings keep the human-review boundary visible

After the product definition, outsourcing permission, and route constraints have
been reviewed, prepare the local quote packet:

```bash
cargo run --quiet --bin gentle_cli -- services project-quote \
  @docs/examples/external_services/geneart_protein_expression_request.json
```

This remains a local handoff artifact. It does not upload data, request pricing
from a provider API, or submit an order.

Use this only as a scaffold. Replace the synthetic tutorial payload with a
reviewed product definition before any real quote or provider handoff.

## Step 6: Ask For Cloning Strategy Only After Expression Constraints

Once the product metric, host, folding/PTM needs, toxicity policy, and endpoint
are known, ask the cloning consultation to rank routine families:

```bash
cargo run --quiet --bin gentle_cli -- planning consult cloning \
  --objective '{"schema":"gentle.planning_objective.v1","biological_intent":"protein_expression_max_yield"}' \
  --format json
```

That report answers a different question: which cloning routine and helper
vector path should be considered. It should not replace the protein-expression
handoff; it should follow it.

## ClawBio Use

For ClawBio or another agent, the safe instruction is:

- call `planning protein-expression-handoff`
- quote `biological_intent`, `product_definition`,
  `product_readiness`, `sequence_context`, `cds_assessment`,
  `tag_assessment`,
  `host_chassis_candidates`, `vector_route_candidates`, `missing_questions`,
  `service_handoff_candidates`, and `suggested_next_actions`
- do not invent a final expression host, promoter, vector, tag, or vendor
  submission from chat context alone

GUI note: the current slice is CLI/shared-shell first. A dedicated Synthetic
Biology or protein-expression handoff inspector for existing
`gentle.protein_expression_handoff.v1` reports is the natural next GUI slice,
not something this tutorial assumes exists today.

## Expected Outcome

At the end of this tutorial you should have:

- a `gentle.protein_expression_handoff.v1` report
- a clear distinction between product definition and expression route choice
- read-only sequence/readiness/CDS/tag context when a sequence id is supplied
- ranked but review-required chassis and route candidates
- explicit missing biological questions for either product boundaries or
  expression-route review
- a local GeneArt preflight scaffold and optional quote-packet handoff
- no created construct and no vendor-side action

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or
missing a useful figure, please open the matching tutorial issue template and
include the context copied from GENtle Help -> Tutorial -> Copy Feedback
Context.

- Tutorial title:
- Tutorial/chapter id:
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
