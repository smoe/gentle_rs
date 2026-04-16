#!/usr/bin/env bash
set -euo pipefail

STATE="${STATE:-/tmp/vkorc1_rs9923231_promoter_reporter.state.json}"

# Instance preflight
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  genomes status "Human GRCh38 Ensembl 116" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes

# Resolve rs9923231 and extract the annotated local context
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"FetchDbSnpRegion":{"rs_id":"rs9923231","genome_id":"Human GRCh38 Ensembl 116","flank_bp":3000,"output_id":"vkorc1_rs9923231_context","annotation_scope":"full","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}' \
  --confirm

# Derive promoter windows and summarize promoter context
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  variant annotate-promoters vkorc1_rs9923231_context \
  --gene-label VKORC1 \
  --upstream-bp 1000 \
  --downstream-bp 200

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  variant promoter-context vkorc1_rs9923231_context \
  --variant rs9923231 \
  --gene-label VKORC1 \
  --path docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/variant_promoter_context.json

# Ask GENtle for the recommended reporter fragment
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  variant reporter-fragments vkorc1_rs9923231_context \
  --variant rs9923231 \
  --gene-label VKORC1 \
  --retain-downstream-from-tss-bp 200 \
  --retain-upstream-beyond-variant-bp 500 \
  --path docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/promoter_reporter_candidates.json

# Materialize the current recommended baseline fragment
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"ExtractRegion":{"input":"vkorc1_rs9923231_context","from":2412,"to":3501,"output_id":"vkorc1_rs9923231_promoter_fragment"}}' \
  --confirm

# Create matched reference and alternate inserts
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  variant materialize-allele vkorc1_rs9923231_promoter_fragment \
  --variant rs9923231 \
  --allele reference \
  --output-id vkorc1_rs9923231_promoter_reference

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  variant materialize-allele vkorc1_rs9923231_promoter_fragment \
  --variant rs9923231 \
  --allele alternate \
  --output-id vkorc1_rs9923231_promoter_alternate

# Load the pinned local mammalian reporter backbone
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"LoadFile":{"path":"data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb","as_id":"gentle_mammalian_luciferase_backbone_v1"}}' \
  --confirm

# Preview the matched reporter pair
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"Ligation":{"inputs":["vkorc1_rs9923231_promoter_reference","gentle_mammalian_luciferase_backbone_v1"],"circularize_if_possible":false,"protocol":"Blunt","output_prefix":"vkorc1_rs9923231_reporter_reference_assembly","unique":false}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"Branch":{"input":"vkorc1_rs9923231_reporter_reference_assembly_1","output_id":"vkorc1_rs9923231_reporter_reference"}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"Ligation":{"inputs":["vkorc1_rs9923231_promoter_alternate","gentle_mammalian_luciferase_backbone_v1"],"circularize_if_possible":false,"protocol":"Blunt","output_prefix":"vkorc1_rs9923231_reporter_alternate_assembly","unique":false}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"Branch":{"input":"vkorc1_rs9923231_reporter_alternate_assembly_1","output_id":"vkorc1_rs9923231_reporter_alternate"}}' \
  --confirm

# Export the promoter-context figure through the same shared render op
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"SetLinearViewport":{"start_bp":2113,"span_bp":1388}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"RenderSequenceSvg":{"seq_id":"vkorc1_rs9923231_context","mode":"Linear","path":"docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg"}}' \
  --confirm

# Export the expected paired construct previews
cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"RenderSequenceSvg":{"seq_id":"vkorc1_rs9923231_reporter_reference","mode":"Circular","path":"docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_reference.svg"}}' \
  --confirm

cargo run --quiet --bin gentle_cli -- \
  --state "$STATE" \
  op '{"RenderSequenceSvg":{"seq_id":"vkorc1_rs9923231_reporter_alternate","mode":"Circular","path":"docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_alternate.svg"}}' \
  --confirm
