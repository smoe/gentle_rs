#!/usr/bin/env bash
set -euo pipefail

# Instance preflight
cargo run --quiet --bin gentle_cli -- \
  genomes status "Human GRCh38 Ensembl 116" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes

# Resolve rs9923231 and extract the annotated local context
cargo run --quiet --bin gentle_cli -- \
  op '{"FetchDbSnpRegion":{"rs_id":"rs9923231","genome_id":"Human GRCh38 Ensembl 116","flank_bp":3000,"output_id":"vkorc1_rs9923231_context","annotation_scope":"full","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}' \
  --confirm

# Extract the default promoter fragment
cargo run --quiet --bin gentle_cli -- \
  op '{"ExtractRegion":{"input":"vkorc1_rs9923231_context","from":2412,"to":3501,"output_id":"vkorc1_rs9923231_promoter_ref"}}' \
  --confirm

# Import the mammalian luciferase reporter backbone
cargo run --quiet --bin gentle_cli -- \
  shell 'genbank fetch AY738222 --as-id promega_luciferase_ay738222'

# Replay the current GUI-parity workflow skeleton
cargo run --quiet --bin gentle_cli -- \
  workflow @docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json

