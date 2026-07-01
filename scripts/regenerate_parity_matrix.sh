#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
cargo run --quiet --bin gentle_examples_docs -- parity-matrix-generate
