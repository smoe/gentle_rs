#!/usr/bin/env bash
set -euo pipefail

show_help() {
    cat <<'EOF'
Usage:
  gentle_local_checkout_cli.sh [gentle_cli args...]

Environment overrides:
  GENTLE_REPO_ROOT             Absolute path to a local GENtle checkout.
                               Required when this script has been copied into a
                               separate ClawBio checkout.
  CARGO_TARGET_DIR             Cargo target dir. Defaults to
                               $GENTLE_REPO_ROOT/target.
  GENTLE_REFERENCE_CACHE_DIR   Prepared reference cache root. Defaults to
                               $GENTLE_REPO_ROOT/data/genomes.
  GENTLE_HELPER_CACHE_DIR      Prepared helper cache root. Defaults to
                               $GENTLE_REPO_ROOT/data/helper_genomes.

This launcher resolves a local GENtle checkout and then executes:

  cargo run --quiet --manifest-path REPO/Cargo.toml --bin gentle_cli -- ...

Typical ClawBio/OpenClaw usage:

  export GENTLE_REPO_ROOT=/home/clawbio/GENtle
  export GENTLE_CLI_CMD=/home/clawbio/ClawBio/skills/gentle-cloning/gentle_local_checkout_cli.sh
  python clawbio.py run gentle-cloning --demo
EOF
}

resolve_repo_root() {
    if [[ -n "${GENTLE_REPO_ROOT:-}" ]]; then
        printf '%s\n' "${GENTLE_REPO_ROOT}"
        return 0
    fi

    local dir
    dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    while [[ "${dir}" != "/" ]]; do
        if [[ -f "${dir}/Cargo.toml" && -f "${dir}/src/bin/gentle_cli.rs" ]]; then
            printf '%s\n' "${dir}"
            return 0
        fi
        dir="$(dirname "${dir}")"
    done

    echo "Could not resolve a local GENtle checkout. Set GENTLE_REPO_ROOT=/absolute/path/to/GENtle." >&2
    return 1
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    show_help
    exit 0
fi

if ! command -v cargo >/dev/null 2>&1; then
    echo "Could not find 'cargo' on PATH" >&2
    exit 1
fi

repo_root="$(resolve_repo_root)"
manifest_path="${repo_root}/Cargo.toml"

export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"
export GENTLE_REFERENCE_CACHE_DIR="${GENTLE_REFERENCE_CACHE_DIR:-${repo_root}/data/genomes}"
export GENTLE_HELPER_CACHE_DIR="${GENTLE_HELPER_CACHE_DIR:-${repo_root}/data/helper_genomes}"

exec cargo run --quiet --manifest-path "${manifest_path}" --bin gentle_cli -- "$@"
