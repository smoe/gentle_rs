#!/usr/bin/env bash
set -euo pipefail

show_help() {
    cat <<'EOF'
Usage:
  gentle_apptainer_cli.sh IMAGE.sif [gentle_cli args...]

Environment overrides:
  GENTLE_APPTAINER_IMAGE       Path to the Apptainer/Singularity image.
  GENTLE_APPTAINER_HOST_WORKDIR
                               Host directory to bind into the container.
                               Defaults to the current working directory.
  GENTLE_APPTAINER_WORKDIR     Container workdir/bind target. Defaults to /work.

This launcher resolves `apptainer` first and falls back to `singularity`.
It then executes:

  <runtime> exec --bind HOSTDIR:WORKDIR --pwd WORKDIR IMAGE.sif gentle_cli ...

Typical ClawBio/OpenClaw usage:

  apptainer pull gentle.sif docker://ghcr.io/smoe/gentle_rs:cli
  export GENTLE_CLI_CMD='./gentle_apptainer_cli.sh /absolute/path/to/gentle.sif'
  python clawbio.py run gentle-cloning --demo
EOF
}

find_runtime() {
    if command -v apptainer >/dev/null 2>&1; then
        printf '%s\n' apptainer
        return 0
    fi
    if command -v singularity >/dev/null 2>&1; then
        printf '%s\n' singularity
        return 0
    fi
    echo "Could not find 'apptainer' or 'singularity' on PATH" >&2
    return 1
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    show_help
    exit 0
fi

image_path="${GENTLE_APPTAINER_IMAGE:-}"
if [[ -z "${image_path}" ]]; then
    if [[ $# -lt 1 ]]; then
        show_help >&2
        exit 2
    fi
    image_path="$1"
    shift
fi

host_workdir="${GENTLE_APPTAINER_HOST_WORKDIR:-$PWD}"
container_workdir="${GENTLE_APPTAINER_WORKDIR:-/work}"
runtime="$(find_runtime)"

exec "${runtime}" exec \
    --bind "${host_workdir}:${container_workdir}" \
    --pwd "${container_workdir}" \
    "${image_path}" \
    gentle_cli \
    "$@"
