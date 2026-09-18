#!/usr/bin/env bash
set -euo pipefail

readonly GENTLE_BIN_DIR="/opt/gentle/bin"

is_apptainer_runtime() {
    [[ -n "${APPTAINER_NAME:-}" || -n "${APPTAINER_CONTAINER:-}" || -n "${SINGULARITY_NAME:-}" || -n "${SINGULARITY_CONTAINER:-}" ]]
}

show_help() {
    cat <<'EOF'
Usage:
  gentle-image cli [gentle_cli args...]
  gentle-image mcp [gentle_mcp args...]
  gentle-image examples-docs [gentle_examples_docs args...]
  gentle-image COMMAND [args...]

Modes:
  cli           Run the shared CLI adapter.
  mcp           Run the MCP server over stdio.
  examples-docs Run documentation/example helper commands.

The default is `cli --help`. This image does not include the GUI or embedded
JavaScript/Lua interfaces; use a native installation for those interfaces.

Under Apptainer/Singularity, `run IMAGE.sif SUBCOMMAND ...` treats unknown
subcommands as `gentle_cli SUBCOMMAND ...`. Docker callers should use the
explicit `cli` prefix, for example `docker run IMAGE cli capabilities`.
EOF
}

mode="${1:-cli}"
if [[ $# -gt 0 ]]; then
    shift
else
    set -- --help
fi

case "${mode}" in
    cli)
        exec "${GENTLE_BIN_DIR}/gentle_cli" "$@"
        ;;
    mcp)
        exec "${GENTLE_BIN_DIR}/gentle_mcp" "$@"
        ;;
    examples-docs)
        exec "${GENTLE_BIN_DIR}/gentle_examples_docs" "$@"
        ;;
    gui-web|gui|gui-x11|js|lua|gentle|gentle_js|gentle_lua)
        echo "This headless image does not include GUI, JavaScript or Lua interfaces. Use CLI/MCP here or a native installation." >&2
        exit 64
        ;;
    help|-h|--help)
        show_help
        ;;
    *)
        if is_apptainer_runtime && [[ -x "${GENTLE_BIN_DIR}/gentle_cli" ]]; then
            exec "${GENTLE_BIN_DIR}/gentle_cli" "${mode}" "$@"
        fi
        exec "${mode}" "$@"
        ;;
esac
