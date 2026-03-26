#!/usr/bin/env bash
set -euo pipefail

readonly GENTLE_BIN_DIR="/opt/gentle/bin"

export DISPLAY="${DISPLAY:-:99}"
export GENTLE_NOVNC_PORT="${GENTLE_NOVNC_PORT:-6080}"
export GENTLE_VNC_PORT="${GENTLE_VNC_PORT:-5900}"
export GENTLE_XVFB_WHD="${GENTLE_XVFB_WHD:-1920x1080x24}"
export XDG_RUNTIME_DIR="${XDG_RUNTIME_DIR:-/tmp/gentle-runtime}"
export LIBGL_ALWAYS_SOFTWARE="${LIBGL_ALWAYS_SOFTWARE:-1}"
export GALLIUM_DRIVER="${GALLIUM_DRIVER:-llvmpipe}"
export MESA_LOADER_DRIVER_OVERRIDE="${MESA_LOADER_DRIVER_OVERRIDE:-llvmpipe}"

_cleanup_pids=()

cleanup() {
    local pid
    for pid in "${_cleanup_pids[@]:-}"; do
        kill "${pid}" >/dev/null 2>&1 || true
    done
}

trap cleanup EXIT INT TERM

find_novnc_proxy() {
    local candidate
    for candidate in \
        /usr/share/novnc/utils/novnc_proxy \
        /usr/share/novnc/utils/launch.sh \
        /opt/novnc/utils/novnc_proxy
    do
        if [[ -x "${candidate}" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
    done
    echo "Could not find noVNC launcher script" >&2
    return 1
}

wait_for_x_socket() {
    local display_num="${DISPLAY#:}"
    local x_socket="/tmp/.X11-unix/X${display_num}"
    local tries=0
    while [[ ! -S "${x_socket}" ]]; do
        tries=$((tries + 1))
        if (( tries > 100 )); then
            echo "Timed out waiting for X socket ${x_socket}" >&2
            return 1
        fi
        sleep 0.1
    done
}

start_gui_web() {
    local novnc_proxy
    mkdir -p "${XDG_RUNTIME_DIR}" /work
    chmod 700 "${XDG_RUNTIME_DIR}"

    Xvfb "${DISPLAY}" -screen 0 "${GENTLE_XVFB_WHD}" -ac +extension GLX +render -noreset \
        >/tmp/gentle-xvfb.log 2>&1 &
    _cleanup_pids+=("$!")
    wait_for_x_socket

    openbox >/tmp/gentle-openbox.log 2>&1 &
    _cleanup_pids+=("$!")

    x11vnc \
        -display "${DISPLAY}" \
        -rfbport "${GENTLE_VNC_PORT}" \
        -localhost \
        -forever \
        -shared \
        -nopw \
        -o /tmp/gentle-x11vnc.log \
        >/dev/null 2>&1 &
    _cleanup_pids+=("$!")

    novnc_proxy="$(find_novnc_proxy)"
    "${novnc_proxy}" \
        --listen "${GENTLE_NOVNC_PORT}" \
        --vnc "localhost:${GENTLE_VNC_PORT}" \
        >/tmp/gentle-novnc.log 2>&1 &
    _cleanup_pids+=("$!")

    echo "GENtle GUI available at http://localhost:${GENTLE_NOVNC_PORT}/vnc.html?autoconnect=1&resize=scale"
    "${GENTLE_BIN_DIR}/gentle" "$@" &
    local gui_pid=$!
    wait "${gui_pid}"
}

show_help() {
    cat <<'EOF'
Usage:
  gentle-image gui-web [gentle args...]
  gentle-image gui [gentle args...]
  gentle-image cli [gentle_cli args...]
  gentle-image mcp [gentle_mcp args...]
  gentle-image js [gentle_js args...]
  gentle-image lua [gentle_lua args...]
  gentle-image examples-docs [gentle_examples_docs args...]
  gentle-image COMMAND [args...]

Modes:
  gui-web       Launch GENtle under Xvfb/openbox and expose it through noVNC.
  gui           Launch GENtle against the current DISPLAY.
  cli           Run the shared CLI adapter.
  mcp           Run the MCP server binary.
  js            Run the embedded JavaScript shell adapter.
  lua           Run the embedded Lua shell adapter.
  examples-docs Run documentation/example helper commands.

The default mode is `gui-web`.
EOF
}

mode="${1:-gui-web}"
if [[ $# -gt 0 ]]; then
    shift
fi

case "${mode}" in
    gui-web)
        start_gui_web "$@"
        ;;
    gui|gui-x11)
        exec "${GENTLE_BIN_DIR}/gentle" "$@"
        ;;
    cli)
        exec "${GENTLE_BIN_DIR}/gentle_cli" "$@"
        ;;
    mcp)
        exec "${GENTLE_BIN_DIR}/gentle_mcp" "$@"
        ;;
    js)
        exec "${GENTLE_BIN_DIR}/gentle_js" "$@"
        ;;
    lua)
        exec "${GENTLE_BIN_DIR}/gentle_lua" "$@"
        ;;
    examples-docs)
        exec "${GENTLE_BIN_DIR}/gentle_examples_docs" "$@"
        ;;
    help|-h|--help)
        show_help
        ;;
    *)
        exec "${mode}" "$@"
        ;;
esac
