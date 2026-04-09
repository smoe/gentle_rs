# syntax=docker/dockerfile:1.7

ARG DEBIAN_SUITE=forky

# Builder stage:
# - uses Debian rust-all as requested
# - keeps compiler/dev headers out of the final runtime images
# - builds every GENtle binary plus the current rnapkin exception
FROM debian:${DEBIAN_SUITE}-slim AS build

ENV DEBIAN_FRONTEND=noninteractive \
    CARGO_HOME=/usr/local/cargo \
    RUSTUP_HOME=/usr/local/rustup \
    PATH=/usr/local/cargo/bin:/usr/local/rustup/bin:/usr/local/bin:/usr/bin:/bin

# Builder-only packages: compilers, headers, and cargo/git support.
# These are intentionally separate from the runtime packages below so the final
# image does not need to ship the full Rust/C/C++ toolchain.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    clang \
    cmake \
    git \
    libegl1-mesa-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libgbm-dev \
    libglib2.0-dev \
    libgl1-mesa-dev \
    libgtk-3-dev \
    libnss3-dev \
    libssl-dev \
    libwayland-dev \
    libx11-dev \
    libx11-xcb-dev \
    libxcb1-dev \
    libxcursor-dev \
    libxinerama-dev \
    libxkbcommon-dev \
    libxi-dev \
    libxrandr-dev \
    perl \
    pkg-config \
    rust-all \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/gentle

COPY Cargo.toml Cargo.lock build.rs ./
COPY vendor ./vendor
COPY crates ./crates
COPY src ./src
COPY assets ./assets
COPY docs ./docs
COPY icons ./icons
COPY integrations/python ./integrations/python
COPY README.md CONTRIBUTING.md copyright ./

ARG GENTLE_CARGO_PROFILE=release-fast

RUN cargo build --locked --profile "${GENTLE_CARGO_PROFILE}" --features script-interfaces --bins
RUN cargo install --locked --root /opt/rnapkin rnapkin

RUN mkdir -p \
        /opt/gentle-dist-cli/bin \
        /opt/gentle-dist-cli/integrations \
        /opt/gentle-dist-gui/bin \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_cli" /opt/gentle-dist-cli/bin/gentle_cli \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_mcp" /opt/gentle-dist-cli/bin/gentle_mcp \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_js" /opt/gentle-dist-cli/bin/gentle_js \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_lua" /opt/gentle-dist-cli/bin/gentle_lua \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_examples_docs" /opt/gentle-dist-cli/bin/gentle_examples_docs \
    && cp -a docs /opt/gentle-dist-cli/docs \
    && cp -a integrations/python /opt/gentle-dist-cli/integrations/python \
    && cp README.md /opt/gentle-dist-cli/README.md \
    && cp CONTRIBUTING.md /opt/gentle-dist-cli/CONTRIBUTING.md \
    && cp copyright /opt/gentle-dist-cli/copyright \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle" /opt/gentle-dist-gui/bin/gentle \
    && cp -a assets /opt/gentle-dist-gui/assets \
    && cp -a icons /opt/gentle-dist-gui/icons

# Headless runtime stage:
# - carries the CLI/MCP/JS/Lua/Python wrapper contract
# - is the preferred runtime for Apptainer/Singularity, ClawBio, and MCP/agent use
FROM debian:${DEBIAN_SUITE}-slim AS runtime-cli

ENV DEBIAN_FRONTEND=noninteractive \
    GENTLE_BIGWIG_TO_BEDGRAPH_BIN=/usr/local/bin/bigWigToBedGraph \
    GENTLE_BLASTN_BIN=/usr/bin/blastn \
    GENTLE_CONTAINER_FLAVOR=cli \
    GENTLE_MAKEBLASTDB_BIN=/usr/bin/makeblastdb \
    GENTLE_RNAPKIN_BIN=/usr/local/bin/rnapkin \
    HOME=/home/gentle \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/gentle/bin:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/opt/gentle/integrations/python \
    PYTHONUNBUFFERED=1 \
    XDG_RUNTIME_DIR=/tmp/gentle-runtime

# Runtime packages only: helper tools and headless dependencies. Keep the CLI
# image narrow so Apptainer/Singularity and MCP use do not carry the GUI stack.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    ncbi-blast+ \
    passwd \
    primer3 \
    python3 \
    python3-pybigwig \
    && rm -rf /var/lib/apt/lists/*

RUN groupadd --gid 1000 gentle \
    && useradd --uid 1000 --gid 1000 --create-home --shell /bin/bash gentle \
    && mkdir -p /opt/gentle /work /tmp/gentle-runtime \
    && chown -R gentle:gentle /home/gentle /opt/gentle /work /tmp/gentle-runtime

WORKDIR /opt/gentle

COPY --from=build /opt/gentle-dist-cli/ /opt/gentle/

# Container-facing launchers and compatibility shims live in /usr/local/bin so
# they behave like normal commands on PATH. /opt/gentle is reserved for the
# application payload itself (binaries, docs, Python wrapper source).
COPY --from=build /opt/rnapkin/bin/rnapkin /usr/local/bin/rnapkin
COPY docker/bigWigToBedGraph /usr/local/bin/bigWigToBedGraph
COPY docker/entrypoint.sh /usr/local/bin/gentle-entrypoint

RUN chmod +x /usr/local/bin/bigWigToBedGraph /usr/local/bin/gentle-entrypoint \
    && chown -R gentle:gentle /opt/gentle

USER gentle
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/gentle-entrypoint"]
CMD ["cli", "--help"]

# GUI runtime stage:
# - layers the browser-served GUI stack on top of the headless CLI image
# - keeps the CLI image smaller while preserving the existing Docker GUI route
FROM runtime-cli AS runtime-gui

USER root

ENV DISPLAY=:99 \
    GALLIUM_DRIVER=llvmpipe \
    GENTLE_CONTAINER_FLAVOR=gui \
    GENTLE_NOVNC_PORT=6080 \
    GENTLE_VNC_PORT=5900 \
    GENTLE_XVFB_WHD=1920x1080x24 \
    LIBGL_ALWAYS_SOFTWARE=1 \
    MESA_LOADER_DRIVER_OVERRIDE=llvmpipe

# Runtime GUI packages only: browser/VNC desktop plumbing and dynamic graphics
# libraries needed by the native GUI binary.
RUN apt-get update && apt-get install -y --no-install-recommends \
    dbus-x11 \
    fonts-dejavu-core \
    fonts-noto-core \
    libegl1 \
    libfontconfig1 \
    libfreetype6 \
    libgbm1 \
    libglib2.0-0 \
    libgl1 \
    libgtk-3-0 \
    libnss3 \
    libwayland-client0 \
    libwayland-egl1 \
    libx11-6 \
    libx11-xcb1 \
    libxcb1 \
    libxcursor1 \
    libxext6 \
    libxinerama1 \
    libxkbcommon0 \
    libxkbcommon-x11-0 \
    libxi6 \
    libxrandr2 \
    novnc \
    openbox \
    websockify \
    x11vnc \
    xauth \
    xdg-utils \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/gentle-dist-gui/ /opt/gentle/

RUN chown -R gentle:gentle /opt/gentle

USER gentle
WORKDIR /work

EXPOSE 6080

ENTRYPOINT ["/usr/local/bin/gentle-entrypoint"]
CMD ["gui-web"]
