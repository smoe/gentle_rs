# syntax=docker/dockerfile:1.7

ARG DEBIAN_SUITE=forky

# Build only headless GENtle binaries. Native GUI/JS/Lua distributions are separate.
FROM debian:${DEBIAN_SUITE}-slim AS build

ENV DEBIAN_FRONTEND=noninteractive \
    CARGO_HOME=/usr/local/cargo \
    RUSTUP_HOME=/usr/local/rustup \
    PATH=/usr/local/cargo/bin:/usr/local/rustup/bin:/usr/local/bin:/usr/bin:/bin

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    clang \
    cmake \
    git \
    libssl-dev \
    perl \
    pkg-config \
    rust-all \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/gentle

COPY Cargo.toml Cargo.lock build.rs ./
COPY vendor ./vendor
COPY packages ./packages
COPY crates ./crates
COPY src ./src
COPY assets ./assets
COPY data/resources/affymetrix/platform_registry.json ./data/resources/affymetrix/
COPY docs ./docs
COPY integrations/python ./integrations/python
COPY README.md CONTRIBUTING.md copyright ./

ARG GENTLE_CARGO_PROFILE=release-fast

# Guard the resolved Linux dependency graph, not just the selected binaries.
RUN cargo tree --locked --no-default-features --edges normal,build --prefix none > /tmp/gentle-dependencies.txt \
    && if grep -E '^(arboard|eframe|egui|egui_commonmark|egui_extras|gentle-gui|rfd|winit|deno_core|deno_error|v8|mlua|mlua-sys|lua-src|luajit-src) v' /tmp/gentle-dependencies.txt; then \
        echo "Desktop or embedded scripting dependency leaked into the headless build" >&2; exit 1; \
    fi
RUN cargo build --locked --profile "${GENTLE_CARGO_PROFILE}" --no-default-features \
    --bin gentle_cli --bin gentle_mcp --bin gentle_examples_docs -j1
RUN cargo install --locked --root /opt/rnapkin rnapkin

RUN mkdir -p /opt/gentle-dist-cli/bin /opt/gentle-dist-cli/integrations \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_cli" /opt/gentle-dist-cli/bin/gentle_cli \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_mcp" /opt/gentle-dist-cli/bin/gentle_mcp \
    && install -Dm755 "target/${GENTLE_CARGO_PROFILE}/gentle_examples_docs" /opt/gentle-dist-cli/bin/gentle_examples_docs \
    && cp -a assets /opt/gentle-dist-cli/assets \
    && cp -a docs /opt/gentle-dist-cli/docs \
    && cp -a integrations/python /opt/gentle-dist-cli/integrations/python \
    && cp README.md CONTRIBUTING.md copyright /opt/gentle-dist-cli/

# Runtime helpers and fonts support scientific computation and headless exports.
FROM debian:${DEBIAN_SUITE}-slim AS runtime-cli

ENV DEBIAN_FRONTEND=noninteractive \
    GENTLE_BIGWIG_TO_BEDGRAPH_BIN=/usr/local/bin/bigWigToBedGraph \
    GENTLE_BLASTN_BIN=/usr/bin/blastn \
    GENTLE_CONTAINER_FLAVOR=cli \
    GENTLE_MAKEBLASTDB_BIN=/usr/bin/makeblastdb \
    GENTLE_RNAFOLD_BIN=/usr/bin/RNAfold \
    GENTLE_RNAPKIN_BIN=/usr/local/bin/rnapkin \
    HOME=/home/gentle \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/gentle/bin:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/opt/gentle/integrations/python \
    PYTHONUNBUFFERED=1

RUN sed -i -E 's/^Components: main$/Components: main non-free/' /etc/apt/sources.list.d/debian.sources \
    && apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    fonts-dejavu-core \
    ncbi-blast+ \
    passwd \
    primer3 \
    python3 \
    python3-pybigwig \
    vienna-rna \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN groupadd --gid 1000 gentle \
    && useradd --uid 1000 --gid 1000 --create-home --shell /bin/bash gentle \
    && mkdir -p /opt/gentle /work \
    && chown -R gentle:gentle /home/gentle /opt/gentle /work

WORKDIR /opt/gentle

COPY --from=build /opt/gentle-dist-cli/ /opt/gentle/
COPY --from=build /opt/rnapkin/bin/rnapkin /usr/local/bin/rnapkin
COPY docker/bigWigToBedGraph /usr/local/bin/bigWigToBedGraph
COPY docker/entrypoint.sh /usr/local/bin/gentle-entrypoint

RUN chmod +x /usr/local/bin/bigWigToBedGraph /usr/local/bin/gentle-entrypoint \
    && chown -R gentle:gentle /opt/gentle

USER gentle
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/gentle-entrypoint"]
CMD ["cli", "--help"]
