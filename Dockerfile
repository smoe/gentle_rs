# syntax=docker/dockerfile:1.7

ARG DEBIAN_SUITE=forky

# Build only headless GENtle binaries. Native GUI/JS/Lua distributions are separate.
FROM debian:${DEBIAN_SUITE}-slim AS build

ENV DEBIAN_FRONTEND=noninteractive \
    CARGO_HOME=/usr/local/cargo \
    CARGO_INCREMENTAL=0 \
    CARGO_PROFILE_DEV_DEBUG=0 \
    RUSTUP_HOME=/usr/local/rustup \
    PATH=/usr/local/cargo/bin:/usr/local/rustup/bin:/usr/local/bin:/usr/bin:/bin

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    clang \
    cmake \
    curl \
    fonts-dejavu-core \
    git \
    libfontconfig1-dev \
    libfreetype6-dev \
    libssl-dev \
    perl \
    pkg-config \
    rust-all \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# RNAPKIN's published lockfile selects a bitmap backend with unaligned pointer
# dereferences. Use the scoped corrected lock, not optimization to hide the bug.
COPY docker/rnapkin/Cargo.lock /tmp/rnapkin.Cargo.lock
RUN curl --fail --show-error --silent --location \
        https://static.crates.io/crates/rnapkin/rnapkin-0.3.9.crate -o /tmp/rnapkin.crate \
    && echo "4495690197e1cced9b16234d6a66b40ebf190e9613b9c1c6aea837adeed00f17  /tmp/rnapkin.crate" | sha256sum -c - \
    && mkdir -p /opt/rnapkin-src \
    && tar -xzf /tmp/rnapkin.crate -C /opt/rnapkin-src --strip-components=1 \
    && cp /tmp/rnapkin.Cargo.lock /opt/rnapkin-src/Cargo.lock \
    && cargo install --locked --debug --path /opt/rnapkin-src --root /opt/rnapkin -j1 \
    && rm -rf /opt/rnapkin-src /tmp/rnapkin.crate /tmp/rnapkin.Cargo.lock

# Exercise both renderers before GENtle's long build. The final-image smoke
# still repeats this as the unprivileged runtime user with networking disabled.
RUN smoke_dir="$(mktemp -d)" \
    && printf "%s\n" "GGGAAACCC" "(((...)))" > "$smoke_dir/hairpin.dbn" \
    && timeout 30 /opt/rnapkin/bin/rnapkin --height 128 -o "$smoke_dir/hairpin.svg" "$smoke_dir/hairpin.dbn" \
    && timeout 30 /opt/rnapkin/bin/rnapkin --height 128 -o "$smoke_dir/hairpin.png" "$smoke_dir/hairpin.dbn" \
    && test -s "$smoke_dir/hairpin.svg" \
    && grep -q "<svg" "$smoke_dir/hairpin.svg" \
    && test -s "$smoke_dir/hairpin.png" \
    && rm -rf "$smoke_dir"

WORKDIR /opt/gentle

COPY Cargo.toml Cargo.lock build.rs ./
COPY vendor ./vendor
COPY packages ./packages
COPY crates ./crates
COPY src ./src
COPY assets ./assets
# build.rs validates these resources even when GUI features are disabled.
COPY icons ./icons
COPY data/resources/affymetrix/platform_registry.json ./data/resources/affymetrix/
COPY docs ./docs
COPY integrations/python ./integrations/python
COPY README.md CONTRIBUTING.md copyright ./

# Guard the resolved Linux dependency graph, not just the selected binaries.
RUN cargo tree --locked --no-default-features --edges normal,build --prefix none > /tmp/gentle-dependencies.txt \
    && if grep -E '^(arboard|eframe|egui|egui_commonmark|egui_extras|gentle-gui|rfd|winit|deno_core|deno_error|v8|mlua|mlua-sys|lua-src|luajit-src) v' /tmp/gentle-dependencies.txt; then \
        echo "Desktop or embedded scripting dependency leaked into the headless build" >&2; exit 1; \
    fi
# Same bounded opt-level=1 recipe as native installers; helpers stay dev.
ARG GENTLE_GIT_COMMIT=""
ENV GENTLE_GIT_COMMIT=${GENTLE_GIT_COMMIT}
RUN cargo build --locked --profile package-opt1 --no-default-features \
    --bin gentle_cli --bin gentle_mcp --bin gentle_examples_docs -j1

RUN mkdir -p /opt/gentle-dist-cli/bin /opt/gentle-dist-cli/integrations \
    && install -Dm755 "target/package-opt1/gentle_cli" /opt/gentle-dist-cli/bin/gentle_cli \
    && install -Dm755 "target/package-opt1/gentle_mcp" /opt/gentle-dist-cli/bin/gentle_mcp \
    && install -Dm755 "target/package-opt1/gentle_examples_docs" /opt/gentle-dist-cli/bin/gentle_examples_docs \
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
    libfontconfig1 \
    libfreetype6 \
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
