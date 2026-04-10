# GENtle Container Guide

This document describes the Debian-first GENtle container images and how to use
them from Docker on macOS/Linux and from Apptainer on Linux.

The container goal is pragmatic:

- one image definition (`Dockerfile`) as the maintained source of truth
- two published runtime targets from that one definition:
  - a headless `cli` image for CLI/MCP/JS/Lua/Python wrapper use
  - a `gui` image layered on top for browser-served GUI use
- complete GENtle functionality across those images:
  - GUI
  - CLI
  - MCP server
  - embedded JavaScript shell
  - embedded Lua shell
  - Python wrapper
  - external helper-tool support for BLAST, Primer3, BigWig import, and RNA
    structure
- Debian packages whenever practical
- no separate Singularity/Apptainer definition unless Linux/HPC requirements
  later diverge enough to justify it

## Design Summary

- Base image: Debian `forky` (testing)
- Build toolchain: Debian `rust-all`
- Runtime targets:
  - `runtime-cli`: headless CLI/MCP/JS/Lua/Python wrapper image
  - `runtime-gui`: GUI image built on top of `runtime-cli`
- Runtime GUI delivery: `Xvfb` + `openbox` + `x11vnc` + `noVNC`
- Default published tags:
  - `ghcr.io/smoe/gentle_rs:cli`
  - `ghcr.io/smoe/gentle_rs:gui`
  - `ghcr.io/smoe/gentle_rs:latest` and the bare release tag remain GUI tags
- Linux/Apptainer strategy: consume the same OCI image family, rather than
  maintaining a second packaging language in parallel

This means:

- macOS and Linux both use the same Linux container image family
- macOS users normally run it through Docker Desktop, Colima, or OrbStack
- Linux users can choose Docker/Podman or Apptainer

Why the Dockerfile installs packages twice:

- the builder stage installs compilers and `-dev` headers
- the runtime stage installs only what the final image must execute

That split is deliberate. It keeps the final image from carrying the full Rust
and native build toolchain, even though both stages use Debian packages.

## Tool Coverage

The image is intended to include the helper tools that GENtle already knows how
to use:

- `blastn`
- `makeblastdb`
- `primer3`
- `bigWigToBedGraph`
- `rnapkin`

Container policy:

- prefer Debian packages when available
- keep non-Debian exceptions narrow and explicit

Current implementation detail:

- `bigWigToBedGraph` is provided by a compatibility wrapper script that uses
  Debian's `python3-pybigwig`
- `rnapkin` is installed via `cargo install` because it is not presently part
  of the Debian package set we rely on
- release/distribution follow-up: pin the `rnapkin` source/version explicitly
  once the preferred upstream release policy is fixed

Filesystem/layout note:

- `/opt/gentle`
  - GENtle payload: binaries, assets, docs, Python wrapper source
- `/usr/local/bin`
  - container-facing commands and compatibility shims:
    - `gentle-entrypoint`
    - `bigWigToBedGraph`
    - `rnapkin`

This follows normal container conventions: application payload under `/opt`,
generic executable entrypoints on `PATH`.

## Build the Images

Standard local GUI build:

```sh
docker build --target runtime-gui -t gentle:gui-local .
```

Standard local CLI build:

```sh
docker build --target runtime-cli -t gentle:cli-local .
```

Optional build arguments:

- `DEBIAN_SUITE`
  - defaults to `forky`
- `GENTLE_CARGO_PROFILE`
  - defaults to `release-fast`
  - set to `release` if you prefer fully optimized binaries over faster image
    builds

Example:

```sh
docker build \
  --build-arg DEBIAN_SUITE=forky \
  --build-arg GENTLE_CARGO_PROFILE=release \
  --target runtime-cli \
  -t gentle:cli-release .
```

Why `forky` by default:

- it keeps the image Debian-first while avoiding the current `sid` +
  `qemu-aarch64` failure mode in GitHub Actions multi-arch builds
- Debian testing (`forky`) currently ships a meaningfully newer `rust-all`
  than `trixie`, while still providing the runtime packages this image needs,
  including `novnc` and `python3-pybigwig`
- this keeps the project on Debian-packaged `rust-all` instead of switching
  the container over to `rustup`
- GENtle does not need audio output, so the image intentionally omits ALSA
  runtime packages instead of depending on the current `libasound2` virtual
  package split

## Run the GUI in a Browser

Default mode:

```sh
docker run --rm -it \
  -p 6080:6080 \
  -v "$(pwd)":/work \
  gentle:gui-local
```

Then open:

- <http://localhost:6080/vnc.html?autoconnect=1&resize=scale>

What happens internally:

- `Xvfb` creates a virtual X11 display
- `openbox` provides a lightweight desktop/session shell
- `x11vnc` exports that display locally inside the container
- `noVNC` publishes it to your browser
- GENtle starts inside that virtual desktop

This is the primary cross-platform GUI route because it behaves similarly on
macOS and Linux and avoids host-specific X11 setup on macOS.

## Run CLI / MCP / JS / Lua Modes

CLI:

```sh
docker run --rm -it \
  -v "$(pwd)":/work \
  gentle:cli-local capabilities
```

MCP:

```sh
docker run --rm -it \
  -v "$(pwd)":/work \
  gentle:cli-local mcp
```

JavaScript shell:

```sh
docker run --rm -it \
  -v "$(pwd)":/work \
  gentle:cli-local js
```

Lua shell:

```sh
docker run --rm -it \
  -v "$(pwd)":/work \
  gentle:cli-local lua
```

The GUI image entrypoint supports these modes:

- `gui-web`
- `gui`
- `cli`
- `mcp`
- `js`
- `lua`
- `examples-docs`

The CLI image ships the same headless modes:

- `cli`
- `mcp`
- `js`
- `lua`
- `examples-docs`

Run `--help` to see the dispatch summary:

```sh
docker run --rm gentle:gui-local --help
```

## Python Wrapper Inside the Image

The runtime image ships the thin Python wrapper from
`integrations/python/gentle_py` on `PYTHONPATH`, so you can use it directly
without an extra install step:

```sh
docker run --rm -it \
  -v "$(pwd)":/work \
  gentle:cli-local \
  python3 - <<'PY'
from gentle_py import GentleClient

client = GentleClient(state_path="/work/.gentle_state.json")
print(client.capabilities()["schema"])
PY
```

## Linux: Optional Native X11 Route

If you want GENtle to use your existing Linux X11 session instead of noVNC,
you can pass your display through and run `gui` mode:

```sh
docker run --rm -it \
  --net=host \
  -e DISPLAY \
  -v /tmp/.X11-unix:/tmp/.X11-unix \
  -v "$(pwd)":/work \
  gentle:gui-local gui
```

This route is intentionally secondary. The browser-served GUI is the default
because it is more consistent across macOS and Linux.

## Apptainer / Singularity on Linux

GENtle does not currently maintain a separate `.def` file. The intended Linux
path is to consume the same OCI image.

### From a published GHCR image

Once an OCI image is published, pull it directly. For headless Linux/HPC use,
prefer the CLI image:

```sh
apptainer pull gentle-cli.sif docker://ghcr.io/OWNER/REPO:cli
```

Then run CLI or shell commands as usual:

```sh
apptainer exec gentle-cli.sif gentle_cli capabilities
```

If you use:

```sh
apptainer run gentle-cli.sif capabilities
```

the CLI image entrypoint treats that as `gentle_cli capabilities`, which makes
`run` usable for quick headless smoke tests. For reproducible automation,
ClawBio/OpenClaw still prefers the explicit launcher shown below.

If you explicitly want the browser-served GUI image instead:

```sh
apptainer pull gentle-gui.sif docker://ghcr.io/OWNER/REPO:gui
apptainer run gentle-gui.sif gui-web
apptainer exec gentle-gui.sif gentle
```

The unsuffixed GUI tags (`:latest`, `:<tag>`) remain available for backward
compatibility, but `:gui` is the clearest pull target when you want the
browser-served GUI image explicitly.

Bare `apptainer run gentle-gui.sif` on the GUI image now prints a friendly
guidance message instead of immediately dropping into the Docker-oriented
`gui-web` default. That keeps headless Linux/HPC users from accidentally
starting the GUI path when they really meant `gentle_cli ...`.

### From a local Docker build

If you built the Docker image locally first:

```sh
docker save gentle:cli-local -o gentle-cli-local.tar
apptainer build gentle-cli-local.sif docker-archive://gentle-cli-local.tar
```

This keeps Docker as the only maintained image definition while still giving
Linux users a normal `.sif` artifact.

### ClawBio / OpenClaw launcher path

The copied ClawBio skill scaffold includes:

- `integrations/clawbio/skills/gentle-cloning/gentle_apptainer_cli.sh`

After pulling a `.sif`, set:

```sh
export GENTLE_CLI_CMD='skills/gentle-cloning/gentle_apptainer_cli.sh /absolute/path/to/gentle.sif'
```

That launcher resolves `apptainer` first and falls back to `singularity`,
binds the current working directory into `/work`, and executes `gentle_cli`
inside the image so the ClawBio wrapper can keep using the same
`GENTLE_CLI_CMD` contract it already uses for Docker.

Important platform note:

- Apptainer is a Linux-native runtime
- on macOS, use the Docker image directly unless you intentionally want an
  extra Linux VM layer just for Apptainer

## GitHub Publishing Setup

1. build the OCI image from this `Dockerfile`
2. publish both runtime targets to `ghcr.io`
3. tag GUI and CLI images distinctly for release tags
4. publish stable `linux/amd64` GUI and CLI images from GitHub Actions
5. let Linux/Apptainer users pull the same images through `docker://...`

Suggested GitHub Actions building blocks:

- `docker/login-action`
- `docker/metadata-action`
- `docker/build-push-action`

The repository now includes:

- `.github/workflows/container.yml`

which:

- builds both runtime targets (`runtime-cli` and `runtime-gui`) as a check on
  pull requests and `main`
- publishes `linux/amd64` GHCR images for release tags matching `v*`
- keeps unsuffixed GUI tags (`latest`, `v...`) for backward compatibility
- also publishes explicit GUI tags (`gui`, `v...-gui`) and headless tags
  (`cli`, `v...-cli`) so human GUI use and headless automation can target the
  right image directly
- updates the `latest` image tag only from those release-tag publishes

Current arm64 note:

- the release workflow intentionally does **not** publish `linux/arm64` from
  GitHub Actions right now
- the failure mode is an emulated `qemu-aarch64` crash during Debian package
  configuration in the container build, not a confirmed GENtle application
  failure
- practical next options are:
  - a native arm64 builder/runner
  - a future headless-only/container-specific build split that drops the GUI
    toolchain from the published image

This keeps maintenance low while still covering:

- Docker Desktop / Colima / OrbStack on macOS
- Docker/Podman on Linux
- Apptainer on Linux

## Storage / Mounting Guidance

Inside the container, the working directory is:

- `/work`

Recommended pattern:

- mount your project directory or scratch directory at `/work`
- keep GENtle state/project files there
- keep large genome/helper caches in a mounted subdirectory if you want them
  to persist across container runs

Example:

```sh
docker run --rm -it \
  -p 6080:6080 \
  -v "$(pwd)":/work \
  gentle:gui-local
```

## Current Limits / Follow-up

- the in-tree implementation now covers:
  - the shared builder
  - the `runtime-cli` image
  - the `runtime-gui` image
  - runtime dispatch for Docker and Apptainer/Singularity
- GitHub image publishing is implemented in:
  - `.github/workflows/container.yml`
- Apptainer currently depends on consuming the OCI image; there is no
  dedicated native definition file yet
