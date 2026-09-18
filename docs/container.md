# GENtle Container Guide

The maintained Debian-first OCI image is **headless**: CLI, MCP, the examples
helper and the Python CLI wrapper. It does not compile or redistribute the
GENtle GUI or embedded JavaScript/Lua interfaces. Native installations retain
those interfaces.

## Build And Distribution

One Dockerfile supplies one runtime target, `runtime-cli`:

```sh
docker build --target runtime-cli -t gentle:cli-local .
```

The builder uses Debian `forky` and `rust-all`. It runs:

```sh
cargo build --locked --profile release-fast --no-default-features \
  --bin gentle_cli --bin gentle_mcp --bin gentle_examples_docs -j1
```

An earlier image split removed desktop packages only from the final CLI image;
its shared builder still compiled GUI and scripting code. This build disables
those features before compilation and rejects desktop/V8/Lua dependencies in
the resolved normal/build dependency graph.

Optional build arguments remain `DEBIAN_SUITE` (default `forky`) and
`GENTLE_CARGO_PROFILE` (default `release-fast`). Selecting `release` increases
optimization and can require more build memory. Serial Cargo compilation is
not a guarantee that an individual compiler process fits the runner.

After an owner-approved publication, these tags name the same headless image:

- `ghcr.io/smoe/gentle_rs:cli`
- `ghcr.io/smoe/gentle_rs:<tag>-cli`
- `ghcr.io/smoe/gentle_rs:<tag>`
- `ghcr.io/smoe/gentle_rs:latest`

**Migration:** new bare release tags and `latest` no longer launch a browser
GUI. Historical `gui` / `<tag>-gui` images are not refreshed or deleted by this
workflow. Use native GENtle for interactive windows. Prefer a versioned CLI
tag or digest for reproducible automation rather than assuming an older
floating tag already contains these changes.

Published images target `linux/amd64`. Native arm64 container acceptance
remains pending; the earlier emulated arm64 packaging failure is not resolved
by this change.

## Run CLI And MCP

Mount project files and persistent caches beneath `/work`, the runtime working
directory. The container runs as unprivileged UID/GID 1000; grant that user
appropriate access to the mounted files.

```sh
docker run --rm -v "$(pwd)":/work gentle:cli-local cli capabilities
docker run --rm -i -v "$(pwd)":/work gentle:cli-local \
  mcp --state /work/project.gentle.json
docker run --rm gentle:cli-local examples-docs --help
docker run --rm gentle:cli-local --help
```

MCP uses stdio: keep stdin open with `-i`, but do not allocate a TTY.
With no arguments the image runs `cli --help`. Explicit `gui`, `gui-web`,
`js` and `lua` requests fail with an explanatory message, rather than
attempting to launch an omitted binary.

Docker callers should use `cli capabilities`, not bare `capabilities`.
Other explicit executables on PATH, such as `python3` or `blastn`, can be
invoked directly.

## Python Wrapper

The thin wrapper is available on `PYTHONPATH` without an extra installation:

```sh
docker run --rm -i -v "$(pwd)":/work gentle:cli-local python3 - <<'PY'
from gentle_py import GentleClient

client = GentleClient(state_path="/work/.gentle_state.json")
print(client.capabilities()["schema"])
PY
```

The wrapper delegates to `gentle_cli`; it does not require embedded JS or Lua.

## Scientific Helpers And Resources

The image retains BLAST (`blastn` / `makeblastdb`), Primer3, ViennaRNA
(`RNAfold`), BigWig conversion and `rnapkin`. Assets and documentation live
under `/opt/gentle`, alongside the binaries and Python wrapper. DejaVu fonts
support headless figure exports without installing a desktop stack.

- Debian packages supply the tools where available; ViennaRNA requires the
  explicitly enabled Debian `non-free` component.
- `bigWigToBedGraph` is a compatibility wrapper over `python3-pybigwig`.
- `rnapkin` remains the existing `cargo install --locked` exception. Its
  upstream version is not pinned; pinning it remains a reproducibility follow-up.
- Compilers and development headers remain in the builder, not the runtime.
- Large reference data and repository test fixtures are not bundled. Mount
  inputs needed by particular workflows; shipping the examples helper is not
  a claim that every tutorial can replay without additional files.

## Apptainer / Singularity

Linux/HPC users consume the same OCI image, without a separate definition file:

```sh
apptainer pull gentle-cli.sif docker://ghcr.io/smoe/gentle_rs:cli
apptainer exec gentle-cli.sif gentle_cli capabilities
apptainer run gentle-cli.sif capabilities
```

Under Apptainer/Singularity, unknown entrypoint subcommands retain the existing
`gentle_cli SUBCOMMAND ...` dispatch. With no arguments, the default is now
CLI help, never a GUI.

For a locally built image:

```sh
docker save gentle:cli-local -o gentle-cli-local.tar
apptainer build gentle-cli-local.sif docker-archive://gentle-cli-local.tar
```

The copied ClawBio/OpenClaw skill scaffold provides
`integrations/clawbio/skills/gentle-cloning/gentle_apptainer_cli.sh`. In the
ClawBio checkout, configure it with:

```sh
export GENTLE_CLI_CMD='skills/gentle-cloning/gentle_apptainer_cli.sh /absolute/path/to/gentle-cli.sif'
```

It selects Apptainer or Singularity, binds the working directory at `/work`,
and executes `gentle_cli`. On macOS, use Docker Desktop, Colima or OrbStack
rather than assuming a native Apptainer runtime.

## CI And Publication

`.github/workflows/container.yml` builds `runtime-cli` on `v*` tag pushes,
published-release events and explicit manual dispatch. Ordinary pushes and
pull requests run the lightweight container-contract tests, not Docker builds.

The container check loads the image, runs CLI/MCP/docs entrypoints with
networking disabled, checks missing shared libraries and bundled assets,
imports the Python wrapper, and verifies GUI/JS/Lua binaries are absent.
Its receipt records the candidate SHA, lockfile hash, image identity,
`default_features: false`, empty optional features and explicit binaries.

Only an owner-approved publication event publishes GHCR tags; a tag push or
build-only dispatch does not. Native installer jobs and their feature sets
are independent and unchanged. See [release process](release.md) for immutable
candidate checks and manual validation commands.

The actual reduced-image build, helper/runtime checks and peak-memory behavior
must pass on the CI runner. Offline launcher tests or a local `cargo check`
do not establish Docker release acceptance.
