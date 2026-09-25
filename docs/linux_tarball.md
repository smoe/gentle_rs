# Linux Tarball Quick Start

The current release workflow builds a downloadable
`gentle-<tag>-linux-x64-dev.tar.gz` for internal versions, or
`gentle-<tag>-linux-x64.tar.gz` for other versions, on Ubuntu 24.04 x86-64.
Internal packages use unoptimized Cargo `dev` binaries; their performance is
not equivalent to the optimized final-release profile. This is a relocatable
directory, not a system installer, Debian package, AppImage or static binary.
Availability is confirmed only when the tag's release workflow passes and
publishes the archive.

## Runtime Requirements

Use Ubuntu 24.04 x86-64 or a system with compatible glibc (2.39 or newer) and
runtime libraries. Other distributions are not certified merely because the
archive extracts. The package's `runtime-libraries.txt` records the libraries
resolved on its build host, rather than bundling that host's shared libraries.

On Ubuntu 24.04, the desktop runtime libraries can be installed with:

```bash
sudo apt-get install ca-certificates libegl1 libgl1 libvulkan1 \
  libfontconfig1 libfreetype6 libgbm1 libgtk-3-0t64 libnss3 libssl3t64 \
  libx11-6 libx11-xcb1 libxcb1 libxcb-render0 libxcb-shape0 libxcb-xfixes0 \
  libxcursor1 libxinerama1 libxkbcommon0 libxi6 libxrandr2 \
  libwayland-client0 libwayland-cursor0 libwayland-egl1
```

An X11 or Wayland desktop session is needed for the GUI. BLAST+, Primer3 and
other external scientific tools are not bundled; install them separately when
using operations that need them. Reference genomes and credentials are not
included. Use the GHCR CLI image for the curated container environment instead
of treating this GUI-enabled tarball as a minimal headless distribution.

## Extract And Run

For a `.11` archive produced by the current workflow, verify its receipt
against the exact candidate SHA. Publication is not installed-package or
scientific acceptance. Use the filename from that run; historical archives
keep their original names and recipes:

```bash
tar -xzf gentle-v0.1.0-internal.11-linux-x64-dev.tar.gz
cd gentle-v0.1.0-internal.11-linux-x64
sha256sum -c SHA256SUMS
./bin/gentle --version
./bin/gentle
```

Start GENtle from the extracted directory so relative resource and tutorial
paths resolve. Save personal projects outside the release directory when you
want to replace or remove that directory for an upgrade.

The current packaging recipe contains the GUI, CLI, MCP server,
example-document generator and publication-report entrypoint in `bin/`, plus
tracked resources, tutorial docs/fixtures and the Python adapter source.
Embedded JavaScript/Lua entrypoints are omitted for now; they remain available
as optional source builds. Older archives retain their tagged recipe.
Untracked genome caches and local analysis files are never included.

```bash
./bin/gentle_cli --state /path/to/project.gentle.json capabilities
./bin/gentle_mcp --help
```

`VERSION` and `REVISION` identify the source tag and full commit. The companion
GitHub release metadata records the archive's SHA-256 and size; `SHA256SUMS`
inside the archive covers its regular files. These are integrity/provenance
records, not a code signature. Runtime command smoke checks are not evidence
of live GUI or scientific acceptance.
