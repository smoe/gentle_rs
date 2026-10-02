# Protein reverse-translation GUI evidence

Both PNGs are direct 728×476 captures of the native `Protein Evidence` window
on Linux/X11 with `test_files/tp73.project.gentle.json` open.

- `01-protein-evidence-before-fetch.png`: upper-window orientation before online
  mutation; SHA-256
  `2b05cf0fb13685f91c42315fe3f2c11d9ec017b26de0265bd88a05820d9f9014`
- `02-ensembl-fetch-provider-failure.png`: the dedicated Ensembl subsection
  after its request failed at 45.1 seconds; SHA-256
  `d4771f69182e3674fb63bec2a2fcd414e6482cd1f0eca4513b52b042da482c06`
- GENtle source base: `4e7adf6e9dcba577892b95439af3ec749b272c3d`
- review-branch base at capture: `3be2960210ab02069de4d675b50c7002fd814915`;
  the second capture also includes the then-uncommitted ingress-worker patch
- environment: isolated 1920×1080 Xvfb display; captured child window only

These images prove layout and the observed failure only. They do not prove
provider, import, reverse-translation or lineage success. Initial coordinate
automation twice clicked the upper generic UniProt `Fetch` button with the ENSP
identifier and exposed this separate worker-stack failure:

```text
thread '<unknown>' has overflowed its stack
fatal runtime error: stack overflow, aborting
```

The review patch gives that GUI ingress worker an explicit 16 MiB stack. With
the patch, the correctly targeted `Fetch Ensembl` action stayed alive and then
reported that it could not send the lookup request after 45.1 seconds. Direct
Ensembl REST checks returned the release-116 lookup and 806-aa sequence. GENtle
CLI fetch attempts returned HTTP 500 errors for lookup/sequence and a later
attempt exceeded 30 seconds, so no successful provider result was retained.
The focused reverse-translation engine test passes with
`RUST_MIN_STACK=16777216`. Recapture result and lineage views only after the
provider path succeeds.
