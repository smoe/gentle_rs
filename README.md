# GENtle

GENtle is a DNA and cloning workbench reimplemented for both interactive use and automation.

The cloning process is presented as a workflow. Every biotechnical operation is
mapped to an in silico representation of the same.

Protocols can explain themselves through graphical templates. The same engine
that executes a workflow can also render a deterministic cartoon of the
underlying molecular events.

![GENtle Gibson protocol cartoon](docs/figures/gibson_two_fragment_protocol_cartoon.svg)

Built-in conceptual Gibson Assembly strip: two input fragments, 5' chew-back,
annealing across the homologous overlap, polymerase fill-in, and ligase
sealing. This is meant to show how protocol templates can become readable,
factual visual companions to executable GENtle workflows.

The README figure was generated with:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  gibson.two_fragment \
  docs/figures/gibson_two_fragment_protocol_cartoon.svg
```

To start from a canonical editable protocol-cartoon template instead of the
built-in render directly, use:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon template-export \
  gibson.two_fragment \
  gibson.two_fragment.template.json
```

Then tweak the exported JSON and render it with
`protocol-cartoon render-template-svg ...` or apply overrides with
`protocol-cartoon render-with-bindings ...`.

The protocol-cartoon command surface intentionally stays canonical under
`protocol-cartoon ...` so scripted and AI-guided use does not need to choose
between overlapping alias names.

## Principles

- One engine, many interfaces: GUI, CLI, JavaScript, and Lua all use the same core logic.
- Provenance by default: derived results should be traceable and replayable.
- Structured contracts: operations, results, and errors should be machine-readable.
- Thin adapters: biology logic lives in the engine, not in frontend-specific code.

## Documentation

- Architecture: [`docs/architecture.md`](docs/architecture.md)
- Roadmap: [`docs/roadmap.md`](docs/roadmap.md)
- Protocol: [`docs/protocol.md`](docs/protocol.md)
- GUI manual: [`docs/gui.md`](docs/gui.md)
- CLI manual: [`docs/cli.md`](docs/cli.md)
- Acknowledgements: [`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md)
