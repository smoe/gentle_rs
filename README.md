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
sealing. It shows how executable GENtle workflows can also carry readable,
factual visual explanations of the underlying molecular events.

This README figure is rendered directly by the shared protocol-cartoon engine:

```sh
cargo run --quiet --bin gentle_cli -- \
  protocol-cartoon render-svg \
  gibson.two_fragment \
  docs/figures/gibson_two_fragment_protocol_cartoon.svg
```

Canonical protocol-cartoon templates start from readable defaults and can then
be adapted to the actual parameters of a concrete cloning routine.

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
- Tutorial hub: [`docs/tutorial/generated/README.md`](docs/tutorial/generated/README.md)
- Additional walkthroughs: [`docs/tutorial/`](docs/tutorial/) and [`docs/agent_interfaces_tutorial.md`](docs/agent_interfaces_tutorial.md)
- Acknowledgements: [`ACKNOWLEDGEMENTS.md`](ACKNOWLEDGEMENTS.md)
