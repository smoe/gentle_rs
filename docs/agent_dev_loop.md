# Agent Development Loop

Use in-tree binaries while iterating on adapter capability and error contracts.

```sh
cargo run -q --bin gentle_cli -- capabilities
```

The CLI payload includes `capability_registry`, which is projected from
`gentle_protocol`.

```sh
printf 'Content-Length: 58\r\n\r\n{"jsonrpc":"2.0","id":1,"method":"tools/list","params":{}}' \
  | cargo run -q --bin gentle_mcp --
```

MCP `tools/list` projects tool names, titles, and descriptions from the same
protocol registry while preserving MCP-shaped input schemas.
