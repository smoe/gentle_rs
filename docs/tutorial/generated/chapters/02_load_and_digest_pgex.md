# Load pGEX and digest with BamHI/EcoRI

- Chapter id: `load_and_digest_pgex`
- Tier: `core`
- Example id: `load_and_digest_pgex`
- Source example: `docs/examples/workflows/load_and_digest_pgex.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Demonstrates deterministic restriction digest product generation.

## Run Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_and_digest_pgex.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/load_and_digest_pgex.json'
```

## Checkpoints

- Digest operation completes and creates fragment sequence IDs.
- Fragment IDs are deterministic across repeated runs.

## Retained Outputs

- None for this chapter.

## Canonical Workflow JSON

```json
{
  "run_id": "example_load_and_digest_pgex",
  "ops": [
    {
      "LoadFile": {
        "path": "test_files/pGEX-3X.gb",
        "as_id": "pgex"
      }
    },
    {
      "Digest": {
        "input": "pgex",
        "enzymes": [
          "BamHI",
          "EcoRI"
        ],
        "output_prefix": "frag"
      }
    }
  ]
}
```
