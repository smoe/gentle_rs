# Digest -> Ligation -> ExtractRegion minimal slice

- Chapter id: `digest_ligation_extract_region_minimal`
- Tier: `core`
- Example id: `digest_ligation_extract_region_minimal`
- Source example: `docs/examples/workflows/digest_ligation_extract_region_minimal.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Minimal end-to-end cloning operation chain used as a regression slice.

## Run Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/digest_ligation_extract_region_minimal.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/digest_ligation_extract_region_minimal.json'
```

## Checkpoints

- Digest creates at least two fragments for ligation input.
- Ligation creates deterministic output IDs.
- ExtractRegion creates `lig_extract`.

## Retained Outputs

- None for this chapter.

## Canonical Workflow JSON

```json
{
  "run_id": "example_digest_ligation_extract_region_minimal",
  "ops": [
    {
      "LoadFile": {
        "path": "test_files/pGEX_3X.fa",
        "as_id": "pgex_fasta"
      }
    },
    {
      "Digest": {
        "input": "pgex_fasta",
        "enzymes": [
          "BamHI",
          "EcoRI"
        ],
        "output_prefix": "d"
      }
    },
    {
      "Ligation": {
        "inputs": [
          "pgex_fasta",
          "pgex_fasta"
        ],
        "circularize_if_possible": false,
        "output_id": null,
        "protocol": "Blunt",
        "output_prefix": "lig",
        "unique": false
      }
    },
    {
      "ExtractRegion": {
        "input": "lig_1",
        "from": 0,
        "to": 120,
        "output_id": "lig_extract"
      }
    }
  ]
}
```
