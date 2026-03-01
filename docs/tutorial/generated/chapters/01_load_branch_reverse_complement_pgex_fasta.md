# Load FASTA, branch, and reverse-complement

- Chapter id: `load_branch_reverse_complement_pgex_fasta`
- Tier: `core`
- Example id: `load_branch_reverse_complement_pgex_fasta`
- Source example: `docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Start with a deterministic local sequence operation pipeline.

## Run Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json'
```

## Checkpoints

- Workflow executes without warnings or errors.
- Derived sequence IDs include a branch and reverse-complement product.

## Retained Outputs

- None for this chapter.

## Canonical Workflow JSON

```json
{
  "run_id": "example_load_branch_reverse_complement_pgex_fasta",
  "ops": [
    {
      "LoadFile": {
        "path": "test_files/pGEX_3X.fa",
        "as_id": "pgex_fasta"
      }
    },
    {
      "Branch": {
        "input": "pgex_fasta",
        "output_id": "pgex_fasta_branch"
      }
    },
    {
      "ReverseComplement": {
        "input": "pgex_fasta_branch",
        "output_id": "pgex_fasta_branch_rc"
      }
    }
  ]
}
```
