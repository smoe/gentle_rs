# Prepare a reference genome cache (online)

- Chapter id: `prepare_reference_genome_online`
- Tier: `online`
- Example id: `prepare_reference_genome_online`
- Source example: `docs/examples/workflows/prepare_reference_genome_online.json`
- Example test_mode: `online`
- Executed during generation: `no`
- Execution note: set `GENTLE_TEST_ONLINE=1` before `tutorial-generate` to execute this chapter.

Network/tooling-dependent chapter, executed only with explicit online opt-in.

## Run Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/prepare_reference_genome_online.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/prepare_reference_genome_online.json'
```

## Checkpoints

- Genome preparation runs only when GENTLE_TEST_ONLINE is enabled.
- Offline generation still emits the chapter with execution status noted.

## Retained Outputs

- None for this chapter.

## Canonical Workflow JSON

```json
{
  "run_id": "example_prepare_reference_genome_online",
  "ops": [
    {
      "PrepareGenome": {
        "genome_id": "Human GRCh38 Ensembl 116",
        "catalog_path": "assets/genomes.json",
        "cache_dir": "data/genomes",
        "timeout_seconds": null
      }
    }
  ]
}
```
