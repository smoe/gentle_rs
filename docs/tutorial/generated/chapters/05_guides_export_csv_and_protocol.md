# Guide oligo export (CSV + protocol)

- Chapter id: `guides_export_csv_and_protocol`
- Tier: `advanced`
- Example id: `guides_export_csv_and_protocol`
- Source example: `docs/examples/workflows/guides_export_csv_and_protocol.json`
- Example test_mode: `skip`
- Executed during generation: `yes`

Retains one representative oligo CSV and one protocol text artifact for GitHub inspection.

## Run Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/guides_export_csv_and_protocol.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/guides_export_csv_and_protocol.json'
```

## Checkpoints

- CSV export file exists and contains guide rows.
- Protocol text export exists and contains checklist content.

## Retained Outputs

- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.csv)
- [`artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt`](../artifacts/guides_export_csv_and_protocol/exports/demo_guides.protocol.txt)

## Canonical Workflow JSON

```json
{
  "run_id": "example_guides_export_csv_and_protocol",
  "ops": [
    {
      "UpsertGuideSet": {
        "guide_set_id": "demo_guides",
        "guides": [
          {
            "guide_id": "demo_1",
            "seq_id": "target_demo",
            "start_0based": 10,
            "end_0based_exclusive": 30,
            "strand": "+",
            "protospacer": "GACCTGTTGACGATGTTCCA",
            "pam": "AGG",
            "nuclease": "SpCas9",
            "cut_offset_from_protospacer_start": 17,
            "rank": 1
          }
        ]
      }
    },
    {
      "GenerateGuideOligos": {
        "guide_set_id": "demo_guides",
        "template_id": "lenti_bsmbi_u6_default",
        "apply_5prime_g_extension": true,
        "output_oligo_set_id": "demo_lenti",
        "passed_only": false
      }
    },
    {
      "ExportGuideOligos": {
        "guide_set_id": "demo_guides",
        "oligo_set_id": "demo_lenti",
        "format": "csv_table",
        "path": "exports/demo_guides.csv",
        "plate_format": null
      }
    },
    {
      "ExportGuideProtocolText": {
        "guide_set_id": "demo_guides",
        "oligo_set_id": "demo_lenti",
        "path": "exports/demo_guides.protocol.txt",
        "include_qc_checklist": true
      }
    }
  ]
}
```
