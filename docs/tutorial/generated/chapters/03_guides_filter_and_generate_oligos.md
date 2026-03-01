# Guide practical filtering and oligo generation

- Chapter id: `guides_filter_and_generate_oligos`
- Tier: `core`
- Example id: `guides_filter_and_generate_oligos`
- Source example: `docs/examples/workflows/guides_filter_and_generate_oligos.json`
- Example test_mode: `always`
- Executed during generation: `yes`

Covers guide-set creation, practical filtering, and oligo generation.

## Run Commands

```bash
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/guides_filter_and_generate_oligos.json
cargo run --bin gentle_cli -- shell 'workflow @docs/examples/workflows/guides_filter_and_generate_oligos.json'
```

## Checkpoints

- Guide set and passed guide set are present in metadata.
- Oligo set generation succeeds for passed guides.

## Retained Outputs

- None for this chapter.

## Canonical Workflow JSON

```json
{
  "run_id": "example_guides_filter_and_generate_oligos",
  "ops": [
    {
      "UpsertGuideSet": {
        "guide_set_id": "tp73_guides",
        "guides": [
          {
            "guide_id": "g1",
            "seq_id": "tp73",
            "start_0based": 100,
            "end_0based_exclusive": 120,
            "strand": "+",
            "protospacer": "GACCTGTTGACGATGTTCCA",
            "pam": "AGG",
            "nuclease": "SpCas9",
            "cut_offset_from_protospacer_start": 17,
            "rank": 1
          },
          {
            "guide_id": "g2",
            "seq_id": "tp73",
            "start_0based": 220,
            "end_0based_exclusive": 240,
            "strand": "+",
            "protospacer": "TTTTGCCATGTTGACCTGAA",
            "pam": "TGG",
            "nuclease": "SpCas9",
            "cut_offset_from_protospacer_start": 17,
            "rank": 2
          },
          {
            "guide_id": "g3",
            "seq_id": "tp73",
            "start_0based": 340,
            "end_0based_exclusive": 360,
            "strand": "-",
            "protospacer": "GGTACCGATGTTGCCAGTAA",
            "pam": "CGG",
            "nuclease": "SpCas9",
            "cut_offset_from_protospacer_start": 17,
            "rank": 3
          }
        ]
      }
    },
    {
      "FilterGuidesPractical": {
        "guide_set_id": "tp73_guides",
        "config": {
          "gc_min": 0.3,
          "gc_max": 0.7,
          "max_homopolymer_run": 4,
          "max_homopolymer_run_per_base": {},
          "reject_ambiguous_bases": true,
          "avoid_u6_terminator_tttt": true,
          "u6_terminator_window": "spacer_plus_tail",
          "max_dinucleotide_repeat_units": null,
          "forbidden_motifs": [],
          "required_5prime_base": "G",
          "allow_5prime_g_extension": true
        },
        "output_guide_set_id": "tp73_guides_pass"
      }
    },
    {
      "GenerateGuideOligos": {
        "guide_set_id": "tp73_guides",
        "template_id": "lenti_bsmbi_u6_default",
        "apply_5prime_g_extension": true,
        "output_oligo_set_id": "tp73_lenti",
        "passed_only": true
      }
    }
  ]
}
```
