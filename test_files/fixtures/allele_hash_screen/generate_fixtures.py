#!/usr/bin/env python3
"""Regenerate the deterministic allele-hash-screen fixture files."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent

TX1 = "ATGACCAAGTTGACCTGACCGTACCGATGCTAGCTTACGATCGTACGATCTAGCTAGGCTA"
TX2 = "GCTTAGCAAGTCGATCCGTAAGCTTAGGCTAACCGTATGCAACTTGACCGATTCGGAATCCG"


def write(path: str, text: str) -> None:
    (ROOT / path).write_text(text, encoding="utf-8")


def main() -> None:
    hap1 = list(TX1)
    hap2 = list(TX1)
    hap2[19] = "T"
    hap1[44] = "G"
    hap1 = "".join(hap1)
    hap2 = "".join(hap2)
    reads = {
        "fus_hap1_alt_v2": hap1[40:49],
        "fus_hap2_alt_v1": hap2[15:24],
        "fus_shared_region": TX1[0:9],
        "fus_off_target": "GGGGGGGGG",
        "fus_short_uninformative": "ACGT",
    }
    write(
        "fus_transcripts.fa",
        f">FUS_TX1 gene_symbol:FUS synthetic:allele_hash_screen\n{TX1}\n"
        f">FUS_TX2 gene_symbol:FUS synthetic:allele_hash_screen\n{TX2}\n",
    )
    write(
        "fus_variants.tsv",
        "variant_id\ttranscript_id\tcdna_pos_1based\tref\talt\tgenotype\tphase_set\n"
        "v1\tFUS_TX1\t20\tC\tT\t0|1\tps_fus_demo\n"
        "v2\tFUS_TX1\t45\tA\tG\t1|0\tps_fus_demo\n",
    )
    fastq = "".join(
        f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n" for read_id, seq in reads.items()
    )
    write("fus_reads.fastq", fastq)
    write("read_allowlist.txt", "".join(f"{read_id}\n" for read_id in reads))
    write(
        "salmon_unmapped_names.txt",
        "".join(
            f"{read_id} u\n"
            for read_id in [
                "fus_hap1_alt_v2",
                "fus_hap2_alt_v1",
                "fus_shared_region",
                "fus_off_target",
            ]
        ),
    )
    write(
        "salmon_mappings.sam",
        "@HD\tVN:1.0\tSO:unsorted\n"
        f"@SQ\tSN:FUS_TX1\tLN:{len(TX1)}\n"
        "fus_short_uninformative\t0\tFUS_TX1\t1\t255\t4M\t*\t0\t0\tACGT\t*\n",
    )


if __name__ == "__main__":
    main()
