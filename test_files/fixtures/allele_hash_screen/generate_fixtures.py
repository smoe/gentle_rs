#!/usr/bin/env python3
"""Regenerate the deterministic allele-hash-screen fixture files."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent

TX1 = "ATGACCAAGTTGACCTGACCGTACCGATGCTAGCTTACGATCGTACGATCTAGCTAGGCTA"
TX2 = "GCTTAGCAAGTCGATCCGTAAGCTTAGGCTAACCGTATGCAACTTGACCGATTCGGAATCCG"
REP_TRANSCRIPTS = {
    "FUS_REP_TX1": TX1,
    "FUS_REP_TX2": TX2,
    "FUS_REP_TX3": "CGGAAGGCCACCCGACGGGTGGCAAGGCCGAAAGGCCACATCTACTGAATACGTCGAGGCAAAGTG",
    "FUS_REP_TX4": "ATCCATTAACAGCCAGAGGTATTCCATGGAAGTTGCGCCGCGACATGCGTCTTTATTACCTTTACG",
}
ALT_BASE = {"A": "C", "C": "G", "G": "T", "T": "A"}


def write(path: str, text: str) -> None:
    (ROOT / path).write_text(text, encoding="utf-8")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def canonical(sequence: str) -> str:
    return min(sequence, reverse_complement(sequence))


def marker_read(sequence: str, pos_1based: int, alt: str) -> str:
    mutated = list(sequence)
    mutated[pos_1based - 1] = alt
    start = pos_1based - 1 - 4
    return "".join(mutated[start : start + 9])


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

    representation_variants = []
    representation_reads = []
    marker_kmers = []
    requested_counts = {
        "FUS_REP_TX1": (12, 2),
        "FUS_REP_TX2": (5, 5),
        "FUS_REP_TX3": (2, 1),
        "FUS_REP_TX4": (0, 0),
    }
    for transcript_id, sequence in REP_TRANSCRIPTS.items():
        hap1_pos = 20
        hap2_pos = 45
        hap1_alt = ALT_BASE[sequence[hap1_pos - 1]]
        hap2_alt = ALT_BASE[sequence[hap2_pos - 1]]
        phase_set = f"ps_{transcript_id.lower()}"
        representation_variants.extend(
            [
                f"{transcript_id}_hap1\t{transcript_id}\t{hap1_pos}\t{sequence[hap1_pos - 1]}\t{hap1_alt}\t1|0\t{phase_set}\n",
                f"{transcript_id}_hap2\t{transcript_id}\t{hap2_pos}\t{sequence[hap2_pos - 1]}\t{hap2_alt}\t0|1\t{phase_set}\n",
            ]
        )
        hap1_marker = marker_read(sequence, hap1_pos, hap1_alt)
        hap2_marker = marker_read(sequence, hap2_pos, hap2_alt)
        marker_kmers.extend([canonical(hap1_marker), canonical(hap2_marker)])
        hap1_count, hap2_count = requested_counts[transcript_id]
        representation_reads.extend(
            (f"{transcript_id}_hap1_{index:02d}", hap1_marker)
            for index in range(1, hap1_count + 1)
        )
        representation_reads.extend(
            (f"{transcript_id}_hap2_{index:02d}", hap2_marker)
            for index in range(1, hap2_count + 1)
        )

    reference_kmers = {
        canonical(sequence[start : start + 9])
        for sequence in REP_TRANSCRIPTS.values()
        for start in range(len(sequence) - 8)
    }
    assert len(set(marker_kmers)) == len(marker_kmers)
    assert set(marker_kmers).isdisjoint(reference_kmers)
    write(
        "representation_transcripts.fa",
        "".join(
            f">{transcript_id} gene_symbol:FUS synthetic:representation_verdict\n{sequence}\n"
            for transcript_id, sequence in REP_TRANSCRIPTS.items()
        ),
    )
    write(
        "representation_variants.tsv",
        "variant_id\ttranscript_id\tcdna_pos_1based\tref\talt\tgenotype\tphase_set\n"
        + "".join(representation_variants),
    )
    write(
        "representation_reads.fastq",
        "".join(
            f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n"
            for read_id, sequence in representation_reads
        ),
    )


if __name__ == "__main__":
    main()
