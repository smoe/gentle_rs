import random, gzip, os

OUT = "/Users/u005069/GitHub/gentle_rs/test_files/fixtures/target_rescue"
os.makedirs(OUT, exist_ok=True)
K = 11
rng = random.Random(20260704)  # fixed seed for reproducibility (see PROVENANCE.md)

def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

def kmers(s, k=K):
    return {s[i:i+k] for i in range(0, len(s)-k+1)}

def make_disjoint(existing, n):
    """Draw a random sequence whose k-mers are disjoint from every set in `existing`."""
    for _ in range(10000):
        s = rand_seq(n)
        ks = kmers(s)
        # also require the sequence has no internal repeated k-mer (keeps counts predictable)
        if len(ks) != len(s) - K + 1:
            continue
        if all(ks.isdisjoint(e) for e in existing):
            return s, ks
    raise RuntimeError("could not build disjoint sequence")

sa, ka = make_disjoint([], 80)                 # GENEA transcript
sb, kb = make_disjoint([ka], 80)               # GENEB transcript
sc, kc = make_disjoint([ka, kb], 80)           # GENEC transcript (present, not requested)
snone, knone = make_disjoint([ka, kb, kc], 60) # unrelated read, hits nothing

# Reads derived from the transcripts
r_geneA = sa[10:54]                 # 44 bp, pure GENEA -> hits GENEA only
r_geneB = sb[8:56]                  # 48 bp, pure GENEB -> hits GENEB only
r_ambig = sa[0:30] + sb[0:30]       # 60 bp, spans both -> ambiguous, hits GENEA and GENEB
r_none  = snone                     # hits nothing

def transcript_header(tid, gene):
    return (f">{tid} cdna chromosome:GRCh38:1:1:80:1 gene:ENSG_{gene}.1 "
            f"gene_biotype:protein_coding transcript_biotype:protein_coding "
            f"gene_symbol:{gene} description:synthetic {gene} transcript")

with open(f"{OUT}/mini_transcripts.fasta", "w") as fh:
    for tid, gene, seq in [("ENSTA1.1","GENEA",sa),("ENSTB1.1","GENEB",sb),("ENSTC1.1","GENEC",sc)]:
        fh.write(transcript_header(tid, gene) + "\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")

reads = [
    ("readA_geneA", r_geneA),
    ("readB_geneB", r_geneB),
    ("readAB_ambiguous", r_ambig),
    ("readN_nohit", r_none),
    ("readA2_geneA", r_geneA),  # duplicate GENEA read to exercise counts > 1
]

def write_fastq(path):
    lines = []
    for rid, seq in reads:
        lines.append(f"@{rid}\n{seq}\n+\n{'I'*len(seq)}\n")
    data = "".join(lines).encode()
    with open(path, "wb") as fh:
        fh.write(data)
    with gzip.open(path + ".gz", "wb") as fh:
        fh.write(data)

write_fastq(f"{OUT}/mini_reads.fastq")

# Read-ID allowlist (a subset): only the two GENEA reads + the ambiguous read
with open(f"{OUT}/read_allowlist.txt", "w") as fh:
    fh.write("readA_geneA\nreadAB_ambiguous\nreadA2_geneA\n")

# Salmon unmapped_names.txt: <name> <code> per line (u = both mates unmapped)
with open(f"{OUT}/salmon_unmapped_names.txt", "w") as fh:
    for rid in ["readA_geneA", "readN_nohit", "readAB_ambiguous", "readA2_geneA"]:
        fh.write(f"{rid} u\n")

# Salmon --writeMappings SAM-like file. readB maps to the GENEB target transcript;
# readA maps to GENEC (present in FASTA but NOT a requested target) -> excluded from target set.
with open(f"{OUT}/salmon_mappings.sam", "w") as fh:
    fh.write("@HD\tVN:1.0\tSO:unsorted\n")
    fh.write("@SQ\tSN:ENSTB1.1\tLN:80\n")
    fh.write("@SQ\tSN:ENSTC1.1\tLN:80\n")
    fh.write(f"readB_geneB\t0\tENSTB1.1\t9\t255\t48M\t*\t0\t0\t{r_geneB}\t*\n")
    fh.write(f"readA_geneA\t0\tENSTC1.1\t1\t255\t44M\t*\t0\t0\t{r_geneA}\t*\n")

print("GENEA/GENEB 11-mer overlap:", len(ka & kb))
print("readA hits A:", not kmers(r_geneA).isdisjoint(ka), "hits B:", not kmers(r_geneA).isdisjoint(kb))
print("readB hits A:", not kmers(r_geneB).isdisjoint(ka), "hits B:", not kmers(r_geneB).isdisjoint(kb))
print("readAB hits A:", not kmers(r_ambig).isdisjoint(ka), "hits B:", not kmers(r_ambig).isdisjoint(kb))
print("readN hits A:", not kmers(r_none).isdisjoint(ka), "hits B:", not kmers(r_none).isdisjoint(kb))
print("wrote fixtures to", OUT)
