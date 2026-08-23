#!/usr/bin/env python3
"""
Minimal paired-end Illumina-style read simulator.

Stands in for a sequencer / SRA download so the pipeline can be tested
offline and quickly. For real use, this rule is simply swapped out for a
download rule (e.g. `prefetch` + `fasterq-dump`, or `enaDataGet`) that
pulls FASTQ files for a real run accession -- the rest of the pipeline
(QC, trimming, assembly, mapping, variant calling, AMR screening) is
unaffected either way, which is the point of separating pipeline stages
into Snakemake rules with well-defined file inputs/outputs.

Usage:
    simulate_reads.py <genome.fasta> <out_R1.fastq.gz> <out_R2.fastq.gz> \
        --coverage 40 --read-length 150 --frag-mean 400 --error-rate 0.01 --seed 7
"""
import argparse
import gzip
import random
from Bio import SeqIO

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def add_errors(seq, error_rate, rng):
    bases = list(seq)
    for i, b in enumerate(bases):
        if rng.random() < error_rate:
            bases[i] = rng.choice([x for x in "ACGT" if x != b])
    return "".join(bases)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genome_fasta")
    ap.add_argument("out_r1")
    ap.add_argument("out_r2")
    ap.add_argument("--coverage", type=float, default=40)
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--frag-mean", type=int, default=400)
    ap.add_argument("--frag-sd", type=int, default=40)
    ap.add_argument("--error-rate", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    record = next(SeqIO.parse(args.genome_fasta, "fasta"))
    genome = str(record.seq).upper()
    glen = len(genome)

    n_fragments = int((glen * args.coverage) / (2 * args.read_length))

    with gzip.open(args.out_r1, "wt") as r1, gzip.open(args.out_r2, "wt") as r2:
        for i in range(n_fragments):
            frag_len = max(
                args.read_length * 2,
                int(rng.gauss(args.frag_mean, args.frag_sd)),
            )
            if frag_len >= glen:
                continue
            start = rng.randint(0, glen - frag_len)
            fragment = genome[start:start + frag_len]

            fwd = fragment[: args.read_length]
            rev = revcomp(fragment[-args.read_length:])

            fwd = add_errors(fwd, args.error_rate, rng)
            rev = add_errors(rev, args.error_rate, rng)

            qual = "I" * args.read_length  # flat high-quality Phred score
            read_name = f"sim_read_{i}"

            r1.write(f"@{read_name}/1\n{fwd}\n+\n{qual}\n")
            r2.write(f"@{read_name}/2\n{rev}\n+\n{qual}\n")

    print(f"[simulate_reads] wrote {n_fragments} read pairs "
          f"(~{args.coverage}x coverage of {glen} bp genome)")


if __name__ == "__main__":
    main()
