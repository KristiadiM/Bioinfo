#!/usr/bin/env python3
"""
Build a small SYNTHETIC bacterial genome pair for pipeline testing/demo purposes.

This is NOT real biological data. It exists so the whole Snakemake pipeline
can be run end-to-end, quickly and offline, without needing an SRA/ENA
download. Swap in real data by editing config/samples.tsv and
config/config.yaml -- see README.md, section "Using real data".

We generate:
  1. resources/reference_genome.fasta
       A "wild-type / reference-strain" chromosome carrying a susceptible
       allele of a fluoroquinolone-target-like gene (gyrSYN).
  2. resources/sample_genome.fasta (used only to simulate reads -- the
     pipeline itself never sees this file, exactly as in real life you
     don't know the isolate's genome before sequencing it)
       The "clinical isolate" genome: same chromosome backbone as the
       reference, PLUS two horizontally-acquired resistance-gene-like
       inserts (blaSYN, mefSYN), PLUS a single non-synonymous point
       mutation in gyrSYN mimicking a target-mutation resistance allele.
  3. resources/amr_marker_db.fasta
       A small stand-in for a curated AMR gene database (CARD/ResFinder/
       NCBI AMRFinderPlus in real use), containing wild-type blaSYN,
       mefSYN and a decoy gene (vanSYN) that is deliberately absent from
       the isolate, to prove the screening step correctly reports absence.

Design rationale (why two detection strategies are demoed):
  - blaSYN/mefSYN are inserted sequence NOT present in the reference
    chromosome, so read-mapping-based variant calling cannot see them --
    exactly like real plasmid-borne acquired resistance genes. They are
    only recoverable via assembly + homology screening (BLAST here;
    ABRicate/RGI/AMRFinderPlus against CARD/ResFinder in production).
  - gyrSYN's point mutation IS on the shared backbone, so it is only
    recoverable via read mapping + variant calling -- exactly like real
    chromosomal target-mutation resistance (e.g. gyrA/rpoB-type SNPs).
  This mirrors why real AMR genomics pipelines run both approaches.
"""
import random
from pathlib import Path

RESOURCES = Path(__file__).resolve().parent.parent / "resources"
SIM_GENOMES = RESOURCES / "sim_genomes"
RESOURCES.mkdir(exist_ok=True)
SIM_GENOMES.mkdir(exist_ok=True)

random.seed(42)
BASES = "ACGT"


def rand_seq(n, seed=None):
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(n))


def write_fasta(path, records):
    """records: list of (header, sequence)"""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def mutate_codon(seq, codon_start, new_codon):
    return seq[:codon_start] + new_codon + seq[codon_start + 3:]


def main():
    # --- backbone shared by reference and isolate ---
    backbone = rand_seq(12000, seed=1)

    # --- chromosomal target gene (present in both, WT vs mutant allele) ---
    gyr_wt = rand_seq(600, seed=2)
    # introduce a defined non-synonymous SNP at codon 28 (0-based pos 81),
    # analogous in spirit (not sequence) to real gyrA target mutations
    gyr_mut = mutate_codon(gyr_wt, 81, "TTG")

    # --- horizontally acquired resistance-gene-like inserts (isolate only) ---
    bla_gene = rand_seq(900, seed=3)   # beta-lactamase-like
    mef_gene = rand_seq(750, seed=4)   # efflux-pump-like

    # --- decoy gene: present in the marker DB, absent from the isolate ---
    van_gene = rand_seq(500, seed=5)   # glycopeptide-resistance-like

    # positions (bp into backbone) where genes sit in the isolate
    p1, p2, p3 = 3000, 6500, 9500

    # reference genome: backbone + WT gyr gene only (no acquired genes)
    reference_seq = (
        backbone[:p3] + gyr_wt + backbone[p3:]
    )

    # isolate 1 ("resistant"): backbone + bla insert + mef insert + mutant gyr gene
    isolate01_seq = (
        backbone[:p1] + bla_gene
        + backbone[p1:p2] + mef_gene
        + backbone[p2:p3] + gyr_mut
        + backbone[p3:]
    )

    # isolate 2 ("susceptible control"): same backbone, WT gyr, no acquired genes
    isolate02_seq = reference_seq

    write_fasta(
        RESOURCES / "reference_genome.fasta",
        [("chr1_reference synthetic reference strain, susceptible gyrSYN allele",
          reference_seq)],
    )
    write_fasta(
        SIM_GENOMES / "isolate01_resistant.fasta",
        [("chr1_isolate01 synthetic MDR-like clinical isolate (read simulation only)",
          isolate01_seq)],
    )
    write_fasta(
        SIM_GENOMES / "isolate02_susceptible.fasta",
        [("chr1_isolate02 synthetic susceptible control isolate (read simulation only)",
          isolate02_seq)],
    )
    write_fasta(
        RESOURCES / "amr_marker_db.fasta",
        [
            ("blaSYN_beta-lactamase-like|acquired|expected:PRESENT", bla_gene),
            ("mefSYN_efflux-pump-like|acquired|expected:PRESENT", mef_gene),
            ("vanSYN_glycopeptide-res-like|acquired|expected:ABSENT (decoy)", van_gene),
        ],
    )

    print(f"[build_synthetic_genomes] reference:  {len(reference_seq)} bp")
    print(f"[build_synthetic_genomes] isolate01:  {len(isolate01_seq)} bp (bla+, mef+, gyr-mutant)")
    print(f"[build_synthetic_genomes] isolate02:  {len(isolate02_seq)} bp (bla-, mef-, gyr-WT)")
    print("[build_synthetic_genomes] wrote reference_genome.fasta, "
          "sim_genomes/isolate0{1,2}_*.fasta, amr_marker_db.fasta")


if __name__ == "__main__":
    main()
