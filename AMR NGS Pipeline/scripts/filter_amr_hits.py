#!/usr/bin/env python3
"""
Turn a raw BLASTN outfmt-6 table (query = assembly contigs, subject = AMR
marker gene database) into a clean per-gene presence/absence report,
applying identity and coverage thresholds.

This is a minimal stand-in for what ABRicate / AMRFinderPlus / RGI do
against CARD or ResFinder in a production pipeline -- same idea
(homology search + thresholding), simplified so it runs with nothing
but BLAST+ and pandas.

outfmt 6 columns used here:
    qseqid sseqid pident length mismatch gapopen
    qstart qend sstart send evalue bitscore slen
"""
import argparse
import pandas as pd
from Bio import SeqIO

COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blast_tsv", help="raw blastn -outfmt 6 output (with slen)")
    ap.add_argument("amr_db_fasta", help="AMR marker gene FASTA used as BLAST subject db")
    ap.add_argument("out_tsv", help="output summary TSV")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-identity", type=float, default=90.0)
    ap.add_argument("--min-coverage", type=float, default=80.0)
    args = ap.parse_args()

    all_genes = [rec.id for rec in SeqIO.parse(args.amr_db_fasta, "fasta")]

    try:
        df = pd.read_csv(args.blast_tsv, sep="\t", names=COLUMNS)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=COLUMNS)

    if not df.empty:
        df["coverage_pct"] = 100.0 * df["length"] / df["slen"]
        passed = df[
            (df["pident"] >= args.min_identity)
            & (df["coverage_pct"] >= args.min_coverage)
        ]
        # keep best (highest bitscore) hit per subject gene
        best = (
            passed.sort_values("bitscore", ascending=False)
            .drop_duplicates(subset="sseqid", keep="first")
            .set_index("sseqid")
        )
    else:
        best = pd.DataFrame(columns=["pident", "coverage_pct", "qseqid"]).set_index(
            pd.Index([], name="sseqid")
        )

    rows = []
    for gene in all_genes:
        if gene in best.index:
            hit = best.loc[gene]
            rows.append({
                "sample": args.sample,
                "marker_gene": gene,
                "status": "PRESENT",
                "pct_identity": round(float(hit["pident"]), 2),
                "pct_coverage": round(float(hit["coverage_pct"]), 2),
                "contig": hit["qseqid"],
            })
        else:
            rows.append({
                "sample": args.sample,
                "marker_gene": gene,
                "status": "absent",
                "pct_identity": None,
                "pct_coverage": None,
                "contig": None,
            })

    out = pd.DataFrame(rows)
    out.to_csv(args.out_tsv, sep="\t", index=False)
    n_present = (out["status"] == "PRESENT").sum()
    print(f"[filter_amr_hits] {args.sample}: {n_present}/{len(all_genes)} "
          f"marker genes detected (>= {args.min_identity}% id, "
          f">= {args.min_coverage}% cov)")


if __name__ == "__main__":
    main()
