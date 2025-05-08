"""sgRNA Predictor"""

import csv
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
import numpy as np

# Set your email for NCBI Entrez (required)
Entrez.email = "your_email@example.com"

# Stub: Generate a fake score (0–1) — replace this with a real model later
def get_doench_score(context_30mer):
    if len(context_30mer) != 30:
        return None
    return float(np.random.uniform(0, 1))  # Random score for demonstration

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

def has_poly_t(seq):
    return 'TTTT' in seq

def find_sgRNAs(dna_sequence, genome=None, gc_range=(40, 60), min_score=0.2, allow_duplicates=True):
    dna_sequence = dna_sequence.upper()
    candidates = []

    def is_valid_spacer(spacer):
        return (
            len(spacer) == 20 and
            gc_range[0] <= gc_content(spacer) <= gc_range[1] and
            not has_poly_t(spacer)
        )

    # Forward strand scanning
    for i in range(4, len(dna_sequence) - 26):
        context = dna_sequence[i-4:i+26]
        spacer = context[4:24]
        pam = context[24:27]
        if pam[1:] == 'GG' and is_valid_spacer(spacer):
            if allow_duplicates or genome.count(spacer) == 1:
                score = get_doench_score(context)
                if score is None or score >= min_score:
                    candidates.append((spacer, pam, i, '+', round(score, 3)))

    # Reverse strand scanning
    rev_seq = str(Seq(dna_sequence).reverse_complement())
    for i in range(4, len(rev_seq) - 26):
        context = rev_seq[i-4:i+26]
        spacer = context[4:24]
        pam = context[24:27]
        rev_index = len(dna_sequence) - i - 26
        if pam[1:] == 'GG' and is_valid_spacer(spacer):
            if allow_duplicates or genome.count(spacer) == 1:
                score = get_doench_score(context)
                if score is None or score >= min_score:
                    candidates.append((spacer, pam, rev_index, '-', round(score, 3)))

    return sorted(candidates, key=lambda x: -x[4])

def save_to_csv(candidates, filename="sgRNA_candidates.csv"):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Spacer", "PAM", "Position", "Strand", "Doench Score"])
        writer.writerows(candidates)
    print(f"Saved {len(candidates)} sgRNAs to {filename}")

def plot_scores(candidates):
    labels = [f"{c[2]}:{c[3]}" for c in candidates]
    scores = [c[4] for c in candidates]

    plt.figure(figsize=(10, 5))
    plt.bar(labels, scores, color='teal')
    plt.xticks(rotation=90)
    plt.ylabel("Predicted Activity Score")
    plt.title("sgRNA Activity (Placeholder Scores)")
    plt.tight_layout()
    plt.show()

# Data fetching from NCBI
def fetch_from_ncbi(query):
    try:
        handle = Entrez.efetch(db="nucleotide", id=query, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        return str(record.seq)
    except Exception:
        print(f"'⚠️ {query}' not found as an accession. Trying as a keyword search...")
        search_handle = Entrez.esearch(db="nucleotide", term=query)
        search_results = Entrez.read(search_handle)
        ids = search_results["IdList"]
        if not ids:
            raise ValueError(f"⚠️ No results found for '{query}' on NCBI")

        print(f"Found {len(ids)} matching sequences. Using the first one.")
        handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        return str(record.seq)

# User filter preference
def get_user_filters():
    print("\n--- Filtering Options ---")
    while True:
        min_score_input = input("Enter minimum score [default: 0.5]: ")
        if not min_score_input:
            min_score = 0.5
            break
        try:
            min_score = float(min_score_input)
            if 0 <= min_score <= 1:
                break
            else:
                print("Please enter a value between 0 and 1.")
        except ValueError:
            print("⚠️ Invalid input.")

    while True:
        gc_min = input("Enter minimum GC% [default: 40]: ") or 40
        gc_max = input("Enter maximum GC% [default: 60]: ") or 60
        try:
            gc_min, gc_max = int(gc_min), int(gc_max)
            if 0 <= gc_min < gc_max <= 100:
                gc_range = (gc_min, gc_max)
                break
            else:
                print("GC% must be between 0 and 100, and min < max.")
        except ValueError:
            print("⚠️ Invalid input.")

    dup_filter = input("Remove sgRNAs with possible off-targets? (y/n) [default: n]: ").lower() in ['y', 'yes']

    return {
        'min_score': min_score,
        'gc_range': gc_range,
        'allow_duplicates': not dup_filter
    }

# Execution block
if __name__ == "__main__":
    print("CRISPR sgRNA Designer")
    user_input = input("Enter NCBI Accession (e.g., NM_001301718) or Gene Name/Keyword: ").strip()

    try:
        sequence = fetch_from_ncbi(user_input)
        filters = get_user_filters()
        output_file = "sgRNA_candidates.csv"

        sgRNAs = find_sgRNAs(
            sequence,
            genome=sequence,
            min_score=filters['min_score'],
            gc_range=filters['gc_range'],
            allow_duplicates=filters['allow_duplicates']
        )
        save_to_csv(sgRNAs, output_file)
        plot_scores(sgRNAs)

    except Exception as e:
        print(f"[ERROR] {e}")