# CRISPR sgRNA Predictor
A simple Python-based tool to predict and analyze potential sgRNA sequences for CRISPR-Cas9 applications. 
This script scans DNA sequences for suitable guide RNAs (sgRNAs), filters them based on GC content, off-target potential, and predicted activity scores, and outputs results in CSV format with a basic visualization.

## Features
- NCBI Integration : Fetch genomic sequences using an accession number or keyword search.
- sgRNA Prediction : Identifies PAM sites (NGG) and generates candidate sgRNAs.
- Filtering Options :
- Minimum Doench score (predicted activity)
- GC content range
- Remove possible off-targets via duplicate detection
- Visualization : Bar plot of predicted sgRNA activities.
- Stubbed Scoring Model : Placeholder scoring function (get_doench_score) can be replaced with a real machine learning model later.
  
## Requirements
pip install biopython numpy matplotlib

## How to Run
1. Clone or download the repository.
2. Open terminal in the project directory.
3. Run the script: python "sgRNA predictor.py"
4. Input Prompts :
  - Enter an NCBI accession (e.g., NM_001301718) or a gene name/keyword.
  - Set filtering preferences when prompted:
  - Minimum score threshold
  - GC% range
  - Off-target filtering (remove duplicates)

## Output :
A CSV file: sgRNA_candidates.csv
A bar chart showing predicted sgRNA activity scores

## How It Works
1. The script scans both forward and reverse strands for sequences preceding the PAM motif (NGG).
2. Each candidate is evaluated based on:
   - Length (must be 20 nucleotides)
   - GC content within specified bounds
   - No poly-T tracts (TTTT)
   - Predicted activity score (currently random)
3. Optional genome-wide uniqueness check (simple substring match)

📬 Contact
For questions or suggestions, reach out at:
📧 kristiadimikael@gmail.com
