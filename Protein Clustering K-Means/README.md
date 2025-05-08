# Protein Clustering with K-Means
This script fetches protein sequences from NCBI , performs K-Means clustering using k-mer feature extraction , and exports results to a CSV file. 
It also includes visualization using PCA and UMAP, and evaluates clustering performance using Silhouette Score and Davies-Bouldin Index .

## Features
- Fetch protein sequences by name
- Filter out invalid/short/duplicate sequences
- Extract features using k-mers (default: 3-mers)
- Feature selection and standardization
- Optimal K-Means clustering with evaluation metrics
- Export clustering results to CSV (protein_clusters.csv)
- Visualize clusters in 2D using PCA and UMAP

## Requirements
pip install biopython numpy pandas scikit-learn matplotlib seaborn umap-learn

## How to Run
1. Clone or download the repository.
2. Open terminal in the project directory.
3. Run the script: python "Protein clustering K-means.py"
4. When prompted:
  - Enter a protein name e.g. Hemoglobin
  - If searching by name, specify how many sequences to retrieve (default is 50)

## This script will:
1. Fetch sequences from NCBI
2. Process and cluster them
3. Save results to protein_clusters.csv
4. Show visualizations and print evaluation metrics
   
   ### Evaluation Metrics
     After clustering, the script prints:
      - Best number of clusters (k)
      - Silhouette Scores for tested values of k
      - Davies-Bouldin Scores for tested values of k
   ### Visualization
     Two plots are shown if possible:
      - PCA-based clustering
      - UMAP-based clustering
      These help visualize how well-separated the clusters are.

## Configuration Notes
You can customize:

- Email for NCBI API : Update "your_email@gmail.com" at the top.
- K-mer size : Change k=3 in seq_to_kmer_vector().
- Clustering range : Modify the loop that tries different k values (currently 2–5).
- Max sequences fetched : Input during runtime or change default.

## Contributing
Contributions are welcome! Feel free to open issues or pull requests for:
- Bug fixes
- Improvements to visualization
- Exporting FASTA files per cluster
- Command-line interface support
- GUI version

📬 Contact
For questions or suggestions, reach out at:
📧 kristiadimikael@gmail.com
