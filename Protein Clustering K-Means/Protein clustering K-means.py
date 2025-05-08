"""Protein Clustering with K-Means"""

from Bio import Entrez, SeqIO
from collections import Counter
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, davies_bouldin_score
import matplotlib.pyplot as plt
import seaborn as sns
import umap.umap_ as umap

# Set your email
Entrez.email = "your_email@gmail.com"

def fetch_by_name(protein_name):
    try:
        handle = Entrez.efetch(db="protein", protein=protein_name, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records
    except Exception as e:
        print(f"Error fetching accession {protein_name}: {e}")
        return []

def fetch_by_gene_name(protein_name, max_records):
    search_term = f"{protein_name}[Protein Name]"
    search_handle = Entrez.esearch(db="protein", term=search_term, retmax=max_records)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    ids = search_results["IdList"]

    if not ids:
        print("No records found.")
        return []
    fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    fetch_handle.close()
    return records

# Input Session
print("NCBI Protein Sequence Fetcher")
user_input = input("Enter Protein's Name: ").strip()
is_accession = user_input[:2].isalpha() and "_" in user_input

if is_accession:
    sequences = fetch_by_accession(user_input)
else:
    try:
        max_records = int(input("Enter max number of sequences (e.g., 50): "))
    except ValueError:
        max_records = 50
    sequences = fetch_by_gene_name(user_input, max_records)

# Data Preprocessing
print(f"\nInitial sequences fetched: {len(sequences)}")
valid_seqs = {}
for record in sequences:
    seq_str = str(record.seq).upper()
    if len(seq_str) >= 30 and set(seq_str).issubset(set("ACDEFGHIKLMNPQRSTVWY")):
        valid_seqs[seq_str] = record  # Remove duplicates

sequences = list(valid_seqs.values())
print(f"Filtered valid sequences (unique, valid AAs, ≥30 AA): {len(sequences)}")

if len(sequences) < 2:
    print("⚠️ Not enough sequences to proceed. Exiting.")
    exit()

# K-mer Feature Extraction
def protein_kmers(seq, k=3):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

def seq_to_kmer_vector(seq, k=3):
    return Counter(protein_kmers(str(seq), k=k))

all_kmers = set()
vectors = []

for record in sequences:
    vec = seq_to_kmer_vector(record.seq)
    vectors.append(vec)
    all_kmers.update(vec.keys())

df_kmers = pd.DataFrame(vectors, columns=sorted(all_kmers)).fillna(0)
print("Raw feature matrix shape:", df_kmers.shape)

# Remove low-variance and duplicate features
selector = VarianceThreshold(threshold=0.01)
X_reduced = selector.fit_transform(df_kmers)
selected_columns = df_kmers.columns[selector.get_support()]
df_kmers = pd.DataFrame(X_reduced, columns=selected_columns)
df_kmers = df_kmers.drop_duplicates().reset_index(drop=True)
print("Final feature matrix shape:", df_kmers.shape)

# Standardization
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_kmers)

# Clustering
n_samples = X_scaled.shape[0]
silhouette_scores = {}
db_scores = {}
best_score = -1
best_k = None

if n_samples >= 3:
    for k in range(2, min(6, n_samples)):
        kmeans = KMeans(n_clusters=k, random_state=42)
        clusters = kmeans.fit_predict(X_scaled)

        if len(set(clusters)) > 1:
            sil = silhouette_score(X_scaled, clusters)
            db = davies_bouldin_score(X_scaled, clusters)
            silhouette_scores[k] = sil
            db_scores[k] = db
            if sil > best_score:
                best_score = sil
                best_k = k
        else:
            print(f"⚠️ Only one cluster formed at k={k}, skipping.")

if best_k:
    print(f"\n✅ Best number of clusters: {best_k}")
    kmeans = KMeans(n_clusters=best_k, random_state=42)
    clusters = kmeans.fit_predict(X_scaled)
    df_kmers['Cluster'] = clusters
else:
    print("\n⚠️ Insufficient data for clustering. Marking as 'Unassigned'.")
    df_kmers['Cluster'] = 'Unassigned'

# Export results to CSV
print("\nExporting results to 'protein_clusters.csv'...")
results = []

for idx, record in enumerate(sequences):
    results.append({
        "Protein_ID": record.id,
        "Description": record.description,
        "Sequence": str(record.seq),
        "Cluster": df_kmers.loc[idx, 'Cluster']
    })

results_df = pd.DataFrame(results)
results_df.to_csv("protein_clusters.csv", index=False)
print("✅ Results exported successfully.")

# Visualization
if best_k and n_samples >= 3:
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X_scaled)

    # UMAP: Ensure n_neighbors < n_samples
    n_neighbors = min(10, max(2, n_samples - 1))  # Ensure at least 2 and < n_samples
    try:
        umap_model = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.1, random_state=42)
        umap_result = umap_model.fit_transform(X_scaled)
        use_umap = True
    except Exception as e:
        print(f"⚠️ UMAP failed: {e}")
        umap_result = None
        use_umap = False

    # Prepare for visualization
    viz_df = pd.DataFrame({
        'PCA1': pca_result[:, 0],
        'PCA2': pca_result[:, 1],
        'Cluster': df_kmers['Cluster'].astype(str)
    })
    if use_umap:
        viz_df['UMAP1'] = umap_result[:, 0]
        viz_df['UMAP2'] = umap_result[:, 1]

    # Plot
    fig, axes = plt.subplots(1, 2 if use_umap else 1, figsize=(16, 7))
    if not use_umap:
        axes = [axes]  # Make it a list to keep indexing consistent


    sns.scatterplot(ax=axes[0], x='PCA1', y='PCA2', hue='Cluster', data=viz_df, palette='Set2', s=70)
    axes[0].set_title(f"PCA Clustering (k={best_k})")

    if use_umap:
        sns.scatterplot(ax=axes[1], x='UMAP1', y='UMAP2', hue='Cluster', data=viz_df, palette='Set2', s=70)
        axes[1].set_title(f"UMAP Clustering (k={best_k})")

    plt.tight_layout()
    plt.show()


    print("\n📊 Silhouette Scores:")
    for k, s in silhouette_scores.items():
        print(f"  k={k}: {s:.4f}")

    print("\n📊 Davies-Bouldin Scores:")
    for k, d in db_scores.items():
        print(f"  k={k}: {d:.4f}")
else:
    print("\n⚠️ Skipped clustering visualization due to insufficient or unclusterable data.")