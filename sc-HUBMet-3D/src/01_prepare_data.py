"""
Prepare benchmark pairs, protein sequences, 3D receptor structures, and single-cell
expression matrices from data/iav_36_pairs_benchmark.xlsx.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import requests
import torch
from sklearn.decomposition import PCA
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1


EXCEL_PATH = Path("data/iav_36_pairs_benchmark.xlsx")


def fetch_uniprot_fasta(uid: str, pdb_path: Path | None = None, timeout: int = 15) -> str:
    """Fetch canonical protein sequence from UniProt REST API, or extract from PDB if present."""
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    resp = requests.get(url, timeout=timeout)
    if resp.status_code == 200 and resp.text.startswith(">"):
        lines = [ln.strip() for ln in resp.text.splitlines() if not ln.startswith(">")]
        seq = "".join(lines)
        if len(seq) > 20:
            return seq

    if pdb_path is not None and pdb_path.exists():
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("rec", str(pdb_path))
        residues = []
        for res in struct[0].get_residues():
            if res.id[0] == " ":
                rname = res.get_resname().strip().capitalize()
                residues.append(protein_letters_3to1.get(rname, "X"))
        if residues:
            return "".join(residues)

    raise RuntimeError(f"Failed to retrieve sequence for UniProt ID {uid}")


def fetch_structure_pdb(structure_source: str, pdb_id: str, uid: str, out_dir: Path, timeout: int = 30) -> Path:
    """Download experimental PDB structure from RCSB or predicted structure from AlphaFold DB."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if str(structure_source).strip().upper() == "PDB":
        clean_id = pdb_id.strip().upper()
        out_pdb = out_dir / f"{clean_id}.pdb"
        if out_pdb.exists() and out_pdb.stat().st_size > 1000:
            return out_pdb
        resp = requests.get(f"https://files.rcsb.org/download/{clean_id}.pdb", timeout=timeout)
        resp.raise_for_status()
        out_pdb.write_text(resp.text, encoding="utf-8")
        return out_pdb

    out_pdb = out_dir / f"AF-{uid}-F1.pdb"
    if out_pdb.exists() and out_pdb.stat().st_size > 1000:
        return out_pdb

    for ver in ("v4", "v3", "v2"):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{ver}.pdb"
        resp = requests.get(url, timeout=timeout)
        if resp.status_code == 200 and "ATOM" in resp.text:
            out_pdb.write_text(resp.text, encoding="utf-8")
            return out_pdb

    api_resp = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}", timeout=timeout)
    api_resp.raise_for_status()
    records = api_resp.json()
    if records and "pdbUrl" in records[0]:
        pdb_resp = requests.get(records[0]["pdbUrl"], timeout=timeout)
        pdb_resp.raise_for_status()
        out_pdb.write_text(pdb_resp.text, encoding="utf-8")
        return out_pdb

    raise RuntimeError(f"Unable to fetch structure for {uid} ({pdb_id})")


def main() -> None:
    if not EXCEL_PATH.exists():
        raise FileNotFoundError(f"Missing input workbook: {EXCEL_PATH}")

    for folder in ["data/structures/receptors_raw", "data/structures/ligands", "data/processed", "results"]:
        Path(folder).mkdir(parents=True, exist_ok=True)

    pairs_df = pd.read_excel(EXCEL_PATH, sheet_name="benchmark_pairs")
    expr_df = pd.read_excel(EXCEL_PATH, sheet_name="cell_expression_ntpm", index_col=0).fillna(0.0)

    uid_to_seq = {}
    raw_dir = Path("data/structures/receptors_raw")
    for _, row in pairs_df.drop_duplicates(subset=["uniprot_id", "pdb_id"]).iterrows():
        uid = str(row["uniprot_id"]).strip()
        pdb_path = fetch_structure_pdb(row["structure_source"], str(row["pdb_id"]), uid, raw_dir)
        if uid not in uid_to_seq:
            uid_to_seq[uid] = fetch_uniprot_fasta(uid, pdb_path=pdb_path)

    pairs_df["sequence"] = pairs_df["uniprot_id"].map(uid_to_seq)
    print(f"Loaded {len(pairs_df)} benchmark pairs from {EXCEL_PATH}")

    X = expr_df.values.astype(np.float32)
    row_sums = np.maximum(X.sum(axis=1, keepdims=True), 1e-6)
    X_log1p = np.log1p((X / row_sums) * 1e4)
    log1p_df = pd.DataFrame(X_log1p, index=expr_df.index.astype(str), columns=expr_df.columns.astype(str))

    n_comps = min(32, X_log1p.shape[0], X_log1p.shape[1])
    pca_mat = PCA(n_components=n_comps, random_state=42).fit_transform(X_log1p)
    if pca_mat.shape[1] < 32:
        pca_mat = np.pad(pca_mat, ((0, 0), (0, 32 - pca_mat.shape[1])), mode="constant").astype(np.float32)
    else:
        pca_mat = pca_mat[:, :32].astype(np.float32)

    cell_emb_dict = {ct: torch.tensor(pca_mat[i], dtype=torch.float32) for i, ct in enumerate(log1p_df.index)}
    torch.save(cell_emb_dict, "data/processed/cell_embeddings.pt")

    triplet_rows = []
    for ct in log1p_df.index:
        for _, prow in pairs_df.iterrows():
            gene = prow["gene_name"]
            expr_val = float(log1p_df.loc[ct, gene]) if gene in log1p_df.columns else 0.0
            active_sig = 1 if (int(prow["label"]) == 1 and expr_val > 0.5) else 0
            triplet_rows.append({
                "cell_type": ct,
                "pair_id": prow["pair_id"],
                "metabolite": prow["metabolite"],
                "smiles": prow["smiles"],
                "uniprot_id": prow["uniprot_id"],
                "gene_name": gene,
                "sequence": prow["sequence"],
                "tier": prow["tier"],
                "label": int(prow["label"]),
                "receptor_log1p_expr": round(expr_val, 4),
                "active_signaling": int(active_sig),
                "iav_pathway": prow["iav_pathway"],
                "hubmet_iav_r": float(prow["hubmet_iav_r"]),
            })

    triplets_df = pd.DataFrame(triplet_rows)
    out_csv = Path("data/processed/multiomics_triplets.csv")
    triplets_df.to_csv(out_csv, index=False)
    print(f"Saved {len(triplets_df)} multi-omics triplets to {out_csv}")


if __name__ == "__main__":
    main()
