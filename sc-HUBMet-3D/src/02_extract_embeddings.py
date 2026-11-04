"""
Extract chemical embeddings (DeepChem/ChemBERTa-77M-MTR) and protein sequence
embeddings (facebook/esm2_t12_35M_UR50D) for the benchmark dataset.
"""

from pathlib import Path
import os
import logging
import pandas as pd
import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, AutoModel

os.environ["HF_HUB_DISABLE_IMPLICIT_TOKEN"] = "1"
os.environ["HF_HUB_DISABLE_SYMLINKS_WARNING"] = "1"
logging.getLogger("transformers").setLevel(logging.ERROR)
logging.getLogger("huggingface_hub").setLevel(logging.ERROR)


def main() -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Running embedding extraction on {device}")

    triplets_path = Path("data/processed/multiomics_triplets.csv")
    df = pd.read_csv(triplets_path)

    unique_mets = df[["metabolite", "smiles"]].drop_duplicates("metabolite")
    unique_prots = df[["uniprot_id", "sequence"]].drop_duplicates("uniprot_id")

    chem_model_name = "DeepChem/ChemBERTa-77M-MTR"
    chem_tok = AutoTokenizer.from_pretrained(chem_model_name)
    chem_mod = AutoModel.from_pretrained(chem_model_name, add_pooling_layer=False).to(device).eval()

    met_dict = {}
    with torch.no_grad():
        for _, row in unique_mets.iterrows():
            inputs = chem_tok(
                str(row["smiles"]),
                return_tensors="pt",
                max_length=128,
                truncation=True,
                padding=True,
            ).to(device)
            emb = chem_mod(**inputs).last_hidden_state.mean(dim=1).reshape(-1).cpu().to(torch.float32)
            met_dict[row["metabolite"]] = F.normalize(emb[:384].clone(), dim=0)

    print(f"Extracted {len(met_dict)} ChemBERTa (384-d) metabolite embeddings.")

    prot_model_name = "facebook/esm2_t12_35M_UR50D"
    prot_tok = AutoTokenizer.from_pretrained(prot_model_name)
    prot_mod = AutoModel.from_pretrained(prot_model_name, add_pooling_layer=False).to(device).eval()

    prot_dict = {}
    with torch.no_grad():
        for _, row in unique_prots.iterrows():
            seq = str(row["sequence"])[:1022]
            inputs = prot_tok(
                seq,
                return_tensors="pt",
                max_length=1024,
                truncation=True,
            ).to(device)
            emb = prot_mod(**inputs).last_hidden_state.mean(dim=1).reshape(-1).cpu().to(torch.float32)
            prot_dict[row["uniprot_id"]] = F.normalize(emb[:480].clone(), dim=0)

    print(f"Extracted {len(prot_dict)} ESM-2 (480-d) protein embeddings.")

    out_path = Path("data/processed/foundation_embeddings.pt")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"metabolites": met_dict, "proteins": prot_dict}, out_path)
    print(f"Saved embeddings to {out_path}")


if __name__ == "__main__":
    main()
