"""
Train the 3D-Gated Cross-Attention Network (BiophysicalSignalingNet) with 5-fold
pair-grouped cross-validation comparing sequence/expression baseline vs. 3D-gated model.
"""

from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")


class BiophysicalSignalingNet(nn.Module):
    """
    3D-Gated Multi-Head Cross-Attention Network integrating ChemBERTa (384-d),
    ESM-2 (480-d), single-cell expression state (33-d), and 3D biophysical priors (4-d).
    """

    def __init__(self, met_dim: int = 384, prot_dim: int = 480, cell_dim: int = 33, d_3d: int = 4):
        super().__init__()
        self.proj_met = nn.Sequential(nn.Dropout(0.30), nn.Linear(met_dim, 64), nn.LayerNorm(64), nn.GELU())
        self.proj_prot = nn.Sequential(nn.Dropout(0.30), nn.Linear(prot_dim, 64), nn.LayerNorm(64), nn.GELU())
        self.proj_cell = nn.Sequential(nn.Linear(cell_dim, 64), nn.LayerNorm(64), nn.GELU())

        self.gate_3d = nn.Sequential(
            nn.Linear(d_3d, 64),
            nn.ReLU(),
            nn.Linear(64, 64),
            nn.Sigmoid(),
        )

        self.cross_attn = nn.MultiheadAttention(embed_dim=64, num_heads=4, batch_first=True)
        self.classifier = nn.Sequential(
            nn.Linear(128, 32),
            nn.GELU(),
            nn.Dropout(0.20),
            nn.Linear(32, 1),
        )

    def forward(self, z_met, z_prot, z_cell, x_3d, use_3d_prior: bool = True):
        h_met = self.proj_met(z_met)
        h_prot = self.proj_prot(z_prot)
        h_cell = self.proj_cell(z_cell)

        mol_pair = h_met * h_prot
        if use_3d_prior:
            gate = self.gate_3d(x_3d)
            mol_pair = (mol_pair + 1.0) * gate

        query = h_cell.unsqueeze(1)
        kv = torch.stack([mol_pair, h_prot], dim=1)
        attn_out, _ = self.cross_attn(query, kv, kv)

        combined = torch.cat([attn_out.squeeze(1), mol_pair], dim=-1)
        return self.classifier(combined).view(-1)


def build_feature_tensors(df: pd.DataFrame, fm_dict: dict, cell_dict: dict, prior_dict: dict):
    met_emb = fm_dict["metabolites"]
    prot_emb = fm_dict["proteins"]

    zm_list, zp_list, zc_list, x3d_list = [], [], [], []
    for _, row in df.iterrows():
        zm_list.append(torch.as_tensor(met_emb[row["metabolite"]], dtype=torch.float32))
        zp_list.append(torch.as_tensor(prot_emb[row["uniprot_id"]], dtype=torch.float32))
        c_pca = torch.as_tensor(cell_dict[row["cell_type"]], dtype=torch.float32)
        expr = torch.tensor([float(row["receptor_log1p_expr"])], dtype=torch.float32)
        zc_list.append(torch.cat([c_pca, expr], dim=0))
        x3d_list.append(torch.as_tensor(prior_dict[(row["metabolite"], row["uniprot_id"])], dtype=torch.float32).numpy())

    zm_all = torch.stack(zm_list)
    zp_all = torch.stack(zp_list)
    zc_all = torch.stack(zc_list)
    x3d_all = np.vstack(x3d_list).astype(np.float32)
    y_all = torch.tensor(df["active_signaling"].values, dtype=torch.float32)
    return zm_all, zp_all, zc_all, x3d_all, y_all


def run_5fold_cv(df: pd.DataFrame, fm_dict: dict, cell_dict: dict, prior_dict: dict, use_3d_prior: bool) -> np.ndarray:
    torch.manual_seed(42)
    np.random.seed(42)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    zm_all, zp_all, zc_all, x3d_all, y_all = build_feature_tensors(df, fm_dict, cell_dict, prior_dict)
    groups = df["pair_id"].values if "pair_id" in df.columns else (df["metabolite"] + "_" + df["gene_name"]).values

    sgkf = StratifiedGroupKFold(n_splits=5)
    oof_preds = np.zeros(len(df), dtype=np.float32)

    for fold_idx, (tr_idx, va_idx) in enumerate(sgkf.split(df, df["active_signaling"].values, groups=groups)):
        torch.manual_seed(42 + fold_idx)
        scaler = StandardScaler()
        x3d_tr = torch.tensor(scaler.fit_transform(x3d_all[tr_idx]), dtype=torch.float32)
        x3d_va = torch.tensor(scaler.transform(x3d_all[va_idx]), dtype=torch.float32)

        tr_ds = TensorDataset(zm_all[tr_idx], zp_all[tr_idx], zc_all[tr_idx], x3d_tr, y_all[tr_idx])
        va_ds = TensorDataset(zm_all[va_idx], zp_all[va_idx], zc_all[va_idx], x3d_va, y_all[va_idx])

        tr_loader = DataLoader(tr_ds, batch_size=24, shuffle=True)
        va_loader = DataLoader(va_ds, batch_size=64, shuffle=False)

        model = BiophysicalSignalingNet().to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=2e-2)
        criterion = nn.BCEWithLogitsLoss()

        model.train()
        for _ in range(35):
            for zm, zp, zc, x3d, yb in tr_loader:
                zm, zp, zc, x3d, yb = zm.to(device), zp.to(device), zc.to(device), x3d.to(device), yb.to(device)
                optimizer.zero_grad()
                logits = model(zm, zp, zc, x3d, use_3d_prior=use_3d_prior)
                loss = criterion(logits, yb)
                loss.backward()
                optimizer.step()

        model.eval()
        fold_probs = []
        with torch.no_grad():
            for zm, zp, zc, x3d, _ in va_loader:
                zm, zp, zc, x3d = zm.to(device), zp.to(device), zc.to(device), x3d.to(device)
                probs = torch.sigmoid(model(zm, zp, zc, x3d, use_3d_prior=use_3d_prior)).cpu().numpy()
                fold_probs.extend(probs.tolist())

        oof_preds[va_idx] = np.array(fold_probs, dtype=np.float32)

    return oof_preds


def main() -> None:
    df = pd.read_csv("data/processed/multiomics_triplets.csv")
    fm_dict = torch.load("data/processed/foundation_embeddings.pt", map_location="cpu", weights_only=False)
    cell_dict = torch.load("data/processed/cell_embeddings.pt", map_location="cpu", weights_only=False)
    prior_dict = torch.load("data/processed/structural_priors_3d.pt", map_location="cpu", weights_only=False)

    preds_fm_only = run_5fold_cv(df, fm_dict, cell_dict, prior_dict, use_3d_prior=False)
    preds_structmet_3d = run_5fold_cv(df, fm_dict, cell_dict, prior_dict, use_3d_prior=True)

    y_true = df["active_signaling"].values.astype(int)
    auc_fm = roc_auc_score(y_true, preds_fm_only)
    ap_fm = average_precision_score(y_true, preds_fm_only)
    auc_3d = roc_auc_score(y_true, preds_structmet_3d)
    ap_3d = average_precision_score(y_true, preds_structmet_3d)

    print(f"Baseline (use_3d_prior=False): ROC-AUC = {auc_fm:.4f} | PR-AUC = {ap_fm:.4f}")
    print(f"Full Model (use_3d_prior=True): ROC-AUC = {auc_3d:.4f} | PR-AUC = {ap_3d:.4f}")

    out_df = df.copy()
    out_df["pred_fm_only"] = np.round(preds_fm_only, 4)
    out_df["pred_structmet_3d"] = np.round(preds_structmet_3d, 4)

    Path("results").mkdir(parents=True, exist_ok=True)
    out_csv = Path("results/ablation_predictions.csv")
    out_df.to_csv(out_csv, index=False)
    print(f"Saved predictions to {out_csv}")


if __name__ == "__main__":
    main()
