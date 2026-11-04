"""
Plot 5-fold cross-validation ROC/PR curves and cell-type metabolite signaling heatmaps.
"""

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score


def main() -> None:
    df = pd.read_csv("results/ablation_predictions.csv")
    fig1_path = Path("results/fig1_ablation_roc_pr.png")
    fig2_path = Path("results/fig2_single_cell_signaling_heatmap.png")
    fig1_path.parent.mkdir(parents=True, exist_ok=True)

    y_true = df["active_signaling"].values.astype(int)
    p_fm = df["pred_fm_only"].values
    p_3d = df["pred_structmet_3d"].values

    fpr_fm, tpr_fm, _ = roc_curve(y_true, p_fm)
    fpr_3d, tpr_3d, _ = roc_curve(y_true, p_3d)
    auc_fm = roc_auc_score(y_true, p_fm)
    auc_3d = roc_auc_score(y_true, p_3d)

    prec_fm, rec_fm, _ = precision_recall_curve(y_true, p_fm)
    prec_3d, rec_3d, _ = precision_recall_curve(y_true, p_3d)
    ap_fm = average_precision_score(y_true, p_fm)
    ap_3d = average_precision_score(y_true, p_3d)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), dpi=300)

    axes[0].plot(fpr_fm, tpr_fm, "--", color="#6c757d", lw=2.2, label=f"ChemBERTa + ESM-2 + Expr (AUC = {auc_fm:.3f})")
    axes[0].plot(fpr_3d, tpr_3d, "-", color="#d90429", lw=2.6, label=f"StructMet-scFM 3D Gate (AUC = {auc_3d:.3f})")
    axes[0].plot([0, 1], [0, 1], ":", color="#adb5bd", lw=1.2)
    axes[0].set_xlabel("False Positive Rate")
    axes[0].set_ylabel("True Positive Rate")
    axes[0].set_title("A. 5-Fold CV ROC Curve")
    axes[0].legend(loc="lower right", frameon=True)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(rec_fm, prec_fm, "--", color="#6c757d", lw=2.2, label=f"ChemBERTa + ESM-2 + Expr (AP = {ap_fm:.3f})")
    axes[1].plot(rec_3d, prec_3d, "-", color="#d90429", lw=2.6, label=f"StructMet-scFM 3D Gate (AP = {ap_3d:.3f})")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].set_title("B. 5-Fold CV Precision-Recall Curve")
    axes[1].legend(loc="lower left", frameon=True)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(fig1_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    pos_df = df[df["label"] == 1].copy()
    pos_df["pair_name"] = pos_df["metabolite"] + " - " + pos_df["gene_name"]
    piv = pos_df.pivot_table(index="pair_name", columns="cell_type", values="pred_structmet_3d", aggfunc="mean")

    fig2, ax2 = plt.subplots(figsize=(11.0, 8.0), dpi=300)
    sns.heatmap(
        piv,
        cmap="viridis",
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={"label": "Predicted Active Signaling Score"},
        ax=ax2,
    )
    ax2.set_title("Predicted Metabolite-Receptor Signaling Across Cell Types")
    ax2.set_xlabel("Cell Type")
    ax2.set_ylabel("Metabolite - Receptor Pair")
    plt.setp(ax2.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    plt.tight_layout()
    fig2.savefig(fig2_path, dpi=300, bbox_inches="tight")
    plt.close(fig2)

    print(f"Saved figures to {fig1_path} and {fig2_path}")


if __name__ == "__main__":
    main()
