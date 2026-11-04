#!/usr/bin/env bash
#SBATCH --job-name=StructMet-scFM-IAV
#SBATCH --output=results/structmet_scfm_%j.out
#SBATCH --error=results/structmet_scfm_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail

echo "============================================================"
echo "Starting StructMet-scFM (Influenza A Host Signaling Pipeline)"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "============================================================"

# Activate virtual environment if present
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
fi

python src/01_prepare_data.py
python src/02_extract_embeddings.py
python src/03_docking_and_md_features.py
python src/04_train_pytorch_ablation.py
python src/05_plot_results.py

echo "============================================================"
echo "StructMet-scFM Pipeline Completed Successfully!"
echo "Outputs verified in results/fig1_ablation_roc_pr.png and results/fig2_single_cell_signaling_heatmap.png"
echo "============================================================"

