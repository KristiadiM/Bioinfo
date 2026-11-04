# HUBMet-3D: 3D Structure-Gated Metabolite-Receptor Signaling in Acute Influenza A Infection

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/your-username/HUBMet/blob/main/sc_hubmet_3d_pipeline.ipynb)

During acute Influenza A virus (IAV) infection, systemic hyperinflammation (the "cytokine storm") elevates circulating immunometabolites and inflammatory cytokines simultaneously. Purely statistical co-abundance analysis ($|r|$) frequently correlates metabolites with abundant acute-phase proteins and interferon-stimulated bystanders (`CXCL10`, `IL6`, `MX1`, `CRP`, `B2M`, `IFNB1`) rather than physical transmembrane sensors or transporters.

`HUBMet-3D` addresses this false-positive problem by combining 3D structural pocket characterization (`RDKit`, `Bio.PDB`, and `GNINA` molecular docking) with chemical and protein language models (`ChemBERTa` and `ESM-2`) and single-cell transcriptomics (`scanpy`).

---

## Key Results

1. **Filtering Cytokine Storm Decoys (Figure 1)**: Statistical co-abundance alone yields a deceptive ROC-AUC of **0.351** because bystanders are expressed at extreme levels during infection. Applying the 3D biophysical prior gate ($P_{3D}$) crushes decoy scores to $\approx 0.00$, boosting true receptor discrimination to **0.877 ROC-AUC**.
2. **Cross-Attention Foundation Model (Figure 2)**: 5-fold pair-grouped cross-validation (`StratifiedGroupKFold`) combining ChemBERTa (384-d), ESM-2 (480-d), and cell-type expression achieves **0.903 ROC-AUC** and **0.869 PR-AUC** with the 3D physical gate.
3. **Single-Cell Immune Sensing Maps (Figures 3 & 4)**: Projects 3D-gated metabolite signaling across 12 respiratory/immune cell types and 2,638 single PBMCs, eliminating false-positive cytokine bystander signals and isolating genuine host sensors (`S1PR1`, `PTGER4`, `AHR`, `SLC7A5`).

---

## Repository Structure

```text
├── data/
│   └── iav_36_pairs_benchmark.xlsx       # Input workbook (benchmark_pairs & cell_expression_ntpm)
├── src/
│   ├── 01_prepare_data.py                # Load Excel tables, fetch UniProt sequences & PDBs
│   ├── 02_extract_embeddings.py          # Extract ChemBERTa (384-d) and ESM-2 (480-d) embeddings
│   ├── 03_docking_and_md_features.py     # RDKit conformers, pocket Rg, and 3D structural priors
│   ├── 04_train_pytorch_ablation.py      # 5-fold pair-grouped CV with BiophysicalSignalingNet
│   ├── 05_plot_results.py                # Plot ROC/PR curves and cell-type signaling heatmaps
│   ├── conformer_gen.py                  # RDKit ETKDGv3 / MMFF94s conformer generation
│   ├── structure_pocket.py               # Sequence-matched chain selection & pocket detection
│   ├── gnina_docking.py                  # GNINA CLI wrapper & contact calculator (< 4.0 A)
│   └── single_cell_scoring.py            # Cell-type metabolite sensing matrix computation
├── slurm/
│   └── run_pipeline.sh                   # SLURM batch execution script
├── sc_hubmet_3d_pipeline.ipynb           # Primary Jupyter / Google Colab notebook (with rendered outputs)
├── requirements.txt                      # Python dependencies
└── README.md
```

---

## Input Data (`data/iav_36_pairs_benchmark.xlsx`)

All inputs are packaged into a single Excel workbook:
* **`benchmark_pairs`** (36 curated Influenza A metabolite-protein pairs):
  * **Tier 1 (Validated Binders, 12 pairs)**: Validated host immunometabolic receptors and transporters (`HCAR1`, `SUCNR1`, `AHR`, `S1PR1`, `PTGER4`, `ADORA2A`, `GPBAR1`, `SLC7A5`, `SLC1A5`, `FFAR2`, `FFAR1`, `HRH3`).
  * **Tier 2 (Putative Sensors, 12 pairs)**: Candidate metabolite sensors (`OXGR1`, `GPR132`, `GPR35`, `GPR183`, `SLC22A5`, `SLC16A1`, `SLC6A6`, `SLC29A1`, `GPR55`, `NR1I2`) and controls.
  * **Tier 3 (Cytokine Storm Decoys, 12 pairs)**: High-correlation bystander proteins (`CXCL10`, `IL6`, `IFNB1`, `ANGPT2`, `IL1B`, `MX1`, `CRP`, `TNF`, `S100A8`, `CD4`, `LEPR`, `B2M`).
* **`cell_expression_ntpm`**: Reference expression matrix (nTPM) across 12 respiratory, vascular, hepatic, and immune cell lineages.

---

## Getting Started

### 1. Google Colab (Recommended for GPU Docking)
Open [`sc_hubmet_3d_pipeline.ipynb`](sc_hubmet_3d_pipeline.ipynb) in Google Colab:
1. Set the runtime to GPU: **Runtime $\rightarrow$ Change runtime type $\rightarrow$ T4 GPU**.
2. Run Cell 0 to install dependencies and configure GNINA.
3. When Cell 1 runs, upload `data/iav_36_pairs_benchmark.xlsx` (or place it in `data/`).
4. Run all cells.

### 2. Local Python Environment
Clone the repository and install requirements:

```bash
git clone https://github.com/your-username/HUBMet.git
cd HUBMet

pip install -r requirements.txt
```

Run the pipeline sequentially via the command-line scripts:

```bash
python src/01_prepare_data.py
python src/02_extract_embeddings.py
python src/03_docking_and_md_features.py
python src/04_train_pytorch_ablation.py
python src/05_plot_results.py
```

---

## Dependencies
* Python $\ge$ 3.10
* PyTorch $\ge$ 2.0.0
* Transformers $\ge$ 4.35.0 (`DeepChem/ChemBERTa-77M-MTR`, `facebook/esm2_t12_35M_UR50D`)
* RDKit $\ge$ 2023.9.1
* Biopython $\ge$ 1.81
* MDAnalysis $\ge$ 2.6.0
* Scanpy $\ge$ 1.9.6
* GNINA $\ge$ v1.1

📬 Contact For questions or suggestions, reach out at: 📧 kristiadimikael@gmail.com
