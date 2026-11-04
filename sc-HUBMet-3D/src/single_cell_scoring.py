"""
Single-cell expression loading from data/iav_36_pairs_benchmark.xlsx and
3D-prior-weighted metabolite sensing score calculation.
"""

from pathlib import Path
from typing import List, Tuple
import numpy as np
import pandas as pd


EXCEL_PATH = Path("data/iav_36_pairs_benchmark.xlsx")


def load_cell_expression_table(
    target_genes: List[str] | None = None,
    excel_path: Path = EXCEL_PATH,
) -> pd.DataFrame:
    """
    Load the cell_expression_ntpm sheet from data/iav_36_pairs_benchmark.xlsx
    and return a gene x cell_type matrix of nTPM values.
    """
    if not excel_path.exists():
        raise FileNotFoundError(f"Workbook not found: {excel_path}")

    cell_by_gene = pd.read_excel(excel_path, sheet_name="cell_expression_ntpm", index_col=0).fillna(0.0)
    gene_by_cell = cell_by_gene.T
    gene_by_cell.index.name = "gene_symbol"
    if target_genes is not None:
        present = [g for g in target_genes if g in gene_by_cell.index]
        gene_by_cell = gene_by_cell.loc[present]
    return gene_by_cell


def compute_celltype_metabolite_sensing_matrices(
    docked_df: pd.DataFrame,
    expr_gene_by_cell: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute unweighted statistical sensing and 3D-prior-weighted sensing matrices
    across cell types and metabolites.
    """
    log_expr = np.log2(1.0 + expr_gene_by_cell)
    cell_types = list(log_expr.columns)
    met_col = "metabolite" if "metabolite" in docked_df.columns else "metabolite_name"
    gene_col = "gene_name" if "gene_name" in docked_df.columns else "gene_symbol"
    r_col = "hubmet_iav_r" if "hubmet_iav_r" in docked_df.columns else "hubmet_stat_r"

    metabolites = sorted(docked_df[met_col].unique())
    unweighted_mat = pd.DataFrame(0.0, index=cell_types, columns=metabolites)
    weighted_3d_mat = pd.DataFrame(0.0, index=cell_types, columns=metabolites)
    long_records = []

    for _, row in docked_df.iterrows():
        met = row[met_col]
        gene = row[gene_col]
        stat_r = abs(float(row[r_col]))
        p_3d = float(row["prior_3d_weight"])

        if gene not in log_expr.index:
            continue

        for ct in cell_types:
            expr_val = float(log_expr.loc[gene, ct])
            unw = expr_val * stat_r
            w3d = expr_val * stat_r * p_3d
            unweighted_mat.loc[ct, met] += unw
            weighted_3d_mat.loc[ct, met] += w3d
            long_records.append({
                "cell_type": ct,
                "metabolite": met,
                "metabolite_class": row["metabolite_class"],
                "gene_name": gene,
                "protein_family": row["protein_family"],
                "tier": row["tier"],
                "log2_expr": round(expr_val, 3),
                "abs_hubmet_r": round(stat_r, 3),
                "prior_3d_weight": round(p_3d, 4),
                "unweighted_stat_sense": round(unw, 4),
                "weighted_3d_sense": round(w3d, 4),
            })

    return unweighted_mat, weighted_3d_mat, pd.DataFrame(long_records)
