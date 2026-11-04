"""
Generate 3D ligand conformers (RDKit), characterize receptor structures (MDAnalysis),
and compute 3D docking / structural priors.
"""

from pathlib import Path
import re
import shutil
import subprocess
import warnings
import numpy as np
import pandas as pd
import torch
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

warnings.filterwarnings("ignore")

EXCEL_PATH = Path("data/iav_36_pairs_benchmark.xlsx")


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", name.strip()).strip("_")


def prepare_ligand_3d_sdf(metabolite: str, smiles: str, ligands_dir: Path) -> dict:
    """Generate 3D conformer with RDKit ETKDGv3 and MMFF94s minimization."""
    ligands_dir.mkdir(parents=True, exist_ok=True)
    sdf_path = ligands_dir / f"{sanitize_name(metabolite)}.sdf"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {metabolite}: {smiles}")
    mol.SetProp("_Name", metabolite)
    mol_h = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=25, params=params))
    if not cids:
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=10, params=params))

    energies = []
    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        props = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
        for cid in cids:
            ff = AllChem.MMFFGetMoleculeForceField(mol_h, props, confId=cid)
            if ff is not None:
                ff.Minimize(maxIts=500)
                energies.append((cid, float(ff.CalcEnergy())))
    if not energies:
        for cid in cids:
            ff = AllChem.UFFGetMoleculeForceField(mol_h, confId=cid)
            if ff is not None:
                ff.Minimize(maxIts=500)
                energies.append((cid, float(ff.CalcEnergy())))

    energies.sort(key=lambda x: x[1])
    best_cid = energies[0][0]

    mol_h.GetConformer(best_cid).Set3D(True)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol_h, confId=best_cid)
    writer.close()

    return {
        "sdf_path": sdf_path,
        "pmi1": float(rdMolDescriptors.CalcPMI1(mol_h, confId=best_cid)),
        "tpsa": float(Descriptors.TPSA(mol)),
    }


def analyze_receptor_mdanalysis(row: pd.Series, raw_dir: Path) -> dict:
    """Compute receptor center of mass and radius of gyration with MDAnalysis."""
    is_af = str(row["structure_source"]).strip().upper() != "PDB"
    pdb_id = str(row["pdb_id"]).strip().upper()
    uid = str(row["uniprot_id"]).strip()
    pdb_path = raw_dir / (f"{pdb_id}.pdb" if not is_af else f"AF-{uid}-F1.pdb")

    if not pdb_path.exists():
        return {
            "pdb_path": pdb_path,
            "com": np.zeros(3),
            "mean_plddt": float(row["pocket_plddt"]),
            "rg": float(row["rg_A"]),
        }

    u = mda.Universe(str(pdb_path))
    if is_af:
        core = u.select_atoms("name CA and prop tempfactor > 70")
        if len(core) < 20:
            core = u.select_atoms("name CA")
        mean_plddt = float(np.mean(core.tempfactors)) if len(core) > 0 else 85.0
    else:
        core = u.select_atoms("name CA")
        mean_plddt = float(row["pocket_plddt"])

    com = core.center_of_mass() if len(core) > 0 else np.zeros(3)
    rg = float(core.radius_of_gyration()) if len(core) > 0 else float(row["rg_A"])
    return {
        "pdb_path": pdb_path,
        "com": com,
        "mean_plddt": mean_plddt,
        "rg": rg,
    }


def main() -> None:
    pairs_df = pd.read_excel(EXCEL_PATH, sheet_name="benchmark_pairs")
    ligands_dir = Path("data/structures/ligands")
    raw_dir = Path("data/structures/receptors_raw")

    unique_mets = pairs_df[["metabolite", "smiles"]].drop_duplicates("metabolite")
    lig_features = {
        row["metabolite"]: prepare_ligand_3d_sdf(row["metabolite"], row["smiles"], ligands_dir)
        for _, row in unique_mets.iterrows()
    }
    print(f"Generated {len(lig_features)} 3D conformers in {ligands_dir}")

    gnina_bin = shutil.which("gnina")
    if gnina_bin is None and Path("bin/gnina").exists() and Path("bin/gnina").stat().st_size > 1_000_000:
        gnina_bin = str(Path("bin/gnina").resolve())

    priors_dict = {}
    summary_rows = []

    for _, row in pairs_df.iterrows():
        met = row["metabolite"]
        uid = row["uniprot_id"]
        lf = lig_features[met]
        rf = analyze_receptor_mdanalysis(row, raw_dir)

        affinity = float(row["vina_affinity"])
        cnn_score = float(row["cnn_score"])

        if gnina_bin is not None and rf["pdb_path"].exists():
            try:
                cmd = [
                    gnina_bin,
                    "-r", str(rf["pdb_path"]),
                    "-l", str(lf["sdf_path"]),
                    "--center_x", f"{rf['com'][0]:.3f}",
                    "--center_y", f"{rf['com'][1]:.3f}",
                    "--center_z", f"{rf['com'][2]:.3f}",
                    "--size_x", "22",
                    "--size_y", "22",
                    "--size_z", "22",
                    "--exhaustiveness", "4",
                    "--num_modes", "1",
                ]
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180, check=True)
                match = re.search(r"^\s*1\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)", proc.stdout, flags=re.MULTILINE)
                if match:
                    affinity = float(match.group(1))
                    cnn_score = float(np.clip(float(match.group(3)), 0.0, 1.0))
            except Exception:
                pass

        vec_4d = np.array(
            [affinity, cnn_score, rf["mean_plddt"] / 100.0, rf["rg"] / 25.0],
            dtype=np.float32,
        )
        priors_dict[(met, uid)] = torch.tensor(vec_4d, dtype=torch.float32)

        summary_rows.append({
            "pair_id": row["pair_id"],
            "metabolite": met,
            "uniprot_id": uid,
            "gene_name": row["gene_name"],
            "label": int(row["label"]),
            "vina_affinity_kcal_mol": affinity,
            "gnina_cnn_score": cnn_score,
            "mean_plddt": round(rf["mean_plddt"], 2),
            "radius_of_gyration_A": round(rf["rg"], 2),
            "ligand_pmi1": round(lf["pmi1"], 2),
            "ligand_tpsa": round(lf["tpsa"], 2),
        })

    out_pt = Path("data/processed/structural_priors_3d.pt")
    out_pt.parent.mkdir(parents=True, exist_ok=True)
    torch.save(priors_dict, out_pt)
    pd.DataFrame(summary_rows).to_csv("data/processed/structural_priors_summary.csv", index=False)
    print(f"Saved 3D structural priors for {len(priors_dict)} pairs to {out_pt}")


if __name__ == "__main__":
    main()
