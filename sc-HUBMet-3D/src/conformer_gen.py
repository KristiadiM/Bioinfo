"""
RDKit 3D Conformer Generation & Force-Field Optimization for Circulating Blood Metabolites.
"""

from pathlib import Path
from typing import Dict, Any, Optional
import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
except ImportError:
    Chem = None


def generate_metabolite_3d_sdf(
    metabolite_name: str,
    smiles: str,
    output_sdf_path: Path,
    num_confs: int = 25,
    max_iters: int = 500,
    random_seed: int = 42,
) -> Dict[str, Any]:
    """
    Generate an ensemble of 3D conformers from an isomeric/canonical SMILES string
    using RDKit ETKDGv3 distance geometry, optimize with MMFF94s (or UFF),
    and export the lowest-energy 3D conformer to an SDF file for GNINA docking.
    """
    if Chem is None:
        raise ImportError("RDKit is required for 3D conformer generation. Install via `pip install rdkit`.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string for {metabolite_name}: {smiles}")

    mol.SetProp("_Name", metabolite_name)
    mol_h = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.pruneRmsThresh = 0.5

    conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params))
    if not conf_ids:
        params.useRandomCoords = True
        conf_ids = list(AllChem.EmbedMultipleConfs(mol_h, numConfs=max(5, num_confs // 2), params=params))

    if not conf_ids:
        raise RuntimeError(f"Failed to embed 3D coordinates for {metabolite_name} ({smiles})")

    energies = []
    force_field_used = "MMFF94s"
    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant="MMFF94s")
        for cid in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(mol_h, mmff_props, confId=cid)
            if ff is not None:
                ff.Minimize(maxIts=max_iters)
                energies.append((cid, float(ff.CalcEnergy())))
    if not energies:
        force_field_used = "UFF"
        for cid in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(mol_h, confId=cid)
            if ff is not None:
                ff.Minimize(maxIts=max_iters)
                energies.append((cid, float(ff.CalcEnergy())))

    energies.sort(key=lambda x: x[1])
    best_cid, best_energy = energies[0]

    mol_h.SetProp("Metabolite_Name", metabolite_name)
    mol_h.SetProp("SMILES", smiles)
    mol_h.SetProp("ForceField", force_field_used)
    mol_h.SetProp("MinEnergy_kcal_mol", f"{best_energy:.3f}")

    mol_h.GetConformer(best_cid).Set3D(True)
    output_sdf_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf_path))
    writer.write(mol_h, confId=best_cid)
    writer.close()

    return {
        "metabolite_name": metabolite_name,
        "smiles": smiles,
        "sdf_path": str(output_sdf_path),
        "num_conformers_sampled": len(conf_ids),
        "force_field": force_field_used,
        "min_energy_kcal_mol": round(best_energy, 3),
        "mw": round(Descriptors.MolWt(mol), 2),
        "clogp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "rotatable_bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "hbd": int(rdMolDescriptors.CalcNumHBD(mol)),
        "hba": int(rdMolDescriptors.CalcNumHBA(mol)),
    }


def batch_generate_conformers(benchmark_df: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """
    Generate 3D SDF conformers for all unique metabolites in the benchmark dataframe.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    records = []
    unique_mets = benchmark_df[["metabolite_name", "pubchem_cid", "smiles", "metabolite_class"]].drop_duplicates(
        subset=["metabolite_name"]
    )
    for _, row in unique_mets.iterrows():
        safe_name = (
            row["metabolite_name"]
            .replace(" ", "_")
            .replace(":", "_")
            .replace("/", "_")
            .replace("-", "_")
        )
        sdf_path = out_dir / f"{safe_name}_CID{row['pubchem_cid']}.sdf"
        info = generate_metabolite_3d_sdf(
            metabolite_name=row["metabolite_name"],
            smiles=row["smiles"],
            output_sdf_path=sdf_path,
        )
        info["pubchem_cid"] = row["pubchem_cid"]
        info["metabolite_class"] = row["metabolite_class"]
        records.append(info)
    return pd.DataFrame(records)

