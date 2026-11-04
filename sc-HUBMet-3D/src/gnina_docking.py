"""
GNINA 3D CNN Molecular Docking, Binding Pocket Contact Extraction,
and Composite 3D Structural Binding Prior Calculation.
"""

from pathlib import Path
from typing import Dict, Any, List, Optional
import subprocess
import numpy as np
import pandas as pd

try:
    from rdkit import Chem
    from Bio.PDB import PDBParser, NeighborSearch
except ImportError:
    Chem = None
    PDBParser = None


def run_gnina_docking_pair(
    gnina_bin: str,
    receptor_pdb: Path,
    ligand_sdf: Path,
    output_sdf: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 22.0,
    size_y: float = 22.0,
    size_z: float = 22.0,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    seed: int = 42,
    timeout: int = 300,
) -> Dict[str, Any]:
    """
    Run GNINA molecular docking for a single metabolite-protein pair and parse the
    top-ranked pose by CNNscore / CNNaffinity.
    """
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(gnina_bin),
        "--receptor", str(receptor_pdb),
        "--ligand", str(ligand_sdf),
        "--out", str(output_sdf),
        "--center_x", f"{center_x:.3f}",
        "--center_y", f"{center_y:.3f}",
        "--center_z", f"{center_z:.3f}",
        "--size_x", f"{size_x:.1f}",
        "--size_y", f"{size_y:.1f}",
        "--size_z", f"{size_z:.1f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--cnn_scoring", "rescore",
        "--seed", str(seed),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
    return parse_gnina_sdf_output(output_sdf)


def parse_gnina_sdf_output(docked_sdf_path: Path) -> Dict[str, Any]:
    """
    Parse GNINA output SDF file and extract top pose metrics:
    minimizedAffinity (kcal/mol), CNNscore (0-1), and CNNaffinity (pKd).
    """
    if Chem is None:
        raise ImportError("RDKit is required to parse GNINA SDF outputs.")

    suppl = Chem.SDMolSupplier(str(docked_sdf_path), removeHs=False)
    poses = []
    for idx, mol in enumerate(suppl):
        if mol is None:
            continue
        min_aff = float(mol.GetProp("minimizedAffinity")) if mol.HasProp("minimizedAffinity") else 0.0
        cnn_score = float(mol.GetProp("CNNscore")) if mol.HasProp("CNNscore") else 0.0
        cnn_aff = float(mol.GetProp("CNNaffinity")) if mol.HasProp("CNNaffinity") else 0.0
        poses.append(
            {
                "pose_rank": idx + 1,
                "minimized_affinity": round(min_aff, 3),
                "cnn_score": round(cnn_score, 4),
                "cnn_affinity": round(cnn_aff, 3),
                "mol": mol,
            }
        )
    if not poses:
        raise RuntimeError(f"No valid docked poses parsed from {docked_sdf_path}")

    poses.sort(key=lambda p: (p["cnn_score"], p["cnn_affinity"]), reverse=True)
    best = poses[0]
    return {
        "best_pose_rank": best["pose_rank"],
        "minimized_affinity": best["minimized_affinity"],
        "cnn_score": best["cnn_score"],
        "cnn_affinity": best["cnn_affinity"],
        "num_poses_generated": len(poses),
    }


def extract_binding_pocket_contacts(
    receptor_pdb_path: Path,
    docked_sdf_path: Path,
    curated_pocket_residues: List[int],
    is_alphafold: bool = False,
    contact_cutoff: float = 4.0,
    hbond_cutoff: float = 3.5,
) -> Dict[str, Any]:
    """
    Extract heavy-atom receptor-ligand contacts (< 4.0 A), putative polar/H-bond contacts (< 3.5 A),
    hydrophobic contacts, orthosteric pocket overlap, and mean pLDDT of contacted residues.
    """
    if Chem is None or PDBParser is None:
        raise ImportError("RDKit and Biopython are required for contact analysis.")

    suppl = Chem.SDMolSupplier(str(docked_sdf_path), removeHs=True)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise RuntimeError(f"Cannot load ligand pose from {docked_sdf_path}")

    conf = mol.GetConformer()
    lig_atoms = []
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        sym = atom.GetSymbol().upper()
        if sym == "H":
            continue
        pos = conf.GetAtomPosition(i)
        lig_atoms.append((sym, np.array([pos.x, pos.y, pos.z], dtype=float)))

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", str(receptor_pdb_path))
    rec_atoms = [a for a in structure.get_atoms() if a.element != "H"]
    ns = NeighborSearch(rec_atoms)

    contacted_residues = {}
    polar_hbond_residues = set()
    hydrophobic_residues = set()
    polar_elements = {"N", "O", "S", "P", "F"}

    for lig_sym, lig_xyz in lig_atoms:
        nearby = ns.search(lig_xyz, contact_cutoff, level="A")
        for rec_atom in nearby:
            res = rec_atom.get_parent()
            if res.id[0] != " ":
                continue
            res_num = int(res.id[1])
            res_label = f"{res.get_resname()}{res_num}"
            dist = float(np.linalg.norm(lig_xyz - rec_atom.get_coord()))
            bfac = float(rec_atom.get_bfactor())

            if res_num not in contacted_residues or dist < contacted_residues[res_num]["min_dist"]:
                contacted_residues[res_num] = {
                    "res_label": res_label,
                    "min_dist": round(dist, 2),
                    "bfactor": bfac,
                }

            rec_sym = rec_atom.element.upper()
            if lig_sym in polar_elements and rec_sym in polar_elements and dist <= hbond_cutoff:
                polar_hbond_residues.add(res_label)
            elif lig_sym == "C" and rec_sym == "C" and dist <= contact_cutoff:
                hydrophobic_residues.add(res_label)

    contact_nums = sorted(contacted_residues.keys())
    contact_labels = [contacted_residues[r]["res_label"] for r in contact_nums]
    mean_plddt = (
        float(np.mean([contacted_residues[r]["bfactor"] for r in contact_nums]))
        if (is_alphafold and contact_nums)
        else 91.5
    )

    overlap_count = len(set(contact_nums).intersection(set(curated_pocket_residues)))
    overlap_frac = (
        overlap_count / max(1, len(curated_pocket_residues))
        if curated_pocket_residues
        else 0.0
    )

    return {
        "n_pocket_contacts": len(contact_nums),
        "n_polar_hbond_contacts": len(polar_hbond_residues),
        "n_hydrophobic_contacts": len(hydrophobic_residues),
        "contact_residues": ";".join(contact_labels),
        "orthosteric_overlap_fraction": round(overlap_frac, 3),
        "contact_mean_plddt": round(mean_plddt, 2),
    }


def compute_composite_3d_prior(
    cnn_affinity: float,
    cnn_score: float,
    n_contacts: int,
    contact_plddt: float,
    aff_midpoint: float = 5.0,
    aff_scale: float = 0.75,
) -> Dict[str, float]:
    """
    Compute the Composite 3D Structural Binding Prior P_3D(m, p) in [0, 1]:
      - S_aff: Sigmoid-transformed GNINA CNNaffinity centered at physiological pKd = 5.0 (10 uM)
      - S_pose: GNINA 3D CNN pose probability score in [0, 1]
      - S_contact: Pocket contact density weighted by structural confidence (pLDDT / 100)
      - G_phys: Soft physical gate that penalizes decoys failing pose/affinity/contact checks
    """
    s_aff = float(1.0 / (1.0 + np.exp(-(cnn_affinity - aff_midpoint) / aff_scale)))
    s_pose = float(np.clip(cnn_score, 0.0, 1.0))
    s_contact = float(min(1.0, n_contacts / 8.0) * np.clip(contact_plddt / 100.0, 0.4, 1.0))

    raw_composite = 0.50 * s_aff + 0.30 * s_pose + 0.20 * s_contact

    if cnn_affinity >= 4.8 and cnn_score >= 0.45 and n_contacts >= 5:
        g_phys = 1.0
    else:
        aff_deficit = max(0.0, 4.8 - cnn_affinity)
        pose_deficit = max(0.0, 0.45 - cnn_score)
        contact_deficit = max(0.0, (5 - n_contacts) / 5.0)
        g_phys = float(np.exp(-2.5 * aff_deficit - 3.5 * pose_deficit - 2.0 * contact_deficit))

    p_3d = float(np.clip(g_phys * raw_composite, 0.0, 1.0))
    return {
        "s_affinity_norm": round(s_aff, 4),
        "s_pose_norm": round(s_pose, 4),
        "s_contact_norm": round(s_contact, 4),
        "physical_gate": round(g_phys, 4),
        "prior_3d_weight": round(p_3d, 4),
    }

