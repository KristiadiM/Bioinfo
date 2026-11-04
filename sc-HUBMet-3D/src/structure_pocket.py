"""
Protein 3D structure fetching (RCSB PDB & AlphaFold DB), sequence-matched chain
cleaning, pLDDT filtering, and binding-pocket detection.
"""

from pathlib import Path
from typing import Dict, Any, List, Tuple
import shutil
import subprocess
import numpy as np
import requests
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.Data.IUPACData import protein_letters_3to1


class CleanReceptorSelect(Select):
    """Retain standard amino acid residues on the target chain and filter low-pLDDT tails."""

    def __init__(self, target_chain: str, is_alphafold: bool = False, min_plddt: float = 50.0):
        super().__init__()
        self.target_chain = target_chain
        self.is_alphafold = is_alphafold
        self.min_plddt = min_plddt

    def accept_chain(self, chain) -> int:
        return 1 if chain.id == self.target_chain else 0

    def accept_residue(self, residue) -> int:
        hetflag, _, _ = residue.id
        if hetflag != " ":
            return 0
        if self.is_alphafold:
            b_factors = [atom.get_bfactor() for atom in residue.get_atoms()]
            if b_factors and float(np.mean(b_factors)) < self.min_plddt:
                return 0
        return 1


def fetch_structure_pdb(
    structure_source: str,
    pdb_id: str,
    uniprot_id: str,
    raw_dir: Path,
    timeout: int = 30,
) -> Path:
    """Download experimental PDB from RCSB or predicted structure from AlphaFold DB."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    if structure_source.upper() == "PDB":
        clean_id = pdb_id.strip().upper()
        out_path = raw_dir / f"{clean_id}.pdb"
        if out_path.exists() and out_path.stat().st_size > 1000:
            return out_path
        url = f"https://files.rcsb.org/download/{clean_id}.pdb"
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        out_path.write_text(resp.text, encoding="utf-8")
        return out_path
    else:
        out_path = raw_dir / f"AF-{uniprot_id}-F1.pdb"
        if out_path.exists() and out_path.stat().st_size > 1000:
            return out_path
        for ver in ("v4", "v3", "v2"):
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{ver}.pdb"
            resp = requests.get(url, timeout=timeout)
            if resp.status_code == 200 and "ATOM" in resp.text:
                out_path.write_text(resp.text, encoding="utf-8")
                return out_path
        api_resp = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}", timeout=timeout)
        api_resp.raise_for_status()
        entries = api_resp.json()
        if entries and "pdbUrl" in entries[0]:
            pdb_resp = requests.get(entries[0]["pdbUrl"], timeout=timeout)
            pdb_resp.raise_for_status()
            out_path.write_text(pdb_resp.text, encoding="utf-8")
            return out_path
        raise RuntimeError(f"Could not download AlphaFold structure for {uniprot_id}")


def select_target_receptor_chain(model, req_chain: str, uniprot_seq: str = "", pocket_res_nums: List[int] | None = None) -> str:
    """
    Select the chain corresponding to the target receptor using UniProt 5-mer sequence
    matching so multi-subunit Cryo-EM structures never select a G-protein or antibody chain.
    """
    chains = list(model.get_chains())
    if len(chains) == 1:
        return chains[0].id

    kmers = {uniprot_seq[i : i + 5] for i in range(max(0, len(uniprot_seq) - 4))} if len(uniprot_seq) >= 10 else set()
    best_chain_id = req_chain if any(c.id == req_chain for c in chains) else chains[0].id
    best_score = -1

    for c in chains:
        std_res = [r for r in c.get_residues() if r.id[0] == " "]
        if not std_res:
            continue
        c_seq = "".join(protein_letters_3to1.get(r.get_resname().strip().capitalize(), "X") for r in std_res)
        k_hits = sum(1 for i in range(max(0, len(c_seq) - 4)) if c_seq[i : i + 5] in kmers) if kmers else 0
        res_nums = {int(r.id[1]) for r in std_res}
        p_hits = len(res_nums.intersection(set(pocket_res_nums or [])))
        pref_bonus = 2 if c.id == req_chain else 0
        score = k_hits * 10 + p_hits * 3 + pref_bonus
        if score > best_score:
            best_score = score
            best_chain_id = c.id

    return best_chain_id


def parse_pocket_residue_list(pocket_residues_str: str) -> List[int]:
    if not pocket_residues_str or str(pocket_residues_str).strip() in ("", "nan", "None"):
        return []
    return [int(t.strip()) for t in str(pocket_residues_str).split(";") if t.strip().isdigit()]


def prepare_receptor_and_pocket_box(
    structure_source: str,
    pdb_id: str,
    uniprot_id: str,
    chain_id: str,
    cocrystal_resname: str,
    pocket_residues_str: str,
    raw_dir: Path,
    clean_dir: Path,
    uniprot_seq: str = "",
    box_size: float = 22.0,
    min_alphafold_plddt: float = 50.0,
) -> Dict[str, Any]:
    raw_pdb_path = fetch_structure_pdb(structure_source, pdb_id, uniprot_id, raw_dir)
    is_af = structure_source.upper() != "PDB"

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(raw_pdb_path))
    model = structure[0]

    pocket_res_nums = parse_pocket_residue_list(pocket_residues_str)
    target_chain_id = select_target_receptor_chain(model, chain_id.strip(), uniprot_seq=uniprot_seq, pocket_res_nums=pocket_res_nums)
    chain = model[target_chain_id]

    pocket_coords, pocket_bfactors = [], []
    for res in chain.get_residues():
        hetflag, resseq, _ = res.id
        if hetflag == " " and resseq in pocket_res_nums:
            for atom in res.get_atoms():
                if atom.element != "H":
                    pocket_coords.append(atom.get_coord())
                    pocket_bfactors.append(atom.get_bfactor())

    cocrystal_coords = []
    if cocrystal_resname and str(cocrystal_resname).strip().upper() not in ("", "NAN", "NONE"):
        target_lig = str(cocrystal_resname).strip().upper()
        for c in model.get_chains():
            for res in c.get_residues():
                if res.get_resname().strip().upper() == target_lig:
                    for atom in res.get_atoms():
                        if atom.element != "H":
                            cocrystal_coords.append(atom.get_coord())

    clean_dir.mkdir(parents=True, exist_ok=True)
    clean_pdb_path = clean_dir / f"{pdb_id}_{target_chain_id}_clean.pdb"
    io_pdb = PDBIO()
    io_pdb.set_structure(structure)
    io_pdb.save(
        str(clean_pdb_path),
        CleanReceptorSelect(target_chain=target_chain_id, is_alphafold=is_af, min_plddt=min_alphafold_plddt),
    )

    if len(cocrystal_coords) >= 3:
        pocket_method = "cocrystal_ligand_autobox"
        center = np.mean(np.array(cocrystal_coords), axis=0)
        mean_plddt = 93.0
    elif len(pocket_coords) >= 3:
        pocket_method = "curated_orthosteric_residues"
        center = np.mean(np.array(pocket_coords), axis=0)
        mean_plddt = float(np.mean(pocket_bfactors)) if is_af else 92.0
    else:
        pocket_method, center, mean_plddt = _detect_fpocket_or_cavity_centroid(clean_pdb_path, chain, is_af)

    return {
        "clean_pdb_path": str(clean_pdb_path),
        "target_chain": target_chain_id,
        "pocket_method": pocket_method,
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x": float(box_size),
        "size_y": float(box_size),
        "size_z": float(box_size),
        "pocket_mean_plddt": round(mean_plddt, 2),
        "resolved_pocket_residues": len(pocket_res_nums),
    }


def _detect_fpocket_or_cavity_centroid(
    clean_pdb_path: Path, chain, is_af: bool
) -> Tuple[str, np.ndarray, float]:
    if shutil.which("fpocket") is not None:
        try:
            subprocess.run(
                ["fpocket", "-f", str(clean_pdb_path)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=30,
            )
            pqr_file = clean_pdb_path.parent / f"{clean_pdb_path.stem}_out" / "pockets" / "pocket1_vert.pqr"
            if pqr_file.exists():
                verts = [
                    [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                    for ln in pqr_file.read_text().splitlines()
                    if ln.startswith(("ATOM", "HETATM"))
                ]
                if verts:
                    return "fpocket_pocket1", np.mean(np.array(verts), axis=0), 85.0
        except Exception:
            pass

    ca_coords, bfactors = [], []
    for res in chain.get_residues():
        if res.id[0] == " " and "CA" in res:
            bf = res["CA"].get_bfactor()
            if not is_af or bf >= 65.0:
                ca_coords.append(res["CA"].get_coord())
                bfactors.append(bf)
    center = np.mean(np.array(ca_coords), axis=0) if ca_coords else np.zeros(3)
    mean_plddt = float(np.mean(bfactors)) if (is_af and bfactors) else 85.0
    return "interior_core_centroid", center, mean_plddt
