import os
from rdkit import Chem
from rdkit.Chem import AllChem  # Import AllChem for force fields
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

def process_ligands_in_folder(input_folder, output_folder, charge_model="gasteiger", merge_nonpolar_hydrogens=False):
    """
    Automates hydrogen addition (if missing), charge assignment,
    energy minimization, and conversion from SDF to PDBQT format
    for all SDF files in a given folder, saving the results to a new folder.

    Args:
        input_folder (str): Path to the folder containing input SDF files.
        output_folder (str): Path to the folder where output PDBQT files will be saved.
        charge_model (str): The charge model to use ("gasteiger" is directly supported by Meeko).
        merge_nonpolar_hydrogens (bool): If True, nonpolar hydrogens will be merged.
                                          Set to False for more accurate docking results (recommended).
    """
    if not os.path.exists(input_folder):
        print(f"Error: Input folder not found at {input_folder}")
        return

    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    print(f"Output folder '{output_folder}' ensured.")

    # Configure Meeko's MoleculePreparation
    if merge_nonpolar_hydrogens:
        merge_types = ("H",)
        print("Nonpolar hydrogens will be merged for all ligands.")
    else:
        merge_types = ()
        print("All hydrogens will be kept (not merged) for all ligands for accuracy.")

    preparator = MoleculePreparation(
        merge_these_atom_types=merge_types,
        charge_model=charge_model
    )

    processed_count = 0
    failed_count = 0

    # Iterate through all files in the input folder
    for filename in os.listdir(input_folder):
        if filename.lower().endswith(".sdf"):
            input_sdf_path = os.path.join(input_folder, filename)
            # Construct output filename: replace .sdf with .pdbqt
            output_pdbqt_filename = os.path.splitext(filename)[0] + ".pdbqt"
            output_pdbqt_path = os.path.join(output_folder, output_pdbqt_filename)

            print(f"\nProcessing: {filename}")

            try:
                # Read the SDF file. removeHs=False ensures existing hydrogens are kept.
                # Use Chem.SDMolSupplier to handle potentially multiple molecules or issues
                suppl = Chem.SDMolSupplier(input_sdf_path, removeHs=False)
                mol = None
                for m in suppl:
                    if m is not None:
                        mol = m
                        break # Assuming one ligand per SDF file for this script

                if mol is None:
                    print(f"  Warning: Could not read any molecule from {filename}. Skipping.")
                    failed_count += 1
                    continue

                # Add hydrogens if they are missing and ensure 3D coordinates for minimization
                mol = Chem.AddHs(mol, addCoords=True)
                print("  Hydrogens ensured (added if missing).")

                # --- Energy Minimization Step ---
                # Check if the molecule already has 3D coordinates. If not, generate them.
                # If it already has, we'll optimize those.
                if mol.GetNumConformers() == 0:
                    print("  Generating 3D coordinates...")
                    # Generate 3D coordinates if not present. MaxAttempts can be adjusted.
                    # randomSeed ensures reproducibility if you always use the same seed
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG(), randomSeed=42)
                    if mol.GetNumConformers() == 0:
                        print(f"  Error: Could not generate 3D coordinates for {filename}. Skipping.")
                        failed_count += 1
                        continue
                else:
                    print("  Existing 3D coordinates found.")

                print("  Performing energy minimization...")
                # Optimize the molecule using the Universal Force Field (UFF)
                # maxIters can be adjusted. A return value of 0 indicates successful optimization.
                ff_status = AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
                if ff_status == 0:
                    print("  Energy minimization completed successfully.")
                else:
                    print(f"  Warning: Energy minimization did not converge for {filename} (status: {ff_status}).")
                # --- End Energy Minimization Step ---

                # Prepare the molecule with Meeko
                mol_setups = preparator.prepare(mol)
                if not mol_setups:
                    print(f"  Error: Meeko could not prepare molecule from {filename}. Skipping.")
                    failed_count += 1
                    continue

                mol_setup = mol_setups[0] # Take the first conformer setup
                print(f"  Molecule prepared with {charge_model} charges.")

                # Write the PDBQT file
                pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setup)

                if is_ok:
                    with open(output_pdbqt_path, 'w') as f:
                        f.write(pdbqt_string)
                    print(f"  Successfully converted to: {output_pdbqt_filename}")
                    processed_count += 1
                else:
                    print(f"  Error writing PDBQT for {filename}: {error_msg}. Skipping.")
                    failed_count += 1

            except Exception as e:
                print(f"  An unexpected error occurred while processing {filename}: {e}. Skipping.")
                failed_count += 1

    print(f"\n--- Conversion Summary ---")
    print(f"Total SDF files found: {processed_count + failed_count}")
    print(f"Successfully converted: {processed_count}")
    print(f"Failed conversions: {failed_count}")
    print(f"PDBQT files saved to: {os.path.abspath(output_folder)}")

if __name__ == "__main__":
    # --- Configuration ---
    input_sdf_directory = r"C:\Python\VinaDock\sdf_files"
    output_pdbqt_directory = r"C:\Python\VinaDock\pdbqt_ligands_minim"

    # Meeko defaults to "gasteiger" charges, which is good for general use.
    chosen_charge_model = "gasteiger"

    # Set to False to keep all hydrogens (recommended for accuracy)
    # Set to True to merge nonpolar hydrogens (faster, less accurate)
    should_merge_nonpolar = False

    # Run the batch processing
    process_ligands_in_folder(
        input_sdf_directory,
        output_pdbqt_directory,
        charge_model=chosen_charge_model,
        merge_nonpolar_hydrogens=should_merge_nonpolar
    )

    print("\nBatch processing finished.")
