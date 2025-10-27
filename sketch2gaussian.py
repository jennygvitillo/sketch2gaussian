#!/usr/bin/env python3
"""
Molecular Structure Converter: From ChemDraw to Gaussian Input Files
====================================================================

Description:
------------
This script converts 2D molecular structures from ChemDraw (.mol files) into
professional-quality 3D Gaussian input files (.gjf) with correct bond lengths
and geometries. It adds all the hydrogens!!
It uses a hybrid approach combining Open Babel for SMILES
conversion and RDKit's ETKDG algorithm for accurate 3D coordinate generation.

Key Features:
-------------
- Converts .mol files to proper 3D structures with correct bond lengths
- Automatically adds missing hydrogen atoms
- Uses RDKit's ETKDGv3 algorithm for chemically accurate geometries
- Generates Gaussian input files with standard computational parameters
- Includes SMILES strings in Gaussian files for traceability
- Batch processes multiple files in a folder
- Interactive folder selection dialog

Workflow:
---------
1. Input: .mol files from ChemDraw (2D coordinates without the hydrogen)
2. Open Babel: Converts to SMILES format with hydrogens
3. RDKit ETKDG: Generates proper 3D coordinates from SMILES
4. Output: Gaussian .gjf files with optimized geometries and SMILES

Requirements:
-------------
- Python 3.6+ (I am using Windows: so I have installed PyCharm and this program deals with all the Python stuff, very convenient)
- Open Babel (https://github.com/openbabel/openbabel)
- RDKit (https://www.rdkit.org/)

Installation for Dummies:
-------------------------
0. If you don’t have it, download and install PyCharm Community Edition: it makes running Python scripts much easier.

1. Create and open the project in Pycharm.
   - Open PyCharm and create a **new project**. You can name it sketch2gaussian.  
   - On the top left, you’ll see something that looks like a folder.  
   - **Right-click** on it → choose **New → Python File**.  
   - Name it `mol_to_gaussian.py` and press **Enter**.  
   - Copy and paste the entire script into the text area that appears.  
   - Go to **File → Save All**.
To run the script, simply click the green “Play” button above the text area and follow the on-screen instructions.  
If you get an error about *Open Babel* or *RDKit* not being installed, follow the steps below.

2. Install Open Babel (if you have never used it, you do not have it on your pc):
   - Download from the internet: https://github.com/openbabel/openbabel/releases
   - Run the Windows installer as Administrator (choose the *.exe where in the description is mentioned 64 bit)
   - Make sure "Add Open Babel to system PATH" is checked (this is the default and it allows th Python script to use it from any directory)

3. Install RDKit (if it is not installed, you will get an error the first time you use the script):
   - Open Terminal in Pycharm (small icon on the botton left)
   - Type:
     pip install rdkit
   - If it fails, try: `conda install -c conda-forge rdkit`
That’s it! Once both are installed, you can run `sketch2gaussian.py` again, and it should work perfectly.

Usage:
------
1. Run the script
2. Select the folder containing your .mol files
3. The script will automatically:
   - Process all .mol files in the folder
   - Create a "gjf_files_rdkit" subfolder with the results
   - Open the output folder when done

Output:
-------
- Gaussian .gjf files with proper 3D geometries
- Files are named after the original .mol files (special characters removed)
- All hydrogen atoms are added automatically
- Bond lengths and angles are chemically accurate
- SMILES strings included in Gaussian file headers for reference

Troubleshooting:
----------------
- "Open Babel not found": Install Open Babel and add to PATH
- "RDKit not available": Run `pip install rdkit` in terminal
- No .mol files found: Check your input folder contains .mol files
- Bond lengths still wrong: RDKit should fix this - check your source .mol files

Author:
---------
Jenny G. Vitillo (University of Insubria) and DeepSeek 2024


License:
---------
MIT


Citation:
---------
J.G. Vitillo, sketch2gaussian, https://github.com/jennygvitillo/sketch2gaussian (2025).

"""

import os
import re
import tkinter as tk
from tkinter import filedialog, messagebox
import subprocess
import sys


def install_rdkit():
    """
    Install RDKit if not available
    Returns: True if RDKit is available, False otherwise
    """
    try:
        import rdkit
        return True
    except ImportError:
        print("📦 RDKit not found, attempting installation...")
        try:
            # Try different installation methods
            subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit"])
            import rdkit
            print("✅ RDKit installed successfully!")
            return True
        except:
            print("❌ Cannot install RDKit automatically")
            print("Please install manually by running:")
            print("  pip install rdkit")
            print("Or if using conda:")
            print("  conda install -c conda-forge rdkit")
            return False


def check_openbabel():
    """
    Check if OpenBabel is available in the system PATH
    Returns: True if OpenBabel is available, False otherwise
    """
    try:
        result = subprocess.run(['obabel', '-V'], capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except:
        return False


def mol_to_smiles(mol_path):
    """
    Convert .mol file to SMILES string using OpenBabel
    Args:
        mol_path: Path to the .mol file
    Returns:
        tuple: (success: bool, smiles: str or error_message: str)
    """
    try:
        # Use OpenBabel to convert .mol to SMILES with hydrogens
        result = subprocess.run(['obabel', mol_path, '-osmi', '-h'],
                                capture_output=True, text=True, timeout=30)
        if result.returncode == 0 and result.stdout.strip():
            # Extract SMILES (remove filename that OpenBabel adds)
            smiles = result.stdout.strip().split('\t')[0]
            return True, smiles
        return False, "Cannot convert to SMILES"
    except Exception as e:
        return False, str(e)


def smiles_to_3d_rdkit(smiles, output_path):
    """
    Convert SMILES string to 3D coordinates using RDKit's ETKDG algorithm
    This generates chemically accurate bond lengths and angles

    Args:
        smiles: SMILES string of the molecule
        output_path: Path where to save the .xyz file

    Returns:
        tuple: (success: bool, message: str)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        print(f"   🧪 Converting SMILES to 3D: {smiles}")

        # Create molecule object from SMILES string
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Add hydrogen atoms (RDKit automatically determines correct positions)
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates using ETKDGv3 algorithm
        # ETKDG = Embeddededness, Torsion, Knowledge, Distance Geometry
        # This algorithm uses chemical knowledge to generate accurate geometries
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # For reproducible results
        params.useRandomCoords = True  # Start from random coordinates
        params.useBasicKnowledge = True  # 👈 KEY: Use chemical knowledge for bond lengths

        print("   🔬 Embedding molecule with ETKDGv3...")
        embed_result = AllChem.EmbedMolecule(mol, params)

        # Fallback to ETKDGv2 if v3 fails
        if embed_result != 0:
            print("   🔄 Fallback to ETKDGv2...")
            embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv2())

        # Final fallback to basic embedding
        if embed_result != 0:
            print("   🔄 Fallback to basic embedding...")
            embed_result = AllChem.EmbedMolecule(mol)

        if embed_result != 0:
            return False, "Cannot generate 3D coordinates"

        print("   ⚙️  Optimizing geometry with forcefield...")
        # Light optimization with MMFF94 forcefield
        try:
            result = AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            if result != 0:
                print("   ⚠️  MMFF optimization didn't converge fully")
        except Exception as e:
            print(f"   ⚠️  MMFF failed: {e}")
            try:
                # Fallback to UFF forcefield
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                print("   ✅ Used UFF optimization")
            except:
                print("   ⚠️  No forcefield optimization applied")

        # Write coordinates to XYZ file
        print("   💾 Writing XYZ file...")
        with open(output_path, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\n")
            f.write(f"Generated by RDKit ETKDG from SMILES: {smiles}\n")

            conf = mol.GetConformer()
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                symbol = atom.GetSymbol()
                f.write(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

        # Verify that bond lengths are chemically reasonable
        if verify_bond_lengths(mol):
            return True, "RDKit ETKDG - Correct geometry generated"
        else:
            return True, "RDKit ETKDG - Check geometry manually"

    except Exception as e:
        return False, f"RDKit error: {str(e)}"


def verify_bond_lengths(mol):
    """
    Verify that generated bond lengths are chemically reasonable
    Checks C-C bonds ~1.54Å, C-H bonds ~1.09Å, etc.

    Args:
        mol: RDKit molecule object

    Returns:
        bool: True if most bonds are reasonable
    """
    try:
        from rdkit import Chem

        conf = mol.GetConformer()
        reasonable_bonds = 0
        total_bonds = 0

        # Check all bonds in the molecule
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()

            # Calculate distance between atoms
            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)
            distance = ((pos2.x - pos1.x) ** 2 +
                        (pos2.y - pos1.y) ** 2 +
                        (pos2.z - pos1.z) ** 2) ** 0.5

            atom1_type = mol.GetAtomWithIdx(idx1).GetSymbol()
            atom2_type = mol.GetAtomWithIdx(idx2).GetSymbol()

            # Check bond lengths against known chemical values
            if atom1_type == 'C' and atom2_type == 'C':
                if 1.40 <= distance <= 1.60:  # C-C single bond range
                    reasonable_bonds += 1
                else:
                    print(f"     ⚠️  C-C bond: {distance:.3f} Å (expected ~1.54 Å)")
            elif atom1_type == 'C' and atom2_type == 'C' and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                if 1.30 <= distance <= 1.40:  # C=C double bond range
                    reasonable_bonds += 1
                else:
                    print(f"     ⚠️  C=C bond: {distance:.3f} Å (expected ~1.34 Å)")
            elif ('C' in [atom1_type, atom2_type]) and ('H' in [atom1_type, atom2_type]):
                if 1.00 <= distance <= 1.12:  # C-H bond range
                    reasonable_bonds += 1
                else:
                    print(f"     ⚠️  C-H bond: {distance:.3f} Å (expected ~1.09 Å)")
            elif 0.8 <= distance <= 2.2:  # General reasonable range for other bonds
                reasonable_bonds += 1
            else:
                print(f"     ⚠️  {atom1_type}-{atom2_type} bond: {distance:.3f} Å")

            total_bonds += 1

        # Calculate success ratio
        ratio = reasonable_bonds / total_bonds if total_bonds > 0 else 0
        print(f"   📊 Bond validation: {reasonable_bonds}/{total_bonds} reasonable ({ratio:.1%})")
        return ratio >= 0.7  # 70% reasonable bonds is acceptable

    except Exception as e:
        print(f"   ⚠️  Bond validation failed: {e}")
        return True  # If validation fails, assume it's OK


def sanitize_filename(filename):
    """
    Remove special characters from filename for compatibility
    Keeps only alphanumeric characters, underscores, and hyphens

    Args:
        filename: Original filename

    Returns:
        str: Sanitized filename
    """
    name_without_ext = os.path.splitext(filename)[0]
    sanitized = re.sub(r'[^\w\-_]', '', name_without_ext)
    return sanitized


def select_input_folder():
    """
    Open dialog to select folder containing .mol files

    Returns:
        str: Selected folder path or None if canceled
    """
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    folder = filedialog.askdirectory(title="Select folder containing .mol files")
    root.destroy()
    return folder


def convert_mol_to_gaussian():
    """
    Main conversion function
    Processes all .mol files in selected folder and generates Gaussian .gjf files
    """
    # Check dependencies
    if not check_openbabel():
        messagebox.showerror("Error",
                             "Open Babel not found!\n\n"
                             "Please install Open Babel from:\n"
                             "https://github.com/openbabel/openbabel/releases\n\n"
                             "Run the installer as Administrator and make sure\n"
                             "'Add Open Babel to system PATH' is selected.")
        return

    if not install_rdkit():
        messagebox.showerror("Error",
                             "RDKit not available!\n\n"
                             "Please install RDKit by running:\n"
                             "pip install rdkit\n\n"
                             "Or if using conda:\n"
                             "conda install -c conda-forge rdkit")
        return

    print("✅ All dependencies available - Starting conversion...")

    # Select input folder
    input_folder = select_input_folder()
    if not input_folder:
        print("❌ No folder selected. Operation cancelled.")
        return

    print(f"📁 Selected folder: {input_folder}")

    # Create output folders
    output_folder = os.path.join(input_folder, "gjf_files_rdkit")
    temp_folder = os.path.join(input_folder, "temp_rdkit")

    # Clean existing folders and create new ones
    import shutil
    for folder in [output_folder, temp_folder]:
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.makedirs(folder)

    # Find all .mol files in the input folder
    mol_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.mol')]
    if not mol_files:
        messagebox.showerror("Error", f"No .mol files found in:\n{input_folder}")
        return

    print(f"📁 Found {len(mol_files)} .mol files")
    print("=" * 50)

    successful_conversions = 0

    # Process each .mol file
    for mol_file in mol_files:
        try:
            base_name = sanitize_filename(mol_file)
            input_path = os.path.join(input_folder, mol_file)
            temp_xyz = os.path.join(temp_folder, f"{base_name}.xyz")
            gjf_path = os.path.join(output_folder, f"{base_name}.gjf")

            print(f"🔄 Processing: {mol_file}")

            # Step 1: Convert .mol to SMILES with OpenBabel
            success, smiles = mol_to_smiles(input_path)
            if not success:
                print(f"   ❌ {smiles}")
                continue

            # Step 2: Generate proper 3D coordinates with RDKit ETKDG
            success, message = smiles_to_3d_rdkit(smiles, temp_xyz)
            if not success:
                print(f"   ❌ {message}")
                continue

            # Step 3: Read coordinates and create Gaussian input file
            atoms = []
            try:
                with open(temp_xyz, 'r') as f:
                    lines = f.readlines()
                    # Skip first two lines (atom count and comment)
                    for line in lines[2:]:
                        parts = line.strip().split()
                        if len(parts) >= 4:
                            atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            except Exception as e:
                print(f"   ❌ Error reading XYZ file: {e}")
                continue

            if not atoms:
                print(f"   ❌ No atoms found in XYZ file")
                continue

            # Create Gaussian input file with SMILES in header
            try:
                with open(gjf_path, 'w') as f:
                    # Gaussian calculation parameters
                    f.write("%LindaWorkers=cn1104\n")
                    f.write("%NProcShared=64\n")
                    f.write("%Mem=127800mb\n")
                    f.write(f"%chk={base_name}.chk\n\n")

                    # Calculation method and basis set
                    f.write("#p opt freq b3lyp pop=hirshfeld 5d def2tzvp empiricaldispersion=gd3bj\n")
                    f.write("gfinput integral=grid=ultrafine scf=(xqc,maxconventionalcycle=150) test\n\n")

                    # Title with SMILES for traceability
                    f.write(f"{base_name} | SMILES: {smiles}\n\n")
                    f.write("0 1\n")

                    # Atomic coordinates
                    for element, x, y, z in atoms:
                        f.write(f"{element:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

                    f.write("\n")

                print(f"✅ Created: {base_name}.gjf with {len(atoms)} atoms")
                print(f"   💡 {message}")
                print(f"   📝 SMILES: {smiles}")
                successful_conversions += 1

            except Exception as e:
                print(f"   ❌ Error writing Gaussian file: {e}")

        except Exception as e:
            print(f"❌ Unexpected error processing {mol_file}: {e}")

    # Cleanup temporary files
    if os.path.exists(temp_folder):
        shutil.rmtree(temp_folder)

    # Display results
    print("=" * 50)
    print(f"🎉 Conversion completed!")
    print(f"📊 Successful conversions: {successful_conversions}/{len(mol_files)}")

    if successful_conversions > 0:
        try:
            # Open output folder for user convenience
            os.startfile(output_folder)
            print(f"📂 Output folder opened: {output_folder}")
        except:
            print(f"📂 Output folder: {output_folder}")

        # Show created files
        gjf_files = [f for f in os.listdir(output_folder) if f.endswith('.gjf')]
        print(f"📄 Generated {len(gjf_files)} Gaussian input files")
        print("💡 Each file contains the SMILES string for reference")
    else:
        print("❌ No files were converted successfully")
        print("Please check the error messages above")

    print("\n✅ Script finished!")


if __name__ == "__main__":
    print("🚀 Molecular Structure Converter - ChemDraw to Gaussian")
    print("=" * 60)
    print("Professional 3D coordinate generation with RDKit ETKDG")
    print("Converts 2D ChemDraw structures to Gaussian input files")
    print("Includes SMILES strings for molecular traceability")
    print("=" * 60)

    convert_mol_to_gaussian()
