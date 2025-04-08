import os
import pandas as pd
from Bio import PDB
import mdtraj as md
from Bio.PDB import PDBParser, DSSP

# Lista de aminoácidos estándar en PDB
AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
]
hydrophobic_residues = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}




def clean_pdb(pdb_file, cleaned_pdb):
    """
    Removes unsupported lines from a PDB file (like TER) before DSSP processing.
    """
    with open(pdb_file, "r") as f, open(cleaned_pdb, "w") as out:
        for line in f:
            if not line.startswith("TER"):  # Ignore 'TER' lines
                out.write(line)

    return cleaned_pdb

pdb_path = input("")
pdb_path_cleaned = clean_pdb(pdb_path, pdb_path)


def compute_sasa(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    model = structure[0]  # First model
    dssp = DSSP(model, pdb_file)  # Run DSSP

    sasa_values = {}
    for key in dssp.keys():
        residue_id = key[1]
        sasa = dssp[key][3]  # Extract SASA
        sasa_values[residue_id] = sasa

    return sasa_values  # Dictionary {Residue Position: SASA}


features = compute_sasa(pdb_path_cleaned)
print(features)

def compute_secondary_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]
    dssp = DSSP(model, pdb_file)

    sec_struct = {}
    for key in dssp.keys():
        residue_id = key[1]
        ss = dssp[key][2]  # Extract secondary structure
        sec_struct[residue_id] = ss  # 'H' (Helix), 'S' (Sheet), 'L' (Loop)

    return sec_struct  # Dictionary {Residue Position: Secondary Structure}


features = compute_secondary_structure(pdb_path_cleaned)
print(features)


from Bio.PDB.ResidueDepth import ResidueDepth

def compute_residue_depth(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    model = structure[0]
    depth_calc = ResidueDepth(model)

    depth_values = {res[1]: depth_calc[res] for res in depth_calc.keys()}  # Residue Position: Depth
    return depth_values



features = compute_residue_depth(pdb_path_cleaned)
print(features)

# Hydrophobic amino acids
hydrophobic_residues = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}

# Charge properties
charge_dict = {"ASP": -1, "GLU": -1, "HIS": 1, "LYS": 1, "ARG": 1}

def compute_chem_properties(res_name):
    hydrophobicity = 1 if res_name in hydrophobic_residues else 0
    charge = charge_dict.get(res_name, 0)
    
    return [hydrophobicity, charge]



features = compute_chem_properties(pdb_path_cleaned)
print(features)


from Bio.PDB.NeighborSearch import NeighborSearch

def compute_contact_number(pdb_file, cutoff=5.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    atoms = [atom for atom in structure.get_atoms()]
    ns = NeighborSearch(atoms)

    contact_numbers = {}
    for residue in structure.get_residues():
        res_id = residue.get_id()[1]
        contacts = ns.search(residue["CA"].coord, cutoff)  # Search nearby residues
        contact_numbers[res_id] = len(contacts) - 1  # Exclude itself

    return contact_numbers


features = compute_contact_number(pdb_path_cleaned)
print(features)




import pandas as pd

pdb_file = "refined-set/187l/187l_protein.pdb"
features = extract_residue_features(pdb_file)
labeled_data = label_pockets(features, "187L")

df = pd.DataFrame(labeled_data, columns=[
    "PDB_ID", "Chain", "Residue", "Position", "B-Factor", "Pocket_Label"
])

print(df.head())
