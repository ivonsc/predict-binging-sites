import os
import pandas as pd
from Bio import PDB
import mdtraj as md

# Lista de aminoácidos estándar en PDB
AMINO_ACIDS = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL"
]
def extract_aa_composition(pdb_file):
    """Extracts the amino acid composition of the pocket PDB, ignoring water and non-standard residues."""
    
    aa_counts = {aa: 0 for aa in AMINO_ACIDS}  # Initialize counts for all amino acids
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", pdb_file)

    for model in structure:
        for chain in model:
            for residue in chain:
                # Get residue name and check if it's a standard amino acid
                aa_name = residue.get_resname().strip()
                
                # Biopython may treat HOH as a heteroatom, so filter them out
                if aa_name == "HOH" or residue.id[0] != " ":  # residue.id[0] != " " means it's a heteroatom
                    continue  # Skip water molecules and heteroatoms
                
                if aa_name in AMINO_ACIDS:
                    aa_counts[aa_name] += 1  # Count only standard amino acids

    
    total = sum(aa_counts.values())
    return {aa: (aa_counts[aa] / total if total > 0 else 0) for aa in AMINO_ACIDS}  # Normalize



# Ejemplo de uso
pdb_path = input("")
features = extract_aa_composition(pdb_path)
#print("AA composition:", features)  # {'ALA': 0.1, 'ARG': 0.05, ...}

# Propiedades fisicoquímicas promedio de los AA (Ejemplo)
HYDROPHOBICITY = {"ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5, 
                  "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
                  "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
                  "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2}

def compute_pocket_hydrophobicity(aa_composition):
    """Calcula la hidrofobicidad promedio del pocket."""
    return sum(aa_composition[aa] * HYDROPHOBICITY[aa] for aa in aa_composition if aa_composition[aa] > 0.000000000000)

pocket_hydro = compute_pocket_hydrophobicity(features)
#print("\n","Hidrofobicidad:", pocket_hydro)

def compute_pocket_volume(pdb_file):
    """Calcula el volumen del pocket usando MDTraj."""
    traj = md.load(pdb_file)
    xyz = traj.xyz[0]  # Coordenadas atómicas
    min_coords = xyz.min(axis=0)
    max_coords = xyz.max(axis=0)
    volume = ((max_coords - min_coords).prod())  # Aproximación cúbica
    return volume

pocket_volume = compute_pocket_volume(pdb_path)
#print("\n","Volumen del pocket:", pocket_volume)

CHARGES = {"ARG": +1, "LYS": +1, "ASP": -1, "GLU": -1, "HIS": 0.5}

def compute_pocket_charge(aa_composition):
    """Calcula la carga neta del pocket."""
    return sum(aa_composition[aa] * CHARGES.get(aa, 0) for aa in aa_composition)

pocket_charge = compute_pocket_charge(features)
#print("\n", "Carga neta:", pocket_charge)

def extract_features(pdb_file):
    """Extrae todas las características de un pocket PDB."""
    aa_composition = extract_aa_composition(pdb_file)
    pocket_charge = compute_pocket_charge(aa_composition)
    hydrophobicity = compute_pocket_hydrophobicity(aa_composition)
    volume = compute_pocket_volume(pdb_file)

    features = {"hydrophobicity": hydrophobicity, "pocket_charge": pocket_charge, "volume": volume, **aa_composition}
    return features




# Extraer features de múltiples archivos
#df = pd.DataFrame([extract_features(pdb) for pdb in pocket_files])
df = pd.DataFrame([extract_features(pdb_path)])
final_df = df.to_csv(sep=",", index=False)
final_df= final_df.strip()
print(final_df)
'''
##### RANDOM FOREST #####

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# Supongamos que tenemos etiquetas para cada pocket
labels = ["Pocket_A", "Pocket_B", "Pocket_C"]  # Reemplaza con tus clases
df["label"] = labels

# Dividir datos en entrenamiento y prueba
X = df.drop(columns=["label"])
y = df["label"]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Entrenar Random Forest
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# Evaluación
y_pred = rf.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)'''

