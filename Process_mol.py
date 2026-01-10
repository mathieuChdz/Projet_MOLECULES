import sys
from rdkit import Chem

def export_to_graph(mol_path, out_path):
    mol = Chem.MolFromMolFile(mol_path)
    if not mol: return
    
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    
    with open(out_path, 'w') as f:
        f.write(f"{len(atoms)} {len(bonds)}\n")
        for atom in atoms:
            f.write(f"{atom.GetIdx()} {atom.GetAtomicNum()}\n")
        for bond in bonds:
            f.write(f"{bond.GetBeginAtomIdx()} {bond.GetEndAtomIdx()}\n")

if __name__ == "__main__":
    export_to_graph(sys.argv[1], sys.argv[2])