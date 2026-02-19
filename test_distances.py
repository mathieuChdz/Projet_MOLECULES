import os
from rdkit import Chem
from Distance import Tanimoto2D, ShapeTanimoto, score


def run_test():

    MOL_DIR = r"data/molecules"
    mol_files = [f for f in os.listdir(MOL_DIR) if f.endswith(".mol")]

    ref_mol = os.path.join(MOL_DIR, mol_files[0])
    target_mols = [os.path.join(MOL_DIR, f) for f in mol_files[1:5]]

    print(f"--- Molécule de référence : {mol_files[0]} ---")

    print("\n--- Test de Tanimoto2D ---")
    for target in target_mols:
        val = Tanimoto2D(ref_mol, target)
        print(f"Sim(Ref, {os.path.basename(target)}) = {val:.4f}")

    print("\n--- Test de ShapeTanimoto ---")
    for target in target_mols:
        val = ShapeTanimoto(ref_mol, target, n=10)
        print(f"ShapeSim(Ref, {os.path.basename(target)}) = {val:.4f}")

    print("\n--- Test de score hybride (a=0.5) ---")
    a = 0.5
    for target in target_mols:
        val = score(ref_mol, target, a)
        print(f"Score(Ref, {os.path.basename(target)}) = {val:.4f}")


if __name__ == "__main__":
    run_test()
