import sys
from rdkit import Chem

# les atomes sont colorés par leur numéro atomique (C=6, O=8, N=7 ...)
# les liaisons par leur type (simple=101,double=102,triple=103,aromatique=104) je suis pas sur de la quatite datome jai juste mis 100 apres je corrigerai
# puis en reliant chaque liaison au deux atomes qu’elle connecte
# par exemple, pour la molécule C=O, le graphe contient trois nœuds — un nœud “C” (6), un nœud “O” (8) et un nœud “liaison double” (102) — reliés comme C — (102) — O
# ce qui permet à Nauty de comparer deux molécules en tenant compte à la fois des types d’atomes et des types de liaisons, garantissant que seules les structures chimiquement identiques soient reconnues comme isomorphes.


def export_to_graph(mol_path, out_path):
    mol = Chem.MolFromMolFile(mol_path)
    if not mol:
        return

    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()

    num_atoms = len(atoms)
    num_bonds = len(bonds)
    total_nodes = num_atoms + num_bonds

    total_edges = num_bonds * 2

    with open(out_path, 'w') as f:

        f.write(f"{total_nodes} {total_edges}\n")

        for atom in atoms:
            f.write(f"{atom.GetAtomicNum()}\n")

        for bond in bonds:
            b_type = bond.GetBondTypeAsDouble()

            color = 100 + int(b_type if b_type != 1.5 else 4)
            f.write(f"{color}\n")

        for i, bond in enumerate(bonds):
            bond_node_idx = num_atoms + i
            u = bond.GetBeginAtomIdx()
            v = bond.GetEndAtomIdx()
            f.write(f"{u} {bond_node_idx}\n")
            f.write(f"{v} {bond_node_idx}\n")


if __name__ == "__main__":
    export_to_graph(sys.argv[1], sys.argv[2])
