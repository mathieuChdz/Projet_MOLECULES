import os
import sys
import requests
import time
import subprocess
import re
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from jinja2 import Environment, FileSystemLoader
import argparse
import json

# Import modules
from Distance import vecteur_fingerprint2D, score, genere_shape

# Configuration
NAUTY_DIR = "nauty2_9_3"
DATA_DIR = "data"
MOL_DIR = os.path.join(DATA_DIR, "molecules")
GRAPH_DIR = os.path.join(DATA_DIR, "graphs")
IMG_DIR = os.path.join(DATA_DIR, "images")
REPORT_FILE = "resultats.html"
TEMPLATE_FILE = "template.html"  # Nom du fichier template
EXEC_CHECK_ISO = "./Check_iso"
GROUPEMENT_DIR = os.path.join(DATA_DIR, "groupement")


def setup():
    """Supprime tout le dossier data et recrée l'arborescence à neuf."""
    if os.path.exists(DATA_DIR):
        print(f"[*] Nettoyage complet du dossier {DATA_DIR}...")
        try:
            shutil.rmtree(DATA_DIR)
        except OSError as e:
            print(f"[!] Erreur lors du nettoyage : {e}")

    for folder in [DATA_DIR, MOL_DIR, GRAPH_DIR, IMG_DIR, GROUPEMENT_DIR]:
        os.makedirs(folder, exist_ok=True)


def download_data(url):
    print(f"[*] Téléchargement depuis : {url}")
    path = os.path.join(DATA_DIR, "source.sdf")
    try:
        r = requests.get(url, stream=True)
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))

        with open(path, 'wb') as f, tqdm(
            desc="Téléchargement",
            total=total_size,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for data in r.iter_content(chunk_size=1024):
                size = f.write(data)
                bar.update(size)
        return path
    except Exception as e:
        print(f"[!] Erreur de téléchargement : {e}")
        sys.exit(1)


def safe_filename(s):
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^a-zA-Z0-9_\-]", "", s)
    return s[:80]


def convert_mol_to_graph(mol_path, output_path):
    mol = Chem.MolFromMolFile(mol_path)
    if not mol:
        return
    Chem.RemoveStereochemistry(mol)
    new_order = list(Chem.CanonicalRankAtoms(mol, breakTies=True))
    mol = Chem.RenumberAtoms(mol, new_order)
    # 3. Standardiser les Hydrogènes
    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol, addCoords=False)
    num_atoms = mol.GetNumAtoms()
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()

    with open(output_path, "w") as f:
        f.write(f"{num_atoms} {len(bonds)}\n")
        atom_types = [str(a.GetAtomicNum()) for a in atoms]
        f.write(" ".join(atom_types) + "\n")
        for bond in bonds:
            f.write(f"{bond.GetBeginAtomIdx()} {bond.GetEndAtomIdx()}\n")


def split_molecules(sdf_path):
    print("[*] Découpage des molécules...")
    suppl = Chem.SDMolSupplier(sdf_path)
    conversion_args = []

    for i, mol in enumerate(tqdm(suppl, desc="Extraction SDF", unit="mol")):
        if mol is None:
            continue

        chebi_id = mol.GetProp("ChEBI ID") if mol.HasProp("ChEBI ID") else ""
        chebi_name = mol.GetProp("ChEBI NAME") if mol.HasProp("ChEBI NAME") else ""
        pubchem_cid = mol.GetProp("PUBCHEM_COMPOUND_CID") if mol.HasProp("PUBCHEM_COMPOUND_CID") else ""
        pubchem_iupac = mol.GetProp("PUBCHEM_IUPAC_NAME") if mol.HasProp("PUBCHEM_IUPAC_NAME") else ""
        sdf_name = mol.GetProp('_Name') if mol.HasProp('_Name') else ""

        id_part = safe_filename(chebi_id) or safe_filename(pubchem_cid) or safe_filename(sdf_name) or "UNKNOWN_ID"
        name_part = safe_filename(chebi_name) or safe_filename(pubchem_iupac) or safe_filename(sdf_name) or "UNKNOWN_NAME"

        # Format forcé : mol_<index>_<id>_<name>
        out_name = f"mol_{i}_{id_part}_{name_part}"

        mol_path = os.path.join(MOL_DIR, f"{out_name}.mol")
        img_path = os.path.join(IMG_DIR, f"{out_name}.png")
        graph_path = os.path.join(GRAPH_DIR, f"{out_name}.graph")

        with open(mol_path, "w") as f:
            f.write(Chem.MolToMolBlock(mol))

        try:
            AllChem.Compute2DCoords(mol)
            Draw.MolToFile(mol, img_path, size=(300, 300))
        except:
            pass

        conversion_args.append((mol_path, graph_path))

    return conversion_args


def process_conversions(conversion_args):
    total = len(conversion_args)
    print(f"[*] Conversion de {total} molécules en graphes...")

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(convert_mol_to_graph, *args)
                   for args in conversion_args]
        for _ in tqdm(as_completed(futures), total=total, desc="Génération Graphes", unit="graphe"):
            pass


def find_isomorphs():
    print("[*] Recherche d'isomorphes (Nauty)...")

    if not os.path.exists(EXEC_CHECK_ISO):
        print(f"[!] Erreur : L'exécutable {EXEC_CHECK_ISO} est introuvable.")
        sys.exit(1)

    signatures = {}
    all_molecules = []
    graphs = [f for f in os.listdir(GRAPH_DIR) if f.endswith(".graph")]

    for g in tqdm(graphs, desc="Analyse Nauty", unit="mol"):
        mol_name = g.replace(".graph", "")
        all_molecules.append(mol_name)
        path = os.path.join(GRAPH_DIR, g)

        try:
            res = subprocess.run([EXEC_CHECK_ISO, path],
                                 capture_output=True, text=True)
            sig = res.stdout.strip()
            if sig:
                if sig not in signatures:
                    signatures[sig] = []
                signatures[sig].append(mol_name)
        except Exception:
            pass

    groups = [mols for mols in signatures.values() if len(mols) > 1]
    groups.sort(key=len, reverse=True)
    return groups, sorted(all_molecules), signatures


def save_groups_json(signatures):
    """Sauvegarde les groupes d'isomorphes dans data/groupement/groupes_isomorphes.json."""
    os.makedirs(GROUPEMENT_DIR, exist_ok=True)
    output_path = os.path.join(GROUPEMENT_DIR, "groupes_isomorphes.json")

    multi_groups = [sorted(mols) for mols in signatures.values() if len(mols) > 1]
    multi_groups.sort(key=len, reverse=True)
    singletons = sorted([mols[0] for mols in signatures.values() if len(mols) == 1])

    payload = {
        "date": time.strftime("%d/%m/%Y %H:%M:%S"),
        "molecules_analysees": sum(len(v) for v in signatures.values()),
        "nombre_total_groupes": len(signatures),
        "groupes_isomorphes": [
            {
                "id": gid,
                "taille": len(members),
                "molecules": members,
            }
            for gid, members in enumerate(multi_groups, 1)
        ],
        "singletons": singletons,
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    print(f"[*] Groupes sauvegardés en JSON : {output_path}")


def compute_similarity_matrix(all_molecules):
    n = len(all_molecules)
    print(f"[*] Calcul Similarité pour {n} molécules...")

    precomputed_data = []

    for name in tqdm(all_molecules, desc="Génération 3D & FP", unit="mol"):
        path = os.path.join(MOL_DIR, f"{name}.mol")
        vec = vecteur_fingerprint2D(path)
        mol_3d, ids = genere_shape(path)
        precomputed_data.append((vec, mol_3d, ids))

    print(f"\n[*] Calcul de la matrice croisée...")
    matrix = {}
    total_ops = (n * (n + 1)) // 2

    with tqdm(total=total_ops, desc="Calcul Distances", unit="paire") as pbar:
        for i in range(n):
            name_i = all_molecules[i]
            matrix[name_i] = {}
            data_i = precomputed_data[i]

            for j in range(n):
                name_j = all_molecules[j]

                if i == j:
                    matrix[name_i][name_j] = 1.0
                elif j < i:
                    matrix[name_i][name_j] = matrix[name_j][name_i]
                else:
                    data_j = precomputed_data[j]
                    try:
                        val = score(data_i, data_j)
                        matrix[name_i][name_j] = val
                    except Exception:
                        matrix[name_i][name_j] = 0.0
                    pbar.update(1)

    return matrix


def generate_report(groups, all_molecules, similarity_matrix):
    print("[*] Génération du rapport HTML via Jinja2...")
    now = time.strftime("%d/%m/%Y %H:%M:%S")

    # Configuration de Jinja2 pour charger le fichier template.html
    file_loader = FileSystemLoader('.')
    env = Environment(loader=file_loader)

    try:
        template = env.get_template(TEMPLATE_FILE)
    except Exception as e:
        print(f"[!] Erreur : Impossible de trouver '{
              TEMPLATE_FILE}'. Vérifie qu'il est dans le dossier.")
        sys.exit(1)

    # On prépare les données à envoyer au template
    # Jinja fera le travail de boucle et d'affichage
    output = template.render(
        now=now,
        groups=groups,
        all_molecules=all_molecules,
        matrix=similarity_matrix
    )

    # Écriture du fichier final
    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write(output)

    print(f"\n[*] Rapport généré : {os.path.abspath(REPORT_FILE)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyse structurelle et isomorphisme de molécules à partir d'un fichier SDF."
    )
    parser.add_argument(
        "input",
        type=str,
        help="Chemin vers un fichier SDF local OU lien de téléchargement (http/https)."
    )

    args = parser.parse_args()
    input_arg = args.input

    setup()

    if input_arg.startswith("http://") or input_arg.startswith("https://"):
        source_path = download_data(input_arg)
    elif os.path.isfile(input_arg):
        print(f"[*] Utilisation du fichier local : {input_arg}")
        if not os.path.isdir(DATA_DIR):
            os.mkdir(DATA_DIR)
        source_path = os.path.join(DATA_DIR, "source.sdf")
        shutil.copy(input_arg, source_path)
    else:
        print(f"[!] Erreur : L'argument '{
              input_arg}' n'est ni une URL valide ni un fichier existant.")
        sys.exit(1)

    conversion_tasks = split_molecules(source_path)
    process_conversions(conversion_tasks)

    groups, all_mols, signatures = find_isomorphs()

    sim_matrix = compute_similarity_matrix(all_mols)

    generate_report(groups, all_mols, sim_matrix)

    save_groups_json(signatures=signatures)
