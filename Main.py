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
from jinja2 import Environment, FileSystemLoader  # <--- NOUVEL IMPORT

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


def setup():
    """Supprime tout le dossier data et recrée l'arborescence à neuf."""
    if os.path.exists(DATA_DIR):
        print(f"[*] Nettoyage complet du dossier {DATA_DIR}...")
        try:
            shutil.rmtree(DATA_DIR)
        except OSError as e:
            print(f"[!] Erreur lors du nettoyage : {e}")

    for folder in [DATA_DIR, MOL_DIR, GRAPH_DIR, IMG_DIR]:
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
    """Version avec types d'atomes (Couleurs) pour Nauty."""
    try:
        mol = Chem.MolFromMolFile(mol_path)
        if not mol:
            return

        num_atoms = mol.GetNumAtoms()

        atom_types = []
        for atom in mol.GetAtoms():
            atom_types.append(str(atom.GetAtomicNum()))

        edges = []
        for bond in mol.GetBonds():
            u = bond.GetBeginAtomIdx()
            v = bond.GetEndAtomIdx()
            edges.append(f"{u} {v}")

        with open(output_path, "w") as f:
            f.write(f"{num_atoms}\n")
            f.write(" ".join(atom_types) + "\n")
            f.write("\n".join(edges))
            f.write("\n")

    except Exception:
        pass


def split_molecules(sdf_path):
    print("[*] Découpage des molécules...")
    suppl = Chem.SDMolSupplier(sdf_path)
    conversion_args = []

    for i, mol in enumerate(tqdm(suppl, desc="Extraction SDF", unit="mol")):
        if mol is None:
            continue

        base_name = mol.GetProp('_Name') if mol.HasProp('_Name') else "MOL"
        iupac = mol.GetProp("PUBCHEM_IUPAC_NAME") if mol.HasProp(
            "PUBCHEM_IUPAC_NAME") else ""

        base_name = safe_filename(base_name)
        iupac = safe_filename(iupac)

        if iupac:
            out_name = f"mol_{i}_{base_name}_{iupac}"
        else:
            out_name = f"mol_{i}_{base_name}"

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
    return groups, sorted(all_molecules)


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
    # Pour que cela fonctionne, il faut copier-coller les fonctions manquantes
    # (setup, download, etc.) du message précédent ici.

    if len(sys.argv) < 2:
        print("Usage: python main.py <URL_SDF>")
        sys.exit(1)

    setup()
    source_path = download_data(sys.argv[1])
    conversion_tasks = split_molecules(source_path)
    process_conversions(conversion_tasks)
    groups, all_mols = find_isomorphs()
    sim_matrix = compute_similarity_matrix(all_mols)

    # Appel de la nouvelle fonction
    generate_report(groups, all_mols, sim_matrix)
