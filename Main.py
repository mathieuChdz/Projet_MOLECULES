from scipy.cluster.hierarchy import linkage, fcluster
import plotly.graph_objects as go  # <-- À ajouter en haut de ton script
from Distance import vecteur_fingerprint2D, score, genere_shape
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
import plotly.figure_factory as ff
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform

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


def generate_clustering(all_molecules, sim_matrix, num_clusters=23):
    n = len(all_molecules)

    if num_clusters > n:
        print(f"[!] Attention : Seulement {n} molécules trouvées. "
              f"Réduction du nombre de clusters de {num_clusters} à {n}.")
        num_clusters = n

    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = 1.0 - sim_matrix[all_molecules[i]][all_molecules[j]]

    np.fill_diagonal(D, 0)
    D = (D + D.T) / 2.0
    condensed = squareform(D)
    Z = linkage(condensed, method='average')
    seuil_couleur = Z[-num_clusters + 1, 2] - 1e-5

    fig = ff.create_dendrogram(
        np.zeros((n, 1)),
        labels=all_molecules,
        linkagefun=lambda x: Z,
        color_threshold=seuil_couleur
    )

    tickvals = fig.layout.xaxis.tickvals
    ticktext = fig.layout.xaxis.ticktext

    fig.add_trace(go.Scatter(
        x=tickvals,
        y=[0] * len(tickvals),  # On les place tous au niveau 0
        mode='markers',
        marker=dict(symbol='square', size=12, color='#007bff',
                    line=dict(color='white', width=1)),
        text=ticktext,         # Le nom apparaîtra au survol
        hoverinfo='text',
        name='CarresCliquables'  # Nom interne pour le JavaScript
    ))

    fig.update_layout(
        title=f"Clustering Hiérarchique ({num_clusters} Clusters)",
        xaxis=dict(
            showticklabels=False,  # On cache les textes illisibles
            ticks='',
            title=""
        ),
        yaxis=dict(
            title="Distance moléculaire",
            showgrid=True, gridwidth=1, gridcolor='LightGray',
            zeroline=True, zerolinewidth=2, zerolinecolor='black'
        ),
        width=None,
        height=700,
        hovermode='closest',
        plot_bgcolor='white',
        showlegend=False
    )

    dendrogram_html = fig.to_html(
        full_html=False, include_plotlyjs='cdn', div_id="plotly_dendrogram"
    )

    # ... (Le reste de ton code pour découper et trier les clusters reste identique)
    clusters_assign = fcluster(Z, num_clusters, criterion='maxclust')
    clusters_data = []
    for cluster_id in range(1, num_clusters + 1):
        members = [all_molecules[i]
                   for i in range(n) if clusters_assign[i] == cluster_id]
        if members:
            clusters_data.append(
                {'id': cluster_id, 'taille': len(members), 'molecules': members})

    sorted_clusters = sorted(
        clusters_data, key=lambda x: x['taille'], reverse=True)

    return sorted_clusters, dendrogram_html


def generate_report(groups, all_molecules, similarity_matrix, clusters_data, dendrogram_html):
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

    output = template.render(
        now=now,
        groups=groups,
        all_molecules=all_molecules,
        matrix=similarity_matrix,
        clusters=clusters_data,
        dendrogram_html=dendrogram_html  # <-- Ajoute bien ça ici
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

    groups, all_mols = find_isomorphs()

    sim_matrix = compute_similarity_matrix(all_mols)

    clusters_dict, dendro_html = generate_clustering(
        all_mols, sim_matrix, num_clusters=23)

    generate_report(groups, all_mols, sim_matrix, clusters_dict, dendro_html)
