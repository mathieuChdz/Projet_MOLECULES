from scipy.cluster.hierarchy import linkage, fcluster
import plotly.graph_objects as go
from Distance import vecteur_fingerprint2D, genere_shape, get_similarities
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
import json


# Configuration
TEMPLATE_FILE = "template.html"
EXEC_CHECK_ISO = "./Check_iso"


def setup(force_clean=False):
    """Prépare l'arborescence. Ne supprime les données que si force_clean est True."""
    if force_clean and os.path.exists(DATA_DIR):
        print(f"[*] Nettoyage complet du dossier {DATA_DIR}...")
        try:
            shutil.rmtree(DATA_DIR)
        except OSError as e:
            print(f"[!] Erreur lors du nettoyage : {e}")

    for folder in [DATA_DIR, MOL_DIR, GRAPH_DIR, IMG_DIR, GROUPEMENT_DIR, CACHE_DIR]:
        os.makedirs(folder, exist_ok=True)


def load_cache():
    """Tente de charger les données intermédiaires depuis le cache (universel)."""
    cache_mols = os.path.join(CACHE_DIR, "all_mols.json")
    cache_groups = os.path.join(CACHE_DIR, "groups.json")
    cache_signatures = os.path.join(CACHE_DIR, "signatures.json")
    cache_matrix = os.path.join(CACHE_DIR, "sim_matrix_full.json")

    if os.path.exists(cache_mols) and os.path.exists(cache_groups) and os.path.exists(cache_signatures):
        with open(cache_mols, "r") as f:
            all_mols = json.load(f)
        with open(cache_groups, "r") as f:
            groups = json.load(f)
        with open(cache_signatures, "r") as f:
            signatures = json.load(f)

        sim_matrix = None
        if os.path.exists(cache_matrix):
            with open(cache_matrix, "r") as f:
                sim_matrix = json.load(f)

        return all_mols, groups, signatures, sim_matrix

    return None, None, None, None


def save_cache_base(all_mols, groups, signatures):
    """Sauvegarde les résultats de base."""
    with open(os.path.join(CACHE_DIR, "all_mols.json"), "w") as f:
        json.dump(all_mols, f)
    with open(os.path.join(CACHE_DIR, "groups.json"), "w") as f:
        json.dump(groups, f)
    with open(os.path.join(CACHE_DIR, "signatures.json"), "w") as f:
        json.dump(signatures, f)


def save_cache_matrix(sim_matrix):
    """Sauvegarde la matrice de similarité (contenant 2D et 3D)."""
    with open(os.path.join(CACHE_DIR, "sim_matrix_full.json"), "w") as f:
        json.dump(sim_matrix, f)


def download_data(url):
    print(f"[*] Téléchargement depuis : {url}")
    path = os.path.join(DATA_DIR, "source.sdf")
    try:
        r = requests.get(url, stream=True)
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))

        with open(path, 'wb') as f, tqdm(
            desc="Téléchargement", total=total_size, unit='B', unit_scale=True, unit_divisor=1024
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
        chebi_name = mol.GetProp(
            "ChEBI NAME") if mol.HasProp("ChEBI NAME") else ""
        pubchem_cid = mol.GetProp("PUBCHEM_COMPOUND_CID") if mol.HasProp(
            "PUBCHEM_COMPOUND_CID") else ""
        pubchem_iupac = mol.GetProp("PUBCHEM_IUPAC_NAME") if mol.HasProp(
            "PUBCHEM_IUPAC_NAME") else ""
        sdf_name = mol.GetProp('_Name') if mol.HasProp('_Name') else ""

        id_part = safe_filename(chebi_id) or safe_filename(
            pubchem_cid) or safe_filename(sdf_name) or "UNKNOWN_ID"
        name_part = safe_filename(chebi_name) or safe_filename(
            pubchem_iupac) or safe_filename(sdf_name) or "UNKNOWN_NAME"

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
    os.makedirs(GROUPEMENT_DIR, exist_ok=True)
    output_path = os.path.join(GROUPEMENT_DIR, "groupes_isomorphes.json")
    multi_groups = [sorted(mols)
                    for mols in signatures.values() if len(mols) > 1]
    multi_groups.sort(key=len, reverse=True)
    singletons = sorted([mols[0]
                        for mols in signatures.values() if len(mols) == 1])

    payload = {
        "date": time.strftime("%d/%m/%Y %H:%M:%S"),
        "molecules_analysees": sum(len(v) for v in signatures.values()),
        "nombre_total_groupes": len(signatures),
        "groupes_isomorphes": [{"id": gid, "taille": len(members), "molecules": members} for gid, members in enumerate(multi_groups, 1)],
        "singletons": singletons,
    }
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def compute_similarity_matrix(all_molecules):
    n = len(all_molecules)
    print(f"[*] Calcul Similarité (2D + 3D) pour {n} molécules...")

    precomputed_data = []
    for name in tqdm(all_molecules, desc="Génération 3D & FP", unit="mol"):
        path = os.path.join(MOL_DIR, f"{name}.mol")
        vec = vecteur_fingerprint2D(path)
        # On génère TOUJOURS la 3D pour la matrice universelle
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
                    matrix[name_i][name_j] = {"sim2d": 1.0, "sim3d": 1.0}
                elif j < i:
                    matrix[name_i][name_j] = matrix[name_j][name_i]
                else:
                    data_j = precomputed_data[j]
                    try:
                        val = get_similarities(data_i, data_j)
                        matrix[name_i][name_j] = val
                    except Exception:
                        matrix[name_i][name_j] = {"sim2d": 0.0, "sim3d": 0.0}
                    pbar.update(1)

    return matrix


def generate_clustering(all_molecules, sim_matrix, alpha, num_clusters=23):
    n = len(all_molecules)

    if num_clusters > n:
        print(f"[!] Attention : Seulement {
              n} molécules trouvées. Réduction du nombre de clusters à {n}.")
        num_clusters = n

    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            scores = sim_matrix[all_molecules[i]][all_molecules[j]]
            sim_combine = alpha * scores["sim2d"] + \
                (1 - alpha) * scores["sim3d"]
            D[i, j] = 1.0 - sim_combine

    np.fill_diagonal(D, 0)
    D = (D + D.T) / 2.0
    condensed = squareform(D)
    Z = linkage(condensed, method='average')
    seuil_couleur = Z[-num_clusters + 1, 2] - 1e-5

    fig = ff.create_dendrogram(
        np.zeros((n, 1)), labels=all_molecules, linkagefun=lambda x: Z, color_threshold=seuil_couleur
    )

    tickvals = fig.layout.xaxis.tickvals
    ticktext = fig.layout.xaxis.ticktext

    fig.add_trace(go.Scatter(
        x=tickvals, y=[0] * len(tickvals), mode='markers',
        marker=dict(symbol='square', size=12, color='#007bff',
                    line=dict(color='white', width=1)),
        text=ticktext, hoverinfo='text', name='CarresCliquables'
    ))

    fig.update_layout(
        title=f"Clustering Hiérarchique ({
            num_clusters} Clusters) - Alpha={alpha}",
        xaxis=dict(showticklabels=False, ticks='', title=""),
        yaxis=dict(title="Distance moléculaire", showgrid=True, gridwidth=1,
                   gridcolor='LightGray', zeroline=True, zerolinewidth=2, zerolinecolor='black'),
        width=None, height=700, hovermode='closest', plot_bgcolor='white', showlegend=False
    )

    dendrogram_html = fig.to_html(
        full_html=False, include_plotlyjs='cdn', div_id="plotly_dendrogram")
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


def get_combined_matrix(sim_matrix, alpha, all_molecules):
    """
    Combine les scores 2D et 3D en une matrice simple (flottants) 
    pour éviter de casser le template HTML Jinja.
    """
    combined = {}
    for m1 in all_molecules:
        combined[m1] = {}
        for m2 in all_molecules:
            scores = sim_matrix[m1][m2]
            combined[m1][m2] = alpha * scores["sim2d"] + \
                (1 - alpha) * scores["sim3d"]
    return combined


def generate_report(groups, all_molecules, combined_matrix, clusters_data, dendrogram_html):
    print("[*] Génération du rapport HTML via Jinja2...")
    now = time.strftime("%d/%m/%Y %H:%M:%S")
    file_loader = FileSystemLoader('.')
    env = Environment(loader=file_loader)

    try:
        template = env.get_template(TEMPLATE_FILE)
    except Exception as e:
        print(f"[!] Erreur : Impossible de trouver '{
              TEMPLATE_FILE}'. Vérifie qu'il est dans le dossier.")
        sys.exit(1)

    output = template.render(
        now=now, groups=groups, all_molecules=all_molecules,
        matrix=combined_matrix, clusters=clusters_data, dendrogram_html=dendrogram_html
    )

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write(output)

    print(f"\n[*] Rapport généré : {os.path.abspath(REPORT_FILE)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyse structurelle et isomorphisme de molécules à partir d'un fichier SDF.")
    parser.add_argument(
        "input", type=str, help="Chemin vers un fichier SDF local OU lien de téléchargement (http/https).")
    parser.add_argument("-hc", "--hierarchique_cluster", type=int,
                        help="Nombre de clusters générés par le clustering hiérarchique.")
    parser.add_argument("-a", "--alpha", type=float, default=0.5,
                        help="Ratio (entre 0 et 1) pour le calcul du score (ex: poids du Tanimoto 2D vs 3D).")
    parser.add_argument("-o", "--output", type=str, default="data",
                        help="Nom du dossier de sortie pour les résultats (par défaut: 'data').")

    args = parser.parse_args()
    input_arg = args.input
    alpha_arg = args.alpha
    nbrcluster = args.hierarchique_cluster

    if not (0.0 <= alpha_arg <= 1.0):
        print("[!] Erreur : La valeur de l'argument -a doit être comprise entre 0 et 1.")
        sys.exit(1)

    global DATA_DIR, MOL_DIR, GRAPH_DIR, IMG_DIR, GROUPEMENT_DIR, REPORT_FILE, CACHE_DIR
    DATA_DIR = os.path.join(args.output, "data")
    MOL_DIR = os.path.join(DATA_DIR, "molecules")
    GRAPH_DIR = os.path.join(DATA_DIR, "graphs")
    IMG_DIR = os.path.join(DATA_DIR, "images")
    GROUPEMENT_DIR = os.path.join(DATA_DIR, "groupement")
    CACHE_DIR = os.path.join(DATA_DIR, "cache")
    REPORT_FILE = os.path.join(args.output, "resultats.html")

    all_mols, groups, signatures, sim_matrix = load_cache()

    if all_mols is not None:
        print(f"[*] Cache de base trouvé dans '{
              DATA_DIR}'. Réutilisation de l'extraction SDF et de l'analyse Nauty...")
        setup(force_clean=False)

        if sim_matrix is None:
            print(
                "[*] Matrice de similarité universelle introuvable en cache. Calcul de la matrice (2D et 3D)...")
            sim_matrix = compute_similarity_matrix(all_mols)
            save_cache_matrix(sim_matrix)
        else:
            print(
                "[*] Matrice de similarité complète (2D/3D) trouvée en cache ! Calcul instantané.")

    else:
        print(f"[*] Aucun cache trouvé dans '{DATA_DIR}'. Exécution totale...")
        setup(force_clean=True)

        if input_arg.startswith("http://") or input_arg.startswith("https://"):
            source_path = download_data(input_arg)
        elif os.path.isfile(input_arg):
            print(f"[*] Utilisation du fichier local : {input_arg}")
            source_path = os.path.join(DATA_DIR, "source.sdf")
            shutil.copy(input_arg, source_path)
        else:
            print(f"[!] Erreur : L'argument '{
                  input_arg}' n'est ni une URL valide ni un fichier existant.")
            sys.exit(1)

        conversion_tasks = split_molecules(source_path)
        process_conversions(conversion_tasks)

        groups, all_mols, signatures = find_isomorphs()
        save_cache_base(all_mols, groups, signatures)

        sim_matrix = compute_similarity_matrix(all_mols)
        save_cache_matrix(sim_matrix)

    clusters_dict, dendro_html = generate_clustering(
        all_mols, sim_matrix, alpha=alpha_arg, num_clusters=nbrcluster)

    # On prépare une matrice simple avec l'alpha sélectionné pour ne pas casser ton Template HTML
    matrice_html = get_combined_matrix(sim_matrix, alpha_arg, all_mols)

    save_groups_json(signatures=signatures)
    generate_report(groups, all_mols, matrice_html, clusters_dict, dendro_html)
