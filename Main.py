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

from Distance import vecteur_fingerprint2D, score, genere_shape

# Configuration
NAUTY_DIR = "nauty2_9_3"
DATA_DIR = "data"
MOL_DIR = os.path.join(DATA_DIR, "molecules")
GRAPH_DIR = os.path.join(DATA_DIR, "graphs")
IMG_DIR = os.path.join(DATA_DIR, "images")
REPORT_FILE = "resultats.html"
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
        # On peut ajouter un stream pour une barre de téléchargement,
        # mais pour un SDF texte c'est souvent rapide.
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

    except Exception as e:
        # On évite les prints dans les threads pour ne pas casser la barre TQDM
        pass


def split_molecules(sdf_path):
    print("[*] Découpage des molécules...")
    suppl = Chem.SDMolSupplier(sdf_path)
    conversion_args = []

    # On ne connaît pas la longueur exacte du supplier sans le lire,
    # donc on utilise tqdm sans 'total' ou on fait une estimation.
    # Ici une boucle simple avec tqdm.

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
    """Lance la conversion mol -> graph en parallèle avec barre de progression."""
    total = len(conversion_args)
    print(f"[*] Conversion de {total} molécules en graphes...")

    with ProcessPoolExecutor() as executor:
        # On soumet toutes les tâches
        futures = [executor.submit(convert_mol_to_graph, *args)
                   for args in conversion_args]

        # On utilise as_completed pour mettre à jour la barre au fur et à mesure
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

    # Ajout TQDM ici
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
        except Exception as e:
            pass  # On évite de spammer

    groups = [mols for mols in signatures.values() if len(mols) > 1]
    groups.sort(key=len, reverse=True)
    return groups, sorted(all_molecules)


def compute_similarity_matrix(all_molecules):
    """
    Calcule la matrice de similarité avec cache et objets pré-calculés.
    """
    n = len(all_molecules)
    print(f"[*] Calcul Similarité pour {n} molécules...")

    precomputed_data = []

    # Barre 1 : Pré-calcul (Génération 3D)
    for name in tqdm(all_molecules, desc="Génération 3D & FP", unit="mol"):
        path = os.path.join(MOL_DIR, f"{name}.mol")

        # 1. Fingerprint
        vec = vecteur_fingerprint2D(path)

        # 2. Objet 3D (renvoie (mol, ids))
        mol_3d, ids = genere_shape(path)

        precomputed_data.append((vec, mol_3d, ids))

    print(f"\n[*] Calcul de la matrice croisée...")
    matrix = {}
    total_ops = (n * (n + 1)) // 2

    # Barre 2 : Calcul Matrice
    # On utilise un context manager pour la barre manuelle
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

                    # Mise à jour de la barre d'un pas
                    pbar.update(1)

    return matrix


def generate_report(groups, all_molecules, similarity_matrix):
    print("[*] Génération du rapport HTML...")
    now = time.strftime("%d/%m/%Y %H:%M:%S")

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Rapport d'Isomorphisme & Distances</title>
        <style>
            body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 1400px; margin: auto; background: #f4f7f6; padding: 20px; color: #333; }}
            h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
            h2 {{ color: #34495e; margin-top: 40px; border-left: 5px solid #27ae60; padding-left: 10px; }}
            .summary {{ background: white; padding: 20px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }}
            .mol-container {{ display: flex; flex-wrap: wrap; gap: 15px; justify-content: center; }}
            .mol-item {{ width: 200px; background: white; border-radius: 8px; padding: 10px; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            .mol-item img {{ width: 180px; height: 180px; object-fit: contain; }}
            .mol-link {{ display: block; margin-top: 8px; font-size: 0.8em; color: #2980b9; word-break: break-all; }}
            .group-card {{ background: white; padding: 20px; margin-top: 20px; border-left: 5px solid #27ae60; border-radius: 5px; }}

            .matrix-wrapper {{ overflow-x: auto; background: white; padding: 15px; border-radius: 8px; margin-top: 20px; }}
            table.matrix {{ border-collapse: collapse; width: 100%; font-size: 0.85em; }}
            table.matrix th, table.matrix td {{ padding: 8px; text-align: center; border: 1px solid #dfe6e9; }}
            table.matrix th {{ background-color: #f8f9fa; min-width: 100px; }}
            th.rotate {{ height: auto; white-space: nowrap; padding: 10px; vertical-align: bottom; }}
            th.rotate > div {{ transform: none; width: auto; }}

            .d-high {{ background-color: #27ae60; color: white; }}
            .d-med {{ background-color: #7bed9f; color: black; }}
            .d-low {{ background-color: #dfe6e9; color: #b2bec3; }}
        </style>
    </head>
    <body>
        <h1>Rapport d'analyse moléculaire</h1>
        <div class="summary">
            <p><strong>Date :</strong> {now}</p>
            <p><strong>Molécules :</strong> {len(all_molecules)}</p>
            <p><strong>Groupes d'isomorphes :</strong> {len(groups)}</p>
        </div>
    """

    # Matrice
    html += """<h2>📊 Matrice de Similarité</h2><div class="matrix-wrapper"><table class="matrix"><thead><tr><th></th>"""
    for mol in all_molecules:
        parts = mol.split("_")
        short = parts[-1] if len(parts) > 1 else mol
        html += f'<th class="rotate"><div>{short[:15]}...</div></th>'
    html += "</tr></thead><tbody>"

    for m1 in all_molecules:
        parts = m1.split("_")
        row = parts[-1] if len(parts) > 1 else m1
        html += f"<tr><th>{row[:20]}</th>"
        for m2 in all_molecules:
            val = similarity_matrix[m1][m2]
            cls = "d-high" if val >= 0.95 else "d-med" if val >= 0.7 else "d-low" if val < 0.3 else ""
            html += f'<td class="{cls}">{val:.2f}</td>'
        html += "</tr>"
    html += "</tbody></table></div>"

    # Isomorphes
    if groups:
        html += "<h2>🔍 Groupes d'isomorphes</h2>"
        for i, grp in enumerate(groups):
            html += f'<div class="group-card"><b>Groupe {
                i+1}</b><div class="mol-container">'
            for m in grp:
                img = f"data/images/{m}.png"
                label = " ".join(m.split("_")[3:]) if len(
                    m.split("_")) > 3 else m
                html += f'<div class="mol-item"><img src="{
                    img}"><div class="mol-link">{label}</div></div>'
            html += "</div></div>"

    # Toutes les molécules
    html += "<h2>📂 Toutes les molécules</h2><div class='mol-container'>"
    for m in all_molecules:
        img = f"data/images/{m}.png"
        html += f'<div class="mol-item"><img src="{
            img}"><div class="mol-link">{m}</div></div>'
    html += "</div></body></html>"

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\n[*] Rapport généré : {os.path.abspath(REPORT_FILE)}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python main.py <URL_SDF>")
        sys.exit(1)

    setup()
    source_path = download_data(sys.argv[1])
    conversion_tasks = split_molecules(source_path)
    process_conversions(conversion_tasks)
    groups, all_mols = find_isomorphs()
    sim_matrix = compute_similarity_matrix(all_mols)
    generate_report(groups, all_mols, sim_matrix)
