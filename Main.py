import os
import sys
import requests
import time
import subprocess
import re
from concurrent.futures import ProcessPoolExecutor
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

import shutil
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
            # Supprime tout : dossiers, fichiers, sous-dossiers
            shutil.rmtree(DATA_DIR)
        except OSError as e:
            print(f"[!] Erreur lors du nettoyage : {e}")
            # On continue même si erreur (ex: fichier ouvert), les makedirs feront le reste

    # Recréation des dossiers vides
    for folder in [DATA_DIR, MOL_DIR, GRAPH_DIR, IMG_DIR]:
        os.makedirs(folder, exist_ok=True)


def download_data(url):
    """Télécharge le fichier SDF source."""
    print(f"[*] Téléchargement depuis : {url}")
    path = os.path.join(DATA_DIR, "source.sdf")
    try:
        r = requests.get(url)
        r.raise_for_status()
        with open(path, 'wb') as f:
            f.write(r.content)
        return path
    except Exception as e:
        print(f"[!] Erreur de téléchargement : {e}")
        sys.exit(1)


def safe_filename(s):
    """Nettoie une chaîne pour l'utiliser comme nom de fichier."""
    s = str(s).strip().replace(" ", "_")
    s = re.sub(r"[^a-zA-Z0-9_\-]", "", s)
    return s[:80]


def convert_mol_to_graph(mol_path, output_path):
    """Convertit un fichier .mol en fichier .graph (Format Adjacence)."""
    try:
        mol = Chem.MolFromMolFile(mol_path)
        if not mol:
            return

        num_atoms = mol.GetNumAtoms()
        edges = []
        for bond in mol.GetBonds():
            u = bond.GetBeginAtomIdx()
            v = bond.GetEndAtomIdx()
            edges.append(f"{u} {v}")

        with open(output_path, "w") as f:
            f.write(f"{num_atoms}\n")
            f.write("\n".join(edges))
            f.write("\n")

    except Exception as e:
        print(f"[!] Erreur conversion {mol_path}: {e}")


def split_molecules(sdf_path):
    """Découpe le SDF en fichiers .mol individuels et génère les images."""
    print("[*] Découpage des molécules...")
    suppl = Chem.SDMolSupplier(sdf_path)
    conversion_args = []

    for i, mol in enumerate(suppl):
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
    """Lance la conversion mol -> graph en parallèle."""
    print(f"[*] Conversion de {len(conversion_args)} molécules en graphes...")
    with ProcessPoolExecutor() as executor:
        for args in conversion_args:
            executor.submit(convert_mol_to_graph, *args)


def find_isomorphs():
    """Appelle le programme C check_iso pour chaque graphe."""
    print("[*] Recherche d'isomorphes (Nauty)...")

    if not os.path.exists(EXEC_CHECK_ISO):
        print(f"[!] Erreur : L'exécutable {EXEC_CHECK_ISO} est introuvable.")
        print("    Lancez 'make' pour le compiler.")
        sys.exit(1)

    signatures = {}
    all_molecules = []

    graphs = [f for f in os.listdir(GRAPH_DIR) if f.endswith(".graph")]

    for i, g in enumerate(graphs):
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
            print(f"[!] Erreur sur {g}: {e}")

    groups = [mols for mols in signatures.values() if len(mols) > 1]
    groups.sort(key=len, reverse=True)
    return groups, sorted(all_molecules)


def compute_similarity_matrix(all_molecules):
    """
    Calcule la matrice de similarité en utilisant Distance.score().
    Optimisé pour ne calculer que le triangle supérieur.
    """
    print(
        f"[*] Calcul des distances (Similarité 2D+3D) pour {len(all_molecules)} molécules...")
    matrix = {}
    n = len(all_molecules)
    matrice_vecteur = [(vecteur_fingerprint2D(os.path.join(MOL_DIR, f"{
        all_molecules[i]}.mol")), genere_shape(os.path.join(MOL_DIR, f"{all_molecules[i]}.mol"))) for i in range(n)]

    # Barre de progression simple
    total_ops = (n * (n + 1)) // 2
    current_op = 0

    for i in range(n):
        name_i = all_molecules[i]
        shape_i = matrice_vecteur[i][1]
        vecteur_i = matrice_vecteur[i][0]
        matrix[name_i] = {}

        for j in range(n):
            name_j = all_molecules[j]
            shape_j = matrice_vecteur[j][1]
            vecteur_j = matrice_vecteur[j][0]
            if i == j:
                matrix[name_i][name_j] = 1.0  # Identique à soi-même
            elif j < i:
                # Récupérer la valeur symétrique (déjà calculée)
                matrix[name_i][name_j] = matrix[name_j][name_i]
            else:
                # Calcul réel via Distance.py
                try:
                    # On appelle la fonction score importée
                    val = score(vecteur_i, shape_i,
                                vecteur_j, shape_j)
                    matrix[name_i][name_j] = val
                except Exception as e:
                    print(f"Err {name_i}/{name_j}: {e}")
                    matrix[name_i][name_j] = 0.0

                current_op += 1
                if current_op % 5 == 0:
                    print(f"\r    Progression : {
                          current_op}/{total_ops}", end="")

    print("\n[*] Calcul terminé.")
    return matrix


def generate_report(groups, all_molecules, similarity_matrix):
    """Génère le rapport HTML incluant la matrice de distances."""
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

            /* Container Molécules */
            .mol-container {{ display: flex; flex-wrap: wrap; gap: 15px; justify-content: center; }}
            .mol-item {{ width: 200px; background: white; border-radius: 8px; padding: 10px; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.1); transition: transform 0.2s; }}
            .mol-item:hover {{ transform: translateY(-5px); box-shadow: 0 5px 15px rgba(0,0,0,0.1); }}
            .mol-item img {{ width: 180px; height: 180px; object-fit: contain; }}
            .mol-link {{ display: block; margin-top: 8px; font-size: 0.85em; color: #2980b9; text-decoration: none; word-break: break-all; }}

            .group-card {{ background: white; padding: 20px; margin-top: 20px; border-left: 5px solid #27ae60; border-radius: 5px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }}

            /* Matrice de Distance */
            .matrix-wrapper {{ overflow-x: auto; background: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }}
            table.matrix {{ border-collapse: collapse; width: 100%; font-size: 0.85em; }}
            table.matrix th, table.matrix td {{ padding: 8px; text-align: center; border: 1px solid #dfe6e9; }}
            table.matrix th {{ background-color: #f8f9fa; font-weight: 600; min-width: 100px; }}

            /* Headers horizontaux */
            th.rotate {{ height: auto; white-space: nowrap; padding: 10px; vertical-align: bottom; }}
            th.rotate > div {{ transform: none; width: auto; }}

            /* Couleurs Heatmap */
            .d-high {{ background-color: #27ae60; color: white; }} /* > 0.9 */
            .d-med {{ background-color: #7bed9f; color: black; }} /* > 0.7 */
            .d-low {{ background-color: #dfe6e9; color: #b2bec3; }} /* < 0.5 */
        </style>
    </head>
    <body>
        <h1>Rapport d'analyse moléculaire</h1>

        <div class="summary">
            <p><strong>Date :</strong> {now}</p>
            <p><strong>Molécules analysées :</strong> {len(all_molecules)}</p>
            <p><strong>Groupes d'isomorphes (Graphes) :</strong> {len(groups)}</p>
        </div>
    """

    # --- 1. Matrice de Similarité ---
    html += """
        <h2>📊 Matrice de Similarité (Distance)</h2>
        <p><i>Score combiné (0.0 = différent, 1.0 = identique). Basé sur la structure 2D et la forme 3D.</i></p>
        <div class="matrix-wrapper">
        <table class="matrix">
            <thead>
                <tr>
                    <th></th>
    """
    # En-têtes colonnes
    for mol in all_molecules:
        parts = mol.split("_")
        short_name = parts[-1] if len(parts) > 1 else mol
        html += f'<th class="rotate"><div>{short_name[:15]}...</div></th>'

    html += "</tr></thead><tbody>"

    # Lignes
    for m1 in all_molecules:
        parts = m1.split("_")
        row_name = parts[-1] if len(parts) > 1 else m1
        html += f"<tr><th>{row_name[:20]}</th>"

        for m2 in all_molecules:
            score_val = similarity_matrix[m1][m2]

            # Choix de la classe CSS pour la couleur
            css_class = ""
            if score_val >= 0.95:
                css_class = "d-high"
            elif score_val >= 0.70:
                css_class = "d-med"
            elif score_val < 0.3:
                css_class = "d-low"

            html += f'<td class="{css_class}">{score_val:.2f}</td>'
        html += "</tr>"

    html += "</tbody></table></div>"

    # --- 2. Isomorphes ---
    if groups:
        html += "<h2>🔍 Groupes d'isomorphes (Graphes identiques)</h2>"
        for i, grp in enumerate(groups):
            html += f'<div class="group-card"><b>Groupe {
                i+1}</b> ({len(grp)} molécules)<div class="mol-container">'
            for mol_name in grp:
                img_path = f"data/images/{mol_name}.png"
                parts = mol_name.split("_")
                label = " ".join(parts[3:]) if len(parts) > 3 else mol_name
                html += f"""
                <div class="mol-item">
                    <img src="{img_path}" loading="lazy">
                    <div class="mol-link">{label}</div>
                </div>"""
            html += "</div></div>"
    else:
        html += "<h2>🔍 Isomorphes</h2><p>Aucun groupe de graphes identiques trouvé.</p>"

    # --- 3. Toutes les molécules ---
    html += "<h2>📂 Toutes les molécules</h2><div class='mol-container'>"
    for mol_name in all_molecules:
        img_path = f"data/images/{mol_name}.png"
        html += f"""
        <div class="mol-item">
            <img src="{img_path}" loading="lazy">
            <div class="mol-link">{mol_name}</div>
        </div>
        """
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

    # 1. Traitement initial
    conversion_tasks = split_molecules(source_path)
    process_conversions(conversion_tasks)

    # 2. Analyse Isomorphisme (C)
    groups, all_mols = find_isomorphs()

    # 3. Calcul des Distances (Python - Distance.py)
    sim_matrix = compute_similarity_matrix(all_mols)

    # 4. Rapport Final
    generate_report(groups, all_mols, sim_matrix)
