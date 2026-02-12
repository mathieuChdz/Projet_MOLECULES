import os
import subprocess
import sys
import requests
import time
from rdkit import Chem
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem import Draw
import re

import tanimoto_3d

NAUTY_DIR = "./nauty"
DATA_DIR = "data"
MOL_DIR = f"{DATA_DIR}/molecules"
GRAPH_DIR = f"{DATA_DIR}/graphs"
REPORT_FILE = "resultats.html"
IMG_DIR = f"{DATA_DIR}/images"

def setup():
    os.makedirs(MOL_DIR, exist_ok=True)
    os.makedirs(GRAPH_DIR, exist_ok=True)
    os.makedirs(IMG_DIR, exist_ok=True)
    for folder in [MOL_DIR, GRAPH_DIR, IMG_DIR]:
        for f in os.listdir(folder):
            os.remove(os.path.join(folder, f))

def download_data(url):
    path = os.path.join(DATA_DIR, "source.sdf")
    try:
        r = requests.get(url)
        r.raise_for_status()
        with open(path, 'wb') as f:
            f.write(r.content)
        return path
    except Exception as e:
        sys.exit(1)

def safe_filename(s):
    s = s.strip()
    s = s.replace(" ", "_")
    s = re.sub(r"[^a-zA-Z0-9_\-]", "", s)
    return s[:80]

def split_molecules(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path)
    
    for i, mol in enumerate(suppl):
        if mol is None:
            continue

        # Nom générique
        base_name = mol.GetProp('_Name') if mol.HasProp('_Name') else "MOL"

        # Nom IUPAC PubChem
        iupac = ""
        if mol.HasProp("PUBCHEM_IUPAC_NAME"):
            iupac = mol.GetProp("PUBCHEM_IUPAC_NAME")

        # Nettoyage pour fichiers
        base_name = safe_filename(base_name)
        iupac = safe_filename(iupac)

        # Nom final
        if iupac:
            out_name = f"mol_{i}_{base_name}_{iupac}"
        else:
            out_name = f"mol_{i}_{base_name}"

        mol_path = os.path.join(MOL_DIR, f"{out_name}.mol")
        with open(mol_path, "w") as f:
            f.write(Chem.MolToMolBlock(mol))

        img_path = os.path.join(IMG_DIR, f"{out_name}.png")
        Draw.MolToFile(mol, img_path, size=(300, 300))


def process_conversions():
    tasks = []
    with ProcessPoolExecutor() as executor:
        for f in os.listdir(MOL_DIR):
            if f.endswith(".mol"):
                input_p = os.path.join(MOL_DIR, f)
                output_p = os.path.join(GRAPH_DIR, f.replace(".mol", ".graph"))
                tasks.append(executor.submit(subprocess.run, 
                             ["python", "process_mol.py", input_p, output_p]))
    for t in tasks: t.result()

def compile_engine():
    if not os.path.exists("./check_iso"):
        cmd = ["gcc", "-O3", "-o", "check_iso", "check_iso.c", 
               f"{NAUTY_DIR}/nauty.a", f"-I{NAUTY_DIR}"]
        subprocess.run(cmd, check=True)

def find_isomorphs():
    signatures = {}
    all_molecules = []

    graphs = [f for f in os.listdir(GRAPH_DIR) if f.endswith(".graph")]

    for g in graphs:
        mol_name = g.replace(".graph", "")
        all_molecules.append(mol_name)

        path = os.path.join(GRAPH_DIR, g)
        res = subprocess.run(["./check_iso", path], capture_output=True, text=True)
        sig = res.stdout.strip()

        if sig not in signatures:
            signatures[sig] = []
        signatures[sig].append(mol_name)

    groups = [mols for mols in signatures.values() if len(mols) > 1]
    return groups, sorted(all_molecules)


def generate_report(groups, all_molecules, similarity_scores):
    now = time.strftime("%d/%m/%Y %H:%M:%S")

    html = f"""
    <html>
    <head>
        <meta charset="utf-8">
        <title>Rapport d'Analyse Moléculaire</title>
        <style>
            body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 1100px; margin: auto; background: #f4f7f6; padding: 20px; }}
            h1 {{ color: #2c3e50; }}
            h2 {{ color: #34495e; margin-top: 40px; }}
            .summary {{ background: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; }}
            .mol-container {{ display: flex; flex-wrap: wrap; gap: 20px; }}
            .mol-item {{ width: 220px; background: white; border-radius: 8px; padding: 10px; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,.1); }}
            .mol-item img {{ width: 200px; height: 200px; }}
            .mol-link {{ display: block; margin-top: 5px; font-size: 0.8em; color: #2980b9; word-break: break-word; }}
            .group-card {{ background: white; padding: 20px; margin-top: 20px; border-left: 6px solid #3498db; border-radius: 8px; }}
            
            /* Style pour le nouveau tableau Tanimoto */
            .sim-table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,.1); margin-top: 20px; }}
            .sim-table th {{ background: #34495e; color: white; padding: 12px; text-align: left; }}
            .sim-table td {{ padding: 12px; border-bottom: 1px solid #eee; }}
            .sim-table tr:last-child td {{ border-bottom: none; }}
            .score-badge {{ font-weight: bold; padding: 4px 8px; border-radius: 4px; }}
        </style>
    </head>
    <body>
        <h1>Rapport d'analyse moléculaire</h1>

        <div class="summary">
            <p><b>Date :</b> {now}</p>
            <p><b>Molécules analysées :</b> {len(all_molecules)}</p>
            <p><b>Groupes d'isomorphes :</b> {len(groups)}</p>
        </div>

        <h2>Toutes les molécules analysées</h2>
        <div class="mol-container">
    """

    for mol_name in all_molecules:
        parts = mol_name.split("_")
        iupac = "_".join(parts[3:]) if len(parts) > 3 else ""
        label = iupac.replace("_", " ") if iupac else mol_name
        img_path = f"data/images/{mol_name}.png"

        html += f"""
        <div class="mol-item">
            <img src="{img_path}">
            <div class="mol-link">{label}</div>
        </div>
        """

    html += "</div>"

    if groups:
        html += "<h2>Groupes de structures isomorphes (Identité exacte)</h2>"
        for i, grp in enumerate(groups):
            html += f'<div class="group-card"><b>Groupe {i+1}</b><div class="mol-container">'
            for mol_name in grp:
                parts = mol_name.split("_")
                iupac = "_".join(parts[3:]) if len(parts) > 3 else ""
                label = iupac.replace("_", " ") if iupac else mol_name
                img_path = f"data/images/{mol_name}.png"
                pubchem = f"https://pubchem.ncbi.nlm.nih.gov/#query={label}"
                html += f"""
                <div class="mol-item">
                    <img src="{img_path}">
                    <a class="mol-link" href="{pubchem}" target="_blank">{label}</a>
                </div>
                """
            html += "</div></div>"

    # NOUVELLE SECTION : Similarité Tanimoto 3D
    if similarity_scores:
        html += "<h2>Similarité Spatiale (Tanimoto 3D)</h2>"
        html += """
        <table class="sim-table">
            <thead>
                <tr>
                    <th>Molécule A</th>
                    <th>Molécule B</th>
                    <th>Score de proximité</th>
                </tr>
            </thead>
            <tbody>
        """
        for m1, m2, score in similarity_scores:

            n1 = m1.split("_")[-1].replace("-", " ")
            n2 = m2.split("_")[-1].replace("-", " ")
            
            # Couleurs pour affichage du score
            color = "#27ae60" if score > 0.7 else "#f39c12" if score > 0.4 else "#c0392b"
            
            html += f"""
                <tr>
                    <td>{n1}</td>
                    <td>{n2}</td>
                    <td><span class="score-badge" style="color: {color};">{score:.4f}</span></td>
                </tr>
            """
        html += "</tbody></table>"

    html += """
        <br><br>
        </body>
        </html>
    """

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write(html)


def run_tanimoto_analysis(all_molecules):
    """Calcule les scores pour toutes les paires uniques de molécules."""
    scores = []
    # On récupère les objets RDKit depuis les fichiers .mol générés
    for i in range(len(all_molecules)):
        for j in range(i + 1, len(all_molecules)):
            m1_name = all_molecules[i]
            m2_name = all_molecules[j]
            
            mol1 = Chem.MolFromMolFile(os.path.join(MOL_DIR, f"{m1_name}.mol"))
            mol2 = Chem.MolFromMolFile(os.path.join(MOL_DIR, f"{m2_name}.mol"))
            
            score = tanimoto_3d.calculate_similarity(mol1, mol2)
            scores.append((m1_name, m2_name, score))
    return scores

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: make run URL=\"...\"")
        sys.exit(1)
    
    setup()
    source = download_data(sys.argv[1])
    split_molecules(source)
    process_conversions()
    
    compile_engine()
    groups, all_molecules = find_isomorphs()
    
    print("Calcul de la similarité Tanimoto 3D...")
    similarity_scores = []
    
    for i in range(len(all_molecules)):
        for j in range(i + 1, len(all_molecules)):
            m1_name = all_molecules[i]
            m2_name = all_molecules[j]
            
            mol_a = Chem.MolFromMolFile(os.path.join(MOL_DIR, f"{m1_name}.mol"))
            mol_b = Chem.MolFromMolFile(os.path.join(MOL_DIR, f"{m2_name}.mol"))
            
            if mol_a and mol_b:
                score = tanimoto_3d.calcul_similarite(mol_a, mol_b)
                similarity_scores.append((m1_name, m2_name, score))

    print(f"Génération du rapport : {REPORT_FILE}")
    generate_report(groups, all_molecules, similarity_scores)