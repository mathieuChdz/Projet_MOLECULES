import os
import subprocess
import sys
import requests
import time
from rdkit import Chem
from concurrent.futures import ProcessPoolExecutor

NAUTY_DIR = "./nauty"
DATA_DIR = "data"
MOL_DIR = f"{DATA_DIR}/molecules"
GRAPH_DIR = f"{DATA_DIR}/graphs"
REPORT_FILE = "resultats.html"

def setup():
    os.makedirs(MOL_DIR, exist_ok=True)
    os.makedirs(GRAPH_DIR, exist_ok=True)
    for folder in [MOL_DIR, GRAPH_DIR]:
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

def split_molecules(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path)
    for i, mol in enumerate(suppl):
        if mol:
            name = mol.GetProp('_Name') if mol.HasProp('_Name') else "MOL"
            chebi_id = mol.GetProp('ChEBI ID') if mol.HasProp('ChEBI ID') else ""

            # On met i au début pour être sûr que mol_0_cafeine et mol_1_cafeine soient deux fichiers différents
            out_name = f"mol_{i}_{name}_{chebi_id}".strip("_")
            
            path = os.path.join(MOL_DIR, f"{out_name}.mol")
            with open(path, 'w') as f:
                f.write(Chem.MolToMolBlock(mol))

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
    graphs = [f for f in os.listdir(GRAPH_DIR) if f.endswith(".graph")]
    
    for g in graphs:
        path = os.path.join(GRAPH_DIR, g)
        res = subprocess.run(["./check_iso", path], capture_output=True, text=True)
        sig = res.stdout.strip()
        if sig not in signatures:
            signatures[sig] = []
        signatures[sig].append(g.replace(".graph", ""))
    
    return [mols for mols in signatures.values() if len(mols) > 1]

def generate_report(groups):
    now = time.strftime("%d/%m/%Y %H:%M:%S")
    
    # Début du HTML avec un CSS plus moderne
    html = f"""
    <html>
    <head>
        <meta charset="utf-8">
        <style>
            body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; max-width: 900px; margin: 20px auto; padding: 20px; background-color: #f4f7f6; }}
            h1 {{ color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px; }}
            .summary {{ background: #fff; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 20px; }}
            .group-card {{ background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-bottom: 15px; border-left: 5px solid #3498db; }}
            .mol-link {{ display: inline-block; background: #eef2f7; padding: 5px 10px; border-radius: 4px; margin: 5px; color: #2980b9; text-decoration: none; font-weight: bold; border: 1px solid #dcdfe6; }}
            .mol-link:hover {{ background: #3498db; color: #fff; }}
            .footer {{ font-size: 0.8em; color: #7f8c8d; margin-top: 40px; text-align: center; }}
            .no-dup {{ color: #27ae60; font-weight: bold; }}
        </style>
    </head>
    <body>
        <h1>Rapport d'isomorphisme moléculaire</h1>
        <div class="summary">
            <p><strong>Date d'analyse :</strong> {now}</p>
            <p><strong>Résultat :</strong> {f'<span class="no-dup">Aucun doublon structurel détecté.</span>' if not groups else f'<b>{len(groups)}</b> groupe(s) de doublons identifié(s).'}</p>
        </div>
    """

    if groups:
        for i, grp in enumerate(groups):
            html += f'<div class="group-card"><strong>Groupe {i+1} (Molécules identiques) :</strong><br><br>'
            for mol_name in grp:
                # On essaie d'extraire l'ID ChEBI du nom du fichier s'il existe
                # Rappel : nom formaté comme "mol_0_Nom_CHEBI_12345"
                chebi_url = "#"
                display_name = mol_name
                
                if "CHEBI" in mol_name.upper():
                    # On extrait juste la partie CHEBI_XXXX
                    parts = mol_name.split('_')
                    for p in parts:
                        if p.upper().startswith("CHEBI"):
                            chebi_url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={p}"
                            break
                
                html += f'<a href="{chebi_url}" target="_blank" class="mol-link" title="Voir sur ChEBI">📄 {display_name}</a>'
            html += "</div>"

    html += f"""
        <div class="footer">Généré par le projet M2 AMIS - Algorithme de McKay (Nauty)</div>
    </body>
    </html>
    """

    with open(REPORT_FILE, "w", encoding="utf-8") as f:
        f.write(html)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(1)
    
    setup()
    source = download_data(sys.argv[1])
    split_molecules(source)
    process_conversions()
    compile_engine()
    results = find_isomorphs()
    generate_report(results)