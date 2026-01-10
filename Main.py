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
    html = "<html><body style='font-family:sans-serif;'><h1>Rapport d'isomorphisme</h1>"
    if not groups:
        html += "<p>Aucun doublon structurel détecté.</p>"
    else:
        for grp in groups:
            html += "<div style='border:1px solid #ccc; margin:10px; padding:10px;'>"
            html += "<strong>Groupe identique :</strong> " + " == ".join(grp)
            html += "</div>"
    html += "</body></html>"
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