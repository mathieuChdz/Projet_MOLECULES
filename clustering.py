#!/usr/bin/env python3
"""
Clustering simple de molécules (Tâche 1.2)
"""

import os
import numpy as np
from Distance import Tanimoto2D, score, ShapeTanimoto
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

MOL_DIR = "data/molecules"
OUTPUT_DIR = "clustering_results"

def main():
    # Créer le dossier de sortie
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Charger les molécules
    mol_files = sorted([os.path.join(MOL_DIR, f) for f in os.listdir(MOL_DIR) if f.endswith(".mol")])
    mol_names = [os.path.basename(f).replace(".mol", "") for f in mol_files]
    
    print(f"Clustering de {len(mol_files)} molécules...")
    
    # Calculer la matrice de distances
    n = len(mol_files)
    D = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            sim = Tanimoto2D(mol_files[i], mol_files[j], 0)
            dist = 1.0 - sim
            D[i, j] = dist
            D[j, i] = dist
    
    # Clustering hiérarchique
    condensed = squareform(D)
    Z = linkage(condensed, method='average')
    
    # Dendrogramme
    plt.figure(figsize=(14, 7))
    dendrogram(Z, labels=mol_names, leaf_font_size=8)
    plt.title("Clustering hiérarchique des molécules")
    plt.xlabel("Molécule")
    plt.ylabel("Distance")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/dendrogram.png", dpi=150)
    print(f"→ Dendrogramme sauvegardé: {OUTPUT_DIR}/dendrogram.png")
    
    # Découper en 10 clusters (adapté pour 150 molécules)
    num_clusters = 23
    clusters = fcluster(Z, num_clusters, criterion='maxclust')
    
    # Sauvegarder les résultats
    with open(f"{OUTPUT_DIR}/clusters.txt", "w") as f:
        for cluster_id in range(1, num_clusters + 1):
            members = [mol_names[i] for i in range(n) if clusters[i] == cluster_id]
            f.write(f"Cluster {cluster_id} ({len(members)} molécules):\n")
            for mol in members:
                f.write(f"  - {mol}\n")
            f.write("\n")
    
    print(f"→ Assignations sauvegardées: {OUTPUT_DIR}/clusters.txt")
    print("\nTerminé !")

if __name__ == "__main__":
    main()
