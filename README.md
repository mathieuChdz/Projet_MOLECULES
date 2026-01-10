# Projet : Détecteur d'Isomorphisme Moléculaire (M2 AMIS)

Ce projet permet de récupérer des bases de données chimiques (SDF), d'en extraire les molécules individuelles, et d'identifier les doublons structurels en utilisant l'algorithme de **McKay (Nauty)**.

---

## Installation

### 1. Prérequis communs

* **Python 3.10+**
* **GCC** 
* **Nauty** : La bibliothèque doit être présente dans le dossier `./nauty` (avec le fichier `nauty.a`).

### 2. Installation des dépendances Python

Ouvrez votre terminal (MINGW64 sur Windows ou Terminal sur Mac) :

```bash
pip install -r requirements.txt

```

---

## Guide d'utilisation

Le projet est entièrement automatisé via un `Makefile`.

### Sur Windows (via MINGW64 / MSYS2)

```bash
# 1. Nettoyer les anciens tests
mingw32-make clean

# 2. Compiler et lancer le pipeline avec un url exemple https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1983,1983,2519,5793/SDF
mingw32-make run URL="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1983,1983,2519,5793/SDF"

```

### Sur macOS / Linux

```bash
# 1. Nettoyer les anciens tests
make clean

# 2. Compiler et lancer le pipeline avec un url exemple https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1983,1983,2519,5793/SDF
make run URL="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1983,1983,2519,5793/SDF"

```

---

## Fonctionnement du Pipeline

Le projet suit un flux rigoureux divisé en 5 étapes clés :

1. **Récupération (Python/Requests)** : Téléchargement du fichier `.sdf` source et stockage local dans `data/`.
2. **Parsing (Python/RDKit)** : Découpage du fichier SDF. Chaque molécule est isolée dans un fichier `.mol` unique avec un index (`mol_i_...`) pour éviter tout écrasement.
3. **Conversion (Python/ProcessPool)** : Transformation des fichiers `.mol` en fichiers `.graph` personnalisés. Cette étape est parallélisée pour traiter de gros volumes de données.
4. **Analyse d'Isomorphisme (C/Nauty)** :
* Le programme C `check_iso` lit le fichier `.graph`.
* Il utilise la bibliothèque **nauty** pour générer une **signature canonique** unique basée sur la topologie de la molécule.


5. **Rapport (Python/HTML)** : Les molécules ayant la même signature sont regroupées. Un rapport `resultats.html` est généré avec des liens vers la base de données **ChEBI** pour vérification.

---

## Fondements de la structure

La représentation moléculaire est basée sur la théorie des graphes. Une molécule est modélisée par un graphe  où :

* **Sommets ()** : Représentent les atomes. Chaque sommet est "coloré" (étiqueté) avec son numéro atomique (ex: 6 pour le Carbone, 8 pour l'Oxygène).
* **Arêtes ()** : Représentent les liaisons chimiques entre les atomes.

### Pourquoi l'isomorphisme de McKay ?

Deux fichiers SDF peuvent décrire la même molécule avec une numérotation d'atomes différente. Pour prouver que , nous cherchons si les graphes  et  sont isomorphes.

L'algorithme de McKay est utilisé ici pour transformer chaque graphe en une **forme canonique** .



Cette approche permet une comparaison ultra-rapide ( après calcul de la signature) plutôt que de comparer chaque paire de molécules ().

---

## Structure du projet

* `main.py` 
* `process_mol.py` : Convertisseur Molécule -> Graphe.
* `check_iso.c` : Calcul utilisant nauty.
* `Makefile` : Scripts d'automatisation.
* `data/` : Stockage des molécules et des graphes.
* `nauty/` : Bibliothèque nauty de Brendan McKay.

---
