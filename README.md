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

## Choix de la bibliothèque RDKit

L'utilisation de **RDKit** pour le parsing des fichiers SDF et la manipulation moléculaire a été privilégiée pour plusieurs raisons critiques :

### Robustesse du parsing

RDKit gère nativement les spécifications complexes des formats SDF (V2000 et V3000), incluant la gestion des valences et la détection d'aromaticité.

### Extraction de la table de connexion

L'API permet d'accéder directement à la matrice d'adjacence chimique. Contrairement à une simple lecture de texte, RDKit valide la structure (sanitization), garantissant que seuls des graphes chimiquement cohérents sont exportés vers le moteur C.

### Performance

La bibliothèque est écrite en C++, offrant une interface Python rapide capable de traiter des flux de données importants sans goulot d'étranglement lors de l'étape de fragmentation.

---

## Génération du format intermédiaire `.graph`

Le passage du fichier `.mol` au fichier `.graph` constitue l'étape de modélisation topologique. Une molécule est représentée par un graphe : G = (V, E)

### Codage des Sommets (V)

Chaque atome est mappé à un sommet. L'invariant choisi pour la partition initiale est le numéro atomique ( Z ).

**Justification :**
Pour que l'isomorphisme soit chimiquement exact, il est impératif de distinguer les types atomiques. Dans le fichier `.graph`, chaque sommet reçoit une étiquette numérique (ex. : `6` pour le carbone, `8` pour l’oxygène). Nauty utilise ensuite cette partition pour restreindre les permutations possibles aux seuls atomes de même nature.

### Codage des Arêtes (E)

Les liaisons covalentes sont traduites en arêtes non orientées. La topologie est ainsi extraite de ses coordonnées spatiales (2D ou 3D) pour ne conserver que l'information de connectivité pure.

---

## Mécanisme de comparaison par Signature Canonique

L'enjeu de l'isomorphisme est de s'affranchir de l'indexation arbitraire des atomes dans le fichier source. Pour prouver que $G_1 \simeq G_2$, le moteur calcule une forme canonique.

### Indépendance de l'ordre

Quelle que soit la manière dont la molécule est dessinée ou numérotée, Nauty génère une signature unique (certificat).

### Efficacité algorithmique

Cette approche transforme un problème de comparaison combinatoire complexe en une simple égalité de chaînes de caractères. Le coût de recherche de doublons dans le rapport final est ainsi réduit à une complexité quasi instantanée de ( O(1) ) par accès disque/mémoire.

---

## Structure du projet

* `main.py` 
* `process_mol.py` : Convertisseur Molécule -> Graphe.
* `check_iso.c` : Calcul utilisant nauty.
* `Makefile` : Scripts d'automatisation.
* `data/` : Stockage des molécules et des graphes.
* `nauty/` : Bibliothèque nauty de Brendan McKay.

---

## A essayer

```bash
mingw32-make run URL="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50/SDF"

```
---