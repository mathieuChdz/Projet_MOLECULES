from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdMolAlign, rdShapeHelpers, Descriptors, Fragments, rdFingerprintGenerator
import inspect


def get_extended_fingerprint(mol, radius=2, nBits=2048):
    """
    Génère un fingerprint hybride : Morgan + Groupes Fonctionnels (Ertl) + Polarité.
    """

    mfgen = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=nBits)
    morgan_fp = mfgen.GetFingerprint(mol)

    fragment_functions = [func for name, func in inspect.getmembers(
        Fragments, inspect.isfunction) if name.startswith('fr_')]

    num_frag_bits = len(fragment_functions)

    total_size = nBits + num_frag_bits + 1

    combined_fp = DataStructs.ExplicitBitVect(total_size)

    on_bits = morgan_fp.GetOnBits()
    for bit in on_bits:
        combined_fp.SetBit(bit)

    offset = nBits
    for i, func in enumerate(fragment_functions):
        if func(mol) > 0:
            combined_fp.SetBit(offset + i)

    offset_polar = nBits + num_frag_bits
    tpsa = Descriptors.TPSA(mol)

    is_polar = tpsa > 0
    if is_polar:
        combined_fp.SetBit(offset_polar)

    return combined_fp


def Tanimoto2D(data1, data2):
    return DataStructs.TanimotoSimilarity(data1[0], data2[0])


def vecteur_fingerprint2D(mol_path):
    mol = Chem.MolFromMolFile(mol_path)
    return get_extended_fingerprint(mol)


def genere_shape(mol_path, n=20):
    # n est le nombre de conformères générés
    mol1 = Chem.MolFromMolFile(mol_path)

    # Ajout des Hydrogènes (essentiel pour la 3D et le volume)
    m1 = Chem.AddHs(mol1)

    # Paramètres de génération 3D
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True  # Aide parfois à la convergence
    params.maxIterations = 1000

    # Génération des conformères
    # embedMolecule retourne -1 si échec, on gère ça avec try/except ou check liste
    ids1 = AllChem.EmbedMultipleConfs(m1, numConfs=n, params=params)

    # Optimisation UFF (Force Field) pour relaxer la géométrie
    # On capture les erreurs potentielles d'optimisation (non-convergentes)
    try:
        for cid in ids1:
            AllChem.UFFOptimizeMolecule(m1, confId=cid)
    except ValueError:
        pass  # Continue avec la géométrie non optimisée si l'UFF échoue

    return (m1, ids1)


def ShapeTanimoto(data1, data2):
    m1 = data1[1]
    m2 = data2[1]

    best_sim = 0.0

    # Comparaison N x N conformères (peut être lent si n est grand)
    for id1 in data1[2]:
        for id2 in data2[2]:
            # O3A : Open3D Align. Aligne m2 sur m1.
            pyO3A = rdMolAlign.GetO3A(m2, m1, prbCid=id2, refCid=id1)
            pyO3A.Align()

            # Calcul de la similarité de forme (Shape Tanimoto Distance)
            # Attention: ShapeTanimotoDist renvoie une DISTANCE (0 = identique, 1 = différent)
            dist = rdShapeHelpers.ShapeTanimotoDist(
                m1, m2, confId1=id1, confId2=id2)
            sim = 1.0 - dist

            if sim > best_sim:
                best_sim = sim

    return best_sim


def score(data1, data2, a=0.5):
    """
    Score combiné.
    a : Poids de la 2D (0.0 à 1.0).
    """
    if a == 0:
        return ShapeTanimoto(data1, data2)
    elif a == 1:
        return Tanimoto2D(data1, data2)
    sim2d = Tanimoto2D(data1, data2)
    sim3d = ShapeTanimoto(data1, data2)

    return a * sim2d + (1 - a) * sim3d
