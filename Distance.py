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
    """
    Génère la 3D uniquement si la molécule n'en possède pas déjà une.
    RENVOIE L'OBJET MOLÉCULE + LES IDS.
    """
    mol = Chem.MolFromMolFile(mol_path, removeHs=False)
    if mol is None:
        return None, []

    has_3d = False
    if mol.GetNumConformers() > 0:
        conf = mol.GetConformer()
        if conf.Is3D():
            for i in range(mol.GetNumAtoms()):
                if conf.GetAtomPosition(i).z != 0.0:
                    has_3d = True
                    break

    m_h = Chem.AddHs(mol, addCoords=True)

    if has_3d:
        ids = [0]
    else:
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        params.maxIterations = 1000
        params.pruneRmsThresh = 0.3
        ids = list(AllChem.EmbedMultipleConfs(
            m_h, numConfs=n, params=params))

        try:
            nb_conf_keep = 6
            res = AllChem.UFFOptimizeMoleculeConfs(m_h, numThreads=0)
            ids = sorted(ids, key=lambda cid: res[cid][1])[
                :min(nb_conf_keep, len(ids))]
            valid_ids = []
            for cid in ids:
                try:
                    AllChem.UFFOptimizeMolecule(m_h, confId=cid)
                    valid_ids.append(cid)
                except (ValueError, RuntimeError, Exception):
                    pass
            ids = valid_ids
        except Exception:
            pass

    return m_h, ids


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


def get_similarities(data1, data2):
    """
    Calcule et renvoie séparément la similarité 2D et 3D dans un dictionnaire.
    """
    sim2d = Tanimoto2D(data1, data2)
    sim3d = 0.0

    # On calcule la 3D uniquement si les deux molécules ont pu générer des coordonnées 3D
    if data1[1] is not None and data2[1] is not None:
        sim3d = ShapeTanimoto(data1, data2)

    return {
        "sim2d": float(sim2d),
        "sim3d": float(sim3d)
    }
