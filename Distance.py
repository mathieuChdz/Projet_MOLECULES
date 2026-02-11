from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs,rdMolAlign,rdShapeHelpers



def Tanimoto2D(m1, m2):
    mol1 = Chem.MolFromMolFile(m1)
    mol2 = Chem.MolFromMolFile(m2)

    if mol1 is None or mol2 is None: return 0.0

    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

    fp1 = gen.GetFingerprint(mol1)
    fp2 = gen.GetFingerprint(mol2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def ShapeTanimoto(m1, m2, n=20):
    #n cest le nombre de generation de la structure 3D

    mol1 = Chem.MolFromMolFile(m1)
    mol2 = Chem.MolFromMolFile(m2)

    if mol1 is None or mol2 is None: return 0.0

    #Ajout des Hydrogene pour la structure 3D
    m1 = Chem.AddHs(mol1)
    m2 = Chem.AddHs(mol2)  

    #algo de gen de structure 3D
    params = AllChem.ETKDGv3()
    #Generation de la structure 3D , n fois et on recupere les ids de chaque structure 3D generer pour les deux molecule
    ids1 = AllChem.EmbedMultipleConfs(m1, numConfs=n, params=params)
    ids2 = AllChem.EmbedMultipleConfs(m2, numConfs=n, params=params)

    if len(ids1) == 0 or len(ids2) == 0: return 0.0

    #Universal Force Field optimisation pour les deux molecule pour stabliser la structure 3D
    for id in ids1:
        AllChem.UFFOptimizeMolecule(m1, confId=id)
    for id in ids2:
        AllChem.UFFOptimizeMolecule(m2, confId=id)
    

    best_sim = 0.0
    for id1 in ids1:
        for id2 in ids2:
            #aligner leur strcture , m2 sur m1 on utilise l'algorithme O3A pour aligner les deux molecule en utilisant les id
            pyO3A = rdMolAlign.GetO3A(m2, m1, prbCid=id2, refCid=id1)
            pyO3A.Align()
            #distance volumique
            dist = rdShapeHelpers.ShapeTanimotoDist(m1, m2, confId1=id1, confId2=id2)
            sim=1 - dist
            if sim > best_sim:
                best_sim = sim

    return best_sim
    


def score(m1, m2, a):
    return  a * Tanimoto2D(m1, m2) + (1-a) * ShapeTanimoto(m1, m2)
