from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs



def Tanimoto2D(m1, m2):
    mol1 = Chem.MolFromMolFile(m1)
    mol2 = Chem.MolFromMolFile(m2)

    if mol1 is None or mol2 is None: return 0.0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def ShapeTanimoto(m1, m2):
    return 0.0

def score(m1, m2, a):
    return  a * Tanimoto2D(m1, m2) + (1-a) * ShapeTanimoto(m1, m2)
