from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D, DataStructs

def get_3d_fingerprint(mol, nBits=1024):
    """Génère un vecteur de bits hybride 2D + descripteurs 3D."""
    if mol is None: return None
    
    # Préparation 3D (Conformères)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG()) # passe une molécule de la 2D à la 3D
    AllChem.MMFFOptimizeMolecule(mol)   # secoue la molécule
    
    # Morgan Fingerprint
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=nBits)
    bit_list = list(fp.ToBitString())
    
    # tests sur propriétés 3D
    npr1 = Descriptors3D.NPR1(mol)
    bit_list.append('1' if npr1 > 0.4 else '0')
    
    rg = Descriptors3D.RadiusOfGyration(mol)
    bit_list.append('1' if rg > 3.0 else '0')
    
    return DataStructs.CreateFromBitString("".join(bit_list))

def calcul_similarite(mol1, mol2):
    """Calcule le score Tanimoto entre deux molécules."""
    fp1 = get_3d_fingerprint(mol1)
    fp2 = get_3d_fingerprint(mol2)
    if fp1 is None or fp2 is None: return 0.0
    return DataStructs.TanimotoSimilarity(fp1, fp2)