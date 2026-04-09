#making the resins  


from rdkit import Chem
from rdkit.Chem import AllChem

mols = {
    "dgebf": "C(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1",
    "detda": "CCc1cc(CC)c(N)cc1Cc1cc(CC)cc(CC)c1N"
}

for name, smi in mols.items():
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToXYZFile(mol, f"{name}.xyz")

