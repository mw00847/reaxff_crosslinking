"""
generate 3D molecular structures of polyester chain and HMMA crosslinker

varying the structure of the "polyester" could be really interesting 

OH-NPG-IPA-NPG-IPA-NPG-OH MW~573 g/mol

HMMA

"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


#specify the smiles 

polyester_smiles=(
    "OCC(C)(C)COC(=O)c1cccc(C(=O)OCC(C)(C)COC(=O)c2cccc(C(=O)OCC(C)(C)CO)c2)c1"
)

hmma_smiles="COCN(COC)c1nc(N(COC)COC)nc(N(COC)COC)n1"


#generate the 3D, use a build function with the smiles and name as inputs

def build_3d(smiles, file_name):
    #create an object of the smiles 
    mol=Chem.MolFromSmiles(smiles)
    
    #explicitly add hydrogens
    mol=Chem.AddHs(mol)
    
    #generate starting 3D coordinates with ETKDGv3 
    #include random seed for repeatability
    
    param_ETKDGv3 = AllChem.ETKDGv3()
    param_ETKDGv3.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, param_ETKDGv3)
    
    #if that doesnt work use ETDG
    if result == -1:
        AllChem.EmbedMolecule(mol,AllChem.ETDG())
        

    #minimise the starting coordinates with MMFF
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    
    #extract the 3D coordinates to make the xyz 
    
    conf    = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()
    formula = CalcMolFormula(mol)
    mw      = sum(atom.GetMass() for atom in mol.GetAtoms())
    
    
    #write the xyz file and save it 
    xyz_lines = [str(n_atoms), f"{file_name}"]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz_lines.append(
            f"{atom.GetSymbol():2s}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}"
        )
 
    with open(f"{file_name}.xyz", "w") as f:
        f.write("\n".join(xyz_lines))
 
    print(f"  {file_name}: {n_atoms} atoms | {formula} | MW = {mw:.1f} g/mol -> {file_name}.xyz")
    return mol, mw
    
    
#run the build function for both the polyester and the crosslinker 
pe_mol, pe_mw = build_3d(polyester_smiles, "polyester_ipa_npg")
    
hmma_mol, hmma_mw =build_3d(hmma_smiles,"hmma")
    
    


        
        
        
        
    
    
    
    
    
    
    
    
    
    
    


