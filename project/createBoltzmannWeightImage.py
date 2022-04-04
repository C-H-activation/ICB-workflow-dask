# project/createBoltzmannWeightImage.py

try:
    import os
except ImportError:
    raise ImportError(
         "Module os is required. Please install module os."
    )
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:
    raise ImportError(
         "Module rdkit is required. Please install module rdkit."
    )

def createBoltzmannWeightImage(smiles: str, atoms: list, weights: dict, name: str, wdir=os.getcwd(), print_weights=True, treshold=5):
    """Create a molecular SVG image with RDKit where the calculated
    Boltzmann weights are visible at the C-H positions."""

    # create RDKit molecule
    m = Chem.MolFromSmiles(smiles)
    m2 = Chem.Mol(m)

    # set weights 
    for atom in m2.GetAtoms():
        idx = atom.GetIdx()
        try:
            if weights[idx] >= treshold:
                lbl = '%.0f'%(weights[idx])
                if print_weights:
                    atom.SetProp('atomNote',lbl)
                else:
                    atom.SetProp('atomNote', str(idx))
        except:
            atom.SetProp('atomNote', '')

    # high quality SVG
    d2d = rdMolDraw2D.MolDraw2DSVG(400,400)
    d2d.DrawMolecule(m2)
    d2d.FinishDrawing()
    path = os.path.join(wdir, name)
    with open(path, 'w') as f:
        f.write(d2d.GetDrawingText())

