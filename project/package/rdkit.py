# dags/project/package/rdkit.py

try:
    import os
except ImportError:
    raise ImportError(
        "Module os is required. Please install module os."
    )
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import GetPeriodicTable
except ImportError:
    raise ImportError(
         "Module rdkit is required. Please install module rdkit."
    )
try:
    from kallisto.atom import Atom
    from kallisto.molecule import Molecule
    from kallisto.units import Bohr
except ImportError:
    raise ImportError(
         "Module kallisto is required. Please install module kallisto."
    )

from project.package.bash import changeDirectory

def smilesToXYZ(smiles: str, xyz_name: str, path=os.getcwd()) -> str:
    """Create xmol file from SMILES. Return path to xyz file."""

    source = os.getcwd()
    changeDirectory(path)
    rdkit_molecule = Chem.MolFromSmiles(smiles)
    mh = Chem.AddHs(rdkit_molecule)
    AllChem.EmbedMolecule(mh)
    xyz = Chem.rdmolfiles.MolToXYZBlock(mh).split("\n")
    # remove empty line from list
    xyz = [line for line in xyz if line != ""]
    # remove number of atoms
    xyz = xyz[1:]

    pt = GetPeriodicTable()

    # create atom list
    atoms = []
    for coord in xyz:
        elem, x, y, z = coord.split()[:4]
        position = [float(x) / Bohr, float(y) / Bohr, float(z) / Bohr]
        atom = Atom(symbol=pt.GetAtomicNumber(elem), position=position)
        atoms.append(atom)
    # construct kallisto molecule
    kallisto_molecule = Molecule(symbols=atoms)
    # write xmol file
    kallisto_molecule.writeMolecule(name=xyz_name, path=path)
    changeDirectory(source)
    return os.path.join(path, xyz_name)
