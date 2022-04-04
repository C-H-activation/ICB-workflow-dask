# project/package/kallisto.py

from kallisto.atom import Atom
from kallisto.molecule import Molecule
from kallisto.units import Bohr

def constructMolecule(geometry: str) -> Molecule:
    """Helper function to construct a Molecule."""
    try:
        with open(geometry, "r+") as fileObject:
            # read atoms from input file
            atoms = read_xyz(fileObject)
            # create molecule from atoms
            molecule = Molecule(symbols=atoms)
    except FileNotFoundError:
        raise FileNotFoundError("Input file not found.")
    return molecule


def read_xyz(f):
    """Read xmol file and return list of atoms."""
    atoms = []
    lines = f.readlines()
    nat = int(lines[0])
    for line in lines[2 : nat + 2]:
        atom, x, y, z = line.split()[:4]
        symbol = atom.strip()[0].upper() + atom.strip()[1:].lower()
        position = [float(x) / Bohr, float(y) / Bohr, float(z) / Bohr]
        atoms.append(Atom(symbol=symbol, position=position))
    return atoms
