# project/createGraphDeltas.py

try:
    import numpy as np
except ImportError:
    raise ImportError(
         "Module numpy is required. Please install module numpy."
    )
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError(
         "Module rdkit is required. Please install module rdkit."
    )
try:
    from kallisto.atom import Atom
    from kallisto.molecule import Molecule
    from kallisto.units import Bohr
    from kallisto.sterics import getClassicalSterimol
except ImportError:
	raise ImportError(
		"Module kallisto is required. Please install module kallisto."
	) 
try:
    from rdkit.Chem import GetPeriodicTable
except ImportError:
	raise ImportError(
		"Module rdkit is required. Please install module rdkit."
	) 

from project.package.myhelper import get_atom_neighbors


def stericalPenalty(l: float, bmin: float, bmax: float):
    """Create penalties from Sterimol descriptors via a multivariant
    linear regression model fitted to B3-LYP/def2-TZVP penalties as obtained
    for SiMe3 probe groups in ortho position to the substituent."""

    intercept = 15.24292747
    l_s = -0.23004497
    bmin_s = -7.76699641
    bmax_s = 4.95959058
    penalty = intercept + l_s*l + bmin_s*bmin + bmax_s*bmax
    if penalty < 0.0:
        return float(0)
    return penalty


def kallisto_molecule_from_rdkit_molecule(rdkit_molecule: Chem.rdchem.Mol) -> Molecule:
    """Create a kallisto molecule from RDKit molecule.
    Args:
        rdkit_molecule: RDKit molecule
    Returns:
        A kallisto molecule (kallisto.molecule.Molecule)
    """
    # get all xyz coordinates and split into list of lines
    xyz = Chem.rdmolfiles.MolToXYZBlock(rdkit_molecule).split("\n")
    # remove empty lines from list
    xyz = [line for line in xyz if line != ""]
    # remove number of atoms as given in xmol files (first line)
    xyz = xyz[1:]

    # setup periodic table
    pt = GetPeriodicTable()
    # create list of atoms
    atoms = []
    # create kallisto molecule
    for coord in xyz:
        elem, x, y, z = coord.split()[:4]

        # convert atomic coordinates from Angstrom to Bohr
        position = [float(x) / Bohr, float(y) / Bohr, float(z) / Bohr]
        atom = Atom(symbol=pt.GetAtomicNumber(elem), position=position)
        atoms.append(atom)
    return Molecule(symbols=atoms)


def atomBoundToExcluded(atom: Chem.rdchem.Atom) -> bool:
    """Returns if atoms is bounded to an excluded nuclear charge."""
    # exclude positions next to one of those nuclear charges
    # Nitrogen, Oxygen, and Sulfur
    excluded = [7, 8, 16]

    nbrs = atom.GetNeighbors()
    for nbr in nbrs:
        if nbr.GetAtomicNum() in excluded:
            return True
    return False


def getDistance(origin: list, partner: list):
    """Calculate the distance between origin and partner."""
    return np.sqrt(np.power(origin[0] - partner[0], 2) + np.power(origin[1] - partner[1], 2) + np.power(origin[2] - partner[2], 2))


def createGraphDeltas(smiles: str) -> dict:
    """Calculate the penalty depending on neighbors in molecular graph."""

    # create substrate from SMILES
    substrate = Chem.MolFromSmiles(smiles)

    # identify aromatic Carbon atoms from pattern
    # and unpack tuple to get list of aromatic Carbon atoms
    pattern = Chem.MolFromSmarts("c")
    aromaticCarbons = substrate.GetSubstructMatches(pattern)
    aromaticCarbons = [atom for (atom,) in aromaticCarbons]

    # get ring information and extract rings
    ri = substrate.GetRingInfo()
    rings = ri.AtomRings()

    # C-H count
    count, atoms = 0, []

    # setup initial Hydrogen count and shift
    hcount, shift = -1, 0

    # Check for each atom
    for index, atom in enumerate(substrate.GetAtoms()):
        shift += 1
        hcount += atom.GetNumImplicitHs()
        hcount += atom.GetNumExplicitHs()
        for ring in rings:
            # extract Carbon atoms that
            # 1) belong to a ring system
            # 2) are aromatic
            # 3) occur only once (exclude doubles)
            # 4) have ONE C-H bond.
            if(
                (index in ring)
                and (index in aromaticCarbons)
                and (index not in atoms)
                and (atom.GetNumImplicitHs() == 1)
				):
                    atoms.append(index)
                    count += 1

    # add Hydrogen and embed
    substrate = Chem.AddHs(substrate)
    AllChem.EmbedMolecule(substrate)

    # build a kallisto molecule, get number of atoms (nat)
    # and molecular graph (covbonds)
    from kallisto.rmsd import recursiveGetSubstructures
    mk = kallisto_molecule_from_rdkit_molecule(substrate)
    # get coordinates of molecule
    coords = mk.get_positions()
    # get molecular graph
    covbonds = mk.get_bonds()
    # get atomic number of molecules
    at = mk.get_atomic_numbers()
    # get coordination numbers
    cns = mk.get_cns(cntype='cov')
    # write xmol file (debug only)
    mk.writeMolecule('kallisto.xyz')
    # get number of atoms
    nat = mk.get_number_of_atoms()

    # get atoms
    rdatoms = substrate.GetAtoms()

    # get nuclear charges
    at = [atom.GetAtomicNum() for atom in rdatoms]

    # define atoms that connect rings from dict
    from collections import defaultdict
    ra = defaultdict(lambda: 0)
    for ring in rings:
        for atom in ring:
            ra[atom] += 1
    # sort out atoms that are only in one ring and dont penalize position for excluded atom types
    connecting_ring_atoms = dict((k, v) for k, v in ra.items() if v > 1)

    # get sorted rings with indices
    rings_sorted = defaultdict(lambda: [])
    for index, ring in enumerate(rings):
        rings_sorted[index] = sorted(ring)

    # iterate over all CHs
    sterical_penalty = {}

    for atom in atoms:
        # initialize
        sterical_penalty[atom] = {}

        # get RDKit atom object
        rdatom = rdatoms[atom]

        # extract neighbors of RDKit atom object
        nbrs = rdatom.GetNeighbors()

        # periodic table
        pt = Chem.GetPeriodicTable()

        # iterate over neighbors
        for nbr in nbrs:

            # sort out Hydrogen partners
            if nbr.GetAtomicNum() != 1:

                # extract substructures using kallisto
                substructures = recursiveGetSubstructures(nat, covbonds, nbr.GetIdx())

                for substructure in substructures:
                    # sort out the ring substructures (see below)
                    # keep substituent only (has no 'atom' in substructure)
                    if atom not in substructure:
                        nbr_covalent = [x.GetIdx() for x in nbr.GetNeighbors()]
                        # tuple unpack single element set to int
                        (nbr_neighbor,) = set(substructure).intersection(set(nbr_covalent))
                        # type cast to int neccessay
                        nn = int(nbr_neighbor)
                        # get nearest-nearest atom
                        nn_atom = rdatoms[nn]

                        # Only penalties for heavy atoms (not Hydrogen)
                        if rdatoms[nn].GetAtomicNum() != 1:
                            # create molecule for modified Sterimol
                            subatoms = []
                            subatoms.append(Atom(symbol=6, position=coords[nbr.GetIdx()][:]))
                            for idx in substructure:
                                subatoms.append(Atom(symbol=at[idx], position=coords[idx][:]))
                            submolecule = Molecule(symbols=subatoms)
                            # write substituent (debug only)
                            #submolecule.writeMolecule("sub.xyz")
							# Sterimol descriptors: L, Bmin, Bmax
                            l, bmin, bmax = getClassicalSterimol(submolecule, 0, 1)
                            if cns[nn] > 1.0:
                                # get Sterimol penalty from multivariant regression model
                                sterical_penalty[atom][nbr.GetIdx()] = stericalPenalty(l, bmin, bmax)
                                #print('nbr, neighbor substituent, L, bmin, bmax, penalty', int(nbr.GetIdx()), int(nn), l, bmin, bmax, sterical_penalty[atom][nbr.GetIdx()], cns[nn], cns[nbr.GetIdx()])
                            else:
                                sterical_penalty[atom][nbr.GetIdx()] = 0

                            ## Special cases
                            # increase penalty for atoms by sqrt(atom_mass)
                            if cns[nn] < 1.0:
                                print(sterical_penalty[atom][nbr.GetIdx()], np.sqrt(nn_atom.GetMass()))
                                sterical_penalty[atom][nbr.GetIdx()] += np.sqrt(nn_atom.GetMass())

    # create panalties
    graph_deltas = {}
    for atom in atoms:
        graph_deltas[atom] = 0
        penalties = sterical_penalty[atom]
        for _, penalty in penalties.items():
            graph_deltas[atom] += penalty

    # now we construct penalties for atoms next to connected ring atoms
    atoms_connecting_rings = []
    for atom in atoms:
        sterical_penalty[atom] = {}
        # get RDKit atom object
        rdatom = rdatoms[atom]

        # extract neighbors of RDKit atom object
        nbrs = rdatom.GetNeighbors()

        # iterate over neighbors
        for nbr in nbrs:

            # sort out Hydrogen partners
            if nbr.GetAtomicNum() != 1:

                # neighbor to connecting ring atom
                if nbr.GetIdx() in connecting_ring_atoms:
                    atoms_connecting_rings.append(rdatom.GetIdx())
                    #print('in ring, connecting', nbr.GetIdx(), rdatom.GetIdx())
                    alpha, beta, gamma = get_atom_neighbors(nbr.GetIdx(), covbonds, [rdatom.GetIdx()])
                    subatoms = []
                    subatoms.append(Atom(symbol=6, position=coords[nbr.GetIdx()][:]))
                    for idx in alpha:
                        subatoms.extend([Atom(symbol=at[idx], position=coords[idx][:])])
                    for idx in beta:
                        subatoms.extend([Atom(symbol=at[idx], position=coords[idx][:])])
                    for idx in gamma:
                        subatoms.extend([Atom(symbol=at[idx], position=coords[idx][:])])
                    submolecule = Molecule(symbols=subatoms)
                    # write submolecule (debug only)
                    # submolecule.writeMolecule("submolecule-2.xyz")
                    # Sterimol descriptors: L, Bmin, Bmax
                    l, bmin, bmax = getClassicalSterimol(submolecule, 0, 1)
                    # get a value penalty from multivariant regression model
                    sterical_penalty[atom][nbr.GetIdx()] = stericalPenalty(l, bmin, bmax)
     
    ## extend penalties with ortho to connecting ring atoms
    for atom in atoms_connecting_rings:
        penalties = sterical_penalty[atom]
        for _, penalty in penalties.items():
            graph_deltas[atom] += penalty
    return graph_deltas
