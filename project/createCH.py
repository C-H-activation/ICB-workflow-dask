# project/createCH.py

def createCH(smiles: str, jobid: str) -> list:
    """Identify all aromatic C-H positions from SMILES and 
    cut away the Hydrogen atom attached to aromatic Carbon.
    Write substructures to directories.
    Return list of aromactic C-H atoms."""

    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError(
             "Module rdkit is required. Please install module rdkit."
        )
    try:
    	import os
    except ImportError:
    	raise ImportError(
    		"Module os is required. Please install module os."
    	) 

    from project.package.bash import copyFile, createDirectory
    from project.package.rdkit import smilesToXYZ
    from project.package.myhelper import removeAtomByIndex

    # path specific
    positions = "run_{}".format(jobid)

    # initialize base directory for method
    dirp = createDirectory(name=positions)

    # create substrate from SMILES
    substrate = Chem.MolFromSmiles(smiles)
    rdatoms = substrate.GetAtoms()

    # identify aromatic Carbon atoms from pattern
    # and unpack tuple to get list of aromatic Carbon atoms
    pattern = Chem.MolFromSmarts("c")
    aromaticCarbons = substrate.GetSubstructMatches(pattern)
    aromaticCarbons = [atom for (atom,) in aromaticCarbons]

    # get ring information and extract rings
    ri = substrate.GetRingInfo()
    rings = ri.AtomRings()

    # C-H candidates
    atoms, hydrogenNeighbors = [], []

    # setup initial Hydrogen count and shift
    hcount, shift = -1, 0

    nitrogen_nuclear_charge = 7
    # Check for each atom
    for index, atom in enumerate(substrate.GetAtoms()):
        shift += 1
        hcount += atom.GetNumImplicitHs()
        hcount += atom.GetNumExplicitHs()
        for ring in rings:
            #################################################
            # SPECIAL CASE for Nitrogen
            #################################################
            # sort out ortho to Nitrogen in six-membered ring
            if len(ring) == 6:
                if checkForNeighborAtom(ring, atom, nitrogen_nuclear_charge):
                    continue
            # sort out ortho to Nitrogen in five-membered ring
            # ONLY when more than one Nitrogen is in ring
            elif len(ring) == 5:
                at = [rdatoms[atom].GetAtomicNum() for atom in ring]
                if checkForAtomDuplicates(at, nitrogen_nuclear_charge):
                    if checkForNeighborAtom(ring, atom, nitrogen_nuclear_charge):
                        continue
            #################################################
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
                    hydrogenNeighbors.append(hcount)

    # Hydrogen atoms are added to the end of the xmol file
    # therefore shift all indices with respect to the number
    # of atoms that are not Hydrogen
    hydrogenNeighbors = [x + shift for x in hydrogenNeighbors]

    # prepare xmol structures
    xyz_name = "smiles.xyz"
    cwd = os.getcwd()
    dirp = os.path.join(cwd, positions)
    smilesToXYZ(smiles=smiles,xyz_name=xyz_name,path=dirp)
    # combine path and file name
    source = os.path.join(dirp, xyz_name)
    for i, k in enumerate(atoms):
        path = os.path.join(dirp, str(k))
        createDirectory(name=path)
        xyz_name = str(k) + ".xyz"
        destination = os.path.join(path, xyz_name)
        copyFile(source, destination)
        # remove Hydrogen atom connected to aromatic Carbon
        removeAtomByIndex(index=hydrogenNeighbors[i], inp=destination, path=path)

    return atoms

def checkForNeighborAtom(ring, atom, nuclear_charge) -> bool:
    """Check if atom is next to a specific Atom (True) or not (False)."""
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError(
             "Module rdkit is required. Please install module rdkit."
        )
    # get all neighbors
    nbrs = atom.GetNeighbors()
    # iterate over neighbors

    for nbr in nbrs:
        if nbr.GetAtomicNum() == nuclear_charge:
            return True
    return False

def checkForAtomDuplicates(atoms, nuclear_charge):
    """Check for duplicates of nuclear_charges in atoms list."""
    # check if nulcear charge exists at all
    if nuclear_charge not in atoms:
       return False 
    # set up a hash table to count occurence
    hash_dict={}
    for elem in atoms:
        if elem not in hash_dict:
           hash_dict[elem] = 1
        else:
           hash_dict[elem] += 1
     # nuclear charge is more than once in atom list
    if hash_dict[nuclear_charge] > 1:
        return True
    return False
