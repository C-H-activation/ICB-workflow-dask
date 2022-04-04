# project/molecularSimilarity.py

def molecularSimilarity(smiles: str) -> int:
    """Calculate the similarity to the machine-learning training set."""

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise ImportError(
             "Module rdkit is required. Please install module rdkit."
        )
    try:
        from scipy.spatial.distance import rogerstanimoto
    except ImportError:
    	raise ImportError(
    		"Module scipy is required. Please install module scipy."
    	) 

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

    # dimension of vector
    dim = 256
    size = 2

    # get molecular similarity
    from project.package.db.smiles import DB2
    dbs = [Chem.MolFromSmiles(smi) for smi in DB2]
    fps = [AllChem.GetMorganFingerprintAsBitVect(x, size, dim) for x in dbs]
    target = AllChem.GetMorganFingerprintAsBitVect(substrate, size, dim)
    distances = []
    for i, fp in enumerate(fps):
        distance = 1 - rogerstanimoto(fp, target)
        distances.append(distance)

    # take maximum similarity
    molecular_similarity = max(distances)

    # save results
    return molecular_similarity
