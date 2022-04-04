# dags/createMLDeltas.py

def createMLDeltas(atoms: str, jobid: str, templates: list, modelPath: str, penalty: dict):
    """Create deltas from machine-learning model (random forest) for all C-H positions."""

    try:
        import numpy as np
    except ImportError:
        raise ImportError(
            "Module numpy is required. Please install module numpy."
        )
    try:
        from collections import defaultdict
    except ImportError:
        raise ImportError(
            "Module collections is required. Please install module collections."
        )
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
    except ImportError:
        raise ImportError(
            "Module rdkit is required. Please install module rdkit."
        )
    try:
        from openbabel import pybel
    except ImportError:
        raise ImportError(
            "Module openbabel is required. Please install module openbabel."
        )

    from project.package.ml.data import NUCLEAR_CHARGE

    # define Compound object
    class Compound(object):
        """ The Compound object is used to store data."""
    
        def __init__(self, xyz=None):
        
            empty_array = np.asarray([], dtype = float)
    
            self.molid = float("nan")
            self.name = None
            self.smile = None
            self.subsmile = None
    
            # information about the compound
            self.nat = float("nan")
            self.at = {}
            self.atomtypes = []
            self.atomtype_indices = defaultdict(list)
            self.nuclear_charges = empty_array
            self.coordinates = empty_array
            
            # representation
            self.representation = empty_array
    
            if xyz is not None:
                self.read_xyz(xyz)
    
        def read_xyz(self, filename):
            """(Re-)initializes the Compound-object with data from an xyz-file.
    
               :param filename: Input xyz-filename.
               :type filename: string
            """
    
            f = open(filename, "r")
            lines = f.readlines()
            f.close()
    
            self.natoms = int(lines[0])
            self.atomtypes = []
            self.nuclear_charges = np.empty(self.natoms, dtype=int)
            self.coordinates = np.empty((self.natoms, 3), dtype=float)
    
            self.name = filename
    
            for i, line in enumerate(lines[2:self.natoms+2]):
                tokens = line.split()
    
                if len(tokens) < 4:
                    break
    
                self.atomtypes.append(tokens[0])
                self.atomtype_indices[tokens[0]].append(i)
                self.nuclear_charges[i] = NUCLEAR_CHARGE[tokens[0]]
        
                self.coordinates[i] = np.asarray(tokens[1:4], dtype=float)
       
            self.at = dict([(key, len(value)) for key,value in self.atomtype_indices.items()])
    
    
        def xyz_to_smiles(self):
            """Convert xyz file to smiles."""
            mol = next(pybel.readfile("xyz", self.name))
            smi = mol.write(format="smi")
            self.smile = smi.split()[0].strip()
    
    
        def get_representation(self, dim=128, size=4):
            """Convert xyz to smile and calculate Morgan representation for cutout."""
            # get SMILES from xyz
            self.xyz_to_smiles()
    	
            # extract substructure in complex via RDKit (turn off sanitization)
            rdmol = Chem.MolFromSmiles(self.smile, sanitize=False)

            # set Carbon and Iridium nuclear charge
            carbon, iridium = 6, 77

            for atom in rdmol.GetAtoms():
                # get Iridium atom index
                if atom.GetAtomicNum() == iridium:
                    # get neighbors of Iridum
                    neighbors = [x.GetAtomicNum() for x in atom.GetNeighbors()]
                    # sanity check: Carbon needs to be bound to Iridium
                    if carbon not in neighbors:
                        return
                    # Carbon position in neighbor list
                    carbonPosition = neighbors.index(carbon)
                    # get indices of neighbors in SMILES
                    indices = [x.GetIdx() for x in atom.GetNeighbors()]
    				# get aufpunkt: Carbon atom attached to Iridium
                    carbonIndex = indices[carbonPosition]
    
            # create submol in rdkit
            env = Chem.FindAtomEnvironmentOfRadiusN(rdmol, size, carbonIndex)
            submol = Chem.PathToSubmol(rdmol, env)
            submol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(submol, 
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                catchErrors=True,
            )
            Chem.GetSymmSSSR(submol)
            # set ECFP4 fingerprint as representation
            representation = AllChem.GetMorganFingerprintAsBitVect(submol, 2, nBits=dim)
            self.representation = np.array(representation)
            # set subsmile
            self.subsmile = Chem.MolToSmiles(submol)

    try:
    	import os
    except ImportError:
    	raise ImportError(
    		"Module os is required. Please install module os."
    	) 
    try:
        import json
    except ImportError:
    	raise ImportError(
    		"Module json is required. Please install module json."
    	) 
    try:
        import joblib
    except ImportError:
    	raise ImportError(
    		"Module joblib is required. Please install module joblib."
    	) 
    try:
        import pandas as pd
    except ImportError:
    	raise ImportError(
    		"Module pandas is required. Please install module pandas."
    	) 
    try:
        from scipy.spatial.distance import rogerstanimoto
    except ImportError:
    	raise ImportError(
    		"Module scipy is required. Please install module scipy."
    	) 

    from project.package.bash import copyFile, createDirectory
    from project.package.utility import existFile
    from project.package.myhelper import removeAtomByIndex

    # path specific
    positions = "run_{}".format(jobid)
    mlDir = os.path.join(positions, "predict")

    # initialize base directory for ML prediction
    createDirectory(name=mlDir)

    # load machine-learning model
    modelName = "pls.mod"
    model = joblib.load(os.path.join(modelPath, modelName))

    # load CSV data with binary training fingerprints
    csvName = os.path.join(modelPath, "db.csv")
    df = pd.read_csv(csvName)
    df_list = df.values.tolist()

    # dimension of vector
    dim = 256

    # copy all structures
    deltas = {}
    similarity = {}
    for i, atom in enumerate(atoms):
        # save available templates in array
        # some templates might not be converged properly
        available = []
        atom = str(atom)
        # initialize similarity
        similarity[atom] = 0
        # initialize delta
        deltas[atom] = 0
        mlAtomDir = os.path.join(mlDir, atom)
        createDirectory(name=mlAtomDir)
        for template in templates:
            template = str(template)
            # create path to xmol file
            source = os.path.join(positions, atom)
            source = os.path.join(source, template)
            source = os.path.join(source, "xtbopt.xyz")
            fileName = "{}.xyz".format(template)
            destination = os.path.join(mlAtomDir, fileName)
            # skip if optimization was not successful
            if existFile(source): 
                copyFile(source, destination)
                available.append(template)
            else:
                continue

        # create prediction vector
        X = []
        similarity_list = []
        for template in available:
            template = str(template)
            # create path to xmol file
            source = os.path.join(positions, atom)
            source = os.path.join(source, template)
            source = os.path.join(source, "xtbopt.xyz")
            fileName = "{}.xyz".format(template)
            destination = os.path.join(mlAtomDir, fileName)
            # remove Hydrogen atom bound to Iridium
            # necessary since Iridium should not be hypervalent for openbabel
            removeAtomByIndex(index=19, inp=destination, out=destination)
            # create compound and representation
            cpd = Compound(xyz=destination)
            #logger.info("cpd.at", cpd.at)
			# sanity check: compound creation
            if len(cpd.at) == 0:
                continue
            cpd.get_representation(dim=dim)
            rep = cpd.representation
			# sanity check: representation creation
            if len(rep) == 0:
                continue
            X.append(rep)

            # Rogers-Tanimoto similarity between vectors rep and row
            for i, row in enumerate(df_list):
                simTerm = 1 - rogerstanimoto(rep, row)
                similarity_list.append(simTerm)

        # extract similarity
        if similarity_list:
            similarity[atom] = max(similarity_list)

        # predict only when we have data
        if X:
            # create numpy array
            X = np.array(X)
            # predict deltas using randon forest
            yp = model.predict(X)
            # determine mean ML delta
            yp_mean = float(sum(yp) / len(yp))
            # sanity check
            # when predicted negative use graph penalty instead
            if yp_mean < -0.5:
                deltas[atom] = penalty[int(atom)]
            else:
                # save mean delta to atomic position
                deltas[atom] = yp_mean
        else:
            deltas[atom] = 14

    return deltas, similarity
