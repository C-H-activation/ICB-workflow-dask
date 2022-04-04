# project/createBarriers.py


def createWeights(atoms: list, barriers: dict, percent=False, T=298.15, R=0.0083145):
    """Create Boltzmann weights from barriers given in kJ/mol.
    T: temperature is given in Kelvin (default=298.15 K). 
    R: ideal gas constant given in kJ/mol (default=0.0083145 J/(K mol)."""

    try:
        import numpy as np
    except ImportError:
    	raise ImportError(
    		"Module numpy is required. Please install module numpy."
    	) 

    # initialize weights
    weights = {}

    norm = 0
    for i, atom in enumerate(atoms):
        expTerm = np.exp(- barriers[atom] / (R * T))
        weights[atom] = expTerm
        norm += expTerm

    # devide by norm to get weights
    for atom in atoms:
        weights[atom] /= norm

    if percent:
        for key, value in weights.items():
            weights[key] = value * 100

    return weights

def similaritySwitch(similarity: float, steepness: float) -> float:
    """Return float depending on the entered similarity."""
    try:
        import numpy as np
    except ImportError:
    	raise ImportError(
    		"Module numpy is required. Please install module numpy."
    	) 

    return (1 - np.tanh(steepness*(1 - similarity)))

def createBarriers(
    atoms: list,
    minimum: float,
    xtb_deltas: dict,
    ml_deltas: dict,
    graph_deltas: dict,
    similarity: dict,
    molecular_similarity: float):
    """Create barriers for C-H activation for each position."""

    try:
        from statistics import mean
    except ImportError:
    	raise ImportError(
    		"Module statistics is required. Please install module statistics."
    	) 

    barriers = {}

    # use integer keys
    similarity = {int(k):v for k,v in similarity.items()}
    ml_deltas = {int(k):v for k,v in ml_deltas.items()}
    graph_deltas = {int(k):v for k,v in graph_deltas.items()}

    for atom in atoms:
        sim = similarity[atom]
        ml = ml_deltas[atom]
        graph = graph_deltas[atom]

        # for unhindered Hydrogens use ML barrier
        if graph != 0.0:
            theta = similaritySwitch(sim, 8.0)
            omtheta = 1 - theta

            ml *= theta
            graph *= omtheta

        mixed = ml + graph
        barriers[atom] = mixed + minimum

    #  T = 298.15 K and R = 0.0083145 J/(K mol)
    weights = createWeights(atoms=atoms, barriers=barriers, percent=True) 

    return barriers, weights
