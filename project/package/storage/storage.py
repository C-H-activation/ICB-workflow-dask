# dags/project/storage/storage.py

class Storage:
    """This class saves results from the Iridium target airflow."""
    
    def __init__(
        self,
        runid=-1,
        hashid='',
        smiles='',
        atoms=[],
        molecular_similarity=0,
        minimum_barrier=0,
        positional_similarity={},
        xtb_deltas={},
        ml_deltas={},
        barriers={},
        weights={},
    ):
        self.runid = runid
        self.hashid = hashid
        self.smiles = smiles
        self.atoms = atoms
        self.molecular_similarity = molecular_similarity
        self.minimum_barrier = minimum_barrier
        self.positional_similarity = positional_similarity
        self.xtb_deltas = xtb_deltas
        self.ml_deltas = ml_deltas
        self.barriers = barriers
        self.weights = weights
