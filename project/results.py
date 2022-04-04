# project/results.py

class Results(object):
    """Results objects saves all results from the graph."""

    def __init__(self):
        self.atoms = None
        self.molecular_similarity = None
        self.substrate_energy = None
        self.substrate_graph = None
        self.sort_ch = None
        self.templates = None
        self.create_ts = None
        self.xtb_energy = None
        self.xtb_minimum = None
        self.xtb_gaps = None
        self.xtb_deltas = None
        self.xtb_ignore = None
        self.similarity = None
        self.ml_deltas = None
        self.graph_deltas = None
        self.barriers = None
        self.weights = None
