# project/header.py

def createHeader(smiles: str, name: str, run_id: str, cpus: str) -> str:
    """Create header of program."""
    header = """
        Regioselectivity determination for the Iridium catalyzed borylation
    
    	Compound       : {0}
        Name           : {1}
        SLURM ID       : {2}
    	Number of cpus : {3}
    
    """.format(smiles, name, run_id, cpus)
    return header
