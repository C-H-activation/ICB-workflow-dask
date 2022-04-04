# project/sortCH.py

from project.package.kallisto import constructMolecule


def sortCH(atoms: str, index:int, jobid: str):
    """Sort xmol structures according to breadth first search algorithm.
    List 'atoms' contains all C-H positions for which structures should be
    sorted. The number tells us which atom should be on zeroth position when
    The structure is sorted.

    Example:
     atoms = [2, 4]

     There exist two C-H positions in the molecule (2 and 4).
     There should exist already a directory 'positions', which incorporates one 
     directory for each C-H position.

     Example from above:
      ./positions/2/2.xyz and ./positions/4/4.xyz

     Those structures MUST exist since they come from the first task in the DAGs definition.
     For structure ./position/2/2.xyz this present function will take the 2nd atom of the xmol file and
     will sort all entries such that the order will match a breadth first search order with respect
     to the connectivity of all atoms in the molecule."""

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
    
    # debugging
    #import logging
    #logger = logging.getLogger("sort_ch")
    #logger.info("Value of ch_positions {}. and type {}".format(atoms, type(atoms)))

    from project.package.myhelper import sortXYZ
    
    # extract atom by index
    atom = atoms[index]

    # get current woring directory
    cwd = os.getcwd()
    
	# path specific
    positions = "run_{}".format(jobid)
    
    # build path to files
    path = os.path.join(cwd, positions)
    path = os.path.join(path, str(atom))
    inp = os.path.join(path, 'wo.xyz')
    out = str(atom) + '_wo_sorted.xyz'
    sortXYZ(inp=inp, start_atom=atom, out=out, path=path)

    # construct kallisto Molecule from sorted structure
    fdir = os.path.join(path, out)
    kallisto_molecule = constructMolecule(geometry=fdir)
    graph = kallisto_molecule.get_bonds()
    return [sorted(x) for x in graph]
