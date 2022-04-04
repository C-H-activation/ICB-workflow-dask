# project/createStericPenalties.py

try:
    import numpy as np
except ImportError:
    raise ImportError(
         "Module numpy is required. Please install module numpy."
    )
try:
    import os
except ImportError:
	raise ImportError(
		"Module os is required. Please install module os."
	) 
try:
    import tempfile
except ImportError:
	raise ImportError(
		"Module tempfile is required. Please install module tempfile."
	) 

def extractCrestEnergy(path:str) -> float:
    """Extract the electronic energy from CREST geometry."""

    from project.package.utility import existFile
    if not existFile(path):
        return 1

    with open(path) as f:
        lines = f.readlines()

    return float(lines[1].split()[0])

def createConformationEnsemble(inp: str, out: str, dirp=os.getcwd()):
    """Create conformational ensemble with CREST using GFN-FF/ALPB(THF)."""

    try:
    	import subprocess
    except ImportError:
    	raise ImportError(
    		"The module subprocess is required. Please install the module subprocess."
    	)

    env = os.environ.copy()
    args = ["crest", inp, "--gfnff", "--T", "16", "--alpb", "thf"]
    output = os.path.join(dirp, out)
    with open(output, "w", newline=None) as outputFile:
        subprocess.call(
            args,
            shell=False,
            stdin=None,
            stderr=subprocess.STDOUT,
            universal_newlines=False,
            cwd=dirp,
            stdout=outputFile,
            env=env,
        )

def createStericalPenalty(maxProx: float, prox: float) -> float:
    """Input a proximity shell and return a sterical hindrance value in kJ/mol."""
    damping = 0.5 * (1 - np.tanh((3.5 - prox) / 2)) 
    return damping * 0.5 * (maxProx - maxProx * np.tanh((maxProx * 0.5 - prox) / 2))


def createStericPenalties(wdir: str, jobid: str, atoms: str):
    """Create proximity shells for substrate using kallisto."""

    try:
        import kallisto.reader.strucreader as ksr
    except ImportError:
        raise ImportError(
             "Module kallisto is required. Please install module kallisto."
        )
    try:
        import click
    except ImportError:
        raise ImportError(
             "Module click is required. Please install module click."
        )
    try:
        import json
    except ImportError:
    	raise ImportError(
    		"Module json is required. Please install module json."
    	) 
    try:
        import shutil
    except ImportError:
    	raise ImportError(
    		"Module shutil is required. Please install module shutil."
    	) 
    try:
        from collections import defaultdict
    except ImportError:
    	raise ImportError(
    		"Module collections is required. Please install module collections."
    	) 

    # debugging
    #import logging
    #logger = logging.getLogger("create_steric_penalties")

    from project.package.bash import changeDirectory, copyFile, createDirectory
    from project.package.utility import existDirectory
    from project.package.myhelper import grouper

    # path specific
    path = "run_{}".format(jobid)
    path = os.path.join(wdir, path)
    changeDirectory(destination=path)
    
    # create molecule
    dummy = click.File
    mol = ksr.constructMolecule(geometry="xtbopt.xyz", out=dummy)
    nat = mol.get_number_of_atoms()

    destination = os.path.join(path, 'crest')
    # check for previous CREST run
    if existDirectory(destination):
        # change to job directory
        changeDirectory(destination=destination)
    else:
	    # create CREST directory
        createDirectory(name='crest', path=path)

        # copy substrate file
        inp = os.path.join(path, "xtbopt.xyz")
        copyFile(source=inp, destination=destination)

        # change to job directory
        changeDirectory(destination=destination)
        
        # CREST GFN-FF/ALPB(THF) conformer generation
        createConformationEnsemble(inp=inp, out="crest.out", dirp=destination)

    # split CREST conformers
    prox = {}
    energies = {}
    delimiter = nat + 2
    with open('crest_conformers.xyz') as f:
        for i, g in enumerate(grouper(delimiter, f, fillvalue=None), 1):
            with tempfile.NamedTemporaryFile('w', delete=False) as fout:
                for j, line in enumerate(g, 1): # count number of lines in group
                    if line is None:
                        j -= 1 # don't count this line
                        break
                    fout.write(line)

            inp = os.path.join(destination, "conf_{0}.xyz".format(i))
            shutil.move(fout.name, inp)

            # initialize
            conf = str(i)
            prox[conf] = []
            energies[conf] = 0

            # extract conformer energy and proximity shells
            #logger.info("Input from CREST", inp)
            energies[conf] = extractCrestEnergy(path=inp)
            mol = ksr.constructMolecule(geometry=inp, out=dummy)
            nat = mol.get_number_of_atoms()
            # size is optimize for getting sterics of ortho-substituents
            size = (1, 1.8)
            prox[conf] = mol.get_prox(size=size)
            
    # find minumum conformer energy
    argmin = min(energies, key=energies.get)
    minimum = float(energies[argmin])
    
    # change to working directory
    changeDirectory(destination=wdir)

    # create conformer deltas in kJ/mol
    hartree2kj = 2625.6
    conformer = {}
    for k, v in energies.items():
        conf = str(k)
        conformer[conf] = (v - minimum) * hartree2kj
    #logger.info("conformer", conformer)

    # extract maximum proximity shell
    prox_list = []
    # over all conformers
    for k, v in prox.items():
        conf = str(k)
        shift = conformer[conf]
        # ignore relative barriers bigger than 10 kJ/mol
        if shift > 10:
            continue
        # over all atoms in each conformer
        for i, j in enumerate(v):
            if i in atoms:
                prox_list.append(j)
    limit = max(prox_list)

    # incorporate all conformers
    penalties = defaultdict(list)
    for k, v in prox.items():
        conf = str(k)
        shift = conformer[conf]
        # ignore relative barriers bigger than 10 kJ/mol
        if shift > 10:
            continue
        penalties[conf] = [createStericalPenalty(limit, x) + shift for x in v]
        #logger.info("conf, penalties", conf, penalties[conf])

    # extract minimum penalty for each atom
    penalty = {}
    for atom in atoms:
        tmp = []
        for k, v in penalties.items():
            tmp.append(v[atom])
        penalty[str(atom)] = min(tmp)

    # set penalties
    return penalty
