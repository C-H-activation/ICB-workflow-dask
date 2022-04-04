# project/createTS.py

def createTS(atoms: str, index: int, jobid: str, templates: list):
    """Create transition-states and exchange substrates with kallisto"""

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
        import click
    except ImportError:
    	raise ImportError(
    		"Module click is required. Please install module click."
    	) 
    try:
        from kallisto.rmsd import exchangeSubstructure
        import kallisto.reader.strucreader as ksr
    except ImportError:
    	raise ImportError(
    		"Module kallisto is required. Please install module kallisto."
    	) 
    
    from project.package.bash import changeDirectory, copyFile, createDirectory
    from project.package.xtb import writeXTBConstraints
    from project.package.ts.store import writeTS

    # debugging
    #import logging
    #logger = logging.getLogger("create_TS")
    #logger.info("Value of ch_positions {}. and type {}".format(atoms, type(atoms)))
 
    # extract atom by index
    atom = atoms[index]

    # get current woring directory
    cwd = os.getcwd()
    
    # path specific
    positions = "run_{}".format(jobid)
    
    # dummy click File needed by kallisto API
    dummy = click.File

    # build path to file
    path = os.path.join(cwd, positions)
    path = os.path.join(path, str(atom))
    inp = os.path.join(path, str(atom) +'_wo_sorted.xyz')

    # define name of newstructure
    name = "newstructure.xyz"
    for template in templates:
        # create directory for template
        createDirectory(name=str(template), path=path)
        # directory name for rotated substrates
        rotated = str(template) + "r"
        # create directory for rotated substrate
        createDirectory(name=rotated, path=path)
        destination = os.path.join(path, str(template))
        rotatedDestination = os.path.join(path, rotated)
        ts_name = "ts{}.xyz".format(template)
        writeTS(n=template, name=ts_name, path=destination)
        writeTS(n=template, name=ts_name, path=rotatedDestination)

        # copy sorted substrate to transition-state directory
        copyFile(source=inp, destination=destination)
        copyFile(source=inp, destination=rotatedDestination)

        # change to transition-state directory
        changeDirectory(destination=destination)
        # create substrate object
        substrate = ksr.constructMolecule(geometry=inp, out=dummy)
        # number of atoms in substrate
        substrate_nat = substrate.get_number_of_atoms()
        # covalent bonding partners in substrate
        substrate_bonds = substrate.get_bonds()
        #logger.info("{} {} {}".format(substrate, inp, substrate_nat))
        # create transition-state template object
        ts = ksr.constructMolecule(geometry=ts_name, out=dummy)
        # number of atoms in template
        ts_nat = ts.get_number_of_atoms()
        # covalent bonding partners in template
        ts_bonds = ts.get_bonds()
        #logger.info("{} {} {}".format(ts, ts_name, ts_nat))

        # exchange benzene with substrate
        # benzene has substructure number 2
        benzene = 2
        # Iridium has index 18 in xmol
        iridium = 18
        mol = exchangeSubstructure(
            n=ts_nat,
            center=iridium,
            subnr=benzene,
            bonds=ts_bonds,
            ref=ts,
            newsub=substrate,
            newSubBonds=substrate_bonds,
            name="newstructure",
            rotate=0,
            exclude=False,
        )
        mol.writeMolecule(name=name, path=destination)
        writeXTBConstraints(index=template)
        #logger.info("{}".format(mol.get_number_of_atoms()))

        # change to transition-state directory with rotation
        changeDirectory(destination=rotatedDestination)

        # exchange benzene with substrate 
        # and rotate substrate by 180 degrees
        # benzene has substructure number 2
        benzene = 2
        # Iridium has index 18 in xmol
        iridium = 18
        mol = exchangeSubstructure(
            n=ts_nat,
            center=iridium,
            subnr=benzene,
            bonds=ts_bonds,
            ref=ts,
            newsub=substrate,
            newSubBonds=substrate_bonds,
            name="newstructure",
            rotate=180,
            exclude=False,
        )
        mol.writeMolecule(name=name, path=rotatedDestination)
        writeXTBConstraints(index=template)

        #logger.info("{}".format(mol.get_number_of_atoms()))
        # get back to working directory
        changeDirectory(destination=cwd) 

        # remove kallisto generated structure in wrong directory
        if os.path.exists(name):
            os.remove(name)
