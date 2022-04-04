# project/optimizeStructure.py


def extractInputName(inp: str) -> str:
    """Extract the input name from a path."""

    if "/" in inp:
        data = inp.split("/")
        inp = data[-1]

    return inp

def extractEnergy(inp: str) -> float:
    """Extract energy from optimized xtb structure."""
    from project.package.utility import existFile
    if existFile(inp):
        with open(inp, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if i == 1:
                    fields = line.split()
                    energy = fields[1]
                    break
        return float(energy)
    else:
        return float(14)


def performXTBCalculation(inp: str, args: list, dirp: str) -> str:
    """Perform an xTB calculation.
	
    Inputs:	
    inp: String including the input (e.g., structures/input.xyz).
    args: Arguments defining the xTB calculation.
    dirp: Name of temporary directory.

	"""
    try:
    	import subprocess
    except ImportError:
    	raise ImportError(
    		"workflow requires the module subprocess. Please install the module subprocess."
    	)
    try:
    	import os
    except ImportError:
    	raise ImportError(
    		"workflow requires the module os. Please install the module os."
    	)

    # save input name
    inputName = extractInputName(inp)

    # insert elements into arguments
    args.insert(0, "xtb")
    args.insert(1, inputName) 

    # get environment settings
    env = os.environ.copy()

    outputName = "xtb.out"
    outputPath = os.path.join(dirp, outputName)
    with open(os.path.join(dirp, outputName), "w", newline=None) as outputFile:
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

    path = os.path.join(dirp, "xtbopt.xyz")
    return extractEnergy(path)
