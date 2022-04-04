# dags/project/package/xtb.py

try:
    import os
except ImportError:
    raise ImportError(
        "Module os is required. Please install module os."		
    )

# define globale line seperator
s = os.linesep

def fixAtoms(atoms: list, name="fix.inp", path=os.getcwd()):
    """Write an xtb fixing file for passed in atoms."""
    
    # build string
    string = ', '
    string = string.join(map(str, atoms))

    with open(name, 'w') as f:
        f.write('$fix' + s)
        f.write(" atoms: {}".format(string) + s)
        f.write('$end' + s)


def writeXTBConstraints(index: int, path=os.getcwd()):
    """Write xtb contraints depending on transition-state template index."""

    if index in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
	    # fix dihedral: B(63)-Ir(19)-H(20)-C(86)
        fixAtoms(atoms=[63, 19, 20, 86])
        return

    # fix dihedral: B(42)-Ir(19)-H(20)-C(86)
    fixAtoms(atoms=[42, 19, 20, 86])
    return


def calculateBarrier(ts: float, cat: float, sub: float, conv=2625.5) -> float:
    """Calculate all barriers for C-H activation using energies for
    ts: transition-state energy in Hartree
    cat: catalyst energy in Hartree
    sub: substrate energy in Hartree.
    Returns by default a barrier in kJ/mol."""

    # default conversion factor for  Hartree -> kJ/mol
    return (ts - cat - sub) * conv


def extractHomoLumoGap(inp: str) -> float:
    """Extract the homo-lumo gap from an xtb output."""
    try:
        with open(inp, 'r+') as f:
            contents = f.readlines()
        target = 'HOMO-LUMO GAP'
        targets = [s for s in contents if target in s]
        # return homo-lumo gap
        return float(targets[0].split()[3])
    except:
        return float(0)
