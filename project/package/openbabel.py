# dags/project/package/openbabel.py

try:
    import os
except ImportError:
    raise ImportError(
        "Module os is required. Please install module os."
    )
try:
    from openbabel import pybel
except ImportError:
    raise ImportError(
         "Module openbabel is required. Please install module openbabel."
    )

from project.package.bash import changeDirectory

def smilesToXYZ(smiles: str, xyz_name: str, path=os.getcwd()) -> str:
    """Create xmol file from SMILES. Return path to xyz file."""

    source = os.getcwd()
    changeDirectory(path)
    inp = pybel.readstring("smi", smiles)
    inp.addh()
    inp.make3D()
    inp.write("xyz", xyz_name, overwrite=True)
    changeDirectory(source)

    return os.path.join(path, xyz_name)

