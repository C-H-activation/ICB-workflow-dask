# dags/project/package/myhelper.py

try:
    import os
except ImportError:
    raise ImportError(
        "Module os is required. Please install module os."		
    )
try:
    import click
except ImportError:
    raise ImportError(
        "Module click is required. Please install module click."		
    )
try:
    from collections import defaultdict
except ImportError:
    raise ImportError(
        "Module collections is required. Please install module collections."		
    )
try:
	from itertools import zip_longest
except ImportError:
	raise ImportError(
		"The module itertools is required. Please install the module itertools."
	)
try:
    import kallisto.reader.strucreader as ksr
except ImportError:
    raise ImportError(
        "Module kallisto is required. Please install module kallisto."		
    )
try:
    from kallisto.data import chemical_symbols
except ImportError:
    raise ImportError(
        "Module kallisto is required. Please install module kallisto."		
    )
try:
    from kallisto.units import Bohr
except ImportError:
    raise ImportError(
        "Module kallisto is required. Please install module kallisto."		
    )


class MolecularGraph():
    """Define a molecular graph.
	
    Methods:

     addEdge - adds edge to graph
     BFS - breadth first search sorting
	 
    Variables:

     inp (str) - input xmol structure
     out (str) - output xmol structure (BFS sorted)
    """

    def __init__(self, inp: str, out: str):
        self.inp = inp
        self.out = out
        self.graph = defaultdict(list)
        self.molecule = ksr.constructMolecule(geometry=self.inp, out=click.File)
        self.nat = self.molecule.get_number_of_atoms()
        self.at = self.molecule.get_atomic_numbers()
        self.coordinates = self.molecule.get_positions()

    # method to add an edge to graph
    def addEdge(self, u: int, v: int):
        self.graph[u].append(v)

    # method to perform breadth first search
    def BFS(self, s: int):
    
        # mark all vertices as not visited
        visited = [False] * (len(self.graph))

        # create a queue
        q = []

        # mark source node as visited and enqueue it
        q.append(s)
        visited[s] = True

        seperator = os.linesep
        with open(self.out, 'w') as outFile:
            outFile.write(f'{self.nat}' + seperator)
            outFile.write('Created by Airflow' + seperator)
            while q:
                
                # dequeue a vertex from queue
                s = q.pop(0)
                outFile.write(
					"{:3} {:9.4f} {:9.4f} {:9.4f}".format(
                        chemical_symbols[self.at[s]],
                        self.coordinates[s][0] * Bohr,
                        self.coordinates[s][1] * Bohr,
                        self.coordinates[s][2] * Bohr,
                    ) + seperator		
                )

                # get adjacent vertices of dequeued vertex s
                # If an adjacent has not been visited, then mark
                # it as visited and enqueue it
                for i in self.graph[s]:
                    if visited[i] is False:
                        q.append(i)
                        visited[i] = True
    

def sortXYZ(inp: str, start_atom: int, out: str, path=os.getcwd()):
    """Sort xyz structure according to connectivity matrix from kallisto program.
    Atom declared by 'start' will be the zeroth atom in the sorted structure."""

    # we need to create a dummy click output file to serve
    # the implemented kallisto API
    output = click.File(mode='r', lazy=True)
    kmol = ksr.constructMolecule(geometry=inp, out=output)
    
    # get number of atoms
    nat = kmol.get_number_of_atoms()
    bonds = kmol.get_bonds()

    # initialize molecular graph
    output = os.path.join(path, out)
    g = MolecularGraph(inp, output)

    # travel through the graph
    for i in range(nat):
        partners = bonds[i]
        for j in partners:
            # add edge (j, i) to graph
            g.addEdge(j, i)
    
	# breath first search sorting: 
    #  start_atom will be zeroth atom in new structure
    #  output name defined by self.out
    g.BFS(start_atom)


def removeAtomByIndex(index: int, inp: str, out="wo.xyz", path=os.getcwd()):
    """Remove atom with index from structure."""

    from project.package.utility import existFile

    # check if input is available
    if not existFile(inp):
        return

    location = os.path.join(path, inp)
    f = open(location, 'r')
    lines = f.readlines()
    f.close()

    # reduce atom count
    nat = int(lines[0]) - 1
    # drop number of atoms line
    lines.pop(0)
    # increase index by 1
    index += 1
    # drop atom
    lines.pop(index)
    
    location = os.path.join(path, out)
    f = open(location, 'w')
    s = os.linesep
    f.write("{:>5}".format(nat) + s)
    for line in lines:
        f.write(line)
    f.close()


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks."
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def extractSubstrate(inp: str, start: int, path=os.getcwd()):
    """Cut away atoms and extract substrate with nat atoms from catalyst structure."""
    try:
        with open(inp, "r+") as f:
            lines = f.readlines()
	    
        nat = int(lines[0]) - start 
        xyz = lines[2 + start:]
        dirp = os.path.join(path, "submolecule.xyz")
        s = os.linesep
        with open(dirp, "w") as f:
            f.write(str(nat) + s)
            f.write(s)
            for coord in xyz:
                f.write(coord)
    except:
        return


def get_atom_neighbors(atom_idx: int, molecule_covalent_nbrs: list, alpha=None):
    """Get all neighbors for atom_idx.

	Extract all covalent bonding partner (alpha), all nearest neighbours (beta),
    and all nearest-nearest neighbors (gamma) for atom with index 'atom_idx'.

    Args:
    atom_idx: index of atom to extract neighbors for
    molecule_covalent_list: list of covalent partner of each atom

    Returns:
    alpha: list of all covalent bonding atom indices of atom_idx
    beta: list of nearest neighbor atom indices of atom_idx
    gamma: list of nearest-nearest neighbor atom indices of atom_idx

	"""
    # extract alpha neighbors
    if alpha is None:
        alpha = molecule_covalent_nbrs[atom_idx]

    # extract beta neighbors
    beta = list()
    for _, a in enumerate(alpha):
        b = molecule_covalent_nbrs[a]
        diff = list(set([atom_idx]) ^ set(b))
        if len(diff) > 0:
            beta.extend(diff)

    # extract gamma neighbors
    gamma = list()
    for _, b in enumerate(beta):
        c = molecule_covalent_nbrs[b]
        inter = list(set(alpha).intersection(set(c)))
        diff = list(set(inter) ^ set(c))
        gamma.extend(diff)
    gamma = list(dict.fromkeys(gamma))
    return alpha, beta, gamma
