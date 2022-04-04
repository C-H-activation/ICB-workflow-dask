import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from datetime import datetime

# dask speficic for parallelism
from dask.distributed import Client
from dask.multiprocessing import get

# kallisto molecule
from project.package.kallisto import constructMolecule
# xtb gaps
from project.package.xtb import extractHomoLumoGap
# helper function for structure validation
from project.package.myhelper import extractSubstrate

# workflow specific tasks
from project.molecularSimilarity import molecularSimilarity
from project.createCH import createCH
from project.createStericPenalties import createStericPenalties
from project.sortCH import sortCH
from project.createTS import createTS
from project.createMLDeltas import createMLDeltas
from project.createGraphDeltas import createGraphDeltas
from project.createBarriers import createBarriers
from project.optimizeStructure import performXTBCalculation
from project.createBoltzmannWeightImage import createBoltzmannWeightImage

# results object
from project.results import Results

# helper functions
from project.header import createHeader
from project.package.bash import changeDirectory
from project.package.utility import existFile

# rich output formatter
from rich.console import Console
from rich.markdown import Markdown
from rich.live import Live
from rich.table import Table

# arguments via ArgumentParser
from argparse import ArgumentParser

# configs
from config import ml_path, xtb_ncpus

# define the workflow
def main():

    #######################################################################
    # general setup
    #######################################################################

    # xtb optimization level
    xtb_opt_level = "lax"

	# input arguments: SMILES, run_id, and number of CPUs
    parser = ArgumentParser()
    parser.add_argument('--smiles', dest='smiles', type=str)
    parser.add_argument('--name', dest='name', type=str)
    parser.add_argument('--run_id', dest='run_id', type=int)
    parser.add_argument('--cpus', dest='cpus', type=int)
    args = parser.parse_args()

    # sanity check for SMILES string for SMILES string
    if args.smiles == '':
        exit()

    # number of CPUs for xtb optimizations
    xtb_cpus = str(xtb_ncpus)

    # define dask Client
    client = Client(threads_per_worker=1, n_workers=args.cpus)

    # set start time
    start_time = datetime.now()

    #######################################################################
    # construct a rich console and create header
    #######################################################################
    console = Console()

    header = createHeader(args.smiles, args.name, args.run_id, args.cpus)
    md = Markdown(header)
    console.print(md)

    #######################################################################
    # worklfow begin
    #######################################################################

    # initialize a results object
    res = Results()

    #######################################################################
    # create all CH positions
    #######################################################################
    md = Markdown(" Generate all C-H positions")
    console.print(md)
    step_begin = datetime.now()
    create_ch = {'run': (createCH, args.smiles, args.run_id)}
    res.atoms = get(create_ch, 'run')
    step_end = datetime.now()
    console.print('time: ', step_end - step_begin)

    #######################################################################
    # optimize substrate with GFN2-xtb/ALPB(THF)
    #######################################################################
    md = Markdown(" Optimize compound using GFN2-xtb/ALP(THF)")
    console.print(md)
    step_begin = datetime.now()
    path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
    fileName = os.path.join(path, "smiles.xyz")
    xtb_args = ["--opt", xtb_opt_level, "--alpb", "thf", "-P", xtb_cpus]
    optimize_substrate = {'run': (performXTBCalculation, fileName, xtb_args, path)}
    res.substrate_energy = get(optimize_substrate, 'run')
    step_end = datetime.now()
    console.print('time: ', step_end - step_begin)
    
    #######################################################################
    # sort CH positions
    #######################################################################
    md = Markdown(" Prepare all structures for docking (sorting)")
    console.print(md)
    step_begin = datetime.now()
    sort_ch = {}
    for i, atom in enumerate(res.atoms):
        key = 'run-{}'.format(atom)
        sort_ch[key] = (sortCH, res.atoms, i, args.run_id)

    # save sorted substrate molecular graph
    # important to check if TS optimization retrains the structure
    res.substrate_graph = {}
    for i, atom in enumerate(res.atoms):
        key = 'run-{}'.format(atom)
        res.substrate_graph[atom] = get(sort_ch, key)
    step_end = datetime.now()
    console.print('substrate graphs: ', res.substrate_graph)
    console.print('time: ', step_end - step_begin)

    #######################################################################
    # create all transition states
    #######################################################################

    md = Markdown(" Dock compound into transition-state templates")
    console.print(md)

    # define templates of the trainsition states
    #res.templates = [1, 4, 7, 19, 22, 25]
    res.templates = [1]

    step_begin = datetime.now()
    create_ts = {}
    for i, atom in enumerate(res.atoms):
        key = 'run-{}'.format(atom)
        create_ts[key] = (createTS, res.atoms, i, args.run_id, res.templates)

    res.create_ts = {}
    for i, atom in enumerate(res.atoms):
        key = 'run-{}'.format(atom)
        res.create_ts[atom] = get(create_ts, key)
    step_end = datetime.now()
    console.print('time: ', step_end - step_begin)

    #######################################################################
    # GFN2-xtb/ALPB(THF) optimizations of trainistion states
    #######################################################################

    md = Markdown(" Optimize transition states with GFN2-xtb/ALPB(THF)")
    console.print(md)

    step_begin = datetime.now()
    # optimize substrate and fix catalyst
    res.xtb_energy = {}
    # extract HOMO-LUMO gaps
    res.xtb_gaps = {}
    # xtb structures to be ignored
    res.xtb_ignore = {}

    optimize_1 = {}
    xtb_args = ["--opt", xtb_opt_level, "--alpb", "thf", "--input", "constrain.inp", "-P", xtb_cpus]
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        for template in res.templates:
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, "newstructure.xyz")
            key = 'run-{}-{}-1'.format(atom, template)
            optimize_1[key] = (performXTBCalculation, fileName, xtb_args, tpath)

    for i, atom in enumerate(res.atoms):
        for template in res.templates:
            xtb_begin = datetime.now()
            key = 'run-{}-{}-1'.format(atom, template)
            tmp = get(optimize_1, key)
            xtb_end = datetime.now()
            if tmp != float(14):
                print(f' - opt 1: {atom} - {template} {xtb_end - xtb_begin}')
            else:
                print(f' - opt 1: {atom} - {template} failed')

    # optimize everything and fix a dihedral
    optimize_2 = {}
    xtb_args = ["--opt", xtb_opt_level, "--alpb", "thf", "--input", "fix.inp", "-P", xtb_cpus]
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        for template in res.templates:
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, "xtbopt.xyz")
            key = 'run-{}-{}-2'.format(atom, template)
            optimize_2[key] = (performXTBCalculation, fileName, xtb_args, tpath)

    for i, atom in enumerate(res.atoms):
        for template in res.templates:
            xtb_begin = datetime.now()
            key = 'run-{}-{}-2'.format(atom, template)
            res.xtb_energy[key] = get(optimize_2, key)
            xtb_end = datetime.now()
            if res.xtb_energy[key] != float(14):
                print(f' - opt 2: {atom} - {template} {xtb_end - xtb_begin}')
            else:
                print(f' - opt 2: {atom} - {template} failed')

    # HL gaps 1
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        for template in res.templates:
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, "xtb.out")
            key = 'run-{}-{}-2'.format(atom, template)
            res.xtb_gaps[key] = extractHomoLumoGap(fileName)
            
    # check if GFN2-xtb optimized structure is okay
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        res.xtb_ignore[atom] = []
        for template in res.templates:
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, 'xtbopt.xyz')
            extractSubstrate(fileName, 85, tpath)
            ttpath = os.path.join(tpath, 'submolecule.xyz')
            try:
                kallisto_molecule = constructMolecule(geometry=ttpath)
                graph = kallisto_molecule.get_bonds()
                for i, node in enumerate(graph):
                    if sorted(node) == res.substrate_graph[atom][i]:
                        continue
                    else:
                        res.xtb_ignore[atom].append(template)
            except:
                res.xtb_ignore[atom].append(template)

    # optimize substrate and fix catalyst (rotated)
    optimize_3 = {}
    xtb_args = ["--opt", xtb_opt_level, "--alpb", "thf", "--input", "constrain.inp", "-P", xtb_cpus]
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        for template in res.templates:
            template = str(template) + 'r'
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, "newstructure.xyz")
            key = 'run-{}-{}-3'.format(atom, template)
            optimize_3[key] = (performXTBCalculation, fileName, xtb_args, tpath)

    for i, atom in enumerate(res.atoms):
        for template in res.templates:
            xtb_begin = datetime.now()
            template = str(template) + 'r'
            key = 'run-{}-{}-3'.format(atom, template)
            tmp = get(optimize_3, key)
            xtb_end = datetime.now()
            if tmp != float(14):
                print(f' - opt 3: {atom} - {template} {xtb_end - xtb_begin}')
            else:
                print(f' - opt 3: {atom} - {template} failed')

    # optimize everything and fix a dihedral (rotated)
    optimize_4 = {}
    xtb_args = ["--opt", xtb_opt_level, "--alpb", "thf", "--input", "fix.inp", "-P", xtb_cpus]
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        for template in res.templates:
            template = str(template) + 'r'
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, "xtbopt.xyz")
            key = 'run-{}-{}-4'.format(atom, template)
            optimize_4[key] = (performXTBCalculation, fileName, xtb_args, tpath)

    for i, atom in enumerate(res.atoms):
        for template in res.templates:
            xtb_begin = datetime.now()
            template = str(template) + 'r'
            key = 'run-{}-{}-4'.format(atom, template)
            res.xtb_energy[key] = get(optimize_4, key)
            xtb_end = datetime.now()
            if res.xtb_energy[key] != float(14):
                print(f' - opt 4: {atom} - {template} {xtb_end - xtb_begin}')
            else:
                print(f' - opt 4: {atom} - {template} failed')

    # HL gaps 2
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        for template in res.templates:
            template = str(template) + 'r'
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, "xtb.out")
            key = 'run-{}-{}-4'.format(atom, template)
            res.xtb_gaps[key] = extractHomoLumoGap(fileName)

    # check if GFN2-xtb optimized structure is okay
    for i, atom in enumerate(res.atoms):
        path = os.path.join(os.getcwd(), 'run_{}'.format(args.run_id))
        path = os.path.join(path, str(atom))
        res.xtb_ignore[atom] = []
        for template in res.templates:
            template = str(template) + 'r'
            tpath = os.path.join(path, str(template))
            fileName = os.path.join(tpath, 'xtbopt.xyz')
            extractSubstrate(fileName, 85, tpath)
            ttpath = os.path.join(tpath, 'submolecule.xyz')
            if existFile(ttpath): 
                kallisto_molecule = constructMolecule(geometry=ttpath)
            else:
                res.xtb_ignore[atom].append(template)
            graph = kallisto_molecule.get_bonds()
            for i, node in enumerate(graph):
                if sorted(node) == res.substrate_graph[atom][i]:
                    continue
                else:
                    res.xtb_ignore[atom].append(template)

    # remove duplicates from ignore list
    for key, val in res.xtb_ignore.items():
        res.xtb_ignore[key] = sorted(set(val))
        
    step_end = datetime.now()
    console.print('Ignore those strctures: ', res.xtb_ignore)
    console.print('time: ', step_end - step_begin)

    #######################################################################
    # create xtb deltas from energies
    #######################################################################

    md = Markdown(" GFN2-xtb/ALPB(THF) minimum and penalties")
    console.print(md)

    step_begin = datetime.now()
    res.xtb_minimum = 0
    res.xtb_deltas = {}

    create_xtb_deltas = {'run': (createXTBDeltas, res.atoms, res.templates, res.substrate_energy, res.xtb_energy, res.xtb_gaps, res.xtb_ignore)}
    res.xtb_minimum, res.xtb_deltas = get(create_xtb_deltas, 'run')
    step_end = datetime.now()
    console.print('time: ', step_end - step_begin)

    #######################################################################
    # create graph deltas from molecular graph
    #######################################################################

    md = Markdown(" Molecular graph penalties")
    console.print(md)

    step_begin = datetime.now()
    res.graph_deltas = {}

    create_graph_deltas = {'run': (createGraphDeltas, args.smiles)}
    res.graph_deltas = get(create_graph_deltas, 'run')
    step_end = datetime.now()
    console.print('time: ', step_end - step_begin)

    #######################################################################
    # create machine-learning deltas using TS structures
    #######################################################################

    md = Markdown(" PLS-P2 machine-learning penalties")
    console.print(md)

    step_begin = datetime.now()
    res.ml_deltas = {}
    res.similarity = {}

    modelPath = ml_path
    create_ml_deltas = {'run': (createMLDeltas, res.atoms, args.run_id, res.templates, modelPath, res.graph_deltas)}
    res.ml_deltas, res.similarity = get(create_ml_deltas, 'run')
    step_end = datetime.now()
    console.print('time: ', step_end - step_begin)


    #######################################################################
    # create barriers in kJ/mol and weights in percent
    #######################################################################

    md = Markdown(" Final barriers and Boltzmann weights")
    console.print(md)

    step_begin = datetime.now()
    res.barriers = {}
    res.weights = {}

    res.molecular_similarity = 0
    create_barriers = {'run': (createBarriers, res.atoms, res.xtb_minimum, res.xtb_deltas, res.ml_deltas, res.graph_deltas, res.similarity, res.molecular_similarity)}
    res.barriers, res.weights = get(create_barriers, 'run')
    step_end = datetime.now()
    console.print('time: ', step_end - step_begin)

    #######################################################################
    # workflow end
    #######################################################################

    console.print()
    end_time = datetime.now()

    #######################################################################
    # Print all properties
    #######################################################################

    # create a results table including allC-H positions and details
    table = Table(show_header=True, header_style="bold magenta")
    table.add_column("[bold]C-H[/bold]", style="dim", width=12, justify="center")
    table.add_column("[bold]Similarity[/bold]", style="dim", width=16, justify="center")
    table.add_column("[bold]Neighbor penalty[/bold]", style="dim", width=16, justify="center")
    table.add_column("[bold]PLS-PF penalty[/bold]", style="dim", width=16, justify="center")
    table.add_column("[bold]Barrier[/bold]", style="dim", width=14, justify="center")
    table.add_column("[bold]Weight[/bold]", style="dim", width=14, justify="center")

    # reformat dictionaries and fill table with C-H information (rows)
    res.similarity = {int(k):v for k,v in res.similarity.items()}
    res.graph_deltas = {int(k):v for k,v in res.graph_deltas.items()}
    res.ml_deltas = {int(k):v for k,v in res.ml_deltas.items()}
    res.barriers = {int(k):v for k,v in res.barriers.items()}
    res.weights = {int(k):v for k,v in res.weights.items()}
    for atom in res.atoms:
        atom = int(atom)
        table.add_row(
            f"{atom}",
            f"{res.similarity[atom]:.2f}",
            f"{res.graph_deltas[atom]:.2f}",
            f"{res.ml_deltas[atom]:.2f}",
            f"{res.barriers[atom]:.2f}",
            f"{res.weights[atom]:.2f}",
        )

    console.print(table)
    console.print(f"GFN2-xtb minimum energy: {res.xtb_minimum:.2f} kJ/mol")
    console.print("All energies in kJ/mol. Boltzmann weights in percent.")

    #######################################################################
    # Save molecular image with calculated Boltzmann weights
    #######################################################################

    imgName = os.path.join(os.getcwd(), f'{args.name}.svg')
    createBoltzmannWeightImage(args.smiles, res.atoms, res.weights, imgName)

    console.print()
    console.print("Created molecular image with Boltzmann weights.")

    console.print()
    console.print('Overall calculation time: ', end_time - start_time)

if __name__ == '__main__':
    main()
