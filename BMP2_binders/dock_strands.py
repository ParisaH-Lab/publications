'''
Script from:
Danny D. Sahtoe
dsahtoe@uw.edu, danny.sahtoe@gmail.com
University of Washington, Institute for Protein Design
Baker Lab
Created Aug 2016
'''
#try:
import sys
import subprocess
import argparse
import glob
import itertools
import os
from operator import itemgetter
from itertools import groupby
from Bio.PDB import *
import pyrosetta
import numpy as np
from pprint import pprint
#except ImportError:
#    print('\nOne or more critical modules are not detected on your system. Exiting.')
#    sys.exit(1)

# Setup commandline options
parser = argparse.ArgumentParser()
parser.add_argument('-s', type=str, metavar='\b', default=None,
                    help='Specify input pdb', required=True)
parser.add_argument('-edge', type=str, metavar='\b', default=None,
                    help='Specify edgestrand range e.g. 12-14 to select resi 12,13,14. Needs to be at least 3 residues')
parser.add_argument('-database', type=str, metavar='\b', default=None, required=False,
                    help='Specify full path to database containing beta motifs')
parser.add_argument('-autodetect', action='store_true',
                    help='Autodetect edgestrands; cannot be used in conjunction with -edge option')
#parser.add_argument('-autodetect_only', action='store_true',
#                    help='Only output pml script with autodetected edgestrands highlighted')
parser.add_argument('-bb_only', action='store_true',
                    help='Ignore sidechains when detecting edgestrand for example when flexible side chain covers otherwise accesible edge. Useful in autodetect_mode')
parser.add_argument('-fastrelax', action='store_true',
                    help='Minimization at torsion level, NOT cartesian.')
parser.add_argument('-coordev', type=float, metavar='\b', default=None, required=False,
                    help='Strength of constraints in fastrelax, the lower the float the stronger the constraints. By default 1.0 if not set' )
args = parser.parse_args()


def fastrelax(pose,coordev):
    '''Kinematic minimization with Rosetta FastRelax. Min with coord csts'''
    #Setup fast_relax
    scorefxn = pyrosetta.get_fa_scorefxn()
    scorefxn(pose)
    fast_relax = pyrosetta.rosetta.protocols.relax.FastRelax(2)
    scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.coordinate_constraint,1.0)
    fast_relax.set_scorefxn(scorefxn)
    # setup coord cst
    addcsts = pyrosetta.rosetta.protocols.relax.AtomCoordinateCstMover()
    addcsts.cst_sd(coordev)
    addcsts.apply(pose)
    #setup task ops
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
    ex1_ex2 = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGeneric()
    ex1_ex2 = pyrosetta.rosetta.core.pack.task.operation.ExtraRotamersGeneric()
    ex1_ex2.ex1(True)
    ex1_ex2.ex2(True)
    tf.push_back(ex1_ex2)
    limitchi2aro = pyrosetta.rosetta.protocols.task_operations.LimitAromaChi2Operation()
    limitchi2aro.include_trp(True)
    limitchi2aro.chi2max(110)
    limitchi2aro.chi2min(70)
    tf.push_back(limitchi2aro)
    task = tf.create_task_and_apply_taskoperations(pose)
    #Setup common movemap elements i.e. jump=True and always alllow bb and sidechain min.
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    # Set jump to true if jump(s) present --> sets all jumps to True
    mm.set_chi(True)
    mm.set_bb(True)
    fast_relax.set_movemap(mm)
    fast_relax.set_task_factory(tf)
    fast_relax.apply(pose)
    print(('total_score',scorefxn(pose)))
    pyrosetta.dump_pdb(pose,'relax_renumbered_'+pdb.split('/')[-1])


def check_for_strands(pose):
    sheetsel = pyrosetta.rosetta.core.select.residue_selector.SecondaryStructureSelector(
        'E')
    sheetsel.set_minE(3)  # minimum length of a strand
    betaresidues = pyrosetta.rosetta.core.select.get_residues_from_subset(
        sheetsel.apply(pose))
    return betaresidues


def switch_renumber(pose, pdbname):
    switch = pyrosetta.rosetta.protocols.simple_moves.SwitchChainOrderMover()
    switch.chain_order('1')
    switch.apply(pose)
    pyrosetta.dump_pdb(pose, pdbname)


def get_abego(betastretch):
    abegos = pyrosetta.rosetta.core.sequence.get_abego(pose, 1)
    abegostretch = []
    for i in betastretch:
        abegostretch.append(abegos[i])
    return abegostretch


def check_longest_beta(beta_stretch):
    long_enough = False
    group = [list(j) for i, j in groupby(beta_stretch)]
    for i in group:
        if 'B' in i and len(i) >= 4:
            long_enough = True
    return long_enough


def new_get_edgestrands(pose, pdbname):
    # Get all beta residues in a list
    betaresidues = check_for_strands(pose)
    # Make dict where you for each beta residue store per atom sasa of NH 'H' and CO 'O'
    sasa = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
    sasa.calculate(pose)
    atomsasa = {}
    atoms = ['H', 'O']
    prolines = []
    for resn in betaresidues:
        resobj = pose.residue(resn)
        atomsasa[resn] = []
        if resobj.annotated_name() == 'P':
            prolines.append(resn)
            print(('Residue '+str(resn)+' is a PRO\n'))
            atomsasa[resn].append(sasa.get_atom_sasa()[
                                  resn][resobj.atom_index('O')])
            atomsasa[resn].append(0)
        else:
            for atomn in atoms:
                atomsasa[resn].append(sasa.get_atom_sasa()[
                                      resn][resobj.atom_index(atomn)])
    print('Atomic sasa:')
    print('betaresidue - NH - CO')
    pprint(atomsasa)
    # Put each stretch of at least 3 resi in dictionary
    beta_stretches = {}
    beta_stretches[pdbname.split('/')[-1]] = []
    i = 1
    for k, g in groupby(enumerate(betaresidues), lambda i_x: i_x[0]-i_x[1]):
        group = list(map(itemgetter(1), g))
        if not len(group) < 3:
            beta_stretches[pdbname.split('/')[-1]].append(group)
            i = i + 1

    # Calculate average sasa (NH and O) of each residue in each edge and put stretch above threshold in new dict
    all_long_enough_exposed_edgestrands = {}
    all_long_enough_exposed_edgestrands[pdbname.split('/')[-1]] = []
    # Seperate proline stretches and non proline stretches. Doing this makes it easier later
    proline_stretches_raw = {}
    remove_from_list = []
    for i in beta_stretches:
        for stretch in beta_stretches[i]:
            for pro in prolines:
                if pro in stretch:
                    proline_stretches_raw[pro] = []
                    proline_stretches_raw[pro].append(stretch)
                    remove_from_list.append(stretch)
    for i in beta_stretches:
        for j in remove_from_list:
            try:
                beta_stretches[i].remove(j)
            except:
                pass
    # Remove duplicates from proline stretches
    vals = []
    proline_stretches = {}
    for i in proline_stretches_raw:
        if not proline_stretches_raw[i][0] in vals:
            vals.append(proline_stretches_raw[i][0])
            proline_stretches[i] = []
            proline_stretches[i].append(proline_stretches_raw[i][0])
    print('\nBeta strands lacking prolines')
    pprint(beta_stretches)
    # Check average atomic sasa edgestrand to decide if it is exposed or not, threshold is...
    for i in beta_stretches:
        for j in beta_stretches[i]:
            avg = []
            abego_stretch = get_abego(j)
            for res in j:
                avg.append(atomsasa[res])
            if np.mean(avg) >= exposed_threshold:
                # print(('Exposed edge avg',np.mean(avg),i,j))
                j.append('{0:.4g}'.format(np.mean(avg)))
                if 'A' in abego_stretch or 'E' in abego_stretch or 'G' in abego_stretch:
                    if check_longest_beta(abego_stretch) == True:
                        all_long_enough_exposed_edgestrands[pdbname.split(
                            '/')[-1]].append(j)
                        # print(('Bulge present, long enough',i,j))
                    # else:
                    #     # print(('Bulge present, not long enough',i,j))
                else:
                    all_long_enough_exposed_edgestrands[pdbname.split(
                        '/')[-1]].append(j)
                    # print(('No bulge appending',i,j))
            elif np.mean(avg) <= exposed_threshold:
                j.append('{0:.4g}'.format(np.mean(avg)))
                # print(('This is a non exposed strand',np.mean(avg),i,j))

    print('\nBeta strands with prolines')
    pprint(proline_stretches)
    # Check average atomic sasa for edgestrands with PROLINE to decide if it is exposed or not, threshold is...
    for i in proline_stretches:
        for j in proline_stretches[i]:
            avg = []
            abego_stretch = get_abego(j)
            for res in j:
                avg.append(atomsasa[res])
            if np.mean(avg) >= exposed_threshold:
                # print(('Exposed edge avg',np.mean(avg),i,j))
                j.append('{0:.4g}'.format(np.mean(avg)))
                if 'A' in abego_stretch or 'E' in abego_stretch or 'G' in abego_stretch:
                    if check_longest_beta(abego_stretch) == True:
                        all_long_enough_exposed_edgestrands[pdbname.split(
                            '/')[-1]].append(j)
                        # print(('Pro and Bulge present,  long enough',i,j))
                    # else:
                    #     # print(('Pro and Bulge present,  not long enough',i,j))
                else:
                    if check_longest_beta(abego_stretch) == True:
                        all_long_enough_exposed_edgestrands[pdbname.split(
                            '/')[-1]].append(j)
                        # print(('Pro but No bulge and EEE appending',i,j))
                    # else:
                    #     # print(('Pro and No bulge and no EEE appending',i,j))
            elif np.mean(avg) <= exposed_threshold:
                j.append('{0:.4g}'.format(np.mean(avg)))
                # print(('This is a non exposed strand',np.mean(avg),i,j))

    global edges
    edges = all_long_enough_exposed_edgestrands
    print('\nExposed edgestrands:')
    pprint(edges)
    visualize_edgestrands_pymol()
    return all_long_enough_exposed_edgestrands


def visualize_edgestrands_pymol():
    '''Writes pymol script to visualize detected egdestrands'''
    colours = ['red', 'blue', 'yellow', 'violet', 'cyan',
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine',
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',
               'wheat', 'white', 'grey']
    count = 0
    for i in edges:
        with open(i+'.pml', 'w') as fout:
            fout.write('load '+name)
            fout.write('\nshow cartoon')
            fout.write('\nhide (h.)')
            fout.write('\ncolor gray90')
            for j, color in zip(edges[i], colours):
                avgatomsasa = edges[i][count][-1]
                fout.write('\nselect edge_'+str(count+1) +
                           '_AVGatomsasa'+avgatomsasa+'_'+i[:-4]+', i. ')
                for k in edges[i][count][:-1]:
                    fout.write(str(k)+'+')
                fout.write(' and '+i[:-4])
                count = count + 1
                fout.write('\ncolor '+color+', edge_'+str(count) +
                           '_AVGatomsasa'+avgatomsasa+'_'+i[:-4]+' and '+i[:-4])
            fout.write(
                '\ncolor blue, elem N\ncolor red, elem O\ncolor yellow, elem S')


def edgestrand_aligner(args, pose, name):
    '''Make edgestrand docks by aligning library of 2 stranded parallel and anti
    parallel sheets'''
    def aligner(target, edge_dictionary, path):
        # Main alignement loop, worst code ever, need to clarify/rewrite
        for model in edge_dictionary:
            for i in range(len(edge_dictionary[model])):
                print('edge '+str(i+1), edge_dictionary[model][i][:-1])
                edge_dictionary[model][i] = edge_dictionary[model][i][:-1]
                window = 3
                step = 1
                num_frags_target = (
                    ((max(edge_dictionary[model][i]) - min(edge_dictionary[model][i])) + 1) - window) / step + 1
                for lib in os.listdir(path):
                    strand_model = parser.get_structure(lib, path+'/'+lib)
                    print(path+'/'+lib)
                    for target_fragment in range(0, int(num_frags_target*step), step):
                        resi_to_be_aligned_target = list(range(min(
                            edge_dictionary[model][i])+target_fragment, min(edge_dictionary[model][i])+target_fragment+window))
                        # print resi_to_be_aligned_target
                        ca_T = []
                        for resi in resi_to_be_aligned_target:
                            ca_T.append(target[0]["A"][resi]["CA"])
                        # Start loop over strand fragments
                        # Make list with residues per chain
                        chs = ["A", "B"]
                        for ch in chs:
                            chain_atom_list = [int(str(res).split(" ")[4].split("=")[
                                                   1]) for model in strand_model for chain in model if chain.get_id() == ch for res in chain]
                            num_frags_strand = (
                                chain_atom_list[-1] - chain_atom_list[0] + 1 - window) / step + 1
                            for strand_fragment in range(0, int(num_frags_strand*step), step):
                                resi_to_be_aligned_strand = list(range(
                                    chain_atom_list[0]+strand_fragment, chain_atom_list[0]+strand_fragment+window))
                                # print resi_to_be_aligned_strand
                                ca_S = []  # Make list in which CA atoms of to be aligned residue are stored
                                for resi in resi_to_be_aligned_strand:
                                    ca_S.append(strand_model[0][ch][resi]["CA"])

                                # Make superimposer object
                                super_imposer = Superimposer()
                                super_imposer.set_atoms(ca_T, ca_S)
                                super_imposer.apply(strand_model.get_atoms())
                                # save
                                io = PDBIO()
                                io.set_structure(strand_model)
                                outpdb = 'aln_'+name.split('/')[-1]+'-edge'+str(i+1)+'-'+str(resi_to_be_aligned_target[0])+'-'+str(
                                    resi_to_be_aligned_target[-1])+'_'+lib[:-8]+ch+'-'+str(resi_to_be_aligned_strand[0])+'-'+str(resi_to_be_aligned_strand[-1])+'.pdb'
                                io.save(outpdb)
                                # Process outpdb so that only the non-aligned chain is in the pdb, remove TER/END and rename new strand to chain A for downstream processing
                                subprocess.check_output(
                                    "grep -v ' "+ch+" ' "+outpdb+" > temp.aln && mv temp.aln "+outpdb, shell=True)
                                subprocess.check_output(
                                    "sed '/TER/d' "+outpdb+" > temp && mv temp "+outpdb, shell=True)
                                subprocess.check_output(
                                    "sed '/END/d' "+outpdb+" > temp && mv temp "+outpdb, shell=True)
                                subprocess.check_output(
                                    "sed 's/ B / A /g' "+outpdb+" > temp && mv temp "+outpdb, shell=True)

    parser = PDBParser()
    if not args.database:
        print('\n location of database not specified, exiting\n')
        sys.exit(1)
    path = args.database
    if args.autodetect:
        edge_dictionary = new_get_edgestrands(pose, name)
        # get_beta_edgestrands() return edgedict and Pose in tuple, use only index 0 --> this is the dictionary
        if args.bb_only:
            target_in = name
        else:
            target_in = args.s
        target = parser.get_structure(target_in, target_in)
        aligner(target, edge_dictionary, path)
    if args.edge:
        if args.bb_only:
            target_in = name
        else:
            target_in = args.s
        target = parser.get_structure(target_in, target_in)
        edge_dictionary = {}
        target_strand = args.edge.split('-')
        target_strand = [int(i) for i in range(int(target_strand[0]),int(target_strand[1])+1)]
        assert int(target_strand[1]) > int(
            target_strand[0]), 'Start residue edge-strand must be smaller than stop residue edge-strand'
        target_strand.append('dummy_value')
        edge_dictionary[target_in] = [target_strand]
        aligner(target, edge_dictionary, path)


def strip_side_chains(pose, name):
    '''Mutate all side chains to alanine'''
#    pose = pyrosetta.io.pose_from_pdb(args.s)
    for i in range(1, len(pose.sequence())+1):
        mut = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(i, 'ALA')
        mut.apply(pose)
    pyrosetta.dump_pdb(pose, name)
    return pose


# main
pyrosetta.init('-ignore_unrecognized_res -mute all -beta_nov16 -renumber_pdb')
exposed_threshold = 2
pdb = args.s
line = '='
logname = len('Processing pdb: '+pdb)
print('\n')
print((logname*line))
print(('Processing pdb: '+pdb))
print((logname*line))
pose = pyrosetta.io.pose_from_file(pdb)

if args.autodetect and args.edge:
    print('\nUse either -autodetect or -edge flag but not both')
    sys.exit(1)
if args.autodetect:
    # if there are no strands at all in pose we can continue
    if len(check_for_strands(pose)) == 0:
        print(('no strands in pdb '+pdb))
    else:
        # If there are beta residues go through individual chains to find edges
        print('\nFinding edgestrands in '+pdb)
        print('____________________________________')
        if len(check_for_strands(pose)) == 0:
            print(('no strands in pose '+pdb))
        else:
            if args.bb_only:
                print(
                    '\n-bb_only option passed: Converting pose to polyalanine before finding edgestrands\n')
                name = 'polyala_'+pdb
                strip_side_chains(pose, name)
                switch_renumber(pose, name)
                edgestrand_aligner(args, pose, name)
            else:
                name = pdb
                switch_renumber(pose, name)
                edgestrand_aligner(args, pose, name)
if args.edge:
    if args.bb_only:
        print(
            '\n-bb_only option passed: Converting pose to polyalanine before finding edgestrands\n')
        name = 'polyala_'+pdb
        strip_side_chains(pose, name)
        edgestrand_aligner(args, pose, name)
    else:
        name = pdb
        edgestrand_aligner(args, pose, name)

if args.fastrelax:
    if args.coordev:
        fastrelax(pose,args.coordev)
    else:
        fastrelax(pose,1.0)
