'''
Script from:
Danny D. Sahtoe
dsahtoe@uw.edu, danny.sahtoe@gmail.com
University of Washington, Institute for Protein Design
Created Aug 2016
'''
import subprocess
import argparse
import sys
import pyrosetta

# Setup commandline options
parser = argparse.ArgumentParser()
parser.add_argument('-s',type=str, metavar='\b', default=None, help='Specify target pdb (required)',required=True)
parser.add_argument('-strand',type=str, metavar='\b', default=None, help='Strand pdb (required)',required=True)
parser.add_argument('-fix_target',action='store_true', help='Do not repack target interface')
parser.add_argument('-repack_target_residues',type=str, metavar='\b', default=None, help='comma separated list of residues to be repacked on target interface')
parser.add_argument('-threshold',type=str, metavar='\b', default=None, help='Edge-to-edge energies higher than this float will be filtered out (required)',required=True)
args = parser.parse_args()


def fastrelax(target_pdb,pose,args,prefix):
    # somehow if not dumping pose but passing it, residueselectors remain empty
    pdbinfo = pyrosetta.rosetta.core.pose.PDBInfo(pose)
    pose.pdb_info(pdbinfo)
    # dump_pdb(pose,args.strand+'_temp.pdb')
    # pose = pose_from_file(args.strand+'_temp.pdb')
    if len(pose.chain_sequence(2)) > 2:
        #Setup fast_relax
        scorefxn = pyrosetta.get_fa_scorefxn()
        scorefxn(pose)
        hbond_lr_bb_target = scorefxn.score_by_scoretype(target_pdb, pyrosetta.rosetta.core.scoring.hbond_lr_bb)
        fast_relax = pyrosetta.rosetta.protocols.relax.FastRelax(2)
        fast_relax.set_scorefxn(scorefxn)
        #Setup residue selectors interface, Neighbors near chain_B, within 8 A, but not chain_B
        chain_B = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")
        chain_A = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
        interfaceres_A = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(chain_B, 8, False )
        interfaceres_B = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(chain_A, 8, False )
        interface = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(interfaceres_A,interfaceres_B)
        fixnotinterface = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(interface)
        scorefxn(pose)
        #Setup common min elements i.e. jump=True and always alllow bb and sidechain min of strand.
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        mm = pyrosetta.rosetta.core.kinematics.MoveMap()
        mm.set_bb_true_range(pose.conformation().chain_begin(2),pose.conformation().chain_end(2))
        mm.set_jump(True)
        # make pack task according to flags
        try:
            if args.fix_target:
                tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), chain_A, False ))
                tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), chain_B, False ))
            elif args.repack_target_residues:
                repack_residues = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
                for i in args.repack_target_residues.split(','):
                    repack_residues.append_index(int(i))
                chB_and_repackres = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(chain_B,repack_residues)
                chA_fix = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(chB_and_repackres)
                tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), chA_fix, False ))
                tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), chain_B, False ))
                tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), repack_residues, False ))
            else:
                tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), fixnotinterface, False ))
                tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), interface, False ))
        except:
            pass

        # apply settings and start min
        task = tf.create_task_and_apply_taskoperations( pose )
        fast_relax.set_movemap(mm)
        fast_relax.set_task_factory(tf)
        fast_relax.apply(pose)

        # Check whether hbond_lr_bb of complex is better than target alone by certain amount (need to find objective number)
        # If it is dump the pdb after setting denovo strand to chain A and target to chain B
        hbond_lr_bb_complex = scorefxn.score_by_scoretype(pose, pyrosetta.rosetta.core.scoring.hbond_lr_bb)
        if hbond_lr_bb_complex-hbond_lr_bb_target <= float(args.threshold):
            #Switch chain order
            switch = pyrosetta.rosetta.protocols.simple_moves.SwitchChainOrderMover()
            switch.chain_order("21")
            switch.apply(pose)
            pyrosetta.dump_pdb(pose, prefix+"_"+args.strand)
            with open(prefix+"_"+args.strand,'a') as f:
                f.write('delta_bb_hbond '+str(hbond_lr_bb_complex-hbond_lr_bb_target))
            #setup directory and extract motif for grafting step
            subprocess.check_output('grep " A " '+prefix+"_"+args.strand+' > motif_'+prefix+"_"+args.strand,shell=True)
            subprocess.check_output('mkdir '+prefix+'_'+args.strand[:-4]+'; mv '+prefix+'_'+args.strand+' '+prefix+'_'+args.strand[:-4]+'; mv motif_'+prefix+'_'+args.strand+' '+prefix+'_'+args.strand[:-4],shell=True)
            subprocess.check_output('echo '+str(hbond_lr_bb_complex-hbond_lr_bb_target)+'\t'+prefix+'_'+args.strand+' >> hbond_score.list',shell=True)

def n1_cut():
    pose.conformation().delete_residue_slow(pose.conformation().chain_begin(2))

def c1_cut():
    pose.conformation().delete_residue_slow(len(pose.sequence()))

def n2_cut():
    pose.conformation().delete_residue_range_slow(pose.conformation().chain_begin(2),pose.conformation().chain_begin(2)+1)

def n1c1_cut():
    pose.conformation().delete_residue_slow(pose.conformation().chain_begin(2))
    pose.conformation().delete_residue_slow(len(pose.sequence()))

def c2_cut():
    pose.conformation().delete_residue_range_slow(pose.conformation().chain_end(2)-1,pose.conformation().chain_end(2))

pyrosetta.init()
target_pdb = pyrosetta.io.pose_from_file(args.s)
working_pose = pyrosetta.io.pose_from_file(args.s)
strand = pyrosetta.io.pose_from_file(args.strand)
pose = pyrosetta.Pose()
sfx_rep = pyrosetta.ScoreFunction()
sfx_rep.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_rep,0.55)

termini_cut = {'n1': [n1_cut],
        'c1': [c1_cut],
    'n2': [n2_cut],
    'n1c1': [n1c1_cut],
    'c2': [c2_cut]}

# For 3 residue strands, seperate dict, need to find better solution
termini_cut_mod = {'n1': [n1_cut],
        'c1': [c1_cut]}

subprocess.check_output('touch hbond_score.list',shell=True)

if args.strand.startswith('aln_'):
    endresi = subprocess.check_output('tail -1 '+args.strand,shell=True)
    endresi = int(endresi.split()[5])
    startresi = subprocess.check_output('head -1 '+args.strand,shell=True)
    startresi = int(startresi.split()[5])
    length = endresi - startresi + 1
    prefix='rlx_aln_untrimmed_'
    working_pose.append_pose_by_jump(strand,1)
    pose = working_pose.clone()
    # dont continue if the clash is severe
    if sfx_rep(pose) >= 5000:
        print('*** Stopped! Severe clash found!! ***')
        sys.exit(1)
    fastrelax(target_pdb,pose,args,prefix)
    if length > 3:
        for i in termini_cut:
            pose = working_pose.clone()
            prefix = 'rlx_aln_'+i
            for j in termini_cut[i]:
                j()
            fastrelax(target_pdb,pose,args,prefix)
    else:
        for i in termini_cut_mod:
            pose = working_pose.clone()
            prefix = 'rlx_aln_'+i
            for j in termini_cut_mod[i]:
                j()
            fastrelax(target_pdb,pose,args,prefix)    
