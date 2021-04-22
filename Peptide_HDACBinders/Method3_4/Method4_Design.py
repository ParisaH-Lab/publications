import argparse
import glob
import math
import os
import sys
import itertools
from itertools import chain
from typing import Union, List

from pyrosetta import Pose, get_score_function, pose_from_file
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.protocols.protein_interface_design.filters import (
     HbondsToResidueFilter,
     ShapeComplementarityFilter,
)
from pyrosetta.rosetta.core.select.residue_selector import (
     ResidueIndexSelector,
     ChainSelector,
     AndResidueSelector,
     NeighborhoodResidueSelector,
     ResiduePDBInfoHasLabelSelector,
)
from pyrosetta.rosetta.protocols.hbnet import UnsatSelector
from pyrosetta.rosetta.protocols.simple_ddg import DdgFilter
from pyrosetta.rosetta.core.chemical.ChemicalManager import get_instance
from pyrosetta.rosetta.core.conformation.ResidueFactory import create_residue
from pyrosetta.rosetta.protocols.cyclic_peptide import PeptideCyclizeMover, DeclareBond
from pyrosetta.rosetta.protocols.moves import NullMover, MoverStatus
from pyrosetta.rosetta.protocols.rosetta_scripts import ParsedProtocol
from pyrosetta.rosetta.protocols.generalized_kinematic_closure import GeneralizedKIC
from pyrosetta.rosetta.core.chemical import (
     UPPER_TERMINUS_VARIANT,
     LOWER_TERMINUS_VARIANT,
     CUTPOINT_LOWER,
     CUTPOINT_UPPER,
)
from pyrosetta.rosetta.core.pose import remove_variant_type_from_pose_residue
from pyrosetta.rosetta.protocols.filters import CombinedFilter
from pyrosetta.rosetta.core.scoring import ScoreType


#global variable set to turn on debug mode
_DEBUG = True

def gen_kic_mover(
    p: Pose,
    res_to_n: int,
    res_to_c: int,
    n_term: int,
    c_term: int,
    scorefxn: ,
    flank: bool
    ) -> bool:
    """Define requirements for genKIC and applies to the pose
    with pre-defined user metrics.
        
    Args:
        p: pose to apply generalized kic to
        res_to_n: number of residues added to N-terminus
        res_to_c: number of residues added to C-terminus
        n_term: number of N-terminal
        c_term: number of c_terminal residue
        scorefxn: score function you want to run with
        flank: should the last residues at C, N be part of the loop
        
    Returns:
        Whether closure was successful or not
    """
    # Defining 3 pivot residues
    pivot_res = list(chain.from_iterable(
                                         [(p.size() - (res_to_c + flank) + 1,'CA'),
                                          (n_term, 'CA'),
                                          (res_to_n + flank , 'CA')]
                                         )
                     )
    # defines the regions to close
    terminus_ranges = ((1, res_to_n + flank + 1),
                        (p.size() - (res_to_c + flank) + 1,
                        p.size() + 1)
                    )
    # defines loop to close
    loop=[]
    for terminus_range in reversed(terminus_ranges):
        for res_num in range(*terminus_range):
            loop.append(res_num)

    # adds hbond based preselction mover
    hbonds_to_residue_filters = [HbondsToResidueFilter(resnum=i,
                                                   partners=0,
                                                   energy_cutoff=-0.25,
                                                   backbone=True,
                                                   sidechain=False,
                                                   bb_bb=True,
                                                   from_other_chains=False,
                                                   from_same_chain=True) for i in loop]
    combined_hbond_filter = CombinedFilter()
    t=-1*(len(loop)/3) # hbond threshold. Can be different.
    combined_hbond_filter.set_threshold(t)
    for filter_ in hbonds_to_residue_filters:
        combined_hbond_filter.add_filter(filter_, -0.5)
    null_mover = NullMover()
    parsed_protocol = ParsedProtocol()
    parsed_protocol.add_mover_filter_pair(null_mover, 'preselection_mover', combined_hbond_filter)

    # sets up the mover
    gk = GeneralizedKIC()
    gk.set_selector_type('lowest_energy_selector')
    gk.set_selector_scorefunction(scorefxn)
    gk.set_closure_attempts(int(1E3))
    gk.set_min_solution_count(10)
    
    gk.add_perturber('randomize_alpha_backbone_by_rama')
    gk.set_perturber_custom_rama_table('flat_symm_dl_aa_ramatable')
    
    gk.set_preselection_mover(parsed_protocol)
    
    for terminus_range in reversed(terminus_ranges):
        for res_num in range(*terminus_range):
            gk.add_loop_residue(res_num)
            if res_num not in (n_term, c_term):
                # add residue to perturber
                gk.add_residue_to_perturber_residue_list(res_num)

    for res_num in pivot_res:
        if type(res_num) != int or res_num in (n_cys, c_cys):
            continue
        # The pivot residues are not necessarily in good regions of
        # Ramachandran space, so we should filter by the rama_prepro energy
        # of pivot positions.
        gk.add_filter('rama_prepro_check')
        gk.set_filter_resnum(res_num)
        gk.set_filter_rama_cutoff_energy(2.0)

    gk.add_filter('loop_bump_check')
    gk.close_bond(n_term, 'N',
                  c_term, 'C',
                  0, '', 0, '',  # optional params -- use default values
                  1.32,
                  114,
                  123,
                  180.,
                  False,
                  False)

    gk.set_pivot_atoms(*pivot_res)
    # applies the mover
    gk.apply(p)
    # returns the status
    return gk.get_last_move_status()


def remove_variants(p: Pose):
    """Remove variant types from the pose
        necessary for bond closure and append/prepend residues.
        
    Args:
        p: pose to perform function on.
    """
    for ir in range(1,p.size()+1):
        if ( p.residue(ir).has_variant_type(UPPER_TERMINUS_VARIANT)):
            remove_variant_type_from_pose_residue( p, UPPER_TERMINUS_VARIANT, ir)
        if ( p.residue(ir).has_variant_type(LOWER_TERMINUS_VARIANT)):
            remove_variant_type_from_pose_residue( p, LOWER_TERMINUS_VARIANT, ir)
        if ( p.residue(ir).has_variant_type(CUTPOINT_LOWER)):
            remove_variant_type_from_pose_residue( p, CUTPOINT_LOWER, ir)
        if ( p.residue(ir).has_variant_type(CUTPOINT_UPPER)):
            remove_variant_type_from_pose_residue( p, CUTPOINT_UPPER, ir)


#this function is basically a filter that counts number of hbonds with some criteria given by the user
def hbond_checker(
    p: Pose,
    chain: int,
    t: float,
    consider_other: bool,
    side_chain: bool,
    consider_same: bool,
    mode: str,
    ) -> Union(bool, float):
    """Count the number of hbonds and returns the number
    or whether it passes the threshold based on mode.
    
    Args:
        p: pose
        chain: what chain we're looking at
        t: threshold of hbonds
        consider_other: consider H-bond from other chains?
        side_chain: Include H-bonds from sidechains?
        consider_same: Consider H-bonds from same chain?
        mode: check or count
        
    Returns:
        True/False if mode is check
        Number of hbonds if mode is count
    """
    loop=[]
    for resNo in range(1,p.size()+1):
        if (p.residue(resNo).chain() == chain):
            loop.append(resNo)
    hbonds_to_residue_filters = []
    number_h=0
    for i in loop:
        indeces=[]
        for j in loop:
            if (abs(j-i) != 2 and abs(j+i) != len(loop)):
                indeces.append(j)
        myList = ','.join(map(str, indeces))
        index_sel=ResidueIndexSelector(myList)

        myFilt=HbondsToResidueFilter()
        if (not consider_other):
            myFilt.set_selector(index_sel)
        myFilt.set_backbone(True)
        myFilt.set_resnum(i)
        myFilt.set_sidechain(side_chain)
        myFilt.set_energy_cutoff(-0.5)
        myFilt.set_from_other_chains(consider_other)
        myFilt.set_from_same_chain(consider_same)
        myFilt.set_partners(0)
        number_h += myFilt.compute(p,i)
        hbonds_to_residue_filters.append(myFilt)

    combined_hbond_filter = CombinedFilter()
    combined_hbond_filter.set_threshold(t)
    for filter_ in hbonds_to_residue_filters:
        combined_hbond_filter.add_filter(filter_, -0.5)
    result = combined_hbond_filter.apply(p)
    if (mode == "check"):
        if (_DEBUG and not(result)):
            print ("the number of hbods",number_h, "is not enough to satisfy")
        return result
    if (mode == "count"):
        if (_DEBUG and not(result)):
            print ("the number of hbods is",number_h)
        return number_h


def bonder(p: Pose, n: int, c: int):
    """Declare amino acid bond.
    Args:
        p: pose
        n: residue to add bond at N-terminus
        c: residue to add bond at C-terminus
    """
    dbm = DeclareBond()
    dbm.set(n, 'N', c, 'C', False)
    dbm.apply(p)


def cyclizer(p: Pose, close_list: List):
    """Close bond and apply cyclic bond constraints.
    
    Args:
        p: pose
        close_list: [n_term, c_term]
    """
    to_close = ','.join(map(str, close_list))
    index_sel = ResidueIndexSelector(to_close)
    pcm = PeptideCyclizeMover()
    pcm.set_selector(index_sel)
    pcm.apply(p)


def relax(p: Pose,chain: int):
    """Relax chain of pose."""
        
    i=0
    indeces=[]
    for resNo in range(1,pose.size()+1):
        if (pose.residue(resNo).chain() == chain):
            indeces.append(resNo)
            if (i ==0):
                nb=resNo
            i=i+1
            if (pose.residue(resNo).name3() == "CYS" or pose.residue(resNo).name3() == "SHA"):
                ce=resNo
    cyclizer(pose,indeces)
    my_score=get_score_function()
    my_score.set_weight(ScoreType.atom_pair_constraint, 1.0)
    my_score.set_weight(ScoreType.angle_constraint, 1.0)
    my_score.set_weight(ScoreType.dihedral_constraint, 1.0)
    my_score.set_weight(ScoreType.hbond_sr_bb, 10.0)
    my_score.set_weight(ScoreType.hbond_lr_bb, 10.0)
    frm = pyrosetta.rosetta.protocols.relax.FastRelax(my_score)
    frm.ramp_down_constraints(False)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_jump(1,False)
    mm.set_bb_true_range(nb, ce-1)
    mm.set_bb_true_range(ce+1, nb+i-1)
    mm.set_chi_true_range(nb, ce-1)
    mm.set_chi_true_range(ce+1, nb+i-1)
#I have to stop SHA from moving or it comes out of pocket
#mm.set_bb_true_range(nb, nb+i-1)
#mm.set_chi_true_range(nb, nb+i-1)
    frm.set_movemap(mm)
    frm.apply(pose)
    bonder(pose,nb,nb+i-1)
    pyrosetta.rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(pose, 1)
    if (_DEBUG):
        pose.dump_pdb('relax.pdb')


def satisfy(pose: Pose,chain: int) -> bool:
    """Check unsatified hbond donor/acceptor at interface.
        
    Returns:
        True if few unsatified hbond donor/acceptor are present
    """
    checker=0
    indeces=[]
    for resNo in range(1,pose.size()+1):
        if (pose.residue(resNo).chain() == chain):
            indeces.append(resNo)
    to_check = ','.join(map(str, indeces))
    index_sel=ResidueIndexSelector(to_check)
    unsats=UnsatSelector()
    unsats.set_consider_mainchain_only(False)
    unsats.set_scorefxn(pyrosetta.get_score_function())
    unsat_chain=AndResidueSelector(index_sel,unsats)
    subset=unsat_chain.apply(pose)
    for i in range(1,pose.size()+1):
        nearby=[]
        if (subset[i]):
            r1=pyrosetta.rosetta.core.conformation.Residue(pose.residue(i))
            for j in range(1,pose.size()+1):
                if (abs(j-i) > 1):
                    #r2=pyrosetta.rosetta.core.conformation.Residue(pose.residue(j))
                    for at in range (1,pose.residue(j).natoms()):
                        d_sq=r1.xyz('O').distance_squared(pose.xyz(core.id.AtomID(at,j)))
                        if (d_sq < 25.01):
                            nearby.append([i,j])
        if (len(nearby) > 30):
            if _DEBUG:
                print ("you have an unsatisfied carbonyl here, so better not continue on this")
            checker=checker+1
            break
    
    if (checker == 0):
        return True
    else:
        return False


def qualifier(pose: Pose,num: int):
    """Define neighborhood properties of a given residue."""

    polar=0
    hydrophobic=0
    r1=pyrosetta.rosetta.core.conformation.Residue(pose.residue(num))
    for i in range(1,pose.size()+1):
        if (pose.residue(i).chain() == 2):
            r2=pyrosetta.rosetta.core.conformation.Residue(pose.residue(i))
            if (r2.name3() == "GLY"):
                d_sq_app=r1.xyz('CA').distance_squared(r2.xyz('CA'))
            else:
                d_sq_app=r1.xyz('CA').distance_squared(r2.xyz('CB'))
            if (d_sq_app < 100.01):
                for at in range (1,pose.residue(i).natoms()):
                    d_sq=r1.xyz('CA').distance_squared(pose.xyz(core.id.AtomID(at,i)))
                    if (d_sq < 64.01):
                        if (("O" in r2.atom_name(at)) or ("N" in r2.atom_name(at)) or (r2.atom_type(at).is_polar_hydrogen())):
                            polar+=1
                        elif (("C" in r2.atom_name(at)) or ("S" in r2.atom_name(at))): #ignoring nonpolar hydrogens
                            hydrophobic+=1
    if _DEBUG:
        print ("number of polar atoms is",polar," and number of hydrophobic atoms is",hydrophobic)

    if (polar > hydrophobic):
        pose.pdb_info().add_reslabel(num,"POLAR")
    elif (polar < hydrophobic):
        pose.pdb_info().add_reslabel(num,"HYDROPH")
    else:
        pose.pdb_info().add_reslabel(num,"BOTH")

    if _DEBUG:
        pose.dump_pdb('modified.pdb')


def counter(pose: Pose,remark: str) -> int:
    """Checks residues that are labeled with a certain remark."""

    count=0
    mySel=ResiduePDBInfoHasLabelSelector(remark)
    SelSet=mySel.apply(pose)
    for i in range(1,pose.size()+1):
        if (SelSet[i]):
            count+=1

    return count


def designer(pose):
    """Perform design."""
    
    #step0: initial step to get rid of all the Gly:(no check to accept. all Gly are removed here. phew)
    normal_design1= protocols.rosetta_scripts.XmlObjects.create_from_file("design_normal.xml")
    whole_normal1= normal_design1.get_mover("ParsedProtocol")
    whole_normal1.apply(pose)

    #step1: running the polars
    polar_design= protocols.rosetta_scripts.XmlObjects.create_from_file("design_polar.xml")
    whole_polar= polar_design.get_mover("ParsedProtocol")
    #defining number of residues around polar areas
    count=counter(pose,"POLAR")
    #minimum hbonds
    hbond_num=-1*(count/2)
    #first round of design for getting hbond across the interface
    whole_polar.apply(pose)
    p_curr=pose.clone()
    p_old=pose.clone()
    old_hbond=hbond_checker(p_old,1,10,1,1,"count",0)
    for i in range(5):
        whole_polar.apply(p_curr)#p_curr is being modified here
        if (_DEBUG):
            print ("the old hbond is", old_hbond, "and new value is",hbond_checker(p_curr,1,10,1,1,"count",0) )
        if (hbond_checker(p_curr,1,10,1,1,"count",0) > old_hbond):#if hbonds are getting better. hbond count is positive
            old_hbond=hbond_checker(p_curr,1,10,1,1,"count",0)
            p_old=p_curr.clone() #a history is being kept of the good p_curr
        else:#if hbonds remained about the same, keep the old pose
            p_curr=p_old.clone()
#at the end of this cycle, p_old and p_curr are basically same things

    #step2: designing the hydrophobic area
    hydro_design= protocols.rosetta_scripts.XmlObjects.create_from_file("design_hydrophobic.xml")
    whole_hydro= hydro_design.get_mover("ParsedProtocol")
    #first round of design and getting initial shape complementarity value
    whole_hydro.apply(pose)
    sc=ShapeComplementarityFilter()
    sc.jump_id(1)
    old_sc=sc.report_sm(p_old)
    for i in range(5):
        whole_hydro.apply(p_curr)#p_curr is being modified here
        if (_DEBUG):
            print ("the old sc is", old_sc, "and new value is",sc.report_sm(p_curr) )
        if (sc.report_sm(p_curr) > old_sc):#if sc is getting better. Remember, good sc is more positive
            old_sc=sc.report_sm(p_curr)
            p_old=p_curr.clone() #a history is being kept of the good p_curr
        else:#if sc remained about the same or worse, keep the old pose
            p_curr=p_old.clone()

    #at the end of this cycle, p_old and p_curr are basically same things
    #step3: final design
    normal_design= protocols.rosetta_scripts.XmlObjects.create_from_file("design_normal.xml")
    whole_normal= normal_design.get_mover("ParsedProtocol")
    #first round of design and getting initial shape complementarity value
    whole_normal.apply(pose)
    ddg=DdgFilter(0,pyrosetta.get_score_function(),1)
    ddg.repack(0)
    old_ddg=ddg.report_sm(p_old)
    for i in range(3):
        whole_normal.apply(p_curr)#p_curr is being modified here
        if (_DEBUG):
            print ("the old ddg is", old_ddg, "and new value is",ddg.report_sm(p_curr) )
        if (ddg.report_sm(p_curr) < old_ddg):#if ddg is getting better
            old_ddg=ddg.report_sm(p_curr)
            p_old=p_curr.clone() #a history is being kept of the good p_curr
        else:#if ddg remained about the same or worse, keep the old pose
            p_curr=p_old.clone()

    return p_curr

##############################################################################################
########################## MAIN FUNCTION THAT CALLS EVERYTHING ELSE ##########################
##############################################################################################
def main(argv):
    
    parser = argparse.ArgumentParser(description='Program')
    parser.add_argument('-i', '--input', action='store', type=str,
                    required=True,
                    help='input target pdb')
    parser.add_argument('-l', '--length', action='store', type=int,
                    required=True,
                    help='max length of the loop, here peptide')
    parser.add_argument('-m', '--min', action='store', type=int,
                    default=4,
                    help='min length of the loop, be aware of terminal inclusion')
    parser.add_argument('-r', '--resn', action='store', type=str,
                    default='GLY',
                    help='residue type to append')
    parser.add_argument('-n', '--nstruct', action='store', type=int,
                    default=1,
                    help='how many times to run each KIC run')
    parser.add_argument('-c', '--chain', action='store', type=int,
                        default=1,
                        help='what chain number is the one I am extending')

    parser.add_argument('-t', '--testnum', action='store', type=str,
                        default='t1',
                        help='name of directory')
    args = parser.parse_args()

    #adding additional params files of interest
    params = glob.glob('params/*.params')
    #initiaiting Rosetta
    init(extra_options='-in:file:fullatom true -mute all -write_all_connect_info -extra_res_fa SHA.params -ignore_waters false -auto_setup_metals true -extra_res_fa {}'.format(' '.join(params)))
    scorefxn = get_score_function()
    
    # get the pose from target and scaffold
    p_in=rosetta.core.import_pose.pose_from_file(args.input)

    in_fname = args.input
    include_initial_termini_in_loop = True #changed this for purpose of new model 
    base_fname = in_fname.split('.')[0]
    out_put_fname = '{dir}/{base}_N-{n_add}_C-{c_add}_{numrun}.pdb'
    
    if (args.length + int(include_initial_termini_in_loop) < 3):
        print('the loop needs to be at least 3 residues', file=sys.stderr)
        sys.exit(1)

    #defining the residue I want to append
    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set( 'fa_standard' )
    res = rosetta.core.conformation.ResidueFactory.create_residue(rts.name_map(args.resn))

    for loop in range (args.min-int(include_initial_termini_in_loop),args.length+1):
        if (_CHECK):
            break
        #This whole part bewlo is just to find a closed solution
        out_dir = '{}_genKIC_{}'.format(args.testnum,loop)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        for runs in range (0,args.nstruct):
            p=Pose()
            for resNo in range(1,p_in.size()+1):
                if (p_in.residue(resNo).chain() == args.chain):
                    p.append_residue_by_bond(p_in.residue(resNo),False)
            num=p.size()+1
            for ir in range(1,p.size()+1):
                if ( p.residue(ir).has_variant_type(core.chemical.UPPER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.  UPPER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.LOWER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.LOWER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_LOWER)):
                        core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_LOWER, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_UPPER)):
                        core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_UPPER, ir)
            if _DEBUG:
                p.dump_pdb('pose_clone.pdb')
        
            # setting residues to prepend to N-term and setting omega to 180
            N_add=loop//2
            for i in range(0, N_add):
                p.prepend_polymer_residue_before_seqpos(res, 1, True)
            if _DEBUG:
                p.dump_pdb('prepend_test.pdb')
            for res_no in range(1,N_add+1):
                p.set_omega(res_no,180.)
            # setting residues to append to C-term and setting omega to 180
            C_add=loop-N_add
            for ir in range(1,p.size()+1):
                if ( p.residue(ir).has_variant_type(core.chemical.UPPER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.UPPER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.LOWER_TERMINUS_VARIANT)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.LOWER_TERMINUS_VARIANT, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_LOWER)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_LOWER, ir)
                if ( p.residue(ir).has_variant_type(core.chemical.CUTPOINT_UPPER)):
                    core.pose.remove_variant_type_from_pose_residue( p, core.chemical.CUTPOINT_UPPER, ir)
            for i in range(0,C_add):
                p.append_residue_by_bond(res, True)
            if _DEBUG:
                p.dump_pdb('append_test.pdb')
            for res_no in range((p.size()-C_add)-1, p.size()+1):
                p.set_omega(res_no,180.)
            # declaring the bond before genKIC call
            to_close=(1,p.size())
            pcm=protocols.cyclic_peptide.PeptideCyclizeMover()
            pcm.apply(p)
            if _DEBUG:
                p.dump_pdb('bonded.pdb')
            # calling genKIC mover
            st=gen_kic_mover(p,N_add,C_add,to_close[0],to_close[1],scorefxn,int(include_initial_termini_in_loop))
            if st == protocols.moves.MoverStatus.MS_SUCCESS:
                p_fin=Pose()
                for resi in range(1,p.size()+1):
                    p_fin.append_residue_by_bond(p.residue(resi),False)
                old_chain=2 #assumption is peptide is always chain 1.
                for resNo in range(1,p_in.size()+1):
                    new_chain=p_in.residue(resNo).chain()
                    if (p_in.residue(resNo).chain() != args.chain):
                        if (resNo == num):
                            p_fin.append_residue_by_jump(p_in.residue(resNo),p.size(),'','',True)
                        elif (new_chain != old_chain):
                            p_fin.append_residue_by_jump(p_in.residue(resNo),p.size(),'','',True)
                            old_chain=new_chain
                        else:
                            p_fin.append_residue_by_bond(p_in.residue(resNo),False)
                db = protocols.cyclic_peptide.DeclareBond()
                db.set(to_close[0],'N',to_close[1],'C',False,False,0,0,True)
                length=(to_close[1]-to_close[0])+1
                hbond_lim=-1*(length/3)
                if (_DEBUG):
                    print ("size of your peptide is", length, "and hbond limit is", hbond_lim)
                db.apply(p_fin)
                if (_DEBUG) :
                    p_fin.dump_pdb('genKIC.pdb')
                checker=0
                #from this moment we will be focusing on the other parts (not just the backbone sample)
                #Step1: Backbone shear or relax with high hbond score to see if I can help with the hbonds
                relax(p_fin,args.chain)
                #Step2: are all the carbonlys that are buried satisified?
                if _DEBUG:
                    print ("----Initial backbone unsat check-----")
                if (satisfy(p_fin,args.chain)):
                  #Step3: Check if we have enough hbonds internally
                    if _DEBUG:
                        print ("----Initial backbone hbond check-----")
                    if (hbond_checker(p_fin,args.chain,hbond_lim,0,0,"check",1)):
                #if (hbond_checker(p_fin,args.chain,2,0,0,"check",1)):
                        checker=1
                    #if either 2 and 3 are correct, move on with design
                    if (checker == 1):
                        #Step5: for each position, check where its C-a/C-b vector is pointing and what kind of site it is in
                        #based on the surrounding, tag it with a group (polar, hydrophobic, charged), add to PDB REMARKs
                        hacky_counter=0 #this is just to speed up because I know peptide is in chain A
                        p_fin.dump_pdb('essential.pdb')
                        p_c=rosetta.core.import_pose.pose_from_file('essential.pdb')
                        for number in range(1,p_c.size()+1):
                            if (p_c.residue(number).chain() == args.chain):
                                if (_DEBUG):
                                    print ("residue number under evaluation is", number)
                                qualifier(p_c,number)
                                hacky_counter+=1
                            if (hacky_counter == length): #when chain A is done!
                                break
                        #Step6: design
                        outcome=designer(p_c)
                        #only print if I have some hbonds
                        if _DEBUG:
                            print ("----Final backbone unsat check-----")
                        if (hbond_checker(outcome,args.chain,hbond_lim,0,0,"check",1)):
                            outcome.dump_pdb(out_put_fname.format(dir=out_dir,
                                                          base=base_fname,
                                                          n_add=N_add,
                                                          c_add=C_add,
                                                          numrun=runs+1))

            else:
                print ("no successful solutions found for genKIC")


if __name__ == '__main__':
    main(sys.argv)
