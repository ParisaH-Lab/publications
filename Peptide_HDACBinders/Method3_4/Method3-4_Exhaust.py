import argparse
import glob

import math
import os
import sys
import csv
import itertools
from itertools import chain
import numpy as np
import rmsd
import re

from pyrosetta.rosetta.core.chemical import (
     UPPER_TERMINUS_VARIANT,
     LOWER_TERMINUS_VARIANT,
     CUTPOINT_LOWER,
     CUTPOINT_UPPER,
)
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pose import remove_variant_type_from_pose_residue
from pyrosetta.rosetta.protocols.protein_interface_design.filters import HbondsToResidueFilter
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import NeighborhoodResidueSelector
from pyrosetta.rosetta.protocols.residue_selectors import UnsatSelector
from pyrosetta.rosetta.core.select.residue_selector import ResiduePDBInfoHasLabelSelector
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.protocols.simple_ddg import DdgFilter
from rosetta import *
from pyrosetta.bindings.utility import bind_method
from pyrosetta.rosetta.core.pack.rotamer_set import bb_independent_rotamers
from pyrosetta.rosetta.core.conformation import Residue


#global variable to turn on debug mode
_DEBUG = True
# Defining global lists
L_res=["GLY","ALA","CYS","MET","VAL","LEU","ILE","ASP","GLU","ASN","GLN","THR","SER","TYR","TRP","PHE","LYS","ARG","HIS","PRO"]
D_res=["GLY","DALA","DCYS","DMET","DVAL","DLEU","DILE","DASP","DGLU","DASN","DGLN","DTHR","DSER","DTYR","DTRP","DPHE","DLYS","DARG","DHIS","DPRO"]


def variant_remove(p):
    """Remove variant types from pose."""
    for ir in rrange(p.size()):
        if ( p.residue(ir).has_variant_type(UPPER_TERMINUS_VARIANT)):
            remove_variant_type_from_pose_residue( p,UPPER_TERMINUS_VARIANT, ir)
        if ( p.residue(ir).has_variant_type(LOWER_TERMINUS_VARIANT)):
            remove_variant_type_from_pose_residue( p,LOWER_TERMINUS_VARIANT, ir)
        if ( p.residue(ir).has_variant_type(CUTPOINT_LOWER)):
            remove_variant_type_from_pose_residue( p, CUTPOINT_LOWER, ir)
        if ( p.residue(ir).has_variant_type(CUTPOINT_UPPER)):
            remove_variant_type_from_pose_residue( p, CUTPOINT_UPPER, ir)


def extend(n, p, a, c_, n_):
    """Extend anchr residue.
    Assumes one chain.
    
    Args:
        n: number to add to pose
        p: pose
        a: residue to add
        c_: boolean, whether to add to c_terminal
        n_: boolean, whether to add to n_terminal
    """
    variant_remove(p)

    if (n_):
        for i in range(0, n):
            p.prepend_polymer_residue_before_seqpos(a, 1, True)
        if _DEBUG:
            p.dump_pdb('prepend_test.pdb')
            for res_no in range(1,n+1):
                p.set_omega(res_no,180.)
    if (c_):
        for i in range(0, n):
            p.append_residue_by_bond(a, True)
            if _DEBUG:
                p.dump_pdb('append_test.pdb')
        for res_no in range(1, p.size()+1):
            p.set_omega(res_no,180.)

    if (_DEBUG):
        p.dump_pdb('extended.pdb')


def bin_sample(p,res, phi, psi):
    """Sample phi and psi torsion for residue res in pose p."""

    p.set_phi(res, phi)
    p.set_psi(res, psi)

    if (_DEBUG):
        p.dump_pdb('set_{}_{}.pdb'.format(phi,psi))


def relax(pose,c):
    """Relax chain c of pose using cartesian minimization."""

    indeces=[]
    for resNo in rrange(pose.size()):
        if (pose.residue(resNo).chain() == c):
            indeces.append(resNo)
    my_score=pyrosetta.get_score_function()
    #set cart weigths
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_angle, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_length, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_ring, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_torsion, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_proper, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_improper, 1.0)
    #set metal constraints
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.metalbinding_constraint, 1.0)
    frm = pyrosetta.rosetta.protocols.relax.FastRelax(my_score)
    frm.ramp_down_constraints(False)
    #set cart to be true
    frm.cartesian(True)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_jump(1,False)
    mm.set_chi_true_range(indeces[0], indeces[-1])
    frm.set_movemap(mm)
    frm.apply(pose)
    pyrosetta.rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(pose, 1)
    if (_DEBUG):
        pose.dump_pdb('relax.pdb')


def unsat_count(pose,chain):
    """Very rough function that checks for unsatisfied hbond donors
       or acceptors at interface."""
    
    checker=0
    indeces=[]
    for resNo in rrange(pose.size()):
        if (pose.residue(resNo).chain() == chain):
            indeces.append(resNo)
    to_check = ','.join(map(str, indeces))
    index_sel=ResidueIndexSelector(to_check)
    unsats=UnsatSelector()
    unsats.set_consider_mainchain_only(False)
    unsats.set_scorefxn(pyrosetta.get_score_function())
    unsat_chain=AndResidueSelector(index_sel,unsats)
    subset=unsat_chain.apply(pose)
    all=[]
    for i in rrange(pose.size()):
        nearby=[]
        if (subset[i]):
            r1=pyrosetta.rosetta.core.conformation.Residue(pose.residue(i))
            for j in rrange(pose.size()):
                if (abs(j-i) > 1):
                    for at in range (1,pose.residue(j).natoms()):
                        d_sq=r1.xyz('O').distance_squared(pose.xyz(core.id.AtomID(at,j)))
                        if (d_sq < 25.01):
                            nearby.append([i,j])
        if (len(nearby) > 30):
            all.append(i)
    
    return len(all)


def metrics(p,c):
    """Score the pose and output different interface metric."""

    #you need to score so that Rosetta can give you the scores you want.
    scf=pyrosetta.get_score_function()
    scf(p)
    # getting score of the whole pose
    tot_energy_value=p.energies().total_energies()[core.scoring.score_type_from_name("total_score")]
    # getting score of the peptide chain
    pep_energy_value=0
    for resNo in rrange(p.size()):
        if (p.residue(resNo).chain() == c):
            pep_energy_value+=p.energies().residue_total_energies(p.residue(resNo).seqpos())[core.scoring.score_type_from_name("total_score")]
    # getting shape complemetarity
    sc=ShapeComplementarityFilter()
    sc.jump_id(1)
    sc_value=sc.report_sm(p)
    #getting ddg_norepack
    ddg1=DdgFilter(0,pyrosetta.get_score_function(),1)
    ddg1.repack(0)
    ddg_value1=ddg1.report_sm(p)
    # getting ddg_repack
    ddg2=DdgFilter(0,pyrosetta.get_score_function(),1)
    ddg2.repack(0)
    ddg_value2=ddg2.report_sm(p)
    #buried unsat
    unsat_value=unsat_count(p,c)

    metrics=[("total energy",tot_energy_value),("peptide energy",pep_energy_value),("shape complementarity",sc_value),("ddg no_repack",ddg_value1), ("ddg repack",ddg_value2), ("number unsats",unsat_value)]

    return metrics


def add_to_score(metric,res_inf,res):
    """Write array to csv."""
    
    everything=[]
    for info in res_inf:
        everything.append(info)
    for energy in metric:
        everything.append(energy[1])

    with open(r'scores_{}.csv'.format(res), 'a') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerow(everything)


def coord_find(p,ir,ia):
    """Find coordinate of atom ia in residue ir of pose p."""

    name=re.compile('.HA ')
    if (len(re.findall(name, ia)) == 0):
        coord_xyz=p.xyz(core.id.AtomID(p.residue(ir).atom_index(ia),ir))
        coord_arr=[]
        x=coord_xyz[0]
        y=coord_xyz[1]
        z=coord_xyz[2]
        coord_arr.append(x)
        coord_arr.append(y)
        coord_arr.append(z)
    
        if (_DEBUG):
            print (ia, coord_arr)
    
        return coord_arr


def find_cent(A):
    """Find the center of mass of coordinates A."""
    
    sumA=[0,0,0]
    for i in range(len(A)):
        sumA[0]=sumA[0]+A[i][0]
        sumA[1]=sumA[1]+A[i][1]
        sumA[2]=sumA[2]+A[i][2]
    
    for i in range(3):
        sumA[i]=sumA[i]/len(A)
    
    if (_DEBUG):
        print ("========= find center function==========")
        print ("the elements are", '\n'
               ,A,'\n'
               ,"and the center is:",'\n'
               ,sumA)
    
    return sumA


def transform(orig,rotd,resi, r_end):
    """Apply transformation from rotd to origin.
       The transformation matrix is calculated such that
       r_end and resi overlat.
    """

    p_scaff=[]
    p_targ=[]
    for atom in rrange(orig.residue(resi).natoms()):
        if (not (orig.residue(resi).atom_is_backbone(atom))):
            if (not orig.residue(resi).atom_is_hydrogen(atom)):
                p_scaff.append(coord_find(rotd,r_end,rotd.residue(r_end).atom_name(atom)))
                p_targ.append(coord_find(orig,resi,orig.residue(resi).atom_name(atom)))

    #step1: moving scaffold to the center
    T=find_cent(p_scaff)
    plusv=numeric.xyzVector_double_t(-1*T[0],-1*T[1],-1*T[2])
    #does not rotate
    noR=numeric.xyzMatrix_double_t.cols(1,0,0,0,1,0,0,0,1)
    rotd.apply_transform_Rx_plus_v(noR,plusv)
    if (_DEBUG):
        print ("============scaffold translation, step 1===============")
        print (T,noR)
        rotd.dump_pdb('translate_scaffold.pdb')

    #Step1': get the coordinates of target at the center
    T_targ=find_cent(p_targ)
    v_targ=numeric.xyzVector_double_t(-1*T_targ[0],-1*T_targ[1],-1*T_targ[2])
    orig.apply_transform_Rx_plus_v(noR,v_targ)
    if (_DEBUG):
        print ("============target translation, step 1===============")
        print (T_targ,noR)
        orig.dump_pdb('translate_target.pdb')

    #need to re-load the matrix now because the pose has changed
    p_scaff_new=[]
    p_targ_new=[]
    for atom in rrange(orig.residue(resi).natoms()):
        if (not (orig.residue(resi).atom_is_backbone(atom))):
            if (not orig.residue(resi).atom_is_hydrogen(atom)):
                p_scaff_new.append(coord_find(rotd,r_end,rotd.residue(r_end).atom_name(atom)))
                p_targ_new.append(coord_find(orig,resi,orig.residue(resi).atom_name(atom)))


    #Step 2: get the rotation matrix
    #the magic of libraries
    semi_V=rmsd.kabsch(p_scaff_new,p_targ_new)
    V=np.linalg.inv(semi_V)
    if (_DEBUG):
        print ("============scaffold rotation, step 2===============")
        print ("the transformation matrix is",'\n', V)

    #Rotate the pose
    Rx=numeric.xyzMatrix_double_t.cols(V[0][0],V[1][0],V[2][0],V[0][1],V[1][1],V[2][1],V[0][2],V[1][2],V[2][2])
    noT=numeric.xyzVector_double_t(0,0,0)
    if (_DEBUG):
        print ("the old rmsd is:", rmsd.kabsch_rmsd(p_scaff, p_targ))

    #moving the pose
    rotd.apply_transform_Rx_plus_v(Rx,noT)
    if (_DEBUG):
        rotd.dump_pdb('rotate_scaffold.pdb')

    #Step3: translate the pose back to target (both the new and the original)
    scaff_trans=numeric.xyzVector_double_t(T_targ[0],T_targ[1],T_targ[2])
    rotd.apply_transform_Rx_plus_v(noR,scaff_trans)
    orig.apply_transform_Rx_plus_v(noR,scaff_trans)
    if (_DEBUG):
        rotd.dump_pdb('back_to_orig.pdb')
        print ("8888888", scaff_trans)

    #generating final set
    p_scaff_final=[]
    for atom in rrange(orig.residue(resi).natoms()):
        if (not (orig.residue(resi).atom_is_backbone(atom))):
            if (not orig.residue(resi).atom_is_hydrogen(atom)):
                p_scaff_final.append(coord_find(rotd,r_end,rotd.residue(r_end).atom_name(atom)))
    if (_DEBUG):
        print ("the new rmsd is:", rmsd.kabsch_rmsd(p_scaff_final, p_targ))


def packer (p,resn,resi,c):
    """Perform mutation and minimization."""

    indeces=[]
    for resNo in rrange(p.size()):
        if (p.residue(resNo).chain() == c):
            indeces.append(resNo)
    #step one is to add the residue you want
    mut=protocols.simple_moves.MutateResidue()
    mut.set_res_name(resn)
    mut.set_target(resi)
    mut.set_preserve_atom_coords(True)
    mut.apply(p)

    #step two is to minimize with a cart_min probably
    my_score=pyrosetta.get_score_function()
    #set cart weigths
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_angle, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_length, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_ring, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_torsion, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_proper, 1.0)
    my_score.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded_improper, 1.0)
    min=protocols.minimization_packing.MinMover()
    min.score_function(my_score)
    min.cartesian(True)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_jump(1,False)
    mm.set_bb_true_range(indeces[0], indeces[-1])
    mm.set_chi_true_range(indeces[0], indeces[-1])
    min.set_movemap(mm)
    min.apply(p)


def res_setter(p,resi,resn,phi,psi):
    """Set the torsion and identity of peptides."""

    #mutate the residue
    mut=protocols.simple_moves.MutateResidue()
    mut.set_res_name(resn)
    mut.set_target(resi)
    mut.set_preserve_atom_coords(True)
    mut.apply(p)

    #setting up the torsions
    p.set_phi(resi, float(phi))
    p.set_psi(resi, float(psi))

##############################################################################################
########################## MAIN FUNCTION THAT CALLS EVERYTHING ELSE ##########################
##############################################################################################
def main(argv):
    
    parser = argparse.ArgumentParser(description='Program')
    parser.add_argument('-i', '--input', action='store', type=str,
                    required=True,
                    help='input target pdb')
    parser.add_argument('-n', '--number', action='store', type=int,
                    required=True,
                    help='number of extension at each end')
    parser.add_argument('-m', '--mut', action='store', type=int,
                    required=True,
                    help='residue to run all the samples on')
    parser.add_argument('-a', '--aa', action='store', type=str,
                    default='GLY',
                    help='aa type to append')
    parser.add_argument('-c', '--chain', action='store', type=int,
                    default=1,
                    help='what chain number is the one I am extending')
    parser.add_argument('-r', '--relax', action='store', type=bool,
                    default=False,
                    help='should I run relax')
                        

    args = parser.parse_args()
    #list of additional param files one wants to add.
    params = glob.glob('*.params')
    #initiaiting Rosetta
    init(extra_options='-beta_nov16_cart -in:file:fullatom true -write_all_connect_info -extra_res_fa SHA.params -ignore_waters false -mute all -auto_setup_metals true -extra_res_fa {}'.format(' '.join(params)))
    scorefxn = get_score_function()
    
    # get the pose from target and scaffold
    p_in=rosetta.core.import_pose.pose_from_file(args.input)

    #defining the residue I want to append
    chm = rosetta.core.chemical.ChemicalManager.get_instance()
    rts = chm.residue_type_set( 'fa_standard' )
    res = rosetta.core.conformation.ResidueFactory.create_residue(
                                                rts.name_map(args.aa))
    
    p=Pose()
    for resNo in rrange(p_in.size()):
        if (p_in.residue(resNo).chain() == args.chain):
            p.append_residue_by_bond(p_in.residue(resNo),False)
    if (_DEBUG):
        p.dump_pdb('just_pep.pdb')

    extend(args.number, p, res, True,True)

    #starting SHA
    pep_s=Pose()
    for resNo in rrange(p.size()):
        pep_s.append_residue_by_bond(p.residue(resNo),False)

    if (_DEBUG):
        pep_s.dump_pdb('clone_check.pdb')
    

    with open('scores.csv', 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerow(["residue1_phi_psi","total_score","peptide_total_score","sc","ddg_no_repack","ddg_repack","unsats"])

    # going through all the residues, all changes slowly.
    for phi in range(30,181,10):
        if (phi > -30 and phi < 30):
            continue
        for psi in range(-180,181,10):
            bin_sample(p,args.mut, phi, psi)
            #after sampling making sure SHA is still there
            transform(pep_s,p,3,3)

            #adding things back as one chain
            variant_remove(p_in)
            p_fin=Pose()
            old_chain=2 #assumption is peptide is always chain 1.
            for resi in rrange(p.size()):
                p_fin.append_residue_by_bond(p.residue(resi),False)
                counter=0
            for resNo in rrange(p_in.size()):
                new_chain=p_in.residue(resNo).chain()
                if (p_in.residue(resNo).chain() != args.chain):
                    if (counter == 0):
                        p_fin.append_residue_by_jump(p_in.residue(resNo),
                        .size(),'','',True)
                        counter +=1
                    elif (new_chain != old_chain):
                        p_fin.append_residue_by_jump(p_in.residue(resNo),
                        .size(),'','',True)
                        old_chain=new_chain
                    else:
                        p_fin.append_residue_by_bond(p_in.residue(resNo),False)
            if (_DEBUG):
                p_fin.dump_pdb('moved_{}_{}.pdb'.format(phi,psi))

            if (phi > 0):
                for resNa in D_res:
                    packer(p_fin,resNa, args.mut,args.chain)
                    add_to_score(metrics(p_fin,1),[(args.mut,resNa,phi,psi),
                        (args.mut,resNa, p.phi(args.mut),p.psi(args.mut))],args.mut)
            else:
                for resNa in L_res:
                    packer(p_fin,resNa, args.mut,args.chain)
                    add_to_score(metrics(p_fin,1),[(args.mut,resNa,phi,psi),
                        (args.mut,resNa, p.phi(args.mut),p.psi(args.mut))],args.mut)
    if (args.relax):
        relax(p_fin,args.chain)

    if (_DEBUG):
        p_fin.dump_pdb('final.pdb')

if __name__ == '__main__':
    freeze_support()
    main(sys.argv)
