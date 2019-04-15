#!/usr/bin/env python
# @file: make_LK_peptide.py
# @brief: Make an LK peptide with the sequence pattern "LKKLLKKL..."
# @author: ralford

import sys, os
from pyrosetta import *
from rosetta import *
init()

nrepeat = 6
seq = "AAAW"
for i in range(1, nrepeat):
    seq = seq + "AAA"
seq = seq + "WA"

pose = pose_from_sequence( seq );

for i in range( 1, pose.total_residue()+1 ):
    pose.set_phi( i, -57.8 )
    pose.set_psi( i, -47.0 )
    pose.set_omega( i, 180.0 )

# Repack residues
scorefxn = get_fa_scorefxn()
task_pack = standard_packer_task(pose)
task_pack.restrict_to_repacking()
pack = protocols.simple_moves.PackRotamersMover( scorefxn, task_pack )
pack.apply(pose)

pose.dump_pdb( "LWA_peptide_n4.pdb" )
