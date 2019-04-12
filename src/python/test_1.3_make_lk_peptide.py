#!/usr/bin/env python
# @file: make_LK_peptide.py
# @brief: Make an LK peptide with the sequence pattern "LKKLLKKL..."
# @author: ralford

import sys, os
from pyrosetta import *
from rosetta.core.conformation import *
from rosetta.core.chemical import *

from rosetta import *

init()

seq = "LKKLLKLLKKLLKLLKKLLKLLKKL"
pose = pose_from_sequence( seq );

for i in range( 1, pose.total_residue()+1 ):
    pose.set_phi( i, -57.8 )
    pose.set_psi( i, -47.0 )
    pose.set_omega( i, 180.0 )

pose.dump_pdb( "LK_peptide_n6.pdb" )
