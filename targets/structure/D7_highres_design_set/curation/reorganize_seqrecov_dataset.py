#!/usr/bin/env python
# @file: reorganize_seqrecov_dataset.py
# @brief: Reorganize the sequence recovery dataset into folders, etc. 
# @author: Rebecca Faye Alford (ralford3@jhu.edu)

import sys, os

with open( 'natives.list', 'rb' ) as f:
    content = f.readlines()
    content = [ x.strip() for x in content ]

for i in range(1, len(content)):
    entry = content[i].split('\t')
    pdbcode = entry[0][0:4].upper()
    pdbfile = str(pdbcode) + "_tr_ignorechain.pdb"

    # Make directory
   # mkdir_command = "mkdir " + pdbcode
   # os.system( mkdir_command )

    # Move pdbfile 
    mv_pdb_command = "mv " + pdbfile + " " + pdbcode + "/"
    os.system( mv_pdb_command )

    # Move spanfile
    spanfile_command = "mv " + pdbcode + "_tr_ignorechain.span " + pdbcode + "/" + pdbcode + "_tr.span"
    os.system( spanfile_command )

    # Move alternatively named spanfile
    spanfile_command = "mv " + pdbcode + "_tr_ignorechai.span " + pdbcode + "/" + pdbcode + "_tr.span" 
    os.system( spanfile_command ) 


