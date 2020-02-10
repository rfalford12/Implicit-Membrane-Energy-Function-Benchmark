#!/usr/bin/env python
# @file: download_all_pdb_chains.py
# @brief: Download all pdb chains
# @author: Rebecca Faye Alford (ralford3@jhu.edu)

import sys, os

with open( 'pdblist_2018.txt', 'rb' ) as f:
    content = f.readlines()
content = [ x.strip() for x in content ]

for i in range(1, len(content)):
    entry = content[i].split('\t')
    pdbcode = entry[0][0:4].upper()
    chain = entry[0][4:5]
    pdbfile = str(pdbcode) + ".pdb"
    if ( not os.path.exists( pdbfile ) ):
#        wget_command = "./get_transformed_pdb.pl " + str(pdbcode) 
#        os.system( wget_command )

        spanfile_command = "~/apps/Rosetta/main/source/bin/spanfile_from_pdb.macosclangdebug -in:ignore_unrecognized_res -in:file:s " + pdbcode + "_tr.pdb"
        os.system( spanfile_command ) 

        # Clean the PDB and download by chains
    #clean_pdb = "python clean_pdb.py " + pdbfile + " " + chain
    #os.system( clean_pdb )
