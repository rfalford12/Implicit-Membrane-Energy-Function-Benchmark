#!/usr/bin/env python
# @file: download_all_pdb_chains.py
# @brief: Download all pdb chains
# @author: Rebecca Faye Alford (ralford3@jhu.edu)

import sys, os

with open( 'docking_bound_set.dat', 'rb' ) as f:
    content = f.readlines()
    content = [ x.strip() for x in content ]

    for i in range(1, len(content)):
        entry = content[i].split('\t')
        pdbcode = entry[0][0:4].upper()
        chains = entry[0][5:7]
        pdbfile = str(pdbcode) + ".pdb"
        if ( not os.path.exists( pdbfile ) ):
            #wget_command = "./get_transformed_pdb.pl " + str(pdbcode) 
            #os.system( wget_command )

            #spanfile_command = "~/apps/Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -in:ignore_unrecognized_res -in:file:s " + pdbcode + "_" + chains + "/" + pdbcode + "_" + chains + ".pdb"
            #os.system( spanfile_command ) 

            # Clean the PDB and download by chains
            #clean_pdb = "python ~/apps/Rosetta/tools/protein_tools/scripts/clean_pdb.py " + pdbcode + " " + chains
            #os.system( clean_pdb )
            #os.system( "mkdir " + pdbcode + "_" + chains )
            #os.system( "mv " + pdbcode + "* " + pdbcode + "_" + chains )
            
            condor_file = pdbcode + "_" + chains + "/" + pdbcode + "_" + chains + ".condor" 
            #os.system( "condor_submit " + condor_file ) 
            entry_name = pdbcode + "_" + chains 
            cleaned_pdb = pdbcode + "_" + chains + "/" + pdbcode + "_" + chains + ".pdb"
            cleaned_spanfile = pdbcode + "_" + chains + "/" + pdbcode + "_" + chains + ".span"
            with open( condor_file, 'w' ) as f: 
                f.write( "#!/bin/bash\n" )
                f.write( "# Automatically generated condor submission file for membrane benchmark\n" )
                f.write( "Universe = vanilla\n" )
                f.write( "output = " + entry_name + "/" + entry_name + ".out\n" )
                f.write( "error = " + entry_name + "/" + entry_name + ".error\n" )
                f.write( "notify_user = rfalford12@gmail.com\n" )
                f.write( "Executable = /home/ralford/apps/Rosetta-stable/main/source/bin/docking_prepack_protocol.linuxgccdebug\n" )  
                f.write( "Arguments = -in:file:s " + cleaned_pdb + " -nstruct 10 -score:weights ref2015.wts -mp:setup:spanfiles " + cleaned_spanfile + " -packing:pack_missing_sidechains 0 -docking:partners " + chains[0] + "_" + chains[1] + " -out:path:all " + entry_name + "\n" )
                f.write( "Queue 10\n" )
            os.system( "condor_submit " + condor_file )

