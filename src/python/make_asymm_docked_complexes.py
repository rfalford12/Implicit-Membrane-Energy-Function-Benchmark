#!/usr/bin/env python
###################################################################
#@file:         make_asymm_docked_complexes.py                                                                                   
#@description:  Assemble membrane protein complexes                                             
#@args:         --energy_fxn (wts) --targets (list of targets)                                                          
#@author:       Rebecca F. Alford                   
#@email:        rfalford12@gmail.com                                          
###################################################################

import random, os, sys
import hpc_util, read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def run_docking_calc( energy_fxn, config, targets, test_name ): 
  """
    A function for assembling symmetric membrane protein complexes

    Arguments: 
      energy_fxn = energy function to use for calculations (typically, name of the weights file)
      config = path to benchmark, Rosetta executables
      targets = list of targets
      test_name = Name of test
  """

  print( "Predicting protein-protein complex structures for set", targets ) 

  # Read list of energy landscape test cases
  list_of_targets = config.benchmark_path + "targets/" + targets + "/targets.list"
  with open( list_of_targets, 'rt' ) as f: 
      temp = f.readlines()
      temp = [ x.strip() for x in temp ]
      target_ids = [ x.split("/")[0] for x in temp ]
      target_pdbs = [ x.split("/")[1] for x in temp]

  # Generate path to executable
  executable = config.rosetta_path + "mp_dock" + "." + config.platform + config.compiler + config.buildenv

  # Change directories to a data analysis dir
  outdir = config.benchmark_path + "/data/" + energy_fxn + "/" + test_name 
  if ( not os.path.isdir( outdir ) ): 
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

  for i in range(0, len(target_ids)): 

    # Setup case-specific variables (pdbfile, spanfile)
    native = targets + "/" + target_ids[i] + "/" + target_pdbs[i]
    spanfile = targets + "/" + target_ids[i] + "/" + target_ids[i] + ".span"
    prepacked = targets + "/" + target_ids[i] + "/" + target_ids + ".prepack.pdb"
    partners = target_pdbs[i].split("_")[1][0] + "_" + target_pdbs[i].split("_")[1][1]

    casedir = outdir + "/" + target_pdbs[i]
    if ( not os.path.isdir( casedir ) ): 
      os.system( "mkdir " + casedir )
      os.chdir( casedir )

    # Generate a string of arguments from the case-specific variables
    s = Template( " -in:file:s $prepacked -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners $partners -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 1000 -out:path:all $outdir -mp:lipids:composition DLPC -mp:pore:accomodate_pore true" )
    arguments = s.substitute( native=native, prepacked=prepacked, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir, partners=partners )

    # Write jobfile and submit to the HPC
    print("Submitting docking calculations for case:", target_pdbs[i] ) 
    jobname = case + "_docking"
    hpc_util.submit_condor_job( casedir, jobname, executable, arguments )
