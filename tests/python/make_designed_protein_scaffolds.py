#!/usr/bin/env python
""" Generate redesigned protein scaffolds

This script runs fixed backbone design on native membrane
protein structures, required data for the sequence recovery
and side chain distribution tests (Test 8 and 9). 

Authors: 
  Rebecca Alford <ralford3@jhu.edu> 

Example: 
  $ import make_designed_protein_scaffolds
  $ make_designed_protein_scaffolds.run_fixed_backbone_design_calc( 
    config, energy_fxn,  targets, test_name, restore )

Arguments: 
  - config: Container with path to benchmark and rosetta files
  - energy_fxn: Weights file for energy function of interest
  - targets: Location of design targets
  - test_name: Name of test
  - restore: Use talaris2014 energy function parameters

Requirements: 
  - Rosetta release 246 or greater
"""

import random, os, sys
import hpc_util, read_config
from string import Template

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def run_fixed_backbone_design_calc( config, energy_fxn, targets, test_name, restore ): 
  """
  A function to perform fixed backbone design on a set of targets
  
  Arguments: 
    energy_fxn = energy function to use for calculations (typically, name of the weights file)
    config = path to benchmark, Rosetta executables
    targets = list of targets
    xml = Path to RosettaScript defining the peptide landscape search protocol
  """

  print( "Performing fixed backbone design on set", targets ) 

  # Read list of energy landscape test cases
  list_of_targets = config.benchmark_path + "targets/design/monomer_chains.list"
  with open( list_of_targets, 'rt' ) as f: 
    test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

  # Generate path to executable
  executable = config.rosetta_path + "fixbb" + "." + config.platform + config.compiler + config.buildenv

  # Make a directory for the subset and lipid composition
  outdir = config.benchmark_path + "data/" + energy_fxn + "/sequence-recovery" 
  os.system( "mkdir " + outdir )
  os.chdir( outdir )

  # For each test case, generate specific arguments, condor files, and then run
  for case in test_cases:

    # Make one directory per case
    casedir = outdir + "/" + case
    os.system( "mkdir " + casedir )
    os.chdir( casedir )

    # Setup arguments by substitution
    pdbfile = config.benchmark_path + "/targets/design/" + case + "/" + case + "_tr_ignorechain.pdb"
    spanfile = config.benchmark_path + "/targets/design/" + case + "/" + case + "_tr.span"
    s = Template( " -in:file:s $pdbfile -mp:setup:spanfiles $spanfile -score:weights $sfxn -in:membrane -out:path:all $outdir -in:file:load_PDB_components false -in:ignore_unrecognized_res" )
    arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir )
    if ( restore == True ): 
        arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior"
    else: 
      arguments = arguments + " -mp:lipids:composition DLPC -mp:lipids:temperature 37"

    # Write arguments and executable to a separate file
    jobfile = casedir + "/" + case + "_seqrecov.sh"
    with open( jobfile, 'a' ) as f: 
        f.write( "#!/bin/bash\n" )
        f.write( executable + " " + arguments + "\n" )
        f.close()
    os.system( "chmod +x " + jobfile )

    # Generate a condor submission file and submit the job to Jazz
    print("Submitting fixed backbone design calculation for sequence recovery case:", case) 
    queue_no = 1
    high_mem = True 
    hpc_util.submit_condor_job( casedir, case, executable, arguments, queue_no, high_mem )


