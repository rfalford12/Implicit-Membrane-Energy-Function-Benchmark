#!/usr/bin/env python
###################################################################
#@file:         make_designed_protein_scaffold.py                                                                                   
#@description:  Generate redesigned protein scaffold                                               
#@args:     --energy_fxn (wts) --targets (list of targets)                                                          
#@author:     Rebecca F. Alford                   
#@email:    rfalford12@gmail.com                                          
###################################################################

import random
import hpc_util, read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )


def run_fixed_backbone_design_calc( energy_fxn, config, targets, test_name ): 
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
  list_of_targets = config.benchmark_path + "targets/" + targets + "/targets.list"
  with open( list_of_targets, 'rt' ) as f: 
    test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

  # Generate path to executable
  executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv

  # Change directories to a data analysis dir
  outdir = config.benchmark_path + "/data/" + energy_fxn + "/" + test_name 
  if ( not os.path.isdir( outdir ) ): 
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

      # Make a directory for this design case
      casedir = outdir + "/" + case
      if ( not os.path.isdir( casedir ) ): 
        os.system( "mkdir " + casedir )
        os.chdir( casedir )

      # Setup case-specific variables (pdbfile, spanfile, xmlargs)
      pdbfile = targets + "/" + case + "/" + case + ".pdb"
      spanfile = targets + "/" + case + "/" + case + ".span"

      # Default composition is DLPC
      s = Template( " -in:file:s $pdbfile -mp:setup:spanfiles $spanfile -score:weights $sfxn -in:membrane -out:path:all $outdir -in:file:load_PDB_components false -in:ignore_unrecognized_res" )
      arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir )

      # Write jobfile and submit to the HPC
      print("Submitting fixed backbone design case:", case ) 
      jobname = case + "_fixbb_design"
      jobfile = hpc_util.make_jobfile( casedir, case, executable, arguments )

